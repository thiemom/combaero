import contextlib
import math
import sys
import warnings
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterator, MutableMapping
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, Literal

import combaero as cb
from combaero import _solver_tools

# Physics configuration types for network introspectability
CompressibilityLiteral = Literal["incompressible", "compressible"]
FrictionModelLiteral = Literal["haaland", "colebrook", "serghides", "petukhov"]
HeatTransferModelLiteral = Literal[
    "none", "gnielinski", "dittus_boelter", "sieder_tate", "petukhov"
]

if TYPE_CHECKING:
    from .graph import FlowNetwork

CombustionMethodLiteral = Literal["complete", "equilibrium"]


# ============================================================================
# Convective Heat Transfer Models and Surface (Phase 3)
# ============================================================================


def _safe_rho(rho_raw: float, rho_min: float = 0.01) -> tuple[float, float]:
    """Smooth lower bound on density for numerical robustness.

    Returns (rho_safe, d_rho_safe / d_rho_raw) for Jacobian chain rule.

    Uses softplus: rho_safe = rho_min + rho_min * ln(1 + exp((rho_raw - rho_min) / rho_min))
    - When rho_raw >> rho_min: rho_safe ~= rho_raw  (transparent)
    - When rho_raw -> 0:        rho_safe -> rho_min * (1 + ln2)  (bounded)
    - When rho_raw -> -inf:     rho_safe -> rho_min  (floor)
    - Derivative is always in (0, 1]: solver always sees a gradient.
    """
    z = (rho_raw - rho_min) / rho_min
    if z > 20.0:  # overflow guard
        return rho_raw, 1.0
    if z < -20.0:  # underflow guard
        return rho_min, 0.0
    exp_z = math.exp(z)
    rho_safe = rho_min + rho_min * math.log1p(exp_z)
    d_safe = exp_z / (1.0 + exp_z)  # sigmoid - always in (0, 1)
    return rho_safe, d_safe


@dataclass
class SmoothModel:
    """Parameters for channel_smooth."""

    correlation: str = "gnielinski"  # "gnielinski" | "dittus_boelter" | "sieder_tate" | "petukhov"
    mu_ratio: float = 1.0  # mu_bulk / mu_wall (Sieder-Tate)
    roughness: float = 0.0  # absolute roughness [m]


@dataclass
class RibbedModel:
    """Rib-roughened walls, evaluated from a provenanced parameter set.

    The correlation gives the RIBBED-SIDE heat transfer. Combining it with the
    smooth walls is this element's job, because how many walls are ribbed is a
    design choice rather than a property of the correlation -- the channel is
    the internal wall by definition, and the other side of it is a different
    channel.

    Parameters
    ----------
    correlation_set : object
        A ``combaero.RibCorrelationSet``. Defaults to Han (1988) for 90 deg
        orthogonal ribs. Supply your own to match measured data: the set
        records its own source and provenance class, so a tuned coefficient
        cannot pass for a published one.
    e_D, p_e, alpha_deg : float
        Rib height / hydraulic diameter, pitch / height, and angle to the flow.
    W_H : float
        Channel width / height. It lives here rather than on the element
        because a channel is defined by its hydraulic diameter, which does not
        determine an aspect ratio -- a 2:1 duct and a square one can share a
        Dh. The correlation needs the ratio for both the wall-law geometry
        group and the ribbed/smooth area split, so the rib model carries it.
    n_ribbed_walls : int
        How many of the four walls carry ribs. 2 means two OPPOSITE walls,
        which is the configuration Han measured.
    smooth_wall_Nu_multiplier : float
        A user knob on the SMOOTH walls' contribution to the channel average.
        It encodes nothing and defaults to 1.0.

        It exists because the ribbed side already has a better knob -- change
        ``C_G`` on the parameter set, which records what you did -- while the
        smooth walls come from the base Gnielinski correlation and have none.

        Its practical use is the documented gap below. A plain smooth wall
        gives ``h_s/h_r`` around 0.42, where Han's own channel average implies
        0.70, because ribs enhance the adjacent smooth wall by 10-50% as well.
        Setting this to about **1.67** reproduces Han's measured average for
        the two-ribbed-wall square-channel case.
    """

    correlation_set: object = None
    e_D: float = 0.0
    p_e: float = 0.0
    alpha_deg: float = 90.0
    W_H: float = 1.0
    n_ribbed_walls: int = 2
    smooth_wall_Nu_multiplier: float = 1.0


@dataclass
class ImpingementModel:
    """Jet array impingement cooling for ONE spanwise row, evaluated from a
    provenanced parameter set (Florschuetz, Truman and Metzger 1981).

    A real array's rows see progressively more crossflow as spent air from
    upstream rows accumulates -- there is no single "channel Nu" for a whole
    array the way there is for a smooth or ribbed duct. Model the array by
    chaining ``n_rows`` separate elements, each with its own ``ImpingementModel``
    and ``row``, the same way a real duct is built from segments elsewhere in
    this codebase. This element does not do that chaining or the row-to-row
    mass-flow bookkeeping for you: see ``row`` below for what it assumes
    instead.

    Parameters
    ----------
    correlation_set : object
        A ``combaero.JetArrayCorrelationSet``. Defaults to
        ``florschuetz_1981_inline()``. Supply ``florschuetz_1981_staggered()``
        for a staggered hole pattern, or your own tuned set.
    d_jet : float
        Jet hole diameter [m]. This row's characteristic length for both
        ``Re_j`` and ``h = Nu * k / d_jet``.
    xn_d, yn_d, z_d : float
        Streamwise hole spacing, spanwise hole spacing, and channel height
        (jet-plate-to-target-plate gap), each normalised by ``d_jet``.
    row : int
        This row's position, 1-indexed counting from upstream. Row 1 sees
        zero crossflow by definition (``Gc/Gj = 0``, matching Florschuetz's
        own ``Nu1``); ``crossflow_to_jet_ratio_at_row`` computes the rest
        from geometry alone, so no other row's state is needed here.
    C_D : float
        Jet-plate discharge coefficient, feeding the ``Gc/Gj`` closed form.
        Defaults to the source's own recommendation absent a measured value
        (``FLORSCHUETZ_1981_DEFAULT_CD = 0.79``) -- combaero has no
        discharge-coefficient correlation of its own for a jet-plate array
        yet, see issue #375.

    What this element assumes about mass flow. The ELEMENT's mass flow is
    this row's jet flow, and ``self.area`` (on the enclosing
    ``ConvectiveSurface``) is THIS ROW's target-plate footprint, which counts
    its holes: ``n = area / (xn_d d_jet * yn_d d_jet)``, so each hole carries
    ``m_dot / n`` and ``Re_j = 4 (m_dot/n) / (pi d_jet mu)`` -- Florschuetz's
    jet mass velocity on the hole area. (Until #460 the flow was taken as
    ``rho v * area``, which made it scale with the target area.) Whether every
    row gets
    the same total flow (a uniform-supply approximation) or a row-dependent
    one (Florschuetz's own Eq. 7, deliberately not implemented -- see
    ``impingement_correlation.h``'s module comment) is up to how the caller
    assembles the chain of elements, not this model.

    What this element does NOT model. The jet plate's own orifice pressure
    loss: model that with a proper ``OrificeElement`` upstream, using a real
    discharge-coefficient correlation -- duplicating it here with a
    hand-rolled formula would risk double-counting it. Friction/pressure
    drop reported here (``f``, ``dP``) is the CROSSFLOW's own plain
    smooth-duct value: Han's book gives no friction correlation for
    impingement the way it does for ribs (there is no ``R``/``f``
    relationship in Eq. 4.9), so this is the best available proxy for the
    spent-air flow along the channel, not a claim that impingement leaves
    friction unchanged.
    """

    correlation_set: object = None
    d_jet: float = 0.0
    xn_d: float = 0.0
    yn_d: float = 0.0
    z_d: float = 0.0
    row: int = 1
    C_D: float = cb.FLORSCHUETZ_1981_DEFAULT_CD
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0


@dataclass
class SingleJetImpingementModel:
    """A single free round jet impinging on a flat plate (Goldstein,
    Behbahani and Heppelmann 1986), at a representative radial position.

    Unlike ``ImpingementModel``, there is no array and no crossflow -- this
    is one jet, one target patch. The correlation's Nu is the AVERAGE over the
    disc of radius ``R`` around the stagnation point (Han: "an average
    heat-transfer coefficient correlation"; Eq. 4.1's ``Nu_bar``), so ``R_D``
    sets the patch it describes. A ``ChannelElement`` left without a
    convective area gets exactly that disc, ``pi (R_D d_jet)^2`` (#462).

    Parameters
    ----------
    correlation_set : object
        A ``combaero.SingleJetImpingementSet``. Defaults to
        ``goldstein_1986_single_jet()``.
    bc : object
        A ``combaero.ImpingementThermalBC``. The correlation's ``(R/D)``
        exponent switches on which surface boundary condition it was fitted
        under; this does not change what the element computes for, only
        which of the source's two curves it reads from.
    d_jet : float
        Jet diameter [m].
    L_D : float
        Jet-to-target-plate spacing / ``d_jet``. Defaults to 7.75, the
        source's own optimum spacing -- a defensible default, not a claim
        about any particular rig.
    R_D : float
        Radius of the averaging disc / ``d_jet``: Nu is Goldstein's average
        over the disc of radius ``R`` around the stagnation point. No
        default: the source's closed-form check point uses 5.0, but that is
        a worked example, not a universal choice.

    Mass flow. The ELEMENT's mass flow is the jet's flow, through its one
    hole: ``Re = 4 m_dot / (pi d_jet mu)``, Goldstein's nozzle Reynolds
    number. ``self.area`` is the target patch the reported h applies to; it
    does not set the flow. (Until #460 it did, as ``rho v * area``.)

    Pressure drop is not modelled here for the same reason as
    ``ImpingementModel``: model the nozzle's own loss with a proper
    ``OrificeElement`` upstream, not a formula duplicated inside this
    heat-transfer model. ``f``/``dP`` reported here are the target-side
    channel's own plain smooth-duct values, borrowed for a T_aw computation
    Han's Eq. 4.2/4.3 recovery-factor model would otherwise have to supply
    -- see ``han_impingement.md`` items 6-7, not yet wired in.
    """

    correlation_set: object = None
    bc: object = cb.ImpingementThermalBC.ConstantHeatFlux
    d_jet: float = 0.0
    L_D: float = 7.75
    R_D: float = 0.0
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0


@dataclass
class PinFinModel:
    """A pin-fin array, evaluated from provenanced parameter sets (#335).

    Heat transfer and friction are separate sets, each in its source's own
    form, evaluated in one canonical basis: ``Re_D`` on the pin diameter and
    the velocity at the minimum flow area, ``Nu_D = h D / k``, and
    ``f = dP / (2 rho Vmax^2 N)``. See ``pin_fin_correlation.h``.

    Parameters
    ----------
    nu_set : object
        A ``combaero.PinFinNuSet`` on the Total surface. Defaults to
        ``metzger_1986_staggered_nu()`` for staggered arrays and
        ``chyu_1998_nu(Inline, Total)`` for inline ones.
    f_set : object
        A ``combaero.PinFinFrictionSet``. Defaults to
        ``metzger_1982_staggered_friction()`` for staggered arrays and
        ``chyu_1990_friction(Inline)`` for inline ones -- the only inline
        friction source in hand, fitted to its digitised Fig. 6 at one
        geometry (H/D 1, S/D = X/D = 2.5). Staggered friction is never
        substituted for inline.
    modifier : object
        An optional ``combaero.PinFinRatioModifier`` applied to ``nu_set``
        (and to ``f_set`` when it carries a friction ratio), e.g.
        ``chyu_1998_inline_over_staggered()`` to transfer a staggered set to
        an inline array. Its ``from_arrangement`` must match ``nu_set``'s and
        its ``to_arrangement`` this array's.
    pin_diameter : float
        Pin diameter D [m].
    S_D, X_D, H_D : float
        Transverse pitch, streamwise pitch and pin height (= channel height),
        each over D.
    N_rows : int
        Pin rows along the flow. Sets the pressure drop and is checked against
        the sets' own row counts (no row correction is applied yet).
    arrangement : object
        ``combaero.PinArrangement``.
    k_pin : float
        Pin thermal conductivity [W/(m K)] for the fin efficiency. 20 is a
        nickel superalloy at turbine metal temperatures; ``math.inf`` treats
        the pins as isothermal.

    Area and velocity conventions. ``ConvectiveSurface.area`` is the BASE
    (planform) area of the endwall this surface couples to; the effective
    coefficient returned is referenced to it,
    ``h_eff = h (A_endwall_exposed + eta_fin A_pin_half) / A_base``, with each
    wall feeding its pins to mid-height. The element's flow area is the
    UNOBSTRUCTED channel cross-section; Vmax follows from it exactly through
    the array geometry.

    The defaults are a working example inside the Metzger sets' box
    (5 mm pins, H/D 1, S/D = X/D = 2.5, 10 rows, Inconel-class pins), not a
    generic "typical" pin fin.
    """

    nu_set: object = None
    f_set: object = None
    modifier: object = None
    pin_diameter: float = 0.005
    S_D: float = 2.5
    X_D: float = 2.5
    H_D: float = 1.0
    N_rows: int = 10
    arrangement: object = cb.PinArrangement.Staggered
    k_pin: float = 20.0

    def __post_init__(self) -> None:
        if not (self.pin_diameter > 0.0):
            raise ValueError("PinFinModel: pin_diameter must be positive [m]")
        cb.validate_pin_fin_geometry(self.geometry())
        f_arr = self.resolved_f_set().arrangement
        if f_arr != self.arrangement:
            raise ValueError(
                f"PinFinModel: f_set {self.resolved_f_set().name!r} is "
                f"{f_arr.name}, but this array is {self.arrangement.name}. "
                "Staggered and inline friction differ by about 1.6x (Chyu "
                "1990, Fig. 6) and are not substituted for each other."
            )
        if self.modifier is not None:
            nu_arr = self.resolved_nu_set().arrangement
            if self.modifier.from_arrangement != nu_arr:
                raise ValueError(
                    f"PinFinModel: modifier {self.modifier.name!r} transfers from "
                    f"{self.modifier.from_arrangement.name}, but nu_set is "
                    f"{nu_arr.name}."
                )
            if self.modifier.to_arrangement != self.arrangement:
                raise ValueError(
                    f"PinFinModel: modifier {self.modifier.name!r} transfers to "
                    f"{self.modifier.to_arrangement.name}, but this array is "
                    f"{self.arrangement.name}."
                )

    def geometry(self):
        return cb.PinFinGeometry(
            S_D=self.S_D,
            X_D=self.X_D,
            H_D=self.H_D,
            N_rows=self.N_rows,
            arrangement=self.arrangement,
        )

    def resolved_nu_set(self):
        if self.nu_set is not None:
            return self.nu_set
        if self.modifier is not None:
            # A transfer modifier starts from its own source arrangement.
            if self.modifier.from_arrangement == cb.PinArrangement.Staggered:
                return cb.metzger_1986_staggered_nu()
            return cb.chyu_1998_nu(self.modifier.from_arrangement, cb.PinNuSurface.Total)
        if self.arrangement == cb.PinArrangement.Inline:
            return cb.chyu_1998_nu(cb.PinArrangement.Inline, cb.PinNuSurface.Total)
        return cb.metzger_1986_staggered_nu()

    def resolved_f_set(self):
        if self.f_set is not None:
            return self.f_set
        if self.arrangement == cb.PinArrangement.Inline:
            return cb.chyu_1990_friction(cb.PinArrangement.Inline)
        return cb.metzger_1982_staggered_friction()


@dataclass
class _ImpingementChannelResult:
    """What either impingement model returns.

    Deliberately NOT a ``ChannelResult``: same rationale as
    ``_RibbedChannelResult`` -- this path computes its own derivatives and
    (for the array case) an ``extrapolated`` flag the C++ correlation
    already carries.
    """

    h: float
    Nu: float
    Re: float
    Pr: float
    f: float
    dP: float
    T_aw: float
    extrapolated: bool = False
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 0.0
    dT_aw_dP: float = 0.0


@dataclass
class _PinFinChannelResult:
    """What a pin-fin array returns.

    ``h`` is the effective coefficient on the BASE area (fin efficiency
    included); ``h_array`` is the correlation's own coefficient on the total
    pin + endwall area, before the fin model, so both halves of that
    modelling step are visible.
    """

    h: float
    h_array: float
    eta_fin: float
    Nu: float
    Re: float
    Pr: float
    f: float
    dP: float
    T_aw: float
    extrapolated: bool = False
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 0.0
    dT_aw_dP: float = 0.0


@dataclass
class _RibbedChannelResult:
    """What a ribbed channel returns.

    Deliberately NOT a ``ChannelResult``: that type is the C++ correlation's
    output and carries derivative fields this path computes differently. It
    also exposes ``h_ribbed`` and ``h_smooth`` separately, because the channel
    average hides a 20% modelling choice and a user should be able to see both
    halves of it.
    """

    h: float
    h_ribbed: float
    h_smooth: float
    Nu: float
    Re: float
    Pr: float
    f: float
    dP: float
    T_aw: float
    e_plus: float
    extrapolated: bool
    # Derivatives the solver's wall-coupling path reads. They are not
    # decoration: without them a ribbed channel joined by a ThermalWall raises
    # AttributeError mid-solve, which is how the coupled combustor example
    # found this after the unit tests missed it -- none of them coupled a wall.
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 0.0
    dT_aw_dP: float = 0.0


ChannelModel = (
    SmoothModel | RibbedModel | ImpingementModel | SingleJetImpingementModel | PinFinModel
)


@dataclass
class ConvectiveSurface:
    """Convective surface description attached to a NetworkElement.

    Attributes
    ----------
    area : float
        Wetted convective area [m^2]. Default 0.0 (disabled).
    model : ChannelModel
        Model-specific subclass holding geometry parameters.
    heating : bool | None
        True if fluid is being heated, False if cooled, None = auto-detect
        from sign of (T_hot - T_fluid).
    Nu_multiplier : float
        Empirical correction factor on Nusselt number (default 1.0).
    f_multiplier : float
        Empirical correction factor on friction factor (default 1.0).
    """

    area: float = 0.0  # A_conv [m^2] - 0 disables
    model: ChannelModel = field(default_factory=SmoothModel)
    heating: bool | None = None  # None = auto-detect
    Nu_multiplier: float = 1.0  # empirical correction on Nu
    f_multiplier: float = 1.0  # empirical correction on f
    # The calling element's flow area for the current evaluation (see
    # htc_and_T); not a parameter of the surface.
    _flow_area: float | None = field(default=None, init=False, repr=False, compare=False)

    def _channel_mdot(self, rho: float, velocity: float, diameter: float) -> float:
        """The solver's mass flow for this surface: through the CHANNEL
        cross-section -- the element's own flow area when it passed one
        (``flow_area`` on ``htc_and_T``), else ``pi Dh^2 / 4`` as
        ``channel_smooth`` assumes. A non-circular channel given a separate
        ``Dh`` needs the real area: ``pi Dh^2/4`` is then off by
        ``(Dh/diameter)^2`` (#460).

        Every ``dh_dmdot`` this class returns is a derivative with respect to
        the element's own ``m_dot`` -- that is what the solver's wall coupling
        multiplies it by -- and ``ChannelElement.htc_and_T`` derives the
        velocity from that mass flow and the channel area, the same convention
        ``channel_smooth`` uses. The convective ``area`` is a heat-transfer
        area, not a flow area. Using it here scaled the ribbed and impingement
        derivatives by ``A_cross / A_surface`` (#456).
        """
        area = self._flow_area if self._flow_area else math.pi / 4.0 * diameter * diameter
        return rho * abs(velocity) * area

    def _flow_area_arg(self) -> float:
        """The element's flow area for channel_smooth's mass-flow derivatives,
        NaN when it passed none (channel_smooth then assumes pi D^2/4, #463)."""
        return self._flow_area if getattr(self, "_flow_area", None) else math.nan

    def _ribbed_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """Ribbed-channel heat transfer and pressure drop.

        Two things are asymmetric here and both follow from the source rather
        than from convenience.

        FRICTION NEEDS NO WALL WEIGHTING. The ``f`` in the roughness function's
        definition is already the equivalent four-sided channel friction
        factor, so the correlation returns a channel-level value directly. The
        element does not multiply pipe friction by anything -- correlations own
        their ``f``, which is what #331 established after a round-trip through
        a restated friction factor moved a drop by -24.7%.

        HEAT TRANSFER DOES. The correlation gives the ribbed side; the smooth
        walls come from the base correlation, and the channel average is the
        area-weighted combination.

        A DOCUMENTED GAP. Using the plain smooth correlation for the smooth
        walls gives h_s/h_r around 0.42, where Han's own reported channel
        average implies 0.70 -- because ribs enhance the adjacent smooth wall
        by 10-50% as well, which no correlation here covers. The channel
        average is therefore about 20% below Han's measurement for the
        two-ribbed-wall square case. That is a knowable, quantified
        under-prediction rather than an invented constant, and
        ``smooth_wall_Nu_multiplier`` is how a user closes it: about 1.67
        reproduces Han.
        """
        model = self.model
        rib_set = model.correlation_set or cb.han_1988_orthogonal()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        Re = rho * velocity * diameter / mu if mu > 0.0 else 0.0
        Pr = cs.transport.Pr

        geom = cb.RibGeometry(
            e_D=model.e_D,
            p_e=model.p_e,
            W_H=model.W_H,
            alpha_deg=model.alpha_deg,
        )
        rib = cb.evaluate_rib(rib_set, geom, Re)

        # Ribbed side, from the correlation's Stanton number.
        k = cs.transport.k
        cp = cs.thermo.cp
        h_ribbed = rib.St_r * rho * abs(velocity) * cp

        # Smooth walls, from the base correlation, with the user's knob.
        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=model.smooth_wall_Nu_multiplier,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        frac_ribbed, frac_smooth = self._ribbed_wall_fractions()
        h_avg = frac_ribbed * h_ribbed + frac_smooth * smooth.h

        # The correlation's f is already the four-sided channel value.
        dP = rib.f * (length / diameter) * 0.5 * rho * velocity * abs(velocity)

        # Wall-coupling derivatives. The ribbed side is h_r = St(Re) rho|v| cp,
        # and BOTH Re and |v| are proportional to the mass flow, so
        #     dh_r/dmdot = (Re dSt/dRe + St) rho |v| cp / mdot.
        # The St term dominates (h_r grows roughly as mdot^0.8); it was missing,
        # and the mass flow was taken on the convective area (#456). The
        # smooth side brings its own derivative from the base correlation;
        # both are area-weighted exactly as h is.
        mdot = self._channel_mdot(rho, velocity, diameter)
        dh_ribbed_dmdot = (
            (Re * rib.dSt_dRe + rib.St_r) * rho * abs(velocity) * cp / mdot if mdot else 0.0
        )
        dh_dmdot = (
            frac_ribbed * dh_ribbed_dmdot + frac_smooth * smooth.dh_dmdot
        ) * self.Nu_multiplier
        # Temperature sensitivity of the ribbed side is not exposed by the
        # correlation, so only the smooth side contributes. Recorded rather
        # than approximated from a stand-in, the same gap the array path had.
        dh_dT = frac_smooth * smooth.dh_dT * self.Nu_multiplier

        return _RibbedChannelResult(
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
            h=h_avg * self.Nu_multiplier,
            h_ribbed=h_ribbed,
            h_smooth=smooth.h,
            Nu=h_avg * diameter / k if k > 0.0 else 0.0,
            Re=Re,
            Pr=Pr,
            f=rib.f * self.f_multiplier,
            dP=dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            e_plus=rib.e_plus,
            extrapolated=rib.extrapolated,
        )

    def _ribbed_wall_fractions(self) -> tuple[float, float]:
        """Area fractions of the ribbed and smooth walls, for a W x H duct.

        Han's configuration is two OPPOSITE walls, which for a duct of width W
        and height H are the two of width W. The ribbed fraction is then
        W/(W+H) -- the same area weighting the four-sided friction conversion
        uses, so friction and heat transfer are treated alike.

        One and four ribbed walls follow from the same picture. Anything else
        is rejected rather than interpolated: "three ribbed walls" has no
        unambiguous geometry.
        """
        model = self.model
        W_H = model.W_H
        n = model.n_ribbed_walls
        if n == 4:
            return 1.0, 0.0
        if n == 2:
            ribbed = W_H / (W_H + 1.0)
            return ribbed, 1.0 - ribbed
        if n == 1:
            ribbed = W_H / (2.0 * (W_H + 1.0))
            return ribbed, 1.0 - ribbed
        raise ValueError(
            f"n_ribbed_walls must be 1, 2 or 4, got {n}. 2 means two opposite "
            "walls, which is the configuration the correlations were measured "
            "on."
        )

    def _impingement_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """Jet array impingement, one row (Florschuetz, Truman and Metzger 1981).

        See ``ImpingementModel``'s docstring for the mass-flow and
        pressure-drop conventions this follows.
        """
        model = self.model
        jet_set = model.correlation_set or cb.florschuetz_1981_inline()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        k = cs.transport.k
        Pr = cs.transport.Pr

        # The element's mass flow IS this row's jet flow (Florschuetz: G_j is
        # the jet mass velocity on the hole area). The convective area only
        # counts the holes. It used to supply the flow as well, as
        # rho v * area -- wrong by A_surface/A_cross (#460).
        mdot_total = self._channel_mdot(rho, velocity, diameter)
        hole_footprint = model.xn_d * model.d_jet * model.yn_d * model.d_jet
        n_holes = self.area / hole_footprint if hole_footprint > 0.0 else 0.0
        mdot_per_hole = mdot_total / n_holes if n_holes > 0.0 else 0.0
        hole_area = math.pi / 4.0 * model.d_jet**2
        v_jet = mdot_per_hole / (rho * hole_area) if rho * hole_area > 0.0 else 0.0
        Re_j = rho * v_jet * model.d_jet / mu if mu > 0.0 else 0.0

        Gc_Gj = cb.crossflow_to_jet_ratio_at_row(model.yn_d, model.z_d, model.C_D, model.row)
        jet = cb.jet_array_impingement_nu(
            jet_set, Re_j, Gc_Gj, Pr, model.xn_d, model.yn_d, model.z_d
        )
        h = jet.Nu * k / model.d_jet if model.d_jet > 0.0 else 0.0

        # T_aw, f and dP are not covered by Eq. 4.9 at all -- Han's book gives
        # no impingement-specific friction/recovery model, unlike ribs' own
        # R/f relationship. Borrowed wholesale from the plain smooth
        # correlation applied to the crossflow's own bulk state, per
        # ImpingementModel's docstring.
        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=1.0,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        # Re_j is proportional to the element's mass flow; the per-hole split
        # above is this row's own bookkeeping, not the solver's m_dot (#456).
        mdot = self._channel_mdot(rho, velocity, diameter)
        dRe_j_dmdot = abs(Re_j / mdot) if mdot else 0.0
        dh_dmdot = jet.dNu_dRe_j * dRe_j_dmdot * (k / model.d_jet if model.d_jet > 0.0 else 0.0)
        dh_dmdot *= self.Nu_multiplier
        # Temperature sensitivity of the correlation itself is not exposed;
        # only the borrowed smooth-side contributes, the same documented gap
        # ribs accepted for their own T-sensitivity.
        dh_dT = smooth.dh_dT * self.Nu_multiplier

        return _ImpingementChannelResult(
            h=h * self.Nu_multiplier,
            Nu=jet.Nu,
            Re=Re_j,
            Pr=Pr,
            f=smooth.f * self.f_multiplier,
            dP=smooth.dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            extrapolated=jet.extrapolated,
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
        )

    def _single_jet_impingement_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """A single free jet (Goldstein, Behbahani and Heppelmann 1986).

        See ``SingleJetImpingementModel``'s docstring for the mass-flow and
        pressure-drop conventions this follows.
        """
        model = self.model
        jet_set = model.correlation_set or cb.goldstein_1986_single_jet()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        k = cs.transport.k

        # The element's mass flow IS the jet's flow (Goldstein: nozzle Re on
        # the jet's own exit area), not rho v * the target-patch area (#460).
        mdot_total = self._channel_mdot(rho, velocity, diameter)
        hole_area = math.pi / 4.0 * model.d_jet**2
        v_jet = mdot_total / (rho * hole_area) if rho * hole_area > 0.0 else 0.0
        Re = rho * v_jet * model.d_jet / mu if mu > 0.0 else 0.0

        jet = cb.single_jet_impingement(jet_set, model.bc, Re, model.L_D, model.R_D)
        h = jet.Nu * k / model.d_jet if model.d_jet > 0.0 else 0.0

        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=1.0,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        mdot = self._channel_mdot(rho, velocity, diameter)
        dRe_dmdot = abs(Re / mdot) if mdot else 0.0
        dh_dmdot = jet.dNu_dRe * dRe_dmdot * (k / model.d_jet if model.d_jet > 0.0 else 0.0)
        dh_dmdot *= self.Nu_multiplier
        dh_dT = smooth.dh_dT * self.Nu_multiplier

        return _ImpingementChannelResult(
            h=h * self.Nu_multiplier,
            Nu=jet.Nu,
            Re=Re,
            Pr=cs.transport.Pr,
            f=smooth.f * self.f_multiplier,
            dP=smooth.dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
        )

    def _pin_fin_terms(self, Re: float, Pr: float):
        """Nu and f (with modifier) at canonical Re_D, and their Re-slopes."""
        model = self.model
        geom = model.geometry()
        nu = cb.evaluate_pin_fin_nu(model.resolved_nu_set(), geom, Re, Pr)
        fr = cb.evaluate_pin_fin_friction(model.resolved_f_set(), geom, Re)
        Nu, dNu = nu.Nu, nu.dNu_dRe
        f, df = fr.f, fr.df_dRe
        extrapolated = nu.extrapolated or fr.extrapolated
        if model.modifier is not None:
            m = cb.evaluate_pin_fin_modifier(model.modifier, geom, Re)
            Nu, dNu = Nu * m.ratio_Nu, dNu * m.ratio_Nu + Nu * m.dratio_Nu_dRe
            if m.has_f:
                f, df = f * m.ratio_f, df * m.ratio_f + f * m.dratio_f_dRe
            extrapolated = extrapolated or m.extrapolated
        return Nu, dNu, f, df, extrapolated

    def _pin_fin_result(self, T, P, X, velocity, diameter, length, T_hot, heating):
        """Pin-fin array: Metzger et al. (and others) via pin_fin_correlation.h.

        ``velocity`` is the bulk velocity through the UNOBSTRUCTED channel
        cross-section; Vmax = velocity / (A_min/A_frontal), exactly. The
        correlation's h lives on the total pin + endwall area; the fin model
        turns it into ``h`` on this surface's BASE area:

            h_eff = h (f_endwall + eta_fin(h) f_pin)

        with the per-wall fractions of ``pin_fin_area_fractions``.

        ``dh_dmdot`` uses the element's mass flow through the channel
        cross-section, ``rho v pi Dh^2 / 4`` -- the same convention as
        ``channel_smooth``. T_aw and its derivatives are borrowed from the
        smooth correlation, as the impingement paths do; ``dh_dT`` is the
        smooth result's relative property sensitivity scaled to ``h_eff``, a
        documented proxy: the pin correlations expose no temperature
        derivative of their own.
        """
        model = self.model
        D = model.pin_diameter
        geom = model.geometry()

        rho, _ = _safe_rho(cb.density(T, P, X))
        cs = cb.complete_state(T, P, X)
        mu = cs.transport.mu
        k = cs.transport.k
        Pr = cs.transport.Pr

        amin = cb.pin_fin_amin_over_afrontal(geom)
        v_max = velocity / amin
        Re = rho * v_max * D / mu if mu > 0.0 else 0.0

        Nu, dNu, f, _, extrapolated = self._pin_fin_terms(Re, Pr)
        h_array = Nu * k / D
        dh_array_dRe = dNu * k / D

        frac = cb.pin_fin_area_fractions(geom)
        eff = cb.pin_fin_array_efficiency(
            h_array, model.k_pin, D, model.H_D * D, frac.pin_over_total
        )
        h_eff = h_array * (frac.endwall_exposed + eff.eta_fin * frac.pin)
        dh_eff_dh = frac.endwall_exposed + frac.pin * (eff.eta_fin + h_array * eff.deta_fin_dh)

        smooth = cb.channel_smooth(
            T,
            P,
            X,
            velocity,
            diameter,
            length,
            T_hot=T_hot,
            heating=heating,
            Nu_multiplier=1.0,
            f_multiplier=1.0,
            flow_area=self._flow_area_arg(),
        )

        mdot = self._channel_mdot(rho, velocity, diameter)
        dRe_dmdot = abs(Re / mdot) if mdot else 0.0
        dh_dmdot = dh_eff_dh * dh_array_dRe * dRe_dmdot * self.Nu_multiplier
        dh_dT = smooth.dh_dT * (h_eff / smooth.h) * self.Nu_multiplier if smooth.h > 0.0 else 0.0

        dP = 2.0 * rho * v_max * abs(v_max) * model.N_rows * f

        return _PinFinChannelResult(
            h=h_eff * self.Nu_multiplier,
            h_array=h_array,
            eta_fin=eff.eta_fin,
            Nu=Nu,
            Re=Re,
            Pr=Pr,
            f=f * self.f_multiplier,
            dP=dP * self.f_multiplier,
            T_aw=smooth.T_aw,
            extrapolated=extrapolated,
            dh_dmdot=dh_dmdot,
            dh_dT=dh_dT,
            dT_aw_dmdot=smooth.dT_aw_dmdot,
            dT_aw_dT=smooth.dT_aw_dT,
            dT_aw_dP=smooth.dT_aw_dP,
        )

    def htc_and_T(
        self,
        T: float,
        P: float,
        X: list[float],
        velocity: float,
        diameter: float,
        length: float,
        T_hot: float = math.nan,
        flow_area: float | None = None,
    ):
        """Compute heat transfer coefficient and adiabatic wall temperature.

        Parameters
        ----------
        T : float
            Bulk static temperature [K].
        P : float
            Bulk static pressure [Pa].
        X : list[float]
            Mole fractions [-].
        velocity : float
            Bulk flow velocity [m/s].
        diameter : float
            Hydraulic diameter [m].
        length : float
            Channel length [m].
        T_hot : float, optional
            Wall temperature [K]. Used for auto-detection of heating/cooling.
        flow_area : float, optional
            The element's flow cross-section [m^2], from which it computed
            ``velocity``. Lets the surface recover the element's mass flow
            exactly; without it ``pi diameter^2 / 4`` is assumed.

        Returns
        -------
        ChannelResult | None
            Full ChannelResult with h, T_aw, and Jacobians (dh_dmdot, dh_dT, etc.),
            or None if area=0. Access convective area via ``self.area``.
        """

        if self.area == 0.0 or abs(velocity) < 1e-12:
            return None
        self._flow_area = flow_area

        # Auto-detect heating direction from T_hot - T sign
        if self.heating is not None:
            heating = self.heating
        elif math.isfinite(T_hot):
            heating = T_hot >= T
        else:
            heating = True  # default when T_hot unknown

        # Dispatch to appropriate C++ channel_* function based on model type
        if isinstance(self.model, SmoothModel):
            result = cb.channel_smooth(
                T,
                P,
                X,
                velocity,
                diameter,
                length,
                T_hot=T_hot,
                correlation=self.model.correlation,
                heating=heating,
                mu_ratio=self.model.mu_ratio,
                roughness=self.model.roughness,
                Nu_multiplier=self.Nu_multiplier,
                f_multiplier=self.f_multiplier,
                flow_area=self._flow_area_arg(),
            )
        elif isinstance(self.model, RibbedModel):
            result = self._ribbed_result(T, P, X, velocity, diameter, length, T_hot, heating)
        elif isinstance(self.model, ImpingementModel):
            result = self._impingement_result(T, P, X, velocity, diameter, length, T_hot, heating)
        elif isinstance(self.model, SingleJetImpingementModel):
            result = self._single_jet_impingement_result(
                T, P, X, velocity, diameter, length, T_hot, heating
            )
        elif isinstance(self.model, PinFinModel):
            result = self._pin_fin_result(T, P, X, velocity, diameter, length, T_hot, heating)
        else:
            raise TypeError(
                f"Unsupported channel model {type(self.model).__name__}. "
                "Enhanced-surface correlations were removed in 0.7.0 pending "
                "provenanced replacements; see issue #339."
            )

        return result


# ============================================================================
# WallConnection (Phase 4/5)
# ============================================================================


@dataclass
class WallLayer:
    """Description of a single physical layer in a multi-layer ThermalWall.

    Attributes
    ----------
    thickness : float
        Thickness of the layer [m].
    conductivity : float
        Thermal conductivity of the material [W/(m*K)].
    material : str
        Material name (e.g., 'inconel718', 'haynes230', or 'custom').
    """

    thickness: float
    conductivity: float
    material: str = "custom"

    def update_conductivity(self, T: float) -> None:
        """Update conductivity based on temperature if a database material is selected."""

        # 'generic' is legacy, 'custom' is new name for manual input
        if self.material.lower() not in ("generic", "custom"):
            with contextlib.suppress(RuntimeError, ValueError):
                # The C++ binding handles T-clamping and k-clamping internally
                self.conductivity = cb.get_material_conductivity(self.material, T)

    @property
    def r_val(self) -> float:
        """Thermal resistance per unit area [m^2*K/W]."""
        # Defensive k-clamping to prevent division by zero/near-zero
        k_safe = max(1e-3, self.conductivity)
        return self.thickness / k_safe


@dataclass
class ThermalWall:
    """Thermal coupling between two elements through a shared multi-layer wall.

    Attributes
    ----------
    id : str
        Unique identifier for this wall connection.
    element_a : str
        Element ID for side A.
    element_b : str
        Element ID for side B.
    layers : list[WallLayer]
        Stack of wall layers ordered from side A to side B.
    contact_area : float | None
        Override effective area [m^2]. If None, uses min(A_conv_a, A_conv_b).
    R_fouling : float
        Additional fouling resistance [m^2*K/W] (applied to cold side).
    """

    id: str
    element_a: str
    element_b: str
    layers: list[WallLayer] = field(default_factory=list)
    contact_area: float | None = None
    R_fouling: float = 0.0

    # Internal cache for temperature-dependent property iteration
    _last_profile: list[float] | None = field(default=None, init=False, repr=False)

    def compute_coupling(
        self,
        h_a: float,
        T_aw_a: float,
        A_conv_a: float,
        h_b: float,
        T_aw_b: float,
        A_conv_b: float,
    ) -> cb.WallCouplingResult:
        """Compute heat transfer rate and wall temperature using analytical Jacobians.

        Parameters
        ----------
        h_a, h_b : float
            Convective HTCs on sides A and B [W/(m^2*K)].
        T_aw_a, T_aw_b : float
            Adiabatic wall temperatures on sides A and B [K].
        A_conv_a, A_conv_b : float
            Convective areas on both sides [m^2].

        Returns
        -------
        WallCouplingResult
            Object containing Q [W], T_hot [K], and analytical Jacobians.
        """
        A_eff = self.contact_area if self.contact_area is not None else min(A_conv_a, A_conv_b)

        # The profile only places the layer temperatures; it refuses h <= 0,
        # which a correlation past its range can return at a solver iterate,
        # so it sees the coupling's own floor (WALL_HTC_KNEE).
        h_pa = max(float(h_a), cb.WALL_HTC_KNEE)
        h_pb = max(float(h_b), cb.WALL_HTC_KNEE)

        def profile(t_over_k: list[float]) -> list[float]:
            prof, _q = cb.wall_temperature_profile(
                T_aw_a, T_aw_b, h_pa, h_pb, t_over_k, self.R_fouling
            )
            return [float(tp) for tp in prof]

        t_over_k_layers = [L.r_val for L in self.layers]
        prof = profile(t_over_k_layers)
        if any(L.material.lower() not in ("generic", "custom") for L in self.layers):
            # k(T) at THIS call's wall temperatures, iterated to its fixed
            # point, so Q depends on the inputs alone. It used to be lagged
            # from the previous call -- any trial point the solver probed --
            # which made the residual depend on evaluation history (#481).
            for _ in range(30):
                for i, layer in enumerate(self.layers):
                    layer.update_conductivity(0.5 * (prof[i] + prof[i + 1]))
                new = [L.r_val for L in self.layers]
                change = max(
                    abs(a - b) / max(abs(b), 1e-300)
                    for a, b in zip(new, t_over_k_layers, strict=True)
                )
                t_over_k_layers = new
                prof = profile(t_over_k_layers)
                if change < 1e-12:
                    break

        res = _solver_tools.wall_coupling_and_jacobian_multilayer(
            h_a, T_aw_a, h_b, T_aw_b, t_over_k_layers, A_eff, self.R_fouling
        )
        self._last_profile = prof
        return res


class EnergyBoundary:
    """Energy source/sink that attaches to mixing nodes.

    Convention: Q > 0 / fraction > 0 means heat **into** the fluid (heating),
                Q < 0 / fraction < 0 means heat **out of** the fluid (cooling).

    Two additive modes:
      - Q [W]: absolute heat transfer rate
      - fraction [-]: relative to total enthalpy flow (H_tot = sum(mdot_i * h_i))
        e.g. fraction=-0.05 means 5% enthalpy loss

    Both can be set simultaneously; effective Q = Q + fraction * H_tot.
    C++ functions handle the conversion to delta_h and all Jacobian corrections.
    """

    def __init__(self, id: str, Q: float = 0.0, fraction: float = 0.0) -> None:
        self.id = id
        self.Q = Q  # [W]
        self.fraction = fraction  # [-]


@dataclass
class NetworkMixtureState:
    P: float
    Pt: float
    T: float
    Tt: float
    m_dot: float
    Y: list[float]

    @property
    def X(self) -> list[float]:
        """Cached mole fractions converted from mass fractions with defensive clamping."""
        try:
            return self._X  # type: ignore[return-value]
        except AttributeError:
            # Defensive guard against unphysical intermediate solver states
            Y_safe = [max(1e-12, min(1.0, float(yi))) for yi in self.Y]
            y_sum = sum(Y_safe)
            if y_sum > 0:
                Y_safe = [yi / y_sum for yi in Y_safe]
            x = list(cb.mass_to_mole(Y_safe))
            object.__setattr__(self, "_X", x)
            return x

    def density(self) -> float:
        """Static density [kg/m^3]."""
        return cb.density(self.T, self.P, self.X)

    def enthalpy(self) -> float:
        """Static specific enthalpy [J/kg]."""
        return cb.h_mass(self.T, self.X)

    def total_enthalpy(self) -> float:
        """Total (stagnation) specific enthalpy [J/kg]."""
        return cb.h_mass(self.Tt, self.X)

    def cp(self) -> float:
        """Specific heat capacity at constant pressure [J/(kg*K)]."""
        return cb.cp_mass(self.T, self.X)

    def speed_of_sound(self) -> float:
        """Speed of sound [m/s]."""
        return cb.speed_of_sound(self.T, self.X)

    def gamma(self) -> float:
        """Ratio of specific heats (cp/cv) [-]."""
        cp = self.cp()
        cv = cb.cv_mass(self.T, self.X)
        return cp / cv if cv > 0 else 1.4


#: Nominal cross-section for a chamber-like node whose area could not be
#: inferred and which nothing in the network measures against. It is a
#: bookkeeping placeholder, not a geometry: components for which the area IS
#: the physics (AreaChangeElement, ChannelElement, TeeJunctionElement) refuse
#: to default and raise instead. Any consumer of a defaulted area reports
#: ``area_source = "default"`` in its diagnostics (issue #262).
DEFAULT_CHAMBER_AREA = 0.1


# ---------------------------------------------------------------------------
# Flow-area inference
# ---------------------------------------------------------------------------
# Elements and nodes that need a reference flow area infer one from the
# topology when the user did not supply it.  Historically every such site
# searched for a neighbouring ChannelElement only and, failing that, fell back
# to a hard-coded constant (0.01 / 0.02 / 0.1 m^2).  Those constants trace to
# "something had to be there", not to any geometry: with a non-channel
# neighbour -- a momentum chamber, a combustor, an ejector outlet -- the
# fallback fabricated an area change that never existed, and the resulting
# residual is unsatisfiable at the real mass flow, so the solver stalls with a
# message that blames the solver rather than the geometry (issue #262).
#
# _infer_flow_area consults every area-bearing neighbour instead.  Channels
# keep the highest priority so that networks which resolved before resolve to
# the same area now.


def _known_node_area(node: "NetworkNode | None") -> float | None:
    """Flow area of a node, but only when it is genuinely known.

    A node whose area was auto-sized rather than user-supplied
    (``_auto_area``) carries a placeholder, not a measurement; inferring from
    it would launder one guess into another.
    """
    if node is None or getattr(node, "_auto_area", True):
        return None
    area = getattr(node, "area", None)
    return float(area) if area is not None and area > 0.0 else None


def _element_area_facing(elem: "NetworkElement", node_id: str) -> float | None:
    """Flow area the element presents to ``node_id``, or None if it has none.

    Orifices are deliberately excluded: a bore is a restriction, not the
    channel area a neighbour should inherit.
    """
    if isinstance(elem, ChannelElement):
        if elem.diameter is not None:
            return math.pi * (elem.diameter / 2.0) ** 2
        return None
    if isinstance(elem, AreaChangeElement):
        # The element changes area across itself, so the face matters.
        if elem.to_node == node_id:
            return float(elem.F1) if elem.F1 is not None else None
        if elem.from_node == node_id:
            return float(elem.F0) if elem.F0 is not None else None
        return None
    if isinstance(elem, TeeJunctionElement):
        if elem.F_C is None:
            return None
        # Branch arm carries F_C/psi; straight and common arms carry F_C.
        if node_id == elem.branch_node:
            return float(elem.F_C) / elem.psi
        return float(elem.F_C)
    if isinstance(elem, (BorderCarnotLossElement, PressureLossElement)):
        return float(elem.area) if elem.area is not None else None
    return None


def _infer_flow_area(
    graph: "FlowNetwork",
    node_ids: "list[str]",
    exclude: "NetworkElement | None" = None,
) -> "tuple[float | None, str]":
    """Infer a flow area from the neighbourhood of ``node_ids``.

    Returns ``(area, provenance)``; ``area`` is None when nothing in the
    neighbourhood knows one, and ``provenance`` names the source so callers
    can report where a resolved area came from.

    Channels are searched first across all of ``node_ids`` so that a network
    which previously resolved from a channel still resolves to that channel's
    area regardless of what else is attached.
    """
    neighbours: list[tuple[str, NetworkElement]] = []
    for node_id in node_ids:
        for elem in graph.get_upstream_elements(node_id) + graph.get_downstream_elements(node_id):
            if elem is exclude:
                continue
            neighbours.append((node_id, elem))

    for node_id, elem in neighbours:
        if isinstance(elem, ChannelElement):
            area = _element_area_facing(elem, node_id)
            if area is not None:
                return area, f"channel '{elem.id}'"

    for node_id, elem in neighbours:
        area = _element_area_facing(elem, node_id)
        if area is not None:
            return area, f"{type(elem).__name__} '{elem.id}'"

    for node_id in node_ids:
        area = _known_node_area(graph.nodes.get(node_id))
        if area is not None:
            return area, f"node '{node_id}'"

    for node_id, elem in neighbours:
        for other in set(elem.all_source_nodes()) | set(elem.all_sink_nodes()):
            if other == node_id:
                continue
            area = _known_node_area(graph.nodes.get(other))
            if area is not None:
                return area, f"node '{other}'"

    return None, ""


def _unresolved_area_error(
    component: "NetworkNode | NetworkElement",
    parameter: str,
    node_ids: "list[str]",
) -> ValueError:
    """The error raised when a load-bearing flow area cannot be inferred.

    Names the component, the parameter that clears the error, and where the
    search looked -- a defaulted area here fabricates geometry, so failing
    loudly is the only honest option (issue #262).
    """
    where = ", ".join(f"'{n}'" for n in node_ids)
    return ValueError(
        f"{type(component).__name__} '{component.id}': cannot determine "
        f"{parameter} -- no neighbouring channel, area-bearing element or "
        f"node with a known area was found at {where}. Set {parameter} "
        f"explicitly on this component (or give a neighbour a known area). "
        f"A default would fabricate a geometry that does not exist and the "
        f"solver would stall on it."
    )


class NetworkNode(ABC):
    #: True when this node supplies a meaningful temperature-rise ratio theta
    #: (T_burned / T_unburned - 1).  Read by ``PressureLossElement`` with
    #: theta-aware correlations.  Overridden on combustion-capable nodes
    #: (``CombustorNode``).  Nodes advertising ``has_theta = True`` must
    #: expose ``_T_unburned`` and ``_T_burned`` float attributes after
    #: ``compute_derived_state`` has been called.
    has_theta: bool = False

    def __init__(self, id: str):
        self.id = id

    @property
    def has_convective_surface(self) -> bool:
        """True if this node has a convective heat transfer surface."""
        return False

    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this node contributes to the solver."""
        pass

    @abstractmethod
    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Returns (residuals, local_jacobian)."""
        pass

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """
        Computes (Tt, Y, Jac) for nodes based on upstream conditions.
        Jac contains derivatives d(Tt)/d(stream_i) and d(Y)/d(stream_i).
        By default, it just takes the first upstream state (pass-through).
        """
        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        up = upstream_states[0]
        return up.Tt, up.Y, None

    @abstractmethod
    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Called automatically by FlowNetwork to resolve neighbors."""
        pass

    def mach(self, state: NetworkMixtureState) -> float:
        """Returns the local Mach number at this node. Default 0.0 (stagnation)."""
        return 0.0

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float]:
        """Compute generalized node-level diagnostics (e.g., thermo properties)."""
        if state.P <= 0:
            return {}

        try:
            cs = cb.complete_state(state.T, state.P, state.X)
            return {
                "Tt": float(state.Tt),
                "Pt": float(state.Pt),
                "h": cs.thermo.h,
                "s": cs.thermo.s,
                "u": cs.thermo.u,
                "rho": cs.thermo.rho,
                "gamma": cs.thermo.gamma,
                "a": cs.thermo.a,
                "cp": cs.thermo.cp,
                "cv": cs.thermo.cv,
                "mw": cs.thermo.mw,
                "mu": cs.transport.mu,
                "k": cs.transport.k,
                "Pr": cs.transport.Pr,
                "nu": cs.transport.nu,
            }
        except Exception:
            ts = cb.thermo_state(state.T, state.P, state.X)
            return {
                "Tt": float(state.Tt),
                "Pt": float(state.Pt),
                "h": ts.h,
                "s": ts.s,
                "u": ts.u,
                "rho": ts.rho,
                "gamma": ts.gamma,
                "a": ts.a,
                "cp": ts.cp,
                "cv": ts.cv,
                "mw": ts.mw,
            }


# ---------------------------------------------------------------------------
# Shared diagnostic helpers
# ---------------------------------------------------------------------------


def _element_pressure_block(
    state_in: "NetworkMixtureState",
    state_out: "NetworkMixtureState",
    mach_in: float,
    mach_out: float,
) -> dict[str, float]:
    """Return the universal 14-field pressure/Mach block for any element."""
    p_ratio_total = state_out.Pt / state_in.Pt if state_in.Pt > 0 else 0.0
    p_ratio_static = state_out.P / state_in.P if state_in.P > 0 else 0.0
    return {
        "P_in": float(state_in.P),
        "P_out": float(state_out.P),
        "T_in": float(state_in.T),
        "T_out": float(state_out.T),
        "Pt_in": float(state_in.Pt),
        "Pt_out": float(state_out.Pt),
        "Tt_in": float(state_in.Tt),
        "Tt_out": float(state_out.Tt),
        "mach_in": float(mach_in),
        "mach_out": float(mach_out),
        "dPt": float(state_in.Pt - state_out.Pt),
        "dP": float(state_in.P - state_out.P),
        "pr_total": float(p_ratio_total),
        "pr_static": float(p_ratio_static),
    }


def _element_reference_block(
    cs: Any,
    *,
    m_dot: float,
    area: float,
    Dh: float | None = None,
    location: str = "inlet",
) -> dict[str, float | str]:
    """Return the thermo/transport reference-state block evaluated at the inlet.

    All keys use the ``_ref`` suffix so it is unambiguous which state they
    refer to.  ``velocity`` and ``Re`` are always included; Re is 0 when
    geometry is insufficient to compute it.
    """
    rho = cs.thermo.rho
    v = (abs(m_dot)) / ((rho) * (area)) if area > 0 and rho > 0 else 0.0
    mu = cs.transport.mu
    re = (v) * (Dh) / (mu) if (Dh and Dh > 0 and v > 0 and mu > 0) else 0.0
    return {
        "ref_location": location,
        "velocity": v,
        "Re": re,
        "rho_ref": float(rho),
        "h_ref": float(cs.thermo.h),
        "s_ref": float(cs.thermo.s),
        "u_ref": float(cs.thermo.u),
        "cp_ref": float(cs.thermo.cp),
        "cv_ref": float(cs.thermo.cv),
        "gamma_ref": float(cs.thermo.gamma),
        "a_ref": float(cs.thermo.a),
        "mw_ref": float(cs.thermo.mw),
        "mu_ref": float(mu),
        "k_ref": float(cs.transport.k),
        "nu_ref": float(cs.transport.nu),
        "Pr_ref": float(cs.transport.Pr),
    }


class LazyJacobian(MutableMapping):
    """An element's local Jacobian, computed on first access (#489).

    The solver's residual-only evaluations never read an element's Jacobian;
    it is assembled only when the root finder asks for J (hybr: at the start
    and on restarts). An element whose derivatives are costly -- the
    compressible channel's implicit Fanno derivatives -- returns this from
    ``residuals`` instead of a dict, and pays for them only when they are
    read. To any reader it is the ``{eq_idx: {unknown: d}}`` dict.
    """

    def __init__(self, compute: Callable[[], dict[int, dict[str, float]]]) -> None:
        self._compute: Callable[[], dict[int, dict[str, float]]] | None = compute
        self._jac: dict[int, dict[str, float]] = {}

    def _resolved(self) -> dict[int, dict[str, float]]:
        if self._compute is not None:
            compute, self._compute = self._compute, None
            self._jac = compute()
        return self._jac

    def __getitem__(self, key: int) -> dict[str, float]:
        return self._resolved()[key]

    def __setitem__(self, key: int, value: dict[str, float]) -> None:
        self._resolved()[key] = value

    def __delitem__(self, key: int) -> None:
        del self._resolved()[key]

    def __iter__(self) -> Iterator[int]:
        return iter(self._resolved())

    def __len__(self) -> int:
        return len(self._resolved())


class NetworkElement(ABC):
    def __init__(self, id: str, from_node: str, to_node: str):
        self.id = id
        self.from_node = from_node
        self.to_node = to_node

    @property
    def has_convective_surface(self) -> bool:
        """True if this element has a convective heat transfer surface."""
        return False

    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this element contributes to the solver."""
        pass

    @abstractmethod
    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """
        Returns (residuals, local_jacobian).
        local_jacobian is a dict mapping equation_index to a dict of {unknown_name: partial_derivative},
        or a LazyJacobian that computes it when first read.
        """
        pass

    @abstractmethod
    def n_equations(self) -> int:
        """Number of residual equations contributed."""
        pass

    @abstractmethod
    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Called automatically by FlowNetwork to resolve neighbors."""
        pass

    def all_source_nodes(self) -> list[str]:
        """All node IDs this element draws flow FROM."""
        return [self.from_node]

    def all_sink_nodes(self) -> list[str]:
        """All node IDs this element delivers flow TO."""
        return [self.to_node]

    def flow_at_node(self, node_id: str, x: Any, indices: list[int]) -> float:
        """Mass flow this element contributes at the given connected node."""
        return float(x[indices[0]])

    def flow_jac_at_node(self, node_id: str, indices: list[int]) -> dict[int, float]:
        """d(mass_flow_at_node)/d(solver_unknown): {global_index: coefficient}."""
        return {indices[0]: 1.0}

    def validate(self) -> None:
        """Perform element-specific validation checks.
        Should raise ValueError with a clear message if validation fails.
        """
        return None

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient and adiabatic wall temperature."""
        return None

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        """Compute diagnostic properties for this element (e.g. Mach, P_ratio)."""
        return {}


class PlenumNode(NetworkNode):
    """
    A stagnation volume where v ~ 0, so Pt = P_static.
    For Phase 1: pressure-only node.
    For Phase 2+: automatically handles mixing when multiple upstream connections exist.
    """

    def __init__(self, id: str):
        super().__init__(id)
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this plenum."""
        self.energy_boundaries.append(eb)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        return [f"{self.id}.P", f"{self.id}.Pt"]

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a plenum (simple mixing + energy boundaries)."""

        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        mix_res = _solver_tools.mixer_from_streams_and_jacobians(
            streams, Q=Q_total, fraction=fraction_total
        )

        return mix_res.T_mix, mix_res.Y_mix, mix_res

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Base residual: Pt = P for plenum
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        return res, jac

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float]:
        # Plenums: v ~ 0, so stagnation = static
        ts = cb.thermo_state(state.T, state.P, state.X)
        return {
            "Tt": float(state.T),
            "Pt": float(state.P),
            "h": ts.h,
            "s": ts.s,
            "u": ts.u,
            "rho": ts.rho,
            "gamma": ts.gamma,
            "a": ts.a,
            "cp": ts.cp,
            "cv": ts.cv,
            "mw": ts.mw,
        }

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Store upstream elements for mixing calculations
        self.upstream_elements = graph.get_upstream_elements(self.id)


class MomentumChamberNode(NetworkNode):
    """
    Similar to a plenum, but preserves dynamic pressure.
    Enforces Pt = P_static + 0.5 * rho * v^2 instead of Pt = P_static.
    Supports directional flow vectors (port angles) for momentum conservation.
    Automatically handles mixing when multiple upstream connections exist.
    """

    def __init__(
        self,
        id: str,
        area: float | None = None,
        port_angles_deg: dict[str, float] | None = None,
        length: float | None = None,
        surface: ConvectiveSurface | None = None,
        t_hot: float | None = None,
        Dh: float | None = None,
        main_inlet: str | None = None,
    ):
        super().__init__(id)
        # MERGE CHAMBER (#471). With main_inlet declared the chamber accepts
        # further inflows as SIDE STREAMS, still one outlet. Its axis is the
        # outlet direction; the main stream arrives along it, each side stream
        # brings m u_jet cos(theta) of axial momentum (the injecting element's
        # injection_momentum) and discharges at the chamber static pressure.
        # The node's (P, Pt) is the chamber/outlet state; the main-inlet
        # element sees the main-face state from chamber_merge_face_state.
        # Splits stay refused: one momentum equation cannot set two outlet
        # pressures.
        self.main_inlet = main_inlet
        self._auto_area = area is None
        self.area = area if area is not None else 0.0
        #: Where self.area came from once resolved (issue #262).
        self._area_source = "user" if area is not None else ""
        self.port_angles_deg: dict[str, float] = port_angles_deg or {}
        self.length = length
        self.Dh = Dh
        self.surface = surface or ConvectiveSurface()
        self.t_hot = t_hot
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []

    @property
    def has_convective_surface(self) -> bool:
        return self.surface.area > 0.0

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this momentum chamber."""
        self.energy_boundaries.append(eb)

    def set_port_angle(self, element_id: str, angle_deg: float) -> None:
        """Set the angle of the flow for a specific connected element."""
        self.port_angles_deg[element_id] = angle_deg

    def get_port_angle(self, element_id: str) -> float:
        """Get the angle of the flow for a specific connected element (defaults to 0.0)."""
        return self.port_angles_deg.get(element_id, 0.0)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        return [f"{self.id}.P", f"{self.id}.Pt"]

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a momentum chamber (simple mixing + energy boundaries)."""

        # Store total mass flow for use in residuals
        self._total_m_dot = sum(s.m_dot for s in upstream_states) if upstream_states else 1.0

        # Build named m_dot Jacobian: {var_name: d(m_total)/d(var)}.
        # flow_jac_at_node(nid) from the solver already gives d(flow_into_nid)/d(each unknown),
        # stored on each upstream state as _m_dot_jac_names. Accumulate across elements,
        # deduplicating by element_id so multi-source elements (tee) are not double-counted.
        self._upstream_element_ids = []
        self._upstream_m_dot_jac: dict[str, float] = {}
        _seen_elem_ids: set[str] = set()
        for s in upstream_states:
            if hasattr(s, "_element_id"):
                self._upstream_element_ids.append(s._element_id)
                eid = s._element_id
                if eid not in _seen_elem_ids and hasattr(s, "_m_dot_jac_names"):
                    _seen_elem_ids.add(eid)
                    for var, coeff in s._m_dot_jac_names.items():
                        self._upstream_m_dot_jac[var] = (
                            self._upstream_m_dot_jac.get(var, 0.0) + coeff
                        )

        # Merge chamber: keep each inflow's supply state for the main-face
        # offset, which needs the chamber's own state (main_face_offset).
        self._inflows = [
            (getattr(s, "_element_id", None), s, dict(getattr(s, "_m_dot_jac_names", {})))
            for s in upstream_states
        ]

        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        mix_res = _solver_tools.mixer_from_streams_and_jacobians(
            streams, Q=Q_total, fraction=fraction_total
        )

        return mix_res.T_mix, mix_res.Y_mix, mix_res

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient for the momentum chamber."""
        if self.surface.area == 0.0:
            return None

        T_hot = self.t_hot if self.t_hot is not None else math.nan
        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho, _ = _safe_rho(state.density())
        u = m_dot_total / (rho * self.area) if self.area > 0 else 1.0

        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)
        length = self.length if self.length is not None else diameter

        return self.surface.htc_and_T(
            T=state.T,
            P=state.P,
            X=state.X,
            velocity=u,
            diameter=diameter,
            length=length,
            T_hot=T_hot,
            flow_area=self.area,
        )

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:

        # Momentum chamber: Pt = P_static + 0.5 * rho * v^2
        # Use total mass flow computed during compute_derived_state
        m_dot_total = getattr(self, "_total_m_dot", 0.0)

        # Use C++ function for residual and analytical Jacobian
        result = _solver_tools.momentum_chamber_residual_and_jacobian(
            state.P, state.Pt, m_dot_total, state.T, state.Y, self.area
        )

        res = [result.residual]
        jac: dict[int, dict[str, float]] = {
            0: {
                f"{self.id}.P": result.d_res_dP,
                f"{self.id}.Pt": result.d_res_dP_total,
                # T is not an unknown at this node; the solver relays the
                # entry through the propagation chain (see the "." branch in
                # _residuals_and_jacobian). The compressible closure depends on
                # T through T0, a(T) and s(T), not just through rho, so the
                # term is no longer small enough to omit.
                f"{self.id}.T": result.d_res_dT,
            }
        }

        # Add Jacobian entries for upstream element mass flows using the
        # named Jacobian built in compute_derived_state. This correctly handles
        # multi-unknown elements (e.g. TeeJunctionElement with m_dot_com /
        # m_dot_branch) by mapping d(residual)/d(m_total) through the chain
        # d(m_total)/d(each_unknown) provided by flow_jac_at_node.
        for var_name, coeff in getattr(self, "_upstream_m_dot_jac", {}).items():
            if coeff != 0.0:
                jac[0][var_name] = jac[0].get(var_name, 0.0) + result.d_res_dmdot * coeff

        return res, jac

    def has_side_streams(self) -> bool:
        """A merge: a declared main inlet and at least one OTHER inflow. The
        main itself need not be flowing in -- a reversed main keeps its face
        (#493)."""
        if self.main_inlet is None:
            return False
        return any(eid != self.main_inlet for eid, _, _ in getattr(self, "_inflows", []))

    def main_face_state(
        self, state: NetworkMixtureState, m_main: float | None = None
    ) -> tuple[float, float, dict[str, float], dict[str, float]]:
        """(P_face, Pt_face, dP_face/d(name), dPt_face/d(name)) for the main inlet.

        Exact (compressible) impulse over the constant-area chamber, C++'s
        ``chamber_merge_face_state``: the face static pressure on the subsonic
        root, the face Pt from the same stagnation closure as the chamber.
        Derivatives are keyed by unknown name; ``"<node>.T"`` keys are relayed
        by the solver.

        ``m_main`` is the main element's signed flow into the chamber, for a
        main that is not among the inflows -- REVERSED (#493). The impulse
        holds for either sign (the x-momentum flux through the face is
        m_main^2 / (rho A) both ways), so the face is kept: the fluid at the
        face is then the chamber's, and the outlet carries the side streams
        less the main's outflow. It used to switch off, a 245 Pa step in the
        main's exit pressure as its flow crossed zero.
        """
        elems = {e.id: e for e in self.upstream_elements}
        main_state, main_names = None, {}
        S, dS = 0.0, {}
        for eid, st, names in getattr(self, "_inflows", []):
            if eid == self.main_inlet:
                main_state, main_names = st, names
                continue
            hook = getattr(elems.get(eid), "injection_momentum", None)
            if hook is None:
                continue
            J, dJ = hook(st, state)
            S += J
            for k, v in dJ.items():
                dS[k] = dS.get(k, 0.0) + v
        out_names = dict(getattr(self, "_upstream_m_dot_jac", {}))
        m_out = float(getattr(self, "_total_m_dot", 0.0))
        main_T_key = None
        if main_state is not None:
            m_main = float(main_state.m_dot)
            T_main = float(main_state.T)
            X_main = main_state.X
        else:
            # Reversed (or idle) main: chamber fluid leaves through the face.
            m_main = float(m_main) if m_main is not None else 0.0
            T_main = float(state.T)
            X_main = state.X
            main_names = {f"{self.main_inlet}.m_dot": 1.0}
            m_out += m_main
            out_names[f"{self.main_inlet}.m_dot"] = (
                out_names.get(f"{self.main_inlet}.m_dot", 0.0) + 1.0
            )
            main_T_key = f"{self.id}.T"
        r = cb.chamber_merge_face_state(
            m_main, T_main, X_main, m_out, state.P, state.T, state.X, S, self.area
        )
        self._face_choked = bool(r.choked)
        main_src = elems[self.main_inlet].from_node if self.main_inlet in elems else None

        def chain(d_mm, d_mo, d_J, d_P, d_T, d_Tm) -> dict[str, float]:
            out: dict[str, float] = {}

            def add(k: str, v: float) -> None:
                out[k] = out.get(k, 0.0) + v

            for k, v in main_names.items():
                add(k, d_mm * v)
            for k, v in out_names.items():
                add(k, d_mo * v)
            for k, v in dS.items():
                add(k, d_J * v)
            add(f"{self.id}.P", d_P)
            add(f"{self.id}.T", d_T)
            if main_T_key is not None:
                add(main_T_key, d_Tm)
            elif main_src is not None:
                add(f"{main_src}.T", d_Tm)
            return out

        return (
            float(r.P_face),
            float(r.Pt_face),
            chain(r.dPf_dm_main, r.dPf_dm_out, r.dPf_dJ, r.dPf_dP, r.dPf_dT, r.dPf_dT_main),
            chain(r.dPtf_dm_main, r.dPtf_dm_out, r.dPtf_dJ, r.dPtf_dP, r.dPtf_dT, r.dPtf_dT_main),
        )

    def mach(self, state: NetworkMixtureState) -> float:
        """Computes Mach number using internal total mass flow and area."""

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        if self.area <= 0 or m_dot_total <= 1e-12:
            return 0.0

        rho = state.density()
        velocity = m_dot_total / (rho * self.area)
        return float(cb.mach_number(velocity, state.T, state.X))

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float]:
        cs = cb.complete_state(state.T, state.P, state.X)

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho = cs.thermo.rho
        u = (abs(m_dot_total)) / ((rho) * (self.area)) if rho > 0 and self.area > 0 else 0.0
        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)
        re = (u) * (diameter) / (cs.transport.mu) if u > 0 and cs.transport.mu > 0 else 0.0
        mach = u / cs.thermo.a if cs.thermo.a > 0 else 0.0
        Tt, Pt = (state.Tt, state.Pt)

        nu_val = 0.0
        htc_val = 0.0
        t_aw_val = float(state.T)
        if self.has_convective_surface:
            h_res = self.htc_and_T(state)
            if h_res is not None:
                nu_val = h_res.Nu
                htc_val = h_res.h
                t_aw_val = h_res.T_aw

        return {
            "Tt": float(Tt),
            "Pt": float(Pt),
            "h": cs.thermo.h,
            "s": cs.thermo.s,
            "u": cs.thermo.u,
            "rho": rho,
            "gamma": cs.thermo.gamma,
            "a": cs.thermo.a,
            "cp": cs.thermo.cp,
            "cv": cs.thermo.cv,
            "mw": cs.thermo.mw,
            "mu": cs.transport.mu,
            "k": cs.transport.k,
            "Pr": cs.transport.Pr,
            "nu": cs.transport.nu,
            "Re": re,
            "Dh": diameter,
            "velocity": u,
            "mach": mach,
            "Nu": nu_val,
            "htc": htc_val,
            "T_aw": t_aw_val,
        }

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        self.upstream_elements = graph.get_upstream_elements(self.id)
        if self.main_inlet is not None:
            # Merge chamber: one area, and the main inlet's port area IS it.
            # Inherit it from the main-inlet element; an explicit area that
            # disagrees is refused in validate(), never silently reconciled.
            main = graph.elements.get(self.main_inlet)
            main_area = getattr(main, "area", None) if main is not None else None
            self._main_inlet_area = main_area if main_area else None
            if self._auto_area and self._main_inlet_area:
                self.area = self._main_inlet_area
                self.Dh = getattr(main, "Dh", None) or 2.0 * math.sqrt(self.area / math.pi)
                self.surface.area = self.area
                self._area_source = f"inherited from {self.main_inlet}"
                return
        # Inherit hydraulic diameter from upstream channel when not user-specified
        if self.Dh is None:
            for elem in self.upstream_elements:
                if isinstance(elem, ChannelElement) and elem.diameter is not None:
                    self.Dh = elem.diameter
                    break
        # Derive cross-section area from Dh when area was not user-specified
        if self._auto_area and self.Dh is not None:
            self.area = math.pi * (self.Dh / 2.0) ** 2
            self.surface.area = self.area
        elif self._auto_area:
            # No channel to inherit from: widen the search to any area-bearing
            # neighbour before falling back (issue #262).
            area, source = _infer_flow_area(graph, [self.id])
            if area is not None:
                self.area = area
                self.Dh = 2.0 * math.sqrt(area / math.pi)
                self.surface.area = area
                self._area_source = source
            else:
                # Unlike AreaChangeElement, this area is a dynamic-head
                # reference rather than the geometry being modelled, and a
                # network can legitimately never read it -- so a nominal
                # default stands. It is recorded, not silent: anything that
                # consumes it reports the provenance in diagnostics.
                self.area = DEFAULT_CHAMBER_AREA
                self.surface.area = self.area
                self._area_source = "default"

    def validate(self) -> None:
        if self.main_inlet is None:
            return
        ids = [e.id for e in self.upstream_elements]
        if self.main_inlet not in ids:
            raise ValueError(
                f"MomentumChamberNode {self.id!r}: main_inlet {self.main_inlet!r} "
                f"does not flow into it (inflows: {sorted(ids)})."
            )
        main_area = getattr(self, "_main_inlet_area", None)
        if not self._auto_area and main_area and abs(self.area / main_area - 1.0) > 1e-9:
            raise ValueError(
                f"MomentumChamberNode {self.id!r}: area {self.area:.6g} m^2 differs "
                f"from its main inlet {self.main_inlet!r} ({main_area:.6g} m^2). "
                "A merge chamber is a constant-area control volume; leave the "
                "area unset to inherit it."
            )


class PressureBoundary(NetworkNode):
    """
    Supplies fixed stagnation pressure and temperature. Contributes no unknowns and no residuals.
    Essential as a reference pressure and absolute mass sink/source for the network.
    """

    def __init__(
        self,
        id: str,
        Pt: float = 101325.0,
        Tt: float = 300.0,
        Y: list[float] | None = None,
        coupling: str = "auto",
    ) -> None:
        super().__init__(id)
        self.Pt = Pt
        self.Tt = Tt
        self.Y = Y
        #: How the supplied pressure couples to a duct meeting this boundary.
        #:
        #: ``"total"``   -- it is the stagnation pressure at the connection.
        #:                  Correct for an INFLOW: a reservoir supplying the
        #:                  network, where Pt is what the reservoir holds.
        #: ``"static"``  -- it is the static pressure at the duct face, and the
        #:                  duct's exit dynamic head is dissipated in the
        #:                  expansion. The standard model for a bare duct
        #:                  discharging into a plenum or to atmosphere.
        #: ``"auto"``    -- infer from flow direction: inflow takes total,
        #:                  outflow takes static.
        #:
        #: ``auto`` is well defined wherever a boundary has one role. It cannot
        #: decide for a boundary whose flow REVERSES during a solve, or that
        #: serves both roles -- inflow wants total, outflow wants static, and
        #: there is no single right answer. That case is what the explicit
        #: setting is for; whoever builds such a network knows whether the exit
        #: is a plain opening or a diffuser, and the solver does not.
        #:
        #: Pinning stagnation pressure at an outflow caps the mass flux at the
        #: sonic value FOR THAT PRESSURE, a constraint the physical problem
        #: never imposed. Measured on the combustor of #351: 1.34x an
        #: impossible ceiling under total coupling, M = 0.752 and comfortable
        #: under static. See issue #360.
        self.coupling: str = coupling

    def exit_head_lost(self, is_outflow: bool) -> bool:
        """Whether a duct meeting this boundary loses its exit dynamic head."""
        if self.coupling == "static":
            return True
        if self.coupling == "total":
            return False
        return bool(is_outflow)

    def unknowns(self) -> list[str]:
        return []

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        return [], {}

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        # Boundary nodes define their own state

        Y = self.Y if self.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        return self.Tt, Y, None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class MassFlowBoundary(NetworkNode):
    """
    Supplies fixed mass flow and stagnation temperature. Its pressure floats to satisfy the flow equations.
    Cannot exclusively define a network, a pressure reference is also required.
    """

    def __init__(
        self,
        id: str,
        m_dot: float = 0.1,
        Tt: float = 300.0,
        Y: list[float] | None = None,
    ) -> None:
        super().__init__(id)
        self.m_dot = m_dot
        self.Tt = Tt
        self.Y = Y

    def unknowns(self) -> list[str]:
        # Pressure floats to whatever is required to push the defined mass flow
        return [f"{self.id}.P", f"{self.id}.Pt"]

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Stagnation assumptions: Pt = P_static
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        return res, jac

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """
        Computes (Tt, Y, Jac) for a MassFlowBoundary.
        - If upstream_states exists, it acts as a SINK (Outlet): inherits state via mixing.
        - If upstream_states is empty, it acts as a SOURCE (Inlet): uses user settings.
        """

        if upstream_states:
            # SINK behavior: Automatically mix upstream streams using C++ core logic
            # This ensures energy conservation and correct species transport at outlets.
            streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
            mix_res = _solver_tools.mixer_from_streams_and_jacobians(streams)
            return mix_res.T_mix, mix_res.Y_mix, mix_res

        # SOURCE behavior: Use developer/user-defined boundary constants
        Y = self.Y if self.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        return self.Tt, Y, None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class WallNode(NetworkNode):
    """Closed-end boundary that imposes zero mass flow.

    The Pt = P residual enforces zero velocity at the dead end.
    The solver's automatic mass-balance equation (applied to every non-PressureBoundary
    node) forces m_dot = 0 on all connected elements without any additional logic here.
    """

    def __init__(self, id: str) -> None:
        super().__init__(id)

    def unknowns(self) -> list[str]:
        return [f"{self.id}.P", f"{self.id}.Pt"]

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        return res, jac

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


def _zero_loss(ctx: object) -> tuple[float, float]:
    return 0.0, 0.0


class CombustorNode(NetworkNode):
    """
    Combustion chamber adding energy (and optionally mass) to the flow.
    Evaluates adiabatic combustion using the C++ core backend based on
    the selected CombustionMethodLiteral.
    Automatically handles mixing when multiple upstream connections exist.

    Pressure loss is **not** a property of the combustor node.  Attach a
    :class:`PressureLossElement` on an adjacent edge (upstream diffuser or
    downstream liner exit) to apply combustor pressure loss; theta-aware
    correlations (``LinearThetaFractionLoss`` / ``LinearThetaHeadLoss``)
    automatically pick up ``_T_unburned`` / ``_T_burned`` from this node.
    """

    # Supplies theta = T_burned / T_unburned - 1 to adjacent PressureLossElement.
    has_theta: bool = True

    def __init__(
        self,
        id: str,
        method: CombustionMethodLiteral = "complete",
        area: float | None = None,
        Dh: float | None = None,
        surface: ConvectiveSurface | None = None,
        t_hot: float | None = None,
    ):
        super().__init__(id)
        self.method = method
        self._auto_area = area is None
        self.area = area if area is not None else 0.0
        #: Where self.area came from once resolved (issue #262).
        self._area_source = "user" if area is not None else ""
        self.Dh = Dh
        self.surface = surface or ConvectiveSurface()
        self.t_hot = t_hot
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []
        # Populated by compute_derived_state so PressureLossElement can read theta.
        self._T_unburned: float = 300.0
        self._T_burned: float = 300.0

    @property
    def has_convective_surface(self) -> bool:
        return self.surface.area > 0.0

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this combustor (post-combustion)."""
        self.energy_boundaries.append(eb)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        return [f"{self.id}.P", f"{self.id}.Pt"]

    def compute_derived_state(
        self, upstream_states: list[NetworkMixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a combustor (Reaction + Mixing + energy boundaries)."""

        # Store total mass flow and unburned temperature for use in diagnostics
        self._total_m_dot = sum(s.m_dot for s in upstream_states) if upstream_states else 1.0
        self._T_unburned = (
            sum(s.m_dot * s.Tt for s in upstream_states) / self._total_m_dot
            if self._total_m_dot > 0
            else 300.0
        )
        # The streams behind T_unburned, for its Jacobian in a theta-sourced
        # PressureLossElement: (m, Tt, {unknown: dm/d(unknown)}, source node).
        # T_u = sum(m Tt) / M is not a relayed node property (#481 C1).
        self._unburned_streams = [
            (
                float(s.m_dot),
                float(s.Tt),
                dict(getattr(s, "_m_dot_jac_names", {})),
                getattr(s, "_src_node", None),
            )
            for s in upstream_states
        ]

        # Store upstream element IDs for momentum-chamber Jacobian (d_res_dmdot entries)
        self._upstream_element_ids = [
            s._element_id for s in upstream_states if hasattr(s, "_element_id")
        ]

        if not upstream_states:
            # Default fallback
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.Tt, s.Pt, s.Y) for s in upstream_states]
        if self._total_m_dot > 0.0:
            P_ref = sum(s.m_dot * s.Pt for s in upstream_states) / self._total_m_dot
        else:
            P_ref = upstream_states[0].Pt

        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        # Combustor only: mixing + combustion, NO pressure loss (that's on the edge).
        try:
            mix_res = _solver_tools.combustor_residuals_and_jacobians(
                streams,
                P_ref,
                Q=Q_total,
                fraction=fraction_total,
                pressure_loss=_zero_loss,
                use_equilibrium=(self.method == "equilibrium"),
            )
        except Exception as exc:
            # Unphysical intermediate state during Newton iteration (e.g. extreme phi).
            # Re-raise so the solver penalty path can guide the step back to physics.
            raise RuntimeError(
                f"CombustorNode '{self.id}': combustion call failed "
                f"(T_unburned={self._T_unburned:.1f} K)"
            ) from exc
        self._last_mix_res = mix_res
        # Expose burned temperature for adjacent PressureLossElement (theta source).
        self._T_burned = float(mix_res.T_mix)
        return mix_res.T_mix, mix_res.Y_mix, mix_res

    def residuals(
        self, state: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Stagnation constraint: Pt = P (combustor is a large-area, low-velocity volume).
        # The dynamic pressure 0.5*rho*v^2 is negligible compared to P at combustor scales,
        # so P ~= Pt is an excellent approximation (same as PlenumNode).
        # Pressure loss is applied by an adjacent PressureLossElement, not here.
        res = [state.Pt - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.Pt": 1.0}}
        return res, jac

    def diagnostics(self, state: NetworkMixtureState) -> dict[str, float]:
        if state.P <= 0 or state.T <= 0:
            return {}
        cs = cb.complete_state(state.T, state.P, state.X)

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho = cs.thermo.rho
        u = (abs(m_dot_total)) / ((rho) * (self.area)) if rho > 0 and self.area > 0 else 0.0
        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)
        re = (u) * (diameter) / (cs.transport.mu) if u > 0 and cs.transport.mu > 0 else 0.0
        mach = u / cs.thermo.a if cs.thermo.a > 0 else 0.0
        Tt, Pt = (state.Tt, state.Pt)

        nu_val = 0.0
        htc_val = 0.0
        t_aw_val = float(state.T)
        if self.has_convective_surface:
            h_res = self.htc_and_T(state)
            if h_res is not None:
                nu_val = h_res.Nu
                htc_val = h_res.h
                t_aw_val = h_res.T_aw

        try:
            phi = float(cb.equivalence_ratio(state.X))
        except Exception as e:
            print(f"DEBUG: CombustorNode.diagnostics phi calculation failed: {e}", file=sys.stderr)
            phi = 0.0

        T_u = getattr(self, "_T_unburned", state.Tt)
        theta = (state.Tt / T_u) - 1.0 if T_u > 0 else 0.0

        return {
            "Tt": float(Tt),
            "Pt": float(Pt),
            "m_dot": m_dot_total,
            "h": cs.thermo.h,
            "s": cs.thermo.s,
            "u": cs.thermo.u,
            "rho": rho,
            "gamma": cs.thermo.gamma,
            "a": cs.thermo.a,
            "cp": cs.thermo.cp,
            "cv": cs.thermo.cv,
            "mw": cs.thermo.mw,
            "mu": cs.transport.mu,
            "k": cs.transport.k,
            "Pr": cs.transport.Pr,
            "nu": cs.transport.nu,
            "Re": re,
            "Dh": diameter,
            "velocity": u,
            "mach": mach,
            "Nu": nu_val,
            "htc": htc_val,
            "T_aw": t_aw_val,
            "phi": phi,
            "theta": theta,
        }

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient for the combustor wall."""
        if self.surface.area == 0.0:
            return None

        T_hot = self.t_hot if self.t_hot is not None else math.nan
        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho, _ = _safe_rho(state.density())
        u = m_dot_total / (rho * self.area) if self.area > 0 else 1.0

        diameter = self.Dh if self.Dh is not None else 0.01
        length = 0.1  # Combustors are modelled as a single zone

        return self.surface.htc_and_T(
            T=state.T,
            P=state.P,
            X=state.X,
            velocity=u,
            diameter=diameter,
            length=length,
            T_hot=T_hot,
            flow_area=self.area,
        )

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        self.upstream_elements = graph.get_upstream_elements(self.id)
        # Inherit hydraulic diameter from upstream channel when not user-specified
        if self.Dh is None:
            for elem in self.upstream_elements:
                if isinstance(elem, ChannelElement) and elem.diameter is not None:
                    self.Dh = elem.diameter
                    break
        # Derive cross-section area from Dh when area was not user-specified
        if self._auto_area and self.Dh is not None:
            self.area = math.pi * (self.Dh / 2.0) ** 2
            self.surface.area = self.area
        elif self._auto_area:
            # No channel to inherit from: widen the search to any area-bearing
            # neighbour before falling back (issue #262).
            area, source = _infer_flow_area(graph, [self.id])
            if area is not None:
                self.area = area
                self.Dh = 2.0 * math.sqrt(area / math.pi)
                self.surface.area = area
                self._area_source = source
            else:
                # Unlike AreaChangeElement, this area is a dynamic-head
                # reference rather than the geometry being modelled, and a
                # network can legitimately never read it -- so a nominal
                # default stands. It is recorded, not silent: anything that
                # consumes it reports the provenance in diagnostics.
                self.area = DEFAULT_CHAMBER_AREA
                self.surface.area = self.area
                self._area_source = "default"


# The discharge-hole family: a hole in a wall, no pipe, no beta. 'fixed' and
# the ISO 5167 metering correlations are deliberately not in this map.
_DISCHARGE_SELECTORS = {
    "IdelchikThick": cb.DischargeCdCorrelation.Idelchik1966Thick,
    "IdelchikBeveled": cb.DischargeCdCorrelation.Idelchik1966Beveled,
    "IdelchikRounded": cb.DischargeCdCorrelation.Idelchik1966Rounded,
    "McGreehanSchotsch": cb.DischargeCdCorrelation.McGreehanSchotsch1988,
    "Lichtarowicz": cb.DischargeCdCorrelation.Lichtarowicz1965,
}


class OrificeElement(NetworkElement):
    """
    Orifice flow element with incompressible or compressible formulation.

    - regime='incompressible': m_dot = Cd * A * sqrt(2 * rho * dP)
    - regime='compressible': Uses isentropic nozzle flow with smooth choked transition

    The discharge coefficient Cd is determined by the 'correlation' parameter:
      - 'fixed': Uses the user-supplied Cd value directly.
      - 'ReaderHarrisGallagher': Sharp thin-plate (ISO 5167-2 / RHG).
      - 'Stolz': ISO 5167:1980 (Corner taps).
      - 'Miller': Miller (1996) simplified correlation.
      - 'IdelchikThick': Deep hole in a wall, Idelchik diagram 4-18a
        (requires plate_thickness). Valid Re 25 to 1e6.
      - 'IdelchikBeveled': Beveled-edge hole, diagram 4-18b
        (requires bevel_depth).
      - 'IdelchikRounded': Rounded-edge hole, diagram 4-18c
        (requires edge_radius).
      - 'McGreehanSchotsch': Cooling hole with inlet crossflow (1988).
      - 'Lichtarowicz': LONG orifice, Lichtarowicz, Duggins and Markland
        (1965), l/d 2-10 and Re 10 to 2e4 (requires plate_thickness).
        Refused below l/d = 1.5, where the source reports hysteresis.

    The first four are NORMED metering correlations: Cd is referenced to the
    tapping differential and is a function of beta = d/D. The last five are
    DISCHARGE correlations for a hole in a wall, where zeta is referenced to
    the hole velocity and there is no beta. They are not interchangeable, and
    there is deliberately no 'Auto' arm choosing between them from geometry.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        Cd: float = 0.6,
        diameter: float | None = None,
        regime: Literal["incompressible", "compressible"] = "incompressible",
        correlation: str = "ReaderHarrisGallagher",
        plate_thickness: float = 0.0,
        edge_radius: float = 0.0,
        bevel_depth: float = 0.0,
        area: float | None = None,
        injection_angle_deg: float = 90.0,
    ):
        super().__init__(id, from_node, to_node)
        # Angle of the discharged jet to the axis of a merge chamber it feeds
        # as a side stream (#471). 90 = normal injection, no axial momentum.
        self.injection_angle_deg = float(injection_angle_deg)

        if diameter is not None:
            self.diameter: float | None = diameter
            self.area: float | None = math.pi * (diameter / 2.0) ** 2
        elif area is not None:
            self.area = area
            self.diameter = math.sqrt(4.0 * area / math.pi)
        else:
            # Both None: inherit from upstream channel in resolve_topology
            self.diameter = None
            self.area = None

        self.Cd = Cd
        self.regime = regime
        self.correlation = correlation
        self.use_correlation = correlation != "fixed"
        self.plate_thickness = plate_thickness
        self.edge_radius = edge_radius
        # Bevel DEPTH along the hole axis, Idelchik diagram 4-18b's l/Dh at a
        # bevel angle of 40-60 deg. Not the wall thickness.
        self.bevel_depth = bevel_depth
        self.upstream_diameter: float | None = None
        self.downstream_diameter: float | None = None
        # OrificeGeometry built in resolve_topology
        self._orifice_geom: object | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Evaluate both sides of the node for geometry discovery
        upstream_elements = graph.get_upstream_elements(self.from_node)
        downstream_elements = graph.get_downstream_elements(self.to_node)

        # Ensure we don't pick ourselves if it's a tight chain
        upstream_channels = [
            e for e in upstream_elements if isinstance(e, ChannelElement) and e.id != self.id
        ]
        downstream_channels = [
            e for e in downstream_elements if isinstance(e, ChannelElement) and e.id != self.id
        ]

        if len(upstream_channels) == 1:
            self.upstream_diameter = upstream_channels[0].diameter

        if len(downstream_channels) == 1:
            self.downstream_diameter = downstream_channels[0].diameter

        # Bore must be explicitly set by the user; fall back to 0.08 when unspecified.
        # (upstream_diameter is kept only for the velocity-of-approach Beta correction.)
        if self.diameter is None:
            self.diameter = 0.08
            self.area = math.pi * (self.diameter / 2.0) ** 2

        # Pre-compute beta for velocity-of-approach factor

        d_bore = math.sqrt(4.0 * self.area / math.pi)
        self.beta = 0.0
        self._orifice_geom = None

        if self.upstream_diameter and self.upstream_diameter > 0:
            if d_bore < self.upstream_diameter:
                self.beta = d_bore / self.upstream_diameter
            else:
                warnings.warn(
                    f"OrificeElement '{self.id}': inferred bore diameter "
                    f"({d_bore:.4f} m) >= channel diameter "
                    f"({self.upstream_diameter:.4f} m). "
                    f"Skipping velocity-of-approach correction (E=1).",
                    stacklevel=2,
                )

        if self.use_correlation:
            # Build geometry descriptor for Cd correlation.
            # D=0 when no upstream channel known (plenum connection) -> beta=0 -> RHG extrapolates to ~0.597.
            D_up = (
                self.upstream_diameter
                if self.upstream_diameter and self.upstream_diameter > 0
                else 1.0
            )
            if D_up <= d_bore:
                D_up = d_bore * 10.0
            geom = cb.OrificeGeometry()
            geom.d = d_bore
            geom.D = D_up
            geom.t = self.plate_thickness
            geom.r = self.edge_radius
            self._orifice_geom = geom

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def _hole_reynolds(self, state_in: "NetworkMixtureState") -> float:
        """Reynolds number of the flow through ONE hole, based on its diameter.

        The discharge-hole correlations (Idelchik, McGreehan-Schotsch) are
        functions of this, not of the pipe Reynolds number. For a plain
        orifice there is one hole carrying all the flow; EffusionPlateElement
        overrides this to divide by the hole count.
        """
        d = self._orifice_geom.d if self._orifice_geom is not None else 0.0
        if d <= 0.0:
            return 1.0e5
        mu = cb.transport_state(state_in.T, state_in.P, state_in.X).mu
        if mu <= 0.0:
            mu = 1.8e-5
        m_dot_hole = abs(state_in.m_dot) / self._hole_count()
        if m_dot_hole <= 1e-12:
            # No flow yet: a mid-range seed, not the top of the table, so the
            # first Newton step starts on a responsive part of the curve.
            return 1.0e4
        return (4.0 * m_dot_hole) / (math.pi * d * mu)

    def _hole_count(self) -> float:
        """Holes carrying the element's mass flow in parallel. One, here."""
        return 1.0

    def _effective_Cd(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> float:
        """Return effective Cd: from correlation or user-supplied fixed value."""
        if not self.use_correlation or self._orifice_geom is None:
            # Clamp manual Cd strictly <= 1.0 to strictly preserve vena contracta physics (A_eff <= A_geom)
            return max(1e-4, min(self.Cd, 1.0))

        # Estimate Re_D from current m_dot and upstream viscosity.
        D_up = self._orifice_geom.D
        if D_up > 0.0:
            ts = cb.transport_state(state_in.T, state_in.P, state_in.X)
            mu = ts.mu if ts.mu > 0.0 else 1.8e-5
            m_dot_abs = abs(state_in.m_dot)
            Re_D = (4.0 * m_dot_abs) / (math.pi * D_up * mu) if m_dot_abs > 1e-12 else 1e4
        else:
            Re_D = 1e5  # typical reference when no upstream channel

        flow_state = cb.OrificeState()
        flow_state.Re_D = Re_D
        flow_state.dP = max(state_in.Pt - state_out.P, 1.0)
        flow_state.rho = cb.density(state_in.T, state_in.P, state_in.X)
        flow_state.mu = cb.transport_state(state_in.T, state_in.P, state_in.X).mu

        # The ISO 5167 correlations below are functions of the PIPE Reynolds
        # number Re_D, but the discharge-hole family (Idelchik,
        # McGreehan-Schotsch) is a function of the HOLE Reynolds number.
        # Feeding Re_D to those was wrong by 1-8% in a pipe, and by up to 38%
        # for a plenum-fed hole where D_up = 0 freezes Re_D at the 1e5
        # fallback and the correlation stops responding to flow entirely.
        Re_hole = self._hole_reynolds(state_in)

        # Correlation selection logic
        if self.correlation == "ReaderHarrisGallagher":
            return float(cb.Cd_sharp_thin_plate(self._orifice_geom, flow_state))
        elif self.correlation == "Stolz":
            b = self._orifice_geom.beta
            b2 = b * b
            b8 = b2 * b2 * b2 * b2
            cd = (
                0.5961
                + 0.0261 * b2
                - 0.216 * b8
                + 0.000521 * math.pow(1e6 * b / max(flow_state.Re_D, 1.0), 0.7)
            )
            return float(cd)
        elif self.correlation == "Miller":
            b = self._orifice_geom.beta
            b2 = b * b
            b8 = b2 * b2 * b2 * b2
            cd = (
                0.5959
                + 0.0312 * math.pow(b, 2.1)
                - 0.184 * b8
                + 91.71 * math.pow(b, 2.5) * math.pow(max(flow_state.Re_D, 1.0), -0.75)
            )
            return float(cd)
        elif self.correlation in _DISCHARGE_SELECTORS:
            # The discharge-hole family: a hole in a wall, not a plate in a
            # pipe, so DischargeHoleGeometry and no pipe diameter -- the
            # reason the old ThickPlate/RoundedEntry arms were wrong: they
            # multiplied an ISO 5167 metering Cd by a correction factor. The
            # supply-side crossflow U1/Vi (0 for a plenum-fed hole) only moves
            # McGreehan-Schotsch; the others have no crossflow term.
            u1_vi, _ = self._supply_crossflow(state_in, state_out, flows)
            return float(
                cb.discharge_cd(
                    _DISCHARGE_SELECTORS[self.correlation],
                    self._discharge_hole(),
                    cb.DischargeHoleState(Re=Re_hole, U1_over_Vi=u1_vi),
                )
            )
        else:
            raise ValueError(
                f"OrificeElement: unknown correlation {self.correlation!r}. "
                "The 'Auto' arm was removed: it picked a correlation from the "
                "geometry behind the caller's back, which is how a "
                "rounded-entry request came back as Stolz. Name one of "
                "'fixed', 'ReaderHarrisGallagher', 'Stolz', 'Miller', "
                "'IdelchikThick', 'IdelchikBeveled', 'IdelchikRounded', "
                "'McGreehanSchotsch', 'Lichtarowicz'."
            )

    def injection_momentum(
        self, state_in: "NetworkMixtureState", state_chamber: "NetworkMixtureState"
    ) -> tuple[float, dict[str, float]]:
        """Streamwise momentum this jet brings into a merge chamber (#471).

        J = m w cos(theta), with w from C++'s ``jet_impulse``: the isentropic
        velocity from the supply's (Pt, Tt) to the chamber static pressure,
        or, choked, the sonic momentum plus its pressure thrust. Cd sets the
        jet's area, not its velocity, so it does not enter. Returns (J,
        dJ/d(name)) keyed by unknown name. Normal injection (90 deg) brings
        none.
        """
        cos_t = math.cos(math.radians(self.injection_angle_deg))
        if abs(cos_t) < 1e-12:
            return 0.0, {}
        jet = cb.jet_impulse(
            float(state_in.m_dot), state_in.Pt, state_in.Tt, state_chamber.P, state_in.X
        )
        return jet.J * cos_t, {
            f"{self.id}.m_dot": jet.dJ_dm * cos_t,
            f"{self.from_node}.Pt": jet.dJ_dPt * cos_t,
            f"{self.from_node}.T": jet.dJ_dTt * cos_t,
            f"{self.to_node}.P": jet.dJ_dP * cos_t,
        }

    def _discharge_hole(self) -> "cb.DischargeHoleGeometry":
        """One real hole, as the discharge-hole correlations see it."""
        hole = cb.DischargeHoleGeometry(
            d=self._orifice_geom.d,
            L=self.plate_thickness,
            r=self.edge_radius,
        )
        hole.bevel = self.bevel_depth
        return hole

    def validate(self) -> None:
        # Lichtarowicz refuses short holes (l/d < 1.5, the source's own
        # hysteresis warning). Ask it once here, so the refusal names this
        # element at set-up rather than surfacing mid-solve. The limit lives
        # in C++ only.
        if self.correlation == "Lichtarowicz" and self._orifice_geom is not None:
            try:
                cb.discharge_cd(
                    _DISCHARGE_SELECTORS[self.correlation],
                    self._discharge_hole(),
                    cb.DischargeHoleState(Re=1.0e4),
                )
            except ValueError as exc:
                raise ValueError(
                    f"{type(self).__name__} {self.id!r}: plate_thickness/d = "
                    f"{self.plate_thickness / self._orifice_geom.d:.3g}. {exc}"
                ) from exc

    def _supply_crossflow(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None,
    ) -> tuple[float, dict[str, float]]:
        """McGreehan-Schotsch's supply-side U1/Vi and d(U1/Vi)/d(name).

        0 for a plain orifice: a plenum feeds it. A duct-fed hole overrides
        this (EffusionPlateElement with crossflow segments).
        """
        return 0.0, {}

    def _dCd_dnames(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> dict[str, float]:
        """d(Cd)/d(unknown), keyed by name, analytic.

        Only the discharge-hole family carries it (C++'s
        discharge_cd_and_derivatives): through the hole Reynolds number
        (linear in m_dot, so dRe/dm = Re/m) and, for a duct-fed hole, through
        U1/Vi. 'fixed' has none, and the normed metering correlations keep
        their documented gap.
        """
        selector = _DISCHARGE_SELECTORS.get(self.correlation)
        if selector is None or self._orifice_geom is None or not self.use_correlation:
            return {}
        u1_vi, du = self._supply_crossflow(state_in, state_out, flows)
        Re = self._hole_reynolds(state_in)
        _, dCd_dRe, dCd_dU = cb.discharge_cd_and_derivatives(
            selector, self._discharge_hole(), cb.DischargeHoleState(Re=Re, U1_over_Vi=u1_vi)
        )
        out: dict[str, float] = {}
        m = float(state_in.m_dot)
        if abs(m) > 1e-12:
            out[f"{self.id}.m_dot"] = float(dCd_dRe) * Re / abs(m) * (1.0 if m > 0 else -1.0)
        for name, d in du.items():
            out[name] = out.get(name, 0.0) + float(dCd_dU) * d
        return out

    def _cd_in_range(self, state_in: "NetworkMixtureState") -> bool:
        """Is the discharge-hole Cd inside its source's Re and l/d range?

        A flag for diagnostics, never a selector. True for 'fixed' (the value
        is the caller's) and for the normed metering family, whose ranges are
        ISO 5167's and are not checked here.
        """
        selector = _DISCHARGE_SELECTORS.get(self.correlation)
        if selector is None or self._orifice_geom is None:
            return True
        return bool(
            cb.discharge_cd_in_range(
                selector,
                self._discharge_hole(),
                cb.DischargeHoleState(Re=self._hole_reynolds(state_in)),
            )
        )

    def residuals(
        self,
        state_in: "NetworkMixtureState",
        state_out: "NetworkMixtureState",
        flows: dict[str, float] | None = None,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:

        m_dot = state_in.m_dot
        effective_cd = self._effective_Cd(state_in, state_out, flows)

        if self.regime == "compressible":
            res_cpp = _solver_tools.orifice_compressible_residuals_and_jacobian(
                m_dot,
                state_in.Pt,
                state_in.Tt,
                state_in.Y,
                state_out.P,
                effective_cd,
                self.area,
                getattr(self, "beta", 0.0),
            )
        else:
            # Use incompressible Bernoulli formulation. The density
            # reference is the upstream static by default; the solver's
            # compressible-seed proxy switches it to the downstream static
            # (_incompressible_p_ref = "outlet"), which tracks the
            # compressible solution far better on blow-down networks
            # (6-7x closer junction pressures on the 2026-07-05 outlet-ref
            # seed experiment, tmp/outlet_ref_seed_experiment.py).
            _p_ref_outlet = getattr(self, "_incompressible_p_ref", "inlet") == "outlet"
            res_cpp = _solver_tools.orifice_residuals_and_jacobian(
                m_dot,
                state_in.Pt,
                state_out.P if _p_ref_outlet else state_in.P,
                state_in.T,
                state_in.Y,
                state_out.P,
                effective_cd,
                self.area,
                beta=getattr(self, "beta", 0.0),
            )

        res = [m_dot - res_cpp.m_dot_calc]

        # Assemble Jacobian with respect to all node and element unknowns
        jac = {
            0: {
                f"{self.id}.m_dot": 1.0,
                f"{self.from_node}.Pt": -res_cpp.d_mdot_dP_total_up,
                f"{self.from_node}.T": -res_cpp.d_mdot_dT_up,
            }
        }
        if self.regime != "compressible" and (
            getattr(self, "_incompressible_p_ref", "inlet") == "outlet"
        ):
            # Density evaluated at the downstream static: its sensitivity
            # moves onto the downstream node's P.
            jac[0][f"{self.to_node}.P"] = -(
                res_cpp.d_mdot_dP_static_down + res_cpp.d_mdot_dP_static_up
            )
        else:
            jac[0][f"{self.from_node}.P"] = -res_cpp.d_mdot_dP_static_up
            jac[0][f"{self.to_node}.P"] = -res_cpp.d_mdot_dP_static_down
        # Add species sensitivities
        for i, val in enumerate(res_cpp.d_mdot_dY_up):
            jac[0][f"{self.from_node}.Y[{i}]"] = -val

        # m_calc is proportional to Cd, and a discharge-hole Cd moves with the
        # flow (Re_hole) and, duct-fed, with the crossflow (U1/Vi):
        # d(m_calc)/d(x) += (m_calc/Cd) dCd/dx.
        if effective_cd > 0.0:
            scale = res_cpp.m_dot_calc / effective_cd
            for name, dcd in self._dCd_dnames(state_in, state_out, flows).items():
                jac[0][name] = jac[0].get(name, 0.0) - scale * dcd

        return res, jac

    def n_equations(self) -> int:
        return 1

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float]:
        if state_in.P <= 0 or state_in.T <= 0:
            return {}
        cs = cb.complete_state(state_in.T, state_in.P, state_in.X)

        # Throat Mach: use isentropic pressure ratio from total inlet to static outlet
        gamma = cs.thermo.gamma
        pr_crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
        p_ratio_throat = state_out.P / state_in.Pt if state_in.Pt > 0 else 1.0
        mach_throat = (
            1.0
            if p_ratio_throat <= pr_crit
            else float(
                cb.mach_from_pressure_ratio(state_in.Tt, state_in.Pt, state_out.P, state_in.X)
            )
        )

        # Physical velocity / Mach at the geometric bore (not Cd-adjusted)
        Dh = math.sqrt(4.0 * self.area / math.pi) if self.area > 0 else 0.0
        ref = _element_reference_block(
            cs, m_dot=state_in.m_dot, area=self.area, Dh=Dh, location="inlet"
        )
        v_in = ref["velocity"]
        mach_in = v_in / cs.thermo.a if cs.thermo.a > 0 else 0.0

        # Outlet Mach (use outlet state for accuracy)
        try:
            cs_out = cb.complete_state(state_out.T, state_out.P, state_out.X)
            v_out = (
                (abs(state_in.m_dot)) / ((cs_out.thermo.rho) * (self.area))
                if self.area > 0 and cs_out.thermo.rho > 0
                else 0.0
            )
            mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 else 0.0
        except Exception:
            mach_out = 0.0

        return {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            "mach_throat": float(mach_throat),
            "Cd": float(self._effective_Cd(state_in, state_out, flows)),
            "is_correlation": float(self.use_correlation),
            "Cd_in_range": float(self._cd_in_range(state_in)),
        }


def _effusion_wall(
    h_i: float, T_c: float, h_gas: float, T_aw_gas: float, t_over_k: float
) -> tuple[float, float, float]:
    """(T_wall_hot, T_wall_cold, q) of a plate between gas and coolant.

    Per unit area, one conduction layer, through C++'s wall coupling. The one
    wall core both EffusionPlateElement.overall_effectiveness (no conduction)
    and wall_heat_transfer use.
    """
    w = _solver_tools.wall_coupling_and_jacobian_multilayer(
        h_gas, T_aw_gas, h_i, T_c, [t_over_k], 1.0, 0.0
    )
    q = float(w.Q)
    return float(w.T_hot), T_c + q / h_i, q


def _chamber_gas_side(node: Any, state: "NetworkMixtureState") -> dict[str, float] | None:
    """What a discharge node offers a wall as its gas side, or None.

    Only a MomentumChamberNode has flow, hence a coefficient: its own surface
    correlation (the user's choice on that node) is the UNBLOWN h, and its
    velocity the gas velocity. The gas TEMPERATURE is the stream the wall
    sees: a merge chamber's state is the MIXED outlet, already diluted by
    the side streams -- with liner effusion at 20-30% of the flow that is a
    large, unconservative drop -- so when the chamber has a main inlet the
    approaching main stream's temperature is used instead. A plenum (state,
    no flow) returns None.
    """
    if not isinstance(node, MomentumChamberNode) or node.area <= 0.0:
        return None
    gas = node.htc_and_T(state)
    if gas is None:
        return None
    rho_g, _ = _safe_rho(cb.density(state.T, state.P, state.X))
    T_gas, from_main = float(gas.T_aw), False
    if node.main_inlet is not None and node.has_side_streams():
        for eid, st, _names in getattr(node, "_inflows", []):
            if eid == node.main_inlet:
                T_gas, from_main = float(st.T), True
    return {
        "h_unblown": float(gas.h),
        "T_gas": T_gas,
        "T_from_main_inlet": from_main,
        "U_gas": abs(getattr(node, "_total_m_dot", 0.0)) / (rho_g * node.area),
        "rho_gas": rho_g,
    }


class EffusionPlateElement(OrificeElement):
    """A multi-perforated (effusion) wall panel: N holes discharging in parallel.

    Coolant enters from `from_node` and leaves through the wall into
    `to_node`, so mass leaves the coolant circuit -- which is what makes this
    an element rather than a `ConvectiveSurface` model.

    GEOMETRY IS GIVEN THE WAY A PLATE IS DESIGNED, not as an area: pitch,
    hole diameter, wall thickness and inclination angle. The hole count
    follows from the panel area and the pitch, and is rounded to a whole
    number of holes -- `hole_count_exact` and `pitch_actual` record what the
    rounding did, because the rounded count is what the element flows.

        n_holes  = round(panel_area / (pitch_x * pitch_y))
        L_hole   = wall_thickness / sin(angle)     <- the drilled length
        A_total  = n_holes * pi * d^2 / 4
        porosity = A_total / panel_area

    HOMOGENISATION, AND WHEN IT BREAKS. One element is one panel with one
    coolant pressure and one gas pressure, so it cannot represent coolant
    migration WITHIN itself. van de Noort and Ireland (2022) show this is not
    a small effect: with uniform inlet AND outlet pressure, holes at the end
    of an array can pass ~75% of what a central hole passes, and under a
    spanwise pressure gradient one row took ~10% of the total while others
    took ~20% each.

    The answer is to use more panels, not a cleverer one. Because the solver
    balances mass over every element at a node, a panel hung off each segment
    of a coolant channel reproduces their flow network directly:

        coolant:  P0 --[Channel]-- n1 --[Channel]-- n2 --[Channel]-- ...
                                    |               |
                            [EffusionPanel]  [EffusionPanel]
                                    |               |
        gas:                       g1              g2

    At each node, `m_channel_in = m_channel_out + m_effusion`, which is the
    coolant mass flow falling along the wall -- the effusion CHANNEL case,
    with no channel-specific element needed. Resolution is the caller's
    choice of segment count. Compare `dP_drive` across panels to see whether
    one panel is smearing a gradient that deserves several.

    INGESTION. If the gas pressure exceeds the coolant pressure the panel
    ingests hot gas. That is a real failure mode (van de Noort's CMF > 0.5),
    and a homogenised panel would otherwise average it into a healthy net
    outflow, so `diagnostics()` reports `is_ingesting` rather than staying
    silent about it.

    DISCHARGE COEFFICIENT. Default `IdelchikThick`: diagram 4-18a is a
    thick-walled hole in a large wall between two plena, which is exactly a
    plenum-fed effusion plate, and it is valid down to Re = 25. For a panel
    fed by a channel rather than a plenum the approach flow is a CROSSFLOW,
    which is McGreehan-Schotsch's `U1/Vi` term -- use `'McGreehanSchotsch'`
    there. Note that neither carries a velocity-of-approach `beta` factor:
    a plenum has no approach velocity to correct for. van de Noort applies
    one built on a half-pitch square inlet area, which is worth knowing about
    but is a different convention; at their own pitch it is a 0.3% effect.

    THERMAL, INTERNAL SIDE. `internal_heat_transfer()` gives the coolant-side
    coefficient from Andrews 86-GT-225: the hole APPROACH flow over the
    coolant-side surface plus the THROAT, summed. The approach term dominates
    at effusion Reynolds numbers, which is the paper's own headline and the
    reason a throat-only treatment recovers only a fraction of the measured
    coefficient.

    THE PLATE OWNS ITS WALL (#471). `wall_heat_transfer()` -- reported in
    the diagnostics -- solves the wall with the coolant side above and a gas
    side the DISCHARGE NODE decides: a MomentumChamberNode (flow) supplies
    its own unblown coefficient, gas temperature and velocity; a plenum
    (state, no flow) takes the imposed `gas_heat_flux`. No ThermalWall
    connects to this element.

    AN OVERALL EFFECTIVENESS CORRELATION IS NEVER THE CLOSURE. Andrews
    (88-GT-290) states internal and film cooling are NOT additive -- "the
    film cooling reduces the mean gas temperature adjacent to the wall...
    which in turn reduces the heat flux removed by the internal wall
    cooling" -- and an overall correlation already contains the internal
    convection computed here. The overall effectiveness is an OUTPUT. An
    adiabatic film (Baldauf + Sellers, `gas_film='baldauf_sellers'`) is
    selectable, off by default: scored on Andrews, the data refuse it offered
    alone; the missing physics is the gas-side augmentation
    (`gas_augmentation`, the caller's).

    PLENUM-FED. This 2-port plate leaves McGreehan-Schotsch's U1/Vi at 0 and
    uses Andrews' still-plenum coolant side; `coolant_crossflow_ignored`
    flags a supply node a channel runs through.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        hole_diameter: float,
        wall_thickness: float,
        pitch: float | None = None,
        panel_area: float | None = None,
        pitch_x: float | None = None,
        pitch_y: float | None = None,
        panel_length: float | None = None,
        panel_width: float | None = None,
        angle_deg: float = 90.0,
        correlation: str = "IdelchikThick",
        Cd: float = 0.6,
        edge_radius: float = 0.0,
        internal_Nu_multiplier: float = 1.0,
        wall_conductivity: float = 20.0,
        gas_film: Literal["none", "baldauf_sellers"] = "none",
        gas_augmentation: float = 1.0,
        turbulence_intensity: float = 0.05,
        gas_heat_flux: float = 0.0,
        crossflow_segments: tuple[str | None, str | None] | None = None,
        crossflow_area: float | None = None,
    ) -> None:
        if hole_diameter <= 0.0:
            raise ValueError("EffusionPlateElement: hole_diameter must be positive")
        if crossflow_segments is not None and not (crossflow_area and crossflow_area > 0.0):
            raise ValueError(
                "EffusionPlateElement: a crossflow-fed panel needs the backside "
                "duct's crossflow_area"
            )
        # DUCT-FED (#471): the panel's supply node is a station of a backside
        # duct; crossflow_segments = (arriving, leaving) segment ids there. The
        # duct's mean velocity at the station, over the static-referenced ideal
        # jet velocity, is McGreehan-Schotsch's U1/Vi.
        self.crossflow_segments = crossflow_segments
        self.crossflow_area = crossflow_area
        self._xf_elems: list[tuple[NetworkElement, str]] = []
        if wall_conductivity <= 0.0:
            raise ValueError("EffusionPlateElement: wall_conductivity must be positive")
        if gas_augmentation <= 0.0:
            raise ValueError("EffusionPlateElement: gas_augmentation must be positive")
        if gas_film not in ("none", "baldauf_sellers"):
            raise ValueError("EffusionPlateElement: gas_film is 'none' or 'baldauf_sellers'")
        if gas_film == "baldauf_sellers" and panel_length is None:
            raise ValueError(
                "EffusionPlateElement: gas_film='baldauf_sellers' needs panel_length "
                "to count the hole rows the film builds over"
            )
        # The plate's own wall (#471): conduction through it, and the gas
        # side the discharge node decides -- see wall_heat_transfer().
        self.wall_conductivity = float(wall_conductivity)
        self.gas_film = gas_film
        self.gas_augmentation = float(gas_augmentation)
        self.turbulence_intensity = float(turbulence_intensity)
        self.gas_heat_flux = float(gas_heat_flux)
        self.panel_length = panel_length
        if internal_Nu_multiplier <= 0.0:
            raise ValueError("EffusionPlateElement: internal_Nu_multiplier must be positive")
        if wall_thickness <= 0.0:
            raise ValueError("EffusionPlateElement: wall_thickness must be positive")
        if not 0.0 < angle_deg <= 90.0:
            raise ValueError(
                "EffusionPlateElement: angle_deg must be in (0, 90]; it is the "
                "hole inclination to the wall PLANE, 90 being a normal hole"
            )

        # The coolant-side tuner, matching `Nu_multiplier` on ConvectiveSurface.
        # Andrews 86-GT-225's correlation runs -10.4% against his own Fig. 8,
        # so a user matching their own plate needs this; the harness never
        # sets, recommends or scores it. See docs/VALIDATION_POLICY.md.
        self.internal_Nu_multiplier = float(internal_Nu_multiplier)

        px = pitch_x if pitch_x is not None else pitch
        py = pitch_y if pitch_y is not None else pitch
        if px is None or py is None:
            raise ValueError(
                "EffusionPlateElement: give pitch (square array) or both pitch_x and pitch_y"
            )
        if px <= 0.0 or py <= 0.0:
            raise ValueError("EffusionPlateElement: pitch must be positive")

        if panel_area is None:
            if panel_length is None or panel_width is None:
                raise ValueError(
                    "EffusionPlateElement: give panel_area, or both panel_length and panel_width"
                )
            panel_area = panel_length * panel_width
        if panel_area <= 0.0:
            raise ValueError("EffusionPlateElement: panel_area must be positive")

        cell_area = px * py
        n_exact = panel_area / cell_area
        n_holes = int(round(n_exact))
        if n_holes < 1:
            raise ValueError(
                f"EffusionPlateElement: the panel holds {n_exact:.3g} holes at "
                f"this pitch, which rounds to none. Enlarge the panel or "
                f"reduce the pitch."
            )

        self.hole_diameter = hole_diameter
        self.wall_thickness = wall_thickness
        self.pitch_x = px
        self.pitch_y = py
        self.panel_area = panel_area
        self.angle_deg = angle_deg
        self.n_holes = n_holes
        # What the rounding cost: the caller asked for a pitch, and the whole
        # number of holes implies a slightly different one.
        self.hole_count_exact = n_exact
        self.pitch_actual = math.sqrt(panel_area / n_holes)

        # Drilled length along the hole axis. An inclined hole is longer than
        # the wall is thick, which is the whole reason effusion holes are
        # inclined: more internal surface for the same wall.
        alpha = math.radians(angle_deg)
        self.hole_length = wall_thickness / math.sin(alpha)

        hole_area = math.pi * hole_diameter * hole_diameter / 4.0
        self.total_hole_area = n_holes * hole_area
        self.porosity = self.total_hole_area / panel_area

        super().__init__(
            id,
            from_node,
            to_node,
            Cd=Cd,
            area=self.total_hole_area,
            correlation=correlation,
            plate_thickness=self.hole_length,
            edge_radius=edge_radius,
        )
        # OrificeElement derived an equivalent single-bore diameter from the
        # TOTAL area. Keep it for the flow equation, but the correlations must
        # see one real hole -- see _hole_diameter_for_correlation.
        self.equivalent_bore = self.diameter
        # Discharged into a merge chamber, the holes' inclination to the wall
        # is the jets' angle to the gas axis (holes pointing downstream).
        self.injection_angle_deg = float(angle_deg)

    def _hole_count(self) -> float:
        return float(self.n_holes)

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """No upstream-diameter discovery: a wall panel has no pipe, hence no
        beta. The geometry the correlations need is the single hole, which is
        known at construction. The discharge node decides the gas side."""
        self._orifice_geom = cb.OrificeGeometry()
        self._orifice_geom.d = self.hole_diameter
        self._orifice_geom.D = 0.0
        self._orifice_geom.t = self.hole_length
        self._orifice_geom.r = self.edge_radius
        self.beta = 0.0
        self._discharge_node = graph.nodes.get(self.to_node)
        # A supply node a channel runs THROUGH feeds the holes with a
        # crossflow. This 2-port plate is plenum-fed: McGreehan-Schotsch's
        # U1/Vi stays 0 and Andrews' coolant side assumes a still plenum, so
        # say so rather than model it (the 3-port liner segment will).
        through = (
            [
                e
                for e in graph.get_upstream_elements(self.from_node)
                + graph.get_downstream_elements(self.from_node)
                if isinstance(e, ChannelElement)
            ]
            if self.from_node in graph.nodes
            else []
        )
        self._coolant_crossflow_ignored = len(through) >= 2 and self.crossflow_segments is None
        self._xf_elems = []
        if self.crossflow_segments is not None:
            for eid in self.crossflow_segments:
                if eid is not None:
                    self._xf_elems.append((graph.elements[eid], eid))

    def network_flow_inputs(self) -> list[tuple[str, str]]:
        """The backside segments' flows at this panel's supply station."""
        return [(eid, self.from_node) for _, eid in self._xf_elems]

    def _supply_crossflow(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None,
    ) -> tuple[float, dict[str, float]]:
        if self.crossflow_segments is None:
            return 0.0, {}
        flows = flows or {}
        prev_id, next_id = self.crossflow_segments
        m_a = flows.get(prev_id, 0.0) if prev_id else 0.0
        m_b = flows.get(next_id, 0.0) if next_id else 0.0
        r = cb.crossflow_velocity_ratio(
            m_a, m_b, state_in.P, state_in.T, state_in.X, state_out.P, self.crossflow_area
        )
        d: dict[str, float] = {}
        for elem, eid in self._xf_elems:
            coeff_m = r.d_dm_a if eid == prev_id else r.d_dm_b
            names = elem.unknowns()
            for local, c in elem.flow_jac_at_node(self.from_node, list(range(len(names)))).items():
                d[names[local]] = d.get(names[local], 0.0) + coeff_m * c
        for key, v in (
            (f"{self.from_node}.P", r.d_dP),
            (f"{self.from_node}.T", r.d_dT),
            (f"{self.to_node}.P", r.d_dP_down),
        ):
            d[key] = d.get(key, 0.0) + v
        return float(r.U1_over_Vi), d

    def validate(self) -> None:
        super().validate()
        if self.correlation in ("ReaderHarrisGallagher", "Stolz", "Miller"):
            raise ValueError(
                f"EffusionPlateElement {self.id!r}: {self.correlation!r} is a "
                "NORMED metering correlation for a standardised plate in a "
                "pipe, referenced to the tapping differential and a function "
                "of beta = d/D. An effusion panel has no pipe and no beta. "
                "Use 'IdelchikThick' (plenum-fed), 'McGreehanSchotsch' "
                "(channel-fed, with approach crossflow), or 'fixed'."
            )
        if self.porosity >= 1.0:
            raise ValueError(
                f"EffusionPlateElement {self.id!r}: porosity is "
                f"{self.porosity:.3f}; the holes do not fit in the panel."
            )
        if self.crossflow_segments is not None and self.correlation != "McGreehanSchotsch":
            raise ValueError(
                f"EffusionPlateElement {self.id!r}: a duct-fed panel needs the "
                "supply-side crossflow term, which only 'McGreehanSchotsch' has "
                f"({self.correlation!r} would silently treat the duct as a plenum)."
            )

    def internal_heat_transfer(self, state_in: NetworkMixtureState) -> dict[str, float]:
        """Coolant-side heat transfer for this panel, Andrews 86-GT-225.

        The hole approach and the throat, summed (Eq. 19). Returns the two
        Nusselt numbers, the Reynolds number they share, and the coefficient
        expressed BOTH ways -- because which area an `h` belongs to is
        exactly the kind of thing that goes wrong silently:

          `h_hole_area`  : on the hole internal surface, pi d L_hole
          `h_plate_area` : on the approach area per hole, cell - hole, which
                           is Andrews' own definition and what his Fig. 8
                           plots

        The two differ by `area_ratio` -- 3.46 for Andrews' plate C -- so
        using one where the other is meant is a factor-of-three error, not a
        refinement. For a panel energy balance use
        `Q = h_plate_area * area_approach_total * dT`.

        Properties are taken at the COOLANT inlet state. The coolant heats up
        through the hole, so a large temperature rise makes this a first
        approximation; Andrews' own rig was near-ambient.

        NON-SQUARE ARRAYS. Andrews' plates are square and Eq. (18) carries
        the pitch explicitly, so a rectangular array is represented by its
        equivalent square pitch, `pitch_actual`. Fine while the aspect ratio
        is near one; a strongly rectangular array is outside what the
        correlation was fitted on.
        """
        ts = cb.transport_state(state_in.T, state_in.P, state_in.X)
        mu = ts.mu if ts.mu > 0.0 else 1.8e-5
        k = ts.k
        pr = ts.Pr

        m_hole = abs(state_in.m_dot) / self.n_holes
        if m_hole <= 0.0:
            return {}
        Re = 4.0 * m_hole / (math.pi * self.hole_diameter * mu)

        X = self.pitch_actual
        nu_approach = cb.effusion_approach_nusselt(Re, pr, X / self.hole_length)
        nu_throat = cb.effusion_throat_nusselt(Re, pr, self.hole_length / self.hole_diameter)

        area_hole = math.pi * self.hole_diameter * self.hole_length
        area_approach = self.panel_area / self.n_holes - math.pi * (self.hole_diameter**2) / 4.0

        # `internal_Nu_multiplier` is the user's rig-matching knob and is 1.0
        # unless they set it. Applied to the SUM, so the split between the
        # approach and throat terms stays the correlation's.
        nu_internal = (nu_approach + nu_throat) * self.internal_Nu_multiplier
        h_hole = nu_internal * k / self.hole_diameter
        return {
            "Re_hole": float(Re),
            "Nu_approach": float(nu_approach),
            "Nu_throat": float(nu_throat),
            "Nu_internal": float(nu_internal),
            "internal_Nu_multiplier": float(self.internal_Nu_multiplier),
            "h_hole_area": float(h_hole),
            "h_plate_area": float(h_hole * area_hole / area_approach),
            "area_ratio": float(area_approach / area_hole),
            "area_approach_total": float(area_approach * self.n_holes),
            "area_hole_total": float(area_hole * self.n_holes),
        }

    def overall_effectiveness(
        self,
        state_in: NetworkMixtureState,
        h_gas_unblown: float,
        T_gas: float,
        U_gas: float | None = None,
        gas_augmentation: float = 1.0,
        eta_film: float = 0.0,
    ) -> dict[str, float]:
        """Overall cooling effectiveness of this panel, as an OUTPUT.

            h_gas = h_gas_unblown * gas_augmentation
            eta   = (h_i + h_gas eta_film) / (h_i + h_gas)
            T_wall = T_gas - eta (T_gas - T_coolant)

        The gas side is the caller's; `eta` and the wall temperature are
        the outputs. An overall-effectiveness CORRELATION must never be the
        closure here -- it would already contain the internal convection
        this computes, and the two would double-count.

        THE FILM IS TWO NUMBERS, NOT ONE, and the signature says so. The
        standard film-cooling form is `q = h_f (T_aw - T_w)` with
        `T_aw = T_gas - eta_film (T_gas - T_c)`, so a film owes you both a
        driving temperature (`eta_film`) AND a conductance
        (`gas_augmentation = h_f/h_0`). An adiabatic effectiveness
        correlation such as `film_effectiveness_baldauf_2002` supplies only
        the first: an adiabatic wall passes no heat, so the experiment
        fixes `T_aw` and measures no coefficient. Supplying `eta_film`
        while leaving `gas_augmentation` at 1.0 therefore OVER-PREDICTS,
        and both default to their no-film values so that omitting the pair
        is consistent rather than half-right.

        `gas_augmentation` IS THE CALLER'S, NEVER FITTED HERE. It matches
        `Nu_multiplier` on `ConvectiveSurface`: the library states what the
        correlations give, and matching a specific rig is the user's job.
        See docs/VALIDATION_POLICY.md.

        MIND WHAT `h_gas_unblown` IS MEASURED AGAINST. An augmentation
        ratio is meaningless without its baseline, and mixing baselines is
        a factor-level error, not a refinement. For Andrews 88-GT-290's rig
        the required `gas_augmentation` is 3.1-4.3 against a fully
        developed Dittus-Boelter and 1.6-2.5 against the same duct with a
        thermal-entry correction -- a factor of 1.75 from that choice
        alone. Published film-cooling ratios are usually referenced to a
        flat-plate turbulent boundary layer at the same x, which is a third
        baseline again.

        WHAT SCORING IT AGAINST ANDREWS SHOWED. At `gas_augmentation = 1.0`
        and `eta_film = 0.0` the closure scores +3.1% on his effusion plate
        C and +23.7% on plate B -- and plate C's agreement is a
        CANCELLATION, a real film raising `eta` against its augmentation
        lowering it, not evidence that either is absent. The two plates
        need gas-side coefficients differing by 1.7x, and that ratio
        survives any choice of film model (scaling Baldauf's `eta_film`
        from 0 to 1.25x moves it only 1.96 to 1.60). So the missing physics
        is a `gas_augmentation` correlation for full-coverage effusion,
        which needs a HEATED-wall measurement; a better film correlation
        cannot supply it. See
        `validation/cooling/extractions/andrews_effusion_overall_eta.md`.

        Pass `U_gas` to get the jet ratios reported beside the result.
        They are REPORTED AND NEVER APPLIED -- no threshold is offered,
        because the requirement collapses on none of the velocity ratio,
        the blowing ratio or the momentum flux ratio. If `velocity_ratio`
        exceeds about 1 -- Andrews' plate B ejects at 53 m/s into a
        26.8 m/s crossflow -- treat `eta` as an upper bound unless
        `gas_augmentation` already accounts for it.

        Returns `{}` when there is no coolant flow or no gas side. With
        `U_gas` omitted the jet ratios are absent rather than guessed.
        """
        internal = self.internal_heat_transfer(state_in)
        if not internal or h_gas_unblown <= 0.0 or gas_augmentation <= 0.0:
            return {}
        if not 0.0 <= eta_film < 1.0:
            raise ValueError(
                "EffusionPlateElement: eta_film must be in [0, 1); it is an "
                "adiabatic effectiveness, not an overall one"
            )

        h_i = internal["h_plate_area"]
        h_gas = h_gas_unblown * gas_augmentation
        T_c = state_in.T
        # The two-resistance closure eta = (h_i + h_gas eta_film)/(h_i + h_gas)
        # is the wall core below with no conduction resistance.
        T_aw = T_gas - eta_film * (T_gas - T_c)
        t_wall, _, _ = _effusion_wall(h_i, T_c, h_gas, T_aw, 0.0)
        eta = (T_gas - t_wall) / (T_gas - T_c) if T_gas != T_c else 0.0
        out = {
            "eta_overall": float(eta),
            "T_wall": float(t_wall),
            "T_adiabatic_wall": float(T_gas - eta_film * (T_gas - T_c)),
            "h_internal_plate_area": float(h_i),
            "internal_Nu_multiplier": float(self.internal_Nu_multiplier),
            "h_gas_unblown": float(h_gas_unblown),
            "gas_augmentation": float(gas_augmentation),
            "h_gas": float(h_gas),
            "eta_film": float(eta_film),
            "resistance_ratio": float(h_i / h_gas),
            "q_flux": float(h_gas * (T_gas - eta_film * (T_gas - T_c) - t_wall)),
        }
        if U_gas is None or U_gas <= 0.0:
            return out

        rho_c = cb.density(T_c, state_in.P, state_in.X)
        area_hole = math.pi * (self.hole_diameter**2) / 4.0
        mass_flux_c = abs(state_in.m_dot) / self.n_holes / area_hole
        u_jet = mass_flux_c / rho_c
        rho_g = rho_c * T_c / T_gas  # same static pressure, ideal gas
        blowing = mass_flux_c / (rho_g * U_gas)
        out.update(
            {
                "u_jet": float(u_jet),
                "velocity_ratio": float(u_jet / U_gas),
                "blowing_ratio": float(blowing),
                "momentum_flux_ratio": float(blowing * blowing / (rho_c / rho_g)),
                "density_ratio": float(rho_c / rho_g),
            }
        )
        return out

    def wall_heat_transfer(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        """The plate's own wall, with the gas side its discharge node decides.

        The plate IS the wall the coolant flows through, so it owns it: no
        ThermalWall connects to it. Its heat goes from the gas into the
        effusing coolant, which discharges back into the same node, so the
        network's energy balance is unchanged; the wall temperatures and the
        heat flux are outputs.

        COOLANT SIDE: Andrews 86-GT-225 on the plate area (`internal_heat_
        transfer`), at the supply temperature.

        GAS SIDE, by discharge node:

        * a MomentumChamberNode has flow. Its own surface correlation gives
          the UNBLOWN coefficient and gas temperature, its velocity the
          blowing ratio. `gas_augmentation` (the caller's, 1.0 unless set;
          never fitted -- see docs/VALIDATION_POLICY.md) scales the
          coefficient; Andrews 88-GT-290 shows that augmentation, not a
          film correlation, is the missing physics. `gas_film` stays 'none'
          by default for the same reason: Baldauf + Sellers
          ('baldauf_sellers') is selectable and reported with its envelope
          flag, but the data refuse it offered alone.
        * anything else (a plenum) has a state but no flow, so no
          coefficient: the gas side is the imposed `gas_heat_flux` [W/m^2],
          0 by default (adiabatic). Radiation is not modelled.

        WALL: one layer, `wall_thickness / wall_conductivity`, through C++'s
        wall_coupling_and_jacobian. Returns {} without coolant flow.
        """
        internal = self.internal_heat_transfer(state_in)
        if not internal:
            return {}
        h_i = internal["h_plate_area"]
        T_c = float(state_in.T)
        A = self.panel_area
        t_over_k = self.wall_thickness / self.wall_conductivity
        out: dict[str, float] = {
            "h_internal_plate_area": float(h_i),
            "wall_conductivity": float(self.wall_conductivity),
            "coolant_crossflow_ignored": float(getattr(self, "_coolant_crossflow_ignored", False)),
            # Duct-fed: the Cd sees the crossflow, but the coolant-side heat
            # transfer is still Andrews' still-plenum correlation (provisional
            # until a crossflow-supply source; #471 PR4).
            "coolant_ht_plenum_assumed": float(self.crossflow_segments is not None),
        }
        gas = _chamber_gas_side(getattr(self, "_discharge_node", None), state_out)
        if gas is None:
            q = self.gas_heat_flux
            t_cold = T_c + q / h_i
            out.update(
                {
                    "gas_side_chamber": 0.0,
                    "gas_heat_flux_imposed": float(q),
                    "q_wall": float(q),
                    "Q_wall": float(q * A),
                    "T_wall_cold": float(t_cold),
                    "T_wall_hot": float(t_cold + q * t_over_k),
                }
            )
            return out

        rho_c, _ = _safe_rho(cb.density(T_c, state_out.P, state_in.X))
        U_g, rho_g, T_g = gas["U_gas"], gas["rho_gas"], gas["T_gas"]
        hole_area = math.pi * self.hole_diameter**2 / 4.0
        G_jet = abs(state_in.m_dot) / (self.n_holes * hole_area)
        blowing = G_jet / (rho_g * U_g) if U_g > 0.0 else math.inf
        eta_film, film_extrapolated = 0.0, False
        if self.gas_film == "baldauf_sellers" and math.isfinite(blowing):
            film = cb.effusion_panel_film_effectiveness(
                max(1, round(self.panel_length / self.pitch_x)),
                self.pitch_x / self.hole_diameter,
                self.pitch_y / self.hole_diameter,
                blowing,
                rho_c / rho_g,
                self.angle_deg,
                self.turbulence_intensity,
            )
            eta_film, film_extrapolated = float(film.eta), bool(film.extrapolated)
        T_aw_g = T_g - eta_film * (T_g - T_c)
        h_g = gas["h_unblown"] * self.gas_augmentation
        t_hot, t_cold, q = _effusion_wall(h_i, T_c, h_g, T_aw_g, t_over_k)
        out.update(
            {
                "gas_side_chamber": 1.0,
                "h_gas_unblown": float(gas["h_unblown"]),
                "gas_augmentation": float(self.gas_augmentation),
                "h_gas": float(h_g),
                "T_gas": float(T_g),
                "T_gas_from_main_inlet": float(gas["T_from_main_inlet"]),
                "eta_film": eta_film,
                "film_extrapolated": float(film_extrapolated),
                "T_adiabatic_wall": float(T_aw_g),
                "q_wall": float(q),
                "Q_wall": float(q * A),
                "T_wall_hot": float(t_hot),
                "T_wall_cold": float(t_cold),
                "eta_overall": float((T_g - t_hot) / (T_g - T_c)) if T_g != T_c else 0.0,
                "U_gas": float(U_g),
                "blowing_ratio": float(blowing),
                "density_ratio": float(rho_c / rho_g),
                "velocity_ratio": float(G_jet / rho_c / U_g) if U_g > 0.0 else math.inf,
            }
        )
        return out

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float]:
        base = super().diagnostics(state_in, state_out, flows)
        m_dot = abs(state_in.m_dot)
        dP_drive = state_in.Pt - state_out.P
        if self.crossflow_segments is not None:
            u1_vi, _ = self._supply_crossflow(state_in, state_out, flows)
            # Rohde-scored bounds of McGreehan-Schotsch's crossflow term
            # (orifice_discharge_coefficient.md): about +8% while U1/Vi < 0.2,
            # degrading badly beyond 0.35, worst for sharp edges.
            base.update(
                {
                    "U1_over_Vi": float(u1_vi),
                    "crossflow_cd_beyond_8pct": float(u1_vi > 0.2),
                    "crossflow_cd_degraded": float(u1_vi > 0.35),
                }
            )
        out = {
            "n_holes": float(self.n_holes),
            "hole_count_exact": float(self.hole_count_exact),
            "porosity": float(self.porosity),
            "hole_length": float(self.hole_length),
            "L_over_d": float(self.hole_length / self.hole_diameter),
            "pitch_over_d": float(self.pitch_actual / self.hole_diameter),
            # Andrews et al. (1988) correlate effusion cooling on the coolant
            # mass flow per unit PLATE area, not per hole. Their measurements
            # span G = 0.1 to 1.6 kg/s/m^2.
            "G_coolant": float(m_dot / self.panel_area),
            "dP_drive": float(dP_drive),
            # A homogenised panel would otherwise average an ingesting hole
            # into a healthy net outflow. Say it instead.
            "is_ingesting": float(dP_drive <= 0.0),
        }
        # Coolant-side heat transfer, when there is flow to carry it, and the
        # plate's own wall with the gas side its discharge node decides.
        out.update(self.internal_heat_transfer(state_in))
        out.update(self.wall_heat_transfer(state_in, state_out))
        out.update(base)
        return out


@dataclass
class _ImpingementPlateResult:
    """What ``ImpingementPlateElement.htc_and_T`` returns to the wall code.

    ``dh_dsources`` is the part no other element has: the derivative of ``h``
    with respect to each crossflow source's mass flow at this plate's
    ``to_node``, keyed by element id. The solver relays it through the
    element's ``network_flow_inputs()`` (#465).
    """

    h: float
    Nu: float
    Re: float
    Pr: float
    T_aw: float
    Gc_Gj: float
    extrapolated: bool = False
    dh_dmdot: float = 0.0
    dh_dT: float = 0.0
    dT_aw_dmdot: float = 0.0
    dT_aw_dT: float = 1.0
    dh_dsources: dict[str, float] = field(default_factory=dict)


class ImpingementPlateElement(OrificeElement):
    """One spanwise row of an impingement jet plate, with its target wall (#465).

    The element's mass flow IS the row's jet flow: it leaves the supply plenum
    (``from_node``) through ``n_holes`` holes and joins the crossflow channel
    between jet plate and target (``to_node``). So the same element carries
    the orifice flow AND the Florschuetz, Truman and Metzger (1981) heat
    transfer on the target, which is what ``ImpingementModel`` on a
    ``ChannelElement`` cannot do: there the element's flow is the channel's.

    A full array is a chain, one plate per row:

        supply plenum ---+---------------+---------------+
                         |               |               |
                    [Plate 1]       [Plate 2]       [Plate 3]
                         |               |               |
        crossflow:      c1 --[Channel]-- c2 --[Channel]-- c3 --[Channel]-- exit

    CROSSFLOW FROM THE NETWORK. Florschuetz's crossflow-to-jet mass velocity
    ratio is

        Gc/Gj = (m_c / m_j) (pi/4) / ((yn/d)(z/d))

    with ``m_j`` this row's jet flow and ``m_c`` the crossflow APPROACHING the
    row: every inflow to ``to_node`` except this element's own. In the chain
    above that is the channel from the previous crossflow node, which carries
    rows 1..i-1 -- Florschuetz's Eq. 8 evaluated at x - xn/2. The ratio is
    read from the solved flows, not from the uniform-supply closed form, so a
    non-uniform supply, a row of different geometry, or a bleed shows up in
    the heat transfer. ``h`` therefore depends on a NEIGHBOUR's mass flow; the
    element exposes that as ``dh_dsources`` and the solver relays it into the
    Jacobian.

    An INITIAL crossflow (flow entering upstream of row 1) is representable
    but outside the source: Florschuetz et al. (1981) had none. ``Gc/Gj``
    beyond 0.8 is flagged by the set's own range check.

    THE TARGET. ``surface.area`` is the target footprint of the holes,
    ``n_holes * xn * yn``. ``T_aw`` is the supply plenum temperature, which
    is Florschuetz's own reference temperature for h (plenum-fed, so static
    and total coincide). The heat leaves the wall into ``to_node``: the
    spent air carries it downstream.

    DISCHARGE COEFFICIENT. Default ``'fixed'`` at Florschuetz's own 0.79
    (measured 0.73-0.85 across their plates). The hole Cd does not depend on
    the crossflow: McGreehan-Schotsch's ``U1/Vi`` is the SUPPLY-side
    approach velocity, zero for a plenum, never Gc/Gj. Selectable
    alternatives, never chosen automatically: ``'IdelchikThick'``,
    ``'Lichtarowicz'`` (long hole, l/d 2-10) and ``'McGreehanSchotsch'``.
    ``Cd_in_range`` in the diagnostics says whether the plate is inside the
    chosen correlation's source range.

    NOT MODELLED. The temperature sensitivity of the correlation's
    properties (``dh_dT`` is 0, the same documented gap the channel
    impingement path has), the crossflow's own temperature (Florschuetz
    referenced h to the plenum), and the jet plate's own heat pick-up.

    Parameters
    ----------
    d_jet : float
        Hole diameter [m].
    xn_d, yn_d, z_d : float
        Streamwise pitch, spanwise pitch and plate-to-target gap over d_jet.
    span : float
        Plate span [m]; ``n_holes = round(span / (yn_d d_jet))``.
    plate_thickness : float
        Jet plate thickness [m]: the hole length the Cd correlations read.
    pattern : {'inline', 'staggered'}
        Selects Florschuetz's inline or staggered set when
        ``correlation_set`` is None.
    row : int or None
        Only for diagnostics: when given, ``Gc_Gj_closed_form`` reports Eq. 8
        at this row beside the network value.
    Nu_multiplier : float
        The user's rig-matching knob. 1.0 unless they set it.
    """

    # htc_and_T's dh_dmdot comes from C++ in the signed jet flow (#481).
    _htc_dmdot_is_signed = True

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        d_jet: float,
        xn_d: float,
        yn_d: float,
        z_d: float,
        span: float,
        plate_thickness: float,
        pattern: Literal["inline", "staggered"] = "inline",
        correlation: str = "fixed",
        Cd: float = cb.FLORSCHUETZ_1981_DEFAULT_CD,
        correlation_set: object | None = None,
        row: int | None = None,
        Nu_multiplier: float = 1.0,
        edge_radius: float = 0.0,
    ) -> None:
        for name, val in (
            ("d_jet", d_jet),
            ("xn_d", xn_d),
            ("yn_d", yn_d),
            ("z_d", z_d),
            ("span", span),
            ("plate_thickness", plate_thickness),
            ("Nu_multiplier", Nu_multiplier),
        ):
            if val <= 0.0:
                raise ValueError(f"ImpingementPlateElement: {name} must be positive")
        if pattern not in ("inline", "staggered"):
            raise ValueError("ImpingementPlateElement: pattern is 'inline' or 'staggered'")

        n_exact = span / (yn_d * d_jet)
        n_holes = int(round(n_exact))
        if n_holes < 1:
            raise ValueError(
                f"ImpingementPlateElement: the span holds {n_exact:.3g} holes at this "
                "spanwise pitch, which rounds to none."
            )

        self.d_jet = d_jet
        self.xn_d = xn_d
        self.yn_d = yn_d
        self.z_d = z_d
        self.span = span
        self.pattern = pattern
        self.row = row
        self.Nu_multiplier = float(Nu_multiplier)
        self.n_holes = n_holes
        self.hole_count_exact = n_exact
        if correlation_set is None:
            correlation_set = (
                cb.florschuetz_1981_inline()
                if pattern == "inline"
                else cb.florschuetz_1981_staggered()
            )
        self.correlation_set = correlation_set

        super().__init__(
            id,
            from_node,
            to_node,
            Cd=Cd,
            area=n_holes * math.pi * d_jet * d_jet / 4.0,
            correlation=correlation,
            plate_thickness=plate_thickness,
            edge_radius=edge_radius,
        )
        # The target footprint of the holes this element flows.
        self.surface = ConvectiveSurface(area=n_holes * xn_d * yn_d * d_jet * d_jet)
        self._row_geometry = cb.JetRowGeometry(
            d=d_jet, xn_d=xn_d, yn_d=yn_d, z_d=z_d, n_holes=float(n_holes)
        )
        self._crossflow_sources: list[str] = []

    @property
    def has_convective_surface(self) -> bool:
        return True

    def _hole_count(self) -> float:
        return float(self.n_holes)

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """A plate has no pipe, so no beta; the correlations see one hole.
        The crossflow sources are every other inflow to ``to_node``."""
        self._orifice_geom = cb.OrificeGeometry()
        self._orifice_geom.d = self.d_jet
        self._orifice_geom.D = 0.0
        self._orifice_geom.t = self.plate_thickness
        self._orifice_geom.r = self.edge_radius
        self.beta = 0.0
        self._crossflow_sources = [
            e.id for e in graph.get_upstream_elements(self.to_node) if e.id != self.id
        ]
        self._topology_error = self._configuration_error(graph)

    def _configuration_error(self, graph: "FlowNetwork") -> str | None:
        """Why this wiring is not the configuration the model covers, or None.

        Florschuetz et al. (1981) is ONE configuration: spent air leaving down
        the channel, the crossflow generated by the rows themselves. A
        hand-wired network can describe others that look alike but are not:

        * spent air leaving through the target or between the jets (#468):
          the row's node must drain through one ImpingementCrossflowElement;
        * an external, bypass crossflow (#467): every stream feeding the
          crossflow chain upstream must be another row's jets;
        * two plates merging at one node would read each other's jets as
          crossflow.
        """
        out = graph.get_downstream_elements(self.to_node)
        if len(out) != 1 or not isinstance(out[0], ImpingementCrossflowElement):
            return (
                f"its node {self.to_node!r} must drain through exactly one "
                "ImpingementCrossflowElement. Spent air leaving through the "
                "target (impingement-effusion/film) or between the jets is a "
                "different configuration with no model yet (#468)."
            )
        node, seen = self.to_node, set()
        while node not in seen:
            seen.add(node)
            inflow = [e for e in graph.get_upstream_elements(node) if e.id != self.id]
            plates = [e for e in inflow if isinstance(e, ImpingementPlateElement)]
            segments = [e for e in inflow if isinstance(e, ImpingementCrossflowElement)]
            if len(plates) > (0 if node == self.to_node else 1):
                return (
                    f"node {node!r} has more than one plate merging; each row "
                    "needs its own crossflow node."
                )
            other = [e for e in inflow if e not in plates and e not in segments]
            if other:
                return (
                    f"crossflow node {node!r} is fed by {other[0].id!r}, which is "
                    "not a jet row: an external or bypass crossflow is a "
                    "different configuration with no model yet (#467)."
                )
            if not segments:
                break
            node = segments[0].from_node
        return None

    def validate(self) -> None:
        super().validate()
        if self.correlation in ("ReaderHarrisGallagher", "Stolz", "Miller"):
            raise ValueError(
                f"ImpingementPlateElement {self.id!r}: {self.correlation!r} is a "
                "metering correlation for a plate in a pipe; a jet plate has no "
                "pipe and no beta. Use 'fixed', 'IdelchikThick', 'Lichtarowicz' "
                "or 'McGreehanSchotsch'."
            )
        if getattr(self, "_topology_error", None):
            raise ValueError(
                f"ImpingementPlateElement {self.id!r}: covers jet rows whose spent "
                "air leaves down the crossflow channel (Florschuetz et al. "
                f"1981) only; {self._topology_error}"
            )

    def network_flow_inputs(self) -> list[tuple[str, str]]:
        """(element id, node id) of every flow ``htc_and_T`` reads besides its
        own: the crossflow sources, each read at ``to_node``."""
        return [(eid, self.to_node) for eid in self._crossflow_sources]

    def residuals(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """The orifice flow. ``flows`` is accepted and unused: the hole Cd
        does not depend on the crossflow (see the class docstring)."""
        return super().residuals(state_in, state_out)

    def htc_and_T(
        self, state: NetworkMixtureState, flows: dict[str, float] | None = None
    ) -> _ImpingementPlateResult:
        """Target-side h from this row's jets and the network's crossflow.

        ``state`` is the supply (``from_node``) state with ``m_dot`` set to
        this element's flow; ``flows`` maps each crossflow source to its flow
        into ``to_node``. Without ``flows`` the row sees no crossflow. The
        physics and both flow derivatives are C++'s ``jet_row_heat_transfer``.
        """
        flows = flows or {}
        cs = cb.complete_state(state.T, state.P, state.X)
        tr = cs.transport
        m_c = float(sum(flows.get(eid, 0.0) for eid in self._crossflow_sources))
        row = cb.jet_row_heat_transfer(
            self.correlation_set, self._row_geometry, float(state.m_dot), m_c, tr.mu, tr.k, tr.Pr
        )
        mult = self.Nu_multiplier
        return _ImpingementPlateResult(
            h=row.h * mult,
            Nu=row.Nu,
            Re=row.Re_j,
            Pr=tr.Pr,
            T_aw=float(state.T),
            Gc_Gj=row.Gc_Gj,
            extrapolated=bool(row.extrapolated),
            dh_dmdot=row.dh_dm_jet * mult,
            dh_dsources=dict.fromkeys(self._crossflow_sources, row.dh_dm_crossflow * mult),
        )

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float]:
        base = super().diagnostics(state_in, state_out)
        if not base:
            return base
        r = self.htc_and_T(state_in, flows=flows)
        out = {
            "n_holes": float(self.n_holes),
            "Re_j": float(r.Re),
            "Nu": float(r.Nu),
            "htc": float(r.h),
            "T_aw": float(r.T_aw),
            "Gc_Gj": float(r.Gc_Gj),
            "surface_extrapolated": float(r.extrapolated),
        }
        if self.row is not None:
            # Florschuetz's uniform-supply Eq. 8 at this row, for comparison.
            cd = self._effective_Cd(state_in, state_out)
            out["Gc_Gj_closed_form"] = float(
                cb.crossflow_to_jet_ratio_at_row(self.yn_d, self.z_d, cd, self.row)
            )
        return {**base, **out}


class EffectiveAreaConnectionElement(OrificeElement):
    """
    An orifice element with user-specified effective area (Cd * A product).

    This element calculates orifice flow using the incompressible flow equation:
        m_dot = A_eff * sqrt(2 * rho * dP)

    where A_eff is the effective area (the product Cd * A).

    **Incompressible Formulation**: The solver uses `cb.orifice_mdot_and_jacobian()`,
    which applies the incompressible Bernoulli equation with regularization for numerical
    stability. Density is evaluated at upstream conditions. For compressible flows with
    significant Mach number effects, use a nozzle element instead (future feature).

    **Effective Area Interpretation**: The effective area represents the combined effect
    of discharge coefficient and geometric area. For example:
    - A_eff = 0.01 m^2 could represent Cd=1.0 * A=0.01 m^2
    - Or equivalently: Cd=0.8 * A=0.0125 m^2
    - The product (Cd * A) is what matters for flow calculation

    Unlike LosslessConnectionElement, this element produces pressure drop proportional to flow rate.
    Use this when you need a simple area-based flow restriction without detailed geometry modeling.

    The effective area is user-specified and already accounts for all geometric effects
    (entrance/exit losses, contraction, vena contracta, etc.), so upstream/downstream
    geometry discovery is not performed.
    """

    def __init__(self, id: str, from_node: str, to_node: str, effective_area: float) -> None:
        """
        Initialize with effective area (Cd * A product).

        Args:
            id: Element identifier
            from_node: Upstream node ID
            to_node: Downstream node ID
            effective_area: Effective area for flow calculation (m^2).
                           This is the product Cd * A, where Cd is the discharge coefficient
                           and A is the geometric area. Stored internally as self.area.

        Example:
            >>> conn = EffectiveAreaConnectionElement("conn1", "inlet", "outlet", 0.01)
            >>> conn.area  # 0.01 m^2 (effective area = Cd * A)
            >>> conn.Cd    # 1.0 (normalized discharge coefficient)

        Note:
            The effective area already includes all loss coefficients. For example, if you
            have a sharp-edged orifice with geometric area 0.0125 m^2 and Cd=0.8, you would
            specify effective_area=0.01 m^2 (0.8 * 0.0125).
        """

        diameter = math.sqrt(4.0 * effective_area / math.pi)
        super().__init__(id, from_node, to_node, Cd=1.0, diameter=diameter, correlation="fixed")

        # Force self.area to be precisely the effective area
        # (to remove any floating point math.sqrt / pi precision issues)
        self.area = effective_area

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """
        Skip upstream/downstream geometry discovery.

        Unlike OrificeElement, we do not discover upstream/downstream diameters because
        the effective area is user-specified and already accounts for all geometric effects.
        No Beta ratio correction is needed or applied.
        """
        # Intentionally empty - effective area is pre-computed by user
        pass

    def n_equations(self) -> int:
        """Return the number of equations (1: mass flow balance)."""
        return 1


class DiameterDischargeCoefficientConnectionElement(OrificeElement):
    """
    An orifice element with user-specified physical diameter and discharge coefficient or loss coefficient.

    This element calculates orifice flow using the incompressible flow equation:
        m_dot = A * Cd * sqrt(2 * rho * dP)

    where A is the computed physical area and Cd is the discharge coefficient.

    **Incompressible Formulation**: Uses the Bernoulli equation with density evaluated
    at upstream conditions. Valid for low Mach number flows (M < 0.3 typically).

    **Parameters**: User provides either Cd (discharge coefficient) or zeta (loss coefficient).
    The other parameter is calculated automatically using the relationship:
        zeta = 1/Cd^2 - 1  or  Cd = 1/sqrt(zeta + 1)

    **Effective Area**: The effective area is calculated as A_eff = A * Cd.

    Unlike LosslessConnectionElement, this element produces pressure drop proportional to flow rate.
    Use this when you have a physical diameter measurement and want to specify loss characteristics
    via discharge coefficient or loss coefficient.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        diameter: float,
        Cd: float | None = None,
        zeta: float | None = None,
    ) -> None:
        """
        Initialize with physical diameter and either Cd or zeta.

        Args:
            id: Element identifier
            from_node: Upstream node ID
            to_node: Downstream node ID
            diameter: Physical geometric diameter (m)
            Cd: Discharge coefficient (optional, 0 < Cd <= 1)
            zeta: Loss coefficient (optional, zeta >= 0)

        Note:
            Exactly one of Cd or zeta must be provided. If both are provided, zeta is ignored.
            If Cd is provided, zeta is calculated as zeta = 1/Cd^2 - 1.
            If zeta is provided, Cd is calculated as Cd = 1/sqrt(zeta + 1).

        Example:
            >>> # Using Cd
            >>> conn = DiameterDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.1, Cd=0.8)
            >>> conn.diameter  # 0.1 m (physical diameter)
            >>> conn.Cd    # 0.8 (discharge coefficient)

            >>> # Using zeta
            >>> conn = DiameterDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.1, zeta=0.5625)
            >>> conn.diameter  # 0.1 m (physical diameter)
            >>> conn.Cd    # 0.8 (calculated from zeta)

        Raises:
            ValueError: If neither Cd nor zeta is provided, or if invalid values are given.
        """
        if Cd is not None and zeta is not None:
            # If both are provided, use Cd and ignore zeta
            zeta = None
        elif Cd is None and zeta is None:
            raise ValueError("Either Cd or zeta must be provided")

        if Cd is not None:
            if not (0 < Cd <= 1):
                raise ValueError("Cd must be in range (0, 1]")
        else:  # zeta is provided
            if zeta < 0:
                raise ValueError("zeta must be >= 0")
            # Calculate Cd from zeta: Cd = 1/sqrt(zeta + 1)
            Cd = 1.0 / (zeta + 1.0) ** 0.5

        # Pass physical diameter and Cd directly to parent
        super().__init__(id, from_node, to_node, Cd=Cd, diameter=diameter, correlation="fixed")

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """
        Skip upstream/downstream geometry discovery.

        The physical area and discharge coefficient are user-specified and already account
        for all geometric effects. No Beta ratio correction is needed or applied.
        """
        # Intentionally empty - area and Cd are pre-computed by user
        pass

    def n_equations(self) -> int:
        """Return the number of equations (1: mass flow balance)."""
        return 1


class LosslessConnectionElement(NetworkElement):
    """
    An ideal connection with no friction or momentum loss.
    Inherently preserves total pressure: Pt_in = Pt_out.
    Does not require geometric parameters.
    """

    def __init__(self, id: str, from_node: str, to_node: str):
        super().__init__(id, from_node, to_node)
        self.area: float | None = None
        self.diameter: float | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass  # Zero-loss elements do not need upstream geometry

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Residual: Pt_in - Pt_out = 0
        res = [state_in.Pt - state_out.Pt]
        jac = {0: {f"{self.from_node}.Pt": 1.0, f"{self.to_node}.Pt": -1.0}}
        return res, jac

    def get_spatial_profile(self, n_steps: int = 100):
        """
        Returns a mock spatial profile array representing the zero-loss state.
        This provides structural symmetry with realistic pipes (e.g. Fanno/Rough)
        for GUI plotters.
        """

        # Simplified placeholder for test completion validation
        # Normally would utilize the true state nodes
        profile = []
        for _ in range(n_steps):
            st = cb.IncompressibleStation()
            st.x = 0.0
            st.P = 100000.0
            st.T = 300.0
            st.rho = 1.2
            st.v = 1.0
            st.M = 0.0
            st.h = 300000.0
            profile.append(st)
        return profile

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        if state_in.P <= 0 or state_in.T <= 0:
            return {}

        cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        ref = _element_reference_block(
            cs_in, m_dot=state_in.m_dot, area=self.area or 0.0, Dh=self.diameter, location="inlet"
        )
        v_in = ref["velocity"]
        mach_in = v_in / cs_in.thermo.a if cs_in.thermo.a > 0 and v_in > 0 else 0.0

        return {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_in),
            **ref,
        }

    def n_equations(self) -> int:
        return 1


# ---------------------------------------------------------------------------
# Discrete pressure loss element: combustor-theta aware or plain cold-flow.
# ---------------------------------------------------------------------------


def _correlation_is_linear_theta(correlation: Callable) -> bool:
    """True if the correlation's value depends on ``ctx.theta``.

    Detected by invoking the correlation twice with different theta values on
    a throwaway context and checking whether the returned xi differs, or by
    checking the ``dxi_dtheta`` return.  The callable is also considered
    theta-aware if its returned ``dxi_dtheta`` is nonzero.
    """
    try:
        ctx_zero = SimpleNamespace(
            theta=0.0,
            T_in=300.0,
            P_in=101325.0,
            X_in=[1.0],
            T_ad=300.0,
            phi=0.0,
            mdot_total=1.0,
            mdot_fuel=0.0,
            mdot_air=0.0,
            X_products=[],
            Y_products=[],
        )
        xi0, dxi_dtheta = correlation(ctx_zero)
    except Exception:
        return False
    if dxi_dtheta != 0.0:
        return True
    ctx_one = SimpleNamespace(**vars(ctx_zero))
    ctx_one.theta = 1.0
    try:
        xi1, _ = correlation(ctx_one)
    except Exception:
        return False
    return abs(xi1 - xi0) > 1e-12


class PressureLossElement(NetworkElement):
    """
    Discrete pressure-loss edge.  Enforces ``Pt_out = Pt_in * (1 - xi)``
    where ``xi`` is supplied by a user correlation callable.

    Four standard correlations are provided in :mod:`combaero.network.pressure_loss`:

    - :class:`ConstantFractionLoss` -- fixed fractional total-pressure drop.
    - :class:`ConstantHeadLoss` -- Euler loss coefficient (dynamic head).
    - :class:`LinearThetaFractionLoss` -- ``xi = k * theta + xi0``.
    - :class:`LinearThetaHeadLoss` -- ``xi = (k * theta + zeta0) * q / P_in``.

    Theta-aware correlations need a reference node with ``has_theta = True``
    (e.g. :class:`CombustorNode`).  Discovery precedence:

    1. ``theta_source=<node_id>`` constructor kwarg (explicit override).
    2. ``to_node`` endpoint if it has ``has_theta = True``.
    3. ``from_node`` endpoint if it has ``has_theta = True``.
    4. Cold-flow fallback: ``theta = 0``.  If the correlation is linear-theta
       a :class:`UserWarning` is emitted at resolution time.

    Placing the element between two ``has_theta = True`` nodes without an
    explicit ``theta_source`` raises ``ValueError`` during
    :meth:`FlowNetwork.validate`.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        correlation: Callable[[Any], tuple[float, float]],
        theta_source: str | None = None,
        area: float | None = None,
        surface: ConvectiveSurface | None = None,
    ):
        super().__init__(id, from_node, to_node)
        self.correlation = correlation
        self.theta_source = theta_source
        self.area = area
        #: Where self.area came from once resolved (issue #262).
        self._area_source = "user" if area is not None else ""
        self.surface = surface or ConvectiveSurface()
        # Resolved during resolve_topology.  Empty string sentinel disallowed
        # so callers can cleanly test ``if self._theta_source_resolved is not None``.
        self._theta_source_resolved: str | None = None
        self._both_endpoints_theta: bool = False

    @property
    def has_convective_surface(self) -> bool:
        return self.surface.area > 0.0

    def htc_and_T(self, state: NetworkMixtureState):
        """Convective HTC on the discrete loss element.

        Treats the element as a short duct of length ``L = Dh`` (L/D = 1),
        where ``Dh = sqrt(4*A/pi)`` is derived from the flow area. Returns
        ``None`` if either the convective surface or the flow area is missing.
        """

        if self.surface.area == 0.0 or not self.area or self.area <= 0.0:
            return None

        rho, _ = _safe_rho(state.density())
        velocity = abs(state.m_dot) / (rho * self.area) if self.area > 0 else 0.0
        diameter = math.sqrt(4.0 * self.area / math.pi)
        length = diameter  # Short duct approximation: L/D = 1

        return self.surface.htc_and_T(
            T=state.T,
            P=state.P,
            X=state.X,
            velocity=velocity,
            diameter=diameter,
            length=length,
            T_hot=math.nan,
            flow_area=self.area,
        )

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Stash a graph reference so residuals/diagnostics can read the theta
        # source node's burned/unburned temperatures without a re-pass.
        self._graph_ref = graph

        # 0. Area inheritance: when no user-specified area, inherit from the
        #    nearest upstream ChannelElement diameter.  This makes discrete loss
        #    work after any element, not just combustor/plenum nodes.
        if self.area is None:
            for elem in graph.get_upstream_elements(self.from_node):
                if isinstance(elem, ChannelElement) and elem.diameter is not None:
                    self.area = math.pi / 4.0 * elem.diameter**2
                    break
            if self.area is None:
                # Widen to any area-bearing neighbour (a combustor or momentum
                # chamber upstream is the common case) before defaulting.
                area, source = _infer_flow_area(graph, [self.from_node, self.to_node], exclude=self)
                if area is not None:
                    self.area = area
                    self._area_source = source
            if self.area is None:
                # A head-loss correlation reads this area as its velocity
                # reference, so a wrong value rescales the loss rather than
                # failing. Warn only when a correlation will actually consume
                # it -- for a fraction-based correlation nothing reads it and
                # a nominal value is harmless (issue #262).
                self.area = DEFAULT_CHAMBER_AREA
                self._area_source = "default"
                if hasattr(self.correlation, "area"):
                    warnings.warn(
                        f"PressureLossElement '{self.id}': no flow area could be "
                        f"inferred from the topology, so the nominal "
                        f"{DEFAULT_CHAMBER_AREA} m^2 default is used as the "
                        f"velocity reference for "
                        f"{type(self.correlation).__name__}. The head loss scales "
                        f"as 1/area^2 -- set 'area' explicitly to make it "
                        f"meaningful.",
                        UserWarning,
                        stacklevel=2,
                    )
            # Propagate updated area into head-loss correlations and convective surface.
            if hasattr(self.correlation, "area"):
                self.correlation.area = self.area
            if self.surface.area == 0.0:
                self.surface.area = self.area

        # 1. Explicit override wins.
        if self.theta_source is not None:
            src = graph.nodes.get(self.theta_source)
            if src is None:
                raise ValueError(
                    f"PressureLossElement '{self.id}': theta_source "
                    f"'{self.theta_source}' is not in the network."
                )
            if not getattr(src, "has_theta", False):
                raise ValueError(
                    f"PressureLossElement '{self.id}': theta_source "
                    f"'{self.theta_source}' is a '{type(src).__name__}' which "
                    f"does not provide theta (has_theta=False)."
                )
            self._theta_source_resolved = self.theta_source
            return

        to_has_theta = getattr(graph.nodes[self.to_node], "has_theta", False)
        from_has_theta = getattr(graph.nodes[self.from_node], "has_theta", False)

        if to_has_theta and from_has_theta:
            # Ambiguous: raised at validate() for a cleaner user-facing error.
            self._both_endpoints_theta = True
            return
        if to_has_theta:
            self._theta_source_resolved = self.to_node
            return
        if from_has_theta:
            self._theta_source_resolved = self.from_node
            return

        # 4. Cold-flow fallback.  Warn if the correlation depends on theta.
        if _correlation_is_linear_theta(self.correlation):
            warnings.warn(
                f"PressureLossElement '{self.id}' uses a linear-theta "
                f"correlation ({type(self.correlation).__name__}) but neither "
                f"endpoint ('{self.from_node}' upstream, '{self.to_node}' "
                f"downstream) has has_theta=True. The k-term will have no "
                f"effect; xi falls back to cold-flow (theta=0).",
                UserWarning,
                stacklevel=2,
            )

    def _build_ctx(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        graph: "FlowNetwork | None",
    ) -> SimpleNamespace:
        """Build a duck-typed PressureLossContext from the current states."""
        # Theta sourcing: read the reference node's burned/unburned temperatures.
        if self._theta_source_resolved is not None and graph is not None:
            src = graph.nodes[self._theta_source_resolved]
            T_unburned = float(getattr(src, "_T_unburned", state_in.T))
            T_burned = float(getattr(src, "_T_burned", state_in.T))
            theta = (T_burned / T_unburned - 1.0) if T_unburned > 0 else 0.0
            T_ctx = T_unburned
            T_ad = T_burned
        else:
            theta = 0.0
            T_ctx = float(state_in.T)
            T_ad = float(state_in.T)

        return SimpleNamespace(
            theta=theta,
            T_in=T_ctx,
            P_in=float(state_in.P),
            X_in=list(state_in.X),
            T_ad=T_ad,
            phi=0.0,
            mdot_total=float(state_in.m_dot),
            mdot_fuel=0.0,
            mdot_air=0.0,
            X_products=[],
            Y_products=[],
        )

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Graph reference is stashed on the element by the solver pre-pass;
        # fall back to ctx without theta sourcing if unavailable.
        graph = getattr(self, "_graph_ref", None)
        ctx = self._build_ctx(state_in, state_out, graph)
        try:
            xi, dxi_dtheta = self.correlation(ctx)
        except Exception as exc:
            raise RuntimeError(
                f"PressureLossElement '{self.id}' correlation raised: {exc}"
            ) from exc

        # Residual: Pt_out - Pt_in * (1 - xi) = 0.
        res = [state_out.Pt - state_in.Pt * (1.0 - xi)]

        # The Jacobian is differenced (one or two correlation calls per
        # variable and per species present), so it is computed only when the
        # solver reads it (LazyJacobian, #489).
        def jacobian() -> dict[int, dict[str, float]]:
            return self._jacobian(state_in, graph, ctx, xi, dxi_dtheta)

        return res, LazyJacobian(jacobian)

    def _jacobian(
        self,
        state_in: NetworkMixtureState,
        graph: "FlowNetwork | None",
        ctx: SimpleNamespace,
        xi: float,
        dxi_dtheta: float,
    ) -> dict[int, dict[str, float]]:
        jac: dict[int, dict[str, float]] = {
            0: {
                f"{self.from_node}.Pt": -(1.0 - xi),
                f"{self.to_node}.Pt": 1.0,
            }
        }

        # ----- Analytical theta Jacobian (linear-theta correlations) -----
        # d_res/d_T_burned = -Pt_in * dxi_dtheta / T_in.
        if self._theta_source_resolved is not None and dxi_dtheta != 0.0 and ctx.T_in > 0:
            jac[0][f"{self._theta_source_resolved}.T"] = state_in.Pt * dxi_dtheta / ctx.T_in

        # ----- FD-based Jacobian for correlation dependencies not captured
        # analytically (e.g. zeta-based xi depends on mdot, P_in, T_in via rho).
        # One extra correlation call per variable; cheap and robust.
        def _dxi_dvar(attr: str, val: float, step_scale: float) -> float:
            eps = step_scale if val == 0.0 else max(abs(val) * 1e-6, step_scale)
            saved = getattr(ctx, attr)
            try:
                setattr(ctx, attr, saved + eps)
                xi_p, _ = self.correlation(ctx)
                setattr(ctx, attr, saved - eps)
                xi_m, _ = self.correlation(ctx)
            finally:
                setattr(ctx, attr, saved)
            return (xi_p - xi_m) / (2.0 * eps)

        # Only bother if xi actually varies with the velocity/density variables
        # (skip for constant fraction to avoid unneeded correlation calls).
        dxi_dmdot = _dxi_dvar("mdot_total", ctx.mdot_total, 1e-6)
        if dxi_dmdot != 0.0:
            jac[0][f"{self.id}.m_dot"] = state_in.Pt * dxi_dmdot

        dxi_dP_in = _dxi_dvar("P_in", ctx.P_in, 1.0)
        if dxi_dP_in != 0.0:
            # P_in = state_in.P.  Static P is driven by the node's P unknown.
            jac[0][f"{self.from_node}.P"] = state_in.Pt * dxi_dP_in

        if self._theta_source_resolved is None:
            # When no theta source, ctx.T_in = state_in.T (a derived state).
            # Emit a relay-routable T Jacobian on the upstream node.
            dxi_dT_in = _dxi_dvar("T_in", ctx.T_in, 0.01)
            if dxi_dT_in != 0.0:
                jac[0][f"{self.from_node}.T"] = state_in.Pt * dxi_dT_in
        elif graph is not None and ctx.T_in > 0:
            self._add_unburned_jacobian(jac[0], graph, ctx, state_in.Pt, dxi_dtheta, _dxi_dvar)

        # Composition of the inlet state (a head loss's density): per species
        # present, other mass fractions fixed, relayed as "{node}.Y[k]". A
        # combustor's burned composition moves with its fuel flow; this column
        # was missing (#481 C1).
        Y_in = list(state_in.Y)
        X_saved = ctx.X_in
        try:
            for k, y in enumerate(Y_in):
                if not y > 0.0:
                    continue
                h = min(1e-6, 0.5 * y)
                xi_pm = []
                for sgn in (1.0, -1.0):
                    Yk = list(Y_in)
                    Yk[k] += sgn * h
                    ctx.X_in = list(cb.mass_to_mole(Yk))
                    xi_pm.append(self.correlation(ctx)[0])
                d = (xi_pm[0] - xi_pm[1]) / (2.0 * h)
                if d != 0.0:
                    jac[0][f"{self.from_node}.Y[{k}]"] = state_in.Pt * d
        finally:
            ctx.X_in = X_saved

        return jac

    def _add_unburned_jacobian(
        self,
        row: dict[str, float],
        graph: "FlowNetwork",
        ctx: SimpleNamespace,
        Pt_in: float,
        dxi_dtheta: float,
        dxi_dvar: Callable[[str, float, float], float],
    ) -> None:
        """d(res)/d(unknowns) through the source's UNBURNED temperature (#481 C1).

        T_u enters twice: theta = T_b/T_u - 1 (d theta/dT_u = -T_b/T_u^2) and
        as the correlation's reference temperature ctx.T_in (a head loss's
        density). T_u = sum(m_i Tt_i)/M over the source's inflows, so
        dT_u = sum((Tt_i - T_u)/M dm_i + m_i/M dTt_i): the flows' unknowns
        directly, each Tt_i through its node's relayed "{src}.T". It was
        missing: d/d(fuel flow) read 0 against 95 by differences.
        """
        src = graph.nodes[self._theta_source_resolved]
        streams = getattr(src, "_unburned_streams", None)
        M = sum(m for m, _, _, _ in streams) if streams else 0.0
        if not streams or M <= 0.0:
            return
        T_u, T_b = ctx.T_in, ctx.T_ad
        dres_dTu = Pt_in * (dxi_dtheta * (-T_b / (T_u * T_u)) + dxi_dvar("T_in", T_u, 0.01))
        if dres_dTu == 0.0:
            return
        for m, Tt, names, src_node in streams:
            for var, coeff in names.items():
                row[var] = row.get(var, 0.0) + dres_dTu * (Tt - T_u) / M * coeff
            if src_node is not None:
                key = f"{src_node}.T"
                row[key] = row.get(key, 0.0) + dres_dTu * m / M

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        if state_in.P <= 0 or state_in.T <= 0 or state_out.P <= 0 or state_out.T <= 0:
            return {}

        graph = getattr(self, "_graph_ref", None)
        ctx = self._build_ctx(state_in, state_out, graph)
        try:
            xi, _ = self.correlation(ctx)
        except Exception:
            xi = 0.0

        area = getattr(self, "area", None) or getattr(self.correlation, "area", None)
        Dh = math.sqrt(4.0 * area / math.pi) if area and area > 0 else None

        cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        ref = _element_reference_block(
            cs_in, m_dot=state_in.m_dot, area=area or 0.0, Dh=Dh, location="inlet"
        )
        v_in = ref["velocity"]
        mach_in = v_in / cs_in.thermo.a if cs_in.thermo.a > 0 and v_in > 0 else 0.0

        cs_out = cb.complete_state(state_out.T, state_out.P, state_out.X)
        rho_out = cs_out.thermo.rho
        v_out = (
            (abs(state_out.m_dot)) / ((rho_out) * (area))
            if area and area > 0 and rho_out > 0
            else 0.0
        )
        mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 and v_out > 0 else 0.0

        diag: dict[str, float] = {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            "xi": float(xi),
            "theta": float(ctx.theta),
        }

        if self.has_convective_surface:
            h_res = self.htc_and_T(state_in)
            if h_res is not None:
                diag["Nu"] = float(h_res.Nu)
                diag["htc"] = float(h_res.h)
                diag["T_aw"] = float(h_res.T_aw)
                diag["f"] = float(h_res.f)

        return diag


class ChannelElement(NetworkElement):
    """Duct with frictional pressure drop and optional convective surface.

    Couples STAGNATION pressures: Pt_up - Pt_down = dP_friction. The
    static-versus-stagnation distinction lives in the node constitutive
    relation, not here (see the note in ``residuals``).

    Pressure drop
        ``regime="incompressible"`` uses Darcy-Weisbach with a friction factor
        selected by ``friction_model`` from Haaland, Colebrook, Serghides or
        Petukhov; ``regime="compressible"`` uses Fanno flow. Both are computed
        in C++ -- ``friction.h`` owns the correlations and
        ``solver_interface.h`` the (f, J) entry points. Nothing about the
        friction factor is restated in this file.

    Convective surfaces
        A ``ConvectiveSurface`` adds heat transfer. Ribbed and pin-fin surfaces
        also OWN the drop: their correlations return the channel's friction
        factor (ribs: the four-sided Darcy value; pins: Metzger's per-row
        ``dP / (2 rho Vmax^2 N)``), and the element uses it directly rather
        than as a multiplier on pipe friction -- see ``_ribbed_residuals`` and
        ``_pin_fin_residuals``. Impingement surfaces leave the drop to the
        smooth path (model the jet plate with an ``OrificeElement``).
        Provenance lives with each correlation set: ``rib_correlation.h``,
        ``pin_fin_correlation.h``, ``impingement_correlation.h``.

    Known gaps
        The array correlations expose no dP sensitivity to static pressure or
        composition, so those Jacobian columns are absent rather than
        approximated; ``test_channel_element_jacobian_fd.py`` pins this as a
        strict xfail carrying the magnitude.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        length: float,
        diameter: float | None = None,
        Dh: float | None = None,
        roughness: float = 1e-5,
        regime: CompressibilityLiteral = "incompressible",
        friction_model: FrictionModelLiteral = "haaland",
        htc_model: HeatTransferModelLiteral = "none",
        t_hot: float | None = None,
        surface: ConvectiveSurface | None = None,
    ):
        super().__init__(id, from_node, to_node)
        self.length = length
        self.diameter: float | None = diameter
        # Dh: user-specified hydraulic diameter (for non-circular ducts).
        # None until resolved: defaults to diameter (circular assumption).
        self.Dh: float | None = Dh if Dh is not None else diameter
        self.roughness = roughness
        self.area: float | None = math.pi * (diameter / 2) ** 2 if diameter is not None else None

        self.regime = regime
        self.friction_model = friction_model
        self.htc_model = htc_model
        self.t_hot = t_hot
        self.surface = surface or ConvectiveSurface()

    @property
    def has_convective_surface(self) -> bool:
        return True

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def _reynolds(self, state: NetworkMixtureState) -> float:
        """Reynolds number on the hydraulic diameter, for surface correlations.

        Density cancels: rho * v is m_dot / area, so Re = m_dot * Dh / (A * mu)
        and no density guard is needed. This replaced a block that substituted
        air at STP (rho = 1.2, mu = 1.8e-5) whenever a state looked unphysical
        -- fabricated properties that a mid-Newton iterate would silently pick
        up, in the path where a discontinuity costs most.
        """
        area = self.area or 0.0
        if area <= 0.0:
            return 0.0
        mu = cb.complete_state(state.T, state.Pt, state.X).transport.mu
        if mu <= 0.0:
            return 0.0
        return abs(state.m_dot) * (self.Dh or self.diameter or 1.0) / (area * mu)

    def _ribbed_residuals(self, state_in: NetworkMixtureState, state_out: NetworkMixtureState):
        """Ribbed channel: the correlation owns the friction factor.

        No multiplier on pipe friction. The `f` the roughness function is
        defined against is already the equivalent four-sided channel value, so
        the correlation returns what the channel needs and the element uses it
        directly. Routing it through a locally restated smooth friction factor
        is what #331 removed, after the round-trip moved a drop by -24.7%.

        The Jacobian is simpler than it looks. `R` carries no `e+` term in this
        correlation family, so **f does not depend on Reynolds number** and
        therefore not on mass flow: `df/d(mdot) = 0`. The drop is
        `f (L/D) rho v |v| / 2` with `v = mdot / (rho A)`, so

            dP          = f (L/D) mdot |mdot| / (2 rho A^2)
            d(dP)/dmdot = f (L/D) |mdot| / (rho A^2)

        which is even in `mdot` and so does not flip sign at zero -- the drop
        itself carries the sign. That is the odd-quantity case: magnitude from
        the correlation, direction from the flow.

        `dP` also depends on upstream temperature and static pressure through
        `rho(T, P)` -- `d(dP)/d(rho) = -dP/rho`, chained through
        `density_and_jacobians`' analytic `d(rho)/dT`, `d(rho)/dP` and
        `_safe_rho`'s own floor derivative. This was a real, measured gap
        (issue #378): the analytic `A.T` column was exactly 0.0 against a
        finite-difference value of roughly -486 for a representative case,
        because `_safe_rho(state_in.density())` computed `rho` without ever
        differentiating it. Composition (`Y`) sensitivity through molecular
        weight remains an accepted gap -- no correlation-level Jacobian for
        it exists yet, and no other channel-friction path in this file
        exposes one either.
        """
        m_dot = state_in.m_dot
        model = self.surface.model
        f_mult = self.surface.f_multiplier

        rho_raw, drho_raw_dT, drho_raw_dP = _solver_tools.density_and_jacobians(
            state_in.T, state_in.P, state_in.X
        )
        rho, drho_draw = _safe_rho(rho_raw)
        drho_dT = drho_draw * drho_raw_dT
        drho_dP = drho_draw * drho_raw_dP
        area = self.area or 0.0
        if area <= 0.0:
            return [state_in.Pt - state_out.Pt], {
                0: {
                    f"{self.from_node}.Pt": 1.0,
                    f"{self.to_node}.Pt": -1.0,
                }
            }

        dh = self.Dh or self.diameter or 1.0
        mu = cb.complete_state(state_in.T, state_in.P, state_in.X).transport.mu
        Re = abs(m_dot) * dh / (area * mu) if mu > 0.0 else 0.0
        Re = math.copysign(Re, m_dot)

        rib = cb.evaluate_rib(
            model.correlation_set or cb.han_1988_orthogonal(),
            cb.RibGeometry(
                e_D=model.e_D,
                p_e=model.p_e,
                W_H=model.W_H,
                alpha_deg=model.alpha_deg,
            ),
            Re,
        )
        f = rib.f * f_mult

        coeff = f * (self.length / dh) / (2.0 * rho * area * area)
        dP = coeff * m_dot * abs(m_dot)
        d_dP_d_mdot = 2.0 * coeff * abs(m_dot)
        # dP = K * mdot|mdot| / rho for K independent of rho, so
        # d(dP)/d(rho) = -dP/rho.
        d_dP_drho = -dP / rho if rho > 0.0 else 0.0
        d_dP_dT = d_dP_drho * drho_dT
        d_dP_dP_static = d_dP_drho * drho_dP

        res = [state_in.Pt - state_out.Pt - dP]
        jac = {
            0: {
                f"{self.id}.m_dot": -d_dP_d_mdot,
                f"{self.from_node}.Pt": 1.0,
                f"{self.to_node}.Pt": -1.0,
                f"{self.from_node}.T": -d_dP_dT,
                f"{self.from_node}.P": -d_dP_dP_static,
            }
        }
        return res, jac

    def _pin_fin_residuals(self, state_in: NetworkMixtureState, state_out: NetworkMixtureState):
        """Pin-fin array: the correlation owns the drop, per row.

            dP = 2 rho Vmax^2 N f(Re_D),   Vmax = mdot / (rho A_min)
               = 2 N f mdot |mdot| / (rho A_min^2)
            Re_D = mdot D / (mu A_min)                          (signed)

        with ``A_min = A_channel * (A_min/A_frontal)`` from the array geometry
        and ``f`` the canonical per-row factor of ``pin_fin_correlation.h``.
        Unlike ribs, pin friction DOES depend on Re, so

            d(dP)/dmdot = 2N/(rho A_min^2) (2 |mdot| f + mdot |mdot| f' D/(mu A_min))

        which is even in mdot, as the ribbed path's. Upstream T and P enter
        through rho (``-dP/rho``) and through mu in Re
        (``dP/df * f' * (-Re/mu) * dmu``), both from analytic or tabulated
        property derivatives (``density_and_jacobians``,
        ``viscosity_and_jacobians``). Composition sensitivity is the same
        accepted gap as every other channel-friction path here.
        """
        m_dot = state_in.m_dot
        model = self.surface.model
        f_mult = self.surface.f_multiplier

        rho_raw, drho_raw_dT, drho_raw_dP = _solver_tools.density_and_jacobians(
            state_in.T, state_in.P, state_in.X
        )
        rho, drho_draw = _safe_rho(rho_raw)
        drho_dT = drho_draw * drho_raw_dT
        drho_dP = drho_draw * drho_raw_dP
        area = self.area or 0.0
        if area <= 0.0:
            return [state_in.Pt - state_out.Pt], {
                0: {
                    f"{self.from_node}.Pt": 1.0,
                    f"{self.to_node}.Pt": -1.0,
                }
            }

        mu, dmu_dT, dmu_dP = _solver_tools.viscosity_and_jacobians(
            state_in.T, state_in.P, state_in.X
        )
        D = model.pin_diameter
        a_min = area * cb.pin_fin_amin_over_afrontal(model.geometry())
        Re = m_dot * D / (a_min * mu) if mu > 0.0 else 0.0

        _, _, f, df_dRe, _ = self.surface._pin_fin_terms(Re, 0.7)
        f *= f_mult
        df_dRe *= f_mult

        coeff = 2.0 * model.N_rows / (rho * a_min * a_min)
        dP = coeff * f * m_dot * abs(m_dot)
        dRe_dmdot = D / (a_min * mu) if mu > 0.0 else 0.0
        d_dP_d_mdot = coeff * (2.0 * abs(m_dot) * f + m_dot * abs(m_dot) * df_dRe * dRe_dmdot)

        # Through rho (dP ~ 1/rho at fixed f) and through mu in Re.
        d_dP_drho = -dP / rho if rho > 0.0 else 0.0
        d_dP_dmu = coeff * m_dot * abs(m_dot) * df_dRe * (-Re / mu) if mu > 0.0 else 0.0
        d_dP_dT = d_dP_drho * drho_dT + d_dP_dmu * dmu_dT
        d_dP_dP_static = d_dP_drho * drho_dP + d_dP_dmu * dmu_dP

        res = [state_in.Pt - state_out.Pt - dP]
        jac = {
            0: {
                f"{self.id}.m_dot": -d_dP_d_mdot,
                f"{self.from_node}.Pt": 1.0,
                f"{self.to_node}.Pt": -1.0,
                f"{self.from_node}.T": -d_dP_dT,
                f"{self.from_node}.P": -d_dP_dP_static,
            }
        }
        return res, jac

    @property
    def residual_scale_kind(self) -> str:
        """Units of this element's residual row (RESIDUAL-ROW-SCALING, see
        NetworkSolver._build_residual_scales): "mdot" for the compressible
        mass-flow form m - m_calc (#481), "p" for the pressure-drop forms
        (incompressible, ribbed, pin-fin). Must change with the residual."""
        if self.regime == "compressible" and not (
            self.surface and isinstance(self.surface.model, (RibbedModel, PinFinModel))
        ):
            return "mdot"
        return "p"

    def _exit_head_lost(self, m_dot: float) -> bool:
        """Whether this channel's exit dynamic head is lost where it discharges.

        The exit is the node the flow goes INTO: to_node forward, from_node
        when reversed. Only a PressureBoundary declares a coupling; anything
        else (a plenum, a momentum chamber, a junction port) carries its own
        constitutive relation and receives the stagnation pressure. See #360.
        """
        attr = "_downstream_node" if m_dot >= 0.0 else "_upstream_node"
        node = getattr(self, attr, None)
        if node is None or not hasattr(node, "exit_head_lost"):
            return False
        return bool(node.exit_head_lost(True))

    def _compressible_flow(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        f_mult: float,
        derivatives: bool = True,
    ) -> tuple[float, dict[str, float], Any]:
        """Fanno mass flow m_calc and d(m_calc)/d(node unknowns) (#481).

        The duct is fed isentropically from the stagnation state of the node
        the flow comes FROM and discharges against the other node's
        stagnation pressure, matched to its exit static pressure (head lost)
        or exit stagnation pressure (recovered). The flux rises monotonically
        as that pressure falls and saturates at the choked flux: it exists for
        every state, so there is no infeasible region and no barrier.

        derivatives=False returns m_calc only (empty derivs): the kernel then
        skips its implicit derivatives. With derivatives, d/dY of the feeding
        node is included where its composition moves with the solve: the
        solver sets _composition_directions -- per such node, an orthonormal
        basis of the directions it moves in -- before it reads a Jacobian.
        The kernel differentiates along that basis and the gradient is
        projected back onto per-species keys, exact for every relay column.
        Standalone (unset): per species present, other fractions fixed.
        """
        A = self.area
        forward = state_in.Pt >= state_out.Pt
        src, dst = (state_in, state_out) if forward else (state_out, state_in)
        src_id, dst_id = (
            (self.from_node, self.to_node) if forward else (self.to_node, self.from_node)
        )
        lost = self._exit_head_lost(1.0 if forward else -1.0)
        known = getattr(self, "_composition_directions", None)
        if not derivatives:
            dirs: list[Any] = []
        elif known is None:
            n_sp = len(src.Y)
            dirs = [[float(j == k) for j in range(n_sp)] for k, y in enumerate(src.Y) if y > 0.0]
        else:
            dirs = known.get(src_id, [])
        flow = cb.fanno_channel_flow(
            src.Pt,
            src.Tt,
            src.X,
            dst.Pt,
            not lost,
            self.length,
            self.diameter,
            self.roughness,
            self.friction_model,
            f_mult,
            getattr(self, "_last_M_exit", -1.0),
            [list(v) for v in dirs],
            derivatives,
        )
        # Warm start only: the flow does not depend on it.
        if not flow.choked:
            self._last_M_exit = float(flow.M_exit)
        sign = 1.0 if forward else -1.0
        if not derivatives:
            return sign * A * flow.G, {}, flow
        derivs = {
            f"{src_id}.Pt": sign * A * flow.dG_dPt0,
            f"{src_id}.T": sign * A * flow.dG_dTt0,
            f"{dst_id}.Pt": sign * A * flow.dG_dP_target,
        }
        # Composition of the feeding node: the directional derivatives
        # projected back onto species, g = sum_i (dG/dv_i) v_i, chained by the
        # solver with dY/dx from its mixing relay.
        if dirs:
            g = [0.0] * len(src.Y)
            for d, v in zip(flow.dG_ddir, dirs, strict=True):
                for k, vk in enumerate(v):
                    g[k] += d * float(vk)
            for k, gk in enumerate(g):
                if gk != 0.0:
                    derivs[f"{src_id}.Y[{k}]"] = sign * A * gk
        return sign * A * flow.G, derivs, flow

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> list[float]:

        m_dot = state_in.m_dot

        # User-set empirical correction on f. Not a correlation: it encodes
        # nothing and defaults to 1.0, so it is inert unless a caller reaches
        # for it. Correlation-derived multipliers were removed in 0.7.0; this
        # one is a tuning knob and stays. See issue #339.
        f_mult = self.surface.f_multiplier if self.surface else 1.0

        if self.surface and isinstance(self.surface.model, RibbedModel):
            return self._ribbed_residuals(state_in, state_out)
        if self.surface and isinstance(self.surface.model, PinFinModel):
            return self._pin_fin_residuals(state_in, state_out)

        if self.regime == "compressible":
            # Mass-flow form, like the compressible orifice: m - m_calc, with
            # m_calc the Fanno flow the end states drive (#481). The former
            # drop-given-m form had no physical drop past choke and patched
            # one in with a barrier, which converged choked ducts 5-20% above
            # their choked flow.
            # The residual needs the flow only; its implicit derivatives are
            # computed if and when the Jacobian is read (LazyJacobian, #489).
            m_calc, _, _ = self._compressible_flow(state_in, state_out, f_mult, derivatives=False)

            def jacobian() -> dict[int, dict[str, float]]:
                _, derivs, _ = self._compressible_flow(state_in, state_out, f_mult)
                jac_c: dict[str, float] = {f"{self.id}.m_dot": 1.0}
                for name, d in derivs.items():
                    jac_c[name] = jac_c.get(name, 0.0) - d
                return {0: jac_c}

            return [m_dot - m_calc], LazyJacobian(jacobian)
        else:
            # Use incompressible Darcy-Weisbach formulation. Density
            # reference: upstream static by default, downstream static when
            # the solver's compressible-seed proxy sets
            # _incompressible_p_ref = "outlet" (see OrificeElement).
            _p_ref_outlet = getattr(self, "_incompressible_p_ref", "inlet") == "outlet"
            res_cpp = _solver_tools.channel_residuals_and_jacobian(
                m_dot,
                state_in.Pt,
                state_out.P if _p_ref_outlet else state_in.P,
                state_in.T,
                state_in.Y,
                state_out.P,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                f_mult,
            )

        # Residual of P-node unknowns matching channel drop.
        # C++ returns dP_calc = friction stagnation-pressure loss only (the
        # downstream static pressure is not an input). A channel therefore
        # connects stagnation pressures: Pt_up - Pt_down = dP_friction. The
        # static-vs-stagnation distinction lives entirely in the node
        # constitutive relation (PlenumNode: Pt = P; MomentumChamberNode:
        # Pt = P + 0.5*rho*v^2). Referencing the downstream node's static P
        # instead would let an inline MomentumChamberNode leak its dynamic head
        # q_N = Pt - P as a free pressure gain across the junction. For every
        # plenum- or boundary-terminated channel Pt == P, so this is identical
        # to the previous static-face coupling.
        res = [state_in.Pt - state_out.Pt - res_cpp.dP_calc]

        upstream_id = self.from_node
        downstream_id = self.to_node

        jac = {
            0: {
                f"{self.id}.m_dot": -res_cpp.d_dP_d_mdot,
                f"{upstream_id}.Pt": 1.0,
                f"{upstream_id}.T": -res_cpp.d_dP_dT_up,
                f"{downstream_id}.Pt": -1.0,
            }
        }
        if getattr(self, "_incompressible_p_ref", "inlet") == "outlet":
            # Density evaluated at the downstream static: the friction-loss
            # pressure sensitivity moves onto the downstream node's P.
            jac[0][f"{downstream_id}.P"] = -res_cpp.d_dP_dP_static_up
        else:
            jac[0][f"{upstream_id}.P"] = -res_cpp.d_dP_dP_static_up
        for i, val in enumerate(res_cpp.d_dP_dY_up):
            jac[0][f"{upstream_id}.Y[{i}]"] = -val

        return res, jac

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float | str]:
        cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        ref = _element_reference_block(
            cs_in,
            m_dot=state_in.m_dot,
            area=self.area,
            Dh=self.Dh or self.diameter,
            location="inlet",
        )
        v_in = ref["velocity"]
        re_in = ref["Re"]
        mach_in = v_in / cs_in.thermo.a if cs_in.thermo.a > 0 and v_in > 0 else 0.0

        cs_out = cb.complete_state(state_out.T, state_out.P, state_out.X)
        rho_out = cs_out.thermo.rho
        v_out = (
            (abs(state_out.m_dot)) / ((rho_out) * (self.area))
            if self.area > 0 and rho_out > 0
            else 0.0
        )
        mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 and v_out > 0 else 0.0

        h_res = self.htc_and_T(state_in)
        surface_diag: dict[str, float] = {}
        if h_res is not None:
            Nu, htc, T_aw, f = h_res.Nu, h_res.h, h_res.T_aw, h_res.f
            # The surface correlation's own Re (jet Re_j for impingement, pin
            # Re_D for pin fins) and whether it is outside its validity box.
            surface_diag["Re_surface"] = float(h_res.Re)
            surface_diag["surface_extrapolated"] = float(
                bool(getattr(h_res, "extrapolated", False))
            )
        else:
            Nu, htc, T_aw = 0.0, 0.0, state_in.T
            dh_diag = self.Dh or self.diameter or 1.0
            e_D = self.roughness / dh_diag if dh_diag > 0 else 0.0
            # Apply the same Re floor used by the C++ friction correlations
            # (Re_eff = sqrt(Re^2 + 64^2)) so f stays finite at near-zero flow.
            re_eff = math.sqrt(re_in * re_in + 64.0 * 64.0)
            if re_eff < 2300:
                # Poiseuille, the physical law, not a choice of correlation.
                f = 64.0 / re_eff
            else:
                # Report the correlation the residual actually uses. This
                # branch used to hardcode Haaland (rough) or Petukhov (smooth)
                # regardless of friction_model, so a channel set to petukhov
                # reported a friction factor 7.5% away from the one driving
                # its own pressure drop.
                f = cb._core.friction_and_jacobian(self.friction_model, re_eff, e_D).result[0]

        duct: dict[str, float] = {}
        if self.regime == "compressible":
            # The duct's own Fanno solution: inlet and exit Mach inside the
            # duct (the node states above are the plenum/boundary states), and
            # whether it is choked -- its flow then independent of the back
            # pressure.
            f_mult = self.surface.f_multiplier if self.surface else 1.0
            _, _, flow = self._compressible_flow(state_in, state_out, f_mult, derivatives=False)
            duct = {
                "choked": float(flow.choked),
                "M_in_duct": float(flow.M_in),
                "M_exit_duct": float(flow.M_exit),
            }
        return {
            "m_dot": float(state_in.m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            **duct,
            "Dh": float(self.Dh or self.diameter or 0.0),
            "Nu": float(Nu),
            "htc": float(htc),
            "T_aw": float(T_aw),
            "f": float(f),
            **surface_diag,
        }

    def get_spatial_profile(
        self,
        state_in: NetworkMixtureState,
        n_steps: int = 100,
    ) -> list:
        """
        Compute and return the spatial flow profile array along the channel length.
        Requires the solved inlet NetworkMixtureState.
        """

        rho = cb.density(state_in.T, state_in.P, state_in.X)
        u = cb.channel_velocity(state_in.m_dot, self.diameter, rho)

        if self.regime == "incompressible":
            res = cb.channel_flow_rough(
                state_in.T,
                state_in.P,
                state_in.X,
                u,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                n_steps,
                True,
            )
            return res.profile

        elif self.regime == "compressible":
            # The x-march from the duct's inlet STATIC state, which the Mach
            # solution gives (isentropic from the feeding node's stagnation
            # state). It called fanno_channel with its arguments out of order
            # and raised for every compressible channel (#481).
            G = abs(state_in.m_dot) / self.area
            if G <= 0.0:
                return []
            duct = cb.fanno_duct(
                state_in.Pt,
                state_in.Tt,
                G,
                state_in.X,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                self.surface.f_multiplier if self.surface else 1.0,
            )
            res = cb.fanno_channel_rough(
                duct.inlet.T,
                duct.inlet.P,
                duct.inlet.u,
                self.length,
                self.diameter,
                self.roughness,
                state_in.X,
                self.friction_model,
                self.surface.f_multiplier if self.surface else 1.0,
                n_steps,
                True,
            )
            return res.profile

        return []

    def htc_and_T(self, state: NetworkMixtureState):
        """Compute heat transfer coefficient and adiabatic wall temperature.

        Uses the element's ConvectiveSurface to compute HTC and adiabatic wall temperature.

        Parameters
        ----------
        state : NetworkMixtureState
            Flow state at the element inlet.

        Returns
        -------
        ChannelResult | None
            Full ChannelResult with h, T_aw, and Jacobians (dh_dmdot, dh_dT, etc.),
            or None if surface.area = 0. Access convective area via ``self.surface.area``.
        """
        if self.surface.area == 0.0:
            return None

        rho, _ = _safe_rho(state.density())
        u = abs(state.m_dot) / (rho * self.area) if self.area > 0 else 1.0

        # Use nan for T_hot if not specified (matches C++ default)
        T_hot = self.t_hot if self.t_hot is not None else math.nan

        return self.surface.htc_and_T(
            T=state.T,
            P=state.P,
            X=state.X,
            velocity=u,
            diameter=self.Dh or self.diameter,
            length=self.length,
            T_hot=T_hot,
            flow_area=self.area,
        )

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Cached for the exit-coupling lookup in residuals (#360); done before
        # the early return so it is set even when the diameter is explicit.
        self._downstream_node = graph.nodes.get(self.to_node)
        self._upstream_node = graph.nodes.get(self.from_node)
        if self.diameter is not None:
            self._fill_convective_area()
            return
        # Inherit diameter from the nearest geometry source. The channel's
        # area sets its friction and Mach behaviour, so a defaulted diameter
        # silently rewrites the flow; unresolvable is an error (issue #262).
        area, _source = _infer_flow_area(graph, [self.from_node, self.to_node], exclude=self)
        if area is not None:
            self.diameter = math.sqrt(4.0 * area / math.pi)
        if self.diameter is None:
            if getattr(self, "_topo_pass", 0) >= 1:
                raise _unresolved_area_error(self, "diameter", [self.from_node, self.to_node])
            # Defer: an upstream tee may not have resolved F_C yet this pass.
            self._topo_pass = 1
            return
        self.area = math.pi * (self.diameter / 2) ** 2
        if self.Dh is None:
            self.Dh = self.diameter
        # Update convective surface area (was 0 if diameter was deferred)
        self._fill_convective_area()

    def default_convective_area(self) -> float:
        """The convective area a surface gets when none was set, by model.

        * Impingement array (one row): the row's target footprint matched to
          this channel's crossflow cross-section. The channel is z wide in
          the jet direction, so its span is ``A_flow / z`` and the footprint
          ``xn * span = A_flow * (xn/d) / (z/d)``; the hole count then
          follows from the channel itself (#462).
        * Single jet: the disc of radius ``R`` that Goldstein et al. (1986)
          average Nu over (Han Eq. 4.1 is an AVERAGE Nu out to R/D),
          ``pi (R_D d_jet)^2``.
        * Anything else: the wetted wall, ``pi D L``.
        """
        model = self.surface.model if self.surface else None
        if isinstance(model, ImpingementModel) and self.area and model.z_d > 0.0:
            return self.area * model.xn_d / model.z_d
        if isinstance(model, SingleJetImpingementModel):
            return math.pi * (model.R_D * model.d_jet) ** 2
        return math.pi * (self.diameter or 0.0) * self.length

    def _fill_convective_area(self) -> None:
        if self.surface and self.surface.area == 0.0:
            self.surface.area = self.default_convective_area()


class CrossflowSegmentElement(ChannelElement):
    """A channel segment between side-stream STATIONS: merges and bleeds (#465, #471).

    A ``ChannelElement`` (friction, stagnation-pressure coupling) plus the
    momentum of side streams joining or leaving at its end nodes. Over one
    station, ``m_a`` arriving and ``m_b`` leaving along the channel of area
    A, with the side stream's AXIAL velocity ``kappa * u_a``, momentum gives
    the static-pressure drop

        dP = (m_b|m_b| - m_a|m_a| + kappa m_a (m_a - m_b)) / (rho A^2)

    * ``kappa = 0``: normal injection (impingement jets). Florschuetz's
      P + G^2/rho = const.
    * ``kappa = 0.75`` (``cb.STATION_KAPPA_BLEED_BASSETT``): bleed through the
      wall. Bassett, Winterbone & Pearson's (2001) separating straight-run
      coefficient K2/K5 (Eq. 15), theta- and psi-independent, reproduced for
      every split.

    CENTRED: a segment carries half of the station at each end, each at its
    node's own density (C++ ``station_half_drop``), so every node holds the
    static pressure at the middle of its station. Ends are explicit, wired by
    whoever builds the chain:

    * ``from_kappa`` with ``prev_seg``: a station at ``from_node``; the flow
      arriving there is ``prev_seg``'s (0 if None: a closed end).
    * ``to_kappa`` with ``next_seg``: a station at ``to_node``; the flow
      leaving it is ``next_seg``'s.
    * ``entry_K``: ``from_node`` is a reservoir; the segment carries the
      entry acceleration ``(1 + K_in) m^2/(2 rho A^2)`` at the duct (to-node)
      density instead of a from-station.
    * an exit into a plenum (no to-station) loses the dynamic head: static
      continuity, as any Pt-coupled channel into a plenum.

    Station nodes must be PlenumNodes (Pt = P): the node pressure is the
    channel's static pressure, which is what a side stream sees.

    FRICTION follows ``ChannelElement``: Darcy-Weisbach on the
    equivalent-area diameter (the C++ channel's circular convention, #463).
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        length: float,
        area: float,
        Dh: float | None = None,
        roughness: float = 0.0,
        regime: CompressibilityLiteral = "incompressible",
        friction_model: FrictionModelLiteral = "haaland",
        from_kappa: float | None = None,
        to_kappa: float | None = None,
        prev_seg: str | None = None,
        next_seg: str | None = None,
        entry_K: float | None = None,
    ) -> None:
        if length <= 0.0 or area <= 0.0:
            raise ValueError(f"{type(self).__name__}: length and area must be positive")
        if entry_K is not None and from_kappa is not None:
            raise ValueError(
                f"{type(self).__name__} {id!r}: from_node is either a reservoir "
                "(entry_K) or a station (from_kappa), not both"
            )
        super().__init__(
            id,
            from_node,
            to_node,
            length=length,
            diameter=math.sqrt(4.0 * area / math.pi),
            Dh=Dh,
            roughness=roughness,
            regime=regime,
            friction_model=friction_model,
        )
        self.from_kappa = from_kappa
        self.to_kappa = to_kappa
        self.prev_seg = prev_seg
        self.next_seg = next_seg
        self.entry_K = entry_K
        self._prev: NetworkElement | None = None
        self._next: NetworkElement | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        super().resolve_topology(graph)
        self._prev = graph.elements.get(self.prev_seg) if self.prev_seg else None
        self._next = graph.elements.get(self.next_seg) if self.next_seg else None

    def network_flow_inputs(self) -> list[tuple[str, str]]:
        """The neighbour segments' flows, each read at the station it shares."""
        out = []
        if self._prev is not None:
            out.append((self._prev.id, self.from_node))
        if self._next is not None:
            out.append((self._next.id, self.to_node))
        return out

    def _momentum_terms(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None,
    ) -> list[tuple["cb.StationHalfDrop", str, str | None, str, NetworkElement | None]]:
        """Each momentum term: (C++ result, node it sits at, which derivative is
        this segment's own flow ('a' or 'b' or None), neighbour's node, the
        neighbour element)."""
        flows = flows or {}
        m = float(state_in.m_dot)
        terms = []
        if self.entry_K is not None:
            r = cb.channel_entry_drop(
                m, state_out.P, state_out.T, state_out.X, self.area, self.entry_K
            )
            terms.append((r, self.to_node, "a", self.to_node, None))
        if self.from_kappa is not None:
            m_a = flows.get(self._prev.id, 0.0) if self._prev is not None else 0.0
            r = cb.station_half_drop(
                m_a, m, state_in.P, state_in.T, state_in.X, self.area, self.from_kappa
            )
            terms.append((r, self.from_node, "b", self.from_node, self._prev))
        if self.to_kappa is not None and self._next is not None:
            m_b = flows.get(self._next.id, m)
            r = cb.station_half_drop(
                m, m_b, state_out.P, state_out.T, state_out.X, self.area, self.to_kappa
            )
            terms.append((r, self.to_node, "a", self.to_node, self._next))
        return terms

    def _momentum_drop_value(self, state_in, state_out, flows) -> float:
        return float(sum(t[0].dP for t in self._momentum_terms(state_in, state_out, flows)))

    def residuals(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        res, jac = super().residuals(state_in, state_out)
        row = jac[0]
        own = f"{self.id}.m_dot"
        for r, node, own_side, nb_node, nb in self._momentum_terms(state_in, state_out, flows):
            res[0] -= r.dP
            d_own = r.d_dm_a if own_side == "a" else r.d_dm_b
            row[own] = row.get(own, 0.0) - d_own
            if nb is not None:
                d_nb = r.d_dm_b if own_side == "a" else r.d_dm_a
                names = nb.unknowns()
                for local, coeff in nb.flow_jac_at_node(nb_node, list(range(len(names)))).items():
                    row[names[local]] = row.get(names[local], 0.0) - d_nb * coeff
            for key, d in ((f"{node}.P", r.d_dP), (f"{node}.T", r.d_dT)):
                row[key] = row.get(key, 0.0) - d
        return res, jac

    def htc_and_T(self, state: NetworkMixtureState, flows: dict[str, float] | None = None):
        return super().htc_and_T(state)

    def diagnostics(
        self,
        state_in: NetworkMixtureState,
        state_out: NetworkMixtureState,
        flows: dict[str, float] | None = None,
    ) -> dict[str, float | str]:
        out = super().diagnostics(state_in, state_out)
        out["dP_momentum"] = self._momentum_drop_value(state_in, state_out, flows)
        return out


class ImpingementCrossflowElement(CrossflowSegmentElement):
    """The crossflow channel between a jet plate and its target, one row pitch (#465).

    A ``CrossflowSegmentElement`` with normal-injection stations
    (``kappa = 0``): the jets that merge at each row arrive with NO
    streamwise momentum, so the crossflow must accelerate them -- the
    momentum equation behind Florschuetz, Truman and Metzger's (1981) flow
    distribution (Eqs. 7-8), P + G_c^2/rho = const. Without it every row sees
    the same pressure difference and the supply stays uniform; with it the
    downstream rows draw more. Area ``height * span``.

    Chain it with plenum crossflow nodes, one per row:

        c1 --[Crossflow]-- c2 --[Crossflow]-- c3 --[Crossflow]-- exit
         |                  |                  |
      [Plate 1]          [Plate 2]          [Plate 3]

    The neighbours are found from the chain itself: the crossflow segment
    arriving at ``from_node`` (none at row 1: a closed end) and the one
    leaving ``to_node`` (none into the exit, where the dynamic head is lost).
    Centred, each node holds its row's mid-merge static pressure: putting the
    whole merge downstream instead over-fed the downstream rows (Gc/Gj at row
    10 of Florschuetz's strongest-crossflow geometry 13% below Eq. 8, against
    a few percent centred).
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        length: float,
        height: float,
        span: float,
        roughness: float = 0.0,
        regime: CompressibilityLiteral = "incompressible",
        friction_model: FrictionModelLiteral = "haaland",
    ) -> None:
        if length <= 0.0 or height <= 0.0 or span <= 0.0:
            raise ValueError(
                "ImpingementCrossflowElement: length, height and span must be positive"
            )
        super().__init__(
            id,
            from_node,
            to_node,
            length=length,
            area=height * span,
            Dh=2.0 * height * span / (height + span),
            roughness=roughness,
            regime=regime,
            friction_model=friction_model,
            from_kappa=cb.STATION_KAPPA_MERGE_NORMAL,
            to_kappa=cb.STATION_KAPPA_MERGE_NORMAL,
        )
        self.height = height
        self.span = span

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        prev = [
            e
            for e in graph.get_upstream_elements(self.from_node)
            if isinstance(e, ImpingementCrossflowElement)
        ]
        nxt = [
            e
            for e in graph.get_downstream_elements(self.to_node)
            if isinstance(e, ImpingementCrossflowElement)
        ]
        self.prev_seg = prev[0].id if prev else None
        self.next_seg = nxt[0].id if nxt else None
        super().resolve_topology(graph)


class AreaChangeElement(NetworkElement):
    """
    Area change element applying pressure drop for sharp or conical area changes.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        F0: float | None = None,
        F1: float | None = None,
        model_type: Literal["sharp", "conical"] = "sharp",
        length: float | None = None,
        D_h: float | None = 0.0,
    ):
        super().__init__(id, from_node, to_node)
        self.F0: float | None = F0
        self.F1: float | None = F1
        self.model_type = model_type
        self.length = length if length is not None else 0.0
        self.D_h = D_h
        # Where a resolved area came from, for diagnostics (issue #262).
        self._F0_source = "user" if F0 is not None else ""
        self._F1_source = "user" if F1 is not None else ""

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # F0/F1 ARE the physics here: their ratio is the area change the
        # element exists to model. A defaulted value invents a contraction or
        # expansion the network does not contain, so an unresolvable area is
        # an error rather than a guess (issue #262).
        if self.F0 is None:
            self.F0, self._F0_source = _infer_flow_area(graph, [self.from_node], exclude=self)
            if self.F0 is None:
                raise _unresolved_area_error(self, "upstream area F0", [self.from_node])
        if self.F1 is None:
            self.F1, self._F1_source = _infer_flow_area(graph, [self.to_node], exclude=self)
            if self.F1 is None:
                raise _unresolved_area_error(self, "downstream area F1", [self.to_node])

    def validate(self) -> None:
        if (self.F0 or 0.0) <= 0:
            raise ValueError(
                f"AreaChangeElement '{self.id}' has invalid Upstream Area F0={self.F0}. Must be > 0."
            )
        if (self.F1 or 0.0) <= 0:
            raise ValueError(
                f"AreaChangeElement '{self.id}' has invalid Downstream Area F1={self.F1}. Must be > 0."
            )
        if self.D_h is not None and self.D_h < 0:
            raise ValueError(
                f"AreaChangeElement '{self.id}' has invalid Hydraulic Diameter D_h={self.D_h}. Must be >= 0."
            )

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        m_dot = state_in.m_dot

        if self.model_type == "sharp":
            res_cpp = _solver_tools.area_change_residuals_and_jacobian(
                m_dot=m_dot,
                P_total_up=state_in.Pt,
                P_static_up=state_in.P,
                T_up=state_in.T,
                Y_up=state_in.Y,
                P_static_down=state_out.P,
                F0=self.F0,
                F1=self.F1,
                m_scale=1e-4,
                D_h=self.D_h or 0.0,
            )
        else:
            res_cpp = _solver_tools.conical_area_change_residuals_and_jacobian(
                m_dot=m_dot,
                P_total_up=state_in.Pt,
                P_static_up=state_in.P,
                T_up=state_in.T,
                Y_up=state_in.Y,
                P_static_down=state_out.P,
                F0=self.F0,
                F1=self.F1,
                length=self.length,
                m_scale=1e-4,
            )

        res = [state_in.Pt - state_out.P - res_cpp.dP_calc]

        upstream_id = self.from_node
        downstream_id = self.to_node

        jac = {
            0: {
                f"{self.id}.m_dot": -res_cpp.d_dP_d_mdot,
                f"{upstream_id}.Pt": 1.0,
                f"{upstream_id}.P": -res_cpp.d_dP_dP_static_up,
                f"{upstream_id}.T": -res_cpp.d_dP_dT_up,
                f"{downstream_id}.P": -1.0,
            }
        }
        for i, val in enumerate(res_cpp.d_dP_dY_up):
            jac[0][f"{upstream_id}.Y[{i}]"] = -val

        return res, jac

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        try:
            cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        except Exception:
            return {
                "dP_static": 0.0,
                "area_ratio": self.F1 / self.F0 if self.F0 > 0 else 0.0,
            }

        f0 = self.F0
        f1 = self.F1
        f_small = min(f0, f1)
        area_ratio = f1 / f0 if f0 > 0 else 0.0
        m_dot = state_in.m_dot
        m_dot_abs = abs(m_dot)

        rho_in = cs_in.thermo.rho
        a_in = cs_in.thermo.a

        # Inlet velocity / Mach (at larger section F0)
        v_in = (m_dot_abs) / ((rho_in) * (f0)) if f0 > 0 and rho_in > 0 else 0.0
        mach_in = v_in / a_in if a_in > 0 else 0.0

        # Mach at the smaller section - used for loss coefficient lookup
        v_loss_ref = (m_dot_abs) / ((rho_in) * (f_small)) if f_small > 0 and rho_in > 0 else 0.0
        mach_loss_ref = v_loss_ref / a_in if a_in > 0 else 0.0

        # Outlet Mach (use outlet state for accuracy)
        try:
            cs_out = cb.complete_state(state_out.T, state_out.P, state_out.X)
            v_out = (
                (m_dot_abs) / ((cs_out.thermo.rho) * (f1))
                if f1 > 0 and cs_out.thermo.rho > 0
                else 0.0
            )
            mach_out = v_out / cs_out.thermo.a if cs_out.thermo.a > 0 else 0.0
        except Exception:
            mach_out = 0.0

        # Physics result for dP_static
        if self.model_type == "sharp":
            res = cb.sharp_area_change(
                m_dot, rho_in, cs_in.transport.mu, f0, f1, Mach=mach_loss_ref, D_h=self.D_h or 0.0
            )
        else:
            res = cb.conical_area_change(
                m_dot, rho_in, cs_in.transport.mu, f0, f1, self.length, Mach=mach_loss_ref
            )

        dP_static = abs(res.dP)
        q_small = 0.5 * rho_in * v_loss_ref**2
        zeta = dP_static / q_small if q_small > 1e-6 else 0.0

        Dh_in = self.D_h or math.sqrt(4.0 * f0 / math.pi)
        ref = _element_reference_block(
            cs_in, m_dot=m_dot, area=f0, Dh=Dh_in, location="small_section"
        )

        return {
            "m_dot": float(m_dot),
            **_element_pressure_block(state_in, state_out, mach_in=mach_in, mach_out=mach_out),
            **ref,
            "dP_static": float(dP_static),
            "area_ratio": float(area_ratio),
            "zeta": float(zeta),
            "mach_loss_ref": float(mach_loss_ref),
        }


class TeeJunctionElement(NetworkElement):
    """
    Three-port tee junction element using Bassett 2001 pressure-loss coefficients.

    Superseded by ``MultiPortChamberBase`` / ``MultiPortChamberElement`` (momentum-CV
    junction, #177). Retained as the M -> 0 regression baseline for the
    junction validation suite (see ``validation/junction/models/
    tee_junction_element_network.py``) and for the existing Python tests
    under ``python/tests/test_tee_network.py``. Not exposed through the GUI:
    new user-drawn junctions use ``mpce_tee`` instead.

    Merging (tee_type="merging"): straight_node and branch_node are inlets;
    common_node is the outlet.

    Branching (tee_type="branching"): common_node is the inlet; straight_node
    and branch_node are outlets.

    Unknowns: m_dot_com (total common flow) and m_dot_branch (branch flow).
    Straight flow is implicit: m_dot_straight = m_dot_com - m_dot_branch.

    psi = F_C / F_branch (common-arm area / lateral-branch area, default 1.0).
    Bassett 2001 derives K5/K11 (straight-arm coefficients) under the assumption
    F_straight = F_C; the straight-arm area is therefore not a free parameter.
    """

    def __init__(
        self,
        id: str,
        common_node: str,
        straight_node: str,
        branch_node: str,
        theta: float,
        F_C: float | None = None,
        F_branch: float | None = None,
        psi: float = 1.0,  # fallback ratio when F_branch is None and not inherited
        tee_type: Literal["merging", "branching"] = "merging",
        blend_k: float = 30.0,
    ):
        from_node = straight_node if tee_type == "merging" else common_node
        to_node = common_node if tee_type == "merging" else straight_node
        super().__init__(id, from_node, to_node)
        self.common_node = common_node
        self.straight_node = straight_node
        self.branch_node = branch_node
        self.theta = theta
        self.F_C: float | None = F_C
        self._F_branch: float | None = (
            F_branch  # resolved branch area; psi derived from F_C/_F_branch
        )
        self.psi = psi
        self.tee_type = tee_type
        self.blend_k = blend_k

    def all_source_nodes(self) -> list[str]:
        if self.tee_type == "merging":
            return [self.straight_node, self.branch_node]
        return [self.common_node]

    def all_sink_nodes(self) -> list[str]:
        if self.tee_type == "merging":
            return [self.common_node]
        return [self.straight_node, self.branch_node]

    def flow_at_node(self, node_id: str, x: Any, indices: list[int]) -> float:
        m_com = float(x[indices[0]])
        m_branch = float(x[indices[1]])
        if node_id == self.common_node:
            return m_com
        if node_id == self.branch_node:
            return m_branch
        return m_com - m_branch  # straight_node

    def flow_jac_at_node(self, node_id: str, indices: list[int]) -> dict[int, float]:
        if node_id == self.common_node:
            return {indices[0]: 1.0}
        if node_id == self.branch_node:
            return {indices[1]: 1.0}
        return {indices[0]: 1.0, indices[1]: -1.0}  # straight_node

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot_com", f"{self.id}.m_dot_branch"]

    def n_equations(self) -> int:
        return 2

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Step 1: resolve F_C from common or straight arm channels.
        # Branch arm is only used as last resort (back-calculate F_C = F_branch * psi).
        if self.F_C is None:
            found_fc = False
            for node_id in (self.common_node, self.straight_node):
                for e in graph.get_upstream_elements(node_id) + graph.get_downstream_elements(
                    node_id
                ):
                    if e is self:
                        continue
                    if isinstance(e, ChannelElement) and e.diameter is not None:
                        self.F_C = math.pi * (e.diameter / 2.0) ** 2
                        found_fc = True
                        break
                if found_fc:
                    break
            if not found_fc:
                # Last resort: use branch arm channel to back-calculate F_C = F_branch * psi
                for e in graph.get_upstream_elements(
                    self.branch_node
                ) + graph.get_downstream_elements(self.branch_node):
                    if e is self:
                        continue
                    if isinstance(e, ChannelElement) and e.diameter is not None:
                        f_b = math.pi * (e.diameter / 2.0) ** 2
                        if self._F_branch is None:
                            self._F_branch = f_b
                        self.F_C = f_b * self.psi
                        found_fc = True
                        break
            if not found_fc:
                # Widen to any area-bearing neighbour before giving up: the
                # common/straight arms may be fed by a chamber or an area
                # change rather than a channel (issue #262).
                area, _source = _infer_flow_area(
                    graph, [self.common_node, self.straight_node], exclude=self
                )
                if area is not None:
                    self.F_C = area
                    found_fc = True
            if not found_fc:
                # F_C sets every loss reference in the junction closure, so a
                # default would silently rescale all of them.
                if getattr(self, "_topo_pass", 0) >= 1:
                    raise _unresolved_area_error(
                        self, "common-arm area F_C", [self.common_node, self.straight_node]
                    )
                self._topo_pass = 1

        # Step 2: resolve F_branch independently from the branch arm channel.
        if self._F_branch is None:
            for e in graph.get_upstream_elements(self.branch_node) + graph.get_downstream_elements(
                self.branch_node
            ):
                if e is self:
                    continue
                if isinstance(e, ChannelElement) and e.diameter is not None:
                    self._F_branch = math.pi * (e.diameter / 2.0) ** 2
                    break

        # Step 3: recompute psi from resolved areas.
        if self.F_C is not None and self._F_branch is not None and self._F_branch > 0:
            self.psi = self.F_C / self._F_branch

        # Step 4: propagate duct area to auto-area MCN nodes at each arm.
        # MCN nodes at straight/branch ports cannot see an upstream channel in
        # their own resolve_topology (the upstream is this tee element), so they
        # fall back to the 0.1 m^2 sentinel.  The tee already knows F_C and
        # F_branch, so we inject the correct geometry here.
        if self.F_C is not None and self.F_C > 0:
            D_main = math.sqrt(4.0 * self.F_C / math.pi)
            for nid in (self.common_node, self.straight_node):
                _node = graph.nodes.get(nid)
                if isinstance(_node, MomentumChamberNode) and _node._auto_area and _node.Dh is None:
                    _node.Dh = D_main
                    _node.area = self.F_C
                    _node.surface.area = _node.area
        if self._F_branch is not None and self._F_branch > 0:
            D_branch = math.sqrt(4.0 * self._F_branch / math.pi)
            _bra_node = graph.nodes.get(self.branch_node)
            if (
                isinstance(_bra_node, MomentumChamberNode)
                and _bra_node._auto_area
                and _bra_node.Dh is None
            ):
                _bra_node.Dh = D_branch
                _bra_node.area = self._F_branch
                _bra_node.surface.area = _bra_node.area

    def validate(self) -> None:
        if self.tee_type not in ("merging", "branching"):
            raise ValueError(
                f"TeeJunctionElement '{self.id}': tee_type must be 'merging' or 'branching'."
            )
        if not (abs(self.theta) <= math.pi / 2.0):
            raise ValueError(
                f"TeeJunctionElement '{self.id}': |theta| must be <= pi/2, got {self.theta}."
            )
        if (self.F_C or 0.0) <= 0.0:
            raise ValueError(f"TeeJunctionElement '{self.id}': F_C must be > 0, got {self.F_C}.")
        if self.psi <= 0.0:
            raise ValueError(f"TeeJunctionElement '{self.id}': psi must be > 0, got {self.psi}.")

    def _make_branch_input(
        self,
        state: NetworkMixtureState,
        m_dot: float,
        A: float,
        theta: float,
    ) -> "cb._core.BranchInput":
        bi = cb._core.BranchInput()
        bi.P_static = state.P
        bi.Pt = state.Pt
        bi.T = state.T
        bi.m_dot = m_dot
        bi.A = A
        bi.theta = theta
        bi.gamma_eff = state.gamma()
        bi.R_gas = float(cb.specific_gas_constant(state.X))
        return bi

    def residuals(
        self,
        state_com: NetworkMixtureState,
        state_straight: NetworkMixtureState,
        state_branch: NetworkMixtureState,
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        m_com = state_com.m_dot
        m_branch = state_branch.m_dot
        m_straight = m_com - m_branch
        A_branch = self._F_branch if self._F_branch is not None else self.F_C / max(self.psi, 1e-6)

        if self.tee_type == "merging":
            # Suppliers: straight (m_straight > 0) and branch (m_branch > 0)
            # Collector: common (m_com > 0 is the outlet)
            # BranchInput sign: >0 = supplier. Straight and branch supply; common collects.
            bra_com = self._make_branch_input(state_com, -m_com, self.F_C, 0.0)
            bra_str = self._make_branch_input(state_straight, m_straight, self.F_C, 0.0)
            bra_bra = self._make_branch_input(state_branch, m_branch, A_branch, abs(self.theta))
            res = cb._core.compressible_merging_tee_rj(bra_com, bra_str, bra_bra)
            # R_0 = P_str - P_bra  (R_sp: equal static pressure at merge point)
            # R_1 = p0_dat - Pt_com - K_com * q_ref  (R_K,com: loss to collector)
            jac: dict[int, dict[str, float]] = {
                0: {
                    f"{self.straight_node}.P": res.dR0_dP_str,
                    f"{self.branch_node}.P": res.dR0_dP_bra,
                },
                1: {
                    f"{self.common_node}.Pt": res.dR1_dPt_com,
                    f"{self.straight_node}.P": res.dR1_dP_str,
                    f"{self.straight_node}.T": res.dR1_dT_str,
                    f"{self.branch_node}.P": res.dR1_dP_bra,
                    f"{self.branch_node}.T": res.dR1_dT_bra,
                    f"{self.id}.m_dot_com": res.dR1_dmdot_com,
                    f"{self.id}.m_dot_branch": res.dR1_dmdot_branch,
                },
            }
        else:
            # Supplier: common. Collectors: straight and branch.
            bra_com = self._make_branch_input(state_com, m_com, self.F_C, 0.0)
            bra_str = self._make_branch_input(state_straight, m_straight, self.F_C, 0.0)
            bra_bra = self._make_branch_input(state_branch, m_branch, A_branch, abs(self.theta))
            res = cb._core.compressible_branching_tee_rj(bra_com, bra_str, bra_bra)
            jac = {
                0: {
                    f"{self.id}.m_dot_com": res.dR0_dmdot_com,
                    f"{self.id}.m_dot_branch": res.dR0_dmdot_branch,
                    f"{self.common_node}.Pt": res.dR0_dPt_com,
                    f"{self.common_node}.P": res.dR0_dP_com,
                    f"{self.common_node}.T": res.dR0_dT_com,
                    f"{self.straight_node}.Pt": res.dR0_dPt_str,
                },
                1: {
                    f"{self.id}.m_dot_com": res.dR1_dmdot_com,
                    f"{self.id}.m_dot_branch": res.dR1_dmdot_branch,
                    f"{self.common_node}.Pt": res.dR1_dPt_com,
                    f"{self.common_node}.P": res.dR1_dP_com,
                    f"{self.common_node}.T": res.dR1_dT_com,
                    f"{self.branch_node}.Pt": res.dR1_dPt_bra,
                    f"{self.branch_node}.P": res.dR1_dP_bra,
                    f"{self.branch_node}.T": res.dR1_dT_bra,
                },
            }

        return [res.R_0, res.R_1], jac

    def diagnostics(
        self,
        state_com: NetworkMixtureState,
        state_straight: NetworkMixtureState,
        state_branch: NetworkMixtureState,
    ) -> dict[str, float]:
        m_com = state_com.m_dot
        m_branch = state_branch.m_dot
        if self.tee_type == "merging":
            res = cb._core.merging_tee_residuals_and_jacobian(
                m_dot_com=m_com,
                m_dot_branch=m_branch,
                dP0_straight=state_straight.Pt - state_com.P,
                dP0_branch=state_branch.Pt - state_com.P,
                P_static_com=state_com.P,
                T_com=state_com.T,
                Y_com=state_com.Y,
                theta=abs(self.theta),
                psi=self.psi,
                F_C=self.F_C,
                blend_k=self.blend_k,
            )
        else:
            res = cb._core.branching_tee_residuals_and_jacobian(
                m_dot_com=m_com,
                m_dot_branch=m_branch,
                dP0_straight=state_com.Pt - state_straight.P,
                dP0_branch=state_com.Pt - state_branch.P,
                P_static_com=state_com.P,
                T_com=state_com.T,
                Y_com=state_com.Y,
                theta=abs(self.theta),
                psi=self.psi,
                F_C=self.F_C,
                blend_k=self.blend_k,
            )
        is_extrapolated = res.status != cb._core.CorrelationValidity.VALID
        return {
            "m_dot_com": float(m_com),
            "m_dot_straight": float(m_com - m_branch),
            "m_dot_branch": float(m_branch),
            "mass_flow_ratio": float(res.q),
            "K_straight": float(res.K_straight),
            "K_branch": float(res.K_branch),
            "correlation_extrapolated": 1.0 if is_extrapolated else 0.0,
        }


# ---------------------------------------------------------------------------
# Momentum-CV junction machinery (momentum cv implementation guide.pdf). The
# residuals live in the subclasses (mpce_element.py), whose closures carry the
# port losses; BorderCarnotLossElement below is NOT a companion for them.
# ---------------------------------------------------------------------------


class MultiPortChamberBase(NetworkElement):
    """
    Momentum-CV junction element (N >= 2 ports).

    ABSTRACT SINCE 0.6.0. This class owns the topology and port machinery --
    port ordering, the inlet/outlet sign map, area and angle resolution, the
    connecting-element wiring checks -- and no longer carries a residual of
    its own. Its impulse-function model was deprecated in 0.5.0 and removed
    here; ``MultiPortChamberElement`` supersedes it (issue #271).

    Subclasses own one scalar unknown ``{id}.P_jct`` and emit N per-port
    relations plus a global mass residual ``sum_i mdot_port_i = 0``, with a
    per-port orientation derived from the topology (positive = flow out of
    the junction). The per-port relation, and with it the port losses, is
    the subclass's closure -- Mynard's for ``MultiPortChamberElement``, fixed
    K for ``ConstantKTeeElement``. Do not add a
    :class:`BorderCarnotLossElement` on a port as well: that double-counts
    the turning loss (#272).

    Topology contract (see PDF Section 2 + addendum task #10):

    - Each port must be a :class:`MomentumChamberNode` (carries P, Pt at the
      port face). The MCN is marked ``_is_junction_port = True`` during
      :meth:`resolve_topology` so the solver skips its mass-balance row -- the
      junction's R_mass takes over conservation at the port (otherwise the
      MCN's mass row degenerates to 0 = 0 once the junction declares the
      same flow with opposite sign).
    - Exactly one element (channel, loss element, etc.) must connect to each
      port-MCN on the outside. Its m_dot is the port's mass-flow source; the
      junction reads it via the graph at residual-eval time, and reports it
      back through ``flow_at_node`` so the state propagation gives every
      port-MCN its real face flow (collector ports would otherwise see zero
      and their Pt = P + 0.5*rho*v^2 closure would degenerate to Pt = P).
    - Sign convention per port is fixed by the inlet/outlet declaration in
      the constructor. ``inlet_nodes`` carry sign ``-1`` (canonical flow INTO
      junction); ``outlet_nodes`` carry sign ``+1`` (canonical flow OUT). The
      connecting outer element's ``from_node``/``to_node`` is then checked
      against the declared direction at :meth:`resolve_topology` and a clear
      error is raised on mismatch. Runtime flow may still reverse (mdot < 0
      at any port); only the canonical direction is fixed for graph topology.

    Port areas are inherited from the connecting :class:`ChannelElement`'s
    diameter unless explicitly given.

    See ``docs/junction/momentum cv implementation guide.pdf`` Sections 2-3
    and ``docs/junction/junction_model_v3_addendum.md`` Finding 6 for the
    motivation (the K-closure's structural K_run ~ 0.30 ceiling and its
    inability to represent ejector / merge-split direction natively).
    """

    # Opt-in for the solver's junction split seed (`_junction_split_guess`).
    # True here because a chamber junction's port flows follow from the port
    # pressure differences: through the element's own closure where it offers
    # `closure_split_seed` (#272), else a Bernoulli share. Subclasses whose
    # port flows are set by their own internal physics must turn it off -- see
    # `EjectorElement`.
    seeds_ports_by_pressure_split: bool = True

    def __init__(
        self,
        id: str,
        inlet_nodes: list[str],
        outlet_nodes: list[str],
        inlet_angles_deg: list[float] | None = None,
        outlet_angles_deg: list[float] | None = None,
        port_areas: list[float] | None = None,
    ):
        """Build a momentum-CV junction.

        Args:
            id: Element id.
            inlet_nodes: Port-MCN ids where flow CANONICALLY enters the
                junction (junction is graph-downstream of each). Flow can
                reverse at runtime; this just sets the canonical topology
                direction for the graph validator.
            outlet_nodes: Port-MCN ids where flow CANONICALLY leaves the
                junction.
            inlet_angles_deg: Geometric branch angles for the inlet ports
                (length must match inlet_nodes). Default 0 each.
            outlet_angles_deg: Same, for outlet ports. Default 0 each.
            port_areas: Per-port cross-section areas [m^2], ordered as
                ``inlet_nodes + outlet_nodes``. If None, inherited from
                connecting channels at resolve time.
        """
        if not inlet_nodes:
            raise ValueError(f"MultiPortChamberBase '{id}': need >= 1 inlet port.")
        if not outlet_nodes:
            raise ValueError(f"MultiPortChamberBase '{id}': need >= 1 outlet port.")
        n_in = len(inlet_nodes)
        n_out = len(outlet_nodes)

        # NetworkElement base class requires from_node/to_node; use the first
        # inlet and the first outlet as placeholders. The all_source_nodes /
        # all_sink_nodes overrides below carry the real topology.
        super().__init__(id, inlet_nodes[0], outlet_nodes[0])

        self.inlet_nodes = list(inlet_nodes)
        self.outlet_nodes = list(outlet_nodes)
        self.port_nodes = self.inlet_nodes + self.outlet_nodes  # ordered: inlets first
        self.N = len(self.port_nodes)

        # Canonical sign convention (PDF Section 2.2: positive = out of junction):
        #   inlet ports: flow into junction at canonical orientation -> -1
        #   outlet ports: flow out of junction at canonical orientation -> +1
        self._port_signs: list[float] = [-1.0] * n_in + [+1.0] * n_out

        inlet_angles = list(inlet_angles_deg) if inlet_angles_deg is not None else [0.0] * n_in
        outlet_angles = list(outlet_angles_deg) if outlet_angles_deg is not None else [0.0] * n_out
        if len(inlet_angles) != n_in:
            raise ValueError(
                f"MultiPortChamberBase '{id}': inlet_angles_deg has "
                f"{len(inlet_angles)} entries, need {n_in}."
            )
        if len(outlet_angles) != n_out:
            raise ValueError(
                f"MultiPortChamberBase '{id}': outlet_angles_deg has "
                f"{len(outlet_angles)} entries, need {n_out}."
            )
        self.port_angles_deg = inlet_angles + outlet_angles

        self.port_areas: list[float | None] = (
            list(port_areas) if port_areas is not None else [None] * self.N
        )
        if len(self.port_areas) != self.N:
            raise ValueError(
                f"MultiPortChamberBase '{id}': port_areas has "
                f"{len(self.port_areas)} entries, need {self.N}."
            )

        # Resolved during resolve_topology:
        #   self._port_element_ids[i] = id of the element connected at port i
        self._port_element_ids: list[str] = [""] * self.N

    def unknowns(self) -> list[str]:
        return [f"{self.id}.P_jct"]

    def n_equations(self) -> int:
        return self.N + 1  # N impulse + 1 sum-mass

    def row_scale_kinds(self) -> list[str]:
        """Per-row scale kind ("p" or "mdot") for the solver's row-scaling
        vector (RESIDUAL-ROW-SCALING, see NetworkSolver._build_residual_scales),
        in the same order as `residuals()`'s returned list.

        Default matches this class's own impulse-CV rows: N pressure-
        magnitude rows (P_i + rho_i u_i^2 - P_jct) followed by one mass-
        magnitude row (sum of port mdots). Subclasses whose residual rows
        have a DIFFERENT type pattern (e.g. a row expressed as a mass-flow
        residual rather than a pressure residual) must override this --
        the solver applies the wrong scale silently otherwise, which has
        previously caused severe convergence stalls (see solver.py's
        scaling block for the historical "doomed primary" stall class this
        guards against).
        """
        return ["p"] * self.N + ["mdot"]

    def all_source_nodes(self) -> list[str]:
        # Inlet ports = nodes the junction draws flow FROM (canonical orientation).
        return list(self.inlet_nodes)

    def all_sink_nodes(self) -> list[str]:
        # Outlet ports = nodes the junction delivers flow TO.
        return list(self.outlet_nodes)

    def _port_throughflow_terms(self, node_id: str) -> list[tuple[str, float]]:
        """Outer-element terms whose weighted sum is the junction throughflow
        at ``node_id`` (positive = canonical direction through the port).

        Supplier (inlet) ports and collector ports of a single-supplier
        junction carry exactly their own outer element's flow. A collector
        of a multi-supplier junction (e.g. merge tee outlet) receives the
        sum of the supplier feeds -- consistent with how the solver's state
        propagation composes the collector's total flow from one stream per
        supplier port.
        """
        if node_id not in self.port_nodes:
            return []
        i = self.port_nodes.index(node_id)
        if self._port_signs[i] < 0 or len(self.inlet_nodes) == 1:
            return [(self._port_element_ids[i], 1.0)]
        return [(self._port_element_ids[j], 1.0) for j in range(self.N) if self._port_signs[j] < 0]

    def flow_at_node(self, node_id: str, x: Any, indices: list[int]) -> float:
        # Port-MCN mass rows are skipped by the solver (_is_junction_port);
        # this is consumed by the state propagation instead, so collector-port
        # MCNs see their real face flow in the Pt = P + 0.5*rho*v^2 closure
        # and mix upstream streams with true mass weights. `indices` are this
        # element's own unknowns (P_jct); the port flows live on the outer
        # elements, resolved via the _port_outer_mdot_idx map the solver
        # stashes when it builds the unknown vector.
        idx_map = getattr(self, "_port_outer_mdot_idx", None)
        if not idx_map:
            return 0.0
        total = 0.0
        for outer_id, coeff in self._port_throughflow_terms(node_id):
            global_idx = idx_map.get(outer_id)
            if global_idx is not None:
                total += coeff * float(x[global_idx])
        return total

    def flow_jac_at_node(self, node_id: str, indices: list[int]) -> dict[int, float]:
        idx_map = getattr(self, "_port_outer_mdot_idx", None)
        if not idx_map:
            return {}
        jac: dict[int, float] = {}
        for outer_id, coeff in self._port_throughflow_terms(node_id):
            global_idx = idx_map.get(outer_id)
            if global_idx is not None:
                jac[global_idx] = jac.get(global_idx, 0.0) + coeff
        return jac

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # 1. For each port: locate the (unique) connecting element on the
        #    OUTSIDE of the port-MCN and determine the orientation sign.
        for i, port_id in enumerate(self.port_nodes):
            port_node = graph.nodes.get(port_id)
            if port_node is None:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port node "
                    f"'{port_id}' is not in the network."
                )
            if not isinstance(port_node, MomentumChamberNode):
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port '{port_id}' "
                    f"must be a MomentumChamberNode, got "
                    f"'{type(port_node).__name__}'."
                )
            # Mark the port-MCN so the solver skips its mass row.
            port_node._is_junction_port = True
            port_node._junction_id = self.id

            outside_elems = [
                e
                for e in (
                    graph.get_upstream_elements(port_id) + graph.get_downstream_elements(port_id)
                )
                if e is not self
            ]
            if len(outside_elems) != 1:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port '{port_id}' "
                    f"must have exactly one outside element (channel, loss "
                    f"element, etc.), found {len(outside_elems)}."
                )
            outer = outside_elems[0]
            self._port_element_ids[i] = outer.id

            # Validate the outer element's orientation matches the declared
            # inlet/outlet for this port. The sign convention was fixed at
            # __init__ from the inlet/outlet split:
            #   inlet port (sign = -1):  outer.to_node should be the port-MCN
            #     (outer feeds flow INTO the port, junction takes it from there).
            #   outlet port (sign = +1): outer.from_node should be the port-MCN
            #     (outer takes flow OUT of port, junction sends flow to there).
            expected_sign = self._port_signs[i]
            if outer.from_node == port_id:
                topo_sign = +1.0  # outer takes flow OUT of port = +1 outflow
            elif outer.to_node == port_id:
                topo_sign = -1.0  # outer feeds flow INTO port = -1 outflow
            else:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': connecting element "
                    f"'{outer.id}' at port '{port_id}' has no direct from/to "
                    f"link to that port (multi-port outer elements are not "
                    f"supported)."
                )
            if topo_sign != expected_sign:
                role = "inlet" if expected_sign < 0 else "outlet"
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': port '{port_id}' is "
                    f"declared as an {role} (sign={expected_sign:+.0f}) but the "
                    f"connecting element '{outer.id}' is wired in the opposite "
                    f"direction (sign={topo_sign:+.0f}). For an inlet port, the "
                    f"outer element should feed INTO the port-MCN (outer.to_node "
                    f"== port). For an outlet port, the outer element should take "
                    f"flow OUT of the port-MCN (outer.from_node == port)."
                )

            # 2. Inherit port area from connecting channel if not set.
            if self.port_areas[i] is None:
                if isinstance(outer, ChannelElement) and outer.diameter is not None:
                    self.port_areas[i] = math.pi * (outer.diameter / 2.0) ** 2
                elif isinstance(outer, BorderCarnotLossElement) and outer.area is not None:
                    self.port_areas[i] = outer.area
            if self.port_areas[i] is None:
                raise ValueError(
                    f"MultiPortChamberBase '{self.id}': could not infer "
                    f"area at port '{port_id}'. Provide port_areas explicitly "
                    f"or attach a ChannelElement with diameter set."
                )

            # 3. Give auto-sized port-MCNs a real flow area. Nodes resolve
            #    before elements (FlowNetwork.resolve_all_topology), and a
            #    collector port has no upstream channel to inherit Dh from,
            #    so it would otherwise keep MomentumChamberNode's 0.1 m^2
            #    fallback and its Pt = P + 0.5*rho*v^2 closure would see a
            #    near-zero face velocity. Mirrors the MCN auto path (Dh ->
            #    area, surface.area) exactly so a later node-resolve pass is
            #    idempotent and collector ports match what inlet ports
            #    already get from upstream-channel Dh inheritance.
            if port_node._auto_area and port_node.Dh is None:
                area_i = float(self.port_areas[i])
                port_node.Dh = 2.0 * math.sqrt(area_i / math.pi)
                port_node.area = area_i
                port_node.surface.area = area_i

    def diagnostics(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float] | None = None,
    ) -> dict[str, float]:
        """Per-port bookkeeping: pressures, temperatures, areas, signs, flows.

        MACHINERY, not model. It reports the port map and the states handed
        in and computes no physics, which is why it survives the 0.6.0
        removal of this class's residual. ``MultiPortChamberElement.diagnostics``
        calls it through ``super()`` and adds its closure quantities on top.
        """
        """DEPRECATED, removal in 0.6.0 -- see ``_warn_v1_junction_model``.
        Note that ``MultiPortChamberElement.diagnostics`` delegates here via ``super()``,
        which is why the warning is guarded on the exact type rather than on
        ``isinstance``.
        """
        diag: dict[str, float] = {"P_jct": float(P_jct), "n_ports": float(self.N)}
        for i in range(self.N):
            diag[f"port_{i}_P"] = float(states[i].P)
            diag[f"port_{i}_T"] = float(states[i].T)
            diag[f"port_{i}_area"] = float(self.port_areas[i] or 0.0)
            diag[f"port_{i}_sign"] = float(self._port_signs[i])
            if port_mdots is not None and i < len(port_mdots):
                diag[f"port_{i}_m_dot"] = float(port_mdots[i])
        return diag


class BorderCarnotLossElement(NetworkElement):
    """
    Border-Carnot turning-loss element: a two-port in-line loss. It was the
    lateral-port companion of the momentum-CV junction's impulse model,
    removed in 0.6.0 (#322).

    NOT FOR JUNCTION PORTS. ``MultiPortChamberElement`` and
    ``ConstantKTeeElement`` carry each port's turning loss in their own
    closure, and adding this element double-counts it: on Bassett's 90 deg,
    psi = 1 dividing tee, K6 at q = 0.5 is 0.867 from the closure alone
    (Bassett 0.867) and 1.249 with this element (+16% at q = 0.3, +78% at
    0.7; #272).

    Residual (PDF Section 3.1):

        Pt_in - Pt_out - L(delta_geom) * 0.5 * rho_in * u_in^2 = 0
        L = 4 * (1 - cos((3/4) * delta_geom))^2

    The (3/4) factor is Hager's effective-angle correction for sharp-edged
    lateral branches. The form is *intended* to reproduce Hager xi_l and
    Bassett K_inc at M -> 0 on a sharp-edged 90-deg lateral. That was never
    demonstrated: its only check, paired with the removed junction model,
    missed Bassett K6 by 11-29%, and the tests went with that model (#322).
    Do not rely on it being anchored to Hager or Bassett. delta_geom = 0
    gives L = 0.

    The element is sign-free in m_dot (mdot^2 in the dynamic head). The
    initial form has no direction asymmetry parameter; an
    ``alpha_contract`` calibration multiplier can be added later if validation
    data forces it (PDF Section 3.3, deferred).
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        delta_geom_deg: float,
        area: float | None = None,
    ):
        super().__init__(id, from_node, to_node)
        self.delta_geom_deg = float(delta_geom_deg)
        self.area = area
        self._area_source = "user" if area is not None else ""

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # This element already refused to guess an area before #262; it now
        # shares the widened inference path, so a chamber or area-change
        # neighbour resolves it where only a channel used to.
        if self.area is None:
            self.area, self._area_source = _infer_flow_area(
                graph, [self.from_node, self.to_node], exclude=self
            )
            if self.area is None:
                raise _unresolved_area_error(self, "area", [self.from_node, self.to_node])

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        delta_rad = math.radians(self.delta_geom_deg)
        cpp = _solver_tools.border_carnot_loss_residual_and_jacobian(
            mdot=state_in.m_dot,
            Pt_in=state_in.Pt,
            Pt_out=state_out.Pt,
            P_in=state_in.P,
            T_in=state_in.T,
            Y_in=list(state_in.Y),
            area=float(self.area),
            delta_geom=delta_rad,
        )
        res = [cpp.residual]
        jac: dict[int, dict[str, float]] = {
            0: {
                f"{self.from_node}.Pt": cpp.d_res_dPt_in,
                f"{self.to_node}.Pt": cpp.d_res_dPt_out,
                f"{self.from_node}.P": cpp.d_res_dP_in,
                f"{self.from_node}.T": cpp.d_res_dT_in,
                f"{self.id}.m_dot": cpp.d_res_dmdot,
            }
        }
        return res, jac

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        return {
            "m_dot": float(state_in.m_dot),
            "Pt_in": float(state_in.Pt),
            "Pt_out": float(state_out.Pt),
            "dPt": float(state_in.Pt - state_out.Pt),
            "delta_geom_deg": float(self.delta_geom_deg),
            "area": float(self.area or 0.0),
        }


# ---------------------------------------------------------------------------
# Vortex element: rotating-cavity centrifugal pressure rise (Vatistas 1991)
# ---------------------------------------------------------------------------


class VortexElement(NetworkElement):
    """
    Models the centrifugal pressure rise of a rotating cavity or disc
    using the Vatistas (1991) n-vortex model.

    The shaft angular velocity omega drives a vortex whose circulation is
    Gamma = 2*pi * r_c^2 * omega * 2^(1/n), i.e. the core rotates as solid
    body at omega.  The pressure rise from r_in to r_out follows from
    integrating the centripetal acceleration:
        dP_rise = vatistas_delta_p(r_out, rho, Gamma, r_c, n)
                - vatistas_delta_p(r_in,  rho, Gamma, r_c, n)
    where rho is taken from the upstream node state.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        r_c: float = 0.02,
        r_out: float = 0.10,
        r_in: float = 0.0,
        omega_rpm: float = 0.0,
        n: float = 2.0,
    ):
        super().__init__(id, from_node, to_node)
        if n < 1.0:
            raise ValueError(f"Vatistas shape parameter n must be >= 1; got {n}")
        if r_c <= 0.0:
            raise ValueError(f"Core radius r_c must be > 0; got {r_c}")
        if r_out <= 0.0:
            raise ValueError(f"Outer radius r_out must be > 0; got {r_out}")
        if r_in < 0.0:
            raise ValueError(f"Inner radius r_in must be >= 0; got {r_in}")
        self.r_c = float(r_c)
        self.r_out = float(r_out)
        self.r_in = float(r_in)
        self.omega_rpm = float(omega_rpm)
        self.n = float(n)

    def _gamma(self) -> float:
        """Circulation [m^2/s] from shaft speed assuming solid-body vortex core."""
        omega = self.omega_rpm * 2.0 * math.pi / 60.0
        return 2.0 * math.pi * self.r_c**2 * omega * 2.0 ** (1.0 / self.n)

    def _dp_rise(self, rho: float) -> float:
        """Total pressure rise [Pa] at inlet density rho."""
        gamma = self._gamma()
        if gamma == 0.0:
            return 0.0
        dp_out = cb.vatistas_delta_p(self.r_out, rho, gamma, self.r_c, self.n)
        dp_in = (
            cb.vatistas_delta_p(self.r_in, rho, gamma, self.r_c, self.n) if self.r_in > 0.0 else 0.0
        )
        return dp_out - dp_in

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def n_equations(self) -> int:
        return 1

    def residuals(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        rho_raw = state_in.density()
        rho, d_rho_safe = _safe_rho(rho_raw)
        dp_rise = self._dp_rise(rho)

        res = [state_out.Pt - state_in.Pt - dp_rise]

        # d(dp_rise)/d(rho) = dp_rise / rho  (dP scales linearly with rho)
        # d(rho)/d(P)  ~  rho / P  for ideal gas
        # d(rho)/d(T)  ~ -rho / T  for ideal gas
        dp_drho = (dp_rise / rho) * d_rho_safe if rho > 0 else 0.0
        drho_dP = rho / state_in.P if state_in.P > 0 else 0.0
        drho_dT = -rho / state_in.T if state_in.T > 0 else 0.0

        jac: dict[int, dict[str, float]] = {
            0: {
                f"{self.to_node}.Pt": 1.0,
                f"{self.from_node}.Pt": -1.0,
                f"{self.from_node}.P": -dp_drho * drho_dP,
                f"{self.from_node}.T": -dp_drho * drho_dT,
            }
        }
        return res, jac

    def diagnostics(
        self, state_in: NetworkMixtureState, state_out: NetworkMixtureState
    ) -> dict[str, float]:
        if state_in.P <= 0 or state_in.T <= 0:
            return {}
        cs_in = cb.complete_state(state_in.T, state_in.P, state_in.X)
        rho = cs_in.thermo.rho
        dp_rise = self._dp_rise(rho)
        return {
            "m_dot": float(state_in.m_dot),
            "omega_rpm": float(self.omega_rpm),
            "dP_vortex": float(dp_rise),
            **_element_pressure_block(state_in, state_out, mach_in=0.0, mach_out=0.0),
        }
