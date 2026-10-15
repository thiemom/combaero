import math as _math

from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    CombustorNode,
    ConstantFractionLoss,
    ConstantHeadLoss,
    ConvectiveSurface,
    EffusionLiner,
    EffusionPlateElement,
    FlowNetwork,
    ImpingementArray,
    ImpingementModel,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberBase,
    OrificeElement,
    PinFinModel,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    RibbedModel,
    SingleJetImpingementModel,
    SmoothModel,
    ThermalWall,
    VortexElement,
    WallLayer,
    WallNode,
)
from combaero.network.components import DEFAULT_CHANNEL_DIAMETER
from combaero.network.ejector_element import EjectorElement
from combaero.network.ideal_tee_element import IdealMomentumTeeElement
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

from .schemas import (
    AreaChangeData,
    ChannelData,
    CombustorData,
    CompositionData,
    ConstantFractionLossData,
    ConstantHeadLossData,
    DiscreteLossData,
    EffusionLinerData,
    EffusionPlateData,
    EjectorData,
    ImpingementArrayData,
    ImpingementModelData,
    LinearThetaFractionLossData,
    LinearThetaHeadLossData,
    MassBoundaryData,
    MomentumChamberData,
    MPCETeeData,
    NetworkGraphSchema,
    OrificeData,
    PinFinModelData,
    PlenumData,
    PressureBoundaryData,
    RibbedModelData,
    SingleJetImpingementModelData,
    SmoothModelData,
    ThermalWallData,
    VortexData,
    WallData,
)

#: |channel area / port area - 1| above which a tee port and its directly
#: connected channel are joined through a sharp step (#536). Below it the bores
#: count as the same size, so round-off in a typed diameter inserts no step.
#: The same bound as the junction's unmodelled-step diagnostic (#539), so a
#: step the builder leaves out is never one that diagnostic reports.
TEE_PORT_STEP_TOL = MultiPortChamberBase.PORT_AREA_MISMATCH_TOL


def build_pressure_loss(data, area: float = 0.1):
    """Converts PressureLossData schema to a pressure-loss callable."""
    if data is None:
        return None
    t = getattr(data, "type", None)
    if isinstance(data, ConstantFractionLossData) or t == "constant_fraction":
        return ConstantFractionLoss(xi=data.xi)
    if isinstance(data, LinearThetaFractionLossData) or t == "linear_theta_fraction":
        return LinearThetaFractionLoss(k=data.k, xi0=data.xi0)
    if isinstance(data, ConstantHeadLossData) or t == "constant_head":
        return ConstantHeadLoss(zeta=data.zeta, area=area)
    if isinstance(data, LinearThetaHeadLossData) or t == "linear_theta_head":
        return LinearThetaHeadLoss(k=data.k, zeta0=data.zeta0, area=area)
    return None


def map_surface_model(data):
    """Maps UI SurfaceModelData to combaero.network models."""
    if isinstance(data, SmoothModelData) or data.type == "smooth":
        return SmoothModel()
    if isinstance(data, RibbedModelData) or data.type == "ribbed":
        # getattr with the schema's own defaults: a network saved before a
        # field existed should load with that field defaulted, not raise an
        # AttributeError the user cannot act on.
        d = RibbedModelData()
        return RibbedModel(
            e_D=getattr(data, "e_D", d.e_D),
            p_e=getattr(data, "p_e", d.p_e),
            alpha_deg=getattr(data, "alpha_deg", d.alpha_deg),
            W_H=getattr(data, "W_H", d.W_H),
            n_ribbed_walls=getattr(data, "n_ribbed_walls", d.n_ribbed_walls),
            smooth_wall_Nu_multiplier=getattr(
                data, "smooth_wall_Nu_multiplier", d.smooth_wall_Nu_multiplier
            ),
        )
    if isinstance(data, ImpingementModelData) or data.type == "impingement":
        d = ImpingementModelData()
        return ImpingementModel(
            d_jet=getattr(data, "d_jet", d.d_jet),
            xn_d=getattr(data, "xn_d", d.xn_d),
            yn_d=getattr(data, "yn_d", d.yn_d),
            z_d=getattr(data, "z_d", d.z_d),
            row=getattr(data, "row", d.row),
            C_D=getattr(data, "C_D", d.C_D),
        )
    if isinstance(data, SingleJetImpingementModelData) or data.type == "single_jet_impingement":
        import combaero as cb

        d = SingleJetImpingementModelData()
        bc_str = getattr(data, "bc", d.bc)
        bc = (
            cb.ImpingementThermalBC.ConstantWallTemperature
            if bc_str == "constant_wall_temperature"
            else cb.ImpingementThermalBC.ConstantHeatFlux
        )
        return SingleJetImpingementModel(
            bc=bc,
            d_jet=getattr(data, "d_jet", d.d_jet),
            L_D=getattr(data, "L_D", d.L_D),
            R_D=getattr(data, "R_D", d.R_D),
        )
    if isinstance(data, PinFinModelData) or data.type == "pin_fin":
        import combaero as cb

        d = PinFinModelData()
        D = getattr(data, "pin_diameter", d.pin_diameter)
        H = getattr(data, "channel_height", d.channel_height)
        return PinFinModel(
            pin_diameter=D,
            H_D=H / D if D else d.channel_height / d.pin_diameter,
            S_D=getattr(data, "S_D", d.S_D),
            X_D=getattr(data, "X_D", d.X_D),
            N_rows=getattr(data, "N_rows", d.N_rows),
            arrangement=(
                cb.PinArrangement.Staggered
                if getattr(data, "is_staggered", d.is_staggered)
                else cb.PinArrangement.Inline
            ),
            k_pin=getattr(data, "k_pin", d.k_pin),
        )

    # Dimpled was removed in 0.7.0. Reject rather than silently
    # falling back to smooth: a saved network asking for one of these and
    # getting an unenhanced surface back is a wrong answer presented as a
    # right one.
    surface_type = getattr(data, "type", "unknown")
    raise ValueError(
        f"Surface type {surface_type!r} was removed in 0.7.0: its correlation "
        "could not be traced to its cited source. See issue #339. Use a "
        "smooth surface."
    )


def resolve_composition(
    comp: CompositionData, T: float, P: float, solver_settings=None
) -> list[float]:
    """
    Converts UI CompositionData into a mass fraction list (Y)
    using combaero.species utilities.
    """
    import combaero as cb

    if comp.source == "dry_air":
        Y = cb.species.dry_air_mass()
    elif comp.source == "humid_air":
        s = solver_settings
        ambient_T = (
            comp.ambient_T
            if comp.ambient_T is not None
            else (s.ambient_T if s is not None and s.ambient_T is not None else 288.15)
        )
        ambient_P = (
            comp.ambient_P
            if comp.ambient_P is not None
            else (s.ambient_P if s is not None and s.ambient_P is not None else 101325.0)
        )
        ambient_RH = (
            comp.relative_humidity
            if comp.relative_humidity is not None
            else (s.ambient_RH if s is not None and s.ambient_RH is not None else 0.6)
        )
        Y = cb.species.humid_air_mass(ambient_T, ambient_P, ambient_RH)
    elif comp.source == "fuel":
        # Default fuel to CH4 for now, or use a specific fuel if provided
        Y = cb.species.pure_species("CH4")
        Y = cb.species.to_mass(Y)
    elif comp.source == "custom" and comp.custom_fractions:
        vec = cb.species.from_mapping(comp.custom_fractions)
        Y = cb.species.to_mass(vec) if comp.mode == "mole" else vec
    else:
        Y = cb.species.dry_air_mass()

    return [float(y) for y in Y]


def _expand_initial_guess(short_guess: dict, prefix: str) -> dict:
    """Translate short GUI guess keys (``P``, ``T``, ``m_dot``) into the fully
    qualified solver unknown names (e.g. ``"plenum1.P"``). ``P`` and ``T`` are
    broadcast to both static and total variants because the solver exposes
    them as separate unknowns on some node types; the solver silently ignores
    entries that do not correspond to an actual unknown.
    Any key already containing a ``.`` is passed through unchanged so that
    hand-authored qualified guesses still work. Explicit short keys
    (e.g. ``Pt``) take precedence over broadcast values from ``P``/``T``.
    """
    if not short_guess:
        return {}
    result: dict = {}
    # Pass 1 — broadcast P and T to both static and total variants.
    for key, val in short_guess.items():
        if "." in key:
            continue
        if key == "P":
            result[f"{prefix}.P"] = val
            result[f"{prefix}.Pt"] = val
        elif key == "T":
            result[f"{prefix}.T"] = val
            result[f"{prefix}.Tt"] = val
    # Pass 2 — explicit short keys (Pt, Tt, m_dot, X[i], etc.)
    # are applied after the broadcast so they win on conflict.
    for key, val in short_guess.items():
        if "." in key:
            result[key] = val
            continue
        if key in ("P", "T"):
            continue  # already handled in pass 1
        result[f"{prefix}.{key}"] = val
    return result


def _guess_from_prior_result(node_data: dict, prefix: str) -> dict:
    """Build an initial-guess dict from a prior solve result stored in node_data.

    The GUI frontend writes the last solve result into node_data["result"].
    When data.initial_guess is empty (no explicit user override), this
    provides a warm-start from the previously found state so that repeated
    solves on a difficult network benefit from accumulated Newton progress.
    Works for both pressure-node results (state.P/T/Pt/Tt) and element
    results (m_dot stored at the top level).

    Only CONVERGED prior results are used: the best iterate of a failed
    solve is stored too, and seeding from it parks Newton at the stall
    point -- one failed solve (e.g. a timeout) would otherwise poison
    every subsequent run of a network that converges fine from cold.
    """
    prior = node_data.get("result")
    if not prior or not prior.get("success"):
        return {}
    short: dict[str, float] = {}
    state = prior.get("state")
    if state:
        if (p := state.get("P")) is not None:
            short["P"] = float(p)
        if (t := state.get("T")) is not None:
            short["T"] = float(t)
        if (pt := state.get("Pt")) is not None:
            short["Pt"] = float(pt)
        if (tt := state.get("Tt")) is not None:
            short["Tt"] = float(tt)
    if (m_dot := prior.get("m_dot")) is not None:
        short["m_dot"] = float(m_dot)
    return _expand_initial_guess(short, prefix)


def _build_discrete_loss_correlation(
    corr_type: str,
    xi: float,
    k: float,
    xi0: float,
    zeta: float,
    zeta0: float,
    area: float,
    m_ref: float | None = None,
):
    """Build a pressure loss correlation object from discrete loss params."""
    if corr_type == "constant_fraction":
        return ConstantFractionLoss(xi=xi, m_ref=m_ref)
    if corr_type == "linear_theta_fraction":
        return LinearThetaFractionLoss(k=k, xi0=xi0, m_ref=m_ref)
    if corr_type == "constant_head":
        return ConstantHeadLoss(zeta=zeta, area=area)
    if corr_type == "linear_theta_head":
        return LinearThetaHeadLoss(k=k, zeta0=zeta0, area=area)
    return ConstantFractionLoss(xi=xi, m_ref=m_ref)


def _channel_D_at_tee_port(
    schema_edges,
    schema_nodes,
    tee_id: str,
    port: str,
) -> float | None:
    """Return the inner diameter [m] of a channel directly connected to a tee port.

    The user-facing schema connects channels straight to tee handles (the
    port-MCN is auto-inserted at network-build time, not at schema level).
    Used by ``mpce_tee`` dispatch to inherit ``F_C`` / ``F_branch``
    from the connected channel geometry when the
    user has not explicitly set them. Returns None when no channel is
    directly connected (e.g., plenum, boundary, multi-hop).
    """
    hit = _channel_at_tee_port(schema_edges, schema_nodes, tee_id, port)
    return None if hit is None else hit[1]


def _channel_at_tee_port(
    schema_edges,
    schema_nodes,
    tee_id: str,
    port: str,
):
    """The edge joining a tee port straight to a channel, and that channel's D [m].

    None when no channel with a positive D is directly connected.
    """
    src_handle = f"port-{port}-source"
    tgt_handle = f"port-{port}-target"
    for edge in schema_edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        neighbour_id: str | None = None
        if edge.source == tee_id and edge.sourceHandle == src_handle:
            neighbour_id = edge.target
        elif edge.target == tee_id and edge.targetHandle == tgt_handle:
            neighbour_id = edge.source
        if neighbour_id is None:
            continue
        neighbour = next((n for n in schema_nodes if n.id == neighbour_id), None)
        if neighbour is None or neighbour.type != "channel":
            continue
        d = neighbour.data.get("D") if isinstance(neighbour.data, dict) else None
        if d is None:
            continue
        try:
            d_val = float(d)
        except (TypeError, ValueError):
            continue
        if d_val > 0.0:
            return edge, d_val
    return None


def _channel_id_at_tee_port(edge, tee_id: str) -> str:
    """The channel at the far end of a tee-port edge."""
    return edge.target if edge.source == tee_id else edge.source


#: Common-arm area [m^2] of a tee with no channel on its common or straight
#: port and no user F_C (reported as "default"; the GUI mirrors it).
DEFAULT_TEE_F_C = 0.01


def _resolve_tee_areas_with_sources(
    schema_edges,
    schema_nodes,
    tee_id: str,
    F_C_explicit: float | None,
    F_branch_explicit: float | None,
    psi_fallback: float,
) -> tuple[float, float, float, str, str]:
    """Resolve (F_C, F_branch, psi) for a 3-port tee, with each area's source.

    Precedence:
      1. Explicit user-set value (non-null) wins.            -> "user"
      2. Channel D at the connected port (auto-inherited).   -> "channel '<id>'"
      3. Fallback: F_C=0.01 m^2                              -> "default";
         F_branch derived from psi_fallback                  -> "derived from F_C / psi"

    The tee's areas are resolved HERE, in the builder, and handed to the
    element as explicit port areas -- so the element itself would report them
    all as "user". The sources returned here are the real provenance, attached
    to the element for ``resolve_preview`` (#559).

    Returns (F_C, F_branch, psi, F_C_source, F_branch_source).
    """
    F_C = F_C_explicit
    F_C_source = "user" if F_C is not None else ""
    if F_C is None:
        # Try common port first, then straight (both are F_C in Bassett/Idelchik
        # geometry where F_s = F_c by definition).
        for port in ("common", "straight"):
            hit = _channel_at_tee_port(schema_edges, schema_nodes, tee_id, port)
            if hit is not None:
                F_C = _math.pi / 4.0 * hit[1] * hit[1]
                F_C_source = f"channel '{_channel_id_at_tee_port(hit[0], tee_id)}'"
                break
    if F_C is None:
        F_C = DEFAULT_TEE_F_C
        F_C_source = "default"

    F_branch = F_branch_explicit
    F_branch_source = "user" if F_branch is not None else ""
    if F_branch is None:
        hit = _channel_at_tee_port(schema_edges, schema_nodes, tee_id, "branch")
        if hit is not None:
            F_branch = _math.pi / 4.0 * hit[1] * hit[1]
            F_branch_source = f"channel '{_channel_id_at_tee_port(hit[0], tee_id)}'"

    if F_branch is not None and F_branch > 0.0:
        psi = F_C / F_branch
    else:
        # No branch channel; use stored psi to derive F_branch
        psi = max(float(psi_fallback), 1e-12)
        F_branch = F_C / psi
        F_branch_source = "derived from F_C / psi"

    return F_C, F_branch, psi, F_C_source, F_branch_source


def _resolve_tee_areas(
    schema_edges,
    schema_nodes,
    tee_id: str,
    F_C_explicit: float | None,
    F_branch_explicit: float | None,
    psi_fallback: float,
) -> tuple[float, float, float]:
    """Resolve (F_C, F_branch, psi) for a 3-port tee, applying inheritance.

    See ``_resolve_tee_areas_with_sources`` for the precedence; this returns
    just the three resolved floats.
    """
    F_C, F_branch, psi, _, _ = _resolve_tee_areas_with_sources(
        schema_edges, schema_nodes, tee_id, F_C_explicit, F_branch_explicit, psi_fallback
    )
    return F_C, F_branch, psi


def _port_from_handle(handle: str | None) -> str | None:
    """Extract port name from a tee handle string.

    'port-common-source' -> 'common', 'port-branch-target' -> 'branch', etc.
    Returns None for non-tee handles.
    """
    if (
        handle
        and handle.startswith("port-")
        and (handle.endswith("-source") or handle.endswith("-target"))
    ):
        return handle.split("-")[1]
    return None


def _find_junction_port(
    edges: list,
    junction_id: str,
    port: str,
    nodes_map: dict,
    tee_port_node_map: dict | None = None,
    wall_edge_remap: dict | None = None,
    all_ports: tuple[str, ...] = ("common", "straight", "branch"),
) -> str:
    """Return the physical-node ID connected to a named junction port.

    Shared by the tee (ports common/straight/branch) and the ejector (ports
    primary/secondary/outlet); ``all_ports`` only shapes the error message.
    Checks both edge directions: junction as source (sourceHandle ==
    port-X-source) and junction as target (targetHandle == port-X-target).
    Falls back to *tee_port_node_map* for auto-created junction nodes when an
    element is connected directly to a port without an explicit plenum.
    """
    src_handle = f"port-{port}-source"
    tgt_handle = f"port-{port}-target"
    for edge in edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        if edge.source == junction_id and edge.sourceHandle == src_handle:
            nid = edge.target
        elif edge.target == junction_id and edge.targetHandle == tgt_handle:
            nid = edge.source
        else:
            continue
        # An auto-inserted port node wins: it is also created for a plenum or
        # boundary wired straight to the port (build pass 1b).
        if tee_port_node_map and (junction_id, port) in tee_port_node_map:
            return tee_port_node_map[(junction_id, port)]
        if nid in nodes_map:
            if wall_edge_remap and edge.id in wall_edge_remap:
                return wall_edge_remap[edge.id]
            return nid
        raise ValueError(
            f"Junction '{junction_id}' port '{port}' connects to '{nid}', which is "
            "not a physical node (plenum or boundary). All junction ports must "
            "connect directly to plena or boundaries."
        )
    raise ValueError(
        f"Junction '{junction_id}' port '{port}' is not connected. All ports "
        f"({', '.join(all_ports)}) must be wired to a node."
    )


def _trace_boundary(
    port_id: str, junction_id: str, edges: list, nodes_map: dict, max_hops: int = 8
):
    """BFS outward from a junction port (never through the junction itself) to
    the first boundary node reached. Returns the boundary node object (a
    PressureBoundary or MassFlowBoundary) or None. Used only to derive a
    warm-start; a miss just means no seed for that port.
    """
    seen = {port_id, junction_id}
    frontier = [port_id]
    # An auto-inserted port MCN (a channel/element on the port) is a synthetic
    # node with no canvas edges: start from whatever is wired to the junction's
    # handle for that port instead, or the outlet of e.g. a channel-fed ejector
    # gets no seed at all (and the stiff cold solve stalls).
    if port_id.startswith("__tee_jct__") and port_id.rsplit("_", 1)[-1]:
        role = port_id.rsplit("_", 1)[-1]
        for edge in edges:
            if edge.source == junction_id:
                handle = edge.sourceHandle
            elif edge.target == junction_id:
                handle = edge.targetHandle
            else:
                handle = None
            if handle and f"port-{role}-" in handle:
                frontier.append(edge.target if edge.source == junction_id else edge.source)
        for nid in frontier[1:]:
            node = nodes_map.get(nid)
            if node is not None and node.is_flow_boundary:
                return node
            seen.add(nid)
    for _ in range(max_hops):
        nxt: list[str] = []
        for nid in frontier:
            for edge in edges:
                if edge.data and edge.data.get("type") == "thermal":
                    continue
                if edge.source == nid and edge.target != junction_id:
                    other = edge.target
                elif edge.target == nid and edge.source != junction_id:
                    other = edge.source
                else:
                    continue
                if other in seen:
                    continue
                seen.add(other)
                node = nodes_map.get(other)
                if node is not None and node.is_flow_boundary:
                    return node
                nxt.append(other)
        frontier = nxt
    return None


def _seed_ejector_warmstart(
    net: FlowNetwork,
    ejector: EjectorElement,
    port_ids: tuple[str, str, str],
    edges: list,
    nodes_map: dict,
) -> None:
    """Warm-start a GUI-built ejector so the (stiff, double-choked) solve
    converges without the user hand-seeding interior nodes.

    The port momentum-chambers' stagnation pressures start far from their
    feeding-boundary values, which puts the cold Newton start outside the
    ejector's narrow convergence basin (confirmed: no init_strategy cracks it
    cold). Here we trace each port to its boundary, run the ejector's own
    closed-form physics to get a consistent operating point, and seed the port
    MCN P/Pt plus the element's P_jct (== P_c*). Only seeds values the user did
    not already provide; best-effort (any failure leaves the seed unset).
    """
    import combaero as _cb
    from combaero.network._ejector_huang1999 import (
        ETA_P,
        EjectorGeometry,
        choked_mass_flow,
        critical_back_pressure,
        entrainment_ratio,
    )
    from combaero.network.ejector_element import DEFAULT_EJECTOR_MIXING_AREA, choke_plane_gamma

    primary_id, secondary_id, outlet_id = port_ids
    pb = _trace_boundary(primary_id, ejector.id, edges, nodes_map)
    sb = _trace_boundary(secondary_id, ejector.id, edges, nodes_map)
    ob = _trace_boundary(outlet_id, ejector.id, edges, nodes_map)

    # Nominal air properties for the SEED only (exact gamma/R not needed here).
    X = _cb.species.dry_air()
    t_g = float(getattr(pb, "Tt", 440.0))
    t_e = float(getattr(sb, "Tt", 280.0))
    try:
        gamma = choke_plane_gamma(t_e, X)
    except Exception:
        gamma = 1.4
    r_gas = float(_cb.specific_gas_constant(X))

    # Primary stagnation: direct for a pressure boundary; inverse choked-flow
    # for a mass boundary (mdot is linear in p0).
    p_g = None
    if isinstance(pb, PressureBoundary):
        p_g = float(pb.Pt)
    elif isinstance(pb, MassFlowBoundary):
        unit = choked_mass_flow(1.0, t_g, ejector.throat_area, gamma, r_gas, eta=ETA_P)
        if unit > 0.0:
            p_g = float(pb.m_dot) / unit
    p_e = float(sb.Pt) if isinstance(sb, PressureBoundary) else None
    p_o = float(ob.Pt) if isinstance(ob, PressureBoundary) else None

    # Regime discriminator for the SEED. The choked-inverse p_g above is the
    # primary stagnation the mass flow WOULD need to reach sonic at the throat.
    # If that is below the outlet back pressure, the primary cannot choke (the
    # forced mass flow is below choked_mass_flow(p_back)); the physical branch
    # is the subsonic jet pump, where the primary Pt floats to ~the back
    # pressure and the mixing pressure P_py sits just below it. Seeding the
    # choked-inverse value there parks Newton in the (spurious) double-choked
    # basin with Pt_primary < Pb, from which it stalls to nan. Detect that and
    # seed the jet-pump operating point instead.
    jetpump_seed = (
        isinstance(pb, MassFlowBoundary) and p_g is not None and p_o is not None and p_g < p_o
    )
    if jetpump_seed:
        p_g = p_o  # primary floats up to the back pressure (dp ~ 0)

    def _seed_node(node_id: str, value: float | None) -> None:
        if value is None or value <= 0.0:
            return
        node = net.nodes.get(node_id)
        if node is None:
            return
        existing = dict(getattr(node, "initial_guess", {}) or {})
        if any(k.endswith(".P") or k.endswith(".Pt") for k in existing):
            return  # user already seeded this port
        existing[f"{node_id}.P"] = float(value)
        existing[f"{node_id}.Pt"] = float(value)
        node.initial_guess = existing

    _seed_node(primary_id, p_g)
    _seed_node(secondary_id, p_e)
    _seed_node(outlet_id, p_o)

    # Geometry for the closed forms. An unset mixing area resolves (later, in
    # resolve_topology) from the outlet pipe, else the default; use the same
    # value here, or every seed below silently fails on area_ratio_mix.
    a3 = ejector.mixing_area
    if a3 is None:
        a3 = ejector._outlet_pipe_area(net)[0] or DEFAULT_EJECTOR_MIXING_AREA
    geom = EjectorGeometry(
        ejector.nozzle_exit_area / ejector.throat_area, float(a3) / ejector.throat_area
    )

    # Seed the owned unknown P_jct (repurposed as the mixing-plane pressure
    # P_py). In the jet-pump regime it sits just below the back pressure; in
    # the critical regime the diagnostic critical back pressure P_c* is the
    # right neighbourhood.
    p_jct_seed: float | None = None
    if jetpump_seed and p_o is not None:
        p_jct_seed = 0.99 * p_o
    elif p_g and p_e:
        try:
            p_jct_seed = float(
                critical_back_pressure(
                    p_g, t_g, p_e, t_e, geom, gamma, r_gas, ejector.recovery_efficiency
                ).p_c_pa
            )
        except Exception:
            p_jct_seed = None
    if p_jct_seed is not None and p_jct_seed > 0.0:
        guess = dict(getattr(ejector, "initial_guess", {}) or {})
        if not any(k.endswith(".P_jct") for k in guess):
            guess[f"{ejector.id}.P_jct"] = p_jct_seed
            ejector.initial_guess = guess

    # Seed the flows on each port's element chain (#567). Left to the generic
    # defaults (tens of kg/s for a ~0.1 kg/s ejector) a channel on a port
    # sends the cold Newton start to Mach-20 port states it never returns
    # from; seeded with the closed-form operating point it converges.
    mp: float | None = None
    if isinstance(pb, MassFlowBoundary):
        mp = float(pb.m_dot)
    elif p_g:
        mp = float(choked_mass_flow(p_g, t_g, ejector.throat_area, gamma, r_gas, eta=ETA_P))
    # Critical regime only. In the jet-pump regime the seed sits at p_g = p_o
    # (no drive), where the closed form returns omega ~ 0 and a zero-entrainment
    # seed parks Newton at that state; the solver's own mass-flow propagation
    # from the primary's mass boundary already converges those cases.
    if jetpump_seed:
        return
    omega: float | None = None
    if mp and p_g and p_e:
        try:
            omega = float(entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).omega)
        except Exception:
            omega = None
    if omega is not None and not (0.0 < omega < 50.0):
        omega = None  # an out-of-envelope closed form is no seed
    if mp:
        _seed_port_chain_flows(net, ejector, primary_id, mp, inflow=True)
        if omega is not None:
            _seed_port_chain_flows(net, ejector, secondary_id, omega * mp, inflow=True)
            _seed_port_chain_flows(net, ejector, outlet_id, (1.0 + omega) * mp, inflow=False)


def _seed_port_chain_flows(
    net: FlowNetwork, junction: object, port_id: str, flow: float, inflow: bool
) -> None:
    """Seed ``<id>.m_dot`` on the linear chain of two-port elements running
    outward from a junction port to its boundary, signed by each element's
    orientation (flow towards the port on an inlet, away on an outlet). Stops
    at a branch, a multi-port element or a boundary; never overrides a guess
    the user gave."""
    node_id, prev = port_id, junction
    for _ in range(32):
        around = net.get_upstream_elements(node_id) + net.get_downstream_elements(node_id)
        elems = [e for e in around if e is not prev and e is not junction]
        if len(elems) != 1:
            return
        elem = elems[0]
        frm, to = getattr(elem, "from_node", None), getattr(elem, "to_node", None)
        if node_id not in (frm, to) or frm is None or to is None:
            return
        # Walking outward from the port: an element pointing back at the port
        # carries +flow on an inlet chain, one pointing away carries -flow.
        points_to_port = to == node_id
        signed = flow if points_to_port == inflow else -flow
        unk = f"{elem.id}.m_dot"
        if unk in elem.unknowns():
            guess = dict(getattr(elem, "initial_guess", {}) or {})
            if unk not in guess:
                guess[unk] = signed
                elem.initial_guess = guess
        node_id = frm if to == node_id else to
        node = net.nodes.get(node_id)
        if node is None or getattr(node, "is_flow_boundary", False):
            return
        prev = elem


def build_network_from_schema(schema: NetworkGraphSchema) -> FlowNetwork:
    net = FlowNetwork()
    nodes_map = {}  # ID -> Physical Node
    plenum_node_ids: set[str] = set()
    element_nodes = []  # List of node schemas that are actually elements

    # 1. First Pass: Create Physical Nodes
    for node_schema in schema.nodes:
        node_id = node_schema.id
        node_data = node_schema.data
        node_type = node_schema.type

        if node_type == "plenum":
            data = PlenumData(**node_data)
            node = PlenumNode(node_id)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
            plenum_node_ids.add(node_id)
        elif node_type == "mass_boundary":
            data = MassBoundaryData(**node_data)
            Y = resolve_composition(data.composition, data.Tt, 101325.0, schema.solver_settings)
            node = MassFlowBoundary(node_id, m_dot=data.m_dot, Tt=data.Tt, Y=Y)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "pressure_boundary":
            data = PressureBoundaryData(**node_data)
            Y = resolve_composition(data.composition, data.Tt, data.Pt, schema.solver_settings)
            node = PressureBoundary(node_id, Pt=data.Pt, Tt=data.Tt, Y=Y)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "combustor":
            data = CombustorData(**node_data)
            node = CombustorNode(
                node_id,
                method=data.method,
                area=data.area,
                Dh=data.Dh,
                surface=ConvectiveSurface(
                    area=data.area or 0.0,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier
                    * (schema.solver_settings.Nu_multiplier or 1.0),
                    f_multiplier=data.f_multiplier * (schema.solver_settings.f_multiplier or 1.0),
                ),
            )
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "momentum_chamber":
            data = MomentumChamberData(**node_data)
            node = MomentumChamberNode(
                node_id,
                area=data.area,
                Dh=data.Dh,
                surface=ConvectiveSurface(
                    area=data.area or 0.0,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier
                    * (schema.solver_settings.Nu_multiplier or 1.0),
                    f_multiplier=data.f_multiplier * (schema.solver_settings.f_multiplier or 1.0),
                ),
            )
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            if not node.initial_guess:
                node.initial_guess = _guess_from_prior_result(node_data, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "probe":
            pass  # Display-only diagnostic nodes -- not part of the solver graph
        elif node_type == "wall":
            data = WallData(**node_data)
            node = WallNode(node_id)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        else:
            # This is an element node (Channel, Orifice, etc.)
            element_nodes.append(node_schema)

    element_ids = {node.id for node in element_nodes}
    # Impingement arrays (#465) expand into a row chain. Their ends fan out to
    # N plates (supply) or collect the last crossflow segment (exit), which a
    # MomentumChamberNode junction cannot represent, so element-to-array links
    # get a PlenumNode instead -- the jet supply plenum the source assumes.
    array_ids = {node.id for node in element_nodes if node.type == "impingement_array"}
    gui_arrays: dict[str, ImpingementArray] = {}
    # Junction elements with named handle ports (common/straight/branch for
    # the tee; primary/secondary/outlet for the ejector). Both get the same
    # per-port MomentumChamberNode auto-insertion below.
    junction_ids = {
        node.id for node in element_nodes if node.type in ("mpce_tee", "ejector", "effusion_liner")
    }
    # An effusion liner's ports are a reservoir, a dump and a discharge that
    # N panels merge into: plenum semantics. A momentum chamber there would
    # refuse the N inflows and mis-state the reservoir (#471).
    liner_ids = {node.id for node in element_nodes if node.type == "effusion_liner"}

    # 1b. Create implicit junction nodes for element -> element visual links.
    # For direct element<->element connections we auto-insert a MomentumChamberNode.
    # For direct element<->tee-port connections we auto-insert a PlenumNode so that
    # users don't need to manually place a plenum on every tee arm.
    edge_junction_map: dict[tuple[str, str], str] = {}
    _ejector_seeds: list[tuple[EjectorElement, tuple[str, str, str]]] = []
    tee_port_node_map: dict[tuple[str, str], str] = {}  # (junction_id, port) -> node_id
    # (from_node, to_node) lossless links joining a directly-wired plenum or
    # boundary to the port node inserted for it.
    port_node_links: list[tuple[str, str]] = []
    for edge in schema.edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        src_elem = edge.source in element_ids
        tgt_elem = edge.target in element_ids
        src_tee = edge.source in junction_ids
        tgt_tee = edge.target in junction_ids

        if src_elem and tgt_elem and not src_tee and not tgt_tee:
            # Standard element-to-element: auto-insert MomentumChamberNode
            junction_id = f"__junction__{edge.source}__{edge.target}"
            if junction_id not in nodes_map:
                to_array = edge.source in array_ids or edge.target in array_ids
                junction_node = (
                    PlenumNode(junction_id) if to_array else MomentumChamberNode(junction_id)
                )
                net.add_node(junction_node)
                nodes_map[junction_id] = junction_node
            edge_junction_map[(edge.source, edge.target)] = junction_id

        elif src_tee and tgt_elem:
            # Junction port → element (e.g. tee.common-source → channel): auto-insert
            # MomentumChamberNode so that Pt != P_static and velocity is accounted for.
            port = _port_from_handle(edge.sourceHandle)
            if port and (edge.source, port) not in tee_port_node_map:
                jid = f"__tee_jct__{edge.source}_{port}"
                if jid not in nodes_map:
                    jn = PlenumNode(jid) if edge.source in liner_ids else MomentumChamberNode(jid)
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(edge.source, port)] = jid
                edge_junction_map[(edge.source, edge.target)] = jid

        elif src_elem and tgt_tee:
            # Element → junction port (e.g. channel → tee.straight-target): auto-insert MomentumChamberNode
            port = _port_from_handle(edge.targetHandle)
            if port and (edge.target, port) not in tee_port_node_map:
                jid = f"__tee_jct__{edge.target}_{port}"
                if jid not in nodes_map:
                    jn = PlenumNode(jid) if edge.target in liner_ids else MomentumChamberNode(jid)
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(edge.target, port)] = jid
                edge_junction_map[(edge.source, edge.target)] = jid

        elif (src_tee and edge.target in nodes_map) or (tgt_tee and edge.source in nodes_map):
            # Junction port <-> plenum / boundary wired directly. The junction
            # needs a MomentumChamberNode on each port (MultiPortChamberBase
            # refuses anything else), so insert one and join it to the node with
            # a lossless link, as a channel on the port would get (#559 audit).
            # A momentum chamber the user placed is used as-is; a liner's ports
            # are plenums by design.
            jct, other = (edge.source, edge.target) if src_tee else (edge.target, edge.source)
            handle = edge.sourceHandle if src_tee else edge.targetHandle
            port = _port_from_handle(handle)
            if (
                port
                and jct not in liner_ids
                and not isinstance(nodes_map[other], MomentumChamberNode)
                and (jct, port) not in tee_port_node_map
            ):
                jid = f"__tee_jct__{jct}_{port}"
                if jid not in nodes_map:
                    jn = MomentumChamberNode(jid)
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(jct, port)] = jid
                # Flow direction follows the canvas edge: node -> port (inlet)
                # or port -> node (outlet).
                port_node_links.append((other, jid) if tgt_tee else (jid, other))

    # 1c. A tee port whose channel has a different area meets it through a
    # sharp step (#536). The tee's ports keep the closures' geometry (F_S =
    # F_C), so a straight channel narrower than the main pipe used to join a
    # port of the main area with nothing in between: the ideal tee, which
    # couples on static pressure, then raised the channel's Pt by its dynamic
    # head, and Mynard, which couples on Pt, contracted losslessly. Bolted
    # together, two bores of different size are a sharp contraction or
    # expansion, so that is what is inserted: a hidden AreaChangeElement
    # (Weisbach / Borda-Carnot) and a momentum chamber on the channel side.
    tee_steps: dict[str, tuple[str, str, float]] = {}
    for tee_node in element_nodes:
        if tee_node.type != "mpce_tee":
            continue
        _td = MPCETeeData(**tee_node.data)
        _F_C, _F_B, _ = _resolve_tee_areas(
            schema.edges, schema.nodes, tee_node.id, _td.F_C, _td.F_branch, _td.psi
        )
        for port, A_port in (("common", _F_C), ("straight", _F_C), ("branch", _F_B)):
            hit = _channel_at_tee_port(schema.edges, schema.nodes, tee_node.id, port)
            if hit is None:
                continue
            edge, D = hit
            A_ch = _math.pi / 4.0 * D * D
            if abs(A_ch / A_port - 1.0) <= TEE_PORT_STEP_TOL:
                continue
            key = (edge.source, edge.target)
            port_node = edge_junction_map.get(key)
            if port_node is None:
                continue
            step_id = f"__tee_step__{tee_node.id}_{port}"
            side_id = f"{step_id}__node"
            side = MomentumChamberNode(side_id)
            net.add_node(side)
            nodes_map[side_id] = side
            if edge.source == tee_node.id:  # drawn tee -> channel
                step = AreaChangeElement(
                    step_id, port_node, side_id, F0=A_port, F1=A_ch, model_type="sharp"
                )
            else:
                step = AreaChangeElement(
                    step_id, side_id, port_node, F0=A_ch, F1=A_port, model_type="sharp"
                )
            net.add_element(step)
            edge_junction_map[key] = side_id
            tee_steps[step_id] = (tee_node.id, port, A_ch / A_port)
    net._gui_tee_steps = tee_steps

    # Pre-scan: for each edge whose target is a WallNode, the 2nd+ connection
    # gets its own clone so the solver enforces m_dot=0 per branch individually
    # (not just sum(m_dot)=0, which allows non-zero flows to cancel).
    # This covers both regular elements and tee junction port connections.
    _wall_edge_remap: dict[str, str] = {}  # edge_id -> effective target node_id
    _wall_target_count: dict[str, int] = {}
    for _edge in schema.edges:
        if _edge.data and _edge.data.get("type") == "thermal":
            continue
        _tgt = _edge.target
        if isinstance(nodes_map.get(_tgt), WallNode):
            _wall_target_count[_tgt] = _wall_target_count.get(_tgt, 0) + 1
            if _wall_target_count[_tgt] > 1:
                _n = _wall_target_count[_tgt]
                _clone_id = f"{_tgt}_w{_n}"
                _clone = WallNode(_clone_id)
                net.add_node(_clone)
                nodes_map[_clone_id] = _clone
                _wall_edge_remap[_edge.id] = _clone_id

    # 2. Second Pass: Create Elements
    for elem_node_schema in element_nodes:
        elem_id = elem_node_schema.id
        elem_data = elem_node_schema.data
        elem_type = elem_node_schema.type

        # 3-port element (Mynard-based). The flow_direction field constrains
        # which side is upstream.
        if elem_type == "mpce_tee":
            _md = MPCETeeData(**elem_data)
            # Resolve areas with inheritance precedence (see _resolve_tee_areas).
            _F_C, _A_branch, _psi, _F_C_src, _A_branch_src = _resolve_tee_areas_with_sources(
                schema.edges,
                schema.nodes,
                elem_id,
                F_C_explicit=_md.F_C,
                F_branch_explicit=_md.F_branch,
                psi_fallback=_md.psi,
            )

            _common = _find_junction_port(
                schema.edges, elem_id, "common", nodes_map, tee_port_node_map, _wall_edge_remap
            )
            _straight = _find_junction_port(
                schema.edges, elem_id, "straight", nodes_map, tee_port_node_map, _wall_edge_remap
            )
            _branch = _find_junction_port(
                schema.edges, elem_id, "branch", nodes_map, tee_port_node_map, _wall_edge_remap
            )

            # Branch (separating): single inlet at common, two outlets.
            # Merge (joining): two inlets at straight+branch, single outlet at common.
            if _md.flow_direction == "branch":
                _inlets = [_common]
                _outlets = [_straight, _branch]
                _inlet_angles = [0.0]
                _outlet_angles = [0.0, _md.theta_deg]
                _port_areas = [_F_C, _F_C, _A_branch]
            else:  # "merge"
                _inlets = [_straight, _branch]
                _outlets = [_common]
                _inlet_angles = [0.0, _md.theta_deg]
                _outlet_angles = [0.0]
                _port_areas = [_F_C, _A_branch, _F_C]

            # strict=False: strict mode raises whenever a Newton ITERATE
            # (not the converged answer) explores a wrong-sign port
            # m_dot, killing GUI solves that would otherwise converge
            # ("Unexpected error during residual evaluation"). Soft mode
            # lets the solver explore; the post-solve
            # verify_solution_consistent guard in NetworkSolver rejects
            # any wrong-direction artifact root at convergence.
            if _md.junction_model == "constant_k":
                # Simplest-model tier: fixed handbook K per non-common
                # port, referenced to the common-leg dynamic head. Port
                # order is inlets first: branch = [common, straight,
                # branch]; merge = [straight, branch, common].
                # A K left blank comes from Bassett's analytical pair at the
                # design split and the resolved geometry.
                _K_str_d, _K_bra_d = ConstantKTeeElement.bassett_design_K(
                    _md.flow_direction, _md.design_branch_share, _psi, _md.theta_deg
                )
                _K_str = _md.K_straight if _md.K_straight is not None else _K_str_d
                _K_bra = _md.K_branch if _md.K_branch is not None else _K_bra_d
                if _md.flow_direction == "branch":
                    _k_ports = {1: _K_str, 2: _K_bra}
                else:
                    _k_ports = {0: _K_str, 1: _K_bra}
                _mpce = ConstantKTeeElement(
                    id=elem_id,
                    inlet_nodes=_inlets,
                    outlet_nodes=_outlets,
                    inlet_angles_deg=_inlet_angles,
                    outlet_angles_deg=_outlet_angles,
                    port_areas=_port_areas,
                    flow_direction=_md.flow_direction,
                    strict=False,
                    K_ports=_k_ports,
                )
            elif _md.junction_model == "ideal_momentum":
                # Parameter-free baseline (#517): no K, no declared
                # direction -- the inlet/outlet split above only wires the
                # ports and sets the angle convention. loss_factor (#526) is
                # the user's knob; 1.0 is the closure itself.
                _mpce = IdealMomentumTeeElement(
                    id=elem_id,
                    inlet_nodes=_inlets,
                    outlet_nodes=_outlets,
                    inlet_angles_deg=_inlet_angles,
                    outlet_angles_deg=_outlet_angles,
                    port_areas=_port_areas,
                    loss_factor=_md.loss_factor,
                )
            else:
                _mpce = MultiPortChamberElement(
                    id=elem_id,
                    inlet_nodes=_inlets,
                    outlet_nodes=_outlets,
                    inlet_angles_deg=_inlet_angles,
                    outlet_angles_deg=_outlet_angles,
                    port_areas=_port_areas,
                    flow_direction=_md.flow_direction,
                    strict=False,
                    joining_etransfer_alpha=_md.joining_etransfer_alpha,
                )
            _mpce.initial_guess = _expand_initial_guess(_md.initial_guess, elem_id)
            if not _mpce.initial_guess:
                _mpce.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            # The element sees these areas as explicit port areas ("user"); the
            # builder knows where they really came from. resolve_preview reports
            # these in place of the per-port sources (#559).
            _mpce._preview_fields = {
                "F_C": (_F_C, _F_C_src),
                "F_branch": (_A_branch, _A_branch_src),
            }
            net.add_element(_mpce)
            continue

        # Ejector: 3-port MultiPortChamberBase with fixed roles (primary +
        # secondary inlets, one outlet). Reuses the same named-port resolution
        # and per-port MomentumChamberNode auto-insertion as the tee.
        if elem_type == "ejector":
            _ed = EjectorData(**elem_data)
            _ports = ("primary", "secondary", "outlet")
            _primary = _find_junction_port(
                schema.edges,
                elem_id,
                "primary",
                nodes_map,
                tee_port_node_map,
                _wall_edge_remap,
                all_ports=_ports,
            )
            _secondary = _find_junction_port(
                schema.edges,
                elem_id,
                "secondary",
                nodes_map,
                tee_port_node_map,
                _wall_edge_remap,
                all_ports=_ports,
            )
            _outlet = _find_junction_port(
                schema.edges,
                elem_id,
                "outlet",
                nodes_map,
                tee_port_node_map,
                _wall_edge_remap,
                all_ports=_ports,
            )
            _ejector = EjectorElement(
                id=elem_id,
                primary_node=_primary,
                secondary_node=_secondary,
                outlet_node=_outlet,
                throat_area=_ed.throat_area,
                nozzle_exit_area=_ed.nozzle_exit_area,
                mixing_area=_ed.mixing_area,
                recovery_efficiency=_ed.recovery_efficiency,
            )
            _ejector.initial_guess = _expand_initial_guess(_ed.initial_guess, elem_id)
            if not _ejector.initial_guess:
                _ejector.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(_ejector)
            # Warm-started once the whole network exists (end of the build):
            # the flow seed walks each port's element chain, and an inherited
            # A_3 needs the outlet pipe, which may be built after the ejector.
            _ejector_seeds.append((_ejector, (_primary, _secondary, _outlet)))
            continue

        # Effusion liner (#471): three named ports, expanded into N stations.
        if elem_type == "effusion_liner":
            _ld = EffusionLinerData(**elem_data)
            _lports = ("coolantin", "coolantout", "discharge")
            _ends = [
                _find_junction_port(
                    schema.edges,
                    elem_id,
                    port,
                    nodes_map,
                    tee_port_node_map,
                    _wall_edge_remap,
                    all_ports=_lports,
                )
                for port in _lports
            ]
            liner = EffusionLiner(
                elem_id,
                n_segments=_ld.n_segments,
                length=_ld.length,
                width=_ld.width,
                duct_height=_ld.duct_height,
                hole_diameter=_ld.hole_diameter,
                wall_thickness=_ld.wall_thickness,
                pitch_x=_ld.pitch_x,
                pitch_y=_ld.pitch_y,
                angle_deg=_ld.angle_deg,
                entry_K=_ld.entry_K,
                roughness=_ld.roughness,
                wall_conductivity=_ld.wall_conductivity,
                gas_film=_ld.gas_film,
                gas_augmentation=_ld.gas_augmentation,
                turbulence_intensity=_ld.turbulence_intensity,
                gas_heat_flux=_ld.gas_heat_flux,
            )
            liner.add_to(net, _ends[0], _ends[1], _ends[2])
            gui_arrays[elem_id] = liner
            continue

        source_id = None
        target_id = None
        incoming_count = 0
        outgoing_count = 0

        for edge in schema.edges:
            if edge.data and edge.data.get("type") == "thermal":
                continue
            if edge.target == elem_id:
                incoming_count += 1
                if incoming_count > 1:
                    raise ValueError(f"Element '{elem_id}' has multiple upstream links.")
                if edge.source in nodes_map:
                    source_id = edge.source
                elif edge.source in element_ids:
                    source_id = edge_junction_map.get((edge.source, elem_id))
            if edge.source == elem_id:
                outgoing_count += 1
                if outgoing_count > 1:
                    raise ValueError(f"Element '{elem_id}' has multiple downstream links.")
                if edge.target in nodes_map:
                    target_id = _wall_edge_remap.get(edge.id, edge.target)
                elif edge.target in element_ids:
                    target_id = edge_junction_map.get((elem_id, edge.target))

        if not source_id or not target_id:
            raise ValueError(f"Element '{elem_id}' is not fully connected.")

        if elem_type == "channel":
            data = ChannelData(**elem_data)
            # The convective area is left to the element: it depends on the
            # surface model (#462) and, when D is inherited, on the resolved
            # diameter. ChannelElement.resolve_topology fills it.
            conv_area = 0.0
            elem = ChannelElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                length=data.L,
                diameter=data.D,
                Dh=data.Dh,
                # Unset D inherits; with nothing to inherit from, a marked
                # default keeps a freshly dropped chain solvable (#570).
                default_diameter=DEFAULT_CHANNEL_DIAMETER,
                roughness=data.roughness,
                friction_model=data.friction_model,
                surface=ConvectiveSurface(
                    area=conv_area,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier
                    * (schema.solver_settings.Nu_multiplier or 1.0),
                    f_multiplier=data.f_multiplier * (schema.solver_settings.f_multiplier or 1.0),
                ),
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "orifice":
            data = OrificeData(**elem_data)
            elem = OrificeElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                diameter=data.diameter,
                Cd=data.Cd,
                correlation=data.correlation,
                plate_thickness=data.plate_thickness,
                edge_radius=data.edge_radius,
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "effusion_plate":
            data = EffusionPlateData(**elem_data)
            elem = EffusionPlateElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                hole_diameter=data.hole_diameter,
                wall_thickness=data.wall_thickness,
                pitch_x=data.pitch_x,
                pitch_y=data.pitch_y,
                panel_length=data.panel_length,
                panel_width=data.panel_width,
                angle_deg=data.angle_deg,
                correlation=data.correlation,
                Cd=data.Cd,
                internal_Nu_multiplier=data.internal_Nu_multiplier,
                wall_conductivity=data.wall_conductivity,
                gas_film=data.gas_film,
                gas_augmentation=data.gas_augmentation,
                turbulence_intensity=data.turbulence_intensity,
                gas_heat_flux=data.gas_heat_flux,
            )
            net.add_element(elem)
        elif elem_type == "impingement_array":
            data = ImpingementArrayData(**elem_data)
            arr = ImpingementArray(
                elem_id,
                n_rows=data.n_rows,
                d_jet=data.d_jet,
                xn_d=data.xn_d,
                yn_d=data.yn_d,
                z_d=data.z_d,
                span=data.span,
                plate_thickness=data.plate_thickness,
                pattern=data.pattern,
                correlation=data.correlation,
                Cd=data.Cd,
                Nu_multiplier=data.Nu_multiplier * (schema.solver_settings.Nu_multiplier or 1.0),
                roughness=data.roughness,
            )
            arr.add_to(net, source_id, target_id)
            gui_arrays[elem_id] = arr
        elif elem_type == "area_change":
            data = AreaChangeData(**elem_data)
            elem = AreaChangeElement(
                id=elem_id,
                from_node=source_id,
                to_node=target_id,
                F0=data.F0,
                F1=data.F1,
                model_type=data.model_type,
                length=data.length,
                D_h=data.D_h,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "lossless_connection":
            elem = LosslessConnectionElement(elem_id, from_node=source_id, to_node=target_id)
            net.add_element(elem)
        elif elem_type == "discrete_loss":
            data = DiscreteLossData(**elem_data)
            # area=None signals "inherit from the topology in resolve_topology".
            # The node-area lookup that used to live here is now part of the
            # library's shared inference path, so GUI and headless networks
            # resolve the same area (issue #262).
            area: float | None = data.area
            correlation = _build_discrete_loss_correlation(
                data.correlation_type,
                data.xi,
                data.k,
                data.xi0,
                data.zeta,
                data.zeta0,
                area if area and area > 0 else 0.1,
                m_ref=data.m_ref if data.m_ref and data.m_ref > 0 else None,
            )
            theta_source = data.theta_source if data.theta_source not in (None, "none") else None
            surface = ConvectiveSurface(
                area=area if area and area > 0 else 0.0,
                model=map_surface_model(data.surface),
                Nu_multiplier=data.Nu_multiplier,
                f_multiplier=data.f_multiplier,
            )
            elem = PressureLossElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                correlation=correlation,
                theta_source=theta_source,
                area=area,
                surface=surface,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)
        elif elem_type == "vortex":
            data = VortexData(**elem_data)
            omega = (
                data.omega_rpm
                if data.omega_rpm is not None
                else (schema.solver_settings.omega_rpm or 0.0)
            )
            elem = VortexElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                r_c=data.r_c,
                r_out=data.r_out,
                r_in=data.r_in,
                omega_rpm=float(omega),
                n=data.n,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            if not elem.initial_guess:
                elem.initial_guess = _guess_from_prior_result(elem_data, elem_id)
            net.add_element(elem)

    # 3. Third Pass: Auto-connections
    for edge in schema.edges:
        if (
            edge.source in nodes_map
            and edge.target in nodes_map
            and not (edge.data or {}).get("type")
        ):
            auto_id = f"__auto_link__{edge.source}__{edge.target}"
            if auto_id not in net.elements:
                elem = LosslessConnectionElement(
                    auto_id, from_node=edge.source, to_node=edge.target
                )
                net.add_element(elem)

    for src, dst in port_node_links:
        link_id = f"__auto_link__{src}__{dst}"
        if link_id not in net.elements:
            net.add_element(LosslessConnectionElement(link_id, from_node=src, to_node=dst))

    # Ejector warm-starts: port pressures, P_jct and the port-chain flows from
    # the ejector's closed-form physics (only where the user gave no guess).
    for _ej, _ej_ports in _ejector_seeds:
        _seed_ejector_warmstart(net, _ej, _ej_ports, schema.edges, nodes_map)

    # Result assembly reports each array under its GUI id (runner.py).
    net._gui_arrays = gui_arrays

    # Merge chambers (#471): inflows on a momentum chamber's "side-target"
    # handle are side streams; the one inflow on its flow handle is the main
    # inlet, which sets the chamber's axis and area. Declared, never inferred.
    for node_schema in schema.nodes:
        if node_schema.type != "momentum_chamber":
            continue
        nid = node_schema.id
        inflows = [
            e
            for e in schema.edges
            if e.target == nid and not (e.data and e.data.get("type") == "thermal")
        ]
        side = [e for e in inflows if e.targetHandle == "side-target"]
        if not side:
            continue
        main = [e for e in inflows if e.targetHandle != "side-target"]
        if len(main) != 1:
            raise ValueError(
                f"Momentum chamber '{nid}' has side streams and {len(main)} main "
                "inlets; it needs exactly one main inlet on its flow handle."
            )
        src = main[0].source
        nodes_map[nid].main_inlet = src if src in net.elements else f"__auto_link__{src}__{nid}"

    # 4. Fourth Pass: Thermal Walls
    for edge in schema.edges:
        if edge.data and edge.data.get("type") == "thermal":
            if edge.source in plenum_node_ids or edge.target in plenum_node_ids:
                continue
            wall_id = edge.id
            data = ThermalWallData(**edge.data)

            # Build layers list
            layers = []
            if data.layers:
                # New multi-layer format
                layers = [
                    WallLayer(
                        thickness=L.thickness, conductivity=L.conductivity, material=L.material
                    )
                    for L in data.layers
                ]
            else:
                # Legacy scalar format
                t = data.thickness if data.thickness is not None else 0.003
                k = data.conductivity if data.conductivity is not None else 20.0
                layers = [WallLayer(thickness=t, conductivity=k)]

            if edge.source in gui_arrays or edge.target in gui_arrays:
                # One wall to an array is one wall per row, each on that
                # row's target footprint; the edge's own area does not apply.
                array_is_a = edge.source in gui_arrays
                arr = gui_arrays[edge.source if array_is_a else edge.target]
                if edge.source in gui_arrays and edge.target in gui_arrays:
                    raise ValueError(f"Thermal wall '{wall_id}' joins two impingement arrays.")
                arr.add_wall(
                    net,
                    wall_id,
                    hot_element=edge.target if array_is_a else edge.source,
                    layers=layers,
                    R_fouling=data.R_fouling,
                    array_is_side_a=array_is_a,
                )
                continue
            net.add_wall(
                ThermalWall(
                    id=wall_id,
                    element_a=edge.source,
                    element_b=edge.target,
                    layers=layers,
                    contact_area=data.area,
                    R_fouling=data.R_fouling,
                )
            )
            net.thermal_coupling_enabled = True

    # Every momentum chamber -- drawn, or inserted at a junction port, an edge
    # junction or a merge -- runs the network's regime (#588): an
    # incompressible network keeps no sonic limit hidden in its chambers.
    for node in net.nodes.values():
        if isinstance(node, MomentumChamberNode):
            node.regime = schema.solver_settings.global_regime

    return net
