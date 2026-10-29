"""ISO 5167-2 orifice metering run: tap differential versus permanent loss.

Source
    validation/sources.yaml, ``wanted`` entry "ISO 5167-2:2003" (issues
    #645, #649-#654, #793). The standard itself is not held under sources/;
    the equations below are the ones #649/#650 implemented, typed again here
    independently of include/orifice.h: Eq. (4) Reader-Harris/Gallagher C
    with its tap terms, the flow equation q_m = C E eps (pi d^2/4)
    sqrt(2 dp rho1), the expansibility Eq. (5), and 5.4.1's pressure loss
    (the long form, and the short form 1 - beta^1.9).

What is built
    A metering run: air at 300 K in a 100 mm pipe, a square-edged plate with
    flange taps, and a pipe downstream, between a 2 bar supply and a 1.95 bar
    sink. Three plates, beta = 0.3, 0.5, 0.7. The network solves the flow;
    the plate drives from the approach pipe's static pressure with the
    velocity-of-approach factor E, and between pipes it charges the PERMANENT
    loss to the stagnation pressure rather than the tap differential (#645).
    The flow is evaluated with D and D/2 taps' C (the taps the 5.4.1 loss
    relation belongs to); the ``tap`` setting only picks the differential
    reported as dP_tap, scaled by (C_DD2 / C_tap)^2 at the same flow.

What is reproduced
    1. C: the network's against Eq. (4) at D-D/2 taps and the solved Re_D.
    2. Metering: read the network's flange-tap differential the way a user
       of the standard would -- flange-tap C iterated on Re_D, eps from
       Eq. (5) -- and get the flow back. The remainder is the network's
       compressible orifice kernel against ISO's empirical expansibility;
       ``--incompressible`` runs the network incompressible with eps = 1,
       where the read-back closes to solver precision.
    3. Permanent loss: the network's stagnation drop across the plate against
       5.4.1's long form, both as the network evaluates it (D-D/2 C and
       differential) and as a flange-tap user would (flange C and
       differential), and the short form 1 - beta^1.9.

What is NOT reproduced
    Annex A's tables of C (the standard is not held), corner taps, the
    straight-length requirements, the standard's uncertainty statements, and
    any measured data.

Validation label
    FIDELITY of the implementation against the standard's own equations
    (docs/VALIDATION_POLICY.md). There is no measurement here: agreement
    means the network does what ISO 5167-2 says, not that ISO 5167-2 fits a
    rig.

Runs from an installed package or a repo checkout: it reads no data files.

The same network at beta = 0.5 opens in the GUI from
Examples -> "ISO 5167-2 orifice metering run"
(gui/examples/iso5167_orifice_metering.json).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

D = 0.100  # [m] pipe bore; above 71.12 mm, so no small-pipe term
L_UP = 2.0  # [m] upstream pipe, 20 D
L_DOWN = 1.0  # [m] downstream pipe, 10 D
P_SUPPLY = 2.00e5  # [Pa] total
P_SINK = 1.95e5  # [Pa] total
T_SUPPLY = 300.0  # [K]
BETAS = (0.3, 0.5, 0.7)
KAPPA = 1.4  # isentropic exponent of air, for Eq. (5)
X_AIR = cb.species.dry_air()

# Eq. (4)'s tap terms: (L1, L2') = (l1/D, l2'/D).
FLANGE = (0.0254 / D, 0.0254 / D)
D_D2 = (1.0, 0.47)


def solve_run(beta: float, regime: str) -> dict[str, float]:
    """Pipe -> flange-tapped plate -> pipe; return the plate's diagnostics."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("supply", Pt=P_SUPPLY, Tt=T_SUPPLY))
    net.add_node(MomentumChamberNode("up", regime=regime))
    net.add_node(MomentumChamberNode("down", regime=regime))
    net.add_node(PressureBoundary("sink", Pt=P_SINK, Tt=T_SUPPLY))
    net.add_element(
        ChannelElement("pipe_up", "supply", "up", length=L_UP, diameter=D, regime=regime)
    )
    net.add_element(
        OrificeElement(
            "plate",
            "up",
            "down",
            diameter=beta * D,
            correlation="ReaderHarrisGallagher",
            tap="Flange",
            regime=regime,
        )
    )
    net.add_element(
        ChannelElement("pipe_down", "down", "sink", length=L_DOWN, diameter=D, regime=regime)
    )
    res = NetworkSolver(net).solve()
    if not res["__success__"]:
        raise RuntimeError(f"solve failed at beta {beta}: {res['__message__']}")
    return {k[len("plate.") :]: v for k, v in res.items() if k.startswith("plate.")}


# --- ISO 5167-2:2003, typed from the standard -----------------------------


def iso_C(beta: float, Re_D: float, taps: tuple[float, float]) -> float:
    """Eq. (4), Reader-Harris/Gallagher (D >= 71.12 mm)."""
    L1, L2 = taps
    A = (19000.0 * beta / Re_D) ** 0.8
    M2 = 2.0 * L2 / (1.0 - beta)
    return (
        0.5961
        + 0.0261 * beta**2
        - 0.216 * beta**8
        + 0.000521 * (1e6 * beta / Re_D) ** 0.7
        + (0.0188 + 0.0063 * A) * beta**3.5 * (1e6 / Re_D) ** 0.3
        + (0.043 + 0.080 * math.exp(-10.0 * L1) - 0.123 * math.exp(-7.0 * L1))
        * (1.0 - 0.11 * A)
        * beta**4
        / (1.0 - beta**4)
        - 0.031 * (M2 - 0.8 * M2**1.1) * beta**1.3
    )


def iso_eps(beta: float, dp: float, p1: float) -> float:
    """Eq. (5), expansibility factor of an orifice plate."""
    return 1.0 - (0.351 + 0.256 * beta**4 + 0.93 * beta**8) * (
        1.0 - (1.0 - dp / p1) ** (1.0 / KAPPA)
    )


def iso_loss_ratio(beta: float, C: float) -> float:
    """5.4.1: dw / dp, long form."""
    s = math.sqrt(1.0 - beta**4 * (1.0 - C * C))
    return (s - C * beta**2) / (s + C * beta**2)


def iso_flow(beta: float, dp: float, p1: float, T1: float, eps_on: bool) -> tuple[float, float]:
    """Mass flow from a flange-tap differential, (q_m, C): C iterated on Re_D."""
    rho1 = cb.density(T1, p1, X_AIR)
    mu1 = cb.transport_state(T1, p1, X_AIR).mu
    E = 1.0 / math.sqrt(1.0 - beta**4)
    eps = iso_eps(beta, dp, p1) if eps_on else 1.0
    area = math.pi * (beta * D) ** 2 / 4.0
    C = 0.6
    for _ in range(100):
        q_m = C * E * eps * area * math.sqrt(2.0 * dp * rho1)
        C_new = iso_C(beta, 4.0 * q_m / (math.pi * D * mu1), FLANGE)
        if abs(C_new - C) < 1e-14:
            break
        C = C_new
    return q_m, C


def _plot(rows: list[dict[str, float]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib not installed: pip install combaero[examples] for the plot)")
        return
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from plot_utils import show_or_save

    b = [i / 100 for i in range(10, 76)]
    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    ax.plot(b, [iso_loss_ratio(x, iso_C(x, 1e5, D_D2)) for x in b], "-", label="5.4.1 long form")
    ax.plot(b, [1.0 - x**1.9 for x in b], "--", label="5.4.1 short form, 1 - beta^1.9")
    ax.plot(
        [r["beta"] for r in rows],
        [r["dP_permanent"] / r["dP_tap"] for r in rows],
        "o",
        label="network: dPt across plate / flange dp",
    )
    ax.set_xlabel("beta = d/D")
    ax.set_ylabel("permanent loss / tap differential")
    ax.set_title("ISO 5167-2 orifice plate -- fidelity (the standard's equations)")
    ax.legend(fontsize=8)
    show_or_save(fig, "iso5167_orifice_metering.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--incompressible", action="store_true", help="eps = 1 on both sides")
    parser.add_argument("--no-plot", action="store_true", help="skip the loss-ratio plot")
    args = parser.parse_args()
    regime = "incompressible" if args.incompressible else "compressible"

    print("ISO 5167-2:2003 orifice metering run -- FIDELITY against the standard's equations")
    print(
        f"Air {T_SUPPLY:.0f} K, D = {D * 1000:.0f} mm, flange taps, "
        f"{P_SUPPLY / 1e5:g} -> {P_SINK / 1e5:g} bar, {regime}\n"
    )
    print("     |  C at D-D/2 taps   |  metering (flange taps)            |  permanent loss dw")
    print(
        "beta |  net      Eq. 4    |  dp [Pa]  m_dot net  m_dot ISO  err  |"
        "  net [Pa]  vs 5.4.1 net taps  vs flange user  dw/dp  short"
    )
    rows = []
    worst = {"C": 0.0, "m": 0.0, "k_net": 0.0, "k_user": 0.0}
    for beta in BETAS:
        r = solve_run(beta, regime)
        p1, T1 = r["P_approach"], r["T_in"]  # static, approach pipe
        mu1 = cb.transport_state(T1, p1, X_AIR).mu
        Re_D = 4.0 * r["m_dot"] / (math.pi * D * mu1)
        C_dd2 = iso_C(beta, Re_D, D_D2)
        dp = r["dP_tap"]  # flange taps
        q_iso, C_fl = iso_flow(beta, dp, p1, T1, eps_on=not args.incompressible)
        err_m = q_iso / r["m_dot"] - 1.0
        dw = r["dP_permanent"]
        # The network's own relation: D-D/2 C and differential at the same flow.
        dw_net_taps = iso_loss_ratio(beta, r["Cd"]) * dp * (C_fl / r["Cd"]) ** 2
        dw_user = iso_loss_ratio(beta, C_fl) * dp
        e_net, e_user = dw / dw_net_taps - 1.0, dw / dw_user - 1.0
        worst["C"] = max(worst["C"], abs(r["Cd"] - C_dd2))
        worst["m"] = max(worst["m"], abs(err_m))
        worst["k_net"] = max(worst["k_net"], abs(e_net))
        worst["k_user"] = max(worst["k_user"], abs(e_user))
        rows.append(r)
        print(
            f"{beta:4.1f} |  {r['Cd']:.4f}   {C_dd2:.4f}   |  {dp:7.0f}  {r['m_dot']:9.5f}  "
            f"{q_iso:9.5f}  {err_m:+.2%} |  {dw:8.0f}  {e_net:+17.2%}  {e_user:+14.2%}  "
            f"{dw / dp:.3f}  {1.0 - beta**1.9:.3f}"
        )
    print(f"\nC vs Eq. (4) at the solved Re_D: worst |difference| {worst['C']:.1e}")
    print(f"Flow read back from the flange differential: worst |error| {worst['m']:.2%}")
    print(f"Permanent loss vs 5.4.1 at the network's taps: worst {worst['k_net']:.1e}")
    print(f"Permanent loss vs 5.4.1 applied with flange C: worst {worst['k_user']:.2%}")
    print(
        "The flange differential overstates what the plate costs the system by "
        f"{rows[0]['dP_tap'] / rows[0]['dP_permanent'] - 1:.0%} at beta {BETAS[0]} and "
        f"{rows[-1]['dP_tap'] / rows[-1]['dP_permanent'] - 1:.0%} at beta {BETAS[-1]}."
    )

    if not args.no_plot:
        _plot(rows)


if __name__ == "__main__":
    main()
