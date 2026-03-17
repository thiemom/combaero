"""Combustion result dataclass for network solver."""

from dataclasses import dataclass


@dataclass
class CombustionResult:
    """Complete thermodynamic state of burned/mixed stream."""

    # Composition
    X: list[float]  # mole fractions [14 species], sums to 1
    Y: list[float]  # mass fractions [14 species], sums to 1
    mw: float  # mixture molecular weight [g/mol]

    # Thermodynamic state - STATIC conditions (low-velocity plenum assumption)
    T: float  # static temperature [K]
    P: float  # static pressure [Pa]
    m_dot: float  # total mass flow [kg/s]

    # Derived properties (computed via CombAero)
    h: float  # specific enthalpy [J/mol]
    cp: float  # heat capacity [J/(mol*K)]
    rho: float  # density [kg/m^3]
    gamma: float  # isentropic expansion coefficient [-]
    a: float  # speed of sound [m/s]

    # Combustion diagnostics
    phi: float  # equivalence ratio [-]  (0 if no fuel)
    T_adiabatic: float  # adiabatic flame temperature [K]
    eta: float  # combustion efficiency applied [-]
    Q_released: float  # heat released [J/s] = [W]
