import numpy as np

import combaero as cb


class VatistasVortex:
    """
    Vatistas n-vortex model.
    Reference: Vatistas, Kozel, Mih (1991), Exp. Fluids 11, 73-76.
    DOI: 10.1007/BF00198434

    Parameters
    ----------
    Gamma : float
        Vortex circulation [m^2/s]. Must be > 0.
    r_c : float
        Core radius [m]. Must be > 0.
    n : float
        Shape parameter (default 2.0). Must be >= 1.
        n=2 gives the best fit to most experimental swirl data (Fig. 2a).
    """

    def __init__(self, Gamma: float, r_c: float, n: float = 2.0) -> None:
        if n < 1.0:
            raise ValueError(f"n must be >= 1 (singularities for n < 1); got {n}")
        if Gamma <= 0.0:
            raise ValueError(f"Gamma must be > 0; got {Gamma}")
        if r_c <= 0.0:
            raise ValueError(f"r_c must be > 0; got {r_c}")
        self.Gamma = float(Gamma)
        self.r_c = float(r_c)
        self.n = float(n)

    def V_theta(self, r: np.ndarray | float) -> np.ndarray | float:
        """Tangential velocity [m/s] at radius r [m]. Vectorised."""
        scalar = np.ndim(r) == 0
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        fn = np.frompyfunc(lambda ri: cb.vatistas_v_theta(ri, self.Gamma, self.r_c, self.n), 1, 1)
        result = fn(r_arr).astype(float)
        return float(result[0]) if scalar else result

    def V_theta_max(self) -> float:
        """Peak tangential velocity [m/s] at r = r_c (for all n >= 1)."""
        return float(cb.vatistas_v_theta(self.r_c, self.Gamma, self.r_c, self.n))

    def dV_theta_dr(self, r: np.ndarray | float) -> np.ndarray | float:
        """Analytical d(V_theta)/dr [1/s] at radius r [m]. Vectorised."""
        scalar = np.ndim(r) == 0
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        fn = np.frompyfunc(
            lambda ri: cb.vatistas_v_theta_and_jacobians(ri, self.Gamma, self.r_c, self.n)[1],
            1,
            1,
        )
        result = fn(r_arr).astype(float)
        return float(result[0]) if scalar else result

    def delta_P(self, r: float, rho: float) -> float:
        """
        Static pressure rise P(r) - P(0) [Pa].

        Parameters
        ----------
        r : float
            Radius [m].
        rho : float
            Fluid density [kg/m^3].
        """
        return float(cb.vatistas_delta_p(float(r), float(rho), self.Gamma, self.r_c, self.n))

    def delta_P_bar(self, r: float) -> float:
        """
        Normalised pressure rise (P(r) - P(0)) / (P_inf - P(0)), range [0, 1).

        Exact for n=1 and n=2; numerical for other n.
        The asymptotic total (P_inf - P(0)) is approximated at r_bar=100.
        """
        I_r = cb.vatistas_pressure_integral(float(r) / self.r_c, self.n)
        I_inf = cb.vatistas_pressure_integral(100.0, self.n)
        if I_inf == 0.0:
            return 0.0
        return I_r / I_inf

    def __repr__(self) -> str:
        return f"VatistasVortex(Gamma={self.Gamma}, r_c={self.r_c}, n={self.n})"
