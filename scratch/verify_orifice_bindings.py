import combaero as cb
import math

def test_correlations():
    geom = cb.OrificeGeometry()
    geom.d = 0.05
    geom.D = 0.1

    state = cb.OrificeState()
    state.Re_D = 1e5
    state.dP = 1000
    state.rho = 1.2

    correlations = [
        cb.CdCorrelation.ReaderHarrisGallagher,
        cb.CdCorrelation.Stolz,
        cb.CdCorrelation.Miller
    ]

    for corr in correlations:
        try:
            val = cb.Cd_orifice(geom, state, correlation=corr)
            print(f"{corr.name}: {val:.4f}")
        except Exception as e:
            print(f"{corr.name} FAILED: {e}")

if __name__ == "__main__":
    test_correlations()
