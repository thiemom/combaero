- **A channel's laminar friction is Hagen-Poiseuille down to zero flow.**
  The network friction kernels (`haaland`, `serghides`, `colebrook`) blend
  64/Re into the turbulent fit around Re 2300.
  - **The old floor:** the laminar branch shared the turbulent fits' Re
    floor of 64, so f was 5% low at Re 200, 29% low at Re 64, and saturated at
    1 below that (64x low at Re 1). A channel's drop then went as m^2 at low
    flow instead of linearly, with zero slope at zero flow.
  - **`colebrook` had no laminar branch at all:** it was 63% low at Re 200
    and 54% high at Re 2000.
  - **Now:** the laminar branch has its own floor of Re 1e-6, and Colebrook
    blends like the other two.
  - **Unchanged:** the turbulent fits, to 1e-12 at Re >= 1e4. With the
    default `haaland`, f moves by more than 1% only below Re ~ 300.
  - **One visible side effect:** f reported at exactly zero flow is now
    64/1e-6 rather than 1. The drop there is still zero.
