fix: resolve NASA-9 validation test file-not-found error

The nasa9_gas fixture in test_nasa9_polynomials.py passed a bare filename
"combaero_nasa9.yaml" to Cantera, which only searches cwd (repo root)
and its own data directory. The file lives in cantera_validation_tests/.

Fix: use Path(__file__).parent / "combaero_nasa9.yaml" for an absolute
path that works regardless of working directory.
