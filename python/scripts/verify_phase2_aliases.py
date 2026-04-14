import os
import sys

# Ensure we can import the built module at the START of the path
sys.path.insert(0, os.path.abspath("build"))

# Import the top-level _core from build/
import _core as core

print(f"Loading core from: {core.__file__}")

try:
    print("Testing Channel terminology...")
    # Test ChannelResult (Heat Transfer)
    print(f"Instantiating core.ChannelResult ({type(core.ChannelResult)})...")
    res = core.ChannelResult()
    print("ChannelResult created successfully.")

    # Test ChannelSolverResult (Solver)
    print(f"Instantiating core.ChannelSolverResult ({type(core.ChannelSolverResult)})...")
    res_s = core.ChannelSolverResult()
    print("ChannelSolverResult created successfully.")

    print("\nTesting backward compatibility (Pipe aliases)...")
    # Test PipeResult alias
    print(f"Instantiating core.PipeResult ({type(core.PipeResult)})...")
    res_p = core.PipeResult()
    print("PipeResult alias works.")

    print("\nRefactor Phase 2 Verification: SUCCESS")
except Exception as e:
    print("\nRefactor Phase 2 Verification: FAILED")
    print(f"Error: {e}")
    # Print more
    import traceback

    traceback.print_exc()
    sys.exit(1)
