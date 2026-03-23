import os
import subprocess
import sys
import time
from pathlib import Path
from typing import TypedDict

class ExampleResult(TypedDict):
    name: str
    success: bool
    duration: float
    error: str | None

def run_example(example_path: Path, timeout: int = 120) -> ExampleResult:
    """Run a single example script and return the result."""
    start_time = time.perf_counter()
    try:
        # Run with the same python interpreter and headless plotting
        env = os.environ.copy()
        env["COMBAERO_SAVE_PLOTS"] = "1"

        result = subprocess.run(
            [sys.executable, str(example_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
        duration = time.perf_counter() - start_time
        if result.returncode == 0:
            return {
                "name": example_path.name,
                "success": True,
                "duration": duration,
                "error": None,
            }
        else:
            return {
                "name": example_path.name,
                "success": False,
                "duration": duration,
                "error": result.stderr.strip() or "Process exited with non-zero return code.",
            }
    except subprocess.TimeoutExpired:
        duration = time.perf_counter() - start_time
        return {
            "name": example_path.name,
            "success": False,
            "duration": duration,
            "error": f"Timeout expired after {timeout}s",
        }
    except Exception as e:
        duration = time.perf_counter() - start_time
        return {
            "name": example_path.name,
            "success": False,
            "duration": duration,
            "error": str(e),
        }

def main() -> None:
    """Discover and run all Python examples."""
    examples_dir = Path("python/examples")
    if not examples_dir.exists():
        print(f"Error: {examples_dir} not found.")
        sys.exit(1)

    # Discover all .py files, excluding utilities
    example_files = sorted([
        f for f in examples_dir.glob("*.py")
        if f.name != "plot_utils.py"
    ])

    print(f"Found {len(example_files)} examples in {examples_dir}")
    print("-" * 60)

    results: list[ExampleResult] = []

    for i, example_file in enumerate(example_files, 1):
        print(f"[{i}/{len(example_files)}] Running {example_file.name}...", end=" ", flush=True)
        res = run_example(example_file)
        results.append(res)

        if res["success"]:
            print(f"SUCCESS ({res['duration']:.2f}s)")
        else:
            print("FAILED")

    print("\n" + "=" * 60)
    print("SUMMARY REPORT")
    print("=" * 60)

    successful = [r for r in results if r["success"]]
    failed = [r for r in results if not r["success"]]

    total_time = sum(r["duration"] for r in results)

    print(f"Total Examples: {len(results)}")
    print(f"Successful:     {len(successful)}")
    print(f"Failed:         {len(failed)}")
    print(f"Total Wall Time: {total_time:.2f}s")
    print("-" * 60)

    if failed:
        print("\nFAILED EXAMPLES DIAGNOSIS:")
        for res in failed:
            print(f"\n--- {res['name']} ---")
            print(f"Error:\n{res['error']}")
        print("-" * 60)
        sys.exit(1)
    else:
        print("\nAll examples ran successfully!")
        sys.exit(0)

if __name__ == "__main__":
    main()
