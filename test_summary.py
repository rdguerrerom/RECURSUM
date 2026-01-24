#!/usr/bin/env python3
"""Generate a summary of all test results."""
import subprocess
import os
from pathlib import Path

test_dir = Path("tests/generated")
test_files = sorted(test_dir.glob("test_*.py"))

results = []
for test_file in test_files:
    name = test_file.stem.replace("test_", "")
    result = subprocess.run(
        ["pytest", str(test_file), "-q"],
        capture_output=True,
        text=True,
        cwd="/home/ruben/Research/Science/Projects/RECURSUM"
    )

    # Parse the output
    output = result.stdout + result.stderr
    last_line = output.strip().split('\n')[-1]

    # Extract pass/fail/skip counts
    passed = failed = skipped = 0
    if "passed" in last_line:
        parts = last_line.split()
        for i, part in enumerate(parts):
            if "passed" in part and i > 0:
                passed = int(parts[i-1])
            if "failed" in part and i > 0:
                failed = int(parts[i-1])
            if "skipped" in part and i > 0:
                skipped = int(parts[i-1])

    status = "PASS" if failed == 0 and passed > 0 else ("SKIP" if passed == 0 else "FAIL")
    results.append((name, status, passed, failed, skipped))
    print(f"{name:25} {status:6} ({passed} passed, {failed} failed, {skipped} skipped)")

print("\n" + "=" * 80)
print("SUMMARY:")
total_pass = sum(1 for _, s, _, _, _ in results if s == "PASS")
total_fail = sum(1 for _, s, _, _, _ in results if s == "FAIL")
total_skip = sum(1 for _, s, _, _, _ in results if s == "SKIP")
print(f"Recurrences: {len(results)} total")
print(f"  - {total_pass} passing all tests")
print(f"  - {total_fail} with failing tests")
print(f"  - {total_skip} skipped")
print(f"\nTests: {sum(p for _, _, p, _, _ in results)} passed, "
      f"{sum(f for _, _, _, f, _ in results)} failed, "
      f"{sum(s for _, _, _, _, s in results)} skipped")
