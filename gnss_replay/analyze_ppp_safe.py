"""Check RINEX/TXT PPP trajectories and RINEX-independent safe replay."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from analyze_equiv import compare, read_csv, read_states


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solver_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return [{k: v for k, v in row.items() if k != "processing_ms"}
                for row in csv.DictReader(stream)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    root = args.directory
    reference = read_csv(root / "rinex.csv")
    reference_states = read_states(root / "rinex.stat")
    comparisons = {}
    for name in ("safe_with_rinex", "safe_without_rinex", "full"):
        comparisons[name] = compare(
            reference, read_csv(root / f"{name}.csv"), reference_states,
            read_states(root / f"{name}.stat"))
    a, b = "safe_with_rinex", "safe_without_rinex"
    independent = {
        "solver_csv_excluding_processing_ms_equal":
            solver_rows(root / f"{a}.csv") == solver_rows(root / f"{b}.csv"),
        "state_byte_equal": digest(root / f"{a}.stat") == digest(root / f"{b}.stat"),
        "obs_byte_equal": digest(root / f"{a}_obs.csv") == digest(root / f"{b}_obs.csv"),
        "trace_byte_equal": digest(root / f"{a}.trace") == digest(root / f"{b}.trace"),
    }
    print(json.dumps({"comparisons": comparisons,
                      "rinex_independence": independent},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
