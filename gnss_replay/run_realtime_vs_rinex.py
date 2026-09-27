"""Reproducible TXT_FAST / TXT_REALTIME / RINEX PPP comparison.

The real-time pass lasts as long as the source GNSS interval. Products and
PPP settings are identical across passes; only the observation source and
the monotonic-clock wait differ.
"""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from compare_realtime_vs_rinex import read_solutions

HERE = Path(__file__).resolve().parent
DEFAULT_PRODUCTS = {
    "nav": r"E:\RTKLIB_Data\NAV\BRDC00WRD_S_20262650000_01D_MN.rnx",
    "sp3": r"E:\RTKLIB_Data\SP3\WUM0MGXNRT_20262641000_02D_05M_ORB.SP3",
    "clk": r"E:\RTKLIB_Data\CLK\WUM0MGXNRT_20262641000_02D_05M_CLK.CLK",
    "bia": r"E:\RTKLIB_Data\BIA\WUM0MGXRTS_20262650000_01D_05M_OSB.BIA",
    "ionex": r"E:\RTKLIB_Data\IONEX\COD0OPSP0D_20262650000_01D_01H_GIM.INX",
    "vmf0": r"E:\RTKLIB_Data\tro\VMF3_20260922.H06",
    "vmf6": r"E:\RTKLIB_Data\tro\VMF3_20260922.H12",
    "orog": r"E:\RTKLIB_Data\tro\orography_ell_5x5",
    "atx": r"E:\RTKLIB_Data\Tables\igs20.atx",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=HERE / "x64/Release/gnss_replay.exe")
    parser.add_argument("--txt", type=Path, required=True)
    parser.add_argument("--rinex-obs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-epochs", type=int, default=0,
                        help="Optional short smoke test; 0 runs the whole file")
    for name, default in DEFAULT_PRODUCTS.items():
        parser.add_argument("--" + name, type=Path, default=Path(default))
    args = parser.parse_args()
    if args.max_epochs < 0:
        parser.error("--max-epochs must be nonnegative")
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    products = {name: getattr(args, name).resolve() for name in DEFAULT_PRODUCTS}
    inputs = {"exe": args.exe.resolve(), "txt": args.txt.resolve(),
              "rinex_obs": args.rinex_obs.resolve(), **products}
    for name, path in inputs.items():
        if not path.is_file():
            parser.error(f"{name} does not exist: {path}")
    manifest = {name: {"path": str(path), "sha256": sha256(path)}
                for name, path in inputs.items()}
    manifest["max_epochs"] = args.max_epochs
    manifest["signal_mode"] = "ppp-safe"
    manifest["adr_unc_max_m"] = 1.0
    (args.output_dir / "inputs.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    def run(name, source):
        command = [str(args.exe), "--source", source,
                   "--result-file", str(args.output_dir / f"{name}_solution.csv"),
                   "--state-dump", str(args.output_dir / f"{name}.stat"),
                   "--trace", str(args.output_dir / f"{name}.trace"),
                   "--trace-level", "2"]
        if source == "rinex":
            command += ["--rinex-obs", str(args.rinex_obs)]
        else:
            command += ["--txt", str(args.txt), "--signal-mode", "ppp-safe",
                        "--adr-unc-max", "1.0"]
        for key, path in products.items():
            command += ["--" + key, str(path)]
        if args.max_epochs:
            command += ["--max-epochs", str(args.max_epochs)]
        print(f"Running {name}: {source}", flush=True)
        with (args.output_dir / f"{name}.stdout.log").open("wb") as stdout, \
             (args.output_dir / f"{name}.stderr.log").open("wb") as stderr:
            subprocess.run(command, stdout=stdout, stderr=stderr, check=True,
                           cwd=args.output_dir)

    run("txt_fast", "txt-fast")
    run("rinex", "rinex")
    fast = read_solutions(args.output_dir / "txt_fast_solution.csv")
    rinex = read_solutions(args.output_dir / "rinex_solution.csv")
    if fast.keys() != rinex.keys():
        raise SystemExit("TXT and RINEX epoch sets differ; see output CSV files. "
                         "Real-time run skipped to preserve identical start/end epochs.")
    run("txt_realtime", "txt-realtime")
    subprocess.run([sys.executable, str(HERE / "compare_realtime_vs_rinex.py"),
                    str(args.output_dir)], check=True)


if __name__ == "__main__":
    main()
