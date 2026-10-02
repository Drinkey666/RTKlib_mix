"""Single-factor height-bias experiments in an isolated executable.
No production sources/options/products are overwritten. Reference coordinates
are used only by analysis, never passed to the C solver.
"""
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gnss_replay/height_bias_audit"
EXE = OUT / "bin/Release/height_bias_replay.exe"
DATA = Path("E:/RTKLIB_Data")
CONFIG = (ROOT / "src/smartphone_ppp_config.c").read_text()
import re
options = "".join(re.findall(r'"([^"\n]+)"', CONFIG.split("#define PROFILE_OPTIONS")[1].split("typedef")[0]))

def dataset(day):
    key = "076" if day == "0317" else "080"
    fin = DATA if day == "0317" else DATA / "PPP_Result/diag_20260321"
    return {
        "rinex-obs": str(DATA / "OBS" / ("K80-GNSS00GEO_R_20260760533_20M_01S_MO.rnx" if day=="0317" else "GNSS00GEO_R_20260800210_10M_01S_MO.rnx")),
        "nav": str(DATA / f"NAV/BRDC00IGS_R_2026{key}0000_01D_MN.rnx"),
        "sp3": str((fin / "SP3" if day=="0317" else fin) / f"WUM0MGXFIN_2026{key}0000_01D_05M_ORB.SP3"),
        "clk": str((fin / "CLK" if day=="0317" else fin) / f"WUM0MGXFIN_2026{key}0000_01D_30S_CLK.CLK"),
        "bia": str((fin / "BIA" if day=="0317" else fin) / f"WUM0MGXFIN_2026{key}0000_01D_01D_OSB.BIA"),
        "ionex": str(DATA / ("IONEX/COD0OPSFIN_20260760000_01D_01H_GIM.INX" if day=="0317" else "IONEX/COD0OPSPRD_20260800000_01D_01H_GIM.INX")),
        "vmf0": str(DATA / f"tro/VMF3_202603{17 if day=='0317' else 21}.H00"),
        "vmf6": str(DATA / f"tro/VMF3_202603{17 if day=='0317' else 21}.H06"),
        "orog": str(DATA / "tro/orography_ell_5x5"),
        "atx": str(DATA / "Tables/igs20.atx"),
    }

def run(day, case):
    products = dataset(day)
    extra = {}
    effective = options
    if case == "NF3": extra["HEIGHT_AUDIT_NF"]="3"
    elif case == "NO_SMOOTH": effective=options.replace("-DOPPSM=0.90", "-DOPPSM=0")
    elif case == "NO_ZWD_PRIOR": effective=options.replace("-VMF3ZWDSIG=0.30","-VMF3ZWDSIG=0")
    elif case == "IONEX_300": effective=options.replace("-IONCONSINT=1","-IONCONSINT=300")
    elif case == "GPS_BDS": extra["HEIGHT_AUDIT_NAVSYS"]="33"
    elif case == "GPS_GAL": extra["HEIGHT_AUDIT_NAVSYS"]="9"
    elif case == "GPS": extra["HEIGHT_AUDIT_NAVSYS"]="1"
    elif case == "ATX_2408": products["atx"]=str(DATA / "Tables/igs20_2408.atx")
    elif case == "WHU_GIM": products["ionex"]=str(DATA / "IONEX/whrg0800.26i")
    output = OUT / day / case
    output.mkdir(parents=True, exist_ok=True)
    env = {k:v for k,v in os.environ.items() if not k.startswith("HEIGHT_AUDIT_")}
    env.update(extra, HEIGHT_AUDIT_OPTS=effective, HEIGHT_AUDIT_DIAG=str(output / "height_diag.csv"),
               HEIGHT_AUDIT_IONDIAG=str(output / "ion_diag.csv"))
    cmd = [str(EXE), "--source", "rinex"]
    for name,path in products.items():
        if not Path(path).is_file(): raise FileNotFoundError(path)
        cmd += ["--"+name,path]
    cmd += ["--result-file",str(output/"solution.csv"),"--state-dump",str(output/"state.stat"),
            "--use-dump",str(output/"use.csv"),"--trace",str(output/"ppp.trace"),"--trace-level","3"]
    manifest = {"products": products, "overrides":extra,"pppopt":effective,"command":cmd}
    (output/"inputs.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    with (output/"run.log").open("wb") as stream:
        subprocess.run(cmd,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,cwd=ROOT)
    print(day,case,"complete",flush=True)

if __name__=="__main__":
    cases = ["BASE","NF3","NO_SMOOTH","NO_ZWD_PRIOR","IONEX_300","GPS_BDS","GPS_GAL","GPS","ATX_2408"]
    jobs = [(day,case) for day in ("0317","0321") for case in cases] + [("0321","WHU_GIM")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda args:run(*args), jobs))
