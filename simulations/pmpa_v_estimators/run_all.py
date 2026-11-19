#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run all PMPA protocols in one pass (§7.6 overdetermination, P6 birth,
isotope grounding) with a shared seed.

    python3 run_all.py                 # text only
    python3 run_all.py --plot          # also regenerate all figures
    python3 run_all.py --seed 99       # reproducibility check

The three modules are independent numpy scripts; this wrapper only calls
them with consistent arguments.
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent


def run(script, args):
    cmd = [sys.executable, str(HERE / script)] + args
    print(f"\n{'='*70}\n$ {' '.join(cmd)}\n{'='*70}")
    return subprocess.run(cmd, check=False).returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--N", type=int, default=200000, help="events per point")
    ap.add_argument("--jackknife", action="store_true")
    ap.add_argument("--blocks", type=int, default=24)
    ap.add_argument("--plot", action="store_true", help="regenerate figures")
    a = ap.parse_args()

    code = 0
    code |= run("sim_pmpa.py", ["--all", "--N", str(a.N), "--seed", str(a.seed)]
                + (["--jackknife", "--blocks", str(a.blocks)] if a.jackknife else [])
                + (["--plot", "pmpa_visibility.png"] if a.plot else []))
    code |= run("sim_p6.py", ["--all", "--seed", str(a.seed)]
                + (["--jackknife"] if a.jackknife else [])
                + (["--plot", "p6_birth.png"] if a.plot else []))
    code |= run("sim_isotopes.py", ["--all"]
                + (["--plot", "p6_isotopes.png"] if a.plot else []))
    code |= run("ed_two_slit.py", []
                + (["--plot", "ed_two_slit.png"] if a.plot else []))
    code |= run("charge_channel.py", [])
    code |= run("decoherence_error.py", [])
    code |= run("compton_frequency.py", [])
    print(f"\nall protocols finished (exit={code})")
    sys.exit(code)


if __name__ == "__main__":
    main()
