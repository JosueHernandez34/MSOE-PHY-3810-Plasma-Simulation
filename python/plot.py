"""plot.py -- plot plasma0d output, optionally over measured data.

usage:  python plot.py output.csv [measured.csv] [--save fig.png]

Model CSV columns: time_s, <species>_m3 ..., n_e_m3, Te_eV  (written by plasma0d).
Measured CSV: see data/measured_template.csv. Any column whose name also appears
in the model CSV is overlaid as points; extra columns are ignored.
"""
import argparse
import csv
import sys

import matplotlib.pyplot as plt


def read_csv(path):
    # csv module, not pandas: nothing extra to install.
    with open(path, newline="") as f:
        rows = [r for r in csv.reader(f) if r and not r[0].lstrip().startswith("#")]
    head, body = rows[0], rows[1:]
    cols = {h.strip(): [] for h in head}
    for r in body:
        for h, v in zip(head, r):
            v = v.strip()
            cols[h.strip()].append(float(v) if v else float("nan"))
    return cols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("measured", nargs="?")
    ap.add_argument("--save")
    a = ap.parse_args()

    m = read_csv(a.model)
    d = read_csv(a.measured) if a.measured else None
    t = m["time_s"]
    dens = [k for k in m if k.endswith("_m3")]

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 7))
    for k in dens:
        ax1.loglog(t[1:], m[k][1:], label=k.replace("_m3", ""))   # skip t=0 on a log axis
        if d and k in d and "time_s" in d:
            ax1.loglog(d["time_s"], d[k], "o", mfc="none", label=k.replace("_m3", "") + " (measured)")
    ax1.set_ylabel("density (m$^{-3}$)")
    ax1.legend(ncol=2, fontsize=8)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.set_title("plasma-0d (stand-in parameters unless you changed the config)")

    ax2.semilogx(t[1:], m["Te_eV"][1:], label="Te model")
    if d and "Te_eV" in d and "time_s" in d:
        ax2.semilogx(d["time_s"], d["Te_eV"], "o", mfc="none", label="Te measured")
    ax2.set_ylabel("Te (eV)")
    ax2.set_xlabel("time (s)")
    ax2.legend(fontsize=8)
    ax2.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    if a.save:
        fig.savefig(a.save, dpi=150)
        print("saved", a.save)
    else:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
