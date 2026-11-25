"""plot2d.py -- plot plasma2d output: field maps and the history of the slab averages.

usage:  python plot2d.py fields2d.csv [history2d.csv] [--save fig.png]

fields CSV columns: x_m, y_m, <species>_m3 ..., phi_V   (written by plasma2d)
history CSV columns: time_s, n_e_mean_m3, <species>_mean_m3 ..., charge_C_per_m

STAND-IN numbers only: nothing in the maps is the lab's plasma.
"""
import argparse
import csv

import matplotlib.pyplot as plt


def read_csv(path):
    # csv module, not pandas: nothing extra to install (same choice as plot.py).
    with open(path, newline="") as f:
        rows = [r for r in csv.reader(f) if r]
    head, body = [h.strip() for h in rows[0]], rows[1:]
    cols = {h: [] for h in head}
    for r in body:
        for h, v in zip(head, r):
            cols[h].append(float(v))
    return cols


def grid_of(cols):
    # Cell centers come out row by row (x fastest), so the grid shape is recoverable
    # from how many times y repeats.
    xs = sorted(set(cols["x_m"]))
    ys = sorted(set(cols["y_m"]))
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fields")
    ap.add_argument("history", nargs="?")
    ap.add_argument("--save")
    a = ap.parse_args()

    f = read_csv(a.fields)
    xs, ys = grid_of(f)
    nx, ny = len(xs), len(ys)
    names = [k for k in f if k not in ("x_m", "y_m")]
    panels = len(names) + (1 if a.history else 0)
    fig, axes = plt.subplots(1, panels, figsize=(3.6 * panels, 3.6), squeeze=False)
    axes = axes[0]
    # mm on the axes: the slab is centimeters across.
    ext = [xs[0] * 1e3 - 0.5 * (xs[1] - xs[0]) * 1e3, xs[-1] * 1e3 + 0.5 * (xs[1] - xs[0]) * 1e3,
           ys[0] * 1e3 - 0.5 * (ys[1] - ys[0]) * 1e3, ys[-1] * 1e3 + 0.5 * (ys[1] - ys[0]) * 1e3]
    for ax, k in zip(axes, names):
        z = [f[k][j * nx:(j + 1) * nx] for j in range(ny)]
        im = ax.imshow(z, origin="lower", extent=ext, aspect="equal")
        ax.set_title(k)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        fig.colorbar(im, ax=ax, shrink=0.8)
    if a.history:
        h = read_csv(a.history)
        ax = axes[-1]
        for k in h:
            if k.endswith("_mean_m3"):
                ax.semilogy(h["time_s"], h[k], label=k.replace("_mean_m3", ""))
        ax.set_xlabel("time (s)")
        ax.set_ylabel("slab average (m^-3)")
        ax.legend()
        ax.set_title("history")
    fig.suptitle("plasma-2d stand-in case (placeholders)")
    fig.tight_layout()
    if a.save:
        fig.savefig(a.save, dpi=130)
    else:
        plt.show()


if __name__ == "__main__":
    main()
