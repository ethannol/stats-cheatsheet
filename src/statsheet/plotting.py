"""Shared matplotlib setup. All analysis scripts call setup() once and save via savefig()."""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt

from . import FIG_DIR

matplotlib.use("Agg")

# brand-neutral categorical palette (colour-blind safe, distinct in print)
COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3",
          "#8C8C8C"]


def setup():
    plt.rcParams.update({
        "figure.figsize": (9, 5),
        "figure.dpi": 110,
        "savefig.dpi": 130,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": matplotlib.cycler(color=COLORS),
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "legend.frameon": False,
    })


def savefig(fig, name: str, section: str):
    """Save to figures/<section>/<name>.png and close."""
    out = FIG_DIR / section
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"    [fig] {path.relative_to(FIG_DIR.parent)}")
    return path


def banner(title: str):
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78)


def note(text: str):
    """Print an interpretation note, indented, wrapped at ~90 chars."""
    import textwrap
    for para in text.strip().split("\n\n"):
        print(textwrap.fill(" ".join(para.split()), width=90, initial_indent="  > ",
                            subsequent_indent="    "))
