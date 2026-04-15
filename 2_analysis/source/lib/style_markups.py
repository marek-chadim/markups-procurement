"""style_markups.py — shared matplotlib style for markups-procurement figures.

Python analog of theme_markups.R. Applies Healy-inspired defaults (colorblind
palette, 11pt base, clean axes, 300 DPI, white background) to every figure in
the paper. Reference: Healy (2026), Data Visualization: A Practical
Introduction, Princeton University Press, 2nd ed., Chapter 8.

Usage:
    from style_markups import apply_markups_style, MARKUPS_BLUE, MARKUPS_PINK
    apply_markups_style()
    fig, ax = plt.subplots()
    ax.plot(x, y, color=MARKUPS_BLUE)
    plt.savefig("out.pdf")   # 300 DPI, white bg, tight bbox (from rcParams)

Colors mirror tol_bright in theme_markups.R so R and Python figures share the
exact Paul Tol bright palette.
"""

import matplotlib as _mpl
import matplotlib.pyplot as _plt

# ---- Paul Tol bright palette (colorblind-safe, up to 7 categories) ----
# Source: https://personal.sron.nl/~pault/
MARKUPS_COLORS = {
    "blue":   "#4477AA",
    "red":    "#EE6677",
    "green":  "#228833",
    "yellow": "#CCBB44",
    "cyan":   "#66CCEE",
    "purple": "#AA3377",
    "grey":   "#BBBBBB",
}

# Named shortcuts matching R helper variables
MARKUPS_BLUE = MARKUPS_COLORS["blue"]      # "#4477AA"
MARKUPS_PINK = MARKUPS_COLORS["purple"]    # "#AA3377"
MARKUPS_RED  = MARKUPS_COLORS["red"]       # "#EE6677"
MARKUPS_GREEN = MARKUPS_COLORS["green"]    # "#228833"
MARKUPS_GREY = MARKUPS_COLORS["grey"]      # "#BBBBBB"

# Default color cycle for matplotlib (replaces C0, C1, C2, ...)
MARKUPS_COLOR_CYCLE = [
    MARKUPS_COLORS["blue"],
    MARKUPS_COLORS["purple"],
    MARKUPS_COLORS["green"],
    MARKUPS_COLORS["yellow"],
    MARKUPS_COLORS["cyan"],
    MARKUPS_COLORS["red"],
    MARKUPS_COLORS["grey"],
]

# ---- Single-column journal figure defaults ----
MARKUPS_FIG_WIDTH = 6.5    # inches
MARKUPS_FIG_HEIGHT = 4.0   # inches (golden ratio)
MARKUPS_FIG_DPI = 300


def apply_markups_style():
    """Apply markups-procurement style to matplotlib rcParams.

    Call once at the top of any figure-producing script, after matplotlib is
    imported. Sets base font size, color cycle, white background, tight bbox
    on save, 300 DPI, clean axes (no top/right spines), subtle grid.

    Side effect: mutates matplotlib.rcParams globally for the process.
    """
    _plt.rcParams.update({
        # ---- Font ----
        # STIX Two Text is a Times-clone serif that closely matches LaTeX
        # Computer Modern Roman (used in the paper body). It handles unicode
        # em-dashes, Greek letters, and math symbols without needing
        # text.usetex=True (which would require LaTeX-escaping every label).
        # Falls back through Times and DejaVu Serif if STIX is unavailable.
        "font.family":        "serif",
        "font.serif":         ["STIX Two Text", "Times", "Times New Roman",
                               "DejaVu Serif"],
        # Math text rendered in Computer Modern to match the paper body math
        "mathtext.fontset":   "cm",
        "font.size":          11,
        "axes.titlesize":     12,
        "axes.titleweight":   "bold",
        "axes.labelsize":     11,
        "xtick.labelsize":    10,
        "ytick.labelsize":    10,
        "legend.fontsize":    10,
        "figure.titlesize":   12,
        # ---- Colors (Paul Tol bright) ----
        "axes.prop_cycle":    _mpl.cycler(color=MARKUPS_COLOR_CYCLE),
        # ---- Figure layout ----
        "figure.figsize":     (MARKUPS_FIG_WIDTH, MARKUPS_FIG_HEIGHT),
        "figure.dpi":         100,                     # screen display
        "savefig.dpi":        MARKUPS_FIG_DPI,         # 300 for print
        "savefig.bbox":       "tight",
        "savefig.facecolor":  "white",
        "figure.facecolor":   "white",
        # ---- Axes (clean minimal look) ----
        "axes.grid":          True,
        "axes.grid.axis":     "y",
        "grid.alpha":         0.3,
        "grid.linestyle":     "-",
        "grid.linewidth":     0.5,
        "axes.axisbelow":     True,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.edgecolor":     "#333333",
        # ---- Legend ----
        "legend.frameon":     False,
        "legend.loc":         "best",
    })


def savefig_markups(target, path, **kwargs):
    """Save a figure with markups defaults (dpi=300, bbox_inches='tight',
    facecolor='white'). Passes additional kwargs through to savefig.

    Works with either:
      - a matplotlib Figure object: savefig_markups(fig, path)
      - the pyplot module: savefig_markups(plt, path)
    """
    defaults = {
        "dpi":          MARKUPS_FIG_DPI,
        "bbox_inches":  "tight",
        "facecolor":    "white",
    }
    defaults.update(kwargs)
    target.savefig(path, **defaults)
