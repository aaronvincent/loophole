"""Shared matplotlib styling: render all text in the LaTeX (Computer Modern) font.

Uses matplotlib's *bundled* Computer Modern fonts via the mathtext "cm" font
set, so it needs no system LaTeX/dvipng install -- safe to run headless on a
cluster.  Call ``use_latex_fonts()`` once, after importing pyplot, in each
plotting script.

This deliberately avoids ``text.usetex=True``: that route gives true LaTeX
typesetting but requires a working TeX + dvipng toolchain on every node, which
is fragile in batch environments.
"""

import matplotlib


def use_latex_fonts():
    """Set rcParams so titles, labels, ticks and legends use Computer Modern."""
    matplotlib.rcParams.update({
        "mathtext.fontset": "cm",        # Computer Modern for $...$ math
        "font.family": "serif",
        "font.serif": ["cmr10", "DejaVu Serif"],  # CM Roman for plain text
        "axes.formatter.use_mathtext": True,       # CM digits on offset/sci ticks
        "axes.unicode_minus": False,     # cmr10 has no Unicode minus glyph
    })
