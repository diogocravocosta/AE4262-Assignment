"""
Abel inversion of a PCO B16 image file.

Usage:
    from abel_b16 import abel_invert_b16
    inverted = abel_invert_b16("B0001.b16")

Requirements:
    pip install numpy PyAbel
    readB16.py must be in the same directory.
"""

import numpy as np
import abel
from readB16 import readB16


def abel_invert_b16(filepath):
    """
    Read a B16 file and perform an inverse Abel transform.

    Parameters
    ----------
    filepath : str
        Path to the .b16 file.

    Returns
    -------
    inverted : np.ndarray
        The Abel-inverted image as a 2D float64 array.
    """
    img = readB16(filepath).astype(np.float64)

    result = abel.Transform(
        img,
        direction="inverse",
        method="hansenlaw",
        symmetry_axis=0,
        verbose=False,
    )

    return result.transform