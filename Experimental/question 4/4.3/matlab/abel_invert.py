"""
Abel inversion function callable from MATLAB via py.abel_invert.abel_invert().

Returns a numpy array that MATLAB can convert with double().
"""

import numpy as np
import abel


def abel_invert(img_array):
    """
    Inverse Abel transform with automatic flame-axis centering.

    Parameters
    ----------
    img_array : array-like
        2D image (will be cast to float64).

    Returns
    -------
    np.ndarray
        Abel-inverted 2D array.
    """
    img = np.array(img_array, dtype=np.float64)

    # Center on flame axis (center-of-mass of column intensity sums)
    col_sums = img.sum(axis=0)
    center_col = int(np.round(np.average(np.arange(img.shape[1]), weights=col_sums)))
    half_w = min(center_col, img.shape[1] - center_col - 1)
    img_centered = img[:, center_col - half_w : center_col + half_w + 1]

    result = abel.Transform(
        img_centered,
        direction="inverse",
        method="hansenlaw",
        symmetry_axis=0,
        verbose=False,
    )

    return result.transform
