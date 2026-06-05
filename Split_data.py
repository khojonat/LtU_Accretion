import numpy as np

def split_paired_array(arr, first_is_x: bool = True):
    """
    Generic helper to unpack an interleaved 1D array of the form
    [x0, y0, x1, y1, ...] or [y0, x0, y1, x1, ...].
    This is used to unpack data points from other plots.

    Parameters
    ----------
    arr : array_like
        Flat array with an even number of elements, storing paired values.
    first_is_x : bool, optional
        If True, interpret as [x0, y0, x1, y1, ...] and return (x, y).
        If False, interpret as [y0, x0, y1, x1, ...] and return (x, y).

    Returns
    -------
    x, y : np.ndarray
        Arrays of the same length containing the unpacked x and y values.
    """
    arr = np.asarray(arr)
    if arr.size % 2 != 0:
        raise ValueError("split_paired_array expects an array with an even number of elements.")

    first = arr[0::2]
    second = arr[1::2]

    if first_is_x:
        return first, second
    return second, first