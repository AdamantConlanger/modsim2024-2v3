import numpy as np

# https://stackoverflow.com/a/11146645/18375328
def cartesian_product(*arrays):
    broadcastable = np.ix_(*arrays)
    broadcasted = np.broadcast_arrays(*broadcastable)
    rows, cols = np.prod(broadcasted[0].shape), len(broadcasted)
    dtype = np.result_type(*arrays)

    out = np.empty(rows * cols, dtype=dtype)
    start, end = 0, rows
    for a in broadcasted:
        out[start:end] = a.reshape(-1)
        start, end = end, end + rows
    return out.reshape(cols, rows).T

def cartesian_product_for_variations(variations_initials, variations_coefficients):
    cart_prod = cartesian_product(*(variations_initials + variations_coefficients))
    len_initials = len(variations_initials)
    len_both = len_initials + len(variations_coefficients)
    return [(item[0:len_initials], item[len_initials:len_both]) for item in cart_prod]

