import numpy as np

def get_greville_points(breaks, periodic = False, spline_degree = 3):
    if periodic:
        length = breaks[-1] - breaks[0]
        knots = list(breaks[-spline_degree-1:-1] - length) + list(breaks) + \
                list(breaks[1:spline_degree+1] + length)
    else:
        knots = [breaks[0]]*spline_degree + list(breaks) + [breaks[-1]]*spline_degree

    # Calculate the moving average
    greville_pts = np.convolve(knots[1:-1], np.ones(spline_degree), 'valid') / spline_degree

    if periodic:
        greville_pts = greville_pts[:len(breaks)-1]
        greville_pts = ((greville_pts - breaks[0]) % (length)) + breaks[0]
        greville_pts = np.sort(greville_pts)

    return greville_pts
