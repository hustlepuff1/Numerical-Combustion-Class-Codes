# grid.py

import numpy as np

def make_1d_grid(ni, x_left=0.0, x_right=1.0):
    """
    Uniform 1D grid for shock tube test.
    Returns cell-center coordinates x and spacing dx.
    """
    x = np.linspace(x_left, x_right, ni)
    dx = x[1] - x[0]
    return x, dx

def make_2d_rect_grid(ni, nj, x_left=0.0, x_right=1.0, y_bottom=0.0, y_top=1.0):
    """
    Uniform 2D Cartesian grid for testing 2D Euler.
    Returns:
        x, y : (ni, nj) cell-center coordinates
        dx, dy : spacings
    """
    x1d = np.linspace(x_left, x_right, ni)
    y1d = np.linspace(y_bottom, y_top, nj)
    dx = x1d[1] - x1d[0]
    dy = y1d[1] - y1d[0]
    x, y = np.meshgrid(x1d, y1d, indexing="ij")
    return x, y, dx, dy

def make_wedge_grid(ni, nj,
                    x_left=0.0, x_right=1.0,
                    y_bottom=0.0, y_top=0.2,
                    x_ramp=0.3, theta_deg=15.0):
    """
    Cartesian grid for ODWE.
    Bottom boundary behaves like a wedge:
      - for x < x_ramp: horizontal wall (0 deg)
      - for x >= x_ramp: ramp at theta_deg

    Returns:
        x, y       : (ni, nj) cell centers
        dx, dy     : spacings
        wall_theta : (ni,) wedge angle (radians) at bottom cells
    """
    x1d = np.linspace(x_left, x_right, ni)
    y1d = np.linspace(y_bottom, y_top, nj)
    dx = x1d[1] - x1d[0]
    dy = y1d[1] - y1d[0]
    x, y = np.meshgrid(x1d, y1d, indexing="ij")

    theta = np.zeros(ni)
    theta[x1d >= x_ramp] = np.deg2rad(theta_deg)

    return x, y, dx, dy, theta