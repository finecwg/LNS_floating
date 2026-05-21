"""
Axisymmetric finite-difference operators for cylindrical (r, z) coordinates.

Boundary conventions used by this module's matrices:
  - Axis (r=0, j=0): mirror BC for fields even in r (phi, eta, w^z).
    The (1/r) d_r term collapses by L'Hopital: Delta_H f|_0 = 2 d_rr f|_0,
    which after mirror gives 4(f_1 - f_0)/dr^2.
  - Outer (r=L, j=nr-1): Dirichlet zero, "effectively infinite" — choose L
    large enough that waves don't reach r=L during the simulation.

These operators are exposed for use in the RHS of the monolithic system.
The matrix-side stencils are inlined in the solver's matrix builder.
"""

import numpy as np
from scipy.sparse import lil_matrix


def make_r_grid(nr: int, L: float):
    """Vertex-centered radial grid: r_0 = 0, r_{nr-1} = L."""
    dr = L / (nr - 1)
    r_vec = np.arange(nr) * dr
    return dr, r_vec


def build_Delta_H_axi(nr: int, dr: float, r_vec: np.ndarray):
    """
    1D axisymmetric horizontal Laplacian:  Delta_H f = d_rr f + (1/r) d_r f.

    Returns an (nr, nr) sparse CSR matrix. Boundary handling:
      - j=0:  4(f_1 - f_0)/dr^2  (L'Hopital + mirror)
      - 0 < j < nr-1: standard 3-point stencil with (1/r) d_r centered diff
      - j=nr-1: zero row (Dirichlet zero is enforced by the user; this row
        of the operator returns 0 when applied)
    """
    D = lil_matrix((nr, nr))
    cr = 1.0 / dr**2

    # Axis (j = 0)
    D[0, 0] = -4.0 * cr
    D[0, 1] = +4.0 * cr

    # Interior (0 < j < nr - 1)
    for j in range(1, nr - 1):
        beta = 1.0 / (2.0 * r_vec[j] * dr)
        D[j, j - 1] = cr - beta
        D[j, j] = -2.0 * cr
        D[j, j + 1] = cr + beta

    # Outer Dirichlet: zero row (no contribution)
    return D.tocsr()


def build_Dr(nr: int, dr: float):
    """
    1D centered first derivative d_r.  Mirror BC at axis (so d_r f|_0 = 0
    for even fields).  Zero row at outer boundary.
    """
    D = lil_matrix((nr, nr))
    c = 1.0 / (2.0 * dr)
    for j in range(1, nr - 1):
        D[j, j - 1] = -c
        D[j, j + 1] = +c
    return D.tocsr()
