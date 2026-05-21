"""
Axisymmetric floating-object solver (folder 12).

Adapts 03_Floating_Object/floating_object_v3.py to cylindrical (r, z) coordinates
using the axisymmetric operators from 11_Surface_Dynamics_3dAxiSymm.

Body
----
Flat-bottomed cylinder (disk) of radius a, mass m, constrained to vertical
motion by axisymmetry. Wetted region B = {r : r < a}.

Differences from 2D Cartesian v3 (per 12/STRATEGY.md §8a)
----------------------------------------------------------
* Operators are axisymmetric:
    - scalar Laplacian:  d_rr + (1/r) d_r + d_zz
    - vector Laplacian for w^r:  + extra  -w^r / r^2  on the diagonal
    - mirror BC at r=0 (axis) for even fields; w^r=0 there by regularity
    - Dirichlet zero at r=L ("effectively infinite" — choose L large)

* Body integral weight is the cylindrical area element 2*pi*r dr instead of dx:
    A[v_b, p_s(j)] = -dt * 2*pi * r_j * dr     for body nodes j.

* Mass non-dimensionalisation is genuinely 3D:
    M = m / (rho * L_c^3)     (true 3D mass)
  instead of M = m / (rho * L_c^2) (mass per unit length in 2D Cartesian).

* Hydrostatic equilibrium (3D Archimedes for a flat-bottomed cylinder):
    |z_b^eq| = M / (pi * a^2)

State vector (size 3*nr*nz + 2*nr + 2)
--------------------------------------
    x = [ w^z (nz x nr),
          w^r (nz x nr),
          phi (nz x nr),
          eta (nr),
          p_s (nr),
          z_b (1),
          v_b (1) ]

All v3 corrections from 03_Floating_Object/floating_object_v3.py are preserved:
  - Drop eq (2.15e); use bulk heat eq inclusively at surface with one-sided d_zz.
  - Stress-free condition with factor 2 on phi_rz:
        2 phi_rz + d_r w^z + d_z w^r = 0      at z = 0, free surface.
"""

from dataclasses import dataclass, field
import sys
import time as timer
from pathlib import Path

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu

# Import axisymmetric operators from folder 11
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "11_Surface_Dynamics_3dAxiSymm"))
from operators_axi import build_Delta_H_axi, build_Dr, make_r_grid  # noqa: E402


# =============================================================================
# Configuration
# =============================================================================


@dataclass
class PhysicalParams:
    """Dimensional physical parameters (CGS)."""
    L_cm: float = 50.0    # Domain radius (cm); large for "effectively infinite"
    D_cm: float = 20.0    # Domain depth (cm)
    g: float = 980.0      # Gravity (cm/s^2)
    sigma: float = 174.9  # Surface tension (dyne/cm)
    rho: float = 1.0      # Density (g/cm^3)
    nu: float = 0.0005    # Kinematic viscosity (cm^2/s)


@dataclass
class BodyParams:
    """Axisymmetric disk floating body (flat-bottomed cylinder)."""
    radius_cm: float = 2.5            # disk radius a (cm)
    mass_g: float = 7.0               # total 3D mass m (g)
    initial_z_cm: float = 0.0         # initial bottom position z_b (cm)
    initial_vz_cm_s: float = 0.0      # initial vertical velocity v_b (cm/s)


@dataclass
class NumericalParams:
    nr: int = 200
    nz: int = 60
    nt: int = 2000
    total_time_s: float = 1.0
    unit_length: float = 2.5    # non-dim length scale (cm)
    unit_time: float = 0.05     # non-dim time scale (s)


@dataclass
class SimConfig:
    phys: PhysicalParams = field(default_factory=PhysicalParams)
    body: BodyParams = field(default_factory=BodyParams)
    num: NumericalParams = field(default_factory=NumericalParams)
    save_every: int | None = None
    verbose: bool = True


# =============================================================================
# Offsets
# =============================================================================


def _compute_offsets(nr: int, nz: int) -> dict:
    Nb = nr * nz
    return {
        "w3": 0,
        "w1": Nb,
        "phi": 2 * Nb,
        "eta": 3 * Nb,
        "ps": 3 * Nb + nr,
        "zb": 3 * Nb + 2 * nr,
        "vb": 3 * Nb + 2 * nr + 1,
        "total": 3 * Nb + 2 * nr + 2,
        "Nb": Nb,
        "nr": nr,
        "nz": nz,
    }


# =============================================================================
# Monolithic matrix assembly  (axisymmetric + body)
# =============================================================================


def build_monolithic_matrix_axi_body(
    nr, nz, dr, dz, dt, Fr, We, Re, M, body_mask, r_vec,
):
    """Assemble axisymmetric monolithic matrix with floating-body coupling.

    Block layout (rows == columns, same ordering as the state vector):
        1) w^z   (nr*nz rows)  : heat eq + bottom non-slip + v3 inclusive surface
        2) w^r   (nr*nz rows)  : vector heat eq;
                                 surface = stress-free (free) | horiz no-slip (body)
        3) phi   (nr*nz rows)  : Laplace + bottom ghost + dynamic BC (incl. p_s)
        4) eta   (nr rows)     : free=>kinematic;  body=>eta = z_b
        5) p_s   (nr rows)     : free=>p_s=0;      body=>kinematic match
        6) z_b   (1 row)       : z_b - dt*v_b = z_b^n
        7) v_b   (1 row)       : Newton's 2nd law, integral 2*pi*r_j*dr
    """
    off = _compute_offsets(nr, nz)
    N = off["total"]
    A = lil_matrix((N, N))

    cr = 1.0 / dr**2
    cz = 1.0 / dz**2
    alpha = (dz / dr) ** 2
    dt_Re = dt / Re
    two_dt_Re = 2.0 * dt / Re

    def w3(i, j): return off["w3"] + i * nr + j
    def w1(i, j): return off["w1"] + i * nr + j
    def phi(i, j): return off["phi"] + i * nr + j
    def eta(j): return off["eta"] + j
    def ps(j): return off["ps"] + j
    zb = off["zb"]
    vb = off["vb"]

    # ============================================================
    # Block 1: w^z   (scalar heat eq + v3 inclusive surface)
    # ============================================================
    for i in range(nz):
        for j in range(nr):
            row = w3(i, j)

            # Outer wall: Dirichlet zero (non-slip in far-wall limit)
            if j == nr - 1:
                A[row, row] = 1.0
                continue

            # Bottom (i=0): non-slip implicit
            if i == 0:
                A[row, w3(0, j)] = 1.0
                A[row, phi(1, j)] = 1.0 / dz
                A[row, phi(0, j)] = -1.0 / dz
                continue

            # Radial part of axisymm scalar Laplacian
            if j == 0:
                # Axis (L'Hopital + mirror): Delta_H f|_0 = 4(f_1 - f_0)/dr^2
                rad_diag = -4.0 * cr
                rad_jp = +4.0 * cr
                rad_jm = 0.0
            else:
                beta = 1.0 / (2.0 * r_vec[j] * dr)
                rad_diag = -2.0 * cr
                rad_jp = cr + beta
                rad_jm = cr - beta

            # Vertical part
            if i == nz - 1:
                # v3 inclusive: one-sided d_zz at the surface
                z_diag = +1.0 * cz
                z_extra = [(nz - 2, -2.0 * cz), (nz - 3, +1.0 * cz)]
            else:
                z_diag = -2.0 * cz
                z_extra = [(i + 1, +1.0 * cz), (i - 1, +1.0 * cz)]

            # Implicit Euler: A = I - dt_Re * Delta
            A[row, w3(i, j)] = 1.0 - dt_Re * (rad_diag + z_diag)
            if j == 0:
                A[row, w3(i, 1)] = -dt_Re * rad_jp
            else:
                A[row, w3(i, j + 1)] = -dt_Re * rad_jp
                A[row, w3(i, j - 1)] = -dt_Re * rad_jm
            for (i_other, c_other) in z_extra:
                A[row, w3(i_other, j)] = -dt_Re * c_other

    # ============================================================
    # Block 2: w^r  (= w^1, vector heat eq;
    #               surface row = stress-free (free) | no-slip (body))
    # ============================================================
    for i in range(nz):
        for j in range(nr):
            row = w1(i, j)

            # Axis: w^r = 0 (regularity of an odd field)
            if j == 0:
                A[row, row] = 1.0
                continue

            # Outer wall: Dirichlet zero
            if j == nr - 1:
                A[row, row] = 1.0
                continue

            # Bottom: non-slip implicit (w^r + d_r phi = 0)
            if i == 0:
                A[row, w1(0, j)] = 1.0
                A[row, phi(0, j + 1)] = 1.0 / (2.0 * dr)
                A[row, phi(0, j - 1)] = -1.0 / (2.0 * dr)
                continue

            # Surface row
            if i == nz - 1:
                if body_mask[j]:
                    # Horizontal no-slip under body:  d_r phi + w^r = 0
                    A[row, w1(nz - 1, j)] = 1.0
                    A[row, phi(nz - 1, j + 1)] = 1.0 / (2.0 * dr)
                    A[row, phi(nz - 1, j - 1)] = -1.0 / (2.0 * dr)
                else:
                    # Free-surface stress-free (factor 2):
                    #   2 phi_rz + d_r w^z + d_z w^r = 0   (v3 fix)
                    # discretised: (w^r[-1] - w^r[-2])/dz + Dr w^z[-1]
                    #              + 2 Dr (phi_z^n) = 0
                    # => w^r[-1] - w^r[-2] + dz Dr w^z[-1] = -2 dz Dr (phi_z^n)
                    A[row, w1(nz - 1, j)] = 1.0
                    A[row, w1(nz - 2, j)] = -1.0
                    A[row, w3(nz - 1, j + 1)] = dz / (2.0 * dr)
                    A[row, w3(nz - 1, j - 1)] = -dz / (2.0 * dr)
                continue

            # Interior vector heat equation (axisymm + -w^r/r^2)
            beta = 1.0 / (2.0 * r_vec[j] * dr)
            inv_r2 = 1.0 / r_vec[j] ** 2
            rad_diag = -2.0 * cr - inv_r2
            rad_jp = cr + beta
            rad_jm = cr - beta
            z_diag = -2.0 * cz

            A[row, w1(i, j)] = 1.0 - dt_Re * (rad_diag + z_diag)
            A[row, w1(i, j + 1)] = -dt_Re * rad_jp
            A[row, w1(i, j - 1)] = -dt_Re * rad_jm
            A[row, w1(i + 1, j)] = -dt_Re * cz
            A[row, w1(i - 1, j)] = -dt_Re * cz

    # ============================================================
    # Block 3: phi  (Laplace + bottom ghost + dynamic BC with p_s)
    # ============================================================
    for i in range(nz):
        for j in range(nr):
            row = phi(i, j)

            # Outer wall: Dirichlet zero
            if j == nr - 1:
                A[row, row] = 1.0
                continue

            # Radial Laplace stencil (scaled by dz^2)
            if j == 0:
                rad_diag = -4.0 * alpha
                rad_jp = +4.0 * alpha
                rad_jm = 0.0
            else:
                beta = alpha * dr / (2.0 * r_vec[j])
                rad_diag = -2.0 * alpha
                rad_jp = alpha + beta
                rad_jm = alpha - beta

            if i == 0:
                # Bottom: Laplace with ghost-point non-slip (phi_{-1} = phi_1 - 2 dz w^z_0)
                A[row, phi(0, j)] = rad_diag - 2.0
                if j == 0:
                    A[row, phi(0, 1)] = rad_jp
                else:
                    A[row, phi(0, j + 1)] = rad_jp
                    A[row, phi(0, j - 1)] = rad_jm
                A[row, phi(1, j)] = +2.0
                A[row, w3(0, j)] = -2.0 * dz
                continue

            if i == nz - 1:
                # Surface dynamic BC (applied uniformly with p_s term);
                # the p_s row separately pins p_s = 0 outside the body.
                if j == 0:
                    A[row, phi(nz - 1, 0)] = 1.0 + two_dt_Re * 4.0 * cr
                    A[row, phi(nz - 1, 1)] = -two_dt_Re * 4.0 * cr
                else:
                    beta_c = 1.0 / (2.0 * r_vec[j] * dr)
                    A[row, phi(nz - 1, j)] = 1.0 + two_dt_Re * 2.0 * cr
                    A[row, phi(nz - 1, j + 1)] = -two_dt_Re * (cr + beta_c)
                    A[row, phi(nz - 1, j - 1)] = -two_dt_Re * (cr - beta_c)
                A[row, w3(nz - 1, j)] = +two_dt_Re / dz
                A[row, w3(nz - 2, j)] = -two_dt_Re / dz
                A[row, ps(j)] = dt
                continue

            # Interior 5-point axisymmetric Laplacian
            A[row, phi(i, j)] = rad_diag - 2.0
            if j == 0:
                A[row, phi(i, 1)] = rad_jp
            else:
                A[row, phi(i, j + 1)] = rad_jp
                A[row, phi(i, j - 1)] = rad_jm
            A[row, phi(i + 1, j)] = 1.0
            A[row, phi(i - 1, j)] = 1.0

    # ============================================================
    # Block 4: eta  (free=>kinematic; body=>eta = z_b)
    # ============================================================
    for j in range(nr):
        row = eta(j)
        if j == nr - 1:
            A[row, row] = 1.0  # outer Dirichlet
            continue
        if body_mask[j]:
            A[row, eta(j)] = 1.0
            A[row, zb] = -1.0
        else:
            # Free: eta - dt*(phi_z + w^z[-1]) = eta^n
            A[row, eta(j)] = 1.0
            A[row, phi(nz - 1, j)] = -dt / dz
            A[row, phi(nz - 2, j)] = +dt / dz
            A[row, w3(nz - 1, j)] = -dt

    # ============================================================
    # Block 5: p_s  (free=>p_s=0; body=>kinematic match)
    # ============================================================
    for j in range(nr):
        row = ps(j)
        if j == nr - 1:
            A[row, row] = 1.0  # outer Dirichlet (p_s = 0)
            continue
        if body_mask[j]:
            # Kinematic match: (phi[-1]-phi[-2])/dz + w^z[-1] - v_b = 0
            A[row, phi(nz - 1, j)] = 1.0 / dz
            A[row, phi(nz - 2, j)] = -1.0 / dz
            A[row, w3(nz - 1, j)] = 1.0
            A[row, vb] = -1.0
        else:
            A[row, ps(j)] = 1.0

    # ============================================================
    # Block 6: z_b row
    # ============================================================
    A[zb, zb] = 1.0
    A[zb, vb] = -dt

    # ============================================================
    # Block 7: v_b row  (Newton's 2nd law)
    #   M v_b - dt * 2*pi * Sum_{j in body} p_s[j] * r_j * dr
    #     = M v_b^n - dt * M / Fr
    # ============================================================
    A[vb, vb] = M
    for j in range(nr):
        if body_mask[j]:
            # r_j * 2*pi*dr is the cylindrical area weight; at j=0 it is zero
            # (which is correct: the axis carries zero area).
            A[vb, ps(j)] = -dt * 2.0 * np.pi * r_vec[j] * dr

    return A.tocsc(), off


# =============================================================================
# RHS assembly
# =============================================================================


def build_monolithic_rhs_axi_body(
    w3_n, w1_n, phi_n, eta_n, zb_n, vb_n,
    off, dr, dz, dt, Fr, We, Re, M, body_mask,
    Delta_H_axi, Dr_op,
):
    """Assemble the explicit RHS b(x^n) for the axisymmetric body solver."""
    nr = off["nr"]
    nz = off["nz"]
    N = off["total"]
    b = np.zeros(N)

    # ---- w^z block ----
    # All non-bottom interior + surface rows have RHS = w^z^n;
    # bottom is homogeneous (non-slip BC); outer wall is Dirichlet (RHS=0).
    for i in range(1, nz):
        b[off["w3"] + i * nr : off["w3"] + (i + 1) * nr] = w3_n[i, :]
        b[off["w3"] + i * nr + (nr - 1)] = 0.0  # outer Dirichlet

    # ---- w^r block ----
    # Bulk interior: RHS = w^r^n  (not bottom, not surface, not axis, not outer wall)
    for i in range(1, nz - 1):
        b[off["w1"] + i * nr : off["w1"] + (i + 1) * nr] = w1_n[i, :]
        b[off["w1"] + i * nr + 0] = 0.0           # axis
        b[off["w1"] + i * nr + (nr - 1)] = 0.0    # outer wall

    # Surface row:  free => -2 dz Dr(phi_z^n);  body => 0;  axis/outer => 0
    phi_z_n = (phi_n[-1, :] - phi_n[-2, :]) / dz
    phi_rz_n = Dr_op.dot(phi_z_n)
    rhs_w1s = np.zeros(nr)
    free_mask = ~body_mask
    free_mask[0] = False
    free_mask[-1] = False
    rhs_w1s[free_mask] = -2.0 * dz * phi_rz_n[free_mask]
    b[off["w1"] + (nz - 1) * nr : off["w1"] + nz * nr] = rhs_w1s

    # ---- phi block ----
    # Bottom/interior rows: RHS = 0 (Laplace eq).  Surface row:
    #   phi_s^n - (dt/Fr) eta^n + (dt/We) kappa[eta^n]
    # where kappa[eta] = Delta_H_axi eta.
    kappa_eta = Delta_H_axi.dot(eta_n)
    rhs_phi_s = phi_n[-1, :] - (dt / Fr) * eta_n + (dt / We) * kappa_eta
    rhs_phi_s[-1] = 0.0  # outer Dirichlet
    b[off["phi"] + (nz - 1) * nr : off["phi"] + nz * nr] = rhs_phi_s

    # ---- eta block ----
    # Free non-outer rows: RHS = eta^n;  body rows + outer: RHS = 0.
    rhs_eta = np.zeros(nr)
    rhs_eta[~body_mask] = eta_n[~body_mask]
    rhs_eta[-1] = 0.0
    b[off["eta"] : off["eta"] + nr] = rhs_eta

    # ---- p_s block: RHS = 0 ----  (free => p_s=0;  body => homogeneous match)

    # ---- z_b row ----
    b[off["zb"]] = zb_n

    # ---- v_b row ----
    b[off["vb"]] = M * vb_n - dt * M / Fr

    return b


# =============================================================================
# Pack / unpack
# =============================================================================


def pack_state(w3, w1, phi, eta, ps, zb, vb, off):
    x = np.zeros(off["total"])
    x[off["w3"] : off["w1"]] = w3.ravel()
    x[off["w1"] : off["phi"]] = w1.ravel()
    x[off["phi"] : off["eta"]] = phi.ravel()
    x[off["eta"] : off["ps"]] = eta
    x[off["ps"] : off["zb"]] = ps
    x[off["zb"]] = zb
    x[off["vb"]] = vb
    return x


def unpack_state(x, off):
    nr, nz = off["nr"], off["nz"]
    w3 = x[off["w3"] : off["w1"]].reshape((nz, nr))
    w1 = x[off["w1"] : off["phi"]].reshape((nz, nr))
    phi = x[off["phi"] : off["eta"]].reshape((nz, nr))
    eta = x[off["eta"] : off["ps"]].copy()
    ps = x[off["ps"] : off["zb"]].copy()
    zb = float(x[off["zb"]])
    vb = float(x[off["vb"]])
    return w3, w1, phi, eta, ps, zb, vb


# =============================================================================
# Solver class
# =============================================================================


class FloatingObjectSolverAxi:
    """Axisymmetric LNS solver with self-consistent rigid-body coupling.

    Body: flat-bottomed cylinder (disk of radius a) on the free surface, mass m,
    constrained to vertical motion by axisymmetry.
    """

    def __init__(self, config: SimConfig | None = None):
        if config is None:
            config = SimConfig()
        self.cfg = config
        self._setup()

    def _setup(self):
        c = self.cfg
        p, b, n = c.phys, c.body, c.num

        # Dimensionless numbers
        V = n.unit_length / n.unit_time
        self.Fr = V**2 / (p.g * n.unit_length)
        self.We = p.rho * V**2 * n.unit_length / p.sigma
        self.Re = V * n.unit_length / p.nu

        # Dimensionless grid
        self.L = p.L_cm / n.unit_length
        self.D = p.D_cm / n.unit_length
        self.nr = n.nr
        self.nz = n.nz
        self.dr, self.r_vec = make_r_grid(n.nr, self.L)
        self.dz = self.D / n.nz
        self.T_total = n.total_time_s / n.unit_time
        self.nt = n.nt
        self.dt = self.T_total / n.nt

        # Body dimensionless parameters
        self.a_obj = b.radius_cm / n.unit_length
        # 3D mass scaling:  M = m / (rho * L_c^3)
        self.M = b.mass_g / (p.rho * n.unit_length**3)
        self.zb0 = b.initial_z_cm / n.unit_length
        self.vb0 = b.initial_vz_cm_s * n.unit_time / n.unit_length

        # Body mask (sharp):  r_j < a
        self.body_mask = self.r_vec < self.a_obj
        self.free_mask = ~self.body_mask
        self.n_body = int(self.body_mask.sum())

        # 3D Archimedes equilibrium depth
        if self.a_obj > 0:
            self.zb_eq = -self.M / (np.pi * self.a_obj**2)
        else:
            self.zb_eq = float("nan")

        # Operators (used by the RHS assembly only; the matrix uses inline stencils)
        self.Delta_H_axi = build_Delta_H_axi(n.nr, self.dr, self.r_vec)
        self.Dr_op = build_Dr(n.nr, self.dr)

        # Monolithic matrix + LU factorisation
        t0 = timer.time()
        A, self.off = build_monolithic_matrix_axi_body(
            n.nr, n.nz, self.dr, self.dz, self.dt,
            self.Fr, self.We, self.Re, self.M, self.body_mask, self.r_vec,
        )
        self.A = A
        t_build = timer.time() - t0

        t0 = timer.time()
        self.LU = splu(A)
        t_fac = timer.time() - t0

        if c.verbose:
            print("FloatingObjectSolverAxi initialised")
            print(f"  Domain: r∈[0,{p.L_cm}] cm × z∈[-{p.D_cm},0] cm")
            print(f"  Grid: nr={n.nr}, nz={n.nz}")
            print(f"  Fr={self.Fr:.4f}, We={self.We:.4f}, Re={self.Re:.0f}")
            print(f"  Body: radius={b.radius_cm:.2f} cm, mass={b.mass_g:.3f} g")
            print(f"  M (dim'less 3D)        = {self.M:.4f}")
            print(f"  Predicted z_b^eq       = {self.zb_eq:.6f}"
                  f"  ({self.zb_eq * n.unit_length * 10:+.4f} mm)")
            print(f"  Body nodes: {self.n_body}/{n.nr}  (a_dim = {self.a_obj:.4f})")
            print(f"  System size N = {self.off['total']}, nnz = {A.nnz}")
            print(f"  Matrix build: {t_build:.2f}s,  LU factor: {t_fac:.2f}s")
            print(f"  dt = {self.dt * n.unit_time * 1000:.3f} ms, {n.nt} steps")

    # ------------------------------------------------------------
    # One time step
    # ------------------------------------------------------------
    def _step(self, w3, w1, phi, eta, ps, zb, vb):
        b = build_monolithic_rhs_axi_body(
            w3, w1, phi, eta, zb, vb,
            self.off, self.dr, self.dz, self.dt,
            self.Fr, self.We, self.Re, self.M,
            self.body_mask, self.Delta_H_axi, self.Dr_op,
        )
        x = self.LU.solve(b)
        return unpack_state(x, self.off)

    # ------------------------------------------------------------
    # Main run loop
    # ------------------------------------------------------------
    def run(self) -> dict:
        nr, nz, nt = self.nr, self.nz, self.nt
        save_every = self.cfg.save_every or max(1, nt // 100)
        ut = self.cfg.num.unit_time

        # Initial condition: eta = z_b under the body, zero outside
        eta = np.zeros(nr)
        eta[self.body_mask] = self.zb0
        phi = np.zeros((nz, nr))
        w3 = np.zeros((nz, nr))
        w1 = np.zeros((nz, nr))
        ps = np.zeros(nr)
        zb = self.zb0
        vb = self.vb0

        eta_hist, ps_hist = [], []
        zb_hist, vb_hist, time_hist = [], [], []
        energy_hist, mass_hist, force_hist = [], [], []

        # Cylindrical area weight  (used in diagnostics)
        area_w = 2.0 * np.pi * self.r_vec * self.dr

        if self.cfg.verbose:
            print(f"\nRunning {nt} steps...")

        t_wall = timer.time()

        for n in range(nt):
            w3, w1, phi, eta, ps, zb, vb = self._step(w3, w1, phi, eta, ps, zb, vb)

            if n % save_every == 0 or n == nt - 1:
                eta_hist.append(eta.copy())
                ps_hist.append(ps.copy())
                zb_hist.append(zb)
                vb_hist.append(vb)
                time_hist.append((n + 1) * self.dt * ut)

                # Diagnostics (cylindrical-weighted)
                # Free-surface contribution to fluid mass change
                fluid_mass = np.sum(eta * area_w)
                # Body integrated traction
                body_force = np.sum(ps[self.body_mask] * area_w[self.body_mask])
                KE_body = 0.5 * self.M * vb**2
                PE_body = (self.M / self.Fr) * zb
                fluid_PE = (1.0 / self.Fr) * np.sum(eta**2 * area_w)
                energy_hist.append(KE_body + PE_body + fluid_PE)
                mass_hist.append(fluid_mass)
                force_hist.append(body_force)

            if self.cfg.verbose and (
                n < 5 or (n + 1) % max(1, nt // 20) == 0 or n == nt - 1
            ):
                t_real = (n + 1) * self.dt * ut
                z_b_mm = zb * self.cfg.num.unit_length * 10
                v_b_cms = vb * self.cfg.num.unit_length / self.cfg.num.unit_time
                F_int = np.sum(ps[self.body_mask] * area_w[self.body_mask])
                # max | eta | outside the body
                if np.any(self.free_mask):
                    max_eta_free = np.max(np.abs(eta[self.free_mask]))
                else:
                    max_eta_free = 0.0
                print(
                    f"  Step {n+1:5d} | t={t_real:.4f}s "
                    f"| z_b={z_b_mm:+.3f}mm | v_b={v_b_cms:+.3f}cm/s "
                    f"| F={F_int:+.4e} | max|eta_free|={max_eta_free:.6f}"
                )

        wall_time = timer.time() - t_wall

        if self.cfg.verbose:
            print(f"\nDone in {wall_time:.2f}s ({wall_time/nt*1000:.2f} ms/step)")
            ul = self.cfg.num.unit_length
            print(f"  Final z_b: {zb * ul * 10:+.4f} mm")
            print(f"  Predicted z_b^eq: {self.zb_eq * ul * 10:+.4f} mm")

        return {
            "eta_history": eta_hist,
            "pressure_history": ps_hist,
            "z_b": np.array(zb_hist),
            "v_b": np.array(vb_hist),
            "time": np.array(time_hist),
            "energy": np.array(energy_hist),
            "mass": np.array(mass_hist),
            "body_force": np.array(force_hist),
            "r": self.r_vec.copy(),
            "config": self.cfg,
            "wall_time": wall_time,
            "zb_eq_predicted": self.zb_eq,
            "M": self.M,
            "a_obj": self.a_obj,
        }


# =============================================================================
# CLI
# =============================================================================


def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Axisymmetric floating-body simulation",
    )
    ap.add_argument("--nr", type=int, default=200)
    ap.add_argument("--nz", type=int, default=60)
    ap.add_argument("--nt", type=int, default=2000)
    ap.add_argument("--time", type=float, default=1.0, help="total time (s)")
    ap.add_argument("--radius", type=float, default=2.5, help="body radius (cm)")
    ap.add_argument("--mass", type=float, default=7.0, help="body mass (g)")
    ap.add_argument("--initial-z", type=float, default=0.0)
    ap.add_argument("--initial-vz", type=float, default=0.0)
    ap.add_argument("--L", type=float, default=50.0, help="domain radius L (cm)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    cfg = SimConfig(
        phys=PhysicalParams(L_cm=args.L),
        num=NumericalParams(nr=args.nr, nz=args.nz, nt=args.nt,
                            total_time_s=args.time),
        body=BodyParams(radius_cm=args.radius, mass_g=args.mass,
                        initial_z_cm=args.initial_z,
                        initial_vz_cm_s=args.initial_vz),
        verbose=not args.quiet,
    )
    solver = FloatingObjectSolverAxi(cfg)
    return solver.run()


if __name__ == "__main__":
    main()
