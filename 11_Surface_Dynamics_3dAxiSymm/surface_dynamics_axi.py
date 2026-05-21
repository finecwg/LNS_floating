"""
Axisymmetric LNS free-surface solver (pure fluid, no body).

Adapts the 2D Cartesian monolithic solver to cylindrical (r, z) coordinates
with axisymmetric flow (no swirl, no theta dependence).  Same monolithic-LU
strategy, same v3 corrections, just axisymmetric operators.

State vector  (size 3*nr*nz + nr):
    x = [w^z (nz x nr), w^1 (nz x nr), phi (nz x nr), eta (nr)]
    where w^z is the vertical vortical and w^1 = w^r is the radial vortical.

v3-style corrections:
  - drop eq (2.15e) for w^z surface; use heat eq inclusively with one-sided d_zz
  - factor of 2 on phi_rz in stress-free condition (2.11a):
        2 phi_rz + d_r w^z + d_z w^1 = 0   at z = 0, free surface

Axisymmetric specifics:
  - scalar Laplacian: d_rr + (1/r) d_r + d_zz   (for phi, w^z)
  - vector Laplacian for w^1 (radial component): adds  -w^1 / r^2
  - mirror at r=0 (axis) for even fields; w^1 = 0 there (regularity)
  - Dirichlet zero at r=L ("effectively infinite" — pick L large)
"""

from dataclasses import dataclass, field
import time as timer

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu

from operators_axi import build_Delta_H_axi, build_Dr, make_r_grid


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class PhysicalParams:
    """Dimensional physical parameters (CGS)."""
    L_cm: float = 200.0    # Domain radius (cm) - chosen large for "effectively infinite"
    D_cm: float = 20.0     # Domain depth (cm)
    g: float = 980.0       # Gravity (cm/s^2)
    sigma: float = 174.9   # Surface tension (dyne/cm)
    rho: float = 1.0       # Density (g/cm^3)
    nu: float = 0.0005     # Kinematic viscosity (cm^2/s)


@dataclass
class InitialCondition:
    """Initial perturbation at t = 0 (only Gaussian is implemented)."""
    type: str = "gaussian"
    amplitude_cm: float = 0.05
    width_cm: float = 5.0       # Gaussian half-width in r (cm), centered at r=0


@dataclass
class NumericalParams:
    """Grid and time-stepping parameters."""
    nr: int = 200
    nz: int = 80
    nt: int = 1000
    total_time_s: float = 0.5
    unit_length: float = 2.5    # non-dim length scale (cm)
    unit_time: float = 0.05     # non-dim time scale (s)


@dataclass
class SimConfig:
    phys: PhysicalParams = field(default_factory=PhysicalParams)
    ic: InitialCondition = field(default_factory=InitialCondition)
    num: NumericalParams = field(default_factory=NumericalParams)
    save_every: int | None = None
    verbose: bool = True


# =============================================================================
# Monolithic matrix assembly
# =============================================================================

def _compute_offsets(nr: int, nz: int) -> dict:
    Nb = nr * nz
    return {
        "w3": 0,
        "w1": Nb,
        "phi": 2 * Nb,
        "eta": 3 * Nb,
        "total": 3 * Nb + nr,
        "Nb": Nb,
        "nr": nr,
        "nz": nz,
    }


def build_monolithic_matrix_axi(nr, nz, dr, dz, dt, Fr, We, Re, r_vec):
    """
    Assemble the axisymmetric monolithic matrix.

    Block layout (rows in this order, same as columns):
        1) w^z (size nr*nz): bulk heat eq + bottom non-slip + surface heat eq
                             (inclusively, with one-sided d_zz at i=nz-1)
        2) w^1 (size nr*nz): bulk vector heat eq (extra -w^1/r^2)
                             + bottom non-slip + surface stress-free (factor-2)
        3) phi (size nr*nz): Laplace + bottom non-slip ghost
                             + top dynamic BC (uniform; p_s=0 since no body)
        4) eta (size nr):    kinematic BC

    Boundary conventions:
        - Axis (j=0): mirror for w^z, phi, eta;  w^1 = 0 (Dirichlet, regularity)
        - Outer (j=nr-1): Dirichlet zero for all fields ("effectively infinite")
    """
    off = _compute_offsets(nr, nz)
    N = off["total"]
    A = lil_matrix((N, N))

    cr = 1.0 / dr**2
    cz = 1.0 / dz**2
    alpha = (dz / dr)**2          # = cr * dz^2; useful for Laplace stencils
    dt_Re = dt / Re
    two_dt_Re = 2.0 * dt / Re

    def w3(i, j): return off["w3"] + i * nr + j
    def w1(i, j): return off["w1"] + i * nr + j
    def phi(i, j): return off["phi"] + i * nr + j
    def eta(j): return off["eta"] + j

    # -----------------------------------------------------------------
    # Block 1: w^z (scalar heat eq, with v3-style surface "inclusive" stencil)
    # -----------------------------------------------------------------
    # Wall (j = nr-1) treatment: zero out phi, w^z, w^1 at the wall.  Since
    # phi = 0 along the ENTIRE wall (all i), we automatically get d_z phi = 0
    # there, so u_z = d_z phi + w^z = 0 + 0 = 0.  Combined with u_r = d_r phi
    # + w^r = 0 + 0 = 0, this IS strict non-slip — provided the wall is far
    # enough that the fields are actually small there (decay to zero by
    # geometric spreading + viscous damping before reaching r = L).
    for i in range(nz):
        for j in range(nr):
            row = w3(i, j)

            # Outer wall: non-slip via Dirichlet zero (justified by far-wall limit)
            if j == nr - 1:
                A[row, row] = 1.0
                continue

            # Bottom (i = 0): non-slip implicit, w^z + (phi[1] - phi[0])/dz = 0
            if i == 0:
                A[row, w3(0, j)] = 1.0
                A[row, phi(1, j)] = 1.0 / dz
                A[row, phi(0, j)] = -1.0 / dz
                continue

            # Otherwise: heat equation row.  Build coefficients of
            # (I - dt/Re Delta_axi) where Delta_axi has scalar form.
            #
            # Radial part of Delta_axi:  d_rr + (1/r) d_r
            if j == 0:
                # Axis (L'Hopital + mirror): Delta_H f|_0 = 4(f_1 - f_0)/dr^2
                rad_diag = -4.0 * cr
                rad_jp = +4.0 * cr   # coefficient on f_{j+1} = f_1
                rad_jm = 0.0         # no f_{j-1} term
            else:
                beta = 1.0 / (2.0 * r_vec[j] * dr)
                rad_diag = -2.0 * cr
                rad_jp = cr + beta
                rad_jm = cr - beta

            # Vertical part of Delta_axi
            if i == nz - 1:
                # Surface: one-sided d_zz |_{nz-1} = (f[-1] - 2 f[-2] + f[-3])/dz^2
                z_diag = +1.0 * cz
                z_extra = [(nz - 2, -2.0 * cz), (nz - 3, +1.0 * cz)]
            else:
                # Interior: standard centered
                z_diag = -2.0 * cz
                z_extra = [(i + 1, +1.0 * cz), (i - 1, +1.0 * cz)]

            # Implicit Euler:  A = I - dt_Re * Delta
            A[row, w3(i, j)] = 1.0 - dt_Re * (rad_diag + z_diag)
            if j == 0:
                A[row, w3(i, 1)] = -dt_Re * rad_jp
            else:
                A[row, w3(i, j + 1)] = -dt_Re * rad_jp
                A[row, w3(i, j - 1)] = -dt_Re * rad_jm
            for (i_other, c_other) in z_extra:
                A[row, w3(i_other, j)] = -dt_Re * c_other

    # -----------------------------------------------------------------
    # Block 2: w^1 (= w^r, vector heat eq, surface stress-free factor-2)
    # -----------------------------------------------------------------
    for i in range(nz):
        for j in range(nr):
            row = w1(i, j)

            # Axis (j=0): w^1 = 0 by regularity of an odd field
            if j == 0:
                A[row, row] = 1.0
                continue

            # Outer wall: Dirichlet zero (non-slip in far-wall limit)
            if j == nr - 1:
                A[row, row] = 1.0
                continue

            # Bottom (i = 0): non-slip implicit, w^1 + d_r phi = 0
            if i == 0:
                A[row, w1(0, j)] = 1.0
                A[row, phi(0, j + 1)] = 1.0 / (2.0 * dr)
                A[row, phi(0, j - 1)] = -1.0 / (2.0 * dr)
                continue

            # Surface (i = nz-1): stress-free with factor 2 on phi_rz
            #   2 phi_rz + d_r w^z + d_z w^1 = 0      at z = 0, free surface
            # Discretize d_z w^1 as one-sided, d_r as centred, phi_rz at
            # time n (explicit; goes into RHS).  Coefficient layout:
            #   w^1[-1,j] - w^1[-2,j] + dz Dr(w^z[-1])_j  =  -2 dz (phi_rz)^n_j
            if i == nz - 1:
                A[row, w1(nz - 1, j)] = 1.0
                A[row, w1(nz - 2, j)] = -1.0
                # implicit dz Dr(w^z[-1])_j
                A[row, w3(nz - 1, j + 1)] = dz / (2.0 * dr)
                A[row, w3(nz - 1, j - 1)] = -dz / (2.0 * dr)
                continue

            # Interior (vector heat eq with extra -1/r^2)
            beta = 1.0 / (2.0 * r_vec[j] * dr)
            inv_r2 = 1.0 / r_vec[j]**2
            rad_diag = -2.0 * cr - inv_r2
            rad_jp = cr + beta
            rad_jm = cr - beta
            z_diag = -2.0 * cz

            A[row, w1(i, j)] = 1.0 - dt_Re * (rad_diag + z_diag)
            A[row, w1(i, j + 1)] = -dt_Re * rad_jp
            A[row, w1(i, j - 1)] = -dt_Re * rad_jm
            A[row, w1(i + 1, j)] = -dt_Re * cz
            A[row, w1(i - 1, j)] = -dt_Re * cz

    # -----------------------------------------------------------------
    # Block 3: phi  (Laplace + bottom non-slip ghost + top dynamic BC)
    # Wall: phi = 0 (Dirichlet), justified by "wall is far" — making phi=0
    # along the entire wall implies d_z phi = 0 there too, so non-slip holds.
    # -----------------------------------------------------------------
    for i in range(nz):
        for j in range(nr):
            row = phi(i, j)

            # Outer wall: Dirichlet zero
            if j == nr - 1:
                A[row, row] = 1.0
                continue

            # Radial coeffs for Laplace stencil (all multiplied by dz^2)
            if j == 0:
                # Axis (L'Hopital + mirror): radial part is 4*alpha*(f_1 - f_0)
                rad_diag = -4.0 * alpha
                rad_jp = +4.0 * alpha
                rad_jm = 0.0
            else:
                beta = alpha * dr / (2.0 * r_vec[j])     # = dz^2/(2 r_j dr)
                rad_diag = -2.0 * alpha
                rad_jp = alpha + beta
                rad_jm = alpha - beta

            if i == 0:
                # Bottom: Laplace with ghost-point non-slip
                #   ghost phi_{-1,j} = phi_{1,j} - 2 dz w^z_{0,j}
                # 5-point stencil collapses to:
                #   (rad coefs) + (z: 2 phi_{1,j} - 2 phi_{0,j}) - 2 dz w^z_{0,j} = 0
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
                # Surface: dynamic BC (2.15d) — pure fluid, p_s = 0 always.
                #   (I - 2dt/Re Delta_H_axi) phi_s + (2dt/Re/dz)(w^z[-1] - w^z[-2])
                #     = phi_s^n - (dt/Fr) eta^n + (dt/We) kappa[eta^n]
                # Δ_H_axi here is the SCALAR-valued radial operator (no dz^2).
                if j == 0:
                    A[row, phi(nz - 1, 0)] = 1.0 + two_dt_Re * 4.0 * cr
                    A[row, phi(nz - 1, 1)] = -two_dt_Re * 4.0 * cr
                else:
                    beta_c = 1.0 / (2.0 * r_vec[j] * dr)
                    A[row, phi(nz - 1, j)] = 1.0 + two_dt_Re * 2.0 * cr
                    A[row, phi(nz - 1, j + 1)] = -two_dt_Re * (cr + beta_c)
                    A[row, phi(nz - 1, j - 1)] = -two_dt_Re * (cr - beta_c)
                # w^z_z coupling
                A[row, w3(nz - 1, j)] = +two_dt_Re / dz
                A[row, w3(nz - 2, j)] = -two_dt_Re / dz
                # eta and kappa[eta] go to RHS (explicit)
                continue

            # Interior 5-point Laplace
            A[row, phi(i, j)] = rad_diag - 2.0
            if j == 0:
                A[row, phi(i, 1)] = rad_jp
            else:
                A[row, phi(i, j + 1)] = rad_jp
                A[row, phi(i, j - 1)] = rad_jm
            A[row, phi(i + 1, j)] = 1.0
            A[row, phi(i - 1, j)] = 1.0

    # -----------------------------------------------------------------
    # Block 4: eta (kinematic BC)
    # -----------------------------------------------------------------
    for j in range(nr):
        row = eta(j)
        if j == nr - 1:
            A[row, row] = 1.0
            continue
        # eta^{n+1} - dt (phi_z^{n+1} + w^z[-1]^{n+1}) = eta^n
        A[row, eta(j)] = 1.0
        A[row, phi(nz - 1, j)] = -dt / dz
        A[row, phi(nz - 2, j)] = +dt / dz
        A[row, w3(nz - 1, j)] = -dt

    return A.tocsc(), off


def build_monolithic_rhs_axi(w3_n, w1_n, phi_n, eta_n,
                             off, dr, dz, dt, Fr, We, Re,
                             Delta_H_axi, Dr_op):
    """RHS vector b(x^n) for the axisymmetric monolithic system."""
    nr, nz = off["nr"], off["nz"]
    N = off["total"]
    b = np.zeros(N)

    # -- w^z block --
    # Heat eq RHS = w^z^n at every (i, j), except bottom row (homogeneous BC).
    for i in range(nz):
        if i == 0:
            continue       # bottom: homogeneous
        b[off["w3"] + i * nr : off["w3"] + (i + 1) * nr] = w3_n[i, :]

    # Outer Dirichlet rows have RHS = 0 already.

    # -- w^1 block --
    # Bulk RHS = w^1^n; bottom = 0 (homogeneous); surface = -2 dz phi_rz^n
    for i in range(nz):
        if i == 0 or i == nz - 1:
            continue
        b[off["w1"] + i * nr : off["w1"] + (i + 1) * nr] = w1_n[i, :]
    # Surface row RHS: -2 dz Dr( phi_z^n )
    phi_z_n = (phi_n[-1, :] - phi_n[-2, :]) / dz
    phi_rz_n = Dr_op.dot(phi_z_n)
    rhs_w1s = np.zeros(nr)
    rhs_w1s[1 : nr - 1] = -2.0 * dz * phi_rz_n[1 : nr - 1]
    # j=0 and j=nr-1 rows are Dirichlet (w^1 = 0), RHS already 0.
    b[off["w1"] + (nz - 1) * nr : off["w1"] + nz * nr] = rhs_w1s

    # -- phi block --
    # Bulk Laplace + bottom: RHS = 0 except top row
    # Top row: phi_s^n - (dt/Fr) eta^n + (dt/We) kappa[eta^n]
    kappa_eta = Delta_H_axi.dot(eta_n)
    rhs_phi_s = phi_n[-1, :] - (dt / Fr) * eta_n + (dt / We) * kappa_eta
    rhs_phi_s[-1] = 0.0    # outer Dirichlet
    b[off["phi"] + (nz - 1) * nr : off["phi"] + nz * nr] = rhs_phi_s

    # -- eta block --
    rhs_eta = eta_n.copy()
    rhs_eta[-1] = 0.0
    b[off["eta"] : off["eta"] + nr] = rhs_eta

    return b


def unpack_state(x, off):
    nr, nz = off["nr"], off["nz"]
    w3 = x[off["w3"] : off["w1"]].reshape((nz, nr))
    w1 = x[off["w1"] : off["phi"]].reshape((nz, nr))
    phi = x[off["phi"] : off["eta"]].reshape((nz, nr))
    eta = x[off["eta"] :].copy()
    return w3, w1, phi, eta


# =============================================================================
# Solver class
# =============================================================================

class SurfaceDynamicsSolverAxi:
    """Axisymmetric LNS free-surface solver, no floating body."""

    def __init__(self, config: SimConfig | None = None):
        if config is None:
            config = SimConfig()
        self.cfg = config
        self._setup()

    def _setup(self):
        c = self.cfg
        p, n = c.phys, c.num

        V = n.unit_length / n.unit_time
        self.Fr = V**2 / (p.g * n.unit_length)
        self.We = p.rho * V**2 * n.unit_length / p.sigma
        self.Re = V * n.unit_length / p.nu

        self.L = p.L_cm / n.unit_length
        self.D = p.D_cm / n.unit_length
        self.nr, self.nz = n.nr, n.nz
        self.dr, self.r = make_r_grid(n.nr, self.L)
        self.dz = self.D / n.nz
        self.T_total = n.total_time_s / n.unit_time
        self.nt = n.nt
        self.dt = self.T_total / n.nt
        self.z = -self.D + np.arange(n.nz) * self.dz

        # Operators (used in RHS)
        self.Delta_H_axi = build_Delta_H_axi(n.nr, self.dr, self.r)
        self.Dr = build_Dr(n.nr, self.dr)

        # Build matrix and factor
        t0 = timer.time()
        A, self.off = build_monolithic_matrix_axi(
            n.nr, n.nz, self.dr, self.dz, self.dt,
            self.Fr, self.We, self.Re, self.r,
        )
        self.A = A
        t_build = timer.time() - t0

        t0 = timer.time()
        self.LU = splu(A)
        t_fac = timer.time() - t0

        if c.verbose:
            print(f"SurfaceDynamicsSolverAxi initialised")
            print(f"  Domain: r ∈ [0, {p.L_cm}] cm,  z ∈ [-{p.D_cm}, 0] cm")
            print(f"  Grid: nr={n.nr}, nz={n.nz}")
            print(f"  Fr={self.Fr:.4f}, We={self.We:.4f}, Re={self.Re:.0f}")
            print(f"  Monolithic size N = {self.off['total']}, nnz = {A.nnz}")
            print(f"  Matrix build: {t_build:.2f}s,  LU factor: {t_fac:.2f}s")
            print(f"  dt = {self.dt * n.unit_time * 1000:.3f} ms, {n.nt} steps")

    def _initial_eta(self):
        ic = self.cfg.ic
        amp = ic.amplitude_cm / self.cfg.num.unit_length
        width = ic.width_cm / self.cfg.num.unit_length
        eta = amp * np.exp(-(self.r / width)**2 / 2.0)
        eta[-1] = 0.0     # enforce Dirichlet at r = L
        return eta

    def _step(self, w3, w1, phi, eta):
        b = build_monolithic_rhs_axi(
            w3, w1, phi, eta,
            self.off, self.dr, self.dz, self.dt,
            self.Fr, self.We, self.Re,
            self.Delta_H_axi, self.Dr,
        )
        x = self.LU.solve(b)
        return unpack_state(x, self.off)

    def run(self):
        nr, nz, nt = self.nr, self.nz, self.nt
        save_every = self.cfg.save_every or max(1, nt // 100)
        ut = self.cfg.num.unit_time

        eta = self._initial_eta()
        phi = np.zeros((nz, nr))
        w3 = np.zeros((nz, nr))
        w1 = np.zeros((nz, nr))

        eta_hist, phi_hist, time_hist = [], [], []
        energy_hist, mass_hist = [], []

        if self.cfg.verbose:
            print(f"\nRunning {nt} steps...")
        t_wall = timer.time()

        for n in range(nt):
            w3, w1, phi, eta = self._step(w3, w1, phi, eta)

            if n % save_every == 0 or n == nt - 1:
                eta_hist.append(eta.copy())
                phi_hist.append(phi.copy())
                time_hist.append((n + 1) * self.dt * ut)
                # Axisymmetric integrals: weight by 2 pi r
                weight = 2.0 * np.pi * self.r * self.dr
                mass_hist.append(np.sum(eta * weight))
                KE = 0.5 * np.sum((phi**2).sum(axis=0) * weight) * self.dz
                PE = (1.0 / self.Fr) * np.sum(eta**2 * weight)
                energy_hist.append(KE + PE)

            if self.cfg.verbose and (n < 3 or (n + 1) % max(1, nt // 20) == 0
                                     or n == nt - 1):
                t_real = (n + 1) * self.dt * ut
                print(f"  Step {n+1:5d} | t={t_real:.4f}s "
                      f"| max|eta|={np.max(np.abs(eta)):.6f} "
                      f"| max|phi|={np.max(np.abs(phi)):.6f}")

        wall_time = timer.time() - t_wall
        if self.cfg.verbose:
            print(f"\nDone in {wall_time:.2f}s ({wall_time/nt*1000:.2f} ms/step)")

        return {
            "eta_history": eta_hist,
            "phi_history": phi_hist,
            "time": np.array(time_hist),
            "energy": np.array(energy_hist),
            "mass": np.array(mass_hist),
            "r": self.r.copy(),
            "z": self.z.copy(),
            "config": self.cfg,
            "wall_time": wall_time,
        }


# =============================================================================
# CLI / smoke test
# =============================================================================

def main():
    # Small grid smoke test
    cfg = SimConfig(
        num=NumericalParams(nr=80, nz=20, nt=50, total_time_s=0.1),
        ic=InitialCondition(amplitude_cm=0.05, width_cm=5.0),
        verbose=True,
    )
    solver = SurfaceDynamicsSolverAxi(cfg)
    result = solver.run()
    import numpy as np
    print()
    print("--- Smoke test sanity ---")
    print(f"max|eta(t=0)| = {np.max(np.abs(result['eta_history'][0])):.4e}")
    print(f"max|eta(t_end)| = {np.max(np.abs(result['eta_history'][-1])):.4e}")
    print(f"mass(t=0) = {result['mass'][0]:.4e}, mass(t_end) = {result['mass'][-1]:.4e}")
    print(f"energy(t=0) = {result['energy'][0]:.4e}, energy(t_end) = {result['energy'][-1]:.4e}")
    print(f"Wall time: {result['wall_time']:.2f} s")


if __name__ == "__main__":
    main()
