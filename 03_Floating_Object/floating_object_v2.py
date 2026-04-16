"""
Floating object v2 — self-consistent body dynamics via monolithic constraint solve.

Extends the LNS free-surface solver of Galeano-Rios, Milewski & Vanden-Broeck
(2017) to include a floating rigid body whose motion is determined self-
consistently by Newton's second law driven by gravity and the integrated fluid
pressure (not prescribed, as in v1).

Physics:
    At t=0, a flat-bottomed rigid body of half-width a and mass M (per unit
    length) is in contact with the fluid surface.  Under the body, the standard
    free-surface boundary conditions are replaced by:
        - Kinematic match:   phi_z + w^3 = v_b(t)   (no-penetration)
        - Horizontal no-slip: phi_x + w^1 = 0       (replaces stress-free 2.11a)
        - Surface follows body: eta = z_b(t)        (geometric constraint)
    The surface pressure p_s becomes a Lagrange multiplier (unknown under the
    body, zero elsewhere) enforcing the kinematic match.  The body dynamics are
        dot{z_b} = v_b,
        M dot{v_b} = -M/Fr + integral(p_s dx) over the wetted region.

Numerical method:
    Fully implicit Euler.  All equations (fluid bulk + surface + body) are
    assembled into a single monolithic linear system
        A x^{n+1} = b(x^n),
    where x = [w^3, w^1, phi, eta, p_s, z_b, v_b].  The matrix A is constant
    (depends only on grid, dt, dimensionless numbers, mass M), so it is
    factored once at setup via scipy.sparse.linalg.splu.  Each time step is
    a single back-substitution at O(N) cost.

Usage:
    from floating_object_v2 import FloatingObjectSolverV2, SimConfig, BodyParams

    config = SimConfig(body=BodyParams(mass_per_length_g_per_cm=0.5,
                                       initial_vz_cm_s=-10.0))
    solver = FloatingObjectSolverV2(config)
    result = solver.run()
    # result['z_b'], result['v_b'], result['pressure'], result['eta_history']
"""

from dataclasses import dataclass, field
import time as timer

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu


# =============================================================================
# Configuration dataclasses
# =============================================================================

@dataclass
class PhysicalParams:
    """Dimensional physical parameters (CGS)."""
    L_cm: float = 50.0     # domain length (cm)
    D_cm: float = 20.0     # domain depth (cm)
    g: float = 980.0       # gravity (cm/s^2)
    sigma: float = 174.9   # surface tension (dyne/cm)
    rho: float = 1.0       # density (g/cm^3)
    nu: float = 0.0005     # kinematic viscosity (cm^2/s)


@dataclass
class BodyParams:
    """Floating body parameters (CGS, 2D: mass per unit length)."""
    half_width_cm: float = 5.0           # half-width a (cm); total width = 2a
    center_x_cm: float = 25.0            # x-position of the body centre (cm)
    mass_per_length_g_per_cm: float = 0.5  # mass per unit length in the
                                         # out-of-plane direction (g/cm)
    initial_z_cm: float = 0.0            # initial body bottom position (cm);
                                         # 0 = at the undisturbed fluid surface
    initial_vz_cm_s: float = 0.0         # initial body vertical velocity (cm/s);
                                         # negative = downward (drop-impact)


@dataclass
class NumericalParams:
    """Grid and time-stepping parameters."""
    nx: int = 400
    nz: int = 80
    nt: int = 4000
    total_time_s: float = 2.0
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
# Operator builders (copied from v1 — unchanged)
# =============================================================================

def build_Delta_H(nx: int, dx: float):
    """Horizontal Laplacian (periodic), nx x nx sparse."""
    D = lil_matrix((nx, nx))
    c = 1.0 / dx**2
    for j in range(nx):
        D[j, j] = -2.0 * c
        D[j, (j + 1) % nx] = c
        D[j, (j - 1) % nx] = c
    return D.tocsr()


def build_Dx(nx: int, dx: float):
    """Centred first derivative in x (periodic), nx x nx sparse."""
    D = lil_matrix((nx, nx))
    c = 1.0 / (2.0 * dx)
    for j in range(nx):
        D[j, (j + 1) % nx] = c
        D[j, (j - 1) % nx] = -c
    return D.tocsr()


# =============================================================================
# Monolithic matrix assembly
# =============================================================================

def _compute_offsets(nx: int, nz: int) -> dict:
    Nb = nx * nz
    return {
        "w3": 0,
        "w1": Nb,
        "phi": 2 * Nb,
        "eta": 3 * Nb,
        "ps": 3 * Nb + nx,
        "zb": 3 * Nb + 2 * nx,
        "vb": 3 * Nb + 2 * nx + 1,
        "total": 3 * Nb + 2 * nx + 2,
        "Nb": Nb,
        "nx": nx,
        "nz": nz,
    }


def build_monolithic_matrix_v2(
    nx: int,
    nz: int,
    dx: float,
    dz: float,
    dt: float,
    Fr: float,
    We: float,
    Re: float,
    M: float,
    body_mask: np.ndarray,
):
    """
    Assemble the full monolithic matrix A for v2.

    Row blocks (in order):
        1) w^3 bulk + BCs       : nx*nz rows
        2) w^1 bulk + BCs       : nx*nz rows
        3) phi (Laplace + dyn)  : nx*nz rows
        4) eta                  : nx rows (free: 2.15c; body: eta = z_b)
        5) p_s                  : nx rows (free: 0; body: kinematic match)
        6) z_b                  : 1 row
        7) v_b                  : 1 row (Newton's 2nd law)

    Returns
    -------
    A : scipy.sparse.csc_matrix
    offsets : dict
    """
    off = _compute_offsets(nx, nz)
    N = off["total"]
    A = lil_matrix((N, N))

    # Coefficients
    cx = 1.0 / dx**2
    cz = 1.0 / dz**2
    alpha = (dz / dx) ** 2
    dt_Re = dt / Re
    two_dt_Re = 2.0 * dt / Re

    # Index helpers
    def w3(i, j):
        return off["w3"] + i * nx + j

    def w1(i, j):
        return off["w1"] + i * nx + j

    def phi(i, j):
        return off["phi"] + i * nx + j

    def eta(j):
        return off["eta"] + j

    def ps(j):
        return off["ps"] + j

    zb = off["zb"]
    vb = off["vb"]

    jl = lambda j: (j - 1) % nx
    jr = lambda j: (j + 1) % nx

    # ============================================================
    # Block 1: w^3
    # ============================================================
    for i in range(nz):
        for j in range(nx):
            row = w3(i, j)

            if i == 0:
                # Bottom non-slip: w^3[0,j] + (phi[1,j] - phi[0,j])/dz = 0
                A[row, w3(0, j)] = 1.0
                A[row, phi(1, j)] = 1.0 / dz
                A[row, phi(0, j)] = -1.0 / dz

            elif i == nz - 1:
                # Surface eq (2.15e): (I - 2dt/Re Delta_H) w^3_s = RHS (explicit phi)
                A[row, w3(nz - 1, j)] = 1.0 + two_dt_Re * 2.0 * cx
                A[row, w3(nz - 1, jr(j))] = -two_dt_Re * cx
                A[row, w3(nz - 1, jl(j))] = -two_dt_Re * cx

            else:
                # Interior heat equation
                A[row, w3(i, j)] = 1.0 + 2.0 * dt_Re * (cx + cz)
                A[row, w3(i, jr(j))] = -dt_Re * cx
                A[row, w3(i, jl(j))] = -dt_Re * cx
                A[row, w3(i + 1, j)] = -dt_Re * cz
                A[row, w3(i - 1, j)] = -dt_Re * cz

    # ============================================================
    # Block 2: w^1
    # ============================================================
    for i in range(nz):
        for j in range(nx):
            row = w1(i, j)

            if i == 0:
                # Bottom non-slip: w^1[0,j] + Dx(phi[0,:])[j] = 0
                A[row, w1(0, j)] = 1.0
                A[row, phi(0, jr(j))] = 1.0 / (2.0 * dx)
                A[row, phi(0, jl(j))] = -1.0 / (2.0 * dx)

            elif i == nz - 1:
                # Surface: stress-free (free) or no-slip (body)
                if body_mask[j]:
                    # No-slip: w^1[-1,j] + Dx(phi[-1,:])[j] = 0
                    A[row, w1(nz - 1, j)] = 1.0
                    A[row, phi(nz - 1, jr(j))] = 1.0 / (2.0 * dx)
                    A[row, phi(nz - 1, jl(j))] = -1.0 / (2.0 * dx)
                else:
                    # Stress-free (2.11a): w^1[-1] - w^1[-2] + dz Dx(w^3[-1]) = -dz phi_xz^n
                    A[row, w1(nz - 1, j)] = 1.0
                    A[row, w1(nz - 2, j)] = -1.0
                    A[row, w3(nz - 1, jr(j))] = dz / (2.0 * dx)
                    A[row, w3(nz - 1, jl(j))] = -dz / (2.0 * dx)

            else:
                # Interior heat equation
                A[row, w1(i, j)] = 1.0 + 2.0 * dt_Re * (cx + cz)
                A[row, w1(i, jr(j))] = -dt_Re * cx
                A[row, w1(i, jl(j))] = -dt_Re * cx
                A[row, w1(i + 1, j)] = -dt_Re * cz
                A[row, w1(i - 1, j)] = -dt_Re * cz

    # ============================================================
    # Block 3: phi (Laplace + surface dynamic BC)
    # ============================================================
    for i in range(nz):
        for j in range(nx):
            row = phi(i, j)

            if i == 0:
                # Bottom: Laplace with ghost point -> 2*phi(1,j) - 2*phi(0,j) - 2 dz w^3(0,j)
                A[row, phi(0, j)] = -2.0 * alpha - 2.0
                A[row, phi(0, jr(j))] = alpha
                A[row, phi(0, jl(j))] = alpha
                A[row, phi(1, j)] = 2.0
                A[row, w3(0, j)] = -2.0 * dz

            elif i == nz - 1:
                # Surface: dynamic BC (2.15d)
                #   (1 + 2 dt/Re (2/dx^2)) phi_s + 2 dt/(Re dx^2) * neighbors
                #   + (2 dt/(Re dz)) (w^3[-1] - w^3[-2]) + dt * p_s = RHS (eta^n terms)
                A[row, phi(nz - 1, j)] = 1.0 + two_dt_Re * 2.0 * cx
                A[row, phi(nz - 1, jr(j))] = -two_dt_Re * cx
                A[row, phi(nz - 1, jl(j))] = -two_dt_Re * cx
                A[row, w3(nz - 1, j)] = two_dt_Re / dz
                A[row, w3(nz - 2, j)] = -two_dt_Re / dz
                A[row, ps(j)] = dt

            else:
                # Interior: 5-point Laplacian
                A[row, phi(i, j)] = -2.0 * alpha - 2.0
                A[row, phi(i, jr(j))] = alpha
                A[row, phi(i, jl(j))] = alpha
                A[row, phi(i + 1, j)] = 1.0
                A[row, phi(i - 1, j)] = 1.0

    # ============================================================
    # Block 4: eta
    # ============================================================
    for j in range(nx):
        row = eta(j)
        if body_mask[j]:
            # Body: eta[j] - z_b = 0
            A[row, eta(j)] = 1.0
            A[row, zb] = -1.0
        else:
            # Free: (2.15c) eta - dt*(phi_z + w^3_s) = eta^n
            A[row, eta(j)] = 1.0
            A[row, phi(nz - 1, j)] = -dt / dz
            A[row, phi(nz - 2, j)] = dt / dz
            A[row, w3(nz - 1, j)] = -dt

    # ============================================================
    # Block 5: p_s
    # ============================================================
    for j in range(nx):
        row = ps(j)
        if body_mask[j]:
            # Body: kinematic match (phi[-1] - phi[-2])/dz + w^3[-1] - v_b = 0
            A[row, phi(nz - 1, j)] = 1.0 / dz
            A[row, phi(nz - 2, j)] = -1.0 / dz
            A[row, w3(nz - 1, j)] = 1.0
            A[row, vb] = -1.0
        else:
            # Free: p_s = 0
            A[row, ps(j)] = 1.0

    # ============================================================
    # Block 6: z_b
    # ============================================================
    A[zb, zb] = 1.0
    A[zb, vb] = -dt

    # ============================================================
    # Block 7: v_b (Newton's 2nd law)
    # ============================================================
    A[vb, vb] = M
    for j in range(nx):
        if body_mask[j]:
            A[vb, ps(j)] = -dt * dx

    return A.tocsc(), off


# =============================================================================
# RHS assembly
# =============================================================================

def build_monolithic_rhs_v2(
    w3_n, w1_n, phi_n, eta_n, zb_n, vb_n,
    off, dx, dz, dt, Fr, We, Re, M, body_mask, Delta_H, Dx,
):
    """Assemble the explicit RHS vector b(x^n)."""
    nx = off["nx"]
    nz = off["nz"]
    N = off["total"]
    b = np.zeros(N)
    two_dt_Re = 2.0 * dt / Re

    # ---- w^3 block ----
    # Interior and bottom/surface have different forms
    for i in range(nz):
        if i == 0:
            # Bottom: homogeneous (no explicit part)
            continue
        elif i == nz - 1:
            # Surface: w^3[-1]^n + (2dt/Re) Delta_H (phi_z^n)
            phi_z_n = (phi_n[-1, :] - phi_n[-2, :]) / dz
            rhs_w3s = w3_n[-1, :] + two_dt_Re * Delta_H.dot(phi_z_n)
            b[off["w3"] + (nz - 1) * nx : off["w3"] + nz * nx] = rhs_w3s
        else:
            # Interior: w^3^n
            b[off["w3"] + i * nx : off["w3"] + (i + 1) * nx] = w3_n[i, :]

    # ---- w^1 block ----
    for i in range(nz):
        if i == 0:
            # Bottom: homogeneous
            continue
        elif i == nz - 1:
            # Surface free: -dz * phi_xz^n; body: 0
            phi_z_n = (phi_n[-1, :] - phi_n[-2, :]) / dz
            phi_xz_n = Dx.dot(phi_z_n)
            rhs_w1s = np.zeros(nx)
            rhs_w1s[~body_mask] = -dz * phi_xz_n[~body_mask]
            b[off["w1"] + (nz - 1) * nx : off["w1"] + nz * nx] = rhs_w1s
        else:
            b[off["w1"] + i * nx : off["w1"] + (i + 1) * nx] = w1_n[i, :]

    # ---- phi block ----
    # Bottom and interior: 0; surface: phi^n - dt/Fr * eta^n + dt/We * kappa[eta^n]
    kappa = Delta_H.dot(eta_n)
    rhs_phi_s = phi_n[-1, :] - (dt / Fr) * eta_n + (dt / We) * kappa
    b[off["phi"] + (nz - 1) * nx : off["phi"] + nz * nx] = rhs_phi_s

    # ---- eta block ----
    # Free: eta^n; body: 0
    rhs_eta = np.zeros(nx)
    rhs_eta[~body_mask] = eta_n[~body_mask]
    b[off["eta"] : off["eta"] + nx] = rhs_eta

    # ---- p_s block: all zero ----

    # ---- z_b row ----
    b[off["zb"]] = zb_n

    # ---- v_b row: M*v_b^n - dt*M/Fr ----
    b[off["vb"]] = M * vb_n - dt * M / Fr

    return b


# =============================================================================
# Pack / unpack state
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
    nx, nz = off["nx"], off["nz"]
    w3 = x[off["w3"] : off["w1"]].reshape((nz, nx))
    w1 = x[off["w1"] : off["phi"]].reshape((nz, nx))
    phi = x[off["phi"] : off["eta"]].reshape((nz, nx))
    eta = x[off["eta"] : off["ps"]].copy()
    ps = x[off["ps"] : off["zb"]].copy()
    zb = float(x[off["zb"]])
    vb = float(x[off["vb"]])
    return w3, w1, phi, eta, ps, zb, vb


# =============================================================================
# Solver class
# =============================================================================

class FloatingObjectSolverV2:
    """Monolithic LNS solver with self-consistent body dynamics."""

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
        self.nx = n.nx
        self.nz = n.nz
        self.dx = self.L / n.nx
        self.dz = self.D / n.nz
        self.T_total = n.total_time_s / n.unit_time
        self.nt = n.nt
        self.dt = self.T_total / n.nt

        self.x = np.linspace(0, self.L, n.nx, endpoint=False)

        # Body dimensionless parameters
        self.a_obj = b.half_width_cm / n.unit_length
        self.xc_obj = b.center_x_cm / n.unit_length
        # Dimensionless mass M = m / (rho * L^2)
        self.M = b.mass_per_length_g_per_cm / (p.rho * n.unit_length**2)
        self.zb0 = b.initial_z_cm / n.unit_length
        self.vb0 = b.initial_vz_cm_s * n.unit_time / n.unit_length

        # Body mask
        self.body_mask = np.abs(self.x - self.xc_obj) < self.a_obj
        self.free_mask = ~self.body_mask
        self.n_body = int(self.body_mask.sum())

        # Operators (for RHS)
        self.Delta_H = build_Delta_H(n.nx, self.dx)
        self.Dx = build_Dx(n.nx, self.dx)

        # Monolithic matrix + factorisation
        t0 = timer.time()
        A, self.off = build_monolithic_matrix_v2(
            n.nx, n.nz, self.dx, self.dz, self.dt,
            self.Fr, self.We, self.Re, self.M, self.body_mask,
        )
        self.A = A
        t_build = timer.time() - t0

        t0 = timer.time()
        self.LU = splu(A)
        t_fac = timer.time() - t0

        if c.verbose:
            print("FloatingObjectSolverV2 initialised")
            print(f"  Domain: {p.L_cm}x{p.D_cm} cm, Grid: {n.nx}x{n.nz}")
            print(f"  Fr={self.Fr:.4f}, We={self.We:.4f}, Re={self.Re:.0f}")
            print(f"  Body: 2a={2*b.half_width_cm:.1f} cm, M_dim={self.M:.4f}")
            print(f"  Initial z_b={b.initial_z_cm:.3f} cm, v_b={b.initial_vz_cm_s:.2f} cm/s")
            print(f"  Body nodes: {self.n_body}/{n.nx}")
            print(f"  Monolithic size N = {self.off['total']}, nnz = {A.nnz}")
            print(f"  Matrix build: {t_build:.2f}s,  LU factor: {t_fac:.2f}s")
            print(f"  dt = {self.dt * n.unit_time * 1000:.3f} ms, {n.nt} steps")

    # --------------------------------------------------------------
    # Single time step
    # --------------------------------------------------------------
    def _step(self, w3, w1, phi, eta, ps, zb, vb):
        b = build_monolithic_rhs_v2(
            w3, w1, phi, eta, zb, vb,
            self.off, self.dx, self.dz, self.dt,
            self.Fr, self.We, self.Re, self.M,
            self.body_mask, self.Delta_H, self.Dx,
        )
        x = self.LU.solve(b)
        return unpack_state(x, self.off)

    # --------------------------------------------------------------
    # Main run
    # --------------------------------------------------------------
    def run(self) -> dict:
        nx, nz, nt = self.nx, self.nz, self.nt
        save_every = self.cfg.save_every or max(1, nt // 100)
        ut = self.cfg.num.unit_time

        # Initialise
        eta = np.zeros(nx)
        eta[self.body_mask] = self.zb0
        phi = np.zeros((nz, nx))
        w3 = np.zeros((nz, nx))
        w1 = np.zeros((nz, nx))
        ps = np.zeros(nx)
        zb = self.zb0
        vb = self.vb0

        # Histories
        eta_hist, phi_hist, ps_hist = [], [], []
        zb_hist, vb_hist, time_hist = [], [], []
        energy_hist, mass_hist, force_hist = [], [], []

        if self.cfg.verbose:
            print(f"\nRunning {nt} steps...")

        t_wall = timer.time()

        for n in range(nt):
            w3, w1, phi, eta, ps, zb, vb = self._step(w3, w1, phi, eta, ps, zb, vb)

            if n % save_every == 0 or n == nt - 1:
                eta_hist.append(eta.copy())
                phi_hist.append(phi.copy())
                ps_hist.append(ps.copy())
                zb_hist.append(zb)
                vb_hist.append(vb)
                time_hist.append((n + 1) * self.dt * ut)

                # Diagnostics
                KE_fluid = 0.5 * np.sum(phi**2) * self.dx * self.dz
                PE_fluid = (1.0 / self.Fr) * np.sum(eta**2) * self.dx
                KE_body = 0.5 * self.M * vb**2
                PE_body = (self.M / self.Fr) * zb
                energy_hist.append(KE_fluid + PE_fluid + KE_body + PE_body)
                mass_hist.append(np.sum(eta) * self.dx)
                force_hist.append(np.sum(ps[self.body_mask]) * self.dx)

            if self.cfg.verbose and (n < 5 or (n + 1) % max(1, nt // 20) == 0 or n == nt - 1):
                t_real = (n + 1) * self.dt * ut
                z_b_cm = zb * self.cfg.num.unit_length
                v_b_cms = vb * self.cfg.num.unit_length / self.cfg.num.unit_time
                F_int = np.sum(ps[self.body_mask]) * self.dx
                print(
                    f"  Step {n+1:5d} | t={t_real:.4f}s "
                    f"| z_b={z_b_cm*10:+.3f}mm | v_b={v_b_cms:+.3f}cm/s "
                    f"| F={F_int:+.4e} | max|eta_free|={np.max(np.abs(eta[self.free_mask])):.6f}"
                )

        wall_time = timer.time() - t_wall

        if self.cfg.verbose:
            print(f"\nDone in {wall_time:.2f}s ({wall_time/nt*1000:.2f} ms/step)")

        return {
            "eta_history": eta_hist,
            "phi_history": phi_hist,
            "pressure_history": ps_hist,
            "z_b": np.array(zb_hist),
            "v_b": np.array(vb_hist),
            "time": np.array(time_hist),
            "energy": np.array(energy_hist),
            "mass": np.array(mass_hist),
            "body_force": np.array(force_hist),
            "x": self.x.copy(),
            "config": self.cfg,
            "wall_time": wall_time,
        }


# =============================================================================
# CLI
# =============================================================================

def main():
    import argparse

    ap = argparse.ArgumentParser(description="Floating object v2 simulation")
    ap.add_argument("--nx", type=int, default=400)
    ap.add_argument("--nz", type=int, default=80)
    ap.add_argument("--nt", type=int, default=4000)
    ap.add_argument("--time", type=float, default=2.0, help="total time (s)")
    ap.add_argument("--mass", type=float, default=0.5,
                    help="body mass per unit length (g/cm)")
    ap.add_argument("--initial-z", type=float, default=0.0,
                    help="initial body position (cm)")
    ap.add_argument("--initial-vz", type=float, default=0.0,
                    help="initial body velocity (cm/s); negative = downward")
    ap.add_argument("--half-width", type=float, default=5.0,
                    help="body half-width (cm)")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--save", type=str, default=None)
    args = ap.parse_args()

    cfg = SimConfig(
        num=NumericalParams(nx=args.nx, nz=args.nz, nt=args.nt, total_time_s=args.time),
        body=BodyParams(
            mass_per_length_g_per_cm=args.mass,
            initial_z_cm=args.initial_z,
            initial_vz_cm_s=args.initial_vz,
            half_width_cm=args.half_width,
        ),
        verbose=not args.quiet,
    )

    solver = FloatingObjectSolverV2(cfg)
    result = solver.run()

    z_b_final_cm = result["z_b"][-1] * cfg.num.unit_length
    v_b_final_cms = result["v_b"][-1] * cfg.num.unit_length / cfg.num.unit_time

    print("\n===== Summary =====")
    print(f"  Wall time: {result['wall_time']:.2f}s")
    print(f"  Final z_b: {z_b_final_cm*10:.3f} mm")
    print(f"  Final v_b: {v_b_final_cms:.3f} cm/s")
    print(f"  Min z_b: {np.min(result['z_b']) * cfg.num.unit_length * 10:.3f} mm")
    print(f"  Max z_b: {np.max(result['z_b']) * cfg.num.unit_length * 10:.3f} mm")
    print(f"  Max |v_b|: "
          f"{np.max(np.abs(result['v_b'])) * cfg.num.unit_length / cfg.num.unit_time:.3f} cm/s")

    if args.save:
        np.savez(
            args.save,
            eta_history=np.array(result["eta_history"]),
            pressure_history=np.array(result["pressure_history"]),
            z_b=result["z_b"],
            v_b=result["v_b"],
            time=result["time"],
            energy=result["energy"],
            body_force=result["body_force"],
            x=result["x"],
        )
        print(f"  Saved to {args.save}")

    return result


if __name__ == "__main__":
    main()
