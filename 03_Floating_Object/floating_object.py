"""
Floating object on a free surface — LNS solver via Helmholtz decomposition.

Based on Galeano-Rios, Milewski & Vanden-Broeck (2017), JFM 826, 97-127.
Extends the 02_Surface_Dynamics_updated solver to include a vertically
oscillating flat-bottomed object sitting on the free surface.

Physics:
    The object is a flat-bottomed body of half-width R, oscillating vertically
    with prescribed motion z_obj(t) = -draft + A sin(ωt). Under the object,
    the kinematic BC changes: instead of η being free, the surface velocity
    is prescribed to match the object's velocity. The pressure p_s under the
    object is then whatever is needed to enforce this (computed from the
    dynamic BC). Outside the object, p_s = 0 (free surface).

Two approaches are implemented:
    1. "velocity" (default): Under the object, prescribe the surface vertical
       velocity (φ_z + w³ = v_obj). η under the object follows from integrating
       this. The dynamic BC (2.15d) then yields the pressure p_s implicitly.
    2. "displacement": Directly set η = z_obj under the object after each step
       (the original notebook approach — less physical but simpler).

Usage:
    from floating_object import FloatingObjectSolver, SimConfig
    config = SimConfig()  # or modify parameters
    solver = FloatingObjectSolver(config)
    result = solver.run()
    # result is a dict with histories: eta, phi, time, object_z, pressure, energy, mass
"""

from dataclasses import dataclass, field
import numpy as np
from scipy.sparse import lil_matrix, eye
from scipy.sparse.linalg import splu
import time as timer


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class PhysicalParams:
    """Dimensional physical parameters (CGS)."""
    L_cm: float = 50.0        # domain length (cm)
    D_cm: float = 20.0        # domain depth (cm)
    g: float = 980.0           # gravity (cm/s²)
    sigma: float = 174.9       # surface tension (dyne/cm)
    rho: float = 1.0           # density (g/cm³)
    nu: float = 0.0005         # kinematic viscosity (cm²/s)


@dataclass
class ObjectParams:
    """Floating object parameters (CGS)."""
    half_width_cm: float = 5.0       # half-width R (cm)
    center_x_cm: float = 25.0        # x-position (cm), default = domain center
    amplitude_cm: float = 0.3        # vertical oscillation amplitude A (cm)
    frequency: float = 25.0          # angular frequency ω (rad/s)
    draft_cm: float = 0.3            # equilibrium submergence (cm)


@dataclass
class NumericalParams:
    """Grid and time-stepping parameters."""
    nx: int = 400
    nz: int = 80
    nt: int = 4000
    total_time_s: float = 2.0
    # Non-dimensionalisation scales
    unit_length: float = 2.5     # cm
    unit_time: float = 0.05      # s


@dataclass
class SimConfig:
    """Full simulation configuration."""
    phys: PhysicalParams = field(default_factory=PhysicalParams)
    obj: ObjectParams = field(default_factory=ObjectParams)
    num: NumericalParams = field(default_factory=NumericalParams)
    method: str = "velocity"     # "velocity" or "displacement"
    save_every: int | None = None  # save every N steps; None → auto (nt//100)
    verbose: bool = True


# ---------------------------------------------------------------------------
# Operator builders
# ---------------------------------------------------------------------------

def build_Delta_H(nx: int, dx: float):
    """Horizontal Laplacian (periodic BC), nx × nx sparse matrix."""
    D = lil_matrix((nx, nx))
    c = 1.0 / dx**2
    for j in range(nx):
        D[j, j] = -2.0 * c
        D[j, (j + 1) % nx] = c
        D[j, (j - 1) % nx] = c
    return D.tocsr()


def build_Dx(nx: int, dx: float):
    """Centred first derivative in x (periodic BC), nx × nx sparse matrix."""
    D = lil_matrix((nx, nx))
    c = 1.0 / (2.0 * dx)
    for j in range(nx):
        D[j, (j + 1) % nx] = c
        D[j, (j - 1) % nx] = -c
    return D.tocsr()


def build_Laplace_matrix(nx: int, nz: int, dx: float, dz: float):
    """
    Laplace operator Δφ = 0 on (nx × nz) grid.
    Bottom (i=0): non-slip ghost-point (∂φ/∂z handled via RHS).
    Top (i=nz-1): Dirichlet (handled via RHS).
    Periodic in x.
    """
    N = nx * nz
    A = lil_matrix((N, N))
    alpha = (dz / dx) ** 2

    for i in range(nz):
        for j in range(nx):
            idx = i * nx + j
            A[idx, idx] = -2.0 * alpha - 2.0
            A[idx, i * nx + (j + 1) % nx] = alpha
            A[idx, i * nx + (j - 1) % nx] = alpha

            if i > 0:
                A[idx, (i - 1) * nx + j] = 1.0
            else:
                # bottom ghost: φ_{-1} = φ_1 - 2dz·w³_bot → doubles coeff for i+1
                A[idx, (i + 1) * nx + j] = 2.0

            if i < nz - 1:
                A[idx, (i + 1) * nx + j] = 1.0
            # top: Dirichlet goes into RHS

    return A.tocsr()


def build_heat_matrix_2d(nx: int, nz: int, dx: float, dz: float):
    """
    2D Laplacian for heat equation on (nx × nz) grid.
    Neumann at top/bottom (ghost-point), periodic in x.
    """
    N = nx * nz
    A = lil_matrix((N, N))
    cx = 1.0 / dx**2
    cz = 1.0 / dz**2

    for i in range(nz):
        for j in range(nx):
            idx = i * nx + j
            A[idx, idx] = -2.0 * cx - 2.0 * cz
            A[idx, i * nx + (j + 1) % nx] = cx
            A[idx, i * nx + (j - 1) % nx] = cx

            if i > 0:
                A[idx, (i - 1) * nx + j] = cz
            else:
                A[idx, (i + 1) * nx + j] += cz

            if i < nz - 1:
                A[idx, (i + 1) * nx + j] = cz
            else:
                A[idx, (i - 1) * nx + j] += cz

    return A.tocsr()


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

class FloatingObjectSolver:
    """
    LNS free-surface solver with a floating object.

    Fields solved: η (surface elevation), φ (velocity potential, 2D),
    w³ and w¹ (vortical components, 2D).
    """

    def __init__(self, config: SimConfig | None = None):
        if config is None:
            config = SimConfig()
        self.cfg = config
        self._setup()

    # ---- setup ----

    def _setup(self):
        c = self.cfg
        p, o, n = c.phys, c.obj, c.num

        # Dimensionless parameters
        V = n.unit_length / n.unit_time
        self.Fr = V**2 / (p.g * n.unit_length)
        self.We = p.rho * V**2 * n.unit_length / p.sigma
        self.Re = V * n.unit_length / p.nu

        # Dimensionless grid
        self.L = p.L_cm / n.unit_length
        self.D = p.D_cm / n.unit_length
        self.nx, self.nz = n.nx, n.nz
        self.dx = self.L / n.nx
        self.dz = self.D / n.nz
        self.T_total = n.total_time_s / n.unit_time
        self.nt = n.nt
        self.dt = self.T_total / n.nt

        # Grid vectors
        self.x = np.linspace(0, self.L, n.nx, endpoint=False)

        # Object (dimensionless)
        self.R_obj = o.half_width_cm / n.unit_length
        self.xc_obj = o.center_x_cm / n.unit_length
        self.A_obj = o.amplitude_cm / n.unit_length
        self.omega_obj = o.frequency * n.unit_time
        self.z_draft = o.draft_cm / n.unit_length

        # Object mask: boolean array, True where object sits
        self.obj_mask = np.abs(self.x - self.xc_obj) < self.R_obj
        self.free_mask = ~self.obj_mask

        # Operators
        self.Delta_H = build_Delta_H(n.nx, self.dx)
        self.Dx = build_Dx(n.nx, self.dx)
        self.A_laplace = build_Laplace_matrix(n.nx, n.nz, self.dx, self.dz)
        self.Lap_2d = build_heat_matrix_2d(n.nx, n.nz, self.dx, self.dz)

        I1 = eye(n.nx, format='csr')
        I2 = eye(n.nx * n.nz, format='csr')

        self.A_phi_surf = I1 - (2.0 * self.dt / self.Re) * self.Delta_H
        self.A_w3_surf = I1 - (2.0 * self.dt / self.Re) * self.Delta_H
        self.A_heat = I2 - (self.dt / self.Re) * self.Lap_2d

        # Prefactorise all constant matrices (factor once, backsolve each step)
        t_fac = timer.time()
        self.LU_heat = splu(self.A_heat.tocsc())       # for w3 and w1 bulk
        self.LU_w3_surf = splu(self.A_w3_surf.tocsc()) # for w3 surface BC
        self.LU_phi_surf = splu(self.A_phi_surf.tocsc()) # for phi surface BC
        self.LU_laplace = splu(self.A_laplace.tocsc())  # for bulk Laplace
        t_fac = timer.time() - t_fac

        # Precompute index arrays for vectorised Laplace RHS
        self._lap_top_idx = np.arange(n.nx) + (n.nz - 1) * n.nx
        self._lap_bot_idx = np.arange(n.nx)

        # Precompute the gravity wavelength for info
        omega_dim = o.frequency
        self.wavelength_cm = 2.0 * np.pi * p.g / omega_dim**2

        if c.verbose:
            print(f"FloatingObjectSolver initialised")
            print(f"  Domain: {p.L_cm}×{p.D_cm} cm, Grid: {n.nx}×{n.nz}")
            print(f"  Fr={self.Fr:.4f}, We={self.We:.4f}, Re={self.Re:.0f}")
            print(f"  Object: {2*o.half_width_cm}cm wide, A={o.amplitude_cm*10:.1f}mm, "
                  f"ω={o.frequency} rad/s")
            print(f"  Gravity wavelength λ ≈ {self.wavelength_cm:.1f} cm")
            print(f"  Method: {c.method}")
            print(f"  dt = {self.dt * n.unit_time * 1000:.3f} ms, {n.nt} steps")
            print(f"  LU prefactorisation: {t_fac:.3f}s")

    # ---- object kinematics ----

    def object_z(self, t: float) -> float:
        """Object bottom position (dimensionless)."""
        return -self.z_draft + self.A_obj * np.sin(self.omega_obj * t)

    def object_vz(self, t: float) -> float:
        """Object vertical velocity (dimensionless)."""
        return self.A_obj * self.omega_obj * np.cos(self.omega_obj * t)

    # ---- curvature ----

    def curvature(self, eta: np.ndarray) -> np.ndarray:
        """Linearised curvature κ[η] = Δ_H η, using the prebuilt matrix."""
        return self.Delta_H.dot(eta)

    # ---- time stepping ----

    def _step(
        self,
        eta: np.ndarray,
        phi: np.ndarray,
        w3: np.ndarray,
        w1: np.ndarray,
        t: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Advance one time step.  Returns (eta_new, phi_new, w3_new, w1_new, p_s).

        p_s is the surface pressure (nonzero under the object, zero outside).
        """
        dt = self.dt
        dx, dz = self.dx, self.dz
        nx, nz = self.nx, self.nz
        Fr, We, Re = self.Fr, self.We, self.Re
        obj = self.obj_mask
        free = self.free_mask

        t_next = t + dt
        vz_obj = self.object_vz(t)
        z_obj_next = self.object_z(t_next)

        # ============================================================
        # Step 1: Evolve w³, w¹ in bulk via heat equation (2.15b)
        # ============================================================
        w3_new = self.LU_heat.solve(w3.ravel()).reshape((nz, nx))
        w1_new = self.LU_heat.solve(w1.ravel()).reshape((nz, nx))

        # -- w³ surface BC (eq 2.15e) --
        phi_z_surf_old = (phi[-1, :] - phi[-2, :]) / dz
        rhs_w3s = w3[-1, :] + (2.0 * dt / Re) * self.Delta_H.dot(phi_z_surf_old)
        w3_new[-1, :] = self.LU_w3_surf.solve(rhs_w3s)

        # -- w³ bottom BC (non-slip) --
        w3_new[0, :] = -(phi[1, :] - phi[0, :]) / dz

        # -- w¹ bottom BC (non-slip) --
        w1_new[0, :] = -self.Dx.dot(phi[0, :])

        # -- w¹ surface BC (stress-free, eq 2.11a) --
        phi_xz_s = self.Dx.dot((phi[-1, :] - phi[-2, :]) / dz)
        w3_x_s = self.Dx.dot(w3_new[-1, :])
        w1_new[-1, :] = w1_new[-2, :] - dz * (phi_xz_s + w3_x_s)

        # ============================================================
        # Step 2: Update φ at surface (eq 2.15d)
        # ============================================================
        kappa = self.curvature(eta)
        phi_surf_old = phi[-1, :]
        w3z_s = (w3_new[-1, :] - w3_new[-2, :]) / dz

        # RHS of dynamic BC with p_s = 0 everywhere first
        rhs_phi = (
            phi_surf_old
            - (dt / Fr) * eta
            + (dt / We) * kappa
            - (2.0 * dt / Re) * w3z_s
        )
        phi_surf_new = self.LU_phi_surf.solve(rhs_phi)

        # ============================================================
        # Step 3: Solve Laplace for bulk φ
        # ============================================================
        b_lap = np.zeros(nz * nx)
        b_lap[self._lap_top_idx] = -phi_surf_new
        b_lap[self._lap_bot_idx] += -2.0 * dz * w3_new[0, :]

        phi_new = self.LU_laplace.solve(b_lap).reshape((nz, nx))

        # ============================================================
        # Step 4: Update η (eq 2.15c)
        # ============================================================
        phi_z_surf_new = (phi_new[-1, :] - phi_new[-2, :]) / dz
        w3_surf_new = w3_new[-1, :]

        if self.cfg.method == "velocity":
            # Free surface: normal evolution
            eta_new = eta.copy()
            eta_new[free] = eta[free] + dt * (phi_z_surf_new[free] + w3_surf_new[free])

            # Under object: prescribe vertical velocity.
            # The fluid surface velocity must match the object velocity.
            # η_new = η + dt * v_obj  (instead of η + dt*(φ_z + w³))
            eta_new[obj] = eta[obj] + dt * vz_obj

            # Compute the pressure that was "needed" under the object.
            # From (2.15d) rearranged: p_s = -(φ_t) - (1/Fr)η + (1/We)κ + (2/Re)Δ_H φ - (2/Re)w³_z
            # The φ_t we computed assumed p_s=0.  The actual φ_t under the object would
            # differ by p_s*dt in the implicit scheme.  So:
            # p_s ≈ (phi_surf_new_free - phi_surf_if_constrained) / dt
            # For now, we approximate p_s from the mismatch in vertical velocity:
            # Under the object, the "free" φ_z+w³ differs from v_obj.
            # The pressure adjusts φ to make the velocities match.
            free_vz = phi_z_surf_new[obj] + w3_surf_new[obj]
            p_s = np.zeros(nx)
            # p_s ~ (free_vz - v_obj) * (something involving the Laplace response)
            # This is a first-order diagnostic — the actual pressure requires
            # solving the constrained problem. Stored for inspection.
            p_s[obj] = (free_vz - vz_obj) / dt

        elif self.cfg.method == "displacement":
            # Original notebook approach: free surface evolves normally,
            # then overwrite η under the object.
            eta_new = eta + dt * (phi_z_surf_new + w3_surf_new)
            eta_new[obj] = z_obj_next
            p_s = np.zeros(nx)

        else:
            raise ValueError(f"Unknown method: {self.cfg.method!r}")

        return eta_new, phi_new, w3_new, w1_new, p_s

    # ---- main run loop ----

    def run(self) -> dict:
        """
        Run the full simulation.

        Returns a dict with keys:
            eta_history, phi_history, time, object_z, pressure,
            energy, mass, max_eta, config, wall_time
        """
        nx, nz, nt = self.nx, self.nz, self.nt
        dt = self.dt
        ut = self.cfg.num.unit_time

        save_every = self.cfg.save_every or max(1, nt // 100)

        # Initialise fields
        eta = np.zeros(nx)
        phi = np.zeros((nz, nx))
        w3 = np.zeros((nz, nx))
        w1 = np.zeros((nz, nx))

        # Set initial η under object
        z0 = self.object_z(0.0)
        eta[self.obj_mask] = z0

        # Storage
        eta_hist: list[np.ndarray] = []
        phi_hist: list[np.ndarray] = []
        time_hist: list[float] = []
        objz_hist: list[float] = []
        ps_hist: list[np.ndarray] = []
        energy_hist: list[float] = []
        mass_hist: list[float] = []
        max_eta_hist: list[float] = []

        if self.cfg.verbose:
            print(f"\nRunning {nt} steps...")

        t_wall_start = timer.time()

        for n in range(nt):
            t = n * dt
            eta, phi, w3, w1, p_s = self._step(eta, phi, w3, w1, t)

            if n % save_every == 0 or n == nt - 1:
                t_real = (n + 1) * dt * ut
                eta_hist.append(eta.copy())
                phi_hist.append(phi.copy())
                time_hist.append(t_real)
                objz_hist.append(self.object_z((n + 1) * dt))
                ps_hist.append(p_s.copy())

                # Diagnostics
                KE = 0.5 * np.sum(phi**2) * self.dx * self.dz
                PE = (1.0 / self.Fr) * np.sum(eta**2) * self.dx
                energy_hist.append(KE + PE)
                mass_hist.append(np.sum(eta) * self.dx)

                # Max |η| outside the object region
                max_eta_hist.append(np.max(np.abs(eta[self.free_mask])))

            if self.cfg.verbose:
                if n < 10 or (n + 1) % (nt // 20) == 0 or n == nt - 1:
                    t_real = (n + 1) * dt * ut
                    print(
                        f"  Step {n+1:5d} | t={t_real:.4f}s "
                        f"| max|η_free|={np.max(np.abs(eta[self.free_mask])):.6f} "
                        f"| max|φ|={np.max(np.abs(phi)):.6f}"
                    )

        wall_time = timer.time() - t_wall_start

        if self.cfg.verbose:
            print(f"\nDone in {wall_time:.1f}s ({wall_time/nt*1000:.1f} ms/step)")

        return {
            "eta_history": eta_hist,
            "phi_history": phi_hist,
            "time": np.array(time_hist),
            "object_z": np.array(objz_hist),
            "pressure": ps_hist,
            "energy": np.array(energy_hist),
            "mass": np.array(mass_hist),
            "max_eta": np.array(max_eta_hist),
            "x": self.x.copy(),
            "config": self.cfg,
            "wall_time": wall_time,
        }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    """Run with default parameters and print summary."""
    import argparse

    parser = argparse.ArgumentParser(description="Floating object LNS simulation")
    parser.add_argument("--nx", type=int, default=400)
    parser.add_argument("--nz", type=int, default=80)
    parser.add_argument("--nt", type=int, default=4000)
    parser.add_argument("--time", type=float, default=2.0, help="total time (s)")
    parser.add_argument("--amplitude", type=float, default=0.3, help="oscillation amplitude (cm)")
    parser.add_argument("--frequency", type=float, default=25.0, help="oscillation freq (rad/s)")
    parser.add_argument("--method", choices=["velocity", "displacement"], default="velocity")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--save", type=str, default=None, help="save results to .npz")
    args = parser.parse_args()

    config = SimConfig(
        num=NumericalParams(nx=args.nx, nz=args.nz, nt=args.nt, total_time_s=args.time),
        obj=ObjectParams(amplitude_cm=args.amplitude, frequency=args.frequency),
        method=args.method,
        verbose=not args.quiet,
    )

    solver = FloatingObjectSolver(config)
    result = solver.run()

    # Print summary
    print(f"\n===== Summary =====")
    print(f"  Method: {args.method}")
    print(f"  Wall time: {result['wall_time']:.1f}s")
    print(f"  Max |η| (free surface): {np.max(result['max_eta']):.6f} "
          f"({np.max(result['max_eta']) * config.num.unit_length * 10:.3f} mm)")
    print(f"  Final mass ∫η dx: {result['mass'][-1]:.6e}")
    print(f"  Energy ratio (final/initial): "
          f"{result['energy'][-1] / max(result['energy'][0], 1e-30):.4f}")

    if args.save:
        np.savez(
            args.save,
            eta_history=np.array(result["eta_history"]),
            time=result["time"],
            object_z=result["object_z"],
            energy=result["energy"],
            mass=result["mass"],
            x=result["x"],
        )
        print(f"  Saved to {args.save}")

    return result


if __name__ == "__main__":
    main()
