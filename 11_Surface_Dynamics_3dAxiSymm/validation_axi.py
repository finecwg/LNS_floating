"""
Quantitative validation for the axisymmetric LNS solver.

Two tests:
  (1) Mass-drift convergence: integrate eta * 2*pi*r dr over time, check that
      drift goes as O(dr^2) for grid refinement.
  (2) Bessel-mode benchmark: initialise eta(r,0) = J_0(k r) where k is chosen
      to be a Bessel zero of J_0 at r=L (so the IC respects the Dirichlet BC).
      Measure oscillation frequency at the axis and compare with the linearised
      gravity-capillary dispersion relation.  Measure decay rate and compare
      with the classical Lamb-1932 result gamma ≈ 2 k^2 / Re.
"""

import numpy as np
from scipy.special import j0, jn_zeros
from scipy.optimize import curve_fit

from surface_dynamics_axi import (
    SurfaceDynamicsSolverAxi, SimConfig,
    PhysicalParams, NumericalParams, InitialCondition,
)


# =============================================================================
# Test 1: Mass-drift convergence
# =============================================================================

def mass_drift(nr, nz, nt=200, total_time_s=0.2, verbose=False):
    """Run a Gaussian-bump simulation and return |Delta mass / mass(0)|."""
    cfg = SimConfig(
        num=NumericalParams(nr=nr, nz=nz, nt=nt, total_time_s=total_time_s),
        ic=InitialCondition(amplitude_cm=0.05, width_cm=5.0),
        verbose=False,
    )
    solver = SurfaceDynamicsSolverAxi(cfg)
    result = solver.run()
    m0 = result["mass"][0]
    m_final = result["mass"][-1]
    drift = abs(m_final - m0) / abs(m0)
    if verbose:
        print(f"  nr={nr}, nz={nz}: mass(0)={m0:.4e}, mass(end)={m_final:.4e}, "
              f"drift={drift:.2%}")
    return drift, m0, m_final


def run_grid_convergence():
    print("=" * 60)
    print("Test 1: Mass-drift convergence under grid refinement")
    print("=" * 60)
    print(f"{'nr':>6} {'nz':>6}   {'mass(0)':>14}   {'mass(end)':>14}   {'drift':>10}")
    print("-" * 60)
    drifts = []
    for nr, nz in [(80, 20), (160, 40), (320, 80)]:
        d, m0, mf = mass_drift(nr, nz)
        drifts.append((nr, nz, d, m0, mf))
        print(f"{nr:>6} {nz:>6}   {m0:>14.4e}   {mf:>14.4e}   {d:>10.2%}")
    print()
    if drifts[2][2] < 0.01:
        print(f"  ✓ Drift at nr=320: {drifts[2][2]:.3%}  (below 1%)")
    else:
        print(f"  ⚠ Drift at nr=320 still {drifts[2][2]:.3%}")
    # Check approximate O(dr^2) convergence
    if len(drifts) >= 2:
        ratio = drifts[0][2] / drifts[2][2]
        print(f"  Drift ratio (coarse/fine) = {ratio:.1f} "
              f"(expected ~16 for O(dr^2), with grid refinement 4x)")
    return drifts


# =============================================================================
# Test 2: Bessel-mode oscillation + damping
# =============================================================================

def damped_cosine(t, A, gamma, omega, phi0, offset):
    return offset + A * np.exp(-gamma * t) * np.cos(omega * t + phi0)


def bessel_mode_test(
    nr=400, nz=80,
    nu_cm2_s=0.1,         # higher than water => visible damping
    L_cm_target=25.0,     # tuned so first J0 zero matches a sensible k
    bessel_zero_idx=1,    # 1 = first zero (longest wavelength), 2 = second, ...
    amplitude_cm=0.01,
    n_periods=4.0,
    nt_per_period=400,
):
    """
    Initialise eta(r, 0) = A * J_0(k r) with k = alpha_n / L_dim, where alpha_n
    is the n-th zero of J_0.  Run, fit damped cosine to eta(r=0, t), compare
    measured (omega, gamma) with analytical predictions.

    Analytical predictions (linearised theory):
      omega^2 = (k/Fr + k^3/We) * tanh(k D)              (deep-water with depth)
      gamma  ≈ 2 k^2 / Re                                (Lamb 1932, high-Re limit)
    """
    # Tune L to put J0's n-th zero at r = L
    alpha_n = jn_zeros(0, bessel_zero_idx)[bessel_zero_idx - 1]

    # Estimate the dimensionless k and pick run time
    # We non-dimensionalise with L_c = unit_length = 2.5 cm.
    L_c = 2.5
    L_dim = L_cm_target / L_c          # dimensionless domain radius
    k = alpha_n / L_dim                # dimensionless wavenumber

    phys = PhysicalParams(L_cm=L_cm_target, D_cm=20.0, nu=nu_cm2_s)
    # Estimate omega for sizing nt
    # In dimensionless form with default scales:
    V_c = L_c / 0.05
    Fr = V_c**2 / (phys.g * L_c)
    We = phys.rho * V_c**2 * L_c / phys.sigma
    Re = V_c * L_c / phys.nu
    D = phys.D_cm / L_c
    omega2 = (k / Fr + k**3 / We) * np.tanh(k * D)
    omega = float(np.sqrt(max(omega2, 1e-12)))
    period = 2.0 * np.pi / omega       # dimensionless period
    period_s = period * 0.05
    total_time_s = n_periods * period_s
    nt = int(n_periods * nt_per_period)

    print(f"Bessel zero α_{bessel_zero_idx} = {alpha_n:.4f}")
    print(f"  L_dim = {L_dim:.4f}  →  k = α/L = {k:.4f}")
    print(f"  Fr = {Fr:.4f}, We = {We:.4f}, Re = {Re:.0f}")
    print(f"  predicted omega = {omega:.4f},  period = {period_s*1000:.1f} ms")
    print(f"  predicted gamma = 2k²/Re = {2*k**2/Re:.4e}")
    print(f"  total_time_s = {total_time_s:.4f},  nt = {nt}")

    # Build the solver, then OVERRIDE the initial eta with the Bessel mode.
    cfg = SimConfig(
        phys=phys,
        num=NumericalParams(nr=nr, nz=nz, nt=nt, total_time_s=total_time_s),
        ic=InitialCondition(amplitude_cm=amplitude_cm, width_cm=5.0),  # placeholder
        verbose=False,
    )
    solver = SurfaceDynamicsSolverAxi(cfg)

    # Manual run loop, replacing the initial eta with the Bessel mode
    nr_, nz_ = solver.nr, solver.nz
    save_every = max(1, nt // 200)

    amp_dim = amplitude_cm / L_c
    eta = amp_dim * j0(k * solver.r)
    eta[-1] = 0.0   # outer Dirichlet
    phi = np.zeros((nz_, nr_))
    w3 = np.zeros((nz_, nr_))
    w1 = np.zeros((nz_, nr_))

    eta_axis = []
    time_arr = []
    for n in range(nt):
        w3, w1, phi, eta = solver._step(w3, w1, phi, eta)
        if n % save_every == 0 or n == nt - 1:
            eta_axis.append(eta[0])
            time_arr.append((n + 1) * solver.dt)

    eta_axis = np.array(eta_axis)
    time_arr = np.array(time_arr)

    # Convert eta_axis back to cm for nicer numbers
    eta_axis_cm = eta_axis * L_c

    # Fit damped cosine on the axis trace
    # Skip the first ~5% of the trace to avoid initial transient effects
    skip = max(1, len(eta_axis) // 20)
    t_fit = time_arr[skip:]
    y_fit = eta_axis[skip:]

    p0 = [eta_axis[0], 2 * k**2 / Re, omega, 0.0, 0.0]
    try:
        popt, _ = curve_fit(damped_cosine, t_fit, y_fit, p0=p0, maxfev=5000)
        A_fit, gamma_fit, omega_fit, phi0_fit, offset_fit = popt
        gamma_fit = abs(gamma_fit)
        omega_fit = abs(omega_fit)
    except Exception as e:
        print(f"  Fit failed: {e}")
        return None

    # Compare with analytical
    omega_pred = omega
    gamma_pred = 2 * k**2 / Re

    print()
    print("Measured vs. predicted:")
    print(f"  omega:  measured = {omega_fit:.4f},  predicted = {omega_pred:.4f},  "
          f"rel err = {abs(omega_fit - omega_pred)/omega_pred:.2%}")
    print(f"  gamma:  measured = {gamma_fit:.4e},  predicted = {gamma_pred:.4e},  "
          f"rel err = {abs(gamma_fit - gamma_pred)/max(gamma_pred,1e-30):.2%}")

    return {
        "k": k, "L": L_dim, "Fr": Fr, "We": We, "Re": Re,
        "omega_pred": omega_pred, "gamma_pred": gamma_pred,
        "omega_fit": omega_fit, "gamma_fit": gamma_fit,
        "time": time_arr, "eta_axis": eta_axis,
        "fit_params": popt,
    }


def run_bessel_validation():
    print("=" * 60)
    print("Test 2a: Bessel-mode validation (high Re — water)")
    print("=" * 60)
    res_high = bessel_mode_test(nu_cm2_s=0.0005, n_periods=3.0, nt_per_period=300)
    print()
    print("=" * 60)
    print("Test 2b: Bessel-mode validation (moderate Re — light oil)")
    print("=" * 60)
    res_low = bessel_mode_test(nu_cm2_s=0.1, n_periods=3.0, nt_per_period=300)
    return res_high, res_low


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    run_grid_convergence()
    print()
    run_bessel_validation()
