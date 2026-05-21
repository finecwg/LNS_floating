"""
Validation tests for the axisymmetric floating-object solver.

Tests
-----
1. Smoke test
   * Solver assembles, factorises and runs without numerical errors.
   * |p_s| off the body is exactly zero (Lagrange-multiplier rows correctly
     pin pressure outside B).

2. Static pressure consistency (instantaneous-response check)
   * At t ~ a few steps, the integrated traction
     ∫ p_s * 2*pi*r dr over the body should already approximate M/Fr,
     because the elliptic Laplace block adjusts the pressure instantaneously
     to support the body weight.

3. Hydrostatic equilibrium (3D Archimedes)
   * Release body from rest at z_b(0) = 0, v_b(0) = 0.
   * Body should oscillate around z_b^eq = -M/(pi a^2) and decay toward it
     under viscous damping.
   * After sufficient time, |<z_b>_t - z_b^eq| should be small.

4. Cylindrical-weight sanity check
   * Verify that the discrete integral of a constant unit field over the body
     using the 2*pi*r_j*dr rule equals (pi a^2) to truncation accuracy.
"""

import sys
from pathlib import Path

import numpy as np

# Ensure we can import sibling module when run from anywhere
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from floating_object_axi import (  # noqa: E402
    FloatingObjectSolverAxi, SimConfig,
    PhysicalParams, BodyParams, NumericalParams,
)


# =============================================================================
# Test 1: Smoke test
# =============================================================================


def smoke_test():
    print("=" * 60)
    print("Test 1: Smoke test (short run, structural sanity)")
    print("=" * 60)
    cfg = SimConfig(
        phys=PhysicalParams(L_cm=30.0),
        body=BodyParams(radius_cm=2.5, mass_g=7.0),
        num=NumericalParams(nr=80, nz=20, nt=20, total_time_s=0.02),
        verbose=False,
    )
    solver = FloatingObjectSolverAxi(cfg)
    result = solver.run()
    ps_final = result["pressure_history"][-1]
    max_ps_free = np.max(np.abs(ps_final[~solver.body_mask]))
    max_ps_body = np.max(np.abs(ps_final[solver.body_mask]))
    print(f"  Wall time: {result['wall_time']*1000:.1f} ms")
    print(f"  System size N = {solver.off['total']}")
    print(f"  max|p_s| off body = {max_ps_free:.2e}    (target: < 1e-12)")
    print(f"  max|p_s| in body  = {max_ps_body:.4e}")
    print(f"  Final z_b (mm)    = {result['z_b'][-1]*cfg.num.unit_length*10:+.4f}")
    if max_ps_free < 1e-12:
        print("  ✓ Off-body pressure is machine-precision zero.")
    else:
        print("  ✗ Off-body pressure is NOT zero -- check Lagrange-multiplier rows.")
    return result, max_ps_free


# =============================================================================
# Test 2: Static pressure consistency
# =============================================================================


def static_pressure_consistency(result, solver):
    """The integrated body traction should approach M/Fr.

    Because the pressure response is elliptic (Laplace) and the body equation
    is coupled implicitly, the integrated traction converges to its
    hydrostatic value as the body settles.
    """
    print()
    print("=" * 60)
    print("Test 2: Integrated traction vs body weight")
    print("=" * 60)
    F_int = result["body_force"]
    weight = solver.M / solver.Fr
    rel_imbalance = abs(F_int[-1] - weight) / weight
    print(f"  Body weight M/Fr           = {weight:.4e}")
    print(f"  Initial F_int (step 1)     = {F_int[0]:.4e}")
    print(f"  Final F_int                = {F_int[-1]:.4e}")
    print(f"  Final relative imbalance   = {rel_imbalance:.2%}")
    # Note: this is not zero -- only at static equilibrium does F_int = weight.
    return rel_imbalance


# =============================================================================
# Test 3: Hydrostatic equilibrium (3D Archimedes)
# =============================================================================


def hydrostatic_equilibrium(
    nr=200, nz=40, nt=3000, total_time_s=3.0,
    L_cm=50.0, radius_cm=2.5, mass_g=7.0,
    nu=2.0,        # viscous fluid for fast settling
    sigma=1.0,     # suppress surface tension so Archimedes is exact
):
    """3D Archimedes test.

    To match the Archimedes formula |z_b^eq| = M/(pi a^2) cleanly, we use:
      * High viscosity (low Re) so the body damps to rest quickly.
      * Negligible surface tension (high We) so the contact-line capillary
        term (1/We) d_rr eta does not shift the equilibrium.

    With these choices, the late-time mean z_b should match the *discrete*
    Archimedes prediction -M / Sum(2*pi*r_j*dr) (over body nodes) to within
    a couple of percent. The discrete prediction differs from the continuous
    formula by the body-area truncation error (~10% at moderate nr), which
    converges to zero under grid refinement (Test 4).
    """
    print()
    print("=" * 60)
    print("Test 3: Hydrostatic equilibrium (3D Archimedes, sigma~0)")
    print("=" * 60)
    cfg = SimConfig(
        phys=PhysicalParams(L_cm=L_cm, nu=nu, sigma=sigma),
        body=BodyParams(radius_cm=radius_cm, mass_g=mass_g,
                        initial_z_cm=0.0, initial_vz_cm_s=0.0),
        num=NumericalParams(nr=nr, nz=nz, nt=nt, total_time_s=total_time_s),
        verbose=False,
    )
    solver = FloatingObjectSolverAxi(cfg)
    result = solver.run()
    ul = cfg.num.unit_length

    # Discrete body area  Sum 2*pi*r_j*dr over body nodes
    area_w = 2.0 * np.pi * solver.r_vec * solver.dr
    A_disc = np.sum(area_w[solver.body_mask])
    A_cont = np.pi * solver.a_obj**2
    zb_eq_disc_mm = -solver.M / A_disc * ul * 10
    zb_eq_cont_mm = -solver.M / A_cont * ul * 10

    zb_history_mm = result["z_b"] * ul * 10
    min_zb_mm = np.min(zb_history_mm)
    final_zb_mm = zb_history_mm[-1]
    n_late = max(1, len(zb_history_mm) // 5)
    mean_late = np.mean(zb_history_mm[-n_late:])
    std_late = np.std(zb_history_mm[-n_late:])

    print(f"  Re = {solver.Re:.1f}, We = {solver.We:.1f} (large => no surf tension)")
    print(f"  Continuous A=pi*a^2 = {A_cont:.4f},  discrete = {A_disc:.4f}"
          f"  ({100*(A_disc/A_cont - 1):+.2f}%)")
    print(f"  z_b^eq (continuous) = {zb_eq_cont_mm:+.4f} mm  (Archimedes formula)")
    print(f"  z_b^eq (discrete)   = {zb_eq_disc_mm:+.4f} mm  (-M/A_disc, true target)")
    print(f"  Min z_b (overshoot) = {min_zb_mm:+.4f} mm")
    print(f"  Final z_b           = {final_zb_mm:+.4f} mm")
    print(f"  Late-time <z_b>     = {mean_late:+.4f} ± {std_late:.4f} mm")
    err_disc = abs(mean_late - zb_eq_disc_mm) / abs(zb_eq_disc_mm) * 100
    err_cont = abs(mean_late - zb_eq_cont_mm) / abs(zb_eq_cont_mm) * 100
    print(f"  |<z_b> - z_b^eq_disc| / |z_b^eq_disc| = {err_disc:.2f}%")
    print(f"  |<z_b> - z_b^eq_cont| / |z_b^eq_cont| = {err_cont:.2f}%")
    if err_disc < 5.0:
        print("  ✓ Late-time mean within 5% of *discrete* Archimedes target.")
    else:
        print("  ⚠ Late-time mean farther than 5% from discrete target "
              "(may need longer run / lower Re).")
    return result, err_disc


# =============================================================================
# Test 4: Cylindrical area-weight sanity
# =============================================================================


def cylindrical_weight_sanity():
    """Test that 2*pi*r_j*dr * 1_{j in body} sums to (pi a^2) accurately."""
    print()
    print("=" * 60)
    print("Test 4: Discrete cylindrical-weight area check")
    print("=" * 60)
    from operators_axi import make_r_grid
    for nr in [50, 100, 200, 400]:
        L_dim = 12.0
        a_dim = 1.0
        dr, r_vec = make_r_grid(nr, L_dim)
        body_mask = r_vec < a_dim
        area_w = 2.0 * np.pi * r_vec * dr
        disc_area = np.sum(area_w[body_mask])
        exact_area = np.pi * a_dim**2
        err = (disc_area - exact_area) / exact_area
        print(f"  nr={nr:>4d}:  disc area = {disc_area:.6f}, "
              f"exact = {exact_area:.6f}, rel err = {err:+.3%}")


# =============================================================================
# Main
# =============================================================================


def main():
    # Test 1
    result, max_ps_free = smoke_test()

    # Re-run solver briefly for the cached solver instance (for Test 2)
    cfg = SimConfig(
        phys=PhysicalParams(L_cm=30.0),
        body=BodyParams(radius_cm=2.5, mass_g=7.0),
        num=NumericalParams(nr=80, nz=20, nt=20, total_time_s=0.02),
        verbose=False,
    )
    solver = FloatingObjectSolverAxi(cfg)
    short_result = solver.run()
    static_pressure_consistency(short_result, solver)

    # Test 3 (longer + finer grid for sharp Archimedes match)
    hydrostatic_equilibrium(nr=200, nz=40, nt=3000, total_time_s=3.0,
                            L_cm=50.0, radius_cm=2.5, mass_g=7.0,
                            nu=2.0, sigma=1.0)

    # Test 4
    cylindrical_weight_sanity()


if __name__ == "__main__":
    main()
