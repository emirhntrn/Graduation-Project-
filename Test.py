# ============================================================
# NSGA-II OPTIMIZATION TEST CODE
# RECTANGULAR PLATE-PROFILE PARTIALLY IN CONTACT WITH FLUID
# HYDROELASTIC VIBRATION ANALYSIS
# Objective 1: maximize selected wet frequency
# Objective 2: minimize structural mass
# Added design variable: straight profile height HS
#
# pymoo minimization formulation:
#   f1 = REF_FREQ / freq
#   f2 = mass / REF_MASS
# Usage: For academic and research purposes only.
# ============================================================

import os
import glob
import traceback

from profile_plate_solver_coupled_wet_straight_profiles import evaluate_design, launch_solver_mapdl


def clean_old_files():
    workdir = os.path.join(os.getcwd(), "ansys_workdir")
    os.makedirs(workdir, exist_ok=True)

    for pattern in ["*.err", "*.out", "*.log"]:
        for f in glob.glob(os.path.join(workdir, pattern)):
            try:
                os.remove(f)
            except Exception:
                pass

    return workdir


def run_system_test():
    workdir = clean_old_files()
    mapdl = None

    print("=== COUPLED WET MODAL TEST - STRAIGHT DRY-SIDE PROFILES ONLY ===\n")

    # Reference design.
    t_plate = 0.010       # plate thickness [m]
    t_profile = 0.010     # profile thickness [m]
    NS = 0.250            # profile spacing [m] -> straight profiles at z=[0.25,0.50,0.75]
    nmode = 5             # final model has physical wet modes from mode 1
    freq_tol = 1.0        # removes acoustic/rigid near-zero wet modes

    try:
        mapdl = launch_solver_mapdl(workdir=workdir)

        freq, volume, mass, info = evaluate_design(
            t_plate=t_plate,
            t_profile=t_profile,
            NS=NS,
            nmode=nmode,
            freq_tol=freq_tol,
            nel_factor=40,
            visualize=True,
            mapdl=mapdl,
            workdir=workdir,
            verbose=True,
        )

        print("\n=== TEST FINISHED SUCCESSFULLY ===")
        print(f"Selected wet frequency       = {freq:.6f} Hz")
        print(f"Selected mode number         = {info['selected_mode']}")
        print(f"TOTVOL / structural volume   = {volume:.9f} m^3")
        print(f"Mass                         = {mass:.6f} kg")
        print(f"Number of dry-side profiles  = {info['n_profiles']}")
        print(f"Profile z positions          = {info['profile_z']}")
        print(f"Near-zero mode count (<1 Hz) = {sum(1 for f in info['mode_frequencies'] if abs(f) < 1.0)}")
        print(f"First physical modes         = {[f for f in info['mode_frequencies'] if f > 1.0][:5]}")

        print("\n--- Reference ratios for identical reference design ---")
        print("Frequency / Ref frequency    = 1.000000")
        print("Mass / Ref mass              = 1.000000")
        print("Optimization f1=Ref/Freq     = 1.000000")
        print("Optimization f2=Mass/RefMass = 1.000000")

    except Exception as e:
        print("\n!!! ERROR OCCURRED !!!")
        print("Error detail:", e)
        traceback.print_exc()

    finally:
        if mapdl is not None:
            try:
                mapdl.exit()
            except Exception:
                pass


if __name__ == "__main__":
    run_system_test()
