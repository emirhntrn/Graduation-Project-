# ============================================================
# NSGA-II OPTIMIZATION CODE
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
import time
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ansys.mapdl.core import launch_mapdl

from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.optimize import minimize

# Solver file must be in the same folder
from profile_plate_solver_coupled_wet_straight_profiles import evaluate_design


# ============================================================
# USER SETTINGS
# ============================================================

# Reference design
REF_T_PLATE = 0.010       # mm
REF_T_PROFILE = 0.010     # 10 mm
REF_H_PROFILE = 0.010     # mm profile height
REF_NS = 0.250            # [m]

# Thickness and profile height are expanded to 6 values.
# 10-14 mm range with 6 discrete values:
T_PLATE_VALUES = np.linspace(0.010, 0.014, 4)
T_PROFILE_VALUES = np.linspace(0.010, 0.014, 4)
H_PROFILE_VALUES = np.linspace(0.010, 0.014, 4)
NS_VALUES = np.array([0.200, 0.250, 0.333, 0.500])

# Modal / mesh settings
NMODE = 5
FREQ_TOL = 1.0
NEL_FACTOR = 30

# NSGA-II settings
POP_SIZE = 10
N_OFFSPRINGS = 5
N_EVAL = 50
SEED = 7

CROSSOVER_PROB = 0.90
CROSSOVER_ETA = 20
MUTATION_ETA = 25

JOBNAME = "profile_plate_opt_job"


# ============================================================
# FUNCTIONS
# ============================================================

def get_base_dir():
    """
    Safe folder definition for both .py files and Jupyter notebooks.
    """
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return os.getcwd()


def make_unique_workdir():
    base_dir = get_base_dir()
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    workdir = os.path.join(base_dir, f"ansys_workdir_opt_{timestamp}")
    os.makedirs(workdir, exist_ok=True)
    return workdir


def launch_solver_mapdl(workdir):
    mapdl = launch_mapdl(
        run_location=workdir,
        jobname=JOBNAME,
        override=True,
        clear_on_connect=True,
        loglevel="ERROR",
        additional_switches="-smp",
        license_type="ansys",
    )
    return mapdl


def safe_ratio(numerator, denominator, large_penalty=1.0e6):
    """
    Safe division for objective and normalized ratios.
    """
    try:
        if denominator is None:
            return large_penalty
        if not np.isfinite(denominator):
            return large_penalty
        if abs(float(denominator)) < 1.0e-12:
            return large_penalty
        return float(numerator) / float(denominator)
    except Exception:
        return large_penalty


# ============================================================
# OPTIMIZATION PROBLEM
# ============================================================

class ProfilePlateOptimizationProblem(ElementwiseProblem):
    def __init__(self, mapdl, workdir, ref_freq, ref_mass):
        self.mapdl = mapdl
        self.workdir = workdir
        self.ref_freq = float(ref_freq)
        self.ref_mass = float(ref_mass)

        self.eval_counter = 0
        self.history = []

        # Variables are integer indices:
        # x[0] -> T_PLATE_VALUES index
        # x[1] -> T_PROFILE_VALUES index
        # x[2] -> NS_VALUES index
        # x[3] -> H_PROFILE_VALUES index
        xl = np.array([0, 0, 0, 0])
        xu = np.array([
            len(T_PLATE_VALUES) - 1,
            len(T_PROFILE_VALUES) - 1,
            len(NS_VALUES) - 1,
            len(H_PROFILE_VALUES) - 1,
        ])

        super().__init__(
            n_var=4,
            n_obj=2,
            n_constr=0,
            xl=xl,
            xu=xu,
            vtype=int,
        )

    def _evaluate(self, x, out, *args, **kwargs):
        self.eval_counter += 1

        i_t_plate = int(round(x[0]))
        i_t_profile = int(round(x[1]))
        i_NS = int(round(x[2]))
        i_h_profile = int(round(x[3]))

        i_t_plate = int(np.clip(i_t_plate, 0, len(T_PLATE_VALUES) - 1))
        i_t_profile = int(np.clip(i_t_profile, 0, len(T_PROFILE_VALUES) - 1))
        i_NS = int(np.clip(i_NS, 0, len(NS_VALUES) - 1))
        i_h_profile = int(np.clip(i_h_profile, 0, len(H_PROFILE_VALUES) - 1))

        t_plate = float(T_PLATE_VALUES[i_t_plate])
        t_profile = float(T_PROFILE_VALUES[i_t_profile])
        NS = float(NS_VALUES[i_NS])
        h_profile = float(H_PROFILE_VALUES[i_h_profile])

        print("\n" + "-" * 80)
        print(f"[EVAL {self.eval_counter}] START")
        print(
            f"[EVAL {self.eval_counter}] "
            f"t_plate={t_plate:.4f} m | "
            f"t_profile={t_profile:.4f} m | "
            f"NS={NS:.4f} m | "
            f"h_profile={h_profile:.4f} m"
        )

        try:
            t0 = time.time()

            freq, volume, mass, info = evaluate_design(
                t_plate=t_plate,
                t_profile=t_profile,
                NS=NS,
                HS=h_profile,
                nmode=NMODE,
                freq_tol=FREQ_TOL,
                nel_factor=NEL_FACTOR,
                visualize=False,
                mapdl=self.mapdl,
                workdir=self.workdir,
                jobname=JOBNAME,
                verbose=False,
            )

            elapsed = time.time() - t0

            if freq is None or mass is None:
                raise RuntimeError("freq or mass returned None.")

            if not np.isfinite(freq) or not np.isfinite(mass):
                raise RuntimeError("freq or mass is not finite.")

            if freq <= FREQ_TOL:
                raise RuntimeError(
                    f"Invalid selected frequency: {freq}. "
                    "Increase NMODE or check model constraints/FSI settings."
                )

            # Objective values
            # REF_FREQ / freq is minimized.
            f1 = safe_ratio(self.ref_freq, freq)

            # mass / REF_MASS is minimized.
            f2 = safe_ratio(mass, self.ref_mass)

            freq_ratio = safe_ratio(freq, self.ref_freq)
            mass_ratio = safe_ratio(mass, self.ref_mass)

            selected_mode = info.get("selected_mode", None)
            n_profiles = info.get("n_profiles", None)
            profile_z = info.get("profile_z_positions", None)

            print(
                f"[RESULT {self.eval_counter}] "
                f"freq={freq:.6f} Hz | "
                f"freq/ref={freq_ratio:.6f} | "
                f"mass={mass:.6f} kg | "
                f"mass/ref={mass_ratio:.6f}"
            )
            print(
                f"[RESULT {self.eval_counter}] "
                f"volume={volume:.9f} m^3 | "
                f"mode={selected_mode} | "
                f"profiles={n_profiles} | "
                f"time={elapsed:.2f} s"
            )
            print(
                f"[OBJECTIVE {self.eval_counter}] "
                f"f1=RefFreq/Freq={f1:.6f} | "
                f"f2=Mass/RefMass={f2:.6f}"
            )

            out["F"] = np.array([f1, f2], dtype=float)

            self.history.append({
                "eval": self.eval_counter,
                "status": "OK",
                "t_plate_m": t_plate,
                "t_profile_m": t_profile,
                "h_profile_m": h_profile,
                "NS_m": NS,
                "frequency_Hz": float(freq),
                "frequency_over_ref": float(freq_ratio),
                "volume_m3": float(volume),
                "mass_kg": float(mass),
                "mass_over_ref": float(mass_ratio),
                "objective_f1_ref_over_freq": float(f1),
                "objective_f2_mass_over_ref": float(f2),
                "selected_mode": selected_mode,
                "n_profiles": n_profiles,
                "profile_z_positions": str(profile_z),
                "time_s": float(elapsed),
                "error": "",
            })

        except Exception as e:
            print(f"[FAILED {self.eval_counter}] {e}")
            traceback.print_exc()

            penalty = 1.0e6
            out["F"] = np.array([penalty, penalty], dtype=float)

            self.history.append({
                "eval": self.eval_counter,
                "status": "FAILED",
                "t_plate_m": t_plate,
                "t_profile_m": t_profile,
                "h_profile_m": h_profile,
                "NS_m": NS,
                "frequency_Hz": np.nan,
                "frequency_over_ref": np.nan,
                "volume_m3": np.nan,
                "mass_kg": np.nan,
                "mass_over_ref": np.nan,
                "objective_f1_ref_over_freq": penalty,
                "objective_f2_mass_over_ref": penalty,
                "selected_mode": None,
                "n_profiles": None,
                "profile_z_positions": "",
                "time_s": np.nan,
                "error": str(e),
            })


# ============================================================
# POST-PROCESS HELPERS
# ============================================================

def build_pareto_dataframe(res, ref_freq, ref_mass):
    pareto_rows = []

    if res.X is None or res.F is None:
        return pd.DataFrame(pareto_rows)

    X = np.atleast_2d(res.X)
    F = np.atleast_2d(res.F)

    for i, (x, f) in enumerate(zip(X, F), start=1):
        i_t_plate = int(np.clip(round(x[0]), 0, len(T_PLATE_VALUES) - 1))
        i_t_profile = int(np.clip(round(x[1]), 0, len(T_PROFILE_VALUES) - 1))
        i_NS = int(np.clip(round(x[2]), 0, len(NS_VALUES) - 1))
        i_h_profile = int(np.clip(round(x[3]), 0, len(H_PROFILE_VALUES) - 1))

        t_plate = float(T_PLATE_VALUES[i_t_plate])
        t_profile = float(T_PROFILE_VALUES[i_t_profile])
        NS = float(NS_VALUES[i_NS])
        h_profile = float(H_PROFILE_VALUES[i_h_profile])

        f1 = float(f[0])
        f2 = float(f[1])

        # Back conversion:
        # f1 = REF_FREQ / freq  -> freq = REF_FREQ / f1
        # f2 = mass / REF_MASS  -> mass = f2 * REF_MASS
        freq = safe_ratio(ref_freq, f1)
        mass = float(f2 * ref_mass)

        pareto_rows.append({
            "pareto_id": i,
            "t_plate_m": t_plate,
            "t_profile_m": t_profile,
            "h_profile_m": h_profile,
            "NS_m": NS,
            "frequency_Hz": freq,
            "mass_kg": mass,
            "frequency_over_ref": safe_ratio(freq, ref_freq),
            "mass_over_ref": safe_ratio(mass, ref_mass),
            "objective_f1_ref_over_freq": f1,
            "objective_f2_mass_over_ref": f2,
        })

    df_pareto = pd.DataFrame(pareto_rows)

    if len(df_pareto) > 0:
        df_pareto = df_pareto.drop_duplicates(
            subset=["t_plate_m", "t_profile_m", "h_profile_m", "NS_m", "frequency_Hz", "mass_kg"]
        ).reset_index(drop=True)
        df_pareto["pareto_id"] = np.arange(1, len(df_pareto) + 1)

    return df_pareto


def add_selected_optimum(df_pareto):
    selected_optimum = None

    if len(df_pareto) == 0:
        return df_pareto, selected_optimum

    df_pareto = df_pareto.copy()

    freq_min = df_pareto["frequency_Hz"].min()
    freq_max = df_pareto["frequency_Hz"].max()

    mass_min = df_pareto["mass_kg"].min()
    mass_max = df_pareto["mass_kg"].max()

    # Frequency: higher is better, so ideal score is 0 at maximum frequency.
    if abs(freq_max - freq_min) < 1e-12:
        df_pareto["freq_score"] = 0.0
    else:
        df_pareto["freq_score"] = 1.0 - (
            (df_pareto["frequency_Hz"] - freq_min) / (freq_max - freq_min)
        )

    # Mass: lower is better, so ideal score is 0 at minimum mass.
    if abs(mass_max - mass_min) < 1e-12:
        df_pareto["mass_score"] = 0.0
    else:
        df_pareto["mass_score"] = (
            (df_pareto["mass_kg"] - mass_min) / (mass_max - mass_min)
        )

    df_pareto["distance_to_ideal"] = np.sqrt(
        df_pareto["freq_score"] ** 2 + df_pareto["mass_score"] ** 2
    )

    selected_idx = df_pareto["distance_to_ideal"].idxmin()
    selected_optimum = df_pareto.loc[selected_idx].copy()

    df_pareto["selected_optimum"] = False
    df_pareto.loc[selected_idx, "selected_optimum"] = True

    return df_pareto, selected_optimum


def save_and_plot_results(workdir, df_all, df_pareto, selected_optimum):
    all_eval_path = os.path.join(workdir, "all_evaluations.csv")
    df_all.to_csv(all_eval_path, index=False)
    print(f"[SAVED] All evaluations: {all_eval_path}")
    
    pareto_path = os.path.join(workdir, "pareto_results.csv")
    df_pareto.to_csv(pareto_path, index=False)
    print(f"[SAVED] Pareto results: {pareto_path}")

    if selected_optimum is not None:
        selected_path = os.path.join(workdir, "selected_optimum.csv")
        pd.DataFrame([selected_optimum]).to_csv(selected_path, index=False)
        print(f"[SAVED] Selected optimum: {selected_path}")

    # ------------------------------------------------------------
    # Print all Pareto values
    # ------------------------------------------------------------
    print("\n" + "=" * 80)
    print("--- ALL PARETO FRONT VALUES ---")
    print("=" * 80)

    if len(df_pareto) > 0:
        pareto_print_cols = [
            "pareto_id",
            "t_plate_m",
            "t_profile_m",
            "h_profile_m",
            "NS_m",
            "frequency_Hz",
            "mass_kg",
            "frequency_over_ref",
            "mass_over_ref",
            "objective_f1_ref_over_freq",
            "objective_f2_mass_over_ref",
            "distance_to_ideal",
            "selected_optimum",
        ]
        available_cols = [c for c in pareto_print_cols if c in df_pareto.columns]
        print(
            df_pareto[available_cols]
            .sort_values(by="mass_kg")
            .to_string(index=False)
        )
    else:
        print("No Pareto result returned.")

    # ------------------------------------------------------------
    # Print selected optimum
    # ------------------------------------------------------------
    print("\n" + "=" * 80)
    print("--- SELECTED OPTIMUM DESIGN ---")
    print("=" * 80)

    if selected_optimum is not None:
        print(f"Selected Pareto ID      = {int(selected_optimum['pareto_id'])}")
        print(f"Plate thickness         = {selected_optimum['t_plate_m']:.4f} m")
        print(f"Profile thickness       = {selected_optimum['t_profile_m']:.4f} m")
        print(f"Profile height          = {selected_optimum['h_profile_m']:.4f} m")
        print(f"NS spacing              = {selected_optimum['NS_m']:.4f} m")
        print(f"Frequency               = {selected_optimum['frequency_Hz']:.6f} Hz")
        print(f"Mass                    = {selected_optimum['mass_kg']:.6f} kg")
        print(f"Frequency / Ref         = {selected_optimum['frequency_over_ref']:.6f}")
        print(f"Mass / Ref              = {selected_optimum['mass_over_ref']:.6f}")
        print(f"Objective f1=Ref/Freq   = {selected_optimum['objective_f1_ref_over_freq']:.6f}")
        print(f"Objective f2=Mass/Ref   = {selected_optimum['objective_f2_mass_over_ref']:.6f}")
        print(f"Distance to ideal       = {selected_optimum['distance_to_ideal']:.6f}")
    else:
        print("No selected optimum available.")

    # ------------------------------------------------------------
    # Prepare successful evaluations
    # ------------------------------------------------------------
    df_ok = df_all[df_all["status"] == "OK"].copy()

    if len(df_ok) > 0:
        df_ok = df_ok.sort_values(by="eval").reset_index(drop=True)

        # ------------------------------------------------------------
        # Plot frequency convergence with evaluation number
        # ------------------------------------------------------------
        df_ok["best_so_far_frequency_Hz"] = df_ok["frequency_Hz"].cummax()

        plt.figure(figsize=(9, 6))
        plt.plot(
            df_ok["eval"],
            df_ok["frequency_Hz"],
            marker="o",
            label="Frequency at each evaluation",
        )
        plt.plot(
            df_ok["eval"],
            df_ok["best_so_far_frequency_Hz"],
            marker="s",
            linestyle="--",
            label="Best-so-far frequency",
        )
        plt.xlabel("Evaluation number")
        plt.ylabel("Selected wet frequency [Hz]")
        plt.title("Frequency Convergence During Optimization")
        plt.grid(True)
        plt.legend()

        convergence_plot_path = os.path.join(workdir, "frequency_convergence_by_eval.png")
        plt.savefig(convergence_plot_path, dpi=300, bbox_inches="tight")
        plt.show()
        print(f"[SAVED] Frequency convergence plot: {convergence_plot_path}")

        # ------------------------------------------------------------
        # Plot mass convergence with evaluation number
        # ------------------------------------------------------------
        df_ok["best_so_far_mass_kg"] = df_ok["mass_kg"].cummin()

        plt.figure(figsize=(9, 6))
        plt.plot(
            df_ok["eval"],
            df_ok["mass_kg"],
            marker="o",
            label="Mass at each evaluation",
        )
        plt.plot(
            df_ok["eval"],
            df_ok["best_so_far_mass_kg"],
            marker="s",
            linestyle="--",
            label="Best-so-far minimum mass",
        )
        plt.xlabel("Evaluation number")
        plt.ylabel("Mass [kg]")
        plt.title("Mass Convergence During Optimization")
        plt.grid(True)
        plt.legend()

        mass_convergence_plot_path = os.path.join(workdir, "mass_convergence_by_eval.png")
        plt.savefig(mass_convergence_plot_path, dpi=300, bbox_inches="tight")
        plt.show()
        print(f"[SAVED] Mass convergence plot: {mass_convergence_plot_path}")

        # ------------------------------------------------------------
        # Plot all evaluated designs: Frequency vs Mass
        # ------------------------------------------------------------
        plt.figure(figsize=(8, 6))
        plt.scatter(
            df_ok["mass_kg"],
            df_ok["frequency_Hz"],
            label="All evaluated designs",
        )
        plt.xlabel("Mass [kg]")
        plt.ylabel("Selected wet frequency [Hz]")
        plt.title("All Evaluated Designs - Frequency vs Mass")
        plt.grid(True)
        plt.legend()

        all_plot_path = os.path.join(workdir, "all_designs_frequency_mass.png")
        plt.savefig(all_plot_path, dpi=300, bbox_inches="tight")
        plt.show()
        print(f"[SAVED] All designs frequency-mass plot: {all_plot_path}")

    # ------------------------------------------------------------
    # Plot Pareto front with selected optimum
    # ------------------------------------------------------------
    if len(df_pareto) > 0:
        plt.figure(figsize=(8, 6))
        plt.scatter(
            df_pareto["mass_kg"],
            df_pareto["frequency_Hz"],
            label="Pareto front",
        )

        if selected_optimum is not None:
            plt.scatter(
                selected_optimum["mass_kg"],
                selected_optimum["frequency_Hz"],
                marker="*",
                s=250,
                label="Selected optimum",
            )

        plt.xlabel("Mass [kg]")
        plt.ylabel("Selected wet frequency [Hz]")
        plt.title("Pareto Front - Frequency vs Mass")
        plt.grid(True)
        plt.legend()

        plot_path = os.path.join(workdir, "pareto_frequency_mass.png")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.show()
        print(f"[SAVED] Pareto plot: {plot_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 80)
    print("NSGA-II OPTIMIZATION - COUPLED WET MODAL PROFILE PLATE - PROFILE HEIGHT INCLUDED - STRAIGHT PROFILES + PROFILE HEIGHT")
    print("=" * 80)

    workdir = make_unique_workdir()
    print(f"[INFO] Workdir: {workdir}")

    mapdl = None

    try:
        # ------------------------------------------------------------
        # MAPDL launch
        # ------------------------------------------------------------
        print("[INFO] Launching MAPDL...")
        mapdl = launch_solver_mapdl(workdir=workdir)
        print("[INFO] MAPDL launched.")

        # ------------------------------------------------------------
        # Reference analysis
        # ------------------------------------------------------------
        print("\n" + "=" * 80)
        print("--- REFERENCE ANALYSIS ---")
        print("=" * 80)
        print(
            f"[REF] t_plate={REF_T_PLATE:.4f} m | "
            f"t_profile={REF_T_PROFILE:.4f} m | "
            f"NS={REF_NS:.4f} m | "
            f"h_profile={REF_H_PROFILE:.4f} m"
        )

        t0_ref = time.time()

        ref_freq, ref_volume, ref_mass, ref_info = evaluate_design(
            t_plate=REF_T_PLATE,
            t_profile=REF_T_PROFILE,
            NS=REF_NS,
            HS=REF_H_PROFILE,
            nmode=NMODE,
            freq_tol=FREQ_TOL,
            nel_factor=NEL_FACTOR,
            visualize=False,
            mapdl=mapdl,
            workdir=workdir,
            jobname=JOBNAME,
            verbose=True,
        )

        ref_time = time.time() - t0_ref

        if ref_freq is None or ref_freq <= FREQ_TOL:
            raise RuntimeError(
                "Reference frequency is invalid. "
                "Increase NMODE or check model settings."
            )

        if ref_mass is None or ref_mass <= 0:
            raise RuntimeError("Reference mass is invalid.")

        print("\n--- REFERENCE RESULT ---")
        print(f"Reference wet frequency = {ref_freq:.6f} Hz")
        print(f"Reference volume        = {ref_volume:.9f} m^3")
        print(f"Reference mass          = {ref_mass:.6f} kg")
        print(f"Reference selected mode = {ref_info.get('selected_mode', None)}")
        print(f"Reference profiles      = {ref_info.get('n_profiles', None)}")
        print(f"Reference time          = {ref_time:.2f} s")

        # ------------------------------------------------------------
        # Problem
        # ------------------------------------------------------------
        problem = ProfilePlateOptimizationProblem(
            mapdl=mapdl,
            workdir=workdir,
            ref_freq=ref_freq,
            ref_mass=ref_mass,
        )

        # ------------------------------------------------------------
        # Algorithm
        # ------------------------------------------------------------
        algorithm = NSGA2(
            pop_size=POP_SIZE,
            n_offsprings=N_OFFSPRINGS,
            sampling=IntegerRandomSampling(),
            crossover=SBX(
                prob=CROSSOVER_PROB,
                eta=CROSSOVER_ETA,
                vtype=float,
                repair=None,
            ),
            mutation=PM(
                eta=MUTATION_ETA,
                vtype=float,
                repair=None,
            ),
            eliminate_duplicates=False,
        )

        print("\n" + "=" * 80)
        print("--- OPTIMIZATION STARTED ---")
        print("=" * 80)
        print(f"POP_SIZE     = {POP_SIZE}")
        print(f"N_OFFSPRINGS = {N_OFFSPRINGS}")
        print(f"N_EVAL       = {N_EVAL}")
        print(f"NMODE        = {NMODE}")
        print(f"NEL_FACTOR   = {NEL_FACTOR}")
        print(f"T_PLATE_VALUES    = {T_PLATE_VALUES}")
        print(f"T_PROFILE_VALUES  = {T_PROFILE_VALUES}")
        print(f"H_PROFILE_VALUES  = {H_PROFILE_VALUES}")
        print(f"NS_VALUES         = {NS_VALUES}")
        print("=" * 80)


        res = minimize(
            problem,
            algorithm,
            termination=("n_eval", N_EVAL),
            seed=SEED,
            verbose=True,
            save_history=False,
            copy_algorithm=False,
        )

        print("\n" + "=" * 80)
        print("--- OPTIMIZATION FINISHED ---")
        print("=" * 80)

        # ------------------------------------------------------------
        # Collect results
        # ------------------------------------------------------------
        df_all = pd.DataFrame(problem.history)

        df_pareto = build_pareto_dataframe(
            res=res,
            ref_freq=ref_freq,
            ref_mass=ref_mass,
        )

        df_pareto, selected_optimum = add_selected_optimum(df_pareto)

        # ------------------------------------------------------------
        # Save and plot results
        # ------------------------------------------------------------
        save_and_plot_results(
            workdir=workdir,
            df_all=df_all,
            df_pareto=df_pareto,
            selected_optimum=selected_optimum,
        )

        print("\n" + "=" * 80)
        print("DONE")
        print("=" * 80)
        print("\nOutput folder:")
        print(workdir)

        return {
            "workdir": workdir,
            "all_evaluations": df_all,
            "pareto_results": df_pareto,
            "selected_optimum": selected_optimum,
        }

    finally:
        if mapdl is not None:
            try:
                mapdl.exit()
                print("[INFO] MAPDL closed.")
            except Exception:
                pass


if __name__ == "__main__":
    main()
