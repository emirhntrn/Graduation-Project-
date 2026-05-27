# ============================================================
# NSGA-II SOLVER CODE
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

from __future__ import annotations

import math
import os
import time
import uuid
from typing import Dict, List, Optional, Tuple

from ansys.mapdl.core import launch_mapdl


def _safe_base_dir() -> str:
    """Works both from .py files and from notebooks."""
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return os.getcwd()


def _make_workdir(workdir: Optional[str] = None) -> str:
    if workdir is None:
        workdir = os.path.join(_safe_base_dir(), "ansys_workdir")
    os.makedirs(workdir, exist_ok=True)
    return workdir


def launch_solver_mapdl(
    workdir: Optional[str] = None,
    jobname: Optional[str] = None,
    loglevel: str = "ERROR",
):
    """
    Launch a MAPDL instance with a unique jobname to avoid lock-file conflicts.
    Use this in test/optimization and pass the returned mapdl to evaluate_design(...).
    """
    workdir = _make_workdir(workdir)
    if jobname is None:
        jobname = f"profile_plate_{uuid.uuid4().hex[:8]}"

    mapdl = launch_mapdl(
        run_location=workdir,
        jobname=jobname,
        loglevel=loglevel,
        override=True,
        license_type="ansys",
        additional_switches="-smp",
    )
    mapdl.cwd(workdir)
    return mapdl


def _compute_profile_positions(le: float, NS: float) -> List[float]:
    """
    Places profiles at z = NS, 2NS, ... inside (0, le).
    Example: le=1.0, NS=0.25 -> [0.25, 0.5, 0.75]
    """
    if NS <= 0:
        raise ValueError("NS must be positive.")

    positions = []
    k = 1
    while k * NS < le - 1e-9:
        positions.append(round(k * NS, 10))
        k += 1

    return positions


def _read_modal_frequencies(mapdl, nmode: int) -> List[float]:
    """Reads modal frequencies from the result file."""
    freqs = []
    mapdl.post1()

    for mode_id in range(1, nmode + 1):
        try:
            # SET may fail if fewer modes were expanded.
            mapdl.set(1, mode_id)
            mapdl.get("FREQ_NOW", "MODE", mode_id, "FREQ")
            f = float(mapdl.parameters["FREQ_NOW"])
            freqs.append(f)
        except Exception:
            break

    return freqs


def _select_physical_frequency(freqs: List[float], freq_tol: float = 1.0) -> Tuple[float, Optional[int]]:
    """
    Coupled wet modal models may have near-zero acoustic/rigid-like modes.
    Select the first positive physical frequency above freq_tol.
    """
    for i, f in enumerate(freqs, start=1):
        if math.isfinite(f) and f > freq_tol:
            return float(f), i
    return 0.0, None


def _static_mapdl_png(mapdl, workdir: str, filename_prefix: str = "mapdl_mesh_plot") -> Optional[str]:
    """
    Creates an ANSYS/MAPDL-style PNG without opening the GUI.
    """
    try:
        outbase = os.path.join(workdir, filename_prefix).replace("\\", "/")
        mapdl.allsel("ALL")
        mapdl.esel("S", "TYPE", "", 2)
        mapdl.run("/GRAPHICS,POWER")
        mapdl.run("/ESHAPE,1")
        mapdl.run("/EDGE,1,1")
        mapdl.run("/VIEW,1,1,1,1")
        mapdl.run(f"/SHOW,{outbase},PNG")
        mapdl.eplot()
        mapdl.run("/SHOW,CLOSE")
        mapdl.allsel("ALL")

        candidates = [
            os.path.join(workdir, f"{filename_prefix}.PNG"),
            os.path.join(workdir, f"{filename_prefix}000.png"),
            os.path.join(workdir, f"{filename_prefix}000.PNG"),
            os.path.join(workdir, f"{filename_prefix}.png"),
        ]
        for c in candidates:
            if os.path.exists(c):
                return c
    except Exception:
        return None
    return None


def _show_png_if_notebook(path: Optional[str]):
    if not path:
        return
    try:
        from IPython.display import Image, display
        display(Image(filename=path))
    except Exception:
        pass


def evaluate_design(
    t_plate: float = 0.010,
    t_profile: float = 0.010,
    NS: float = 0.250,
    *,
    A: float = 0.100,
    N_wave: int = 2,
    HS: float = 0.040,
    le: float = 1.0,
    he: float = 1.0,
    nmode: int = 5,
    nel_factor: float = 40.0,
    freq_tol: float = 1.0,
    visualize: bool = False,
    mapdl=None,
    workdir: Optional[str] = None,
    jobname: Optional[str] = None,
    verbose: bool = True,
) -> Tuple[float, float, float, Dict]:

    if t_plate <= 0:
        raise ValueError("t_plate must be positive.")
    if t_profile <= 0:
        raise ValueError("t_profile must be positive.")
    if NS <= 0:
        raise ValueError("NS must be positive.")
    if nmode < 1:
        raise ValueError("nmode must be >= 1.")

    start_time = time.time()

    workdir = _make_workdir(workdir)
    local_mapdl = mapdl is None

    if local_mapdl:
        mapdl = launch_solver_mapdl(workdir=workdir, jobname=jobname)

    # ------------------------------------------------------------
    # PARAMETERS
    # ------------------------------------------------------------
    de = he / 2.0
    ti = float(t_plate)
    TS = float(t_profile)
    PI = math.pi

    RHO_S = 7850.0
    E = 206.8e9
    NU = 0.3
    RHO_F = 1025.0

    nel = he / float(nel_factor)
    profile_z = _compute_profile_positions(le, NS)
    n_profiles = len(profile_z)

    try:
        mapdl.clear()
        mapdl.prep7()

        # ------------------------------------------------------------
        # ELEMENT TYPES AND MATERIALS
        # ------------------------------------------------------------
        mapdl.et(1, "FLUID30")
        mapdl.keyopt(1, 2, 0)
        mapdl.keyopt(1, 4, 0)

        mapdl.run("MPTEMP,,,,,,,,")
        mapdl.mptemp(1, 0)
        mapdl.mpdata("DENS", 2, "", RHO_F)
        mapdl.mpdata("SONC", 2, "", 1507)

        mapdl.et(2, "SHELL181")
        mapdl.run("MPTEMP,,,,,,,,")
        mapdl.mptemp(1, 0)
        mapdl.mpdata("DENS", 1, "", RHO_S)
        mapdl.mpdata("EX", 1, "", E)
        mapdl.mpdata("PRXY", 1, "", NU)

        # ------------------------------------------------------------
        # SECTIONS
        # ------------------------------------------------------------
        mapdl.sectype(1, "SHELL")
        mapdl.secdata(ti, 1, 0.0, 3)
        mapdl.secoffset("MID")
        mapdl.run("SECCONTROL,,,,,,,")

        mapdl.sectype(2, "SHELL")
        mapdl.secdata(TS, 1, 0.0, 3)
        mapdl.secoffset("MID")
        mapdl.run("SECCONTROL,,,,,,,")

        # ------------------------------------------------------------
        # ORIGINAL GEOMETRY - KEYPOINTS
        # ------------------------------------------------------------
        mapdl.k(1, 0, -he, (le + he / 2))
        mapdl.k(2, he, -he, (le + he / 2))
        mapdl.k(3, he, -he, (-he / 2))
        mapdl.k(4, 0, -he, (-he / 2))
        mapdl.k(5, 0, 0, (le + he / 2))
        mapdl.k(6, he, 0, (le + he / 2))
        mapdl.k(7, he, 0, (-he / 2))
        mapdl.k(8, 0, 0, (-he / 2))
        mapdl.k(9, 0, (-he / 2), le)
        mapdl.k(10, 0, (-he / 2), 0)
        mapdl.k(11, 0, (he / 2), 0)
        mapdl.k(12, 0, (he / 2), le)
        mapdl.k(13, 0, 0, le)
        mapdl.k(14, 0, 0, 0)

        # ------------------------------------------------------------
        # ORIGINAL GEOMETRY - LINES
        # ------------------------------------------------------------
        mapdl.l(1, 2)
        mapdl.l(2, 3)
        mapdl.l(3, 4)
        mapdl.l(4, 1)
        mapdl.l(5, 6)
        mapdl.l(6, 7)
        mapdl.l(7, 8)
        mapdl.l(8, 14)
        mapdl.l(14, 13)
        mapdl.l(13, 5)
        mapdl.l(5, 1)
        mapdl.l(6, 2)
        mapdl.l(7, 3)
        mapdl.l(8, 4)
        mapdl.l(9, 10)
        mapdl.l(10, 14)
        mapdl.l(14, 11)
        mapdl.l(11, 12)
        mapdl.l(12, 13)
        mapdl.l(13, 9)

        # Original plate area
        mapdl.al(20, 15, 16, 17, 18, 19)

        added_profile_area_ids = []
        for zp in profile_z:
            mapdl.get("KPMAX", "KP", 0, "NUM", "MAX")
            kp = int(mapdl.parameters["KPMAX"]) + 1

            k1 = kp
            k2 = kp + 1
            k3 = kp + 2
            k4 = kp + 3

            mapdl.k(k1, 0.0, -he / 2.0, zp)
            mapdl.k(k2, 0.0,  he / 2.0, zp)
            mapdl.k(k3, -HS,  he / 2.0, zp)
            mapdl.k(k4, -HS, -he / 2.0, zp)

            aid = mapdl.a(k1, k2, k3, k4)
            try:
                added_profile_area_ids.append(int(aid))
            except Exception:
                pass

        # ------------------------------------------------------------
        # SECTION ATTRIBUTES FOR SHELL AREAS
        # ------------------------------------------------------------
        mapdl.allsel("ALL")

        # Plate surfaces at x=0
        mapdl.asel("S", "LOC", "X", 0, 0)
        mapdl.aatt(1, "", 2, 0, 1)

        # Dry-side/profile surfaces at negative x
        mapdl.asel("S", "LOC", "X", -HS - 1e-6, -1e-6)
        mapdl.aatt(1, "", 2, 0, 2)

        # Explicitly assign added profile areas as section 2.
        for aid in added_profile_area_ids:
            mapdl.asel("S", "AREA", "", aid)
            mapdl.aatt(1, "", 2, 0, 2)

        mapdl.allsel("ALL")

        # ------------------------------------------------------------
        # FLUID VOLUME
        # ------------------------------------------------------------
        mapdl.v(1, 2, 3, 4, 5, 6, 7, 8)
        mapdl.vatt(2, "", 1, 0)
        
        try:
            mapdl.asel("S", "LOC", "X", 0, 0)
            mapdl.asel("U", "LOC", "Y", 0, he)
            mapdl.asbl("ALL", "ALL")
        except Exception:
            pass
            
        mapdl.allsel("ALL")
        mapdl.asel("S", "LOC", "X", 0, 0)
        mapdl.aatt(1, "", 2, 0, 1)

        mapdl.asel("S", "LOC", "X", -HS - 1e-6, -1e-6)
        mapdl.aatt(1, "", 2, 0, 2)

        for aid in added_profile_area_ids:
            mapdl.asel("S", "AREA", "", aid)
            mapdl.aatt(1, "", 2, 0, 2)

        # ------------------------------------------------------------
        # MESH CONTROLS
        # ------------------------------------------------------------
        mapdl.allsel("ALL")
        mapdl.lsel("S", "LOC", "X", 0, 0)
        mapdl.lsel("U", "LOC", "Z", -le, -0.001)
        mapdl.lsel("U", "LOC", "Z", (le + 0.001), 2 * le)
        mapdl.lsel("U", "LOC", "Y", -2 * he, (-0.001 - he / 2))
        mapdl.lesize("ALL", nel, "", "", "", 1, "", "", 1)

        mapdl.allsel("ALL")
        mapdl.lsel("S", "LOC", "X", he, he)
        mapdl.lsel("A", "LOC", "Z", -he / 2, -he / 2)
        mapdl.lsel("A", "LOC", "Z", (le + he / 2), (le + he / 2))
        mapdl.lsel("A", "LOC", "Y", -he, -he)
        mapdl.lesize("ALL", 1.5 * nel, "", "", "", 1, "", "", 1)

        mapdl.allsel("ALL")
        mapdl.lsel("S", "LOC", "Z", (-he / 4), (-he / 4))
        mapdl.lesize("ALL", 1.25 * nel, "", "", 0.333, "", "", "", 1)

        mapdl.allsel("ALL")
        mapdl.lsel("S", "LOC", "Z", (le + he / 4), (le + he / 4))
        mapdl.lesize("ALL", 1.25 * nel, "", "", 3, "", "", "", 1)

        mapdl.allsel("ALL")
        mapdl.asel("S", "LOC", "X", -HS - 1e-6, 1e-6)
        mapdl.asel("R", "LOC", "Y", -he / 2, he / 2)
        mapdl.asel("R", "LOC", "Z", 0, le)

        mapdl.mshape(0, "2D")
        mapdl.mshkey(0)
        mapdl.amesh("ALL")

        # Mesh fluid volume.
        mapdl.allsel("ALL")
        mapdl.mshape(1, "3D")
        mapdl.mshkey(0)
        mapdl.vmesh("ALL")

        # ------------------------------------------------------------
        # PROFILE ROOT CONNECTION
        # Safe CPINTF: only structural DOFs, and avoids boundary nodes.
        # ------------------------------------------------------------
        cp_tol = max(nel * 0.45, 0.008)
        cp_guard = max(nel * 1.10, 0.030)

        y_min_cp = -he / 2.0 + cp_guard
        y_max_cp =  he / 2.0 - cp_guard

        for zp in profile_z:
            if y_min_cp >= y_max_cp:
                continue

            mapdl.allsel("ALL")
            mapdl.nsel("S", "LOC", "X", -1e-6, 1e-6)
            mapdl.nsel("R", "LOC", "Y", y_min_cp, y_max_cp)
            mapdl.nsel("R", "LOC", "Z", zp - cp_tol, zp + cp_tol)

            for dof in ["UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"]:
                try:
                    mapdl.run(f"CPINTF,{dof},{cp_tol}")
                except Exception:
                    pass

        mapdl.allsel("ALL")

        # ------------------------------------------------------------
        # BOUNDARY CONDITIONS
        # ------------------------------------------------------------
        mapdl.allsel("ALL")
        mapdl.lsel("S", "LOC", "X", 0, 0)
        mapdl.lsel("R", "LOC", "Y", -he / 2 - 0.01, he / 2 + 0.01)
        mapdl.lsel("R", "LOC", "Z", -0.01, le + 0.01)
        mapdl.lsel("U", "LOC", "Y", -he / 2 + 0.01, he / 2 - 0.01)
        mapdl.dl("ALL", "", "UX", 0)
        mapdl.dl("ALL", "", "UY", 0)
        mapdl.dl("ALL", "", "UZ", 0)
        mapdl.dl("ALL", "", "ROTX", 0)
        mapdl.dl("ALL", "", "ROTY", 0)
        mapdl.dl("ALL", "", "ROTZ", 0)

        mapdl.allsel("ALL")
        mapdl.lsel("S", "LOC", "X", 0, 0)
        mapdl.lsel("R", "LOC", "Y", -he / 2 - 0.01, he / 2 + 0.01)
        mapdl.lsel("R", "LOC", "Z", -0.01, le + 0.01)
        mapdl.lsel("U", "LOC", "Z", 0.01, le - 0.01)
        mapdl.dl("ALL", "", "UX", 0)
        mapdl.dl("ALL", "", "UY", 0)
        mapdl.dl("ALL", "", "UZ", 0)
        mapdl.dl("ALL", "", "ROTX", 0)
        mapdl.dl("ALL", "", "ROTY", 0)
        mapdl.dl("ALL", "", "ROTZ", 0)

        # FSI on wet structural/fluid interface
        mapdl.allsel("ALL")
        mapdl.nsel("S", "LOC", "X", 0, 0)
        mapdl.nsel("U", "LOC", "Y", 0.000001, he)
        mapdl.nsel("U", "LOC", "Y", (-0.000001 - he / 2), -he)
        mapdl.nsel("U", "LOC", "Z", (le + 0.000001), 2 * le)
        mapdl.nsel("U", "LOC", "Z", -he, -0.000001)
        mapdl.sf("", "FSI")

        # Free surface pressure condition
        mapdl.allsel("ALL")
        mapdl.nsel("S", "LOC", "Y", 0, 0)
        mapdl.nsel("U", "LOC", "X", -HS, -0.001)
        mapdl.d("ALL", "PRES", 0)

        # ------------------------------------------------------------
        # SOLUTION 
        # ------------------------------------------------------------
        mapdl.allsel("ALL")
        mapdl.slashsolu()
        mapdl.antype(2)
        mapdl.modopt("UNSYM", nmode)
        mapdl.mxpand(nmode, "", "", "YES")
        mapdl.lumpm(0)
        mapdl.pstres(0)
        mapdl.modopt("UNSYM", nmode, 0, 10000, "", "OFF")
        mapdl.solve()
        mapdl.finish()

        # ------------------------------------------------------------
        # POSTPROCESS
        # ------------------------------------------------------------
        freqs = _read_modal_frequencies(mapdl, nmode)
        selected_freq, selected_mode = _select_physical_frequency(freqs, freq_tol=freq_tol)

        mapdl.post1()
        if selected_mode is not None:
            try:
                mapdl.set(1, selected_mode)
            except Exception:
                pass
        else:
            try:
                mapdl.set("FIRST")
            except Exception:
                pass

        mapdl.allsel("ALL")
        mapdl.esel("S", "TYPE", "", 2)
        mapdl.run("ETABLE,EVOL,VOLU")
        mapdl.ssum()
        mapdl.get("TOTVOL", "SSUM", 0, "ITEM", "EVOL")
        volume = float(mapdl.parameters["TOTVOL"])
        mass = volume * RHO_S

        plot_file = None
        if visualize:
            print("\n[INFO] MAPDL-style mesh image is being generated. MAPDL GUI will not open.")
            plot_file = _static_mapdl_png(mapdl, workdir)
            _show_png_if_notebook(plot_file)

        elapsed = time.time() - start_time

        info = {
            "mode_frequencies": freqs,
            "selected_mode": selected_mode,
            "freq_tol": freq_tol,
            "n_profiles": n_profiles,
            "profile_z": profile_z,
            "t_plate": t_plate,
            "t_profile": t_profile,
            "NS": NS,
            "volume": volume,
            "mass": mass,
            "elapsed_time": elapsed,
            "workdir": workdir,
            "plot_file": plot_file,
        }

        if verbose:
            print(f"[RESULT] Selected wet frequency = {selected_freq:.6f} Hz")
            print(f"[RESULT] Selected mode          = {selected_mode}")
            print(f"[RESULT] Structural volume     = {volume:.9f} m^3")
            print(f"[RESULT] Mass                  = {mass:.6f} kg")
            print(f"[RESULT] Profiles              = {n_profiles} at z={profile_z}")
            print(f"[RESULT] Time                  = {elapsed:.2f} s")

        return selected_freq, volume, mass, info

    finally:
        if local_mapdl:
            try:
                mapdl.exit()
            except Exception:
                pass
