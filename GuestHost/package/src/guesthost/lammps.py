import os
import re
from pathlib import Path
from typing import List, Sequence, Tuple

from ase.io import read

from guesthost.lattice import LatticeTrajectory
from guesthost.trajectory import Trajectory


def _extract_timesteps_from_dump(trajfile: Path) -> List[int]:
    steps = []
    with trajfile.open() as handle:
        for line in handle:
            if line.startswith("ITEM: TIMESTEP"):
                step_line = next(handle, "").strip()
                try:
                    steps.append(int(step_line))
                except ValueError:
                    continue
    return steps


def _normalize_frame_indices(steps: Sequence | slice | int | None, nframes: int) -> List[int]:
    if steps is None:
        indices = [0]
    elif isinstance(steps, int):
        indices = [steps - 1]
    elif isinstance(steps, range):
        indices = [i - 1 for i in steps]
    elif isinstance(steps, slice):
        start = 0 if steps.start is None else steps.start - 1
        stop = nframes if steps.stop is None else steps.stop
        step = 1 if steps.step is None else steps.step
        indices = list(range(start, stop, step))
    else:
        indices = [int(i) - 1 for i in steps]
    return [i for i in indices if 0 <= i < nframes]


def _parse_log_times(logpath: Path) -> Tuple[List[int], List[float]]:
    if not logpath.exists():
        return [], []

    steps = []
    times = []
    with logpath.open() as handle:
        for line in handle:
            if not line or not line[0].isdigit():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                step = int(parts[0])
                time = float(parts[1])
            except ValueError:
                continue
            steps.append(step)
            times.append(time)
    return steps, times


def _select_times(step_numbers, log_steps, log_times, timestep_fs):
    if log_steps and len(log_steps) == len(log_times):
        mapping = dict(zip(log_steps, log_times))
        if all(step in mapping for step in step_numbers):
            return [mapping[step] for step in step_numbers]

    scale = timestep_fs * 1e-3
    return [step * scale for step in step_numbers]


def _resolve_unitcell_data(unitcell_data=None, unitdata=None):
    if unitcell_data is not None:
        return unitcell_data
    if unitdata is not None:
        return unitdata
    raise TypeError("Pass unitcell_data or unitdata.")


def load_lammps_trajectory(
    trajfile,
    *,
    unitcell_data=None,
    unitdata=None,
    supercell_size=None,
    steps: Sequence | slice | int | None = None,
    timestep_fs=0.1,
    prefer_logtime=True,
):
    """
    Load a LAMMPS dump and return a LatticeTrajectory of HPLattice frames.

    Frame selections use 1-based frame numbers, matching the perovskite_analysis
    API.
    """
    unitcell_data = _resolve_unitcell_data(unitcell_data=unitcell_data, unitdata=unitdata)
    traj_path = Path(trajfile)
    dump_steps = _extract_timesteps_from_dump(traj_path)
    frame_indices = _normalize_frame_indices(steps, len(dump_steps))
    if not frame_indices:
        raise ValueError("No frames selected from trajectory.")

    if len(frame_indices) == 1:
        atoms = read(traj_path, index=frame_indices[0])
        frames = atoms if isinstance(atoms, list) else [atoms]
    else:
        all_frames = read(traj_path, index=":")
        frames = [all_frames[i] for i in frame_indices]

    selected_steps = [dump_steps[i] for i in frame_indices]
    if prefer_logtime:
        log_steps, log_times = _parse_log_times(traj_path.with_name("log.lammps"))
    else:
        log_steps, log_times = [], []
    times = _select_times(selected_steps, log_steps, log_times, timestep_fs)

    trj = Trajectory.from_atoms_list(frames)
    lattices = trj.create_hplattice(
        supercell_size=supercell_size,
        unitcell_data=unitcell_data,
    )
    ltraj = LatticeTrajectory(lattices)
    ltraj.times = times
    ltraj.steps = selected_steps
    ltraj.data = {"stepnumber": selected_steps, "dumppath": str(traj_path)}
    ltraj.frames = frames
    return ltraj


def load_lammps_trajectories(
    maindir,
    *,
    unitcell_data=None,
    unitdata=None,
    supercell_size=None,
    mddir_regex=r"^0[1-6]",
    steps: Sequence | slice | int | None = None,
    dump_regex=r"prod\.mpb$",
    prmatcher=r"[0-9]+_npt([0-9.]+)gpa$",
    append_to_all=None,
    timestep_fs=0.1,
    prefer_logtime=True,
):
    """Load multiple LAMMPS dumps from matching subdirectories."""
    unitcell_data = _resolve_unitcell_data(unitcell_data=unitcell_data, unitdata=unitdata)
    base = Path(maindir)
    dirs = [base / d for d in os.listdir(base) if re.search(mddir_regex, d)]
    if not dirs:
        raise FileNotFoundError(f"No directories matching {mddir_regex} in {maindir}")

    trajectories = {}
    for dpath in dirs:
        workdir = dpath / append_to_all if append_to_all else dpath
        candidates = [
            workdir / name for name in os.listdir(workdir) if re.search(dump_regex, name)
        ]
        if not candidates:
            raise FileNotFoundError(f"No dump matching {dump_regex} in {workdir}")

        match = re.search(prmatcher, str(dpath))
        key = match.group(1) if match else dpath.name
        trajectories[str(key)] = load_lammps_trajectory(
            candidates[0],
            unitcell_data=unitcell_data,
            supercell_size=supercell_size,
            steps=steps,
            timestep_fs=timestep_fs,
            prefer_logtime=prefer_logtime,
        )

    return trajectories
