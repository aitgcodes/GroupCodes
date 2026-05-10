# guesthost

`guesthost` provides tools for analyzing guest-host atomic systems, with helper
classes for molecular fragments, host fragments, unit-cell motifs, lattices, and
lattice trajectories.

## Contents

- [Installation](#installation)
- [Data Model](#data-model)
- [Basic Usage](#basic-usage)
- [Unit-Cell Methods](#unit-cell-methods)
  - [Local Cell Geometry](#local-cell-geometry)
  - [MA Orientation And Order](#ma-orientation-and-order)
  - [Host And Cage Metrics](#host-and-cage-metrics)
  - [Index And Utility Methods](#index-and-utility-methods)
- [Lattice-Wide Methods](#lattice-wide-methods)
- [Layerwise And Global Order Parameters](#layerwise-and-global-order-parameters)
- [Polar Order](#polar-order)
- [Global Orthorhombic Cell Helpers](#global-orthorhombic-cell-helpers)
- [Compute And Results IO](#compute-and-results-io)
- [LAMMPS Trajectories](#lammps-trajectories)
- [Unit-Cell Detection](#unit-cell-detection)
- [Testing](#testing)

## Installation

From the package root:

```bash
pip install -e .
```

## Data Model

The main objects are:

- `Trajectory`: reads ASE-supported trajectory files and stores frames, atomic
  coordinates, cells, symbols, and frame counts.
- `MethylAmmonium`: stores the C, N, H-C, H-N coordinates and the MA `CN`
  orientation vector for one molecule.
- `Host`: stores the central host atom and three in-cell host-environment atoms.
- `Motif`: groups the fragments belonging to one unit cell. For hybrid
  perovskite systems, `motif.ma` accesses the MA fragment and `motif.host`
  accesses the host fragment.
- `HPLattice`: stores motifs on a 3D lattice and provides unit-cell and
  lattice-level analysis methods.
- `LatticeTrajectory`: stores one lattice object per trajectory frame and
  provides trajectory-level compute helpers.

## Basic Usage

```python
from importlib.resources import files
import guesthost as gh

data_path = files("guesthost").joinpath("data", "structures", "mpb_trajectory.xyz")

trj = gh.Trajectory(str(data_path), order=True)

lattices = trj.create_hplattice(
    unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
    supercell_size=(4, 4, 4),
)

lattice = lattices[1]
ind = (0, 0, 0)

theta, phi = lattice.ucell_theta_phi(ind, dir=0)
params = lattice.ucell_localcellparameters(ind)
eta_ma = lattice.ucell_eta_MA(ind, dir=0)
eta_lat = lattice.ucell_eta_lat(ind, dir_coup=0)
```

Individual XYZ structures can be loaded with ASE through the package helper:

```python
from importlib.resources import files
import guesthost as gh

data_path = files("guesthost").joinpath("data", "structures", "mpb_cubic_4x4x4.xyz")
system = gh.load_system_from_xyz(data_path)
```

Unit-cell data can also be inferred from a reference structure:

```python
lattices = trj.create_hplattice(
    reference=trj.atoms_list[0],
    supercell_size=(4, 4, 4),
)
```

Packaged MAPbBr3 reference structures and unit-cell data are available as
constants:

```python
system_4 = gh.MPB_SYS_4x4x4
unitcells_4 = gh.UNITCELL_INDEXDATA_MPB_4x4x4

system_8 = gh.MPB_SYS_8x8x8
unitcells_8 = gh.UNITCELL_INDEXDATA_MPB_8x8x8
```

## Unit-Cell Methods

`HPLattice` methods operate on unit-cell indices of the form `(ix, iy, iz)`.
Directions use zero-based axes: `dir=0, 1, 2` and `dir_coup=0, 1, 2`.

```python
ind = (0, 0, 0)

theta, phi = lattice.ucell_theta_phi(ind, dir=0)
vecs = lattice.ucell_localcellvecs(ind)
params = lattice.ucell_localcellparameters(ind)
angle = lattice.ucell_pb_br_pb_angle(ind, dir=0)
angle_hostrelative = lattice.ucell_pb_br_pb_angle_hostrelative(ind, dir=0, dir_coup=1)
eta_ma = lattice.ucell_eta_MA(ind, dir=0)
eta_lat = lattice.ucell_eta_lat(ind, dir_coup=0)
```

### Local Cell Geometry

| Method | Description |
| --- | --- |
| `ucell_localcellvecs(ind)` | Local Pb-Pb cell vectors for a unit cell. |
| `ucell_localcellparameters(ind)` | Local cell lengths, angles, and volume as a dictionary. |
| `ucell_localcelllengths(ind)` | Local cell vector lengths as a NumPy array. |
| `ucell_localcelllengths_only(ind)` | Local cell vector lengths as a list. |
| `ucell_localcellangles(ind)` | Local cell angles `[alpha, beta, gamma]`. |
| `ucell_localcell_lengths_angles(ind)` | Pair of local cell lengths and angles. |
| `ucell_localcellvecs_ortho(ind, dir, small=False)` | Orthorhombic-like local cell vectors for a coupled direction. |
| `ucell_localcelllengths_ortho(ind, dir, small=False)` | Orthorhombic-like local cell vector lengths. |
| `ucell_localcellangles_ortho(ind, dir, small=False)` | Orthorhombic-like local cell angles. |
| `ucell_localcell_lengths_angles_ortho(ind, dir, small=False)` | Pair of orthorhombic-like lengths and angles. |

### MA Orientation And Order

| Method | Description |
| --- | --- |
| `ucell_theta_phi(ind, dir)` | MA orientation angles `(theta, phi)` relative to the local cell. |
| `ucell_theta_phi_dict(ind, dir)` | MA orientation angles as `{"theta": ..., "phi": ...}`. |
| `ucell_theta_phi_withcutoff(ind, dir, ratio_cut=0.05)` | MA orientation with `phi` set to zero when the in-plane projection is small. |
| `ucell_xyz(theta_deg, phi_deg, theta_dir=2)` | Convert orientation angles to Cartesian-like components. |
| `ucell_eta_MA(ind, dir, **kwargs)` | Two-component MA order parameter. |
| `ucell_eta_MA_from_angles(theta, phi, **kwargs)` | Two-component MA order parameter from angles. |
| `ucell_eta_3d_MA(ind, dir=2, **kwargs)` | Three-component MA order parameter. |
| `ucell_eta_3d_MA_from_angles(theta, phi, **kwargs)` | Three-component MA order parameter from angles. |
| `ucell_MA_dihedrals(ind)` | Nine MA H-C-N-H dihedral angles. |
| `ucell_CN_distance(ind)` | C-N bond distance. |
| `ucell_CH_distances(ind)` | C-H distances for the C-side hydrogens. |
| `ucell_NH_distances(ind)` | N-H distances for the N-side hydrogens. |

### Host And Cage Metrics

| Method | Description |
| --- | --- |
| `ucell_br_br_vectors(ind, dir=0)` | Br(-axis) to Br(+axis) cage vectors. |
| `ucell_br_br_distances(ind, dir=0)` | Lengths of the Br-Br cage vectors. |
| `ucell_tetragonal_lengthening_distortion_br(ind, dir=0)` | Br cage tetragonal lengthening distortion. |
| `ucell_scissoring_distortion_br(ind, dir=0)` | Br cage scissoring distortion. |
| `ucell_scissoring_distortion_angle_br(ind, dir=0)` | Br cage scissoring distortion angle. |
| `ucell_pb_br_pb_angle(ind, dir)` | Pb-Br-Pb angle along a local axis. |
| `ucell_pb_br_pb_angle_hostrelative(ind, dir, dir_coup)` | Pb-Br-Pb angle with host-relative sign convention. |
| `ucell_br_distance_from_plane(ind, dir_pl, dir_br)` | Signed Br distance from a Pb plane. |
| `ucell_pb_br_pb_angle_alongaxis(ind, dir, dir_coup)` | Pb-Br-Pb angle projected along a coupled axis. |
| `ucell_pb_br_pb_angle_alongaxis_hostrelative(ind, dir, dir_coup)` | Projected Pb-Br-Pb angle with host-relative sign convention. |
| `ucell_eta_lat(ind, dir_coup, ref=180.0, N=1.0)` | Two-component lattice order parameter from host-relative Pb-Br-Pb angles. |
| `ucell_eta_lat_alongaxis(ind, dir_coup, ref=180.0, N=1.0)` | Lattice order parameter from projected Pb-Br-Pb angles. |
| `ucell_ma_displacement_vecs_br(ind)` | MA center displacement relative to the Br cage center. |

### Index And Utility Methods

| Method | Description |
| --- | --- |
| `ucell_pb_axis(ind, python_index=False)` | Pb origin and positive-axis Pb indices for the unit cell. |
| `ucell_udata(ind)` | Unit-cell index metadata associated with a motif. |
| `ucell_all_octahedra_indices(ind)` | Three `(Pb, Br, Pb)` axis triplets for the local octahedron. |
| `ucell_ω(eta, R, k_vec=None, orig_cell=(0, 0), MA=False)` | Rotate an order parameter into the omega basis. |

Unicode method names are also available for compatibility:

```python
eta_ma = lattice.ucell_η_MA(ind, dir=0)
eta_lat = lattice.ucell_η_lat(ind, dir_coup=0)
eta_ma_3d = lattice.ucell_η_3d_MA(ind, dir=0)
```

## Lattice-Wide Methods

These methods compute unit-cell quantities for every motif in the lattice.

```python
theta, phi = lattice.all_ucell_theta_phi(dir=0)
eta_ma = lattice.all_ucell_eta_MA(dir=0)
eta_lat = lattice.all_ucell_eta_lat(dir_coup=0)
```

| Method | Description |
| --- | --- |
| `all_ucell_theta_phi(dir)` | MA orientation angles for every unit cell. |
| `all_ucell_eta_MA(dir, **kwargs)` | Two-component MA order parameter for every unit cell. |
| `all_ucell_eta_lat(dir_coup, ref=180.0, N=1.0)` | Two-component lattice order parameter for every unit cell. |

For trajectory-level analysis:

```python
lattraj = gh.LatticeTrajectory(lattices)
theta_phi_per_frame = lattraj.compute(gh.HPLattice.ucell_theta_phi, dir=0)
xi_ma_per_frame = lattraj.compute_layerwise(gh.global_xi_MA, dir=0)
```

## Layerwise And Global Order Parameters

Layerwise functions operate on a complete `HPLattice`. They reduce unit-cell
quantities into layer-resolved or system-wide order parameters.

```python
omega_lat = gh.layerwise_OmegaLat(lattice, dir_coup=0)
eta_lat_layers = gh.layerwise_eta_Lat(lattice, dir_coup=0)
xi_lat_layers = gh.global_layerwise_xi(lattice, dir_coup=0)

omega_ma = gh.layerwise_OmegaMA(lattice, dir=0)
eta_ma_layers = gh.layerwise_eta_MA(lattice, dir=0)
xi_ma = gh.global_xi_MA(lattice, dir=0)

omega_coupling = gh.layerwise_OmegaCoupling(lattice, dir=0)
xi_coupling = gh.global_xi_coupling(lattice, dir=0)
volume = gh.global_volume_local(lattice)
```

Available layerwise/global helpers:

| Function | Description |
| --- | --- |
| `compute_all_ucells(func, lattice, **kwargs)` | Apply a unit-cell callable to every unit cell and return a flat list. |
| `compute_all_ucells_matrix(func, lattice, order=None, cellshape=None, **kwargs)` | Apply a unit-cell callable and return an object matrix. |
| `compute_all_ucells_data(func, lattice, order=None, cellshape=None, **kwargs)` | Apply a unit-cell callable and keep arbitrary return data. |
| `layerwise_OmegaLat(lattice, dir_coup, ...)` | Layer-resolved lattice omega order parameter. |
| `layerwise_eta_Lat(lattice, dir_coup, ...)` | Layer-resolved scalar lattice eta projection. |
| `global_layerwise_xi(lattice, dir_coup, ...)` | Layer-resolved xi from lattice omega values. |
| `global_layerwise_eta_alloctahedra(lattice, dir_coup, ...)` | Layer-resolved octahedral eta values. |
| `global_layerwise_xi_alloctahedra(lattice, dir_coup, ...)` | Layer-resolved octahedral xi matrices. |
| `global_layerwise_xi_new(lattice, dir_coup, ...)` | Mean xi for each layer. |
| `layerwise_OmegaMA(lattice, dir, ...)` | Layer-resolved MA omega order parameter. |
| `layerwise_eta_MA(lattice, dir, ...)` | Layer-resolved scalar MA eta projection. |
| `global_xi_MA(lattice, dir, ...)` | Global MA xi order parameter. |
| `layerwise_OmegaCoupling(lattice, dir, ...)` | Layer-resolved lattice-MA coupling omega. |
| `global_xi_coupling(lattice, dir, ...)` | Global lattice-MA coupling xi. |
| `global_volume_local(lattice)` | Mean local cell volume. |
| `global_xi(lattice, dir_coup, ...)` | Global lattice xi order parameter. |
| `global_S(lattice, dir_coup, ...)` | Layer autocorrelation-like `S` values. |

## Polar Order

The polar order parameter is derived from nearest-neighbor differences in MA
`phi` angles along a coupling direction.

```python
polar = gh.polar_order_parameter(lattice, dir_coup=0)

dphi = polar["dphi"]
z = polar["z"]
chain_order = polar["S_chain"]
```

The result dictionary contains `dphi`, `z`, `z_pi`, `z_pi_by2`,
`z_chain`, `z_chain_pi`, `z_chain_pi_by2`, `S_chain`, `S_chain_pi`, and
`S_chain_pi_by2`.

## Global Orthorhombic Cell Helpers

Global orthorhombic cell helpers describe the whole lattice cell, not a single
unit cell.

```python
global_data = gh.get_global_index_data(lattice)

vecs = gh.global_cell_vectors_OR(lattice, dir=0, global_data=global_data)
lengths, angles = gh.global_cell_lengths_angles_OR(lattice, dir=0)

lengths_from_cell = gh.global_cell_lengths_OR(lattice, dir=0, fromcell=True)
angles_from_cell = gh.global_cell_angles_OR(lattice, dir=0, fromcell=True)
```

## Compute And Results IO

The compute helpers accept `HPLattice` methods or callables shaped like
`func(lattice, ind, **kwargs)`.

```python
values = gh.compute(lattice, gh.HPLattice.ucell_localcelllengths)

lattraj = gh.LatticeTrajectory(lattices[:2])
theta_phi = gh.compute_trajectory(lattraj, gh.HPLattice.ucell_theta_phi, dir=0)

bundle = gh.compute_functions(
    lattraj,
    [
        (gh.HPLattice.ucell_theta_phi, {"dir": 0}, False),
        (gh.global_xi_MA, {"dir": 0}, True),
    ],
)

default_bundle = gh.compute_default_functions(lattice, cellshape=(4, 4, 4))
```

Results can be saved to and loaded from HDF5:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
import guesthost as gh

with TemporaryDirectory() as tmpdir:
    out = Path(tmpdir) / "analysis_results.h5"
    gh.save_results(out, {"xi_ma": bundle["layerwise:global_xi_MA(dir=0)"]})
    loaded = gh.load_results(out)
```

If a `LatticeTrajectory` has `times` or `steps` attributes, `save_results`
writes them as `time` and `steps` datasets.

## LAMMPS Trajectories

LAMMPS dump files can be loaded directly into a `LatticeTrajectory`.
Frame selections use one-based frame numbers.

```python
from pathlib import Path
import guesthost as gh

pkg_root = Path(gh.__file__).resolve().parents[2]
dump_path = pkg_root / "tests" / "structures" / "lammps_trajectory" / "trajectory.lmp"

lattraj = gh.load_lammps_trajectory(
    dump_path,
    unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
    supercell_size=(4, 4, 4),
    steps=range(1, 3),
)

theta_phi = lattraj.compute(gh.HPLattice.ucell_theta_phi, dir=0)
xi_ma = lattraj.compute_layerwise(gh.global_xi_MA, dir=0)
```

Multiple dump directories can be loaded with `load_lammps_trajectories`:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
import guesthost as gh

pkg_root = Path(gh.__file__).resolve().parents[2]
fixture_dir = pkg_root / "tests" / "structures" / "lammps_trajectory"

with TemporaryDirectory() as tmpdir:
    run_dir = Path(tmpdir) / "01_npt1.1gpa"
    run_dir.mkdir()
    for name in ("trajectory.lmp", "log.lammps"):
        (run_dir / name).write_bytes((fixture_dir / name).read_bytes())

    trajs = gh.load_lammps_trajectories(
        tmpdir,
        unitcell_data=gh.UNITCELL_INDEXDATA_MPB_4x4x4,
        supercell_size=(4, 4, 4),
        mddir_regex=r"^01",
        dump_regex=r"trajectory\.lmp$",
        steps=1,
    )
```

For production data, point `load_lammps_trajectories` at the directory that
contains the matching run subdirectories and adjust `mddir_regex` and
`dump_regex` to match the local naming scheme.


## Unit-Cell Detection

`get_unitcell_indexdata` extracts unit-cell index metadata from an ASE system.
By default it returns the fields used by the cached package constants.

```python
unitcell_data = gh.get_unitcell_indexdata(gh.MPB_SYS_4x4x4)
```

Optional orthorhombic Pb-axis metadata can be requested when needed:

```python
unitcell_data = gh.get_unitcell_indexdata(gh.MPB_SYS_4x4x4, ortho_axes=True)

first = unitcell_data[0]
ortho_axes = first["pb_ortho_axes"]
ortho_axes_small = first["pb_ortho_axes_small"]
```

## Testing

From the package root:

```bash
python -m pytest tests
```

The tests include reference-data checks for lattice construction, orientations,
and torsions.
