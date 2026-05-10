import json
from pathlib import Path

from ase.io import read

from guesthost.unitcells import (
    unit_order_from_unitcell_data,
    unit_orders_from_unitcell_data,
)


_PACKAGE_DIR = Path(__file__).resolve().parent
_STRUCTURES_DIR = _PACKAGE_DIR / "data" / "structures"
_UNITCELLS_DIR = _PACKAGE_DIR / "data" / "unitcells"


def _load_json(name):
    with open(_UNITCELLS_DIR / name, "r") as f:
        return json.load(f)


def _convert_entry(data):
    out = dict(data)
    if "pb_axis" in out:
        out["pb_axis"] = (out["pb_axis"][0], out["pb_axis"][1])
    return out


MPB_SYS_4x4x4 = read(_STRUCTURES_DIR / "mpb_cubic_4x4x4.xyz", index=0)
MPB_SYS_8x8x8 = read(_STRUCTURES_DIR / "mpb_cubic_8x8x8.xyz", index=0)

UNITCELL_INDEXDATA_MPB_4x4x4 = [
    _convert_entry(data)
    for data in _load_json("unitcell_indexdata_mpb_4x4x4.json")
]
UNITCELL_INDEXDATA_MPB_8x8x8 = [
    _convert_entry(data)
    for data in _load_json("unitcell_indexdata_mpb_8x8x8.json")
]

UNIT_ORDER_MPB_4x4x4 = unit_order_from_unitcell_data(UNITCELL_INDEXDATA_MPB_4x4x4)
UNIT_ORDERS_MPB_4x4x4 = unit_orders_from_unitcell_data(UNITCELL_INDEXDATA_MPB_4x4x4)
UNIT_ORDER_MPB_8x8x8 = unit_order_from_unitcell_data(UNITCELL_INDEXDATA_MPB_8x8x8)
UNIT_ORDERS_MPB_8x8x8 = unit_orders_from_unitcell_data(UNITCELL_INDEXDATA_MPB_8x8x8)

MA_INDICES_MPB = [0, 1, 2, 3, 4, 5, 6, 7]
HOST_INDICES_MPB = [8, 9, 10, 11]
