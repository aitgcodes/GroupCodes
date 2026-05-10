from importlib.resources import files

import guesthost as gh


def test_load_system_from_xyz():
    data_path = files("guesthost").joinpath("data", "structures", "mpb_cubic_4x4x4.xyz")
    system = gh.load_system_from_xyz(data_path)

    assert len(system) == 768
    assert system.get_chemical_symbols()[0] == "Pb"
