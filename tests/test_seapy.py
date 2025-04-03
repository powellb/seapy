import os
import pytest
import datetime as dt
from netCDF4 import Dataset

from seapy.roms.forcing import gen_bulk_forcing, forcing_data
from seapy.roms.interp import to_clim


tstart = dt.datetime(2019, 8, 13)
tend = dt.datetime(2019, 8, 14)

_here = os.path.dirname(os.path.abspath(__file__))
METEOFILE = os.path.join(_here, "meteo.nc")
OCEANFILE = os.path.join(_here, "ocean.nc")
GRDFILE = os.path.join(_here, "roms_grd.nc")
FRCFILE = os.path.join(_here, "roms_frc.nc")
CLMFILE = os.path.join(_here, "roms_clm.nc")
PMAPFILE = os.path.join(_here, "../ocean_roms_grd_pmap.npz")

OUTPUTS = [
    FRCFILE,
    CLMFILE,
]


@pytest.fixture
def cleanup_outputs():
    """Fixture to delete outputs after the test."""

    # Yield control to the test
    yield  # This is where the test runs

    # Teardown: Delete the file after the test
    for file_path in OUTPUTS:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"Deleted {file_path}")


def test_gen_bulk_forcing():
    """Test the gen_bulk_forcing function."""
    # Define test parameters
    gen_bulk_forcing(
        METEOFILE,
        {
            "pad": 1.0,
            "frc_lat": "latitude",
            "frc_lon": "longitude",
            "frc_time": "time",
            "Tair": forcing_data(field="tmp2m", ratio=1, offset=-273.15),
            "Pair": forcing_data(field="mslp", ratio=0.01, offset=0),
            "Qair": forcing_data(field="rh2m", ratio=1, offset=0),
            "rain": forcing_data(field="pratesfc", ratio=1, offset=0),
            "Uwind": forcing_data(field="ugrd10m", ratio=1, offset=0),
            "Vwind": forcing_data(field="vgrd10m", ratio=1, offset=0),
            "lwrad_down": forcing_data(field="dlwrfsfc", ratio=1, offset=0),
            "swrad": forcing_data(field="dswrfsfc", ratio=1, offset=0),
        },
        FRCFILE,
        GRDFILE,
        tstart,
        tend,
        epoch=tstart,
        clobber=True,
    )

    # check if output file exists
    assert os.path.exists(FRCFILE), f"Output file {FRCFILE} does not exist."

    # check if output file is a valid netCDF file
    with Dataset(FRCFILE, "r") as nc:
        assert nc is not None

    os.remove(FRCFILE)


def test_to_clim():
    """Test the to_clim function."""
    # Define test parameters
    pmap = to_clim(
        OCEANFILE,
        CLMFILE,
        dest_grid=GRDFILE,
        reftime=tstart,
        vmap={"ssh": "zeta", "uo": "u", "vo": "v", "temp": "temp", "salt": "salt"},
        nx=0.16,
        ny=0.16,
        threads=2,
        clobber=True,
    )

    # check if pmap is not None
    assert pmap is not None, "pmap should not be None."

    # check if output file exists
    assert os.path.exists(CLMFILE), f"Output file {CLMFILE} does not exist."

    # check if output file is a valid netCDF file
    with Dataset(CLMFILE, "r") as nc:
        assert nc is not None

    os.remove(CLMFILE)
    os.remove(PMAPFILE)
