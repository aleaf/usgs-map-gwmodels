import os
from pathlib import Path
import shutil
import numpy as np
import pytest
from mfsetup.grid import MFsetupGrid


@pytest.fixture(scope="session")
def project_root_path():
    """Root folder for the project (with setup.py),
    two levels up from the location of this file.
    """
    filepath = os.path.split(os.path.abspath(__file__))[0]
    normpath = os.path.normpath(os.path.join(filepath, '..', '..'))
    return Path(normpath)


@pytest.fixture(scope="session")
def test_data_path(project_root_path):
    """Root folder for the project (with setup.py),
    two levels up from the location of this file.
    """
    return Path(project_root_path, 'mapgwm', 'tests', 'data')


@pytest.fixture(scope="session", autouse=True)
def test_output_folder(project_root_path):
    """(Re)make an output folder for the tests
    at the begining of each test session."""
    folder = os.path.join(project_root_path, 'mapgwm', 'tests', 'output')
    reset = False
    if reset:
        if os.path.isdir(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)
    return Path(folder)


@pytest.fixture(scope='module')
def delta_inset_model_grid():
    ncol=270
    nrow=605
    dxy=500.
    delr = np.ones(ncol) * dxy
    delc = np.ones(nrow) * dxy
    grid = MFsetupGrid(delc, delr, top=None, botm=None,
                       lenuni=2, epsg=5070,
                       xoff=434955, yoff=1040285, angrot=0.0)
    return grid