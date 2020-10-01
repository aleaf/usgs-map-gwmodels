import os
from pathlib import Path
import pandas as pd


def makedirs(path):
    dirs = [os.path.join(path, 'shps/'),
            os.path.join(path, 'figures/')]
    for folder in dirs:
        if not os.path.isdir(folder):
            os.makedirs(folder)
