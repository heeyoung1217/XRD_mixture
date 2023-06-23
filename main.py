import pymatgen as mg
import os
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.ndimage import convolve1d
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.signal import convolve
from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Step 3: Mixture - Generate Mixture Peak
augmented_folder = 'augmented_data'  # Folder containing augmented data
mixture_folder = 'mixture_data'  # Folder to save mixture data

if not os.path.exists(mixture_folder):
    os.makedirs(mixture_folder)

num_selected = 10  # Number of augmented spectra to select for mixture

augmented_files = os.listdir(augmented_folder)

for file_name in augmented_files:
    if file_name.endswith('.xy'):
        base_name = file_name.split('_augmented_')[0]

        # Select num_selected augmented spectra for mixture
        selected_files = [f for f in augmented_files if f.startswith(base_name)][:num_selected]
        spectra = []

        for selected_file in selected_files:
            file_path = os.path.join(augmented_folder, selected_file)
            spectrum = np.loadtxt(file_path, delimiter='\t')
            spectra.append(spectrum)

        # Generate mixture peak
        mixture_spectrum = np.mean(spectra, axis=0)

        # Save mixture peak as a xy.file
        save_path = os.path.join(mixture_folder, f'{base_name}_mixture.xy')
        np.savetxt(save_path, mixture_spectrum, delimiter='\t', fmt='%.2f')