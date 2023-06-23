import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure
from scipy.optimize import curve_fit

def process_cif_file(cif_file):
    # Load CIF file
    structure = Structure.from_file(cif_file)

    # Set up XRD calculator
    xrd_calculator = XRDCalculator()

    # Calculate XRD pattern
    xrd_pattern = xrd_calculator.get_pattern(structure)

    # Define the desired x range and step size
    x_range = np.arange(10, 80, 0.01)

    # Extract the x and y data from the XRD pattern
    x_data = xrd_pattern.x
    y_data = xrd_pattern.y

    # Initialize the y values with zeros
    y_values = np.zeros_like(x_range)

    # Match the x values from the XRD pattern with the desired x range
    for i, x_val in enumerate(x_range):
        closest_idx = np.argmin(np.abs(x_data - x_val))

        # Check if the closest x value is within a certain threshold
        if np.abs(x_data[closest_idx] - x_val) <= 0.01:
            y_values[i] = y_data[closest_idx]

    # Perform further processing or plotting here
    # ...

    return x_range, y_values

# Directory containing CIF files
cif_directory = 'COD(1)'

# Get a list of CIF files in the directory
cif_files = os.listdir(cif_directory)

# Iterate over CIF files
for cif_file in cif_files:
    if cif_file.endswith('.cif'):
        cif_path = os.path.join(cif_directory, cif_file)
        x_range, y_values = process_cif_file(cif_path)
