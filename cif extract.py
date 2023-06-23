import pymatgen as mg
import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os

# Provide your valid Materials Project API key
api_key = "aNniDNfUw1ikp6UBNyLIFjJtci44U0Bv"

# Create an instance of MPRester with your API key
mpr = MPRester(api_key)

# Perform the query and retrieve the data
data = mpr.summary.search(elements=['Fe', 'O', 'Mn', 'Mg'])

# Extract the material IDs
id_list = [entry.task_ids for entry in data]

# Print the list of material IDs
print(id_list)

# Initialize the Materials Project API
material_ids = ["mp-2217438"]

# Create a folder to store the XY files
folder_name = "cif_files"
os.makedirs(folder_name, exist_ok=True)

# Initialize the Materials Project API
# Initialize the Materials Project API
with MPRester(api_key) as mpr:
    for material_id in material_ids:
        try:
            # Retrieve the relevant structure
            structure = mpr.get_structure_by_material_id(material_id)

            # Obtain the formula
            formula = structure.composition.reduced_formula

            # Obtain the conventional structure
            sga = SpacegroupAnalyzer(structure)
            conventional_structure = sga.get_conventional_standard_structure()

            # Calculate the XRD diffraction pattern
            calculator = XRDCalculator(wavelength="CuKa")
            pattern = calculator.get_pattern(conventional_structure)

            # Generate the XY file path
            xy_file_path = os.path.join(folder_name, f"{formula}.xy")

            # Generate the XY file
            with open(xy_file_path, "w") as xy_file:
                xy_file.write("# 2-theta (degrees)\tIntensity\n")
                for angle, intensity in zip(pattern.x, pattern.y):
                    xy_file.write(f"{angle}\t{intensity}\n")

            print(f"Generated XY file for {material_id} ({formula}): {xy_file_path}")

        except Exception as e:
            print(f"Skipping material ID {material_id} due to error: {e}")
            continue
