import pymatgen as mg
import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

import os
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Load the XRD spectrum data (.xy)
xy_folder = 'xy_files'  # Folder containing .xy files
xy_numpy_folder = 'xy_file_numpy'  # Folder to save numpy data

if not os.path.exists(xy_numpy_folder):
    os.makedirs(xy_numpy_folder)

xy_files = os.listdir(xy_folder)

for file_name in xy_files:
    if file_name.endswith('.xy'):
        file_path = os.path.join(xy_folder, file_name)
        data = np.genfromtxt(file_path, skip_header=1)  # Load data skipping the first row

        x_data = data[:, 0]  # Extract x data from the first column
        y_data = data[:, 1]  # Extract y data from the second column

        # Filter x and y data within the specified x range (10-80 with 0.01 step)
        x_range = np.arange(10, 80, 0.01)
        filtered_y_data = np.zeros(len(x_range))

        for i, x in enumerate(x_range):
            index = np.where(np.isclose(x_data, x))[0]
            if len(index) > 0:
                filtered_y_data[i] = y_data[index[0]]

        # Save the numpy data
        save_path = os.path.join(xy_numpy_folder, f'{file_name.split(".")[0]}.npy')
        np.save(save_path, np.vstack((x_range, filtered_y_data)).T)

# Step 2: Data Augmentation
xy_numpy_folder = 'xy_file_numpy'  # Folder containing numpy data
augmented_folder = 'augmented_data'  # Folder to save augmented data

if not os.path.exists(augmented_folder):
    os.makedirs(augmented_folder)

strain_range = np.linspace(-0.04, 0.04, num=70)  # Generate 70 strain values within the range of +/- 4%

xy_files = os.listdir(xy_numpy_folder)

for file_name in xy_files:
    if file_name.endswith('.npy'):
        file_path = os.path.join(xy_numpy_folder, file_name)
        data = np.load(file_path)

        x_data = data[:, 0]  # Extract x data
        y_data = data[:, 1]  # Extract y data

        augmented_data = []

        for strain in strain_range:
            shifted_x_data = x_data * (1 + strain)  # Shift x values by strain percentage
            augmented_data.append(np.column_stack((shifted_x_data, y_data)))  # Combine shifted x and y data

        # Save each augmented spectrum as a xy.file
        base_name = os.path.splitext(file_name)[0]
        for i, spectrum in enumerate(augmented_data):
            save_path = os.path.join(augmented_folder, f'{base_name}_augmented_{i}.xy')
            np.savetxt(save_path, spectrum, delimiter='\t', fmt='%.2f')

# Step 3: Generate Gaussian peaks
augmented_folder = 'augmented_data'  # Folder containing augmented data
gaussian_folder = 'gaussian_peaks'  # Folder to save the generated Gaussian peaks

if not os.path.exists(gaussian_folder):
    os.makedirs(gaussian_folder)

fwhm = 0.015  # FWHM of the Gaussian peaks

xy_files = os.listdir(augmented_folder)

for file_name in xy_files:
    if file_name.endswith('.xy'):
        file_path = os.path.join(augmented_folder, file_name)
        data = np.genfromtxt(file_path)

        x_data = data[:, 0]  # Extract x data
        y_data = data[:, 1]  # Extract y data

        gaussian_peaks = np.zeros_like(x_data)

        # Generate Gaussian peaks for y data values over 1
        peak_indices = np.where(y_data > 1)[0]

        for i in peak_indices:
            peak_center = x_data[i]
            peak_amplitude = y_data[i]
            sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Calculate standard deviation (sigma) from FWHM
            gaussian_peak = peak_amplitude * np.exp(-0.5 * ((x_data - peak_center) / sigma) ** 2)
            gaussian_peaks += gaussian_peak

        # Save the generated Gaussian peaks as a xy.file
        base_name = os.path.splitext(file_name)[0]
        save_path = os.path.join(gaussian_folder, f'{base_name}_gaussian.xy')
        np.savetxt(save_path, np.column_stack((x_data, gaussian_peaks)), delimiter='\t', fmt='%.2f')

        # Optional: Plot the original and generated Gaussian peaks
        plt.plot(x_data, y_data, label='Original Spectrum')
        plt.plot(x_data, gaussian_peaks, label='Generated Gaussian Peaks')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.title(f'{base_name} - Gaussian Peaks')
        plt.show()