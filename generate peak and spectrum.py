import pymatgen as mg
import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# Load CIF file
cif_file = "cif_data/Mg(Fe2O4).cif"
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


# Define the Gaussian function for fitting
def gaussian(x, amplitude, center, width):
    return amplitude * np.exp(-(x - center) ** 2 / (2 * width ** 2))

# Initialize the parameters for the Gaussian peaks
gaussian_params = []

# Fixed FWHM value for all peaks
fwhm = 0.015

# Iterate over the x and y data
for x_val, y_val in zip(x_data, y_data):
    # Check if the y value is over 1
    if y_val > 0.1:
        # Define the initial guesses for the Gaussian parameters
        amplitude_guess = y_val
        center_guess = x_val

        # Perform the curve fitting for the Gaussian peak
        popt, _ = curve_fit(gaussian, x_data, y_data, p0=[amplitude_guess, center_guess, fwhm])

        # Append the fitted parameters to the Gaussian parameters list
        gaussian_params.append(popt)

# Generate the Gaussian peaks using the fitted parameters
gaussian_peaks = []
for params in gaussian_params:
    gaussian_peak = gaussian(x_range, *params)
    gaussian_peaks.append(gaussian_peak)


# Accumulate the Gaussian peaks
accumulated_gaussian_peaks = np.sum(gaussian_peaks, axis=0)

# Plot the generated Gaussian peaks
import matplotlib.pyplot as plt

# Plot the XRD pattern and the Gaussian peaks together
plt.plot(x_range, y_values, label='XRD pattern')
plt.plot(x_range, accumulated_gaussian_peaks, label='Gaussian Peaks')
plt.xlabel('2θ (degrees)')
plt.ylabel('Intensity')
plt.title('XRD Pattern with Gaussian Peaks')
plt.legend()
plt.show()