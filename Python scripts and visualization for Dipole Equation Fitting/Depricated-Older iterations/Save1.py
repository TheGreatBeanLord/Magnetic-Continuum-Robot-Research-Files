import numpy as np
import pandas as pd
from scipy.optimize import least_squares

def load_data(filename):
    # Load data, skipping the first two header lines
    df = pd.read_csv(filename, delim_whitespace=True, skiprows=2, header=None)
    
    # Assign meaningful column names
    df.columns = ["X", "Y", "Z", "Bx", "By", "Bz"]
    
    # Filter data where BOTH X and Z are < 0.5
    df = df[~((df["X"] < 0.005) & (df["Z"] < 0.005))]
    
    return df

def dipole_field(m, x, y, z):
    """ Compute the theoretical dipole field at (x, y, z) given dipole moment m. """
    mu0 = 4 * np.pi * 1e-7  # Permeability of free space
    r = np.sqrt(x**2 + y**2 + z**2)
    if r == 0:
        return np.array([0.0, 0.0, 0.0])  # Avoid division by zero
    r_hat = np.array([x, y, z]) / r
    m_dot_r = np.dot(m, r_hat)
    B = (mu0 / (4 * np.pi * r**3)) * (3 * m_dot_r * r_hat - m)
    return B

def residuals(m, x, y, z, bx, by, bz):
    """ Compute residuals between simulated and theoretical dipole fields. """
    B_theory = np.array([dipole_field(m, xi, yi, zi) for xi, yi, zi in zip(x, y, z)])
    if B_theory.shape != (len(bx), 3):
        return np.zeros(len(bx) * 3)  # Ensure proper shape in case of empty data
    return np.ravel(B_theory - np.column_stack((bx, by, bz)))

def fit_dipole_moment(df):
    """ Fit the dipole moment using least squares. """
    if df.empty:
        print("No valid data points remaining after filtering.")
        return np.array([0.0, 0.0, 0.0])
    
    x, y, z = df["X"].values, df["Y"].values, df["Z"].values
    bx, by, bz = df["Bx"].values, df["By"].values, df["Bz"].values
    
    # Initial guess for m
    m0 = np.array([1e-3, 1e-3, 1e-3])
    
    # Perform least squares fitting
    result = least_squares(residuals, m0, args=(x, y, z, bx, by, bz))
    
    return result.x

# Example usage:
filename = "MagData.fld"
df = load_data(filename)
m_fitted = fit_dipole_moment(df)
print("Fitted dipole moment:", m_fitted)
