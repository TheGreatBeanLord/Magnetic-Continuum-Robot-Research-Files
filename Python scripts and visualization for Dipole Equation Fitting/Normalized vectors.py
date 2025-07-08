import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

def load_data(filename):
    # Load data, skipping the first two header lines
    df = pd.read_csv(filename, delim_whitespace=True, skiprows=2, header=None)
    
    # Assign meaningful column names
    df.columns = ["X", "Y", "Z", "Bx", "By", "Bz"]
    
    # Filter data where BOTH X and Z are < 0.5
    df = df[~((df["X"] < 0.05) & (df["Z"] < 0.05))]
    
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

def  visualize_dipole_field(df, m_fitted):
    """ Visualizes the simulated and fitted dipole field. """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Filter data where BOTH X and Z are < 0.5
    df_filtered = df[~((df["X"] < 0.005) & (df["Z"] < 0.005))]
    
    # Extract X, Z, and field components (only in X-Z plane for 2D visualization)
    x, z = df_filtered["X"].values, df_filtered["Z"].values
    bx, bz = df_filtered["Bx"].values, df_filtered["Bz"].values
    
    # Normalize vectors for visualization (set all vectors to the same size, e.g., length 1)
    max_length = 0.01  # Set the fixed vector length you want for visualization
    bx_norm, bz_norm = bx / np.sqrt(bx**2 + bz**2) * max_length, bz / np.sqrt(bx**2 + bz**2) * max_length
    
    # Plot the original field vectors
    ax.quiver(x, z, bx_norm, bz_norm, color='r', alpha=0.6, label='Simulated Field', angles='xy', scale_units='xy', scale=1)
    
    # Compute theoretical dipole field at the same points
    mu0 = 4 * np.pi * 1e-7
    Bx_theory, Bz_theory = [], []
    
    for xi, zi in zip(x, z):
        r = np.sqrt(xi**2 + zi**2)
        if r == 0:
            Bx_theory.append(0)
            Bz_theory.append(0)
            continue
        r_hat = np.array([xi, 0, zi]) / r
        m_dot_r = np.dot(m_fitted, r_hat)
        B = (mu0 / (4 * np.pi * r**3)) * (3 * m_dot_r * r_hat - m_fitted)
        Bx_theory.append(B[0])
        Bz_theory.append(B[2])
    
    Bx_theory, Bz_theory = np.array(Bx_theory), np.array(Bz_theory)
    
    # Normalize the dipole field vectors to the same length
    Bx_theory_norm, Bz_theory_norm = Bx_theory / np.sqrt(Bx_theory**2 + Bz_theory**2) * max_length, Bz_theory / np.sqrt(Bx_theory**2 + Bz_theory**2) * max_length
    
    # Plot the fitted dipole field vectors (with the same vector size)
    ax.quiver(x, z, Bx_theory_norm, Bz_theory_norm, color='b', alpha=0.6, label='Fitted Dipole Field', angles='xy', scale_units='xy', scale=1)
    
    # Labels and legend
    ax.set_xlabel("X Position (m)")
    ax.set_ylabel("Z Position (m)")
    ax.set_title("Dipole Field Visualization")
    ax.legend()
    plt.show()
    plt.savefig('dipole_field_plot.png')  # Save the figure
    plt.close()  # Close the figure to avoid multiple open windows

# Example usage:
df = load_data("MagData.fld")
m_fitted = fit_dipole_moment(df)
visualize_dipole_field(df, m_fitted)
