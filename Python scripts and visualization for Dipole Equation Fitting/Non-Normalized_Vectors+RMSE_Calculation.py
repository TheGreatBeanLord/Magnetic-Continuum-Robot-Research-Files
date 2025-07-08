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

def calculate_error(sim_bx, sim_bz, fitted_bx, fitted_bz):
    """ Calculate the root mean square error (RMSE) between simulated and fitted field components. """
    # Compute the differences between simulated and fitted field components
    error_bx = sim_bx - fitted_bx
    error_bz = sim_bz - fitted_bz
    
    # Calculate the RMSE for each component
    rmse_bx = np.sqrt(np.mean(error_bx**2))
    rmse_bz = np.sqrt(np.mean(error_bz**2))
    
    # Total RMSE as the combined RMSE for both components
    total_rmse = np.sqrt(rmse_bx**2 + rmse_bz**2)
    return total_rmse, rmse_bx, rmse_bz

def visualize_dipole_field(df, m_fitted):
    """ Visualizes the simulated and fitted dipole field and calculates the error. """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Filter data where BOTH X and Z are < 0.5
    df_filtered = df[~((df["X"] < 0.005) & (df["Z"] < 0.005))]
    
    # Extract X, Z, and field components (only in X-Z plane for 2D visualization)
    x, z = df_filtered["X"].values, df_filtered["Z"].values
    bx, bz = df_filtered["Bx"].values, df_filtered["Bz"].values
    
    # Normalize vectors for visualization (preserve relative size, adjust scaling)
    max_field_strength = np.max(np.sqrt(bx**2 + bz**2))  # Find the maximum field strength
    scale_factor = 0.1 / max_field_strength  # Scale down the vectors to fit the plot better
    
    bx_norm, bz_norm = bx * scale_factor, bz * scale_factor
    
    # Plot the original field vectors with relative size
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
    
    # Normalize the dipole field vectors to preserve relative size
    max_theory_strength = np.max(np.sqrt(Bx_theory**2 + Bz_theory**2))
    theory_scale_factor = 0.1 / max_theory_strength  # Scale down the vectors for the dipole field
    
    Bx_theory_norm, Bz_theory_norm = Bx_theory * theory_scale_factor, Bz_theory * theory_scale_factor
    
    # Plot the fitted dipole field vectors with relative size
    ax.quiver(x, z, Bx_theory_norm, Bz_theory_norm, color='b', alpha=0.6, label='Fitted Dipole Field', angles='xy', scale_units='xy', scale=1)
    
    # Labels and legend
    ax.set_xlabel("X Position (m)")
    ax.set_ylabel("Z Position (m)")
    ax.set_title("Dipole Field Visualization")
    ax.legend()
    plt.show()
    
    # Calculate the error between the simulated and fitted fields
    total_rmse, rmse_bx, rmse_bz = calculate_error(bx, bz, Bx_theory, Bz_theory)
    
    # Print the error metrics
    print(f"Total RMSE (Simulated vs Fitted): {total_rmse:.4f}")
    print(f"RMSE for Bx: {rmse_bx:.4f}")
    print(f"RMSE for Bz: {rmse_bz:.4f}")
    
    return total_rmse, rmse_bx, rmse_bz

# Example usage:
df = load_data("MagData.fld")
m_fitted = fit_dipole_moment(df)
visualize_dipole_field(df, m_fitted)
