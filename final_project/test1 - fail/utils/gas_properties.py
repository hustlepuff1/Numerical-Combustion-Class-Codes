import numpy as np

class GasProperties:
    """
    Real Thermodynamics from NASA CEA (User Run: 64 bar, AP/Al/HTPB).
    Data fits the output for Equilibrium Flow.
    """
    def __init__(self):
        # 1. TEMPERATURE (K)
        # Extracted from your CEA Output (Low T -> High T order)
        # Exit (Ae/At=10) -> Throat -> Chamber
        self.T_data = np.array([
            2431.62, 2520.26, 2634.64, 2797.15, 3090.55, 3548.77, 3729.08
        ])
        
        # 2. GAMMA (Matched to Temperature)
        self.gamma_data = np.array([
            1.1553, 1.1503, 1.1442, 1.1370, 1.1290, 1.1250, 1.1250
        ])
        
        # 3. GAS CONSTANT R
        # R = 8314.46 / Molecular_Weight (M)
        # M from CEA: 27.95, 27.88, 27.77, 27.59, 27.18, 26.46, 26.17
        M_data = np.array([
            27.955, 27.884, 27.777, 27.592, 27.187, 26.464, 26.174
        ])
        self.R_data = 8314.46 / M_data

    def get_properties(self, T_grid):
        """
        Input: T_grid (2D numpy array of Temperatures)
        Output: gamma_grid, R_grid
        """
        # Clamp Temperature to avoid extrapolation errors
        T_safe = np.clip(T_grid, self.T_data[0], self.T_data[-1])
        
        # Linear Interpolation
        gamma = np.interp(T_safe, self.T_data, self.gamma_data)
        R = np.interp(T_safe, self.T_data, self.R_data)
        
        return gamma, R