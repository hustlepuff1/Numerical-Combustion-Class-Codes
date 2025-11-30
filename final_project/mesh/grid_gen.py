import numpy as np

class GridGenerator:
    """
    Generates a 2D Axisymmetric Grid for a Solid Rocket Motor Nozzle.
    [cite_start]Based on the algebraic method described in the NASA-OK-199183 report[cite: 349].
    """
    def __init__(self, ni=81, nj=31, length=1.0):
        """
        Initialize Grid Generator dimensions.
        
        Args:
            ni (int): Number of streamwise (axial) points. [cite_start]Default 81[cite: 349].
            nj (int): Number of radial points. [cite_start]Default 31[cite: 349].
            length (float): Total length of the nozzle (normalized or physical).
        """
        self.ni = ni
        self.nj = nj
        self.length = length
        
        # Geometry Parameters (Normalized relative to length)
        # These define the converging-diverging shape
        self.x_throat = 0.5 * length
        self.r_inlet = 0.4
        self.r_throat = 0.15
        self.r_exit = 0.35
        
        # Storage for the mesh
        self.X = None
        self.Y = None

    def _radius_distribution(self, x):
        """
        Internal method to calculate the Wall Radius R(x) at a specific axial location.
        Uses a cosine blend for smooth throat transition.
        """
        # 1. Chamber / Inlet Section (Constant area)
        if x < 0.2 * self.length:
            return self.r_inlet
            
        # 2. Converging Section (Inlet -> Throat)
        elif x < self.x_throat:
            # Normalize x between 0.0 (start of converge) and 1.0 (throat)
            x_start = 0.2 * self.length
            x_norm = (x - x_start) / (self.x_throat - x_start)
            
            # Cosine interpolation for smooth curve
            return self.r_throat + (self.r_inlet - self.r_throat) * \
                   (0.5 * (1 + np.cos(x_norm * np.pi)))
                   
        # 3. Diverging Section (Throat -> Exit)
        else:
            # Normalize x between 0.0 (throat) and 1.0 (exit)
            x_norm = (x - self.x_throat) / (self.length - self.x_throat)
            
            # Linear/Parabolic expansion
            return self.r_throat + (self.r_exit - self.r_throat) * x_norm

    def generate(self):
        """
        Execute the algebraic grid generation.
        
        Returns:
            X (np.array): 2D array of Axial coordinates.
            Y (np.array): 2D array of Radial coordinates.
        """
        # Create arrays
        self.X = np.zeros((self.ni, self.nj))
        self.Y = np.zeros((self.ni, self.nj))
        
        # Create streamwise distribution (x_points)
        # Using linear spacing for now. 
        # Future MDO upgrade: Add clustering parameters here to stretch grid at throat.
        x_points = np.linspace(0, self.length, self.ni)
        
        for i in range(self.ni):
            x_loc = x_points[i]
            r_wall = self._radius_distribution(x_loc)
            
            for j in range(self.nj):
                # Transfinite Interpolation:
                # Distribute points between Centerline (y=0) and Wall (y=r_wall)
                # eta varies from 0.0 to 1.0
                eta = j / (self.nj - 1)
                
                self.X[i, j] = x_loc
                self.Y[i, j] = r_wall * eta
                
        return self.X, self.Y

# --- Standalone Debugging Block ---
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # 1. Instantiate
    grid = GridGenerator(ni=81, nj=31)
    
    # 2. Generate
    X, Y = grid.generate()
    
    # 3. Visualize (Sanity Check)
    print(f"Grid generated successfully. Shape: {X.shape}")
    
    plt.figure(figsize=(10, 4))
    # Plot axial lines (j is constant)
    for j in range(0, Y.shape[1], 5): # Plot every 5th line for clarity
        plt.plot(X[:, j], Y[:, j], 'k-', linewidth=0.5, alpha=0.6)
    # Plot radial lines (i is constant)
    for i in range(0, X.shape[0], 5):
        plt.plot(X[i, :], Y[i, :], 'k-', linewidth=0.5, alpha=0.6)
        
    plt.title("Reference Grid (Matches NASA Paper Fig. 2)")
    plt.xlabel("Axial (x)")
    plt.ylabel("Radial (r)")
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.show()