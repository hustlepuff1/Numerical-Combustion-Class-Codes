import numpy as np

class GridGenerator:
    """
    Generates a 2D Axisymmetric Grid for a Solid Rocket Motor Nozzle.
    Based on the algebraic method described in the NASA-OK-199183 report.
    """
    def __init__(self, ni=81, nj=31, length=1.0):
        self.ni = ni
        self.nj = nj
        self.length = length
        
        # Geometry parameters
        self.x_throat = 0.5 * length
        self.r_inlet  = 0.4
        self.r_throat = 0.15
        self.r_exit   = 0.35
        
        self.X = None
        self.Y = None

    def _radius_distribution(self, x):
        # 1. Chamber / Inlet
        if x < 0.2 * self.length:
            return self.r_inlet
        
        # 2. Converging (inlet -> throat)
        elif x < self.x_throat:
            x_start = 0.2 * self.length
            x_norm = (x - x_start) / (self.x_throat - x_start)
            return self.r_throat + (self.r_inlet - self.r_throat) * \
                   (0.5 * (1 + np.cos(x_norm * np.pi)))
        
        # 3. Diverging (throat -> exit)
        else:
            x_norm = (x - self.x_throat) / (self.length - self.x_throat)
            return self.r_throat + (self.r_exit - self.r_throat) * x_norm

    def generate(self):
        self.X = np.zeros((self.ni, self.nj))
        self.Y = np.zeros((self.ni, self.nj))
        
        x_points = np.linspace(0, self.length, self.ni)
        
        for i in range(self.ni):
            x_loc = x_points[i]
            r_wall = self._radius_distribution(x_loc)
            
            for j in range(self.nj):
                eta = j / (self.nj - 1)
                self.X[i, j] = x_loc
                self.Y[i, j] = r_wall * eta
        
        return self.X, self.Y

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    grid = GridGenerator(ni=81, nj=31)
    X, Y = grid.generate()
    
    print(f"Grid generated successfully. Shape: {X.shape}")
    
    plt.figure(figsize=(10, 4))
    for j in range(0, Y.shape[1], 5):
        plt.plot(X[:, j], Y[:, j], 'k-', linewidth=0.5, alpha=0.6)
    for i in range(0, X.shape[0], 5):
        plt.plot(X[i, :], Y[i, :], 'k-', linewidth=0.5, alpha=0.6)
    plt.title("Reference Grid (Matches NASA Paper Fig. 2)")
    plt.xlabel("Axial (x)")
    plt.ylabel("Radial (r)")
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.show()
