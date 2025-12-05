import numpy as np

class GridGenerator:
    """
    2D axisymmetric grid for:
      chamber -> converging -> throat -> diverging nozzle -> plume region.

    Provides:
      self.ni, self.nj
      self.length
      self.X, self.Y
    """

    def __init__(self,
                 ni=161,
                 nj=61,
                 chamber_length=0.40,
                 converge_length=0.10,
                 diverge_length=0.50,
                 plume_length=0.60,
                 r_inlet=0.40,
                 r_throat=0.15,
                 r_exit=0.35,
                 r_plume=0.50):

        # Grid sizes
        self.ni = ni
        self.nj = nj

        # Geometry parameters
        self.chamber_length  = chamber_length
        self.converge_length = converge_length
        self.diverge_length  = diverge_length
        self.plume_length    = plume_length

        self.r_inlet  = r_inlet
        self.r_throat = r_throat
        self.r_exit   = r_exit
        self.r_plume  = r_plume

        # Total length
        self.length = (
            self.chamber_length
            + self.converge_length
            + self.diverge_length
            + self.plume_length
        )

        # Allocate arrays (filled in generate())
        self.X = np.zeros((self.ni, self.nj))
        self.Y = np.zeros((self.ni, self.nj))

    # -------------------------------------------------------------
    # Radius distribution R(x)
    # -------------------------------------------------------------
    def _radius_distribution(self, x):
        """
        Piecewise radius function:
          0 .. Lch          : chamber (r_inlet)
          Lch .. Lch+Lconv  : converge to r_throat
          ... +Ldiv         : diverge to r_exit
          rest              : plume (r_plume)
        """
        Lch  = self.chamber_length
        Lcv  = self.converge_length
        Ldv  = self.diverge_length

        if x <= Lch:
            # chamber
            return self.r_inlet

        elif x <= Lch + Lcv:
            # converging
            s = (x - Lch) / Lcv
            return (1.0 - s) * self.r_inlet + s * self.r_throat

        elif x <= Lch + Lcv + Ldv:
            # diverging
            s = (x - (Lch + Lcv)) / Ldv
            return (1.0 - s) * self.r_throat + s * self.r_exit

        else:
            # plume region
            return self.r_plume

    # -------------------------------------------------------------
    # Build full 2D grid
    # -------------------------------------------------------------
    def generate(self):
        """
        Fill self.X, self.Y with coordinates for a structured (x,r) grid.
        """
        x_nodes = np.linspace(0.0, self.length, self.ni)

        for i in range(self.ni):
            x_loc = x_nodes[i]
            r_wall = self._radius_distribution(x_loc)

            for j in range(self.nj):
                eta = j / (self.nj - 1)      # 0 at centerline, 1 at wall/plume
                self.X[i, j] = x_loc
                self.Y[i, j] = eta * r_wall
        
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
