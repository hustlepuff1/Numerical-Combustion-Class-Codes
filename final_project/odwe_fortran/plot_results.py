import os
import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4
Rgas  = 287.0

# Match Fortran grid exactly
ni = 241
nj = 81

def load_snapshot_ascii(fname_txt, ni, nj):
    if not os.path.exists(fname_txt):
        alt = os.path.join("build", fname_txt)
        if os.path.exists(alt):
            fname_txt = alt
        else:
            raise FileNotFoundError(fname_txt)

    # Columns: i, j, U1, U2, U3, U4
    data = np.loadtxt(fname_txt, usecols=(2,3,4,5))

    expected_rows = ni * nj
    if data.shape[0] != expected_rows:
        raise ValueError(f"Snapshot has {data.shape[0]} rows, expected {expected_rows}.")

    # Shape: (nj, ni, 4) then transpose to (4, ni, nj)
    U = data.reshape((nj, ni, 4))
    U = np.transpose(U, (2,1,0))
    return U

def cons_to_prim(U, gamma):
    rho = U[0]
    u = U[1] / rho
    v = U[2] / rho
    E = U[3] / rho
    p = (gamma - 1.0) * rho * (E - 0.5*(u*u + v*v))
    return rho, u, v, p

if __name__ == "__main__":
    fname = "snap_00018.txt"

    U = load_snapshot_ascii(fname, ni, nj)
    rho, u, v, p = cons_to_prim(U, gamma)

    x = np.linspace(0.0, 1.0, ni)
    y = np.linspace(0.0, 0.2, nj)
    X, Y = np.meshgrid(x, y, indexing="ij")

    a = np.sqrt(gamma * p / rho)
    M = np.sqrt(u*u + v*v) / a

    plt.figure(figsize=(8,3))
    cp = plt.contourf(X, Y, M, 40, cmap="jet")
    plt.colorbar(cp)
    plt.title("Fortran ODWE - Mach Contour")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("mach_snapshot.png", dpi=200)
    print("Saved: mach_snapshot.png")
