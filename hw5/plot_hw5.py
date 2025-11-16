#==========================================================
# File: plot_hw5.py
# Plot results from Fortran solver for Homework #5
#==========================================================
import numpy as np
import matplotlib.pyplot as plt

mu_list = [1, 5, 50, 500, 1000]

for mu in mu_list:
    filename = f"solution_mu_{mu:04d}.dat"
    print(f"Loading {filename} ...")
    # Skip comment line starting with '#'
    data = np.loadtxt(filename, comments='#')

    x        = data[:, 0]
    u_exact  = data[:, 1]
    u_imp    = data[:, 2]
    u_exp1   = data[:, 3]
    u_exp2   = data[:, 4]

    plt.figure(figsize=(6, 4))
    plt.plot(x, u_exact,  label='Exact',           linewidth=2)
    plt.plot(x, u_imp,    '--', label='Implicit (1st order)')
    plt.plot(x, u_exp1,   ':',  label='Explicit Euler (1st)')
    plt.plot(x, u_exp2,   '-.', label='Explicit RK2 (2nd)')

    plt.xlabel('x')
    plt.ylabel('u(x, t=0.4)')
    plt.title(f'Model reactive flow: μ = {mu}')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'hw5_mu_{mu:04d}.png', dpi=200)

# Show all figures at the end
plt.show()
