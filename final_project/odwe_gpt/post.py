# post.py
import numpy as np

def save_snapshot_1d(x, rho, u, p, it):
    """
    Placeholder: user can add matplotlib plotting here later.
    For now we just print min/max to check sanity.
    """
    print(f"[snapshot {it}] rho in [{rho.min():.3f}, {rho.max():.3f}], "
          f"p in [{p.min():.3f}, {p.max():.3f}]")
