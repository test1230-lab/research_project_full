#!/usr/bin/env python3
"""Plot the _log_dist.csv directly, with NO extra division applied.
The CSV already stores az = ln(Z/((perp+0.5)*cntnorm)), i.e. the Jacobian
division is already baked in, so the values are read and plotted as-is."""
import sys
import numpy as np
import matplotlib.pyplot as plt

name = sys.argv[1] if len(sys.argv) > 1 else "E0"
with open(f"output/{name}_log_dist.csv") as fh:
    nper, npar, dlvper, dlvpar, cntnorm = fh.readline().split(",")
    az = np.loadtxt(fh, delimiter=",")          # 24x24, az[ipar][iper], read as-is
dlvpar, dlvper = float(dlvpar), float(dlvper)

az = np.ma.masked_less_equal(az, -99.0)         # mask the empty-bin floor (-99.99)
az = az - az.max()                              # show ln(f/fmax); no /v_perp applied

def mirror(a, axis):
    flipped = np.flip(a, axis=axis)
    sl = [slice(None)] * a.ndim
    sl[axis] = slice(0, -1)
    return np.ma.concatenate([flipped[tuple(sl)], a], axis=axis)

full = mirror(mirror(az, axis=0), axis=1)
vpar = az.shape[0] * dlvpar / 1e5               # cm/s -> km/s
vper = az.shape[1] * dlvper / 1e5

plt.figure(figsize=(7, 6))
im = plt.imshow(full, origin="lower", aspect="equal", cmap="viridis",
                extent=[-vper, vper, -vpar, vpar])
plt.colorbar(im, label=r"$\ln(f/f_{max})$  (straight from CSV)")
plt.xlabel(r"$v_\perp$ [km/s]")
plt.ylabel(r"$v_\parallel$ [km/s]")
plt.title(f"{name} _log_dist.csv read directly (no division applied)")
plt.tight_layout()
out = f"output/{name}_csv_direct.png"
plt.savefig(out, dpi=150)
print("wrote", out)
print("header: nper,npar,dlvper,dlvpar,cntnorm =", nper, npar, dlvper, dlvpar, cntnorm.strip())
print("az range (excl. floor):", f"{az.min():.2f} .. {az.max():.2f}")
