#!/usr/bin/env python3
"""Plot the Monte Carlo ion velocity distribution histogram.

The simulation only stores the first quadrant (|v_par|, |v_perp|), so we
reflect it across both axes to build the full v_par-v_perp plane.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

name = sys.argv[1] if len(sys.argv) > 1 else "E100"
hist = np.load(f"output/{name}hist.npy")          # shape (npar, nper), [v_par][v_perp]

# bin width (m/s) from the csv header (dlvpar, dlvper)
with open(f"output/{name}_log_dist.csv") as f:
    hdr = f.readline().split(",")
dlvpar, dlvper = float(hdr[2]), float(hdr[3])


npar, nper = hist.shape
# velocity extent in km/s (bin edges), centred on zero
vpar = (npar - 0.5) * dlvpar / 1000.0
vper = (nper - 0.5) * dlvper / 1000.0

# x = perpendicular velocity, y = parallel velocity
plt.figure(figsize=(7, 6))
im = plt.imshow(hist, origin="lower", aspect="auto", cmap="viridis",
                extent=[-vper, vper, -vpar, vpar])
plt.colorbar(im, label="counts")
plt.xlabel(r"$v_\perp$ [km/s]")
plt.ylabel(r"$v_\parallel$ [km/s]")
plt.title(f"Ion velocity distribution ({name})")
plt.tight_layout()
out = f"output/{name}_hist.png"
plt.savefig(out, dpi=150)
plt.show()
print("wrote", out)
