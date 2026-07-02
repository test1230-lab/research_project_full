import sys
import numpy as np
import matplotlib.pyplot as plt

if sys.argc != 3:
	print("incorrect arg count, the arg is the csv output from the mc sim.")
	sys.exit()

filename = sys.argv[2]

with open(f"output/{filename}_log_dist.csv") as fh:
    header = fh.readline().split(",")
	hist = np.loadtxt(fh, delimiter=",")

dlvpar, dlvper = float(header[2]), float(header[3])          # cm/s per bin

npar, nper = hist.shape
perp_center = (np.arange(nper) + 0.5)                  # bin centres in units of dlvper
f = hist / perp_center[None, :]                        # remove 2*pi*v_perp Jacobian -> f(v_par,v_perp)

def mirror(a, axis):                                   # reflect, dropping the duplicated zero bin
    flipped = np.flip(a, axis=axis)
    sl = [slice(None)] * a.ndim
    sl[axis] = slice(0, -1)
    return np.concatenate([flipped[tuple(sl)], a], axis=axis)

#check units!
full = mirror(mirror(f, axis=0), axis=1)
vpar = npar * dlvpar / 1e5                             # cm/s -> km/s, outer edge
vper = nper * dlvper / 1e5

plt.figure(figsize=(7, 6))
im = plt.imshow(full, origin="lower", aspect="equal", cmap="viridis",
                extent=[-vper, vper, -vpar, vpar])
plt.colorbar(im, label="f (linear)")

plt.xlabel(r"$v_\perp$ [km/s]")
plt.ylabel(r"$v_\parallel$ [km/s]")
plt.title("Velocity Distribution Function")
plt.tight_layout()
plt.savefig(f"./plots/{name}_hist_plot.png", dpi=400)
print(f"wrote: ./plots/{name}_hist_plot.png")