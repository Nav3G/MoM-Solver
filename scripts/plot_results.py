import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Paths
BUILD = Path("output")
OUT = Path("output")
OUT.mkdir(exist_ok=True)

# Load
df = pd.read_csv(BUILD / "I_v2.csv")

z = df["z"] / (3e8 / 300e6)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(z, df["Re_I"] * 1000, label=r"$\mathrm{Re}\,I(z)$")
ax.plot(z, df["Im_I"] * 1000, label=r"$\mathrm{Im}\,I(z)$")
ax.plot(z, df["abs_I"] * 1000, label=r"$|I|$", linestyle="--", color="k")
ax.set_xlabel(r"$z / \lambda$")
ax.set_ylabel(r"$I(z)$ [mA]")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig(OUT / "I_v2.png", dpi=200, bbox_inches="tight")
plt.close(fig)

dz = pd.read_csv(BUILD / "Z_v2.csv")

Re_Z = dz.pivot(index="m", columns="n", values="Re_Z")
Im_Z = dz.pivot(index="m", columns="n", values="Im_Z")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

im1 = ax1.imshow(Re_Z.values, origin="lower", aspect="auto", cmap="RdBu_r")
ax1.set_title(r"$\mathrm{Re}\,Z_{mn}$")
ax1.set_xlabel(r"$n$")
ax1.set_ylabel(r"$m$")
fig.colorbar(im1, ax=ax1)

im2 = ax2.imshow(Im_Z.values, origin="lower", aspect="auto", cmap="RdBu_r")
ax2.set_title(r"$\mathrm{Im}\,Z_{mn}$")
ax2.set_xlabel(r"$n$")
ax2.set_ylabel(r"$m$")
fig.colorbar(im2, ax=ax2)

plt.savefig(OUT / "Z_v2.png", dpi=200, bbox_inches="tight")
plt.close(fig)