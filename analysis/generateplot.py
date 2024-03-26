import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from readdata import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    "-f", "--field"
)

parser.add_argument(
    "-o", "--output"
)

args = parser.parse_args()

if __name__ == "__main__":
    field = args.field

    fname = args.output
    ds = load(fname)

    if field == "VelMag":
        u = ds.data["VelX1"]
        v = ds.data["VelX2"]
        parray = np.rot90(np.sqrt(u**2 + v**2))
    else:
        parray = np.rot90(ds.data[field])

    fig = plt.figure(figsize=(5,5))
    ax1 = fig.add_subplot(111)

    im = ax1.imshow(parray, extent=(0,1,0,1), cmap="jet")
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

    ax1.axis("off")

    plt.savefig(f"{field}_{fname.split('.')[-2]}.png")