import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from readdata import *
from glob import glob
from tqdm import tqdm

if __name__ == "__main__":
    CMAP = "cool"
    ROOTDIR = ".."
    SORDIRS = sorted(list(set(glob(f"{ROOTDIR}/sor_omega*")) - set(glob(f"{ROOTDIR}/*.sh"))))

    omega_list = []
    psolve_list = []
    tsolve_list = []
    tmax_list = []

    for dir in SORDIRS:

        one_minus_omega = dir.split("_")[-1].split(".")[1]
        omega = 1.0 + float(one_minus_omega)/10.0

        om_fnames = sorted(glob(f"{dir}/*.h5"))

        om_ptimes = []
        om_itimes = []
        for i in tqdm(range(len(om_fnames)), desc=f"{dir}"):
            fname = om_fnames[i]
            ds = load(fname)
            om_ptimes.append(ds.data["TimeInfo"][-1])
            om_itimes.append(ds.data["TimeInfo"][1])

        om_tmax = load(om_fnames[-1]).data["TimeInfo"][0]

        omega_list.append(omega)
        psolve_list.append(om_ptimes)
        tsolve_list.append(om_itimes)
        tmax_list.append(om_tmax)

    fig = plt.figure(figsize=(8,8))
    cmap = matplotlib.cm.get_cmap(CMAP)

    for i in range(len(omega_list)):
        color = cmap(i/(len(omega_list)-1))
        plt.semilogy(
            np.linspace(0,tmax_list[i],len(tsolve_list[i])), 
            psolve_list[i], color = color, linewidth = 3, 
            label = r"$\omega = $" + f"{omega_list[i]}"
        )

    plt.xlabel("Simulation Time [s]",fontsize=16)
    plt.ylabel("Pressure Solve Time [ms]",fontsize=16)
    plt.grid()
    plt.legend(fontsize=12,loc="upper right")
    plt.title(r"Tolerance = $1\times10^{-7}$",fontsize=16)

    plt.savefig("varying_omega.png")