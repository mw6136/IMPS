import numpy as np
import matplotlib.pyplot as plt
from readdata import *
from glob import glob
from tqdm import tqdm

if __name__ == "__main__":

    RUNTYPE = "hightol_64x64"
    JAC_DIR = f"../{RUNTYPE}_jac"
    GS_DIR = f"../{RUNTYPE}_gs"
    SOR_DIR = f"../{RUNTYPE}_sor"

    jac_fnames = sorted(glob(f"{JAC_DIR}/*.h5"))
    gs_fnames = sorted(glob(f"{GS_DIR}/*.h5"))
    sor_fnames = sorted(glob(f"{SOR_DIR}/*.h5"))

    gs_ptimes = []
    gs_itimes = []
    for i in tqdm(range(len(gs_fnames))):
        fname = gs_fnames[i]
        ds = load(fname)
        gs_ptimes.append(ds.data["TimeInfo"][-1])
        gs_itimes.append(ds.data["TimeInfo"][1])

    gs_tmax = load(gs_fnames[-1]).data["TimeInfo"][0]

    sor_ptimes = []
    sor_itimes = []
    for i in tqdm(range(len(sor_fnames))):
        fname = sor_fnames[i]
        ds = load(fname)
        sor_ptimes.append(ds.data["TimeInfo"][-1])
        sor_itimes.append(ds.data["TimeInfo"][1])

    sor_tmax = load(sor_fnames[-1]).data["TimeInfo"][0]

    jac_ptimes = []
    jac_itimes = []
    for i in tqdm(range(len(jac_fnames))):
        fname = jac_fnames[i]
        ds = load(fname)
        jac_ptimes.append(ds.data["TimeInfo"][-1])
        jac_itimes.append(ds.data["TimeInfo"][1])

    jac_tmax = load(jac_fnames[-1]).data["TimeInfo"][0]

    fig = plt.figure(figsize=(8,8))
    
    plt.semilogy(np.linspace(0,jac_tmax,len(jac_ptimes)), jac_ptimes, color="navy", linewidth=3, label="Jacobi")
    plt.semilogy(np.linspace(0,gs_tmax,len(gs_ptimes)), gs_ptimes, color="firebrick", linewidth=3, label="Gauss-Seidel")
    plt.semilogy(np.linspace(0,sor_tmax,len(sor_ptimes)), sor_ptimes, color="goldenrod", linewidth=3, label="SOR")

    # plt.semilogy(np.linspace(0,100,len(jac_ptimes))[1:], np.array(jac_itimes)[1:], color="navy", linewidth=3, linestyle=":", label="Total Jacobi")
    # plt.semilogy(np.linspace(0,100,len(gs_ptimes))[1:], np.array(gs_itimes)[1:], color="firebrick", linewidth=3, linestyle=":", label="Total Gauss-Seidel")
    # plt.semilogy(np.linspace(0,100,len(sor_ptimes))[1:],  np.array(sor_itimes)[1:], color="goldenrod", linestyle=":", linewidth=3, label="Total SOR")

    plt.xlabel("Simulation Time [s]",fontsize=16)
    plt.ylabel("Pressure Solve Time [ms]",fontsize=16)
    plt.grid()
    plt.legend(fontsize=14)

    if RUNTYPE == "lowtol_64x64":
        plt.title(r"Tolerance = $1\times10^{-7}$",fontsize=16)
        plt.savefig("lowtol_timetest.png")
    else:
        plt.title(r"Tolerance = $1\times10^{-5}$",fontsize=16)
        plt.savefig("hightol_timetest.png")

