# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last updated: Jul 24, 2024

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import multiprocessing as mp
from functools import partial
import argparse

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

def gettingFacets(filename,includeCoat='true'):
    exe = ["./getFacet-threePhase", filename, includeCoat]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2, z2)))
                    segs.append(((-r1, z1),(-r2, z2)))
                    skip = True
    return segs

def gettingfield(filename, zmin, zmax, rmax, nr):
    exe = ["./getData-threePhase", filename, str(zmin), str(0), str(zmax), str(rmax), str(nr)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    # print(temp2) #debugging
    Rtemp, Ztemp, D2temp, veltemp, Utemp, Vtemp = [],[],[],[],[],[]

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            D2temp.append(float(temp3[2]))
            veltemp.append(float(temp3[3]))

    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    D2 = np.asarray(D2temp)
    vel = np.asarray(veltemp)
    nz = int(len(Z)/nr)

    # print("nr is %d %d" % (nr, len(R))) # debugging
    print("nz is %d" % nz)

    R.resize((nz, nr))
    Z.resize((nz, nr))
    D2.resize((nz, nr))
    vel.resize((nz, nr))

    return R, Z, D2, vel, nz
# ----------------------------------------------------------------------------------------------------------------------

def process_timestep(ti, folder, nGFS, GridsPerR, rmin, rmax, zmin, zmax, lw):
    t = 0.01 * ti
    place = f"intermediate/snapshot-{t:.4f}"
    name = f"{folder}/{int(t*1000):08d}.png"

    if not os.path.exists(place):
        print(f"{place} File not found!")
        return

    if os.path.exists(name):
        print(f"{name} Image present!")
        return

    segs1 = gettingFacets(place)
    segs2 = gettingFacets(place, 'false')

    if not segs1 and not segs2:
        print(f"Problem in the available file {place}")
        return

    nr = int(GridsPerR * rmax)
    R, Z, taus, vel, nz = gettingfield(place, zmin, zmax, rmax, nr)
    zminp, zmaxp, rminp, rmaxp = Z.min(), Z.max(), R.min(), R.max()

    # Plotting
    AxesLabel, TickLabel = 50, 20
    fig, ax = plt.subplots()
    fig.set_size_inches(19.20, 10.80)

    ax.plot([0, 0], [zmin, zmax], '-.', color='grey', linewidth=lw)
    ax.plot([rmin, rmin], [zmin, zmax], '-', color='black', linewidth=lw)
    ax.plot([rmin, rmax], [zmin, zmin], '-', color='black', linewidth=lw)
    ax.plot([rmin, rmax], [zmax, zmax], '-', color='black', linewidth=lw)
    ax.plot([rmax, rmax], [zmin, zmax], '-', color='black', linewidth=lw)

    line_segments = LineCollection(segs2, linewidths=4, colors='green', linestyle='solid')
    ax.add_collection(line_segments)
    line_segments = LineCollection(segs1, linewidths=4, colors='blue', linestyle='solid')
    ax.add_collection(line_segments)

    cntrl1 = ax.imshow(taus, cmap="hot_r", interpolation='Bilinear', origin='lower', extent=[-rminp, -rmaxp, zminp, zmaxp], vmax=1.0, vmin=-3.0)

    # TODO: fixme the colorbar bounds for taup must be set manually based on the simulated case.
    cntrl2 = ax.imshow(vel, interpolation='Bilinear', cmap='Blues', origin='lower', extent=[rminp, rmaxp, zminp, zmaxp], vmax=0.0, vmin=1.0)

    ax.set_aspect('equal')
    ax.set_xlim(rmin, rmax)
    ax.set_ylim(zmin, zmax)
    ax.set_title(f'$t$ = {t:4.3f}', fontsize=TickLabel)

    l, b, w, h = ax.get_position().bounds
    cb1 = fig.add_axes([l+0.05*w, b-0.05, 0.40*w, 0.03])
    c1 = plt.colorbar(cntrl1, cax=cb1, orientation='horizontal')
    c1.set_label(r'$\log_{10}\left(\mathcal{D}\right)$', fontsize=TickLabel, labelpad=5)
    c1.ax.tick_params(labelsize=TickLabel)
    c1.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    cb2 = fig.add_axes([l+0.55*w, b-0.05, 0.40*w, 0.03])
    c2 = plt.colorbar(cntrl2, cax=cb2, orientation='horizontal')
    c2.ax.tick_params(labelsize=TickLabel)
    c2.set_label(r'$\|u_i\|$', fontsize=TickLabel)
    c2.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    ax.axis('off')

    plt.savefig(name, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=mp.cpu_count(), help='Number of CPUs to use')
    parser.add_argument('--nGFS', type=int, default=1550, help='Number of restart files to process')
    parser.add_argument('--ZMAX', type=float, default=2e0, help='Maximum Z value')
    parser.add_argument('--RMAX', type=float, default=4.0, help='Maximum R value')
    parser.add_argument('--ZMIN', type=float, default=-2e0, help='Minimum Z value')
    args = parser.parse_args()

    nGFS = args.nGFS
    ZMAX = args.ZMAX
    RMAX = args.RMAX
    ZMIN = args.ZMIN    
    num_processes = args.CPUs

    rmin, rmax, zmin, zmax = [-RMAX, RMAX, ZMIN, ZMAX]
    GridsPerR = 100

    lw = 2
    folder = 'Video'

    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Prepare the partial function with fixed arguments
    process_func = partial(process_timestep, folder=folder, nGFS=nGFS, 
                           GridsPerR=GridsPerR, 
                           rmin=rmin, rmax=rmax, zmin=zmin, zmax=zmax, lw=lw)

    # Use all available CPU cores
    num_processes = mp.cpu_count()
    
    # Create a pool of worker processes
    with mp.Pool(processes=num_processes) as pool:
        # Map the process_func to all timesteps
        pool.map(process_func, range(nGFS))

if __name__ == "__main__":
    main()