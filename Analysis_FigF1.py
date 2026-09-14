import numpy as np
import glob
import h5py 
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.integrate import solve_ivp
#import cmasher as cmr
from matplotlib.lines import Line2D
from matplotlib import colors
from matplotlib.colors import Normalize


import sys    # needed for exit codes
import os      # file specific calls
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

import matplotlib.colors as mcolors
import numpy as np
plt.rcParams.update({
    'font.size': 14,               # base font size
    'axes.titlesize': 16,          # title size
    'axes.labelsize': 14,          # x/y label size
    'xtick.labelsize': 12,         # x tick labels
    'ytick.labelsize': 12,         # y tick labels
    'legend.fontsize': 12,         # legend text
})



def load_stuff(string):
    v = np.load("Simulations/"+string+"/output/velocities.npy")
    n = np.load("Simulations/"+string+"/output/densities.npy")
    M = np.load("Simulations/"+string+"/output/Mdots.npy")
    T = np.load("Simulations/"+string+"/output/Temperatures.npy")
    sg = np.load("Simulations/"+string+"/output/subgrid.npy")
    return v, n, M, T, sg

fid_v, fid_n, fid_M, fid_T, fid_sg = load_stuff("fid")

#one color 
fcl_v, fcl_n, fcl_M, fcl_T, fcl_sg = load_stuff("fcld_l") #0.1
fch_v, fch_n, fch_M,fch_T, fch_sg = load_stuff("fcld_h") #0.9

#one color 
#NoMix_v, NoMix_n, NoMix_M, NoMix_sg = load_stuff("NoMix")
mixl_v, mixl_n, mixl_M,mixl_T, mixl_sg = load_stuff("fc0.03") 
mixh_v, mixh_n, mixh_M,mixh_T, mixh_sg = load_stuff("fc0.3")

#one color 
rcl001_v, rcl001_n, rcl001_M, rcl001_T,rcl001_sg = load_stuff("rcl0.01") 
#rcl002_v, rcl002_n, rcl002_M, rcl002_sg = load_stuff("rcl0.02") 
rcl06_v, rcl06_n, rcl06_M, rcl06_T,rcl06_sg = load_stuff("rcl0.6") 
#rcl02_v, rcl02_n, rcl02_M, rcl02_sg = load_stuff("rcl0.1") 



col = [
    'mediumturquoise',
    '#FF7F50',
    'sandybrown'
]
ls = [
    '-',
    '--',
    '-.',
    ':',
]

lab = ["$f_{cold}$ = 0.1", "$f_{cold}$ = 0.9", "$f_{mix} = 0.03$", "$f_{cool} = 0.3$", "$r_{cl}$ = 0.01 kpc", "$r_{cl}$ = 0.6 kpc"]



#### single plot con le mcrit 
fig, ax = plt.subplots()
plt.plot(fid_sg[0,:], fid_sg[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(fcl_sg[0,:], fcl_sg[2,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(fch_sg[0,:], fch_sg[2,:], linestyle = ls[1], color = col[0], label = lab[1])

#plt.plot(NoMix_sg[0,:], NoMix_sg[2,:], linestyle = ls[0], color = col[1], label = lab[2])
plt.plot(mixl_sg[0,:], mixl_sg[2,:], linestyle = ls[1], color = col[1], label = lab[2])
plt.plot(mixh_sg[0,:], mixh_sg[2,:], linestyle = ls[2], color = col[1], label = lab[3])

plt.plot(rcl001_sg[0,:], rcl001_sg[2,:], linestyle = ls[0], color = col[2], label = lab[4])
plt.plot(rcl06_sg[0,:], rcl06_sg[2,:], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[2,:], linestyle = ls[3], color = col[2], label = lab[8])

plt.legend()
plt.ylim(0.1,1e4)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$m_{crit}$  [$M_\odot$]")
plt.savefig("mcl_impo")

Rmin=0.5
fig, ax = plt.subplots(figsize=(6, 5),  dpi=300)
plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], fid_sg[1,fid_sg[0,:]>Rmin], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(fcl_sg[0,fid_sg[0,:]>Rmin], fcl_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(fch_sg[0,fid_sg[0,:]>Rmin], fch_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[1], color = col[0], label = lab[1])

#plt.plot(mixl_sg[0,:], mixl_sg[1,:], linestyle = ls[1], color = col[1], label = lab[2])
#plt.plot(mixh_sg[0,:], mixh_sg[1,:], linestyle = ls[2], color = col[1], label = lab[3])

plt.plot(rcl001_sg[0,fid_sg[0,:]>Rmin], rcl001_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[0], color = col[2], label = lab[4])
plt.plot(rcl06_sg[0,fid_sg[0,:]>Rmin], rcl06_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[1,:], linestyle = ls[3], color = col[2], label = lab[8])

plt.legend(
    prop={'family': 'Times New Roman', 'size': 10},
)
plt.ylim(1e-4,1)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$r$  [kpc]")
plt.savefig("rcl_impo")



fig, ax = plt.subplots()
plt.plot(fid_sg[0,:], fid_sg[5,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(fcl_sg[0,:], fcl_sg[5,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(fch_sg[0,:], fch_sg[5,:], linestyle = ls[1], color = col[0], label = lab[1])

#plt.plot(NoMix_sg[0,:], NoMix_sg[2,:], linestyle = ls[0], color = col[1], label = lab[2])
plt.plot(mixl_sg[0,:], mixl_sg[5,:], linestyle = ls[1], color = col[1], label = lab[2])
plt.plot(mixh_sg[0,:], mixh_sg[5,:], linestyle = ls[2], color = col[1], label = lab[3])

plt.plot(rcl001_sg[0,:], rcl001_sg[5,:], linestyle = ls[0], color = col[2], label = lab[4])
plt.plot(rcl06_sg[0,:], rcl06_sg[5,:], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[2,:], linestyle = ls[3], color = col[2], label = lab[8])

plt.legend()
#plt.ylim(0.1,1e4)
#plt.yscale('log')
plt.yscale('symlog')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$udotmix$  []")
plt.savefig("udotmix")



fig, ax = plt.subplots()
plt.plot(fid_sg[0,:], fid_sg[6,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(fcl_sg[0,:], fcl_sg[6,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(fch_sg[0,:], fch_sg[6,:], linestyle = ls[1], color = col[0], label = lab[1])

#plt.plot(NoMix_sg[0,:], NoMix_sg[2,:], linestyle = ls[0], color = col[1], label = lab[2])
plt.plot(mixl_sg[0,:], mixl_sg[6,:], linestyle = ls[1], color = col[1], label = lab[2])
plt.plot(mixh_sg[0,:], mixh_sg[6,:], linestyle = ls[2], color = col[1], label = lab[3])

plt.plot(rcl001_sg[0,:], rcl001_sg[6,:], linestyle = ls[0], color = col[2], label = lab[4])
plt.plot(rcl06_sg[0,:], rcl06_sg[6,:], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[2,:], linestyle = ls[3], color = col[2], label = lab[8])

plt.legend()
#plt.ylim(0.1,1e4)
plt.yscale('symlog')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$udotdrarg$  []")
plt.savefig("udotdrag")



fig, ax = plt.subplots()
plt.plot(fid_sg[0,:], fid_sg[7,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(fcl_sg[0,:], fcl_sg[7,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(fch_sg[0,:], fch_sg[7,:], linestyle = ls[1], color = col[0], label = lab[1])

#plt.plot(NoMix_sg[0,:], NoMix_sg[2,:], linestyle = ls[0], color = col[1], label = lab[2])
plt.plot(mixl_sg[0,:], mixl_sg[7,:], linestyle = ls[1], color = col[1], label = lab[2])
plt.plot(mixh_sg[0,:], mixh_sg[7,:], linestyle = ls[2], color = col[1], label = lab[3])

plt.plot(rcl001_sg[0,:], rcl001_sg[7,:], linestyle = ls[0], color = col[2], label = lab[4])
plt.plot(rcl06_sg[0,:], rcl06_sg[7,:], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[2,:], linestyle = ls[3], color = col[2], label = lab[8])

plt.legend()
#plt.ylim(0.1,1e4)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$tcool$  []")
plt.savefig("tcool")




fig, ax = plt.subplots()
plt.plot(fid_T[0,:], fid_T[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(fcl_T[0,:], fcl_T[1,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(fch_T[0,:], fch_T[1,:], linestyle = ls[1], color = col[0], label = lab[1])

#plt.plot(NoMix_sg[0,:], NoMix_sg[2,:], linestyle = ls[0], color = col[1], label = lab[2])
plt.plot(mixl_T[0,:], mixl_T[1,:], linestyle = ls[1], color = col[1], label = lab[2])
plt.plot(mixh_T[0,:], mixh_T[1,:], linestyle = ls[2], color = col[1], label = lab[3])

plt.plot(rcl001_T[0,:], rcl001_T[1,:], linestyle = ls[0], color = col[2], label = lab[4])
plt.plot(rcl06_T[0,:], rcl06_T[1,:], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[2,:], linestyle = ls[3], color = col[2], label = lab[8])

plt.legend()
#plt.ylim(0.1,1e4)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$T$  []")
plt.savefig("T")
#qui

#### cfr Velocities
    
    ## - Hot




#fig, axs = plt.subplots(3, 2, figsize=(10, 12), gridspec_kw={'wspace': 0.3, 'hspace': 0.2})









fig, axs = plt.subplots(1, 3, figsize=(16, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)

axs[0].plot(fid_M[0,:], fid_M[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[0].plot(fcl_M[0,:], fcl_M[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[0].plot(fch_M[0,:], fch_M[2,:], linestyle = ls[1], color = col[0], label = lab[1])


axs[0].plot(rcl001_M[0,:], rcl001_M[2,:], linestyle = ls[0], color = col[2], label = lab[4])
#axs[1].plot(rcl002_v[0,:], rcl002_v[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[0].plot(rcl06_M[0,:], rcl06_M[2,:], linestyle = ls[1], color = col[2], label = lab[5])
#axs[0,1].plot(rcl02_v[0,:], rcl02_v[2,:], linestyle = ls[3], color = col[2], label = lab[8])



#axs[1].set_ylim(100,1200)
#axs[1].set_yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')

axs[0].set_ylabel("$\dot{M}_{cold}$  [M$_\odot$/yr]")
axs[0].set_xlabel('$R$   [kpc]')



    ## - Cold

axs[1].plot(fid_v[0,:], fid_v[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[1].plot(fcl_v[0,:], fcl_v[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(fch_v[0,:], fch_v[2,:], linestyle = ls[1], color = col[0], label = lab[1])

#axs[1].plot(NoMix_v[0,:], NoMix_v[2,:], linestyle = ls[0], color = col[1], label = lab[2])
#axs[1].plot(mixl_v[0,:], mixl_v[2,:], linestyle = ls[1], color = col[1], label = lab[2])
#axs[1].plot(mixh_v[0,:], mixh_v[2,:], linestyle = ls[2], color = col[1], label = lab[3])

axs[1].plot(rcl001_v[0,:], rcl001_v[2,:], linestyle = ls[0], color = col[2], label = lab[4])
#axs[1].plot(rcl002_v[0,:], rcl002_v[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[1].plot(rcl06_v[0,:], rcl06_v[2,:], linestyle = ls[1], color = col[2], label = lab[5])
#axs[0,1].plot(rcl02_v[0,:], rcl02_v[2,:], linestyle = ls[3], color = col[2], label = lab[8])

axs[1].legend(
    prop={'family': 'Times New Roman', 'size': 10},
    loc='upper left',  # You can change to 'upper right', etc.
    ncol=1,
)


axs[1].set_ylim(-100,800)
#plt.yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$v_{cold}$  [km/s]")
axs[1].set_xlabel('$R$   [kpc]')



Rmin=0.5

axs[2].plot(fid_sg[0,fid_sg[0,:]>Rmin], fid_sg[1,fid_sg[0,:]>Rmin], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[2].plot(fcl_sg[0,fid_sg[0,:]>Rmin], fcl_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[0], color = col[0], label = lab[0])
axs[2].plot(fch_sg[0,fid_sg[0,:]>Rmin], fch_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[1], color = col[0], label = lab[1])

#plt.plot(mixl_sg[0,:], mixl_sg[1,:], linestyle = ls[1], color = col[1], label = lab[2])
#plt.plot(mixh_sg[0,:], mixh_sg[1,:], linestyle = ls[2], color = col[1], label = lab[3])

axs[2].plot(rcl001_sg[0,fid_sg[0,:]>Rmin], rcl001_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[0], color = col[2], label = lab[4])
axs[2].plot(rcl06_sg[0,fid_sg[0,:]>Rmin], rcl06_sg[1,fid_sg[0,:]>Rmin], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[1,:], linestyle = ls[3], color = col[2], label = lab[8])





axs[2].set_ylim(1e-4,1)
axs[2].set_yscale('log')
axs[2].set_xlabel('$R$   [kpc]')
axs[2].set_ylabel("$r$  [kpc]")







'''
axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
'''


#plt.show()
plt.savefig("CfrParameters")

































































r_fid = np.mean(fid_sg[4,:]/fid_sg[3,:])
r_fcl = np.mean(fcl_sg[4,:]/fcl_sg[3,:])
r_fch = np.mean(fch_sg[4,:]/fch_sg[3,:])
r_mixl = np.mean(mixl_sg[4,:]/mixl_sg[3,:])
r_mixh = np.nanmean(mixh_sg[4,:]/mixh_sg[3,:])
r_rcl001 = np.mean(rcl001_sg[4,:]/rcl001_sg[3,:])
r_rcl06 = np.mean(rcl06_sg[4,:]/rcl06_sg[3,:])


print(mixh_sg[4,:]/mixh_sg[3,:])




fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)

axs[0].plot(fid_v[0,:], fid_v[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[0].plot(fcl_v[0,:], fcl_v[1,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[0].plot(fch_v[0,:], fch_v[1,:], linestyle = ls[1], color = col[0], label = lab[1])

#axs[0].plot(NoMix_v[0,:], NoMix_v[1,:], linestyle = ls[0], color = col[1], label = lab[2])
#axs[0].plot(mixl_v[0,:], mixl_v[1,:], linestyle = ls[1], color = col[1], label = str(round(r_mixl,2)))
#axs[0].plot(mixh_v[0,:], mixh_v[1,:], linestyle = ls[2], color = col[1], label = str(round(r_mixh,2)))

axs[0].plot(rcl001_v[0,:], rcl001_v[1,:], linestyle = ls[0], color = col[2], label = lab[4])
#axs[0].plot(rcl002_v[0,:], rcl002_v[1,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[0].plot(rcl06_v[0,:], rcl06_v[1,:], linestyle = ls[1], color = col[2], label = lab[5])
#axs[0,0].plot(rcl02_v[0,:], rcl02_v[1,:], linestyle = ls[3], color = col[2], label = lab[8])
axs[0].legend(
    prop={'family': 'Times New Roman', 'size': 10},
    loc='lower right',  # You can change to 'upper right', etc.
    ncol=1,
)

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
axs[0].set_ylim(300,1600)
#plt.yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$v$  [km/s]")
#plt.show()


    ## - Cold

axs[1].plot(fid_v[0,:], fid_v[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[1].plot(fcl_v[0,:], fcl_v[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(fch_v[0,:], fch_v[2,:], linestyle = ls[1], color = col[0], label = lab[1])

#axs[1].plot(NoMix_v[0,:], NoMix_v[2,:], linestyle = ls[0], color = col[1], label = lab[2])
#axs[1].plot(mixl_v[0,:], mixl_v[2,:], linestyle = ls[1], color = col[1], label = lab[2])
#axs[1].plot(mixh_v[0,:], mixh_v[2,:], linestyle = ls[2], color = col[1], label = lab[3])

axs[1].plot(rcl001_v[0,:], rcl001_v[2,:], linestyle = ls[0], color = col[2], label = lab[4])
#axs[1].plot(rcl002_v[0,:], rcl002_v[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[1].plot(rcl06_v[0,:], rcl06_v[2,:], linestyle = ls[1], color = col[2], label = lab[5])
#axs[0,1].plot(rcl02_v[0,:], rcl02_v[2,:], linestyle = ls[3], color = col[2], label = lab[8])



axs[1].set_ylim(-100,800)
#plt.yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$v$  [km/s]")
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')



#plt.show()
plt.savefig("CfrParameters_impo2")














fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)

axs[0].plot(fid_n[0,:], fid_n[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = str(round(r_fid,2)))


axs[0].plot(fcl_n[0,:], fcl_n[1,:], linestyle = ls[0], color = col[0], label = str(round(r_fcl,2)))
axs[0].plot(fch_n[0,:], fch_n[1,:], linestyle = ls[1], color = col[0], label = str(round(r_fch,2)))

#axs[0].plot(NoMix_v[0,:], NoMix_v[1,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[0].plot(mixl_n[0,:], mixl_n[1,:], linestyle = ls[1], color = col[1], label = str(round(r_mixl,2)))
axs[0].plot(mixh_n[0,:], mixh_n[1,:], linestyle = ls[2], color = col[1], label = str(round(r_mixh,2)))

axs[0].plot(rcl001_n[0,:], rcl001_n[1,:], linestyle = ls[0], color = col[2], label = str(round(r_rcl001,2)))
#axs[0].plot(rcl002_v[0,:], rcl002_v[1,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[0].plot(rcl06_n[0,:], rcl06_n[1,:], linestyle = ls[2], color = col[2], label = str(round(r_rcl06,2)))
#axs[0,0].plot(rcl02_v[0,:], rcl02_v[1,:], linestyle = ls[3], color = col[2], label = lab[8])
axs[0].legend(
    title=r'$\langle \dot{p}_{\mathrm{drag}} / \dot{p}_{\mathrm{mix}}\rangle$',
    loc='lower right',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=3,
    title_fontsize=11,
)

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
axs[0].set_yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$n$  []")
#plt.show()


    ## - Cold

axs[1].plot(fid_n[0,:], fid_n[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[1].plot(fcl_n[0,:], fcl_n[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(fch_n[0,:], fch_n[2,:], linestyle = ls[1], color = col[0], label = lab[1])

#axs[1].plot(NoMix_v[0,:], NoMix_v[2,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[1].plot(mixl_n[0,:], mixl_n[2,:], linestyle = ls[1], color = col[1], label = lab[2])
axs[1].plot(mixh_n[0,:], mixh_n[2,:], linestyle = ls[2], color = col[1], label = lab[3])

axs[1].plot(rcl001_n[0,:], rcl001_n[2,:], linestyle = ls[0], color = col[2], label = lab[4])
#axs[1].plot(rcl002_v[0,:], rcl002_v[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[1].plot(rcl06_n[0,:], rcl06_n[2,:], linestyle = ls[2], color = col[2], label = lab[5])
#axs[0,1].plot(rcl02_v[0,:], rcl02_v[2,:], linestyle = ls[3], color = col[2], label = lab[8])
axs[1].legend()


#axs[1].set_ylim(100,1200)
axs[1].set_yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$n$  []")
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')



#plt.show()
plt.savefig("CfrParameters_impo2_n")









fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)

axs[0].plot(fid_M[0,:], fid_M[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = str(round(r_fid,2)))


axs[0].plot(fcl_M[0,:], fcl_M[1,:], linestyle = ls[0], color = col[0], label = str(round(r_fcl,2)))
axs[0].plot(fch_M[0,:], fch_M[1,:], linestyle = ls[1], color = col[0], label = str(round(r_fch,2)))

#axs[0].plot(NoMix_v[0,:], NoMix_v[1,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[0].plot(mixl_M[0,:], mixl_M[1,:], linestyle = ls[1], color = col[1], label = str(round(r_mixl,2)))
axs[0].plot(mixh_M[0,:], mixh_M[1,:], linestyle = ls[2], color = col[1], label = str(round(r_mixh,2)))

axs[0].plot(rcl001_M[0,:], rcl001_M[1,:], linestyle = ls[0], color = col[2], label = str(round(r_rcl001,2)))
#axs[0].plot(rcl002_v[0,:], rcl002_v[1,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[0].plot(rcl06_M[0,:], rcl06_M[1,:], linestyle = ls[2], color = col[2], label = str(round(r_rcl06,2)))
#axs[0,0].plot(rcl02_v[0,:], rcl02_v[1,:], linestyle = ls[3], color = col[2], label = lab[8])
axs[0].legend(
    title=r'$\langle \dot{p}_{\mathrm{drag}} / \dot{p}_{\mathrm{mix}}\rangle$',
    loc='lower right',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=3,
    title_fontsize=11,
)

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
#axs[0].set_yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$M$  []")
#plt.show()


    ## - Cold

axs[1].plot(fid_M[0,:], fid_M[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[1].plot(fcl_M[0,:], fcl_M[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(fch_M[0,:], fch_M[2,:], linestyle = ls[1], color = col[0], label = lab[1])

#axs[1].plot(NoMix_v[0,:], NoMix_v[2,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[1].plot(mixl_M[0,:], mixl_M[2,:], linestyle = ls[1], color = col[1], label = lab[2])
axs[1].plot(mixh_M[0,:], mixh_M[2,:], linestyle = ls[2], color = col[1], label = lab[3])

axs[1].plot(rcl001_M[0,:], rcl001_M[2,:], linestyle = ls[0], color = col[2], label = lab[4])
#axs[1].plot(rcl002_v[0,:], rcl002_v[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[1].plot(rcl06_M[0,:], rcl06_M[2,:], linestyle = ls[2], color = col[2], label = lab[5])
#axs[0,1].plot(rcl02_v[0,:], rcl02_v[2,:], linestyle = ls[3], color = col[2], label = lab[8])
axs[1].legend()


#axs[1].set_ylim(100,1200)
#axs[1].set_yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$M$  []")
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')



#plt.show()
plt.savefig("CfrParameters_impo2_M")


sys.exit(-1)









































#### cfr densities
    
    ## - Hot

axs[1,0].plot(fid_n[0,:], fid_n[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[1,0].plot(fcl_n[0,:], fcl_n[1,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1,0].plot(fch_n[0,:], fch_n[1,:], linestyle = ls[1], color = col[0], label = lab[1])

axs[1,0].plot(NoMix_n[0,:], NoMix_n[1,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[1,0].plot(mixl_n[0,:], mixl_n[1,:], linestyle = ls[1], color = col[1], label = lab[3])
axs[1,0].plot(mixh_n[0,:], mixh_n[1,:], linestyle = ls[2], color = col[1], label = lab[4])

axs[1,0].plot(rcl001_n[0,:], rcl001_n[1,:], linestyle = ls[0], color = col[2], label = lab[5])
axs[1,0].plot(rcl002_n[0,:], rcl002_n[1,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[1,0].plot(rcl004_n[0,:], rcl004_n[1,:], linestyle = ls[2], color = col[2], label = lab[7])
#axs[1,0].plot(rcl0.2_n[0,:], rcl0.2_n[1,:], linestyle = ls[3], color = col[2], label = lab[8])




axs[1,0].legend()
axs[1,0].set_ylim(2e-4,4e-1)
axs[1,0].set_yscale('log')
#plt.xlabel('$R$   [kpc]')
axs[1,0].set_ylabel("$n$  [cm$^{-3}$]")
#plt.show()





    ## - Cold

axs[1,1].plot(fid_n[0,:], fid_n[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[1,1].plot(fcl_n[0,:], fcl_n[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1,1].plot(fch_n[0,:], fch_n[2,:], linestyle = ls[1], color = col[0], label = lab[1])

axs[1,1].plot(NoMix_n[0,:], NoMix_n[2,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[1,1].plot(mixl_n[0,:], mixl_n[2,:], linestyle = ls[1], color = col[1], label = lab[3])
axs[1,1].plot(mixh_n[0,:], mixh_n[2,:], linestyle = ls[2], color = col[1], label = lab[4])

axs[1,1].plot(rcl001_n[0,:], rcl001_n[2,:], linestyle = ls[0], color = col[2], label = lab[5])
axs[1,1].plot(rcl002_n[0,:], rcl002_n[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[1,1].plot(rcl004_n[0,:], rcl004_n[2,:], linestyle = ls[2], color = col[2], label = lab[7])
#axs[1,1].plot(rcl02_n[0,:], rcl02_n[2,:], linestyle = ls[3], color = col[2], label = lab[8])

#plt.legend()
axs[1,1].set_ylim(1e-2,2e1)
axs[1,1].set_yscale('log')
#plt.xlabel('$R$   [kpc]')
axs[1,1].set_ylabel("$n$  [cm$^{-3}$]")
#plt.show()







#### cfr Mdots
    
    ## - Hot

axs[2,0].plot(fid_M[0,:], fid_M[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[2,0].plot(fcl_M[0,:], fcl_M[1,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[2,0].plot(fch_M[0,:], fch_M[1,:], linestyle = ls[1], color = col[0], label = lab[1])

axs[2,0].plot(NoMix_M[0,:], NoMix_M[1,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[2,0].plot(mixl_M[0,:], mixl_M[1,:], linestyle = ls[1], color = col[1], label = lab[3])
axs[2,0].plot(mixh_M[0,:], mixh_M[1,:], linestyle = ls[2], color = col[1], label = lab[4])

axs[2,0].plot(rcl001_M[0,:], rcl001_M[1,:], linestyle = ls[0], color = col[2], label = lab[5])
axs[2,0].plot(rcl002_M[0,:], rcl002_M[1,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[2,0].plot(rcl004_M[0,:], rcl004_M[1,:], linestyle = ls[2], color = col[2], label = lab[7])
#axs[2,0].plot(rcl02_M[0,:], rcl02_M[1,:], linestyle = ls[3], color = col[2], label = lab[8])







#plt.legend()

#axs[2,0].set_yscale('log')
axs[2,0].set_xlabel('$R$   [kpc]')
axs[2,0].set_ylim( 5e-2, 1.5)
axs[2,0].set_ylabel("$\dot{M}$   [M$_\odot$/yr]")
##lt.show()





    ## - Cold

axs[2,1].plot(fid_M[0,:], fid_M[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[2,1].plot(fcl_M[0,:], fcl_M[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[2,1].plot(fch_M[0,:], fch_M[2,:], linestyle = ls[1], color = col[0], label = lab[1])

axs[2,1].plot(NoMix_M[0,:], NoMix_M[2,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[2,1].plot(mixl_M[0,:], mixl_M[2,:], linestyle = ls[1], color = col[1], label = lab[3])
axs[2,1].plot(mixh_M[0,:], mixh_M[2,:], linestyle = ls[2], color = col[1], label = lab[4])

axs[2,1].plot(rcl001_M[0,:], rcl001_M[2,:], linestyle = ls[0], color = col[2], label = lab[5])
axs[2,1].plot(rcl002_M[0,:], rcl002_M[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[2,1].plot(rcl004_M[0,:], rcl004_M[2,:], linestyle = ls[2], color = col[2], label = lab[7])
#axs[2,1].plot(rcl02_M[0,:], rcl02_M[2,:], linestyle = ls[3], color = col[2], label = lab[8])



#plt.legend()
#plt.ylim(30,1000)
#axs[2,1].set_yscale('log')
axs[2,1].set_xlabel('$R$   [kpc]')
axs[2,1].set_ylim( 5e-2, 1.5)
axs[2,1].set_ylabel("$\dot{M}$   [M$_\odot$/yr]")
#plt.show()


plt.tight_layout()
#plt.show()
plt.savefig("CfrParameters_impo")



sys.exit(-1)



#### cfr Mdots
    
    ## - Hot

plt.plot(fid_M[0,:], fid_M[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(al_M[0,:], al_M[1,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(ah_M[0,:], ah_M[1,:], linestyle = ls[0], color = col[1], label = lab[1])

plt.plot(ml_M[0,:], ml_M[1,:], linestyle = ls[1], color = col[0], label = lab[2])
plt.plot(mh_M[0,:], mh_M[1,:], linestyle = ls[1], color = col[1], label = lab[3])

plt.plot(Ml_M[0,:], Ml_M[1,:], linestyle = ls[2], color = col[0], label = lab[4])
plt.plot(Mh_M[0,:], Mh_M[1,:], linestyle = ls[2], color = col[1], label = lab[5])

plt.legend()
#plt.ylim(300,2500)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylim( 1e-2, 2)
plt.ylabel("$\dot{M}_{hot}$   [M$_\odot$/yr]")
plt.show()





    ## - Cold

plt.plot(fid_M[0,:], fid_M[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(al_M[0,:], al_M[2,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(ah_M[0,:], ah_M[2,:], linestyle = ls[0], color = col[1], label = lab[1])

plt.plot(ml_M[0,:], ml_M[2,:], linestyle = ls[1], color = col[0], label = lab[2])
plt.plot(mh_M[0,:], mh_M[2,:], linestyle = ls[1], color = col[1], label = lab[3])

plt.plot(Ml_M[0,:], Ml_M[2,:], linestyle = ls[2], color = col[0], label = lab[4])
plt.plot(Mh_M[0,:], Mh_M[2,:], linestyle = ls[2], color = col[1], label = lab[5])

plt.legend()
#plt.ylim(30,1000)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylim( 1e-2, 2)
plt.ylabel("$\dot{M}_{cold}$   [M$_\odot$/yr]")
plt.show()





#### cfr densities
    
    ## - Hot

plt.plot(fid_n[0,:], fid_n[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(al_n[0,:], al_n[1,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(ah_n[0,:], ah_n[1,:], linestyle = ls[0], color = col[1], label = lab[1])

plt.plot(ml_n[0,:], ml_n[1,:], linestyle = ls[1], color = col[0], label = lab[2])
plt.plot(mh_n[0,:], mh_n[1,:], linestyle = ls[1], color = col[1], label = lab[3])

plt.plot(Ml_n[0,:], Ml_n[1,:], linestyle = ls[2], color = col[0], label = lab[4])
plt.plot(Mh_n[0,:], Mh_n[1,:], linestyle = ls[2], color = col[1], label = lab[5])

plt.legend()
#plt.ylim(300,2500)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$n_{hot}$  [cm$^{-3}$]")
plt.show()





    ## - Cold

plt.plot(fid_n[0,:], fid_n[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(al_n[0,:], al_n[2,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(ah_n[0,:], ah_n[2,:], linestyle = ls[0], color = col[1], label = lab[1])

plt.plot(ml_n[0,:], ml_n[2,:], linestyle = ls[1], color = col[0], label = lab[2])
plt.plot(mh_n[0,:], mh_n[2,:], linestyle = ls[1], color = col[1], label = lab[3])

plt.plot(Ml_n[0,:], Ml_n[2,:], linestyle = ls[2], color = col[0], label = lab[4])
plt.plot(Mh_n[0,:], Mh_n[2,:], linestyle = ls[2], color = col[1], label = lab[5])

plt.legend()
#plt.ylim(30,1000)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$n_{cold}$  [cm$^{-3}$]")
plt.show()