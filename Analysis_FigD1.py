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


###### ======= RAINER STUFF CHECK.PY ======== ######## ======== ######## ======== ######## ======== ######## ======== ########
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
    P = np.load("Simulations/"+string+"/output/Resolution.npy")
    sg = np.load("Simulations/"+string+"/output/subgrid.npy")
    return v, n, M, T, P, sg

fid_v, fid_n, fid_M, fid_T, fid_P, fid_sg = load_stuff("fid")
V001_v, V001_n, V001_M, V001_T, V001_P, V001_sg = load_stuff("HR1")
V03_v, V03_n, V03_M, V03_T, V03_P, V03_sg= load_stuff("HR_Vc0.3")
V003_v, V003_n, V003_M, V003_T, V003_P, V003_sg = load_stuff("HR_Vc0.03")
V1_v, V1_n, V1_M, V1_T, V1_P, V1_sg = load_stuff("HR_Vc1")
R_v, R_n, R_M, R_T, R_P, R_sg = load_stuff("high_res_rainer")
#one color 

'''
col = [
    '#87CEEB',
    '#FF7F50',
    '#228B22'
]
'''
col = [
    "#87CEEB",  # Sky Blue
    "#FF7F50",  # Coral
    "#3CB371",  # Medium Sea Green
    "#DA70D6"   # Orchid
]
ls = [
    '-',
    '--',
    '-.',
    ':',
]

lab = ["$0.001$","$0.01$", "$0.03$", "$0.1$ (fid)", "$0.3$", "$1$" ]

'''
import matplotlib.cm as cm
cmap = cm.get_cmap('gist_earth')
positions = np.linspace(0.1, 0.9, 5)
col = [cmap(x) for x in positions]
'''


from matplotlib.colors import ListedColormap
import cmcrameri.cm as cm
#import cmocean
#cmap = cmocean.cm.berlin

# Number of lines
N = 7

# Choose a colormap (e.g., 'viridis', 'plasma', 'tab10', 'jet', etc.)
#cmap = cm.get_cmap('tab20b', N)  # Discretize into N colors
cmapc = ListedColormap(cm.berlin(np.linspace(0, 1, N))) #acton, #lapaz 
col = [cmapc(i) for i in range(N)]
col = ['#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177']



from matplotlib.font_manager import FontProperties

#fig, axs = plt.subplots(3, 2, figsize=(10, 12), gridspec_kw={'wspace': 0.3, 'hspace': 0.2})
fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)
axs[0].plot(R_v[0,:], R_v[1,:], color = col[0], label = lab[0])
axs[0].plot(V001_v[0,:], V001_v[1,:], color = col[1], label = lab[1])
axs[0].plot(V003_v[0,:], V003_v[1,:], color = col[2], label = lab[2])
axs[0].plot(fid_v[0,:], fid_v[1,:], color = col[3], linewidth = 3, label = lab[3])
axs[0].plot(V03_v[0,:], V03_v[1,:], color = col[4], label = lab[4])
axs[0].plot(V1_v[0,:], V1_v[1,:], color = col[5], label = lab[5])


#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
axs[0].set_ylim(250,1600)
#plt.yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$v$  [km/s]")
#plt.show()

axs[0].legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title='$V_{crit}$  [$\\text{kpc}^3$]',
    title_fontproperties=FontProperties(family='Times New Roman', size=10)
)
    ## - Cold

axs[1].plot(R_v[0,:], R_v[2,:], color = col[0], label = lab[0])
axs[1].plot(V001_v[0,:], V001_v[2,:], color = col[1], label = lab[1])
axs[1].plot(V003_v[0,:], V003_v[2,:], color = col[2], label = lab[2])
axs[1].plot(fid_v[0,:], fid_v[2,:], color = col[3], linewidth = 3, label = lab[3])
axs[1].plot(V03_v[0,:], V03_v[2,:], color = col[4], label = lab[4])
axs[1].plot(V1_v[0,:], V1_v[2,:], color = col[5], label = lab[5])




#axs[1].set_ylim(0,1200)
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
plt.savefig("Cfr_hr_vel")







from matplotlib.font_manager import FontProperties

#fig, axs = plt.subplots(3, 2, figsize=(10, 12), gridspec_kw={'wspace': 0.3, 'hspace': 0.2})
fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)
axs[0].plot(R_M[0,:], R_M[1,:], color = col[0], label = lab[0])
axs[0].plot(V001_M[0,:], V001_M[1,:], color = col[1], label = lab[1])
axs[0].plot(V003_M[0,:], V003_M[1,:], color = col[2], label = lab[2])
axs[0].plot(fid_M[0,:], fid_M[1,:], color = col[3], linewidth = 3, label = lab[3])
axs[0].plot(V03_M[0,:], V03_M[1,:], color = col[4], label = lab[4])
axs[0].plot(V1_M[0,:], V1_M[1,:], color = col[5], label = lab[5])


#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#plt.yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$Mdot$")
#plt.show()

axs[0].legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title='$V_{crit}$  [$\\text{kpc}^3$]',
    title_fontproperties=FontProperties(family='Times New Roman', size=10)
)
    ## - Cold

axs[1].plot(R_M[0,:], R_M[2,:], color = col[0], label = lab[0])
axs[1].plot(V001_M[0,:], V001_M[2,:], color = col[1], label = lab[1])
axs[1].plot(V003_M[0,:], V003_M[2,:], color = col[2], label = lab[2])
axs[1].plot(fid_M[0,:], fid_M[2,:], color = col[3], linewidth = 3, label = lab[3])
axs[1].plot(V03_M[0,:], V03_M[2,:], color = col[4], label = lab[4])
axs[1].plot(V1_M[0,:], V1_M[2,:], color = col[5], label = lab[5])




#axs[1].set_ylim(0,1200)
#plt.yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$Mdot$ ")
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')

#plt.show()
plt.savefig("Cfr_hr_M")







from matplotlib.font_manager import FontProperties

#fig, axs = plt.subplots(3, 2, figsize=(10, 12), gridspec_kw={'wspace': 0.3, 'hspace': 0.2})
fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)
axs[0].plot(R_n[0,:], R_n[1,:], color = col[0], label = lab[0])
axs[0].plot(V001_n[0,:], V001_n[1,:], color = col[1], label = lab[1])
axs[0].plot(V003_n[0,:], V003_n[1,:], color = col[2], label = lab[2])
axs[0].plot(fid_n[0,:], fid_n[1,:], color = col[3], linewidth = 3, label = lab[3])
axs[0].plot(V03_n[0,:], V03_n[1,:], color = col[4], label = lab[4])
axs[0].plot(V1_n[0,:], V1_n[1,:], color = col[5], label = lab[5])


#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#plt.yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$n$")
#plt.show()

axs[0].legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title='$V_{crit}$  [$\\text{kpc}^3$]',
    title_fontproperties=FontProperties(family='Times New Roman', size=10)
)
    ## - Cold

axs[1].plot(R_n[0,:], R_n[2,:], color = col[0], label = lab[0])
axs[1].plot(V001_n[0,:], V001_n[2,:], color = col[1], label = lab[1])
axs[1].plot(V003_n[0,:], V003_n[2,:], color = col[2], label = lab[2])
axs[1].plot(fid_n[0,:], fid_n[2,:], color = col[3], linewidth = 3, label = lab[3])
axs[1].plot(V03_n[0,:], V03_n[2,:], color = col[4], label = lab[4])
axs[1].plot(V1_n[0,:], V1_n[2,:], color = col[5], label = lab[5])


#axs[1].set_ylim(0,1200)
#plt.yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$n$ ")
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')

#plt.show()
plt.savefig("Cfr_hr_n")

































fig, ax = plt.subplots(figsize=(6, 5),  dpi=300)
Rmin = 0.5
plt.plot(R_P[0, R_P[0,:]>Rmin], R_P[1,R_P[0,:]>Rmin], color = col[0], label = lab[0])
plt.plot(V001_P[0,R_P[0,:]>Rmin], V001_P[1,R_P[0,:]>Rmin], color = col[1], label = lab[1])
plt.plot(V003_P[0,R_P[0,:]>Rmin], V003_P[1,R_P[0,:]>Rmin], color = col[2], label = lab[2])
plt.plot(fid_P[0,R_P[0,:]>Rmin], fid_P[1,R_P[0,:]>Rmin], color = col[3], linewidth = 3, label = lab[3])
plt.plot(V03_P[0,R_P[0,:]>Rmin], V03_P[1,R_P[0,:]>Rmin], color = col[4], label = lab[4])
plt.plot(V1_P[0,R_P[0,:]>Rmin], V1_P[1,R_P[0,:]>Rmin], color = col[5], label = lab[5])


plt.plot(R_P[0,R_P[0,:]>Rmin], R_P[2,R_P[0,:]>Rmin],linestyle = '--', color = col[0])
plt.plot(V001_P[0,R_P[0,:]>Rmin], V001_P[2,R_P[0,:]>Rmin], linestyle = '--', color = col[1])
plt.plot(V003_P[0,R_P[0,:]>Rmin], V003_P[2,R_P[0,:]>Rmin], linestyle = '--', color = col[2])
plt.plot(fid_P[0,R_P[0,:]>Rmin], fid_P[2,R_P[0,:]>Rmin],  linestyle = '--',color = col[3], linewidth = 3)
plt.plot(V03_P[0,R_P[0,:]>Rmin], V03_P[2,R_P[0,:]>Rmin],  linestyle = '--',color = col[4])
plt.plot(V1_P[0,R_P[0,:]>Rmin], V1_P[2,R_P[0,:]>Rmin],  linestyle = '--',color = col[5])

plt.legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title='$V_{crit}$  [kpc$^3$]',
    title_fontproperties=FontProperties(family='Times New Roman', size=11),
    ncol=2
)
#axs[0].legend()
plt.xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
plt.yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
plt.ylabel("$l$  [kpc]")
#plt.show()
plt.ylim(0.01,5)
plt.savefig("Cfr_hr_P")

axs[1].legend()
sys.exit(-1)

































#fig, axs = plt.subplots(3, 2, figsize=(10, 12), gridspec_kw={'wspace': 0.3, 'hspace': 0.2})
fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)
axs[0].plot(fid_n[0,:], fid_n[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')
axs[0].plot(hr_n[0,:], hr_n[1,:], linestyle = ls[0], color = col[0], label = lab[0])

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
axs[0].set_yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$v$  [km/s]")
#plt.show()


    ## - Cold

axs[1].plot(fid_n[0,:], fid_n[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[1].plot(hr_n[0,:], hr_n[2,:], linestyle = ls[0], color = col[0], label = lab[0])

axs[1].legend()


#axs[1].set_ylim(0,1200)
axs[1].set_yscale('log')
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
plt.savefig("Cfr_hr_n")




#fig, axs = plt.subplots(3, 2, figsize=(10, 12), gridspec_kw={'wspace': 0.3, 'hspace': 0.2})
fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)
axs[0].plot(fid_M[0,:], fid_M[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')
axs[0].plot(hr_M[0,:], hr_M[1,:], linestyle = ls[0], color = col[0], label = lab[0])

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
#axs[0].yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$v$  [km/s]")
#plt.show()


    ## - Cold

axs[1].plot(fid_M[0,:], fid_M[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[1].plot(hr_M[0,:], hr_M[2,:], linestyle = ls[0], color = col[0], label = lab[0])

axs[1].legend()


#axs[1].set_ylim(0,1200)
#axs[1].yscale('log')
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
plt.savefig("Cfr_hr_M")
sys.exit(-1)









fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)
axs[0].plot(fid_n[0,:], fid_v[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[0].plot(m2_n[0,:], m2_n[1,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[0].plot(m4_n[0,:], m4_n[1,:], linestyle = ls[1], color = col[1], label = lab[1])
axs[0].plot(m6_n[0,:], m6_n[1,:], linestyle = ls[2], color = col[2], label = lab[2])
axs[0].plot(m8_n[0,:], m8_n[1,:], linestyle = ls[3], color = col[3], label = lab[3])

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
axs[0].set_yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$n$  []")
#plt.show()



    ## - Cold

axs[1].plot(fid_n[0,:], fid_n[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[1].plot(m2_n[0,:], m2_n[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(m4_n[0,:], m4_n[2,:], linestyle = ls[1], color = col[1], label = lab[1])
axs[1].plot(m6_n[0,:], m6_n[2,:], linestyle = ls[2], color = col[2], label = lab[2])
axs[1].plot(m8_n[0,:], m8_n[2,:], linestyle = ls[3], color = col[3], label = lab[3])
axs[1].legend()


#axs[1].set_ylim(0,1200)
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
plt.savefig("Cfr_mmin_n")

















fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)
axs[0].plot(fid_sg[0,:], fid_sg[3,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[0].plot(m2_sg[0,:], m2_sg[3,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[0].plot(m4_sg[0,:], m4_sg[3,:], linestyle = ls[1], color = col[1], label = lab[1])
axs[0].plot(m6_sg[0,:], m6_sg[3,:], linestyle = ls[2], color = col[2], label = lab[2])
axs[0].plot(m8_sg[0,:], m8_sg[3,:], linestyle = ls[3], color = col[3], label = lab[3])

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
#axs[0].set_yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$mix$  []")
#plt.show()



    ## - Cold

axs[1].plot(fid_sg[0,:], fid_sg[4,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[1].plot(m2_sg[0,:], m2_sg[4,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(m4_sg[0,:], m4_sg[4,:], linestyle = ls[1], color = col[1], label = lab[1])
axs[1].plot(m6_sg[0,:], m6_sg[4,:], linestyle = ls[2], color = col[2], label = lab[2])
axs[1].plot(m8_sg[0,:], m8_sg[4,:], linestyle = ls[3], color = col[3], label = lab[3])
axs[1].legend()


#axs[1].set_ylim(0,1200)
#axs[1].set_yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$drag$  []")
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')

#plt.show()
plt.savefig("Cfr_mmin_pdot")



fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)
axs[0].plot(fid_sg[0,:], fid_sg[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[0].plot(m2_sg[0,:], m2_sg[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[0].plot(m4_sg[0,:], m4_sg[2,:], linestyle = ls[1], color = col[1], label = lab[1])
axs[0].plot(m6_sg[0,:], m6_sg[2,:], linestyle = ls[2], color = col[2], label = lab[2])
axs[0].plot(m8_sg[0,:], m8_sg[2,:], linestyle = ls[3], color = col[3], label = lab[3])

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
#axs[0].set_ylim(250,1750)
#axs[0].set_yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$mcl$  []")
#plt.show()



    ## - Cold

axs[1].plot(fid_sg[0,:], fid_sg[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')


axs[1].plot(m2_sg[0,:], m2_sg[1,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(m4_sg[0,:], m4_sg[1,:], linestyle = ls[1], color = col[1], label = lab[1])
axs[1].plot(m6_sg[0,:], m6_sg[1,:], linestyle = ls[2], color = col[2], label = lab[2])
axs[1].plot(m8_sg[0,:], m8_sg[1,:], linestyle = ls[3], color = col[3], label = lab[3])
axs[1].legend()


#axs[1].set_ylim(0,1200)
#axs[1].set_yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[1].set_ylabel("$rcl$  []")
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')

#plt.show()
plt.savefig("Cfr_mmin_mclrcl")
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