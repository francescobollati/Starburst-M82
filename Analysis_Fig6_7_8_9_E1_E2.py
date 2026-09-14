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
    sg = np.load("Simulations/"+string+"/output/subgrid.npy")
    sc = np.load("Simulations/"+string+"/output/scalars.npy")
    a = np.load("Simulations/"+string+"/output/accelerations.npy")
    return v, n, M, T, sg, sc, a

fid_v, fid_n, fid_M, fid_T, fid_sg, fid_sc, fid_a = load_stuff("fid")
mix003_v, mix003_n, mix003_M,mix003_T, mix003_sg, mix003_sc, mix003_a = load_stuff("fc0.03") 
mix01_v, mix01_n, mix01_M,mix01_T, mix01_sg, mix01_sc, mix01_a = load_stuff("fc0.1") 
mix03_v, mix03_n, mix03_M,mix03_T, mix03_sg, mix03_sc, mix03_a = load_stuff("fc0.3") 
mix1_v, mix1_n, mix1_M,mix1_T, mix1_sg, mix1_sc, mix1_a = load_stuff("fc1")
nd_v, nd_n, nd_M, nd_T, nd_sg, nd_sc, nd_a = load_stuff("fc0.03_nd")

def load_T_NoS():
    T = np.load("Simulations/15J_NoSources/output/Temperatures.npy")
    return T 

NoS_T = load_T_NoS()




import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import cmcrameri.cm as cm
#import cmocean
#cmap = cmocean.cm.berlin

# Number of lines
N = 5

# Choose a colormap (e.g., 'viridis', 'plasma', 'tab10', 'jet', etc.)
#cmap = cm.get_cmap('tab20b', N)  # Discretize into N colors
cmapc = ListedColormap(cm.batlow(np.linspace(0, 1, N))) #acton, #lapaz 
cmaph = ListedColormap(cm.lajolla(np.linspace(0, 1, N)))
colorsc = [cmapc(i) for i in range(N)]
colorsh = [cmaph(i) for i in range(N)]


colorsc = ['#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494']#'#ffffcc',
colorsc = colorsc[::-1]
#colorsc = ['#fcc5c0','#fa9fb5','#f768a1','#c51b8a','#7a0177']
lab = ["0.003", "0.03", "0.1", "0.3", "1"]




################### ===== VELOCITIES ====== ########################


r_fid = np.mean(fid_sg[4,:]/fid_sg[3,:])
r_003 = np.mean(mix003_sg[4,:]/mix003_sg[3,:])
r_01 = np.mean(mix01_sg[4,:]/mix01_sg[3,:])
r_03 = np.mean(mix03_sg[4,:]/mix03_sg[3,:])
r_1 = np.mean(mix1_sg[4,:]/mix1_sg[3,:])




fig, axs = plt.subplots(2, 2, figsize=(12, 10), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)


axs[0,0].plot(fid_v[0,:], fid_v[1,:], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
#axs[0].plot(mix001_v[0,:], mix001_v[1,:], linestyle = 'solid', color = colors[1], label = str(round(r_001,2)))
axs[0,0].plot(mix003_v[0,:], mix003_v[1,:], linestyle = 'solid', color = colorsc[1], label = lab[1])
#axs[0].plot(mix005_v[0,:], mix005_v[1,:], linestyle = 'solid', color = colors[3], label = str(round(r_005,2)))
#axs[0].plot(mix006_v[0,:], mix006_v[1,:], linestyle = 'solid', color = colors[4], label = str(round(r_006,2)))
#axs[0].plot(mix008_v[0,:], mix008_v[1,:], linestyle = 'solid', color = colors[5], label = str(round(r_008,2)))
axs[0,0].plot(mix01_v[0,:], mix01_v[1,:], linestyle = 'solid', color = colorsc[2], label = lab[2])
axs[0,0].plot(mix03_v[0,:], mix03_v[1,:], linestyle = 'solid', color = colorsc[3], label = lab[3])
axs[0,0].plot(mix1_v[0,:], mix1_v[1,:], linestyle = 'solid', color = colorsc[4], label = lab[4])
#axs[0].plot(mix003nd_v[0,:], mix003nd_v[1,:], linestyle = '--', color = colors[2], label = str(round(r_003_nd,2)))



axs[0,0].legend(
    title=r'$f_{mix}$',
    loc='lower right',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)

#axs[0].legend()
#axs[0].set_xlabel('$R$   [kpc]')
axs[0,0].set_ylim(70,1400)
#plt.yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0,0].set_ylabel("$v$  [km/s]")
#plt.show()


    ## - Cold

axs[0,1].plot(fid_v[0,:], fid_v[2,:], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = str(round(r_fid,2)))
#axs[1].plot(mix001_v[0,:], mix001_v[2,:], linestyle = 'solid', color = colors[1], label = lab[1])
axs[0,1].plot(mix003_v[0,:], mix003_v[2,:], linestyle = 'solid', color = colorsc[1], label = str(round(r_003,2)))
#axs[1].plot(mix005_v[0,:], mix005_v[2,:], linestyle = 'solid', color = colors[3], label = lab[3])
#axs[1].plot(mix006_v[0,:], mix006_v[2,:], linestyle = 'solid', color = colors[4], label = lab[4])
#axs[1].plot(mix008_v[0,:], mix008_v[2,:], linestyle = 'solid', color = colors[5], label = lab[5])
axs[0,1].plot(mix01_v[0,:], mix01_v[2,:], linestyle = 'solid', color = colorsc[2], label = str(round(r_01,2)))
axs[0,1].plot(mix03_v[0,:], mix03_v[2,:], linestyle = 'solid', color = colorsc[3], label = str(round(r_03,2)))
axs[0,1].plot(mix1_v[0,:], mix1_v[2,:], linestyle = 'solid', color = colorsc[4], label = str(round(r_1,2)))
#axs[1].plot(mix003nd_v[0,:], mix003nd_v[2,:], linestyle = '--', color = colors[2], label = lab[8])

'''
axs[0,1].legend(
    title=r'$\langle \dot{p}_{\mathrm{drag}} / \dot{p}_{\mathrm{mix}}\rangle$',
    loc='upper left',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)
'''


axs[0,1].set_ylim(100,800)
#plt.yscale('log')
#axs[0,1].set_xlabel('$R$   [kpc]')
axs[0,1].set_ylabel("$v$  [km/s]")
#axs[1].set_xlabel('$R$   [kpc]')


axs[0,0].text(0.5, 1.05, 'Hot phase', transform=axs[0,0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[0,1].text(0.5, 1.05
, 'Cold phase', transform=axs[0,1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')





m_fid = np.mean(np.abs(fid_M[3,fid_M[0,:]>2])/fid_M[1,fid_M[0,:]>2])
m_003 = np.mean(np.abs(mix003_M[3,fid_M[0,:]>2])/mix003_M[1,fid_M[0,:]>2])
m_01 = np.mean(np.abs(mix01_M[3,fid_M[0,:]>2])/mix01_M[1,fid_M[0,:]>2])
m_03 = np.mean(np.abs(mix03_M[3,fid_M[0,:]>2])/mix03_M[1,fid_M[0,:]>2])
m_1 = np.mean(np.abs(mix1_M[3,fid_M[0,:]>2])/mix1_M[1,fid_M[0,:]>2])

axs[1,0].plot(fid_M[0,:], fid_M[1,:], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = str(round(m_fid,2)))
#axs[0].plot(mix001_v[0,:], mix001_v[1,:], linestyle = 'solid', color = colors[1], label = str(round(r_001,2)))
axs[1,0].plot(mix003_M[0,:], mix003_M[1,:], linestyle = 'solid', color = colorsc[1], label = str(round(m_003,2)))
#axs[0].plot(mix005_v[0,:], mix005_v[1,:], linestyle = 'solid', color = colors[3], label = str(round(r_005,2)))
#axs[0].plot(mix006_v[0,:], mix006_v[1,:], linestyle = 'solid', color = colors[4], label = str(round(r_006,2)))
#axs[0].plot(mix008_v[0,:], mix008_v[1,:], linestyle = 'solid', color = colors[5], label = str(round(r_008,2)))
axs[1,0].plot(mix01_M[0,:], mix01_M[1,:], linestyle = 'solid', color = colorsc[2], label = str(round(m_01,2)))
axs[1,0].plot(mix03_M[0,:], mix03_M[1,:], linestyle = 'solid', color = colorsc[3], label = str(round(m_03,2)))
axs[1,0].plot(mix1_M[0,:], mix1_M[1,:], linestyle = 'solid', color = colorsc[4], label = str(round(m_1,2)))
axs[1,0].set_ylabel('$\dot{M}$   [$M_\odot$/yr]')
axs[1,0].set_xlabel('$R$   [kpc]')

'''
axs[1,0].legend(
    title=r'$\langle \dot{M}_{loss} / \dot{M}_{hot}\rangle$',
    loc='upper left',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)
'''

axs[1,1].plot(fid_M[0,:], fid_M[2,:], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = str(round(r_fid,2)))
#axs[0].plot(mix001_v[0,:], mix001_v[1,:], linestyle = 'solid', color = colors[1], label = str(round(r_001,2)))
axs[1,1].plot(mix003_M[0,:], mix003_M[2,:], linestyle = 'solid', color = colorsc[1], label = str(round(r_003,2)))
#axs[0].plot(mix005_v[0,:], mix005_v[1,:], linestyle = 'solid', color = colors[3], label = str(round(r_005,2)))
#axs[0].plot(mix006_v[0,:], mix006_v[1,:], linestyle = 'solid', color = colors[4], label = str(round(r_006,2)))
#axs[0].plot(mix008_v[0,:], mix008_v[1,:], linestyle = 'solid', color = colors[5], label = str(round(r_008,2)))
axs[1,1].plot(mix01_M[0,:], mix01_M[2,:], linestyle = 'solid', color = colorsc[2], label = str(round(r_01,2)))
axs[1,1].plot(mix03_M[0,:], mix03_M[2,:], linestyle = 'solid', color = colorsc[3], label = str(round(r_03,2)))
axs[1,1].plot(mix1_M[0,:], mix1_M[2,:], linestyle = 'solid', color = colorsc[4], label = str(round(r_1,2)))

axs[1,1].set_ylabel('$\dot{M}$   [$M_\odot$/yr]')
axs[1,1].set_xlabel('$R$   [kpc]')

#plt.show()
plt.savefig("Cfr_mix_vel")







fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)


Rmin = 0
fid_mch = np.min(fid_T[2,fid_T[0,:]>Rmin])
m003_mch = np.min(mix003_T[2,mix003_T[0,:]>Rmin])
m01_mch = np.min(mix01_T[2,mix01_T[0,:]>Rmin])
m03_mch = np.min(mix03_T[2,mix03_T[0,:]>Rmin])
m1_mch = np.min(mix1_T[2,mix1_T[0,:]>Rmin])




axs[0].plot(fid_T[0,fid_T[0,:]>Rmin], fid_T[1,fid_T[0,:]>Rmin], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = str(round(fid_mch,2)))
#axs[0].plot(mix001_T[0,:], mix001_T[1,:], linestyle = 'solid', color = colors[1], label = lab[1])
axs[0].plot(mix003_T[0,fid_T[0,:]>Rmin], mix003_T[1,fid_T[0,:]>Rmin], linestyle = 'solid', color = colorsc[1], label = str(round(m003_mch,2)))
#axs[0].plot(mix005_T[0,:], mix005_T[1,:], linestyle = 'solid', color = colors[3], label = lab[3])
#axs[0].plot(mix006_T[0,:], mix006_T[1,:], linestyle = 'solid', color = colors[4], label = lab[4])
#axs[0].plot(mix008_T[0,:], mix008_T[1,:], linestyle = 'solid', color = colors[5], label = lab[5])
axs[0].plot(mix01_T[0,fid_T[0,:]>Rmin], mix01_T[1,fid_T[0,:]>Rmin], linestyle = 'solid', color = colorsc[2], label = str(round(m01_mch,2)))
axs[0].plot(mix03_T[0,fid_T[0,:]>Rmin], mix03_T[1,fid_T[0,:]>Rmin], linestyle = 'solid', color = colorsc[3], label = str(round(m03_mch,2)))
axs[0].plot(mix1_T[0,fid_T[0,:]>Rmin], mix1_T[1,fid_T[0,:]>Rmin], linestyle = 'solid', color = colorsc[4], label = str(round(m1_mch,2)))
#axs[0].plot(mix003nd_T[0,:], mix003nd_T[1,:], linestyle = '--', color = colors[2], label = lab[8])

axs[0].plot(NoS_T[0,NoS_T[0,:]>Rmin], NoS_T[1,NoS_T[0,:]>Rmin], linestyle = 'solid', color = 'k')

axs[0].legend(
    title=r'min$(\mathcal{M}_{hot})$',
    loc='lower left',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)
axs[0].set_yscale('log')
axs[0].set_ylabel('$T$   [K]')
axs[0].set_xlabel('$R$   [kpc]')
#plt.show()



axs[1].plot(fid_sg[0,fid_sg[0,:]>Rmin], -fid_sg[11,fid_sg[0,:]>Rmin]+fid_sg[5,fid_sg[0,:]>Rmin]+fid_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
axs[1].plot(mix003_sg[0,fid_sg[0,:]>Rmin], -mix003_sg[11,fid_sg[0,:]>Rmin]+mix003_sg[5,fid_sg[0,:]>Rmin]+mix003_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[1], label = lab[1] )
axs[1].plot(mix01_sg[0,fid_sg[0,:]>Rmin], -mix01_sg[11,fid_sg[0,:]>Rmin]+mix01_sg[5,fid_sg[0,:]>Rmin]+mix01_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[2], label = lab[2])
axs[1].plot(mix03_sg[0,fid_sg[0,:]>Rmin], -mix03_sg[11,fid_sg[0,:]>Rmin]+mix03_sg[5,fid_sg[0,:]>Rmin]+mix03_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[3], label = lab[3])
axs[1].plot(mix1_sg[0,fid_sg[0,:]>Rmin], -mix1_sg[11,fid_sg[0,:]>Rmin]+mix1_sg[5,fid_sg[0,:]>Rmin]+mix1_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[4], label = lab[4])

axs[1].legend(
    title=r'$f_{mix}$',    
    loc='lower left',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)
#axs[1].plot(fid_sg[0,:], fid_sg[6,:], linestyle = '--', color = colors[0], linewidth = 3, label = lab[0])
#axs[1].plot(mix003_sg[0,:], mix003_sg[6,:], linestyle = '--', color = colors[2], label = lab[2])
#axs[1].plot(mix01_sg[0,:], mix01_sg[6,:], linestyle = '--', color = colors[6], label = lab[6])
#axs[1].plot(mix03_sg[0,:], mix03_sg[6,:], linestyle = '--', color = colors[7], label = lab[7])

axs[1].set_yscale('symlog', linthresh=0.01)
axs[1].set_ylabel('$\dot{u}_{tot}$   [(cm/s)$^2$s$^{-1}$]')
axs[1].set_xlabel('$R$   [kpc]')


#plt.show()
plt.savefig("Cfr_mix_Tudot")







fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)


Rmin = 0.5
m_fid = np.mean(np.abs(fid_a[5,fid_a[0,:]>Rmin])/fid_a[4,fid_a[0,:]>Rmin])
m_003 = np.mean(np.abs(mix003_a[5,fid_a[0,:]>Rmin])/mix003_a[4,fid_a[0,:]>Rmin])
m_01 = np.mean(np.abs(mix01_a[5,fid_a[0,:]>Rmin])/mix01_a[4,fid_a[0,:]>Rmin])
m_03 = np.mean(np.abs(mix03_a[5,fid_a[0,:]>Rmin])/mix03_a[4,fid_a[0,:]>Rmin])
m_1 = np.mean(np.abs(mix1_a[5,fid_a[0,:]>Rmin])/mix1_a[4,fid_a[0,:]>Rmin])




l1, = axs[0].plot(fid_a[0,fid_a[0,:]>Rmin], fid_a[2,fid_a[0,:]>Rmin] + fid_a[3,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
l2, =axs[0].plot(mix003_a[0,fid_a[0,:]>Rmin], mix003_a[2,fid_a[0,:]>Rmin] + mix003_a[3,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[1], label = lab[1])
l3, =axs[0].plot(mix01_a[0,fid_a[0,:]>Rmin], mix01_a[2,fid_a[0,:]>Rmin] + mix01_a[3,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[2], label = lab[2])
l4, =axs[0].plot(mix03_a[0,fid_a[0,:]>Rmin], mix03_a[2,fid_a[0,:]>Rmin]+ mix03_a[3,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[3], label =lab[3] )
l5, =axs[0].plot(mix1_a[0,fid_a[0,:]>Rmin], mix1_a[2,fid_a[0,:]>Rmin] + mix1_a[3,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[4], label = lab[4])



axs[0].plot(fid_a[0,fid_a[0,:]>Rmin], fid_a[1,fid_a[0,:]>Rmin], linestyle = '--', color = colorsc[0], linewidth = 3)
axs[0].plot(mix003_a[0,fid_a[0,:]>Rmin], mix003_a[1,fid_a[0,:]>Rmin], linestyle = '--', color = colorsc[1])
axs[0].plot(mix01_a[0,fid_a[0,:]>Rmin], mix01_a[1,fid_a[0,:]>Rmin], linestyle = '--', color = colorsc[2])
axs[0].plot(mix03_a[0,fid_a[0,:]>Rmin], mix03_a[1,fid_a[0,:]>Rmin], linestyle = '--', color = colorsc[3])
axs[0].plot(mix1_a[0,fid_a[0,:]>Rmin], mix1_a[1,fid_a[0,:]>Rmin], linestyle = '--', color = colorsc[4])




'''
axs[0].legend(
    title=r'$f_{mix}$',    
    loc='lower right',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)
'''
legend1 = axs[0].legend(
    handles=[l1, l2, l3, l4, l5],
    title=r'$f_{mix}$',
    loc='lower right',  # change to 'upper right', etc., if needed
    fontsize=10,
    ncol=2,
    title_fontsize=10
)

# Correct way to add the legend as an artist
axs[0].add_artist(legend1)

from matplotlib.lines import Line2D
line1 = Line2D([], [], linestyle='solid', color='k', linewidth=1, label="$a_{drag}+a_{mix}$")
line2 = Line2D([], [], linestyle='--', color='k', linewidth=1, label="$a_{hydro}$")
legend2 = axs[0].legend(handles=[line1,line2], loc='upper right', fontsize=10)


#axs[0].set_yscale('symlog')#, linthresh=)
axs[0].set_ylabel('$a$   [cm/s$^2$]')
axs[0].set_xlabel('$R$   [kpc]')
#plt.show()

axs[1].plot(fid_a[0,fid_a[0,:]>Rmin], fid_a[4,fid_a[0,:]>Rmin]+fid_a[5,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = str(round(m_fid,2)))
axs[1].plot(mix003_a[0,fid_a[0,:]>Rmin], mix003_a[4,fid_a[0,:]>Rmin]+mix003_a[5,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[1], label = str(round(m_003,2)) )
axs[1].plot(mix01_a[0,fid_a[0,:]>Rmin], mix01_a[4,fid_a[0,:]>Rmin]+mix01_a[5,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[2], label = str(round(m_01,2)))
axs[1].plot(mix03_a[0,fid_a[0,:]>Rmin],mix03_a[4,fid_a[0,:]>Rmin]+mix03_a[5,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[3], label = str(round(m_03,2)))
axs[1].plot(mix1_a[0,fid_a[0,:]>Rmin], mix1_a[4,fid_a[0,:]>Rmin]+mix1_a[5,fid_a[0,:]>Rmin], linestyle = 'solid', color = colorsc[4], label = str(round(m_1,2)))

axs[1].legend(
    title=r'$\langle a_{c,mix}/a_{c,drag}\rangle$',
    loc='upper right',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)

#axs[1].set_yscale('symlog', linthresh=0.01)
axs[1].set_ylabel('$a$   [cm/s$^2$]')
axs[1].set_xlabel('$R$   [kpc]')

axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')


#plt.show()
plt.savefig("Cfr_mix_acc")











fig, ax = plt.subplots(figsize=(6, 5),  dpi=300)


#maybe in legend I show average fraction 

rr1 = np.mean(fid_M[5,:]/fid_M[2,:])
rr2 = np.mean(mix003_M[5,:]/mix003_M[2,:])
rr3 = np.mean(mix01_M[5,:]/mix01_M[2,:])
rr4 = np.mean(mix03_M[5,:]/mix03_M[2,:])
rr5 = np.mean(mix1_M[5,:]/mix1_M[2,:])

plt.plot(fid_M[0,:], fid_M[5,:], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
plt.plot(mix003_M[0,:], mix003_M[5,:], linestyle = 'solid', color = colorsc[1], label = lab[1])
plt.plot(mix01_M[0,:], mix01_M[5,:], linestyle = 'solid', color = colorsc[2], label = lab[2])
plt.plot(mix03_M[0,:], mix03_M[5,:], linestyle = 'solid', color = colorsc[3], label = lab[3])
plt.plot(mix1_M[0,:], mix1_M[5,:], linestyle = 'solid', color = colorsc[4], label = lab[4])

plt.ylabel('$\dot{M}$   [$M_\odot$/yr]')
plt.xlabel('$R$   [kpc]')
legend1 = plt.legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title_fontsize = 10 ,
    title='$f_{mix}$',
    ncol=2,
    loc='upper left',
)


plt.gca().add_artist(legend1)


from matplotlib.lines import Line2D
line1 = Line2D([], [], linewidth=3,  linestyle='solid', color=colorsc[0], label=str(round(rr1,2)))
line2 = Line2D([], [], linestyle='solid', color=colorsc[1], linewidth=1, label=str(round(rr2,2)))
line3 = Line2D([], [], linestyle='solid', color=colorsc[2], linewidth=1, label=str(round(rr3,2)))
line4 = Line2D([], [], linestyle='solid', color=colorsc[3], linewidth=1, label=str(round(rr4,2)))
line5 = Line2D([], [], linestyle='solid', color=colorsc[4], linewidth=1, label=str(round(rr5,2)))
legend2 = plt.legend(handles=[line1, line2, line3, line4, line5], loc='upper right',  prop={'family': 'Times New Roman', 'size': 10}, ncol=2, title="$\langle \dot{M}_{IC,cold}/\dot{M}_{cold} \\rangle$", title_fontsize=10)


plt.ylim(0,1)
#axs[0].legend()
plt.ylabel('$\dot{M}_{IC,cold}$  [M$_\odot$/yr]')
plt.xlabel('$R$   [kpc]')

#plt.show()
plt.savefig("Cfr_mix_Pcoldout")










fig, ax = plt.subplots(figsize=(6, 5),  dpi=300)

plt.plot(fid_sc[0,:], fid_sc[4,:], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
plt.plot(mix003_sc[0,:], mix003_sc[4,:], linestyle = 'solid', color = colorsc[1], label = lab[1])
plt.plot(mix01_sc[0,:], mix01_sc[4,:], linestyle = 'solid', color = colorsc[2], label = lab[2])
plt.plot(mix03_sc[0,:], mix03_sc[4,:], linestyle = 'solid', color = colorsc[3], label = lab[3])
plt.plot(mix1_sc[0,:], mix1_sc[4,:], linestyle = 'solid', color = colorsc[4], label = lab[4])
#plt.plot(nd_sc[0,:], nd_sc[4,:], linestyle = '--', color = 'k', label = "NoDrag_fmix0.03")

#plt.plot(fid_sc[0,:], fid_sc[7,:], linestyle = '--', color = colorsc[0], linewidth = 3, label = "SNe")
#plt.plot(mix003_sc[0,:], mix003_sc[7,:], linestyle = '--', color = colorsc[1])
#plt.plot(mix01_sc[0,:], mix01_sc[7,:], linestyle = '--', color = colorsc[2])
#plt.plot(mix03_sc[0,:], mix03_sc[7,:], linestyle = '--', color = colorsc[3])

#plt.legend(loc='upper right')


legend1 = plt.legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title_fontsize = 10 ,
    title='$f_{mix}$',
    ncol=2,
    loc='lower left',
)

plt.gca().add_artist(legend1)

#line2, = plt.plot(mix1_sc[0,:], mix1_sc[7,:], linestyle = ':', color = colorsc[4], label = '$\mathcal{P}_{cold,2}$')
#line3, = plt.plot(mix1_sc[0,:], mix1_sc[1,:], linestyle = '--', color = colorsc[4],label = '$\mathcal{P}_{cold,0}$')
#from matplotlib.lines import Line2D
#line1 = Line2D([], [], linestyle='solid', color=colorsc[4], linewidth=1, label="$\mathcal{P}_{cold,1}$")
#legend2 = plt.legend(handles=[line3, line1, line2], loc='lower right', fontsize=10)


plt.ylim(0,1)
#axs[0].legend()
plt.ylabel('$\mathcal{P}_{cold,1}$')
plt.xlabel('$R$   [kpc]')

#plt.show()
plt.savefig("Cfr_mix_scalarscold_paper")






fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)

# --- First panel ---
ax[0].plot(fid_sc[0,:], fid_sc[4,:], linestyle='solid', color=colorsc[0], linewidth=3, label=lab[0])
ax[0].plot(mix003_sc[0,:], mix003_sc[4,:], linestyle='solid', color=colorsc[1], label=lab[1])
ax[0].plot(nd_sc[0,:], nd_sc[4,:], linestyle='--', color='k', label='0.03  (NoDrag)')
ax[0].plot(mix01_sc[0,:], mix01_sc[4,:], linestyle='solid', color=colorsc[2], label=lab[2])
ax[0].plot(mix03_sc[0,:], mix03_sc[4,:], linestyle='solid', color=colorsc[3], label=lab[3])
ax[0].plot(mix1_sc[0,:], mix1_sc[4,:], linestyle='solid', color=colorsc[4], label=lab[4])


legend1 = ax[0].legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title='$f_{mix}$',
    ncol=1,
    loc='upper right',
    title_fontsize=10
)
ax[0].add_artist(legend1)  # Keep this legend when adding another

# Lines for second legend
from matplotlib.lines import Line2D
line11 = Line2D([], [], linestyle='solid', color=colorsc[4], linewidth=1, label="$\mathcal{P}_{cold,ic}$")
line22, = ax[0].plot(mix1_sc[0,:], mix1_sc[7,:], linestyle=':', color=colorsc[4], label='$\mathcal{P}_{cold,isn}$')
line33, = ax[0].plot(mix1_sc[0,:], mix1_sc[1,:], linestyle='--', color=colorsc[4], label='$\mathcal{P}_{cold,ih}$')

legend2 = ax[0].legend(handles=[line33, line11, line22], loc='lower right', fontsize=10)
ax[0].add_artist(legend2)

ax[0].set_ylim(0, 1)
ax[0].set_ylabel('$\mathcal{P}_{cold,ic}$')
ax[0].set_xlabel('$R$ [kpc]')

# --- Second panel ---
rr1 = np.mean(fid_M[5,:] / fid_M[2,:])
rr2 = np.mean(mix003_M[5,:] / mix003_M[2,:])
rr3 = np.mean(mix01_M[5,:] / mix01_M[2,:])
rr4 = np.mean(mix03_M[5,:] / mix03_M[2,:])
rr5 = np.mean(mix1_M[5,:] / mix1_M[2,:])
rrnd = np.mean(nd_M[5,:] / nd_M[2,:])

ax[1].plot(fid_M[0,:], fid_M[5,:], linestyle='solid', color=colorsc[0], linewidth=3, label=str(round(rr1, 2)))
ax[1].plot(mix003_M[0,:], mix003_M[5,:], linestyle='solid', color=colorsc[1], label=str(round(rr2, 2)))
ax[1].plot(mix01_M[0,:], mix01_M[5,:], linestyle='solid', color=colorsc[2], label=str(round(rr3, 2)))
ax[1].plot(mix03_M[0,:], mix03_M[5,:], linestyle='solid', color=colorsc[3], label=str(round(rr4, 2)))
ax[1].plot(mix1_M[0,:], mix1_M[5,:], linestyle='solid', color=colorsc[4], label=str(round(rr5, 2)))
#ax[1].plot(nd_M[0,:], nd_M[5,:], linestyle='--', color='k', label=str(round(rrnd, 2)))

ax[1].legend(
    loc='lower right',
    prop={'family': 'Times New Roman', 'size': 10},
    ncol=2,
    title="$\langle \dot{M}_{IC,cold} / \dot{M}_{cold} \\rangle$",
    title_fontsize=10
)
ax[1].set_ylabel('$\dot{M}_{IC,cold}$ [M$_\odot$/yr]')
ax[1].set_xlabel('$R$ [kpc]')

plt.savefig("Cfr_mix_SCALARS", bbox_inches='tight')





fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)

# --- First panel ---
ax[0].plot(fid_sc[0,:], fid_sc[16,:], linestyle='solid', color=colorsc[0], linewidth=3, label=lab[0])
ax[0].plot(mix003_sc[0,:], mix003_sc[16,:], linestyle='solid', color=colorsc[1], label=lab[1])
ax[0].plot(mix01_sc[0,:], mix01_sc[16,:], linestyle='solid', color=colorsc[2], label=lab[2])
ax[0].plot(mix03_sc[0,:], mix03_sc[16,:], linestyle='solid', color=colorsc[3], label=lab[3])
ax[0].plot(mix1_sc[0,:], mix1_sc[16,:], linestyle='solid', color=colorsc[4], label=lab[4])

legend1 = ax[0].legend(
    prop={'family': 'Times New Roman', 'size': 10},
    title='$f_{mix}$',
    ncol=1,
    loc='upper right',
    title_fontsize=10
)
ax[0].add_artist(legend1)  # Keep this legend when adding another

# Lines for second legend
from matplotlib.lines import Line2D
line11 = Line2D([], [], linestyle='solid', color=colorsc[4], linewidth=1, label="$\mathcal{P}_{cold,1}$")
line22, = ax[0].plot(mix1_sc[0,:], mix1_sc[13,:], linestyle=':', color=colorsc[4], label='$\mathcal{P}_{hot,1}$')
line33, = ax[0].plot(mix1_sc[0,:], mix1_sc[10,:], linestyle='--', color=colorsc[4], label='$\mathcal{P}_{hot,0}$')

legend2 = ax[0].legend(handles=[line33, line11, line22], loc='lower right', fontsize=10)
ax[0].add_artist(legend2)

ax[0].set_ylim(0, 1)
ax[0].set_ylabel('$\mathcal{P}_{hot,2}$')
ax[0].set_xlabel('$R$ [kpc]')

# --- Second panel ---
rr1 = np.mean(fid_M[5,:] / fid_M[2,:])
rr2 = np.mean(mix003_M[5,:] / mix003_M[2,:])
rr3 = np.mean(mix01_M[5,:] / mix01_M[2,:])
rr4 = np.mean(mix03_M[5,:] / mix03_M[2,:])
rr5 = np.mean(mix1_M[5,:] / mix1_M[2,:])

ax[1].plot(fid_M[0,:], np.abs(fid_M[6,:]), linestyle='solid', color=colorsc[0], linewidth=3, label=str(round(rr1, 2)))
ax[1].plot(mix003_M[0,:], np.abs(mix003_M[6,:]), linestyle='solid', color=colorsc[1], label=str(round(rr2, 2)))
ax[1].plot(mix01_M[0,:], np.abs(mix01_M[6,:]), linestyle='solid', color=colorsc[2], label=str(round(rr3, 2)))
ax[1].plot(mix03_M[0,:], np.abs(mix03_M[6,:]), linestyle='solid', color=colorsc[3], label=str(round(rr4, 2)))
ax[1].plot(mix1_M[0,:], np.abs(mix1_M[6,:]), linestyle='solid', color=colorsc[4], label=str(round(rr5, 2)))
ax[1].set_yscale('log')
ax[1].legend(
    loc='lower right',
    prop={'family': 'Times New Roman', 'size': 10},
    ncol=2,
    title="$\langle \dot{M}_{IC,cold} / \dot{M}_{cold} \\rangle$",
    title_fontsize=10
)
ax[1].set_ylabel('$\dot{M}_{IC,cold}$ [M$_\odot$/yr]')
ax[1].set_xlabel('$R$ [kpc]')

plt.savefig("Cfr_mix_SCALARS2", bbox_inches='tight')







fig, ax = plt.subplots(figsize=(6.3, 5), dpi=300)


rf = np.mean( np.abs(fid_M[4,:] - fid_M[3,:])/fid_M[4,:] )
r003 = np.mean( np.abs(mix003_M[4,:] - mix003_M[3,:])/mix003_M[4,:] )
r01 = np.mean( np.abs(mix01_M[4,:] - mix01_M[3,:])/mix01_M[4,:] )
r03 = np.mean( np.abs(mix03_M[4,:] - mix03_M[3,:])/mix03_M[4,:] )
r1 = np.mean( np.abs(mix1_M[4,:] - mix1_M[3,:])/mix1_M[4,:] )


# Plot all your lines and capture their handles/labels
line1, = ax.plot(fid_M[0,:], (fid_M[4,:] - fid_M[3,:])/0.3,
                 linestyle='solid', color=colorsc[0], linewidth=3, label=lab[0])
line2, = ax.plot(mix003_M[0,:], (mix003_M[4,:] - mix003_M[3,:])/0.3,
                 linestyle='solid', color=colorsc[1], label=lab[1])
line3, = ax.plot(mix01_M[0,:], (mix01_M[4,:] - mix01_M[3,:])/0.3,
                 linestyle='solid', color=colorsc[2], label=lab[2])
line4, = ax.plot(mix03_M[0,:], (mix03_M[4,:] - mix03_M[3,:])/0.3,
                 linestyle='solid', color=colorsc[3], label=lab[3])
line5, = ax.plot(mix1_M[0,:], (mix1_M[4,:] - mix1_M[3,:])/0.3,
                 linestyle='solid', color=colorsc[4], label=lab[4])

# First legend (based on plotted line labels)
legend1 = ax.legend(
    handles=[line1, line2, line3, line4, line5],
    prop={'family': 'Times New Roman', 'size': 10},
    loc='lower left',
    title=r'$f_{mix}$',
    title_fontsize=10,
    ncol=2
)
ax.add_artist(legend1)  # Add first legend manually

# Second legend: for example, different labels for the same lines
second_labels = [str(round(rf,2)), str(round(r003,2)), str(round(r01,2)), str(round(r03,2)), str(round(r1,2))]
legend2 = ax.legend(
    handles=[line1, line2, line3, line4, line5],
    labels=second_labels,
    title=r'$\langle |\dot{M}_{net}/\dot{M}_{growth} |\rangle$',
    title_fontsize=10,
    loc='lower right',
    ncol=2,
    prop={'family': 'Times New Roman', 'size': 10}
)
# No need to add_artist for second legend — last .legend() call attaches it

# Rest of your formatting
ax.set_yscale('log')
ax.set_ylabel('$\dot{M}_{net}/\Delta r$  [$M_\odot$/yr kpc$^{-1}$]')
ax.set_xlabel('$R$   [kpc]')

plt.savefig("Cfr_mix_Mnet")










fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

plt.plot(fid_T[0,fid_T[0,:]>0.5], fid_T[2,fid_T[0,:]>0.5], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
plt.plot(mix003_T[0,fid_T[0,:]>0.5], mix003_T[2,fid_T[0,:]>0.5], linestyle = 'solid', color = colorsc[1], label = lab[1])
plt.plot(mix01_T[0,fid_T[0,:]>0.5], mix01_T[2,fid_T[0,:]>0.5], linestyle = 'solid', color = colorsc[2], label = lab[2])
plt.plot(mix03_T[0,fid_T[0,:]>0.5], mix03_T[2,fid_T[0,:]>0.5], linestyle = 'solid', color = colorsc[3], label = lab[3])
plt.plot(mix1_T[0,fid_T[0,:]>0.5], mix1_T[2,fid_T[0,:]>0.5], linestyle = 'solid', color = colorsc[4])
plt.plot(fid_T[0,fid_T[0,:]>0.5], np.ones(len(fid_T[0,fid_T[0,:]>0.5])), linestyle = '--', color='k')

plt.ylabel('$\mathcal{M}_{hot}$')
plt.xlabel('$R$   [kpc]')
plt.savefig("Cfr_mix_Mach")





fig, ax = plt.subplots(figsize=(7.3, 5), dpi=300)


Rmin = 0.5
line1, = plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], -fid_sg[11,fid_sg[0,:]>Rmin]+fid_sg[5,fid_sg[0,:]>Rmin]+fid_sg[6,fid_sg[0,:]>Rmin], 
                  linestyle='solid', color=colorsc[0], linewidth=1, label=r'$\dot{u}_{tot}$')
line2, = plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], -fid_sg[11,fid_sg[0,:]>Rmin], 
                  linestyle='-.', color=colorsc[0], linewidth=2, label=r'$\dot{u}_{cool}$')
line3, = plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], fid_sg[5,fid_sg[0,:]>Rmin], 
                  linestyle='--', color=colorsc[0], linewidth=2, label=r'$\dot{u}_{mix}$')
line4, = plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], fid_sg[6,fid_sg[0,:]>Rmin], 
                  linestyle=':', color=colorsc[0], linewidth=2, label=r'$\dot{u}_{drag}$')

# Plot the mix1_sg lines (without labels on all but the desired lines)
mix_line_tot, = plt.plot(mix1_sg[0,fid_sg[0,:]>Rmin], 
                         -mix1_sg[11,fid_sg[0,:]>Rmin] + mix1_sg[5,fid_sg[0,:]>Rmin] + mix1_sg[6,fid_sg[0,:]>Rmin],
                         linestyle='solid', color=colorsc[4], linewidth=1, label=r'$1$')

mix_line2 = plt.plot(mix1_sg[0,fid_sg[0,:]>Rmin], -mix1_sg[11,fid_sg[0,:]>Rmin],
                     linestyle='-.', color=colorsc[4], linewidth=2)
mix_line3 = plt.plot(mix1_sg[0,fid_sg[0,:]>Rmin], mix1_sg[5,fid_sg[0,:]>Rmin],
                     linestyle='--', color=colorsc[4], linewidth=2)
mix_line4 = plt.plot(mix1_sg[0,fid_sg[0,:]>Rmin], mix1_sg[6,fid_sg[0,:]>Rmin],
                     linestyle=':', color=colorsc[4], linewidth=2)

# Create the first legend
legend1 = plt.legend(handles=[line1, line2, line3, line4], loc='upper right', fontsize=10, ncol=2)
plt.gca().add_artist(legend1)  # Keep this legend when adding a second

from matplotlib.lines import Line2D
line1_fmix = Line2D([], [], linestyle='solid', color=colorsc[0], linewidth=1, label=r'$0.003$')
legend2 = plt.legend(handles=[line1_fmix, mix_line_tot], loc='lower right', title=r'$f_{mix}$', fontsize=10,title_fontsize=10)





plt.yscale('symlog', linthresh=0.01)
plt.ylabel('$\dot{u}_{tot}$   [(cm/s)$^2$s$^{-1}$]')
plt.xlabel('$R$   [kpc]')

plt.savefig("Cfr_mix_utot_cfr")




fig, ax = plt.subplots(figsize=(7, 5), dpi=300)


Rmin = 0

line2, = plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], np.abs(fid_sg[12,fid_sg[0,:]>Rmin]), 
                  linestyle='-', color=colorsc[0], linewidth=1, label=r'$|\dot{\bf{P}}_{hydro}|$')
line3, = plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], np.abs(fid_sg[3,fid_sg[0,:]>Rmin]+fid_sg[4,fid_sg[0,:]>Rmin]), 
                  linestyle='--', color=colorsc[0], linewidth=1, label=r'$-|\dot{\bf{P}}_{mix}|-|\dot{\bf{P}}_{drag}|$')
#line4, = plt.plot(fid_sg[0,fid_sg[0,:]>Rmin], np.abs(fid_sg[4,fid_sg[0,:]>Rmin]), 
#                  linestyle=':', color=colorsc[0], linewidth=2, label=r'$\dot{p}_{drag}$')

# Plot the mix1_sg lines (without labels on all but the desired lines)
mix_line2, = plt.plot(mix1_sg[0,fid_sg[0,:]>Rmin], np.abs(mix1_sg[12,fid_sg[0,:]>Rmin]),
                     linestyle='-', color=colorsc[4], linewidth=1, label = lab[4])
mix_line3, = plt.plot(mix1_sg[0,fid_sg[0,:]>Rmin], np.abs(mix1_sg[3,fid_sg[0,:]>Rmin]+mix1_sg[4,fid_sg[0,:]>Rmin]),
                     linestyle='--', color=colorsc[4], linewidth=1)
#mix_line4 = plt.plot(mix1_sg[0,fid_sg[0,:]>Rmin], np.abs(mix1_sg[4,fid_sg[0,:]>Rmin]),
#                     linestyle=':', color=colorsc[4], linewidth=2)


#mix01_line2 = plt.plot(mix01_sg[0,fid_sg[0,:]>Rmin], np.abs(mix01_sg[12,fid_sg[0,:]>Rmin]),
#                     linestyle='-.', color=colorsc[2], linewidth=2, label=r'$0.1$')
mix01_line2, = plt.plot(mix01_sg[0,fid_sg[0,:]>Rmin], np.abs(mix01_sg[12,fid_sg[0,:]>Rmin]),
                     linestyle='-', color=colorsc[2], linewidth=1, label= lab[2])
mix01_line3, = plt.plot(mix01_sg[0,fid_sg[0,:]>Rmin], np.abs(mix01_sg[3,fid_sg[0,:]>Rmin]+mix01_sg[4,fid_sg[0,:]>Rmin]),
                     linestyle='--', color=colorsc[2], linewidth=1)


mix03_line2, = plt.plot(mix03_sg[0,fid_sg[0,:]>Rmin], np.abs(mix03_sg[12,fid_sg[0,:]>Rmin]),
                     linestyle='-', color=colorsc[3], linewidth=1, label=lab[3])
mix03_line3, = plt.plot(mix03_sg[0,fid_sg[0,:]>Rmin], np.abs(mix03_sg[3,fid_sg[0,:]>Rmin]+mix03_sg[4,fid_sg[0,:]>Rmin]),
                     linestyle='--', color=colorsc[3], linewidth=1)

#mix003_line2 = plt.plot(mix003_sg[0,fid_sg[0,:]>Rmin], np.abs(mix003_sg[12,fid_sg[0,:]>Rmin]),
#                     linestyle='-.', color=colorsc[2], linewidth=2, label=r'$0.1$')
mix003_line2, = plt.plot(mix003_sg[0,fid_sg[0,:]>Rmin], np.abs(mix003_sg[12,fid_sg[0,:]>Rmin]),
                     linestyle='-', color=colorsc[1], linewidth=1, label = lab[1])  
mix003_line3, = plt.plot(mix003_sg[0,fid_sg[0,:]>Rmin], np.abs(mix003_sg[3,fid_sg[0,:]>Rmin]+mix003_sg[4,fid_sg[0,:]>Rmin]),
                     linestyle='--', color=colorsc[1], linewidth=1)                     



# Create the first legend
legend1 = plt.legend(handles=[line2, line3], loc='lower right', fontsize=10, ncol=1)
plt.gca().add_artist(legend1)  # Keep this legend when adding a second

from matplotlib.lines import Line2D
line1_fmix = Line2D([], [], linestyle='solid', color=colorsc[0], linewidth=1, label=r'$0.003$')
legend2 = plt.legend(handles=[line1_fmix, mix003_line2, mix01_line2, mix03_line2, mix_line2], loc='lower left', title=r'$f_{mix}$', fontsize=10,title_fontsize=10, ncol=2)





plt.yscale('symlog', linthresh=0.01)
plt.ylabel('$\dot{p}$   [g$\cdot$cm/s]')
plt.xlabel('$R$   [kpc]')

plt.savefig("Cfr_mix_pdots_cfr")








fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'wspace': 0.4, 'hspace': 0.2}, dpi=300)

Rmin = 0
axs[0].plot(fid_sg[0,fid_sg[0,:]>Rmin], fid_sg[11,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
axs[0].plot(mix003_sg[0,fid_sg[0,:]>Rmin], mix003_sg[11,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[1], label = lab[1])
axs[0].plot(mix01_sg[0,fid_sg[0,:]>Rmin], mix01_sg[11,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[2], label = lab[2])
axs[0].plot(mix03_sg[0,fid_sg[0,:]>Rmin], mix03_sg[11,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[3], label = lab[3])
axs[0].plot(mix1_sg[0,fid_sg[0,:]>Rmin], mix1_sg[11,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[4], label = lab[4])


axs[0].plot(fid_sg[0,fid_sg[0,:]>Rmin], fid_sg[6,fid_sg[0,:]>Rmin], linestyle = '--', color = colorsc[0], linewidth = 3)
axs[0].plot(mix003_sg[0,fid_sg[0,:]>Rmin], mix003_sg[6,fid_sg[0,:]>Rmin], linestyle = '--', color = colorsc[1])
axs[0].plot(mix01_sg[0,fid_sg[0,:]>Rmin], mix01_sg[6,fid_sg[0,:]>Rmin], linestyle = '--', color = colorsc[2])
axs[0].plot(mix03_sg[0,fid_sg[0,:]>Rmin], mix03_sg[6,fid_sg[0,:]>Rmin], linestyle = '--', color = colorsc[3])
axs[0].plot(mix1_sg[0,fid_sg[0,:]>Rmin], mix1_sg[6,fid_sg[0,:]>Rmin], linestyle = '--', color = colorsc[4])


axs[0].plot(fid_sg[0,fid_sg[0,:]>Rmin], fid_sg[5,fid_sg[0,:]>Rmin], linestyle = ':', color = colorsc[0], linewidth = 3)
axs[0].plot(mix003_sg[0,fid_sg[0,:]>Rmin], mix003_sg[5,fid_sg[0,:]>Rmin], linestyle = ':', color = colorsc[1])
axs[0].plot(mix01_sg[0,fid_sg[0,:]>Rmin], mix01_sg[5,fid_sg[0,:]>Rmin], linestyle = ':', color = colorsc[2])
axs[0].plot(mix03_sg[0,fid_sg[0,:]>Rmin], mix03_sg[5,fid_sg[0,:]>Rmin], linestyle = ':', color = colorsc[3])
axs[0].plot(mix1_sg[0,fid_sg[0,:]>Rmin], mix1_sg[5,fid_sg[0,:]>Rmin], linestyle = ':', color = colorsc[4])

axs[0].legend(
    title=r'$f_{mix}$',
    loc='upper right',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=2,
    title_fontsize=10,
)
axs[0].set_yscale('symlog', linthresh=0.001)
axs[0].set_ylabel('$ucool$   [(cm/s)^2/s]')
axs[0].set_xlabel('$R$   [kpc]')
#plt.show()


    ## - Cold

axs[1].plot(fid_sg[0,fid_sg[0,:]>Rmin], -fid_sg[11,fid_sg[0,:]>Rmin]+fid_sg[5,fid_sg[0,:]>Rmin]+fid_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[0], linewidth = 3, label = lab[0])
axs[1].plot(mix003_sg[0,fid_sg[0,:]>Rmin], -mix003_sg[11,fid_sg[0,:]>Rmin]+mix003_sg[5,fid_sg[0,:]>Rmin]+mix003_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[1], label = lab[1])
axs[1].plot(mix01_sg[0,fid_sg[0,:]>Rmin], -mix01_sg[11,fid_sg[0,:]>Rmin]+mix01_sg[5,fid_sg[0,:]>Rmin]+mix01_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[2], label = lab[2])
axs[1].plot(mix03_sg[0,fid_sg[0,:]>Rmin], -mix03_sg[11,fid_sg[0,:]>Rmin]+mix03_sg[5,fid_sg[0,:]>Rmin]+mix03_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[3], label = lab[3])
axs[1].plot(mix1_sg[0,fid_sg[0,:]>Rmin], -mix1_sg[11,fid_sg[0,:]>Rmin]+mix1_sg[5,fid_sg[0,:]>Rmin]+mix1_sg[6,fid_sg[0,:]>Rmin], linestyle = 'solid', color = colorsc[4], label = lab[4])


#axs[1].plot(fid_sg[0,:], fid_sg[6,:], linestyle = '--', color = colors[0], linewidth = 3, label = lab[0])
#axs[1].plot(mix003_sg[0,:], mix003_sg[6,:], linestyle = '--', color = colors[2], label = lab[2])
#axs[1].plot(mix01_sg[0,:], mix01_sg[6,:], linestyle = '--', color = colors[6], label = lab[6])
#axs[1].plot(mix03_sg[0,:], mix03_sg[6,:], linestyle = '--', color = colors[7], label = lab[7])

axs[1].set_yscale('symlog', linthresh=0.01)
axs[1].set_ylabel('$\dot{u}_{tot}$   [(cm/s)$^2$s$^{-1}$]')
axs[1].set_xlabel('$R$   [kpc]')


#axs[0].text(0.5, 1.05, 'Temperature', transform=axs[0].transAxes,
#            ha='center', va='bottom', fontsize=12, fontweight='bold',
#            fontname='Times New Roman')
#axs[1].text(0.5, 1.05
#, 'Hot phase heating due to source terms', transform=axs[1].transAxes,
#            ha='center', va='bottom', fontsize=12, fontweight='bold',
#            fontname='Times New Roman')


#plt.show()
plt.savefig("Cfr_mix_ucool")
















################### ===== Mdots ====== ########################


fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)


axs[0].plot(fid_M[0,:], fid_M[1,:], linestyle = 'solid', color = colors[0], linewidth = 3, label = lab[0])
axs[0].plot(mix001_M[0,:], mix001_M[1,:], linestyle = 'solid', color = colors[1], label = lab[1])
axs[0].plot(mix003_M[0,:], mix003_M[1,:], linestyle = 'solid', color = colors[2], label = lab[2])
axs[0].plot(mix005_M[0,:], mix005_M[1,:], linestyle = 'solid', color = colors[3], label = lab[3])
axs[0].plot(mix006_M[0,:], mix006_M[1,:], linestyle = 'solid', color = colors[4], label = lab[4])
axs[0].plot(mix008_M[0,:], mix008_M[1,:], linestyle = 'solid', color = colors[5], label = lab[5])
axs[0].plot(mix01_M[0,:], mix01_M[1,:], linestyle = 'solid', color = colors[6], label = lab[6])
axs[0].plot(mix03_M[0,:], mix03_M[1,:], linestyle = 'solid', color = colors[7], label = lab[7])
axs[0].plot(mix003nd_M[0,:], mix003nd_M[1,:], linestyle = '--', color = colors[2], label = lab[8])


#axs[0].legend()
axs[0].set_ylabel('$\dot{M}$   [$M_\odot$/yr]')
axs[0].set_xlabel('$R$   [kpc]')
#plt.show()


    ## - Cold

axs[1].plot(fid_M[0,:], fid_M[2,:], linestyle = 'solid', color = colors[0], linewidth = 3, label = lab[0])
axs[1].plot(mix001_M[0,:], mix001_M[2,:], linestyle = 'solid', color = colors[1], label = lab[1])
axs[1].plot(mix003_M[0,:], mix003_M[2,:], linestyle = 'solid', color = colors[2], label = lab[2])
axs[1].plot(mix005_M[0,:], mix005_M[2,:], linestyle = 'solid', color = colors[3], label = lab[3])
axs[1].plot(mix006_M[0,:], mix006_M[2,:], linestyle = 'solid', color = colors[4], label = lab[4])
axs[1].plot(mix008_M[0,:], mix008_M[2,:], linestyle = 'solid', color = colors[5], label = lab[5])
axs[1].plot(mix01_M[0,:], mix01_M[2,:], linestyle = 'solid', color = colors[6], label = lab[6])
axs[1].plot(mix03_M[0,:], mix03_M[2,:], linestyle = 'solid', color = colors[7], label = lab[7])
axs[1].plot(mix003nd_M[0,:], mix003nd_M[2,:], linestyle = '--', color = colors[2], label = lab[8])


axs[1].legend()


axs[1].set_ylabel('$\dot{M}$   [$M_\odot$/yr]')
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Hot phase', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Cold phase', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')



#plt.show()
plt.savefig("Cfr_mix_Mdots")








################### ===== umix ====== ########################


fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)


axs[0].plot(fid_T[0,:], fid_T[1,:], linestyle = 'solid', color = colors[0], linewidth = 3, label = lab[0])
axs[0].plot(mix001_T[0,:], mix001_T[1,:], linestyle = 'solid', color = colors[1], label = lab[1])
axs[0].plot(mix003_T[0,:], mix003_T[1,:], linestyle = 'solid', color = colors[2], label = lab[2])
axs[0].plot(mix005_T[0,:], mix005_T[1,:], linestyle = 'solid', color = colors[3], label = lab[3])
axs[0].plot(mix006_T[0,:], mix006_T[1,:], linestyle = 'solid', color = colors[4], label = lab[4])
axs[0].plot(mix008_T[0,:], mix008_T[1,:], linestyle = 'solid', color = colors[5], label = lab[5])
axs[0].plot(mix01_T[0,:], mix01_T[1,:], linestyle = 'solid', color = colors[6], label = lab[6])
axs[0].plot(mix03_T[0,:], mix03_T[1,:], linestyle = 'solid', color = colors[7], label = lab[7])
axs[0].plot(mix003nd_T[0,:], mix003nd_T[1,:], linestyle = '--', color = colors[2], label = lab[8])


#axs[0].legend()
axs[0].set_yscale('log')
axs[0].set_ylabel('$T$   [K]')
axs[0].set_xlabel('$R$   [kpc]')
#plt.show()


    ## - Cold

axs[1].plot(fid_sg[0,:], fid_sg[5,:], linestyle = 'solid', color = colors[0], linewidth = 3, label = lab[0])
#axs[1].plot(mix001_sg[0,:], mix001_sg[5,:], linestyle = 'solid', color = colors[1], label = lab[1])
axs[1].plot(mix003_sg[0,:], mix003_sg[5,:], linestyle = 'solid', color = colors[2], label = lab[2])
#axs[1].plot(mix005_sg[0,:], mix005_sg[5,:], linestyle = 'solid', color = colors[3], label = lab[3])
#axs[1].plot(mix006_sg[0,:], mix006_sg[5,:], linestyle = 'solid', color = colors[4], label = lab[4])
#axs[1].plot(mix008_sg[0,:], mix008_sg[5,:], linestyle = 'solid', color = colors[5], label = lab[5])
axs[1].plot(mix01_sg[0,:], mix01_sg[5,:], linestyle = 'solid', color = colors[6], label = lab[6])
axs[1].plot(mix03_sg[0,:], mix03_sg[5,:], linestyle = 'solid', color = colors[7], label = lab[7])
#axs[1].plot(mix003nd_sg[0,:], mix003nd_sg[5,:], linestyle = '--', color = colors[2], label = lab[8])


axs[1].plot(fid_sg[0,:], fid_sg[6,:], linestyle = '--', color = colors[0], linewidth = 3, label = lab[0])
#axs[1].plot(mix001_sg[0,:], mix001_sg[6,:], linestyle = '--', color = colors[1], label = lab[1])
axs[1].plot(mix003_sg[0,:], mix003_sg[6,:], linestyle = '--', color = colors[2], label = lab[2])
#axs[1].plot(mix005_sg[0,:], mix005_sg[6,:], linestyle = '--', color = colors[3], label = lab[3])
#axs[1].plot(mix006_sg[0,:], mix006_sg[6,:], linestyle = '--', color = colors[4], label = lab[4])
#axs[1].plot(mix008_sg[0,:], mix008_sg[6,:], linestyle = '--', color = colors[5], label = lab[5])
axs[1].plot(mix01_sg[0,:], mix01_sg[6,:], linestyle = '--', color = colors[6], label = lab[6])
axs[1].plot(mix03_sg[0,:], mix03_sg[6,:], linestyle = '--', color = colors[7], label = lab[7])
#axs[1].plot(mix003nd_sg[0,:], mix003nd_sg[6,:], linestyle = ':', color = colors[2], label = lab[8])


axs[1].set_yscale('symlog', linthresh=1e4)
axs[1].set_ylabel('$\dot{u}$   [?]')
axs[1].set_xlabel('$R$   [kpc]')


axs[0].text(0.5, 1.05, 'Temperature', transform=axs[0].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')
axs[1].text(0.5, 1.05
, 'Hot phase heating due to source terms', transform=axs[1].transAxes,
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            fontname='Times New Roman')



#plt.show()
plt.savefig("Cfr_mix_Tudot")





################### ===== Ml/Mdots ====== ########################


fig, ax = plt.subplots()

plt.plot(fid_M[0,:], fid_M[1,:]/fid_M[3,:], linestyle = 'solid', color = colors[0], linewidth = 3, label = lab[0])
plt.plot(mix001_M[0,:], mix001_M[1,:]/mix001_M[3,:], linestyle = 'solid', color = colors[1], label = lab[1])
plt.plot(mix003_M[0,:], mix003_M[1,:]/mix003_M[3,:], linestyle = 'solid', color = colors[2], label = lab[2])
plt.plot(mix005_M[0,:], mix005_M[1,:]/mix005_M[3,:], linestyle = 'solid', color = colors[3], label = lab[3])
plt.plot(mix006_M[0,:], mix006_M[1,:]/mix006_M[3,:], linestyle = 'solid', color = colors[4], label = lab[4])
plt.plot(mix008_M[0,:], mix008_M[1,:]/mix008_M[3,:], linestyle = 'solid', color = colors[5], label = lab[5])
plt.plot(mix01_M[0,:], mix01_M[1,:]/mix01_M[3,:], linestyle = 'solid', color = colors[6], label = lab[6])
plt.plot(mix03_M[0,:], mix03_M[1,:]/mix03_M[3,:], linestyle = 'solid', color = colors[7], label = lab[7])
plt.plot(mix003nd_M[0,:], mix003nd_M[1,:]/mix003nd_M[3,:], linestyle = '--', color = colors[2], label = lab[8])


#axs[0].legend()
plt.yscale('log')
plt.ylabel('$Ml/Mh$')
plt.xlabel('$R$   [kpc]')

#plt.show()
plt.savefig("Cfr_mix_mlmh")








################### ===== Scalars ====== ########################


fig, ax = plt.subplots()

plt.plot(fid_sc[0,:], fid_sc[1,:], linestyle = 'solid', color = colors[0], linewidth = 3, label = lab[0])
plt.plot(mix001_sc[0,:], mix001_sc[1,:], linestyle = 'solid', color = colors[1], label = lab[1])
plt.plot(mix003_sc[0,:], mix003_sc[1,:], linestyle = 'solid', color = colors[2], label = lab[2])
plt.plot(mix01_sc[0,:], mix01_sc[1,:], linestyle = 'solid', color = colors[6], label = lab[6])
plt.plot(mix03_sc[0,:], mix03_sc[1,:], linestyle = 'solid', color = colors[7], label = lab[7])


'''
plt.plot(fid_sc[0,:], fid_sc[4,:], linestyle = '--', color = colors[0], linewidth = 3, label = lab[0])
plt.plot(mix001_sc[0,:], mix001_sc[4,:], linestyle = '--', color = colors[1], label = lab[1])
plt.plot(mix003_sc[0,:], mix003_sc[4,:], linestyle = '--', color = colors[2], label = lab[2])
plt.plot(mix005_sc[0,:], mix005_sc[4,:], linestyle = '--', color = colors[3], label = lab[3])
plt.plot(mix006_sc[0,:], mix006_sc[4,:], linestyle = '--', color = colors[4], label = lab[4])
plt.plot(mix008_sc[0,:], mix008_sc[4,:], linestyle = '--', color = colors[5], label = lab[5])
plt.plot(mix01_sc[0,:], mix01_sc[4,:], linestyle = '--', color = colors[6], label = lab[6])
plt.plot(mix03_sc[0,:], mix03_sc[4,:], linestyle = '--', color = colors[7], label = lab[7])
plt.plot(mix003nd_sc[0,:], mix003nd_sc[4,:], linestyle = '--', color = colors[2], label = lab[8], linewidth = 5)#
'''

plt.plot(fid_sc[0,:], fid_sc[7,:], linestyle = ':', color = colors[0], linewidth = 3, label = lab[0])
plt.plot(mix001_sc[0,:], mix001_sc[7,:], linestyle = ':', color = colors[1], label = lab[1])
plt.plot(mix003_sc[0,:], mix003_sc[7,:], linestyle = ':', color = colors[2], label = lab[2])
plt.plot(mix01_sc[0,:], mix01_sc[7,:], linestyle = ':', color = colors[6], label = lab[6])
plt.plot(mix03_sc[0,:], mix03_sc[7,:], linestyle = ':', color = colors[7], label = lab[7])



#axs[0].legend()
plt.ylabel('$\mathcal{P}_{cold}$')
plt.xlabel('$R$   [kpc]')

#plt.show()
plt.savefig("Cfr_mix_scalarscold")
sys.exit(-1)








fig, ax  = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.3}, dpi=300)

ax[0].plot(scalars[0], scalars[1], label="IC hot", color = 'r')
ax[0].fill_between(scalars[0], scalars[2], scalars[3], alpha=0.2, color='r')

ax[0].plot(scalars[0], scalars[4], label="IC cold", color = 'blue')
ax[0].fill_between(scalars[0], scalars[5], scalars[6], alpha=0.2, color='blue')

ax[0].plot(scalars[0], scalars[7], label="SNe", color = 'k')
ax[0].fill_between(scalars[0], scalars[8], scalars[9], alpha=0.2, color='k')
ax[0].set_xlabel("$R$   [kpc]")
ax[0].set_ylabel("$\mathcal{P}_{cold}$")
ax[0].legend()

ax[1].plot(scalars[0], scalars[10], label="scalar0", color = 'r')
ax[1].fill_between(scalars[0], scalars[11], scalars[12], alpha=0.2, color='r')

ax[1].plot(scalars[0], scalars[13], label="scalar1", color = 'blue')
ax[1].fill_between(scalars[0], scalars[14], scalars[15], alpha=0.2, color='blue')

ax[1].plot(scalars[0], scalars[16], label="scalar2", color = 'k')
ax[1].fill_between(scalars[0], scalars[17], scalars[18], alpha=0.2, color='k')
ax[1].set_ylabel("$\mathcal{P}_{hot}$")

ax[1].set_xlabel("$R$   [kpc]")

plt.savefig("ScalarPlot")











































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


fig, ax = plt.subplots()
plt.plot(fid_sg[0,:], fid_sg[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

plt.plot(fcl_sg[0,:], fcl_sg[1,:], linestyle = ls[0], color = col[0], label = lab[0])
plt.plot(fch_sg[0,:], fch_sg[1,:], linestyle = ls[1], color = col[0], label = lab[1])


plt.plot(mixl_sg[0,:], mixl_sg[1,:], linestyle = ls[1], color = col[1], label = lab[2])
plt.plot(mixh_sg[0,:], mixh_sg[1,:], linestyle = ls[2], color = col[1], label = lab[3])

plt.plot(rcl001_sg[0,:], rcl001_sg[1,:], linestyle = ls[0], color = col[2], label = lab[4])
plt.plot(rcl06_sg[0,:], rcl06_sg[1,:], linestyle = ls[1], color = col[2], label = lab[5])
#plt.plot(rcl02_sg[0,:], rcl02_sg[1,:], linestyle = ls[3], color = col[2], label = lab[8])

plt.legend()
plt.ylim(1e-4,1)
plt.yscale('log')
plt.xlabel('$R$   [kpc]')
plt.ylabel("$r_{crit}$  [$M_\odot$]")
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





r_fid = np.mean(fid_sg[4,:]/fid_sg[3,:])
r_fcl = np.mean(fcl_sg[4,:]/fcl_sg[3,:])
r_fch = np.mean(fch_sg[4,:]/fch_sg[3,:])
r_mixl = np.mean(mixl_sg[4,:]/mixl_sg[3,:])
r_mixh = np.nanmean(mixh_sg[4,:]/mixh_sg[3,:])
r_rcl001 = np.mean(rcl001_sg[4,:]/rcl001_sg[3,:])
r_rcl06 = np.mean(rcl06_sg[4,:]/rcl06_sg[3,:])


print(mixh_sg[4,:]/mixh_sg[3,:])




fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'wspace': 0.3, 'hspace': 0.2}, dpi=300)

axs[0].plot(fid_v[0,:], fid_v[1,:], linestyle = 'solid', color = 'k', linewidth = 3, label = str(round(r_fid,2)))


axs[0].plot(fcl_v[0,:], fcl_v[1,:], linestyle = ls[0], color = col[0], label = str(round(r_fcl,2)))
axs[0].plot(fch_v[0,:], fch_v[1,:], linestyle = ls[1], color = col[0], label = str(round(r_fch,2)))

#axs[0].plot(NoMix_v[0,:], NoMix_v[1,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[0].plot(mixl_v[0,:], mixl_v[1,:], linestyle = ls[1], color = col[1], label = str(round(r_mixl,2)))
axs[0].plot(mixh_v[0,:], mixh_v[1,:], linestyle = ls[2], color = col[1], label = str(round(r_mixh,2)))

axs[0].plot(rcl001_v[0,:], rcl001_v[1,:], linestyle = ls[0], color = col[2], label = str(round(r_rcl001,2)))
#axs[0].plot(rcl002_v[0,:], rcl002_v[1,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[0].plot(rcl06_v[0,:], rcl06_v[1,:], linestyle = ls[2], color = col[2], label = str(round(r_rcl06,2)))
#axs[0,0].plot(rcl02_v[0,:], rcl02_v[1,:], linestyle = ls[3], color = col[2], label = lab[8])
axs[0].legend(
    title=r'$\langle \dot{p}_{\mathrm{drag}} / \dot{p}_{\mathrm{mix}}\rangle$',
    loc='lower right',  # You can change to 'upper right', etc.
    fontsize=10,
    ncol=3,
    title_fontsize=10,
)

#axs[0].legend()
axs[0].set_xlabel('$R$   [kpc]')
axs[0].set_ylim(250,1750)
#plt.yscale('log')
#axs[0,0].set_xlabel('$R$   [kpc]')
axs[0].set_ylabel("$v$  [km/s]")
#plt.show()


    ## - Cold

axs[1].plot(fid_v[0,:], fid_v[2,:], linestyle = 'solid', color = 'k', linewidth = 3, label = 'fid')

axs[1].plot(fcl_v[0,:], fcl_v[2,:], linestyle = ls[0], color = col[0], label = lab[0])
axs[1].plot(fch_v[0,:], fch_v[2,:], linestyle = ls[1], color = col[0], label = lab[1])

#axs[1].plot(NoMix_v[0,:], NoMix_v[2,:], linestyle = ls[0], color = col[1], label = lab[2])
axs[1].plot(mixl_v[0,:], mixl_v[2,:], linestyle = ls[1], color = col[1], label = lab[2])
axs[1].plot(mixh_v[0,:], mixh_v[2,:], linestyle = ls[2], color = col[1], label = lab[3])

axs[1].plot(rcl001_v[0,:], rcl001_v[2,:], linestyle = ls[0], color = col[2], label = lab[4])
#axs[1].plot(rcl002_v[0,:], rcl002_v[2,:], linestyle = ls[1], color = col[2], label = lab[6])
axs[1].plot(rcl06_v[0,:], rcl06_v[2,:], linestyle = ls[2], color = col[2], label = lab[5])
#axs[0,1].plot(rcl02_v[0,:], rcl02_v[2,:], linestyle = ls[3], color = col[2], label = lab[8])
axs[1].legend()


axs[1].set_ylim(100,1200)
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
    title_fontsize=10,
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
    title_fontsize=10,
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