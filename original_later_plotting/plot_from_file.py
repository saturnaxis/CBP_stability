import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import rcParams

rcParams.update({'font.size': 24})
cmap = cm.gist_rainbow

aspect = 1.
width = 16.
ms = 5
lw = 2.
fs = 'x-large'

vmin_t = 0
vmax_t = 5

my_cmap_t=cm.get_cmap(cmap)
norm_t = colors.Normalize(vmin_t,vmax_t)
cmmapable_t =cm.ScalarMappable(norm_t,my_cmap_t)
cmmapable_t.set_array(range(0,1))
cmap.set_under('w')

vmin_e = 0
vmax_e = 1
my_cmap_e=cm.get_cmap(cmap)
norm_e = colors.Normalize(vmin_e,vmax_e)
cmmapable_e =cm.ScalarMappable(norm_e,my_cmap_e)
cmmapable_e.set_array(range(0,1))


fig = plt.figure(figsize=(aspect*width,2*width),dpi=300)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

bins = 100
xi = np.arange(0.,181.,1.)
yi = np.arange(1.,5.01,0.01)

data = np.genfromtxt("MaxEcc_[0.200,0.200].txt",delimiter=',',comments='#')
zi = griddata((data[:,1],data[:,0]),np.log10(np.abs(data[:,-1])),(xi[None,:],yi[:,None]),method = 'linear')
ax1.pcolormesh(xi-0.5,yi-0.05,zi,cmap = cmap,vmin=vmin_t,vmax=vmax_t)

ax1.set_xticks([f for f in np.arange(0.,210.,30.)])
ax1.set_ylabel("$a_p/a_{bin}$",fontsize=fs)
ax1.set_ylim(1,5)
ax1.set_xlim(0,180)
ax1.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)

color_label='$\log_{10}$ $t$(yr)'
cax_t = fig.add_axes([0.93,0.52,0.015,0.35])
cbar_t = plt.colorbar(cmmapable_t,cax=cax_t,orientation='vertical')
cbar_t.set_label(color_label,fontsize=fs)
cbar_t.ax.tick_params(axis='both', direction='out',length = 8.0, width = 8.0,labelsize=fs)


delta_e = data[:,2] - data[:,3]
zi = griddata((data[:,1],data[:,0]),(delta_e),(xi[None,:],yi[:,None]),method = 'linear')
ax2.pcolormesh(xi-0.5,yi-0.05,zi,cmap = cmap,vmin=vmin_e,vmax=vmax_e)

ax2.set_xticks([f for f in np.arange(0.,210.,30.)])
ax2.set_ylabel("$a_p/a_{bin}$",fontsize=fs)
ax2.set_xlabel("Mean Anomaly (deg.)",fontsize=fs)
ax2.set_ylim(1,5)
ax2.set_xlim(0,180)
ax2.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)

color_label='$\Delta e$'
cax_e = fig.add_axes([0.93,0.12,0.015,0.35])
cbar_e = plt.colorbar(cmmapable_e,cax=cax_e,orientation='vertical')
cbar_e.set_label(color_label,fontsize=fs)
cbar_e.ax.tick_params(axis='both', direction='out',length = 8.0, width = 8.0,labelsize=fs)

fig.subplots_adjust(hspace=0.1)

fig.savefig("CBP_map[0.2,0.2].png"  , dpi = 600, bbox_inches = 'tight')
plt.close()
