import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import rcParams

rcParams.update({'font.size': 24})
names = ['16','34','35','38','47','64','413','453','1647']
M_J = 9.54e-4
mu = [0.20255 /(0.20255+0.68970),1.0208/(1.0208+1.0479),0.8094/(0.8094+0.8877),0.249/(0.249+0.949),0.342/(0.342+0.957),0.408/(0.408+1.528),0.5423/(0.5423+0.82),0.1951/(0.1951+0.944),0.9678/(0.9678+1.2207)]
M_T = [(0.20255+0.68970),(1.0208+1.0479),(0.8094+0.8877),(0.249+0.949),(0.342+0.957),(0.408+1.528),(0.5423+0.82),(0.1951+0.944),(0.9678+1.2207)]
e_b = [0.15944,0.52087,0.1421,0.1032,0.0288,0.2117,0.0365,0.0524,0.1602]
a_b = [0.22431,0.22882,0.17617,0.1649,0.08145,0.1744,0.10148,0.185319,0.1276]
m_p = [0.333*M_J,0.22*M_J,0.127*M_J,0.38*M_J,0.1*M_J,0.211*M_J,0.211*M_J,0.03*M_J,1.52*M_J]
cmap = cm.nipy_spectral

vmin = 1.3
vmax = 4.5
my_cmap=cm.get_cmap(cmap)
norm = colors.Normalize(vmin,vmax)
cmmapable =cm.ScalarMappable(norm,my_cmap)
cmmapable.set_array(range(0,1))
cmap.set_under('gray')

fs = 'x-large'
width = 10.
aspect = 16./9.
ms = 5

data = np.genfromtxt("../a_crit.txt",delimiter=',',comments='#')

X = data[:,0]
X[X==0.001] = 0.
Y = data[:,1]
Z = data[:,2]

xi = np.linspace(0,0.5,51)
yi = np.linspace(0,0.8,81)
zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'linear',fill_value=0)

bins = 300


fig = plt.figure(1,figsize = (aspect*width,width))
ax = fig.add_subplot(111)

ax.pcolormesh(xi,yi,zi,cmap = cmap,vmin=vmin,vmax=vmax)
for i in xrange(0,len(mu)):
	ax.plot(mu[i],e_b[i],'.',color='w',ms=10)
	ax.text(mu[i],e_b[i]+0.02,names[i],color='w',fontsize='medium',horizontalalignment='center')

ax.set_xlabel("$\mu$",fontsize=fs)
ax.set_ylabel("$e_{bin}$",fontsize=fs)
ax.set_ylim(0,0.8)
ax.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)

color_label='$a_{c}$'
cax = fig.add_axes([0.92,0.11,0.015,0.77])
cbar= plt.colorbar(cmmapable,cax=cax,orientation='vertical',ticks=np.arange(1.5,4.5+0.5,0.5))
cbar.set_label(color_label,fontsize=fs)
cax.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)
cmmapable.set_clim(vmin,vmax)



fig.savefig("Q18_Fig5.png" , dpi = 300, bbox_inches = 'tight')
plt.close()
