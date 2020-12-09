import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.cm as cm
import sys
from matplotlib import rcParams
import os

rcParams.update({'font.size': 22})

def func(x, a,b,c,d,e,f,g):
    return  a + b*x[1] + c*x[1]**2 + d*x[0] + e*x[0]*x[1] + f*x[0]**2 + g*(x[0]*x[1])**2

def func2(x, a,b,c,d,e,f,g):
    return  a + b*x[1] + c*x[1]**2 + d*(x[0])**(1./3.) + e*(x[0])**(1./3.)*x[1] + f*(x[0])**(2./3.) + g*(x[0])**(2./3.)*x[1]**2

cmap = cm.gnuplot_r
vmin = -3.
vmax = 0.
my_cmap=cm.get_cmap(cmap)
norm = colors.Normalize(vmin,vmax)
cmmapable =cm.ScalarMappable(norm,my_cmap)
cmmapable.set_array(range(0,1))
cmap.set_under('gray')
cmap.set_over('white')

fn = sys.argv[1]

Kep_names = ['16','34','35','38','47','64','413','453']
Kep_idx = Kep_names.index(fn)
Kep_ae = np.zeros((8,2))
Kep_ae[:,0] = [0.70,1.09,0.6,0.46,0.27,0.63,0.36,0.79]
Kep_ae[:,1] = [0.0069,0.18,0.042,0.032,0.035,0.054,0.12,0.038]

mu = [0.20255 /(0.20255+0.68970),1.0208/(1.0208+1.0479),0.8094/(0.8094+0.8877),0.249/(0.249+0.949),0.342/(0.342+0.957),0.408/(0.408+1.528),0.5423/(0.5423+0.82),0.1951/(0.1951+0.944),0.9678/(0.9678+1.2207)]
M_T = [(0.20255+0.68970),(1.0208+1.0479),(0.8094+0.8877),(0.249+0.949),(0.342+0.957),(0.408+1.528),(0.5423+0.82),(0.1951+0.944),(0.9678+1.2207)]
e_b = [0.15944,0.52087,0.1421,0.1032,0.0288,0.2117,0.0365,0.0524,0.1602]
a_b = [0.22431,0.22882,0.17617,0.1649,0.08145,0.1744,0.10148,0.185319,0.1276]

Qa = [0.60,0.81,0.48,0.43,0.18,0.52,0.24,0.42]
popt = [1.60,5.10,-2.22,4.12,-4.27,-5.09,4.61]
poptQ = [1.48,3.92,-1.41,5.14,0.33,-7.95,-4.89]
poptQ2 = [0.93,2.67,-0.25,3.72,2.25,-2.72,-4.17]

fs = 'x-large'
width = 10.
aspect = 16./9.
ms = 5
lw = 4

fig = plt.figure(1,figsize = (aspect*width,width))
ax = fig.add_subplot(111)


home = os.getcwd() + "/"

tscale = 1e5
data = np.genfromtxt(home+"%s_out/MaxEcc_%s.txt" % (fn,fn) ,delimiter = ',',comments = '#')
unstab = np.where(data[:,-1]<tscale)[0]
stab = np.where(data[:,-1]==tscale)[0]
data[unstab,2] = 3.
data[unstab,3] = 3.


X = data[:,0]
Y = data[:,1]
astab = np.where(np.logical_and(Y==0,data[:,-1]==tscale))[0]
astab = np.min(X[astab])

Z = (data[:,2]-data[:,3])
Z[unstab] = 3.


xmin,xmax = np.min(X),np.max(X)
nx = int((xmax-xmin)/0.001)+1

xi = np.linspace(xmin,xmax,nx)
yi = np.linspace(0.,0.5,51)
zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'nearest')

cs = ax.pcolormesh(xi,yi,np.log10(zi),cmap=cmap,vmin=vmin,vmax=vmax,zorder=2)
fact = (1.-0.2)
astab = func([mu[Kep_idx],e_b[Kep_idx]], *poptQ)
ax.plot(xi,fact*(1.-astab*a_b[Kep_idx]/xi),'c-',lw=lw,label='Fit 1',zorder=4)

ax.plot(xi,fact*(1.-Qa[Kep_idx]/xi),'y-.',lw=lw,label='Interpolation',zorder=4)

astab = func([mu[Kep_idx],e_b[Kep_idx]], *popt)
ax.plot(xi,fact*(1.-astab*a_b[Kep_idx]/xi),'--',color='violet',lw=lw,label='HW99',zorder=4)
ax.plot(Kep_ae[Kep_idx,0],Kep_ae[Kep_idx,1],'.',color='limegreen',ms=20,zorder=5)


ax.set_ylim([0.,0.5])
ax.set_xlim(0.15,1.5)
ax.set_xticks([0.2 + k*0.1 for k in range(0,14)])

ax.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)
ax.set_yticks([k*0.1 for k in range(0,6)])
#ax.set_yticklabels(['%1.1f' % (k*0.1) for k in range(0,5)])

ax.set_ylabel("$e_p$",fontsize = fs)
ax.set_xlabel("$a_p$ (AU)", fontsize = fs)


ax.legend(loc='upper right',fontsize='large',markerscale=3,framealpha=0,ncol=3,bbox_to_anchor=(0.75,1.1),handletextpad=0.25)

color_label='$\log_{10}[\Delta e]$'
cax = fig.add_axes([0.92,0.11,0.015,0.77])
cbar = plt.colorbar(cmmapable,cax=cax,orientation='vertical')#,ticks=np.linspace(vmin,vmax,11)) #draw colorbar  ,ticks =[0.0,0.2,0.4,0.6,0.8,1.0]
cbar.set_label(color_label,fontsize=fs)
cmmapable.set_clim(vmin,vmax)
cax.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)


fig.savefig(home+"Kep_%s.png" % fn, dpi = 300, bbox_inches = 'tight')
plt.close()
