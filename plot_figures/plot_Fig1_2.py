import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import griddata
import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.cm as cm
import sys
from matplotlib import rcParams

rcParams.update({'font.size': 22})

fig_num = int(sys.argv[1])

cmap = cm.gnuplot_r
vmin = -3.
vmax = 0.
my_cmap=cm.get_cmap(cmap)
norm = colors.Normalize(vmin,vmax)
cmmapable =cm.ScalarMappable(norm,my_cmap)
cmmapable.set_array(range(0,1))
cmap.set_under('gray')
cmap.set_over('white')

fs = 'x-large'
width = 13.
aspect = 1.
tscale = 1e5

val = [0.0,0.1,0.3,0.5]
home = os.getcwd() + "/"
datalist = [f for f in os.listdir('Fig%i_data' % fig_num) if f.endswith('.txt')]

f, axarr = plt.subplots(4,4)
f.set_size_inches(aspect*width,width)

for dat in datalist:

	ax_mu = val.index(float(dat[8:11]))
	ax_ecc = val.index(float(dat[12:15]))

	ax = axarr[ax_mu,ax_ecc]

	data = np.genfromtxt(home+"Fig1_data/"+dat,delimiter = ',',comments = '#')
	unstab = np.where(data[:,-1]<tscale)[0]
	stab = np.where(data[:,-1]==tscale)[0]

	X = data[:,1]
	Y = data[:,0]
	idx_cr = np.where(np.logical_and(X==90.,data[:,-1]==tscale))[0]

	Z = data[:,2]
	Z[unstab] = 1.5
	y_cr = []
	for i in xrange(0,401):
		yi = 1.01 + i*0.01
		row_x = np.where(np.abs(Y-yi)<1e-6)[0]
		all_stab = np.where(data[row_x,-1]==tscale)[0]
		if len(all_stab)>180:
			y_cr.append(yi)


	xi = np.linspace(0,360,181)
	yi = np.linspace(1.,5.,401)
	zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'nearest')

	bins =200 

	cs = ax.contourf(xi,yi,np.log10(zi),bins, cmap = cmap,vmin=vmin,vmax=vmax)#,locator=ticker.LogLocator() ,norm=colors.LogNorm())

	if len(y_cr)>0:
		a_crit = np.min(y_cr)
		ax.axhline(a_crit,linestyle='-',color='c',lw=2)	
	
	ax.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)
	ax.set_xticks(np.arange(0,360+90,90))
	yticks = np.arange(1.,6.0,1.)

	ax.set_yticks(yticks)
	ax.set_ylim(1,5.)
	mu = val[ax_mu]
	eb = val[ax_ecc]
	if mu == 0:
		ax.text(1.0,0.9,"(%1.3f,%1.1f)" % (0.001,eb), color='k',fontsize='medium',horizontalalignment='right',weight='bold',transform=ax.transAxes)
	else:
		ax.text(1.0,0.9,"(%1.1f,%1.1f)" % (mu,eb), color='k',fontsize='medium',horizontalalignment='right',weight='bold',transform=ax.transAxes)
	ax.text(0.0,0.02,"$a_{c}$ = %1.2f" % a_crit, color='c',fontsize='small',horizontalalignment='left',weight='bold',transform=ax.transAxes)
	if ax_ecc > 0:
		ax.set_yticklabels([])
	if ax_ecc == 0 and ax_mu>0:
		ax.set_yticklabels(["%i" % i for i in xrange(1,5)])
	if ax_mu == 3 and ax_ecc<3:
		ax.set_xticklabels(["%i" % (90*i) for i in xrange(0,4)])


# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0,wspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-4]], visible=False)

color_label='$\log_{10}[e_{max}]$'
cax = f.add_axes([0.92,0.11,0.015,0.77])
cbar = plt.colorbar(cmmapable,cax=cax,orientation='vertical')#,ticks=np.linspace(vmin,vmax,11)) #draw colorbar  ,ticks =[0.0,0.2,0.4,0.6,0.8,1.0]
cbar.set_label(color_label,fontsize=fs)

cax.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)


f.text(0.5, 0.04, "Mean Anomaly (deg.)", ha='center',fontsize = fs)
f.text(0.04, 0.5, "$a_p/a_{bin}$", va='center', rotation='vertical',fontsize = fs)

axarr[0,3].text(1.0,1.02,"$\lambda_{bin}$ = %i$^\circ$" % (0), color='k',fontsize='small',horizontalalignment='right',weight='bold',transform=ax.transAxes)

plt.savefig("Q18_Fig%i.png" % fig_num, dpi = 300, bbox_inches = 'tight')
plt.close()
