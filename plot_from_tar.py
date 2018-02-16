import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import griddata
import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.cm as cm
import sys
from matplotlib import rcParams
import tarfile
from rebound.interruptible_pool import InterruptiblePool
import multiprocessing as mp
import subprocess as sb

rcParams.update({'font.size': 20})
cmap = cm.gnuplot_r
vmin = -2.
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

mu = float(sys.argv[1])
eb = float(sys.argv[2])
fname = 'MaxEcc_[%1.3f,%1.3f].txt' % (mu,eb)
tscale = 1e5

#Need to download the tar archive on Zenodo (http://doi.org/10.5281/zenodo.1174228)
if os.path.isfile("MaxEcc.tar.gz"):
	tar = tarfile.open("MaxEcc.tar.gz", "r:gz")
	filelist = tar.getnames()
	fidx = filelist.index(fname)


	tardata = tar.extractfile(filelist[fidx])
	data = np.genfromtxt(tardata,delimiter = ',',comments = '#')
	unstab = np.where(data[:,-1]<tscale)[0]
	stab = np.where(data[:,-1]==tscale)[0]

	X = data[:,1]
	Y = data[:,0]		
	Z = data[:,2]#-data[:,3]
	Z[unstab] = 1.5
	y_cr = []
	a_crit = -1.
	for i in xrange(0,401):
		yi = 1.01 + i*0.01
		row_x = np.where(np.abs(Y-yi)<1e-6)[0]
		all_stab = np.where(data[row_x,-1]==tscale)[0]
		if len(all_stab)>90:
			y_cr.append(yi)
	if len(y_cr)>0:
		a_crit = np.min(y_cr)

	xi = np.linspace(0,180,91)
	yi = np.linspace(1.,5.,401)
	zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'nearest')

	bins =200


	fig = plt.figure(1,figsize=(aspect*width,width),dpi=300)
	ax = fig.add_subplot(111)

	cs = ax.contourf(xi,yi,np.log10(zi),bins, cmap = cmap,vmin=vmin,vmax=vmax)

	if len(y_cr)>0:
		ax.axhline(a_crit,linestyle='-',color='k',lw=2)	

	ax.tick_params(axis='both', direction='out',length = 4.0, width = 4.0)

	ax.set_xticks(np.arange(0,180+30,30))
	yticks = np.arange(1.,6.0,1.)

	ax.set_yticks(yticks)
	ax.set_ylim(1,5.)
	ax.text(1.0,0.9,"(%1.3f,%1.2f)" % (mu,eb), color='black',fontsize='medium',horizontalalignment='right',weight='bold',transform=ax.transAxes)
	ax.text(0.0,0.02,"$a_{crit}$ = %1.2f" % a_crit, color='black',fontsize='small',horizontalalignment='left',weight='bold',transform=ax.transAxes)


	color_label='$\log_{10}[e_{max}]$'
	cax = fig.add_axes([0.92,0.11,0.015,0.77])
	cbar= plt.colorbar(cmmapable,cax=cax,orientation='vertical')
	cbar.set_label(color_label,fontsize=fs)


	cbar.ax.tick_params(labelsize=fs)

	ax.text(1.0,1.02,"$\lambda_{bin}$ = %i$^\circ$" % (0), color='black',fontsize='small',horizontalalignment='right',weight='bold',transform=ax.transAxes)
	ax.set_ylabel("$a_p/a_{bin}$",fontsize = fs)

	ax.set_xlabel("Mean Anomaly (deg.)", fontsize = fs)

	plt.savefig("MuEcc_[%1.3f,%1.2f]_emax.png" % (mu,eb), dpi = 300, bbox_inches = 'tight')
	plt.close()
else:
	print "Please download the tar archive on Zenodo (http://doi.org/10.5281/zenodo.1174228)"

