import numpy as np
import sys
from scipy.interpolate import griddata
from scipy.interpolate import interp2d

mu = float(sys.argv[1])
e_bin = float(sys.argv[2])

data = np.genfromtxt("a_crit.txt",delimiter=',',comments='#')

X = data[:,0]
Y = data[:,1]
Z = data[:,2]

xi = np.linspace(0,0.5,51)
yi = np.linspace(0,0.8,81)
zi = griddata((X,Y),Z,(xi[None,:],yi[:,None]),method = 'linear',fill_value=0)

f = interp2d(xi, yi, zi, kind='linear')

print "a_c = ",f(mu,e_bin)[0]
