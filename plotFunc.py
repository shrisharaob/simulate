from numpy import *
from pylab import *
x = loadtxt("outputFile.csv")
t = x[:, 0]
v = x[:, 1]
n = x[:, 2]
z = x[:, 3]
h = x[:, 4]
cur = x[:, 5]
#plot(t, v, 'b', t, cur, 'g', t, z, 'k', t, h,'r')
##plot(t, cur, 'g', t, 1-h, 'k', t, n, 'r')
##plot(t[0:1e3], v[0:1e3] ,'k')  
##plot(t, v)
plot(t, v, t, cur, t, z)
show()

