from pylab import *
import math,numpy,scipy,matplotlib


t,m,a,e,i,x,y,z,u,v,w=loadtxt('STAR2.aei',skiprows=4,unpack=True)
t2,m2,a2,e2,i2,x2,y2,z2,u2,v2,w2=loadtxt('elgato_orig/mercury/STAR2.aei',skiprows=4,unpack=True)

print max(x-x2),max(y-y2),max(z-z2),max(u-u2),max(v-v2),max(w-w2)

plot( x,y)

show()
