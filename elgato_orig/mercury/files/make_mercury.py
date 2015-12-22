from pylab import *
import numpy
from numpy.random import rayleigh,uniform
from merc import *

G=2.959122082855911E-4 #AU^3/M_sun/day^2

npl = 10

# Set up central star
m1=0.5; m2=0.5
a_star=0.01; e_star=0.; i_star=0.; om_star=0.; Om_star=0.; M_star=180.
sx,sy,sz,svx,svy,svz=el2xv((m1+m2)*G,a_star*(1.-e_star),e_star,radians(i_star),radians(Om_star),radians(Om_star-om_star),radians(M_star))
exes=[sx,sy,sz]
cen=matrix([0.,0.,0.])
com=(m1*cen+m2*matrix(exes))/(m1+m2)
star2=list(ravel(exes-com))
star1=list(ravel(cen-com))

# set up key values
mykeys=['star1','star2']
for z in range(npl):
    mykeys=mykeys+['pl'+str(z+1)]

# set up holding array
stats=zeros((len(mykeys),7))
stats[0]=[m1,0.,0.,0.,0.,0.,0.]; stats[1]=[m2,a_star,e_star,i_star,om_star,Om_star,M_star]
x=[star1[0],star2[0]]; y=[star1[1],star2[1]]; z=[star1[2],star2[2]]

# randomize planets
for k in range(npl):
    pl_ecc=rayleigh(0.1)
    while pl_ecc>=1.: pl_ecc=rayleigh(0.1)
    pl_i=rayleigh(5.73)
    pl_a=10.**(-1.+3.*uniform())
    pl_m=10.**(-4.+2.*uniform())  # Mj/MS = 9.5492E-4
    pl_om=uniform(0.,360.)
    pl_capom=uniform(0.,360.)
    pl_capm=uniform(0.,360.)
    plx,ply,plz,plvx,plvy,plvz=el2xv((sum(ravel(stats[:,0]))+pl_m)*G,pl_a*(1.-pl_ecc),pl_ecc,radians(pl_i),radians(pl_capom),radians(pl_capom-pl_om),radians(pl_capm))
    stats[k+2]=[pl_m,pl_a,pl_ecc,pl_i,pl_om,pl_capom,pl_capm]
    x +=[plx]; y +=[ply]; z +=[plz]

# sort by distance
dist=sqrt(add(square(x),add(square(y),square(z))))
ind=argsort(dist)
stats=array(stats)[ind]

# shift to center of mass coordinates
data=zeros((len(stats),7))
data[:,0]=stats[:,0]
com=matrix([0.,0.,0.]); comv=matrix([0.,0.,0.])
comt=matrix([0.,0.,0.]); comvt=matrix([0.,0.,0.])
summ=0.
x=0.;y=0.;z=0.;u=0.;v=0.;w=0.
for k in range(len(data)-1):
  zz=k+1
  x,y,z,u,v,w=el2xv(sum(ravel(stats[:zz+1,0]))*G,stats[zz,1]*(1.-stats[zz,2]),stats[zz,2],radians(stats[zz,3]),radians(stats[zz,5]),radians(stats[zz,5]-stats[zz,4]),radians(stats[zz,6])) 
  data[zz,1:]=array([x,y,z,u,v,w])
  newx=add(data[zz,1:4],comt)
  newv=add(data[zz,4:],comvt) 
  data[zz,1:4]=newx
  data[zz,4:]=newv
  com+=multiply(data[zz,0],data[zz,1:4])
  comv+=multiply(data[zz,0],data[zz,4:])
  summ=sum(data[:zz+1,0])
  comt=com/summ 
  comvt=comv/summ

# write big.in
f=open('big.in','w')
f.write(')O+_06 Integration parameters  (WARNING: Do not delete this line!!)\n')
f.write(') Lines beginning with ) are ignored.\n')
f.write(')---------------------------------------------------------------------\n')
f.write(' style (Cartesian, Asteroidal, Cometary) = Cartesian\n')
f.write(' epoch (in days) = 0.\n')
f.write(')---------------------------------------------------------------------\n')
for i in arange(1,len(mykeys)):
    f.write(' '+str(mykeys[i].upper())+' m='+str(data[i,0])+'\n')
    f.write(' '+str(data[i,1])+' '+str(data[i,2])+' '+str(data[i,3])+'\n')
    f.write(' '+str(data[i,4])+' '+str(data[i,5])+' '+str(data[i,6])+'\n')
    f.write(' '+str(0.)+' '+str(0.)+' '+str(0.)+'\n')
f.close()

# write stats file
f=open('stats','w')
f.write('Initial conditions for this run\n')
f.write('\nname\tmass\ta\te\ti\tomega\tOmega\tM\tx\ty\tz\tvx\tvy\tvz  \n')
for i in range(len(mykeys)):
    f.write(mykeys[i].upper()+' '+' '.join(map(str,stats[i]))+' ')
    f.write(str(data[i,1])+' '+str(data[i,2])+' '+str(data[i,3])+' ')
    f.write(str(data[i,4])+' '+str(data[i,5])+' '+str(data[i,6])+'\n')
f.close()

