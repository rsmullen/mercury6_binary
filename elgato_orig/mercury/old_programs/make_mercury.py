from pylab import *
import numpy
from numpy import random

G=2.959122082855911E-4 #AU^3/M_sun/day^2

npl = 10

def jacobi(x,y,z,vx,vy,vz,dic,keys):
    com=[0.,0.,0.]
    comv=[0.,0.,0.]
    summ=0.
    for i in range(len(keys)):
        com+=multiply(dic[keys[i]][0],dic[keys[i]][1:4])
        comv+=multiply(dic[keys[i]][0],dic[keys[i]][4:])
        summ+=dic[keys[i]][0]
    com=com/summ
    comv=comv/summ
    newx=subtract([x,y,z],com)
    newv=subtract([vx,vy,vz],comv)
    return newx[0],newx[1],newx[2],newv[0],newv[1],newv[2]

def EA(M,ecc):
    if ecc<0.8: thisE=M
    else: thisE=pi
    niter = 10
    delta = 1.E-18
    for z in range(niter):
        F=thisE-ecc*sin(thisE)-M
        Fp=1.-ecc*cos(thisE)
        newE=thisE-F/Fp
        if abs((newE-thisE)) < delta:
            thisE=newE
            break
        thisE=newE
    return thisE

def el2xv(a,ecc,i,omega,capom,capm,M,m):
    i=radians(i)
    omega=radians(omega)
    capom=radians(capom)
    capm=radians(capm)
    mu=G*(M+m)
    E=EA(capm,ecc)

    l1=cos(capom)*cos(omega)-sin(capom)*sin(omega)*cos(i)
    m1=sin(capom)*cos(omega)+cos(capom)*sin(omega)*cos(i)
    n1=sin(omega)*sin(i)

    l2=-cos(capom)*sin(omega)-sin(capom)*cos(omega)*cos(i)
    m2=-sin(capom)*sin(omega)+cos(capom)*cos(omega)*cos(i)
    n2=cos(omega)*sin(i)

    n=sqrt(mu)*a**(-3./2)
    b=a*(1.-ecc**2.)**(1./2)
    r=a*(1.-ecc*cos(E))

    x=a*l1*cos(E)+b*l2*sin(E)-a*ecc*l1
    y=a*m1*cos(E)+b*m2*sin(E)-a*ecc*m1
    z=a*n1*cos(E)+b*n2*sin(E)-a*ecc*n1

    vx=n*a/r*(b*l2*cos(E)-a*l1*sin(E))
    vy=n*a/r*(b*m2*cos(E)-a*m1*sin(E))
    vz=n*a/r*(b*n2*cos(E)-a*n1*sin(E))
    return x,y,z,vx,vy,vz


dic=dict()
#      m  x  y  z  vx vy vz
star1=[.5]
star2=[.5]
a_star=0.01
sx,sy,sz,svx,svy,svz=el2xv(a_star,0.,0.,0.,0.,0.,star1[0],star2[0])
star2+=[-sx/2.]+[sy]+[sz]+[svx]+[-svy/2.]+[svz]
star1+=[sx/2.]+[sy]+[sz]+[svx]+[svy/2.]+[svz]
dic['star1']=star1
dic['star2']=star2
mykeys=['star1','star2']
for z in range(npl):
    mykeys=mykeys+['pl'+str(z+1)]
#print(mykeys)
aa=[0.,a_star]
mm=[star1[0],star2[0]]
ee=[0.,0.]
ii=[0.,0.]
om=[0.,0.]
Om=[0.,0.]
M=[0.,180.]

a_array=sorted([10.**(-1.+3.*random.uniform()) for i in range(npl)])

for k in range(npl):
    pl_ecc=random.rayleigh(0.1)
    #pl_ecc=random.rayleigh(0.5)  #added
    while pl_ecc>=1.: pl_ecc=random.rayleigh(0.1)
    pl_i=random.rayleigh(5.73)
    pl_a=a_array[k]
    #pl_a=10.**(-1.+2.*random.uniform())   #added
    pl_m=10.**(-4.+2.*random.uniform())  # Mj/MS = 9.5492E-4
    #pl_m=1.E-3  #added
    pl_om=random.uniform(0.,360.)
    pl_capom=random.uniform(0.,180.)
    pl_capm=random.uniform(0.,360.)

    plx,ply,plz,plvx,plvy,plvz=el2xv(pl_a,pl_ecc,pl_i,pl_om,pl_capom,pl_capm,sum(mm),pl_m)
    px,py,pz,pvx,pvy,pvz=jacobi(plx,ply,plz,plvx,plvy,plvz,dic,mykeys[:(k+2)])
    #print px-plx,py-ply,pz-plz,pvx-plvx,pvy-plvy,pvz-plvz
    dic[mykeys[k+2]]=[pl_m]+[px]+[py]+[pz]+[pvx]+[pvy]+[pvz]

    aa+=[pl_a]
    mm+=[pl_m]
    ee+=[pl_ecc]
    ii+=[pl_i]
    om+=[pl_om]
    Om+=[pl_capom]
    M +=[pl_capm]

move=dict()
for k in range(len(mykeys)):
    move[mykeys[k]]=numpy.concatenate([[dic[mykeys[k]][0]],subtract(dic[mykeys[k]][1:],dic[mykeys[0]][1:])])

f=open('big.in','w')
f.write(')O+_06 Integration parameters  (WARNING: Do not delete this line!!)\n')
f.write(') Lines beginning with ) are ignored.\n')
f.write(')---------------------------------------------------------------------\n')
f.write(' style (Cartesian, Asteroidal, Cometary) = Cartesian\n')
f.write(' epoch (in days) = 0.\n')
f.write(')---------------------------------------------------------------------\n')
index=1
for i in argsort(aa):
  if i!=0:
    f.write(' '+str(mykeys[i].upper())+' m='+str(move[mykeys[i]][0])+'\n')
    f.write(' '+str(move[mykeys[i]][1])+' '+str(move[mykeys[i]][2])+' '+str(move[mykeys[i]][3])+'\n')
    f.write(' '+str(move[mykeys[i]][4])+' '+str(move[mykeys[i]][5])+' '+str(move[mykeys[i]][6])+'\n')
    f.write(' '+str(0.)+' '+str(0.)+' '+str(0.)+'\n')
    #print(' '+str(mykeys[i].upper())+' m='+str(move[mykeys[i]][0])+'\n')
    #print(' '+str(move[mykeys[i]][1])+' '+str(move[mykeys[i]][2])+' '+str(move[mykeys[i]][3])+'\n')
    #print(' '+str(move[mykeys[i]][4])+' '+str(move[mykeys[i]][5])+' '+str(move[mykeys[i]][6])+'\n')
    #print(' '+str(0.)+' '+str(0.)+' '+str(0.)+'\n')
    index+=1
f.close()

f=open('stats','w')
f.write('Initial conditions for this run\n')
f.write('\nname\tmass\ta\te\ti\tomega\tOmega\tM\tx\ty\tz\tvx\tvy\tvz  \n')
for i in argsort(aa):
    f.write(mykeys[i].upper()+' '+str(mm[i])+' '+str(aa[i])+' '+str(ee[i])+' '+str(ii[i])+' '+str(om[i])+' '+str(Om[i])+' '+str(M[i])+' ')
    f.write(str(move[mykeys[i]][1])+' '+str(move[mykeys[i]][2])+' '+str(move[mykeys[i]][3])+' ')
    f.write(str(move[mykeys[i]][4])+' '+str(move[mykeys[i]][5])+' '+str(move[mykeys[i]][6])+'\n')
f.close()

