from pylab import *
import numpy
import math
from subprocess import check_output
from re import split
from copy import  deepcopy

mcen=check_output("(cat param.in | head -n 29 | tail -1 | cut -d= -f2)",shell=True)

G=2.959122082855911E-4 #AU^3/M_sun/day^2

def xv2el(data,mu):
  x=data[1:4]
  v=data[4:]
  # Adapted from swifter routine
  TINY=4.E-15
  ellipse=-1
  parabola=0
  hyperbola=1

  a = 0.0
  ecc = 0.0
  inc = 0.0
  capom = 0.0
  omega = 0.0
  capm = 0.0
  
  r = sqrt(dot(x, x))
  v2 = dot(v, v)
  hx = x[1]*v[2] - x[2]*v[1]
  hy = x[2]*v[0] - x[0]*v[2]
  hz = x[0]*v[1] - x[1]*v[0]
  h2 = hx*hx + hy*hy + hz*hz
  h = sqrt(h2)
  if (h2 == 0.0): return
  rdotv = dot(x,v)
  energy = 0.5*v2 - mu/r
  fac = hz/h
  if (fac < -1.0): inc = pi
  elif (fac < 1.0): inc = math.acos(fac)
  fac = sqrt(hx*hx + hy*hy)/h
  if (fac < TINY): 
    u = math.atan2(x[1], x[0])
    if (hz < 0.0): u = -u
  else:
    capom = math.atan2(hx, -hy)
    if(sin(inc)==0.): u=0. # RAS-- to get rid of error
    else: u = math.atan2(x[2]/sin(inc), x[0]*cos(capom) + x[1]*sin(capom))
  if (capom < 0.0): capom = capom + 2.*pi
  if (u < 0.0): u = u + 2.*pi
  if (abs(energy*r/mu) < sqrt(TINY)): iorbit_type = parabola
  else: 
    a = -0.5*mu/energy
    if (a < 0.0): 
      fac = -h2/(mu*a)
      if (fac > TINY): iorbit_type = hyperbola
      else: iorbit_type = parabola
    else: iorbit_type = ellipse
  if (iorbit_type == ellipse):
    fac = 1.0 - h2/(mu*a)
    if (fac > TINY):
      ecc = sqrt(fac)
      cape = 0.0
      face = (a - r)/(a*ecc)
      if (face < -1.0): cape = pi
      elif (face < 1.0): cape = math.acos(face)
      if (rdotv < 0.0): cape = 2.*pi - cape
      fac = 1.0 - ecc*cos(cape)
      cw = (cos(cape) - ecc)/fac
      sw = sqrt(1.0 - ecc*ecc)*sin(cape)/fac
      w = math.atan2(sw, cw)
      if (w < 0.0): w = w + 2.*pi
    else:
      cape = u
      w = u
    capm = cape - ecc*sin(cape)
  elif (iorbit_type==parabola):
    a = 0.5*h2/mu
    ecc = 1.0
    w = 0.0
    fac = 2.0*a/r - 1.0
    if (fac < -1.0): w = pi
    elif (fac < 1.0): w = math.acos(fac)
    if (rdotv < 0.0): w = 2.*pi - w
    tmpf = tan(0.5*w)
    capm = tmpf*(1.0 + tmpf*tmpf/3.0)
  elif (iorbit_type==hyperbola):
    ecc = sqrt(1.0 + fac)
    tmpf = (a - r)/(a*ecc)
    if (tmpf < 1.0): tmpf = 1.0
    capf = log(tmpf + sqrt(tmpf*tmpf - 1.0))
    if (rdotv < 0.0): capf = -capf
    fac = ecc*cosh(capf) - 1.0
    cw = (ecc - cosh(capf))/fac
    sw = sqrt(ecc*ecc - 1.0)*sinh(capf)/fac
    w = math.atan2(sw, cw)
    if (w < 0.0): w = w + 2.*pi
    capm = ecc*sinh(capf) - capf
  omega = u - w
  if (omega < 0.0): omega = omega + 2.*pi
  return [a, ecc, degrees(inc), degrees(omega), degrees(capom), degrees(capm)]

def cen2jac(data):
  dist=sqrt(add(square(data[:,1]),add(square(data[:,2]),square(data[:,3]))))
  ind=argsort(ravel(dist))
  data2=data[ind]
  data3=deepcopy(data2)
  for i in range(len(data2)):
    z=0
    com=matrix([0.,0.,0.])
    comv=matrix([0.,0.,0.])
    summ=0.
    while z <= i:
      com+=multiply(data2[z,0],data2[z,1:4])
      comv+=multiply(data2[z,0],data2[z,4:])
      summ+=data2[z,0]
      z+=1
    com=com/summ
    comv=comv/summ
    newx=subtract(data2[i,1:4],com)
    newv=subtract(data2[i,4:],comv)
    data3[i,1:4]=newx
    data3[i,4:]=newv
  return data3[argsort(ind)]

def elout(dics,mykeys,low,high):

  dic=dict()
  z=0
  for i in mykeys:
    dic[i]=dics[i][low:high,:]
    if z > 0: 
      if low==0:
        with open(i+'.el','a') as f: f.write(i+' Jacobi Orbital Elements and Coordinates \n Time (years)\t mass\t a\t e\t i\t omega\t Omega\t M\t dist_jac\t dist_cen\t x \t y\t z\t u\t v\t w\n')
    z+=1
  ## Convert to jacobi
  data=zeros((len(mykeys),7))
  for i in range(len(dic[mykeys[0]])): 
    for k in range(len(mykeys)): data[k]=dic[mykeys[k]][i,1:]
    jac=cen2jac(data)
    bin_el= xv2el(data[1],(data[0,0]+data[1,0])*G) 
    dist=sqrt(add(square(jac[1,1]),add(square(jac[1,2]),square(jac[1,3]))))
    dist2=sqrt(add(square(data[1,1]),add(square(data[1,2]),square(data[1,3]))))
    temp= list(ravel(dic[mykeys[1]][i,:2])) + bin_el +[dist]+[dist2]+ list(ravel(jac[1,1:]))
    with open(mykeys[1]+'.el','a') as f: f.write(' '.join(map(str,temp))+'\n')
    for k in arange(2,len(mykeys)):
      mu=sum(ravel(jac[:k+1,0]))*G
      el=xv2el(jac[k],mu)
      dist=sqrt(add(square(jac[k,1]),add(square(jac[k,2]),square(jac[k,3]))))
      dist2=sqrt(add(square(data[k,1]),add(square(data[k,2]),square(data[k,3]))))
      temp= list(ravel(dic[mykeys[k]][i,:2])) + el +[dist]+[dist2]+ list(ravel(jac[k,1:]))
      with open(mykeys[k]+'.el','a') as f: f.write(' '.join(map(str,temp))+'\n')
  return 


names=loadtxt('allout',dtype=('S60'),unpack=True)
dic=dict()
mykeys=['STAR1']

convert = lambda text: int(text) if text.isdigit() else text
alphanum_key = lambda key: [convert(c) for c in split('([0-9]+)', key)]
temp =sorted(names, key = alphanum_key)

names= [temp[-1]]+temp[:len(temp)-1]

for i in range(len(names)):
    tmp=names[i].split('/')[-1]
    mykeys=mykeys+[tmp.split('.')[0]]

for i in range(len(names)): dic[mykeys[i+1]]=matrix(loadtxt(names[i],skiprows=4))
# star 1
dic[mykeys[0]] = matrix(zeros(shape(dic[mykeys[1]])))
dic[mykeys[0]][:,0]=dic[mykeys[1]][:,0]
dic[mykeys[0]][:,1]=mcen

shapes=zeros(len(mykeys))
for i in range(len(mykeys)): shapes[i]=shape(dic[mykeys[i]])[0]
shapes_ind=argsort(shapes)
unique = array(sort(list(set(list(shapes)))))

shortkeys=mykeys
if unique[0]>0: elout(dic,mykeys,0,unique[0])
for i in range(len(unique)-1):
  shorties=list(array(mykeys)[where(shapes==unique[i])[0]])
  p=len(shorties)
  while p>0:
   shortkeys=[x for x in shortkeys if x != shorties[p-1]]
   p-=1
  elout(dic,shortkeys,unique[i],unique[i+1])
