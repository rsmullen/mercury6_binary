from pylab import *
import numpy

# coordinate to elements adapted from mercury
def xv2el(gm,x,y,z,u,v,w):
        q=0.;p=0.;ecc=0.;i=0.;n=0.;l=0.
        PI=pi
        TWOPI=2.*pi
  
        hx = y * w  -  z * v
        hy = z * u  -  x * w
        hz = x * v  -  y * u
        h2 = hx*hx + hy*hy + hz*hz
        v2 = u * u  +  v * v  +  w * w
        rv = x * u  +  y * v  +  z * w
        r = sqrt(x*x + y*y + z*z)
        h = sqrt(h2)
        s = h2 / gm
  
  # Inclination and node
        ci = hz / h
        if (abs(ci)<1):
          i = math.acos(ci)
          n = math.atan2(hx,-hy)
          if (n<0): n = n + TWOPI
        else:
          if (ci>0): i = 0.0
          if (ci<0): i = PI
          n = 0.0
  
  # Eccentricity and perihelion distance
        temp = 1.0  +  s * (v2 / gm  -  2.0 / r)
        if (temp<=0.): ecc = 0.0
        else: ecc = sqrt (temp)
        q = s / (1.0 + ecc)
  
  # True longitude
        if (hy!=0):
          to = -hx/hy
          temp = (1.0 - ci) * to
          tmp2 = to * to
          true = math.atan2((y*(1.0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp))
        else: 
          true = math.atan2(y * ci, x)
        if (ci<0): true = true + PI
  
        if (ecc<3.E-8):
          p = 0.0
          l = true
        else:
          ce = (v2*r - gm) / (ecc*gm)
  
  # Mean anomaly for ellipse
          if (ecc<1.):
            if (abs(ce)>1): ce = fsign(1.0,ce)
            bige = math.acos(ce)
            if (rv<0): bige = TWOPI - bige
            l = bige - ecc*sin(bige)
          else:
  # Mean anomaly for hyperbola
            if (ce<1): ce = 1.0
            bige = math.log( ce + sqrt(ce*ce-1.0) )
            if (rv<0): bige = - bige
            l = ecc*sinh(bige) - bige
  
  # Longitude of perihelion
          cf = (s - r) / (ecc*r)
          if (abs(cf)>1): cf = fsign(1.0,cf)
          f = math.acos(cf)
          if (rv<0): f = TWOPI - f
          p = true - f
          p = (p + TWOPI + TWOPI) % TWOPI
  
        if (l<0): l = l + TWOPI
        if (l>TWOPI): l = l % TWOPI
  
        return q,ecc,i,p,n,l

def orbel_zget(q):
	iflag = 0
	if(q<0.0):
	  iflag = 1
          q = -q

	if (q<1.E-3): orbel_zget = q*(1.0 - (q*q/3.0)*(1.0 -q*q))
        else:
	   x = 0.50*(3.0*q + sqrt(9.0*(q**2) +4.0))
	   tmp = x**(1.0/3.0)
	   orbel_zget = tmp - 1.0/tmp

	if(iflag == 1):
           orbel_zget = -orbel_zget
	   q = -q
	
	return orbel_zget


#                    ORBEL_FLON.F

def orbel_flon(ec,capn):
        TINY=4.E-15
	IMAX = 10
	a11 = 156.0;a9 = 17160.0;a7 = 1235520.0
	a5 = 51891840.0;a3 = 1037836800.0
	b11 = 11.0*a11;b9 = 9.0*a9;b7 = 7.0*a7
	b5 = 5.0*a5; b3 = 3.0*a3

# Function to solve "Kepler's eqn" for F (here called
# x) for given e and CAPN. Only good for smallish CAPN 

	iflag = 0
	if( capn < 0.0) :
	   iflag = 1
	   capn = -capn

	a1 = 6227020800.0 * (1.0 - 1.0/ec)
	a0 = -6227020800.0*capn/ec
	b1 = a1

#  Set iflag nonzero if capn < 0., in which case solve for -capn
#  and change the sign of the final answer for F.
#  Begin with a reasonable guess based on solving the cubic for small F	


	a = 6.0*(ec-1.0)/ec
	b = -6.0*capn/ec
	sq = sqrt(0.25*b*b +a*a*a/27.0)
	biga = (-0.5*b + sq)**0.33333333333333330
	bigb = -(+0.5*b + sq)**0.33333333333333330
	x = biga + bigb
	orbel_flon = x
# If capn is tiny (or zero) no need to go further than cubic even for
# e =1.
        if( capn >= TINY):

	  for i in arange(0,IMAX,1):
	    x2 = x*x
	    f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
	    fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.0*x2)))))   
	    dx = -f/fp
	    orbel_flon = x + dx
#   If we have converged here there's no point in going on
            if(abs(dx) > TINY):
              x = orbel_flon
          
# Abnormal return here - we've gone thru the loop 
# IMAX times without convergence
          if(abs(dx) > TINY):
	    if(iflag == 1) :
	       orbel_flon = -orbel_flon
	       capn = -capn
	    print 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	    diff = ec*sinh(orbel_flon) - orbel_flon - capn
	    print 'N, F, ecc*sinh(F) - F - N : '
	    print capn,orbel_flon,diff
	    return

#  Normal return here, but check if capn was originally negative
 	if(iflag == 1) :
	   orbel_flon = -orbel_flon
	   capn = -capn

	return orbel_flon


#                    ORBEL_FGET.F
#Function to solve "Kepler's eqn" for F (here called
#x) for given e and CAPN. 
def orbel_fget(ec,capn):
        IMAX=10
# begin with a guess proposed by Danby	
        if( capn < 0.0) :
	   tmp = -2.0*capn/ec + 1.80
	   x = -log(tmp)
	else:
	   tmp = 2.0*capn/ec + 1.80
	   x = log( tmp)

	orbel_fget = x

	for i in range(IMAX):
          shx=sinh(x)
          chx=cosh(x)
	  esh = ec*shx
	  ech = ec*chx
	  f = esh - x - capn
	  fp = ech - 1.0  
	  fpp = esh 
	  fppp = ech 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.0)
	  dx = -f/(fp + dx*fpp/2.0 + dx*dx*fppp/6.0)
	  orbel_fget = x + dx
#  If we have converged here there's no point in going on
	  if(abs(dx) < TINY): return orbel_fget
	  x = orbel_fget

        print 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
        return 0


#                    ORBEL_FHYBRID.F

def orbel_fhybrid(ea,n):
	abn = n
	if(n<0.0): abn = -abn

	if(abn < 0.6360*ea -0.60) :
	  orbel_fhybrid = orbel_flon(ea,n)
        else :
	  orbel_fhybrid = orbel_fget(ea,n)

	return orbel_fhybrid

def mco_kep(ek,oldl):

      twopi = 2.0 * pi
      piby2 = .50 * pi

# Reduce mean anomaly to lie in the range 0 < l < pi
      if (oldl>=0) : l = oldl % twopi
      else: l = oldl % twopi + twopi
      sign = 1.0
      if (l>pi) :
        l = twopi - l
        sign = -1.0

      ome = 1.0 - ek

      if ((l>=.450)or(ek<.550)):

# Regions A,B or C in Nijenhuis
# -----------------------------

# Rough starting value for eccentric anomaly
        if (l<ome): 
          u1 = ome
        else:
          if (l>(pi-1.0-ek)): u1 = (l+ek*pi)/(1.0+ek)
          else:  u1 = l + ek

# Improved value using Halley's method
        flag = u1>piby2
        if (flag): x = pi - u1
        else:  x = u1
        x2 = x*x
        sn = x*(1.0 + x2*(-.16605 + x2*.00761) )
        dsn = 1.0 + x2*(-.49815 + x2*.03805)
        if (flag): dsn = -dsn
        f2 = ek*sn
        f0 = u1 - f2 - l
        f1 = 1.0 - ek*dsn
        u2 = u1 - f0/(f1 - .50*f0*f2/f1)
      else:

# Region D in Nijenhuis
# ---------------------

# Rough starting value for eccentric anomaly
        z1 = 4.0*ek + .50
        p = ome / z1
        q = .50 * l / z1
        p2 = p*p
        z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.0*q / ( z2 + p + p2/z2 )

# Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - .0750*u1*z3 / (ome + z1*z2 + .3750*z3)
        u2 = l + ek*u2*( 3.0 - 4.0*u2*u2 )

# Accurate value using 3rd-order version of Newton's method
# N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!

# First get accurate values for u2 - sin(u2) and 1 - cos(u2)
      bigg = (u2>piby2)
      if (bigg): z3 = pi - u2
      else: z3 = u2

      big = (z3>(.50*piby2))
      if (big): x = piby2 - z3
      else:  x = z3

      x2 = x*x
      ss = 1.0
      cc = 1.0

      ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
      cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. - x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. - x2/306.))))))))

      if (big) :
        z1 = cc + z3 - 1.0
        z2 = ss + z3 + 1.0 - piby2
      else:
        z1 = ss
        z2 = cc

      if (bigg):
        z1 = 2.0*u2 + z1 - pi
        z2 = 2.0 - z2

      f0 = l - u2*ome - ek*z1
      f1 = ome + ek*z2
      f2 = .50*ek*(u2-z1)
      f3 = ek/6.0*(1.0-z2)
      z1 = f0/f1
      z2 = f0/(f2*z1+f1)
      return sign*( u2 + f0/((f3*z1+f2)*z2+f1) )


#      MCO_EL2X.FOR    (ErikSoft  7 July 1999)
def el2xv (gm,q,ecc,i,p,n,l):

# Change from longitude of perihelion to argument of perihelion
      g = p - n

# Rotation factors
      si=sin(i)
      ci=cos(i)
      sg=sin(g)
      cg=cos(g)
      sn=sin(n)
      cn=cos(n)
      z1 = cg * cn
      z2 = cg * sn
      z3 = sg * cn
      z4 = sg * sn
      d11 =  z1 - z4*ci
      d12 =  z2 + z3*ci
      d13 = sg * si
      d21 = -z3 - z2*ci
      d22 = -z4 + z1*ci
      d23 = cg * si

# Semi-major axis
      a = q / (1.0 - ecc)
#
# Ellipse
      if (ecc<1.0):
        romes = sqrt(1.0 - ecc*ecc)
        temp = mco_kep (ecc,l)
        se=sin(temp)
        ce = cos(temp)
        z1 = a * (ce - ecc)
        z2 = a * romes * se
        temp = sqrt(gm/a) / (1.0 - ecc*ce)
        z3 = -se * temp
        z4 = romes * ce * temp
      else:
# Parabola
        if (ecc==1.0):
          ce = orbel_zget(l)
          z1 = q * (1.0 - ce*ce)
          z2 = 2.0 * q * ce
          z4 = sqrt(2.0*gm/q) / (1.0 + ce*ce)
          z3 = -ce * z4
        else:
# Hyperbola
          romes = sqrt(ecc*ecc - 1.0)
          temp = orbel_fhybrid(ecc,l)
          se=sinh(temp)
          ce=sinh(temp)
          z1 = a * (ce - ecc)
          z2 = -a * romes * se
          temp = sqrt(gm/abs(a)) / (ecc*ce - 1.0)
          z3 = -se * temp
          z4 = romes * ce * temp

      x = d11 * z1  +  d21 * z2
      y = d12 * z1  +  d22 * z2
      z = d13 * z1  +  d23 * z2
      u = d11 * z3  +  d21 * z4
      v = d12 * z3  +  d22 * z4
      w = d13 * z3  +  d23 * z4
      return  x,y,z,u,v,w


