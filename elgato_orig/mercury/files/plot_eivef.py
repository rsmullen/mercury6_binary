from pylab import *
import matplotlib
import numpy
import math
from subprocess import call,check_call
import matplotlib.cm as cm

n=10

dirs=loadtxt('allruns',dtype='S64',unpack=True)

ei=[]; ef=[]; ii=[]; iif=[]; ai=[]; af=[]; mi=[]; mf=[]; tf=[]
eie=[]; efe=[]; iie=[]; iife=[]; aie=[]; afe=[]; mie=[]
for z in range(len(dirs)):
    files=loadtxt(dirs[z]+'/allelem',dtype='S64',unpack=True)
    print dirs[z]
    for q in range(len(files)-1):
        with open(dirs[z]+'/'+files[q]) as f:
            f.readline()
            f.readline()
            first=[f.readline()]
            for i in range(n-1): first+=[f.readline()]
            k=0
            last=[]
            for line in f: 
                if k >= n: 
                    last=last[1:]+[line]
                else:
                    last+=[line]
                    k+=1
#            print dirs[z]+'/'+files[q]
        
        if ((float(last[-1].split(' ')[0])) > 490000.):
          mi+=[float( first[0].split(' ')[1])]
          these=[float(w.split(' ')[2]) for w in first]; ai+=[mean(these)]; aie+=[std(these)]
          these=[float(w.split(' ')[2]) for w in last];  af+=[mean(these)]; afe+=[std(these)]
          these=[float(w.split(' ')[3]) for w in first]; ei+=[mean(these)]; eie+=[std(these)]
          these=[float(w.split(' ')[3]) for w in last];  ef+=[mean(these)]; efe+=[std(these)]
          these=[float(w.split(' ')[4]) for w in first]; ii+=[mean(these)]; iie+=[std(these)]
          these=[float(w.split(' ')[4]) for w in last]; iif+=[mean(these)];iife+=[std(these)]
#    if z > 3: break
#          tf+=[float(last[-1].split(' ')[0])]

#these = where(ravel(tf) > 490000.)[0]
#ei=array(ei)[these]; ef=array(ef)[these]; ii=array(ii)[these]; iif=array(iif)[these]; ai=array(ai)[these]; af=array(af)[these]; mi=array(mi)[these] #solar
#eie=array(eie)[these]; efe=array(efe)[these]; iie=array(iie)[these]; iife=array(iife)[these]; aie=array(aie)[these]; afe=array(afe)[these]

mi=divide(mi,0.00095426) #jupiter
#mi=mi/3.002458E-6 #earth

fig = figure(figsize=(10,10))
ax=subplot(311)
scatter(ei,ef,marker='o',c=log10(mi),cmap=cm.jet,edgecolors='None',zorder=10)
errorbar(ei,ef,xerr=eie,yerr=efe,color='k',marker=None,capsize=0,capthick=0,ls='', mew=0,zorder=0)
plot(arange(0.,1.1,.1),arange(0.,1.1,.1),'k:',zorder=0)
xlim([0.,1.])
ylim([0.,1.])
xlabel('Initial eccentricity',fontsize=18)
ylabel('Final \neccentricity',fontsize=18)
ax.tick_params(which='major', labelsize=14)

ax=subplot(312)
scatter(ii,iif,marker='o',c=log10(mi),cmap=cm.jet,edgecolors='None',zorder=10)
errorbar(ii,iif,xerr=iie,yerr=iife,color='k',marker=None,capsize=0,capthick=0,ls='', mew=0,zorder=0)
plot(arange(0.,180.1,.1),arange(0.,180.1,.1),'k:',zorder=0)
xlim([0.,90.])
ylim([0.,180.])
xlabel('Initial inclination',fontsize=18)
ylabel('Final \ninclination',fontsize=18)
ax.tick_params(which='major', labelsize=14)

ax=subplot(313)
loglog(arange(0.,1000.1,.1),arange(0.,1000.1,.1),'k:',zorder=0)
zz=scatter(ai,af,marker='o',c=log10(mi),cmap=cm.jet,edgecolors='None',zorder=10)
errorbar(ai,af,xerr=aie,yerr=afe,color='k',marker=None,capsize=0,capthick=0,ls='', mew=0,zorder=0)
xlim([0.,1000.])
ylim([0.,1000.])
xlabel('Initial semi-major axis',fontsize=18)
ylabel('Final \nsemi-major axis',fontsize=18)
ax.tick_params(which='major', labelsize=14)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar=fig.colorbar(zz, cax=cbar_ax)
cbar.set_label(r'log$_{10}$Mass (Jupiter masses)',fontsize=18)

subplots_adjust(bottom=0.1, top=.97, wspace=None, hspace=.3)

savefig('params.eps')

fig = figure(figsize=(10,10))
ax=subplot(311)
weights = ones_like(ei)/len(ei)
hist([ei,ef], bins=arange(0.,1.1,.1),normed=0,color=['r','b'],histtype='step',alpha=1.,weights=[weights,weights])
xlabel('Eccentricity',fontsize=18)
ylabel('Fraction',fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
ax=subplot(312)
hist([ii,iif], bins=arange(0.,180.1,5.),normed=0,color=['r','b'],histtype='step',alpha=1.,weights=[weights,weights])
xlabel('Inclination (degrees)',fontsize=18)
ylabel('Fraction',fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
ax=subplot(313)
hist([ai,af], bins=logspace(-1.5,3.,20.),normed=0,color=['r','b'],histtype='step',alpha=1.,label=['Initial','Final'],weights=[weights,weights])
gca().set_xscale("log")
xlabel('Semi-major Axis (AU)',fontsize=18)
ylabel('Fraction',fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
handles, labels = ax.get_legend_handles_labels()
figlegend(handles,('Initial','Final'),loc='center right')
subplots_adjust(bottom=0.07, top=.97, wspace=None, hspace=0.25)
savefig('hist.eps')
#show()



