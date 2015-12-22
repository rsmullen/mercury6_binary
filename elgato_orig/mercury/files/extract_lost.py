from pylab import *
import matplotlib
import numpy
import math
from subprocess import call,check_call

maxt=500000.

#call('ls -1d '+'run* > allruns',shell=True)
names=loadtxt('allruns',dtype='S64',unpack=True)

T=[]
lost=[]
z=0
zz=0
j=0
timearr=arange(0.,maxt,100.)
npl=zeros((len(names),len(timearr)))
npl2=zeros((len(names),len(timearr)))
for i in names:
    t=[]
    #print i
    for line in open(i+'/info.out'):
        if "PL" in line: t+=[float(line.split('ejected at')[-1].split('years')[0].strip())]
    T+=t
    lost+=[len(t)]
        #if len(t)>7: print i,[len(t)]
    temp=0
    for k in range(len(timearr)):
        if (array(t).any() < timearr[k]): 
            temp=len(where(array(t)<timearr[k])[0])
        npl[j,k]=10-temp
        npl2[j,k]=temp
   
        #print i, array(t)[array(where(array(t)<1000.)[0])]
        #times=array(t)[array(where(array(t)<1000.)[0])]
        #if len(times)>1: zz+=1
        #z+=1
    j+=1

#print npl2

n=[]
cdf=[]
for j in range(len(timearr)):
    #print mean(npl[:,j])
    n+=[mean(npl[:,j])]
    #print sum((npl2[:,:j+1]))
    temp=sum(npl2[:,j])
    cdf+=[temp]

cdf=cdf/(sum(ravel(npl2[:,-1])))

fsize=18
lsize=14

fig=figure(figsize=(10.,13.))

ax=subplot(311)
counts,bins,patches=hist(lost,bins=range(11),normed=True)
xlim([0,10])
bin_centers = 0.5 * np.diff(bins) + bins[:-1]
for bins, x in zip(bins, bin_centers):
    # Label the raw counts
    ax.annotate(str(bins), xy=(x, 0), xycoords=('data', 'axes fraction'),xytext=(0, -8), textcoords='offset points', va='top', ha='center',fontsize=lsize)
xlabel('Total number of planets lost after 500kyr', labelpad=20,fontsize=fsize)
ylabel('Fraction of sample',fontsize=fsize)
ax.set_xticklabels([])
ax.tick_params(axis='both', which='major', labelsize=lsize)

ax=subplot(312)
hist(T,bins=logspace(1.,log10(maxt),15))
gca().set_xscale("log")
xlabel('Time of ejection',fontsize=fsize)
ylabel(r'n$_{pl}$ lost',fontsize=fsize)
ax.tick_params(axis='both', which='major', labelsize=lsize)
#print z,zz

p1=fig.add_subplot(313)
semilogx(timearr,n)
xlabel('Time',fontsize=fsize)
ylabel(r'Average n$_{pl}$ remaining',fontsize=fsize)
ylim([0,10])
p1.tick_params(axis='both', which='major', labelsize=lsize)
p1.spines['right'].set_color('red')
p1.spines['left'].set_color('blue')

p2=fig.add_subplot(313,sharex=p1, frameon=False)
semilogx(timearr,cdf,'r')
p2.yaxis.tick_right()
p2.yaxis.set_label_position("right")
xlabel('Time',fontsize=fsize)
ylabel(r'CDF (n$_{pl}$ lost)',fontsize=fsize)
p2.tick_params(axis='both', which='major', labelsize=lsize)
p2.spines['right'].set_color('red')

subplots_adjust(bottom=0.07, top=.97, wspace=None, hspace=None)
savefig('extract.eps')
show()


