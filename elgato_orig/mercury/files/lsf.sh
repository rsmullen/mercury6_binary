#!/bin/bash
thisdir=JT08_nocollisions_npl10_300kyr
#BSUB -n 1
#BSUB -e err_JT08_%I
#BSUB -o out_JT08_%I
#BSUB -q "windfall"
#BSUB -u rsmullen
#BSUB -J JT08[1-100]
#BSUB -R "span[ptile=1]"

module load python/2.7.3

mkdir /rsgrps/kkratterstudents/rsmullen/$thisdir
cp -rf /rsgrps/kkratterstudents/rsmullen/files/param.in /rsgrps/kkratterstudents/rsmullen/$thisdir
cp -rf /rsgrps/kkratterstudents/rsmullen/files/files.in /rsgrps/kkratterstudents/rsmullen/$thisdir
cp -rf /rsgrps/kkratterstudents/rsmullen/files/element.in /rsgrps/kkratterstudents/rsmullen/$thisdir
cp -rf /rsgrps/kkratterstudents/rsmullen/files/close.in /rsgrps/kkratterstudents/rsmullen/$thisdir
cp -rf /rsgrps/kkratterstudents/rsmullen/files/small.in /rsgrps/kkratterstudents/rsmullen/$thisdir

mkdir /rsgrps/kkratterstudents/rsmullen/$thisdir/run${LSB_JOBINDEX}

cd /rsgrps/kkratterstudents/rsmullen/$thisdir/run${LSB_JOBINDEX}/

ln -s /home/u14/rsmullen/mercury/mercury6 
ln -s /home/u14/rsmullen/mercury/element6 
ln -s /home/u14/rsmullen/mercury/close6 
ln -s /home/u14/rsmullen/mercury/message.in 

ln -s /rsgrps/kkratterstudents/rsmullen/files/param.in
ln -s /rsgrps/kkratterstudents/rsmullen/files/files.in
ln -s /rsgrps/kkratterstudents/rsmullen/files/element.in
ln -s /rsgrps/kkratterstudents/rsmullen/files/close.in
ln -s /rsgrps/kkratterstudents/rsmullen/files/small.in

cp -rf /rsgrps/kkratterstudents/rsmullen/files/make_mercury.py .

python make_mercury.py
./mercury6
./element6
./close6
ls -1 *.aei > allout
ls -1 *.clo > allclose
