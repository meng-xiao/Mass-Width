# Mass-Width

## first check out the combine package under CMSSW_7_1_5, this might be out-synced with combine recommendation, but for a first try it should be OK
git clone https://github.com/meng-xiao/HiggsAnalysis-CombinedLimit HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git checkout slc6-root5.34.17

## compile the pacakge
cd $CMSSW_BASE/src 
scram b -j8

## put clean.c file somewhere in the release
mkdir workspace125_onshell 
## parameters for clean.c 
## chan: 4e, 2e2mu, 4mu , 
## vbf category =0, 1
## onshell=0, 1
root -q -b "clean.c(\"2e2mu\",0,1)"  

#onshell only workspace running
text2workspace.py hzz4l_13TeV_onshell.txt -o hzz4l_13TeV_onshell.root -P HiggsAnalysis.CombinedLimit.HighmassModel:HighmassModel --PO=mwAsPOI -v 3
combine -M MultiDimFit hzz4l_13TeV_onshell.root -n scan -t -1 --saveToys --setPhysicsModelParameters mean_pole=125,sigma_pole=0.004,r=1.,rvbf_ggh=1 --setPhysicsModelParameterRanges mean_pole=120,130:sigma_pole=0.0001,5.0001 --algo=grid --points=10
