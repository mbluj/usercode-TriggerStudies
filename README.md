
Tools for tau triggers validation for 7_4_X

cmsrel CMSSW_7_4_6; #or similar 74X

cd CMSSW_7_4_6/src;

cmsenv;

git clone https://github.com/mbluj/usercode-TriggerStudies.git TriggerStudies/Tau;

(cd TriggerStudies/Tau; git checkout 7_4_X)

scram b


