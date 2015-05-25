
Tools for tau triggers validation for 7_3_X

cmsrel CMSSW_7_3_6_patch1; #or similar 73X
cd CMSSW_7_3_6_patch1/src;
cmsenv;

git clone https://github.com/mbluj/usercode-TriggerStudies.git TriggerStudies/Tau;
(cd TriggerStudies/Tau; git checkout 7_3_X)

scram b


