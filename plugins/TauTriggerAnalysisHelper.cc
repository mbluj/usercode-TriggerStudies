#include "TriggerStudies/Tau/interface/TauTriggerAnalysisHelper.h"

void tautrigger::addElIds(std::map<std::string,CutBasedElId> &ids){

  //PHYS14_25ns WP's: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Working_points_for_PHYS14_sample
  ids["POG_PHYS14_25ns_v1_Veto"] = CutBasedElId("POG_PHYS14_25ns_v1_Veto",
						0.016315, 0.010671,
						0.252044, 0.245263,
						0.011100, 0.033987,
						0.345843, 0.134691,
						0.248070, 0.157160);
  ids["POG_PHYS14_25ns_v1_Loose"] = CutBasedElId("POG_PHYS14_25ns_v1_Loose",
						 0.012442, 0.010654,
						 0.072624, 0.145129,
						 0.010557, 0.032602,
						 0.121476, 0.131862,
						 0.221803, 0.142283);
  ids["POG_PHYS14_25ns_v1_Medium"] = CutBasedElId("POG_PHYS14_25ns_v1_Medium",
						  0.007641, 0.009285,
						  0.032643, 0.042447,
						  0.010399, 0.029524,
						  0.060662, 0.104263,
						  0.153897, 0.137468);
  ids["POG_PHYS14_25ns_v1_Tight"] = CutBasedElId("POG_PHYS14_25ns_v1_Tight",
						 0.006574, 0.005681,
						 0.022868, 0.032046,
						 0.010181, 0.028766,
						 0.037553, 0.081902,
						 0.131191, 0.106055);
  return;
}
