#ifndef TriggerStudies_Tau_Helpers_h
#define TriggerStudies_Tau_Helpers_h

// Helper functions
// M. Bluj

// system include files
#include <memory>
#include <cmath>

//
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <utility>


namespace tautrigger {

  struct CutBasedElId{
    CutBasedElId(){};

    CutBasedElId(std::string name,
		 std::pair<float,float> dEtaInCut,
		 std::pair<float,float> dPhiInCut,
		 std::pair<float,float> sigmaIEtaIEtaCut,
		 std::pair<float,float> HoECut,
		 std::pair<float,float> EoPCut):
      name_(name),
      dEtaInCut_(dEtaInCut),dPhiInCut_(dPhiInCut),sigmaIEtaIEtaCut_(sigmaIEtaIEtaCut),
      HoECut_(HoECut),EoPCut_(EoPCut)
    {};
    CutBasedElId(std::string name,
		 float dEtaInCutEB, float dEtaInCutEE,
		 float dPhiInCutEB, float dPhiInCutEE,
		 float sigmaIEtaIEtaCutEB, float sigmaIEtaIEtaCutEE,
		 float HoECutEB, float HoECutEE,
		 float EoPCutEB, float EoPCutEE):
      name_(name),
      dEtaInCut_(dEtaInCutEB,dEtaInCutEE),dPhiInCut_(dPhiInCutEB,dPhiInCutEE),
      sigmaIEtaIEtaCut_(sigmaIEtaIEtaCutEB,sigmaIEtaIEtaCutEE),
      HoECut_(HoECutEB,HoECutEE),EoPCut_(EoPCutEB,EoPCutEE)
    {};
    ~CutBasedElId(){};
    std::string name_;
    std::pair<float,float> dEtaInCut_,dPhiInCut_,sigmaIEtaIEtaCut_,HoECut_,EoPCut_;
    void setName(std::string name) { name_=name; }
    bool checkIdEB(float dEtaIn,
		   float dPhiIn,
		   float sigmaIEtaIEta,
		   float HoE,
		   float EoP){
      if(dEtaIn >= dEtaInCut_.first) return false;
      if(dPhiIn >= dPhiInCut_.first) return false;
      if(sigmaIEtaIEta >= sigmaIEtaIEtaCut_.first) return false;
      if(HoE >= HoECut_.first) return false;
      if(EoP >= EoPCut_.first) return false;
      return true;
    }
    bool checkIdEE(float dEtaIn,
		   float dPhiIn,
		   float sigmaIEtaIEta,
		   float HoE,
		   float EoP){
      if(dEtaIn >= dEtaInCut_.second) return false;
      if(dPhiIn >= dPhiInCut_.second) return false;
      if(sigmaIEtaIEta >= sigmaIEtaIEtaCut_.second) return false;
      if(HoE >= HoECut_.second) return false;
      if(EoP >= EoPCut_.second) return false;
      return true;
    }
    bool checkId(float dEtaIn,
		 float dPhiIn,
		 float sigmaIEtaIEta,
		 float HoE,
		 float EoP,
		 bool isEB){
      if(isEB) return checkIdEB(dEtaIn,dPhiIn,sigmaIEtaIEta,HoE,EoP);
      else return checkIdEE(dEtaIn,dPhiIn,sigmaIEtaIEta,HoE,EoP);
    }
  };

  void addElIds(std::map<std::string,CutBasedElId> &ids);

}

#endif
