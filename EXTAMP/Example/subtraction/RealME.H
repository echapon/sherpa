#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "ATOOLS/Phys/Flavour.H"

#include "CPPProcess_real.h"

namespace PHASIC {
  class Process_Info;
}


class RealME: public PHASIC::Tree_ME2_Base {

public:
  
  RealME(const PHASIC::Process_Info &);
    
    ~RealME() {}
  
  double Calc(const ATOOLS::Vec4D_Vector& momenta);

  int OrderQCD(const int &id=-1) const { return m_oqcd; }
  int OrderEW (const int &id=-1) const { return m_oew;  }

  double RunningCouplingFactor() const;
  
private:

  size_t m_nin, m_nout;
  size_t m_oew,m_oqcd;

  PHASIC::Process_Info m_pinfo;

  double m_symfac;

  MG_Real::CPPProcess m_mg_proc;
};
