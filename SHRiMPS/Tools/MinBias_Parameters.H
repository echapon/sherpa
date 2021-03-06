#ifndef SHRIMPS_Tools_MinBias_Parameters_H
#define SHRIMPS_Tools_MinBias_Parameters_H

#include "SHRiMPS/Tools/Parameter_Structures.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  class Form_Factor;
  class Omega_ik;
  class Cross_Sections;

  struct XSecs_Container {
    double xs_tot, xs_in, xs_el, xs_SD, xs_DD;
    XSecs_Container():
      xs_tot(0.), xs_in(0.), xs_el(0.), xs_SD(0.), xs_DD(0.) {}
  };

  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class MinBias_Parameters {
  private:
    run_mode::code        m_runmode;
    size_t                m_NGWstates;
    double                m_originalY, m_bmax, m_accu;
    FormFactor_Parameters m_ff_params;
    Eikonal_Parameters    m_eik_params;
    Run_Parameters        m_run_params;
    Ladder_Parameters     m_ladder_params; 
    ShowerLink_Parameters m_showerlink_params; 

    std::list<Form_Factor *> * p_ffs;
    std::list<Omega_ik *>    * p_eikonals;
    XSecs_Container          * p_xsecs;

    void RegisterDefaults() const;
    void FillRunParameters();
    void FillFormFactorParameters();
    void FillEikonalParameters();
    void FillLadderParameters();
    void FillShowerLinkParameters();
  public:
    MinBias_Parameters();
    ~MinBias_Parameters();

    void Init();
    void Reset();
    void AddFormFactor(Form_Factor * ff) { p_ffs->push_back(ff); }
    void AddEikonal(Omega_ik * eikonal)  { p_eikonals->push_back(eikonal); }
    void SetXSecs(Cross_Sections * xsecs);
    void UpdateForNewEnergy(const double & energy);

    const run_mode::code & RunMode()   const { return m_runmode; }
    const size_t         & NGWStates() const { return m_NGWstates; }

    const Run_Parameters & GetRunParameters() const {
      return m_run_params; 
    }
    const FormFactor_Parameters & GetFFParameters() const {
      return m_ff_params; 
    }
    const Eikonal_Parameters & GetEikonalParameters() const {
      return m_eik_params; 
    }
    const Ladder_Parameters & GetLadderParameters() const {
      return m_ladder_params; 
    }
    const ShowerLink_Parameters & GetShowerLinkParameters() const {
      return m_showerlink_params; 
    }
    std::list<Form_Factor *> * GetFormFactors() const {
      return p_ffs;
    }
    std::list<Omega_ik *> * GetEikonals() const {
      return p_eikonals;
    }
    XSecs_Container * GetXSecs() const {
      return p_xsecs;
    }
  };
  
  extern MinBias_Parameters MBpars;
}

#endif
