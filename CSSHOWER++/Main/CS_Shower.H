#ifndef CSSHOWER_Main_CS_Shower_H
#define CSSHOWER_Main_CS_Shower_H

#include "PDF/Main/Shower_Base.H"
#include "CSSHOWER++/Main/CS_Cluster_Definitions.H"
#include "CSSHOWER++/Showers/Shower.H"
#include "ATOOLS/Phys/Blob_List.H"


namespace REMNANTS { class Remnant_Handler; }

namespace CSSHOWER {
  struct partcomp{
    bool operator() (const Parton * part1, const Parton * part2) const
    {
      if (part1->Id() < part2->Id()) return true;
      else return false;
    }
  };

  class Splitting_Function_Base;

  typedef std::map<size_t,std::pair<double,double> > KT2X_Map;

  class CS_Shower : public PDF::Shower_Base {

  private:

    struct Q2_Value {
      double m_q2;
      ATOOLS::Flavour m_fl;
      int m_i, m_j, m_k;
      inline Q2_Value(const double &q2,
		      const ATOOLS::Flavour &fl,
		      int i, int j, int k):
	m_q2(q2), m_fl(fl), m_i(i), m_j(j), m_k(k) {}
      inline bool operator<(const Q2_Value &v) const
      { return m_q2<v.m_q2; }
    };

    PDF::ISR_Handler * p_isr;
    int m_kmode, m_type;
    bool m_respectq2;
    size_t m_maxem, m_nem, m_recocheck;

    Shower          * p_shower;
    All_Singlets m_allsinglets;
    CS_Cluster_Definitions *p_cluster;
    All_Singlets *p_next, *p_refs;

    ATOOLS::Mass_Selector *p_ms;

    ATOOLS::Cluster_Amplitude *p_rampl;

    void   GetKT2Min(ATOOLS::Cluster_Amplitude *const ampl,const size_t &id,
		     KT2X_Map &kt2xmap,std::set<size_t> &aset);
    void   GetKT2Min(ATOOLS::Cluster_Amplitude *const ampl,KT2X_Map &kt2xmap);

    int IsDecay(ATOOLS::Cluster_Amplitude *const ampl,
		ATOOLS::Cluster_Leg *const cl) const;

    double HardScale(const ATOOLS::Cluster_Amplitude *const ampl);

    Singlet *TranslateAmplitude(ATOOLS::Cluster_Amplitude *const ampl,
				std::map<ATOOLS::Cluster_Leg*,Parton*> &pmap,
				std::map<Parton*,ATOOLS::Cluster_Leg*> &lmap,
				const KT2X_Map &kt2xmap);

    bool PrepareStandardShower(ATOOLS::Cluster_Amplitude *const ampl);
    int  PerformShowers(const size_t &maxem,size_t &nem);

    double Qij2(const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
		const ATOOLS::Vec4D &pk,const ATOOLS::Flavour &fi,
		const ATOOLS::Flavour &fj) const;

  public:

    // constructor 
    CS_Shower(PDF::ISR_Handler *const,
              MODEL::Model_Base *const,
              const int type);

    // destructor
    ~CS_Shower();

    //member functions
    int  PerformShowers();
    int  PerformDecayShowers();

    void AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl);

    bool ExtractPartons(ATOOLS::Blob_List *const blist);

    void CleanUp();

    // inline functions
    PDF::Cluster_Definitions_Base * GetClusterDefinitions();
    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl,
		       const bool & soft=false);
    double JetVeto(ATOOLS::Cluster_Amplitude *const ampl,
		   const int mode);

    void SetRBOff();
    inline void SetRemnants(REMNANTS::Remnant_Handler * remnants) {
      p_shower->SetRemnants(remnants);
    }

  };// end of class CS_Shower

}// end of namespace CSSHOWER

#endif
