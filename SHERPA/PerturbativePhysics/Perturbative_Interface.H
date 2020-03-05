#ifndef SHERPA_PerturbativePhysics_Perturbative_Interface_H
#define SHERPA_PerturbativePhysics_Perturbative_Interface_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"
#include "SHERPA/Tools/Variations.H"

namespace ATOOLS { class Cluster_Amplitude; }

namespace SHERPA {

  class Matrix_Element_Handler;
  class MI_Handler;
  class Soft_Collision_Handler;
  class Shower_Handler;

  class Perturbative_Interface {
  protected :

    Matrix_Element_Handler *p_me;
    MI_Handler             *p_mi;
    Soft_Collision_Handler *p_sc;
    Shower_Handler         *p_shower;        

    ATOOLS::Blob *p_hard;

    ATOOLS::Poincare m_cms;

    ATOOLS::Cluster_Amplitude *p_ampl;
    ATOOLS::Mass_Selector *p_ms;


    double m_weight, m_globalkfac, m_maxkfac;
    int    m_cmode, m_bbarmode;

    SHERPA::Variation_Weights *p_localkfactorvarweights;

    bool LocalKFactor(ATOOLS::Cluster_Amplitude* ampl);

  public :

    // constructors
    Perturbative_Interface(Matrix_Element_Handler *const meh,
			   Shower_Handler *const psh);
    Perturbative_Interface(Shower_Handler *const psh);
    Perturbative_Interface(MI_Handler *const mih,
			   Shower_Handler *const psh);
    Perturbative_Interface(Soft_Collision_Handler *const sch,
			   Shower_Handler *const psh);

    // destructor
    ~Perturbative_Interface();

    // member functions
    ATOOLS::Return_Value::code DefineInitialConditions(ATOOLS::Blob *blob);

    bool FillBlobs(ATOOLS::Blob_List *blobs); 

    int  PerformShowers();
    int  PerformDecayShowers();

    void CleanUp();
    
    // inline functions
    inline double Weight() const { return m_weight; }
    
    inline Matrix_Element_Handler *MEHandler() { return p_me; }

    ATOOLS::Cluster_Amplitude *Amplitude() const { return p_ampl; }

    Shower_Handler* Shower() { return p_shower; }

  };// end of class Perturbative_Interface

  typedef std::map<std::string,Perturbative_Interface *> PertInterfaceMap;
  typedef PertInterfaceMap::iterator                     PertInterfaceIter;

}// end of namespace SHERPA

#endif
