#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class QCD_Core_Scale: public Core_Scale_Setter {
  public:

    QCD_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::Cluster_Param Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::Cluster_Param QCD_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double s(2.0*ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom());
  double t(2.0*ampl->Leg(0)->Mom()*ampl->Leg(2)->Mom());
  double u(2.0*ampl->Leg(0)->Mom()*ampl->Leg(3)->Mom());
  double muf2(-1.0/(1.0/s+1.0/t+1.0/u)/sqrt(2.0)), mur2(muf2), q2(muf2);
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"  \\mu_q = "<<sqrt(q2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::Cluster_Param(NULL,q2,muf2,mur2,-1);
}

DECLARE_ND_GETTER(QCD_Core_Scale,"QCD",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,QCD_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new QCD_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    QCD_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"QCD core scale"; 
}
