#include "MODEL/Main/Coupling_Data.H"

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace MODEL;
using namespace ATOOLS;

void Coupling_Data::Calculate()
{
  if (p_scl==NULL) return;
  m_fac=1.;
  if (*p_scl>=0.) m_fac=(*p_cpl)(*p_scl)/m_def;
  msg_Debugging()<<METHOD<<": scl = "<<sqrt(*p_scl)<<" -> "<<*this<<"\n";
}

void Coupling_Map::Calculate() const
{
  for (const_iterator cit(begin());
       cit!=end();++cit) cit->second->Calculate();
}

Coupling_Data *Coupling_Map::Get
(const std::string &tag,const ATOOLS::NLO_subevt *sub) const
{
  for (const_iterator cit(lower_bound(tag)),
	 eit(upper_bound(tag));cit!=eit;++cit)
    if (cit->second->Sub()==sub) return cit->second;
  return NULL;
}

namespace MODEL {

  std::ostream &operator<<(std::ostream &str,const Coupling_Data &cd)
  {
    str<<"'"<<cd.ID()<<"'";
    if (cd.Sub()) str<<"[("<<cd.Sub()->m_i<<","<<cd.Sub()->m_j
		     <<")("<<cd.Sub()->m_k<<")]";
    return str<<"{fac="<<cd.Factor()<<",cpl="<<cd.Default()<<"}";
  }

  std::ostream &operator<<(std::ostream &str,const Coupling_Map &cm)
  {
    str<<&cm<<" -> "<<cm.size()<<" couplings stored"<<std::endl;
    for (Coupling_Map::const_iterator it(cm.begin());it!=cm.end();++it) {
      str<<"  "<<it->first<<": "<<*it->second<<std::endl;
    }
    return str;
  }

}
