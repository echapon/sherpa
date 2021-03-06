#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class List_Merger: public Analysis_Object {
  protected:

    std::vector<std::string> m_inlists;

    std::string m_outlist;

  public:

    List_Merger(const std::vector<std::string> &inlists,
		const std::string &outlist);
    
    void Evaluate(const ATOOLS::Blob_List &blobs,
		  double value,double ncount);

    Analysis_Object *GetCopy() const;
    
  };// List_Merger

}// end of namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(List_Merger,"MergeLists",
	       Analysis_Object,Analysis_Key);

Analysis_Object *ATOOLS::Getter<Analysis_Object,Analysis_Key,List_Merger>::
operator()(const Analysis_Key& key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() < 3)
    THROW(missing_input, "MergeLists expects at least three parameters.");
  auto inlists(parameters);
  inlists.pop_back();
  return new List_Merger(inlists, parameters.back());
}

void ATOOLS::Getter<Analysis_Object,Analysis_Key,List_Merger>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"[inlist1, .., inlistN, outlist]"; }

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

List_Merger::List_Merger(const std::vector<std::string> &inlists,
			 const std::string &outlist):
  m_inlists(inlists), m_outlist(outlist) {}

void List_Merger::Evaluate(const ATOOLS::Blob_List &blobs,
			   double value,double ncount)
{
  Particle_List *outlist(new Particle_List());
  for (size_t i(0);i<m_inlists.size();++i) {
    const Particle_List *inlist(p_ana->GetParticleList(m_inlists[i]));
    if (inlist==NULL) {
      msg_Error()<<METHOD<<"(): List "<<i<<" '"<<m_inlists[i]
		 <<"' not found."<<std::endl;
    }
    else {
      size_t k(outlist->size());
      outlist->resize(outlist->size()+inlist->size());
      for (size_t j(0);j<inlist->size();++j)
	(*outlist)[k++] = new Particle(*(*inlist)[j]);
    }
  }
  p_ana->AddParticleList(m_outlist,outlist);
}

Analysis_Object *List_Merger::GetCopy() const
{
  return new List_Merger(m_inlists,m_outlist);
}

