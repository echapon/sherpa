#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ANALYSIS;

class Blob_Selector: public Analysis_Object {
private:
  ATOOLS::btp::code m_type;
  std::string m_outlist;
  int m_mode;
public:
  Blob_Selector(const int type,const int mode,const std::string &outlist):
    m_type((ATOOLS::btp::code)type), m_outlist(outlist), m_mode(mode) {}
  
  Analysis_Object *GetCopy() const 
  {
    return new Blob_Selector(m_type,m_mode,m_outlist);
  }

  void Evaluate(const ATOOLS::Blob_List &blobs,double value,double ncount)
  {
    ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
    p_ana->AddParticleList(m_outlist,outlist);
    {
      for (ATOOLS::Blob_List::const_iterator bit=blobs.begin();
	   bit!=blobs.end();++bit) {
	if ((*bit)->Type()&m_type) {
	  if (m_mode>1) {
	    for (int i=0;i<(*bit)->NInP();++i) {
	      const ATOOLS::Particle *cur=(*bit)->ConstInParticle(i);
	      bool present=false;
	      for (ATOOLS::Particle_List::const_iterator pit=outlist->begin();
		   pit!=outlist->end();++pit) 
		if (cur==*pit) {
		  present=true;
		  break;
		}
	      if (!present) outlist->push_back(new ATOOLS::Particle(*cur));
	    }
	  }
	  for (int i=0;i<(*bit)->NOutP();++i) {
	    const ATOOLS::Particle *cur=(*bit)->ConstOutParticle(i);
	    if (cur->DecayBlob()!=NULL) if (m_mode<1) continue;
	    bool present=false;
	    for (ATOOLS::Particle_List::const_iterator pit=outlist->begin();
		 pit!=outlist->end();++pit) 
	      if (cur==*pit) {
		present=true;
		break;
	      }
	    if (!present) outlist->push_back(new ATOOLS::Particle(*cur));
	  }
	}
      }
    }
  }

};

DECLARE_GETTER(Blob_Selector,"BlobSel",
	       Analysis_Object,Analysis_Key);

Analysis_Object *ATOOLS::Getter<Analysis_Object,Analysis_Key,Blob_Selector>::
operator()(const Analysis_Key& key) const
{ 
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  const auto outlist = parameters.size() > 2 ? parameters[2] : "Analysed";
  const auto type = s.Interprete<int>(parameters[0]);
  const auto mode = s.Interprete<int>(parameters[1]);
  return new Blob_Selector(type, mode, outlist);
}

void ATOOLS::Getter<Analysis_Object,Analysis_Key,Blob_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[type, mode, outlist (optionally)]";
}

