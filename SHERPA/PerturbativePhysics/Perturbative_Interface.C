#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"

#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Decays/Decay_Handler.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;

bool DecayDefineInitialConditions(Cluster_Amplitude* ampl,
                                  Blob* initial_blob,
                                  PDF::Cluster_Definitions_Base *p_clus);
Cluster_Amplitude* DecayClusterConfiguration(Blob *const bl);
class Had_Mass_Selector : public ATOOLS::Mass_Selector {
public:
  inline double Mass(const ATOOLS::Flavour &fl) const { return fl.HadMass(); }

  inline double Mass2(const Flavour &fl) const
  { double m(Mass(fl)); return m*m; }

};


Perturbative_Interface::Perturbative_Interface
(Matrix_Element_Handler *const meh,Shower_Handler *const psh):
  p_me(meh), p_mi(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL),
  p_ms(NULL)
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.SetInputPath(p_me->Path());
  read.SetInputFile(p_me->File());
  m_cmode=ToType<int>(rpa->gen.Variable("METS_CLUSTER_MODE"));
  m_bbarmode=read.GetValue<int>("METS_BBAR_MODE",1);
  m_globalkfac=read.GetValue<double>("GLOBAL_KFAC",0.);
  m_maxkfac=read.GetValue<double>("MENLOPS_MAX_KFAC",10.0);
}

Perturbative_Interface::Perturbative_Interface
(MI_Handler *const mi,Shower_Handler *const psh):
  p_me(NULL), p_mi(mi), p_sc(NULL), p_shower(psh), p_ms(NULL),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_sc(NULL), p_shower(psh), p_ms(NULL),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Soft_Collision_Handler *const sch,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_sc(sch), p_shower(psh), p_ms(NULL),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL)  {}

Perturbative_Interface::~Perturbative_Interface() 
{
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
  }
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
  }
}

Return_Value::code Perturbative_Interface::
DefineInitialConditions(ATOOLS::Blob *blob) 
{
  if (blob==NULL) {
    msg_Error()<<METHOD<<"(): Signal process not found."<<std::endl;
    return Return_Value::Error;
  }
  p_hard=blob;
  if (!p_shower->On() || (p_me && p_me->Process()->Info().m_nlomode==1)) {
    m_weight=1.0;
    return Return_Value::Success;
  }
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
    p_ampl=NULL;
  }
  p_shower->CleanUp();
  DEBUG_FUNC(blob->Id());
  if (p_mi) {
    p_ampl=p_mi->ClusterConfiguration();
    p_mi->Process()->Generator()->SetMassMode(1);
    int stat(p_mi->Process()->Generator()->ShiftMasses(p_ampl));
    if (stat<0) {
      msg_Tracking()<<METHOD<<"(): MI Mass shift failed. Reject event."<<std::endl;
      return Return_Value::Retry_Event;
    }
    if (stat==1) {
      stat=p_mi->Shower()->GetShower()->
	GetClusterDefinitions()->ReCluster(p_ampl);
      if (stat!=1) {
	msg_Tracking()<<METHOD<<"(): MI Reclustering failed. Reject event.\n";
	return Return_Value::Retry_Event;
      }
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl))
      return Return_Value::New_Event;
    return Return_Value::Success;
  }
  if (blob->Type()==btp::Hadron_Decay) {
    if (p_ms) delete p_ms; p_ms=NULL;
    p_ms=new Had_Mass_Selector(); // TODO find proper location
    p_ampl=DecayClusterConfiguration(blob);
    p_ampl->SetMS(p_ms);
    if (!p_shower->GetShower()->PrepareShower(p_ampl))
      return Return_Value::New_Event;
    p_ampl->SetMS(NULL);
    return Return_Value::Success;
  }
  if (p_sc) {
    p_sc->SetClusterDefinitions(p_shower->GetShower()->GetClusterDefinitions());
    p_ampl=p_sc->ClusterConfiguration(blob);
    if (p_ampl==NULL) {
      msg_Out()<<METHOD<<": Soft_Collision_Handler has no amplitude.\n";
      return Return_Value::New_Event;
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl,true)) {
      msg_Out()<<METHOD<<": could not prepare shower.\n"; 
      return Return_Value::New_Event;
    }
    return Return_Value::Success;
  }
  p_ampl=p_me->Process()->Get<Single_Process>()->Cluster
    (p_me->Process()->Integrator()->Momenta(),m_cmode);
  if (p_ampl==NULL) return Return_Value::New_Event;
  if (p_ampl->MS()==NULL)
    p_ampl=p_me->Process()->Get<Single_Process>()->Cluster
      (p_me->Process()->Integrator()->Momenta(),m_cmode|256);
  if (p_ampl==NULL) return Return_Value::New_Event;
  m_weight=1.0;
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
    p_localkfactorvarweights = NULL;
  }
  if (p_me->Process()->Info().m_ckkw&1) {
    if ((m_bbarmode&1) && p_me->HasNLO() &&
        p_me->Process()->Parent()->Info().m_fi.NLOType()==nlo_type::lo) {
      Cluster_Amplitude *oampl=p_me->Process()->
	Get<Single_Process>()->Cluster
	(p_me->Process()->Integrator()->Momenta(),m_cmode);
      if (!LocalKFactor(oampl)) {
	DEBUG_INFO("didn't find process using original amplitude");
	if (m_bbarmode&4) {
	  Cluster_Amplitude *ampl=p_me->Process()->
	    Get<Single_Process>()->Cluster
	    (p_me->Process()->Integrator()->Momenta(),m_cmode|16|256|512);
	  while (ampl->Prev()) ampl=ampl->Prev();
	  if (!LocalKFactor(ampl))
	    DEBUG_INFO("didn't find process using exclusive clustering");
	  ampl->Delete();
	}
      }
      while (oampl->Prev()) oampl=oampl->Prev();
      oampl->Delete();
    }
  }
  p_me->Process()->Generator()->SetMassMode(1);
  int stat(p_me->Process()->Generator()->ShiftMasses(p_ampl));
  if (stat<0) {
    msg_Tracking()<<METHOD<<"(): ME Mass shift failed. Reject event."<<std::endl;
    return Return_Value::New_Event;
  }
  if (stat==1) {
    stat=p_me->Shower()->GetShower()->
      GetClusterDefinitions()->ReCluster(p_ampl);
    if (stat!=1) {
      msg_Tracking()<<METHOD<<"(): ME Reclustering failed. Reject event."<<std::endl;
      return Return_Value::New_Event;
    }
  }
  size_t cmax(0);
  for (size_t i(0);i<p_ampl->Legs().size();++i)
    cmax=Max(cmax,(size_t)p_ampl->Leg(i)->Col().m_i);
  while (Flow::Counter()<cmax);
  p_me->Process()->Parent()->SetRBMap(p_ampl);
  if (/*p_dec*/true) {
    if (!DecayDefineInitialConditions(p_ampl, blob, p_me->Shower()->GetShower()->GetClusterDefinitions())) {
      msg_Tracking()<<METHOD<<"(): Decay clustering failed. Reject event."<<std::endl;
      return Return_Value::Retry_Event;
    }
    Cluster_Amplitude *ampl(p_ampl);
    while (ampl->Prev()) ampl=ampl->Prev();
    int stat(p_me->Process()->Generator()->ShiftMasses(ampl));
    if (stat<0) {
      msg_Tracking()<<METHOD<<"(): DH Mass shift failed. Reject event."<<std::endl;
      return Return_Value::Retry_Event;
    }
    if (stat==1) {
      stat=p_me->Shower()->GetShower()->
	GetClusterDefinitions()->ReCluster(ampl);
      if (stat!=1) {
	msg_Tracking()<<METHOD<<"(): DH Reclustering failed. Reject event."<<std::endl;
	return Return_Value::Retry_Event;
      }
    }
  }
  while (p_ampl->Prev()) p_ampl=p_ampl->Prev();
  if (p_me->Process()->Info().m_ckkw&1) {
    blob->AddData("Sud_Weight",new Blob_Data<double>(m_weight));
    if (p_me->EventGenerationMode()!=0) {
      const double disc = ran->Get();
      const double abswgt = std::abs(m_weight);
      if (abswgt < disc) {
        return Return_Value::New_Event;
      }
      m_weight /= Min(1.0, abswgt);
      // local kfactor varweights not initialized if weight equal to one
      if (p_localkfactorvarweights && abswgt < 1.)
        *p_localkfactorvarweights *= 1.0 / abswgt;
    }
    Blob_Data_Base *winfo((*blob)["Weight"]);
    if (!winfo) THROW(fatal_error,"No weight information in signal blob");
    double meweight(winfo->Get<double>());
    blob->AddData("Weight",new Blob_Data<double>(meweight*m_weight));
    // also update reweighting weights
    Blob_Data_Base *vws((*blob)["Variation_Weights"]);
    if (vws) {
      if (p_localkfactorvarweights) {
        vws->Get<Variation_Weights>() *= *p_localkfactorvarweights;
      } else {
        vws->Get<Variation_Weights>() *= m_weight;
      }
    }
  }
  if (!p_shower->GetShower()->PrepareShower(p_ampl)) 
    return Return_Value::New_Event;
  return Return_Value::Success;
}

bool Perturbative_Interface::LocalKFactor(ATOOLS::Cluster_Amplitude* ampl)
{
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
    p_localkfactorvarweights = NULL;
  }
  if (m_globalkfac) {
    m_weight*=m_globalkfac;
    return true;
  }
  DEBUG_FUNC(ampl->Legs().size());
  Process_Vector procs(p_me->AllProcesses());
  Process_Base::SortFlavours(ampl);
  if (p_hard) {
    Blob_Data_Base *vws((*p_hard)["Variation_Weights"]);
    if (vws) {
      if (p_localkfactorvarweights) {
        delete p_localkfactorvarweights;
      }
      Variations *variations = vws->Get<Variation_Weights>().GetVariations();
      p_localkfactorvarweights = new SHERPA::Variation_Weights(variations);
    }
  }
  while (ampl->Next()!=NULL) {
    ampl=ampl->Next();
    if (ampl->Next() && (m_bbarmode&2)) continue;
    Process_Base::SortFlavours(ampl);
    for (size_t i=0; i<procs.size(); ++i) {
      if (p_localkfactorvarweights) p_localkfactorvarweights->Reset();
      MCatNLO_Process* mcnloproc=dynamic_cast<MCatNLO_Process*>(procs[i]);
      if (mcnloproc) {
        if (mcnloproc->VariationWeights()) THROW(fatal_error, "Variation weights already set.");
        mcnloproc->SetVariationWeights(p_localkfactorvarweights);
	double K(mcnloproc->LocalKFactor(*ampl));
        mcnloproc->SetVariationWeights(NULL);
	if (K==0.0 || dabs(K)>m_maxkfac) continue;
	m_weight*=K;
	return true;
      }
    }
  }
  // no process found along ampl
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
    p_localkfactorvarweights = NULL;
  }
  return false;
}

bool Perturbative_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
{
  if (p_hard==NULL) return false;
  Blob *sblob = new Blob();
  sblob->SetType(btp::Shower);
  sblob->SetStatus(blob_status::needs_showers);
  sblob->SetId();
  sblob->SetPosition(p_hard->Position());
  if (p_shower->On()) {
    if (p_hard->Type()!=btp::Hadron_Decay) {
      for (int i(0);i<p_hard->NInP();++i)
	sblob->AddToOutParticles(p_hard->InParticle(i));
      for (size_t j(0);j<blobs->size();++j) {
        Blob *cb((*blobs)[j]);
        if (cb->Has(blob_status::needs_showers))
          for (int i(0);i<cb->NOutP();++i)
            if (cb->OutParticle(i)->DecayBlob()==NULL)
              sblob->AddToInParticles(cb->OutParticle(i));
      }
    }
    else {
      for (int i(0);i<p_hard->NOutP();++i) {
	if (!(p_hard->OutParticle(i)->GetFlow(1)==0 &&
	      p_hard->OutParticle(i)->GetFlow(2)==0))
	  sblob->AddToInParticles(p_hard->OutParticle(i));
      }
    }
  }
  blobs->push_back(sblob);
  p_shower->FillBlobs(blobs); 
  return true;
}

int Perturbative_Interface::PerformShowers()
{
  DEBUG_FUNC("");
  // see if the event has any weight
  Blob_Data_Base *winfo((*p_hard)["Weight"]);
  if (!winfo) THROW(fatal_error,"No weight information in signal blob");
  double meweight(winfo->Get<double>());
  if (meweight==0.0) return 0;

  PDF::Shower_Base *csh(p_shower->GetShower());

  // look for reweightings and set up the shower accordingly
  Blob_Data_Base *blob_data_base((*p_hard)["Variation_Weights"]);
  if (blob_data_base) {
    csh->SetVariationWeights(&blob_data_base->Get<Variation_Weights>());
  }

  int stat=csh->PerformShowers();
  double weight=csh->Weight();
  m_weight*=weight;
  p_hard->AddData("Shower_Weight",new Blob_Data<double>(weight));
  p_hard->AddData("Weight",new Blob_Data<double>(meweight*weight));
  if (blob_data_base) {
    blob_data_base->Get<Variation_Weights>() *= weight;
  }
  return stat;
}

int Perturbative_Interface::PerformDecayShowers()
{ 
  DEBUG_FUNC("");
  return p_shower->GetShower()->PerformDecayShowers(); 
}

void Perturbative_Interface::CleanUp()
{
  if (p_me && p_me->Process())
    p_me->Process()->Generator()->SetMassMode(0);
  if (p_mi && p_mi->Process())
    p_mi->Process()->Generator()->SetMassMode(0);
  p_shower->CleanUp();
}













// TODO: move to proper place
#include <algorithm>
using namespace std;
typedef std::pair<ATOOLS::Particle *,ATOOLS::Particle *> ParticlePair;
typedef std::vector<ParticlePair> ParticlePair_Vector;

class ParticlePairFirstEnergySort {
public:
  bool operator()(const ParticlePair& a,const ParticlePair& b)
  { return (a.first->Momentum()[0]<b.first->Momentum()[0]); }
};


void AssignPhotons(const Particle_Vector& daughters,
                                       ParticlePair_Vector& photons)
{
  // for every photon, find charged particle that's closest
  // ignore radiation off charged resonance for now
  if (photons.size()) {
    Particle_Vector cdaughters;
    for (size_t i(0);i<daughters.size();++i)
      if (daughters[i]->Flav().Charge()) cdaughters.push_back(daughters[i]);
    if (cdaughters.size()==1) {
      for (size_t i(0);i<photons.size();++i)
        photons[i].second=cdaughters[0];
    }
    else {
      Vec4D cmom(0.,0.,0.,0.);
      Vec4D_Vector cmoms;
      for (size_t i(0);i<cdaughters.size();++i) {
        cmoms.push_back(cdaughters[i]->Momentum());
        cmom+=cmoms[i];
      }
      Poincare ccms(cmom);
      for (size_t i(0);i<cdaughters.size();++i) ccms.Boost(cmoms[i]);
      for (size_t i(0);i<photons.size();++i){
        Vec4D pmom(photons[i].first->Momentum());
        ccms.Boost(pmom);
        size_t id(0);
        double dR(pmom.DR(cmoms[0]));
        for (size_t j(1);j<cmoms.size();++j) {
          double dRj(pmom.DR(cmoms[j]));
          if (dRj<dR) { id=j; dR=dRj; }
        }
        photons[i].second=cdaughters[id];
      }
    }
    for (size_t i(0);i<photons.size();++i) {
      if (photons[i].first==photons[i].second)
        THROW(fatal_error,"Photon has not been assigned.");
      msg_Debugging()<<photons[i].first->Flav()<<" "
                     <<photons[i].first->Momentum()
                     <<" assigned to "<<photons[i].second->Flav()<<std::endl;
    }
  }
}

Vec4D RecombinedMomentum(const Particle * daughter,
                                             const ParticlePair_Vector& photons,
                                             size_t& stat)
{
  Vec4D mom(0.,0.,0.,0.);
  for (size_t i(0);i<photons.size();++i) {
    if (photons[i].second==daughter) {
      mom+=photons[i].first->Momentum();
      stat|=2|4;
    }
  }
  msg_Debugging()<<daughter->Flav()<<": "<<mom<<" "<<stat<<std::endl;
  return mom+daughter->Momentum();
}


void AddPhotonsClustering(Cluster_Amplitude*& ampl,
                                              const Particle_Vector daughters,
                                              ParticlePair_Vector& photons,
                                              size_t& imax,
                                              const std::vector<size_t>& ids)
{
  DEBUG_FUNC(photons.size()<<" photons to be clustered");
  Particle * photon(photons.back().first);
  Particle * daughter(photons.back().second);
  photons.pop_back();
  size_t idmother(0);
  if      (daughter==daughters[0]) idmother=ids[0];
  else if (daughter==daughters[1]) idmother=ids[1];
  else if (daughter==daughters[2]) idmother=ids[2];
  else THROW(fatal_error,"Did not find id for "+daughter->Flav().IDName());
  msg_Debugging()<<"Cluster "<<photon->Flav()<<" "<<photon->Momentum()
                <<" with "<<daughter->Flav()<<" "<<ID(idmother)<<std::endl;
  Cluster_Amplitude* copy=ampl->InitPrev();
  copy->CopyFrom(ampl);
  copy->SetNLO(0);
  copy->SetFlag(1);
  copy->SetMS(ampl->MS());
  Cluster_Leg *lij(ampl->IdLeg(idmother));
  for (size_t i=0; i<ampl->Legs().size(); ++i)
    ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
  lij->SetStat(1|2|4);
  size_t idk(0);
  for (size_t i=0; i<copy->Legs().size(); ++i) {
    copy->Leg(i)->SetK(0);
    if (copy->Leg(i)->Id()!=idmother &&
        (copy->Leg(i)->Col().m_i==lij->Col().m_j ||
         copy->Leg(i)->Col().m_j==lij->Col().m_i))
      idk=copy->Leg(i)->Id();
  }
  if (lij->Col().m_i==0 && lij->Col().m_j==0) {
    // Ad hoc QED partner, must not be another soft photon
    size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
    if (ampl_nout==1) idk=ampl->Leg(0)->Id();
    else {
      size_t select=ampl->Legs().size();
      size_t nvalid=0;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
        if (!(ampl->Leg(i)->Id()&idmother || i>ampl->Legs().size()-1 ||
              ampl->Leg(i)->Flav().Kfcode()==kf_photon)) {
          nvalid++;
        }
      }
      if (nvalid==0) select=0;
      else {
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother ||
                 select>ampl->Legs().size()-1 ||
                 ampl->Leg(select)->Flav().Kfcode()==kf_photon);
      }
      msg_Debugging()<<"choose ("<<ID(ampl->Leg(select)->Id())<<") "
                     <<ampl->Leg(select)->Flav()<<std::endl;
      idk=ampl->Leg(select)->Id();
    }
  }
  else THROW(fatal_error,"Adding QED to coloured particle.");
  if (idk==0) THROW(fatal_error,"Colour partner not found");
  lij->SetK(idk);
  Cluster_Leg *d1(copy->IdLeg(idmother));
  size_t stat1(0), stat2(0);
  d1->SetMom(RecombinedMomentum(daughter,photons,stat1));
  d1->SetStat(stat1);
  d1->SetFlav(daughter->Flav());
  copy->CreateLeg(photon->Momentum(),photon->RefFlav());
  size_t idnew=1<<(++imax);
  copy->Legs().back()->SetId(idnew);
  copy->Legs().back()->SetStat(stat2);
  Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                copy->IdLeg(idmother),
                                copy->Legs().back());
  DEBUG_VAR(*copy);
  Cluster_Amplitude* tmp=copy;
  while (tmp->Next()) {
    tmp=tmp->Next();
    for (size_t i=0; i<tmp->Legs().size(); ++i) {
      if (tmp->Leg(i)->Id()&idmother) {
        tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
      }
      if (tmp->Leg(i)->K()&idmother) {
        tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
      }
    }
    DEBUG_VAR(*tmp);
  }
  ampl=copy;
}

void AddDecayClustering(ATOOLS::Cluster_Amplitude*& ampl,
                                            Blob* blob,
                                            size_t& imax,
                                            size_t idmother)
{
  DEBUG_FUNC("blob->Id()="<<blob->Id()<<" idmother="<<ID(idmother));
  DEBUG_VAR(*blob);
  Particle_Vector daughters;
  ParticlePair_Vector photons;
  for (size_t i(0);i<blob->GetOutParticles().size();++i) {
    Particle * p(blob->OutParticle(i));
    if (p->Info()=='S') photons.push_back(make_pair(p,p));
    else                daughters.push_back(p);
  }
  std::sort(photons.begin(),photons.end(),ParticlePairFirstEnergySort());
  msg_Debugging()<<"daughters: ";
  for (size_t i(0);i<daughters.size();++i)
    msg_Debugging()<<daughters[i]->Flav().IDName()<<" ";
  msg_Debugging()<<" +  "<<photons.size()<<" soft photon(s)"<<std::endl;
  AssignPhotons(daughters,photons);
  if (daughters.size()==2) {
    msg_Debugging()<<"1 to 2 case"<<std::endl;
    Cluster_Amplitude* copy=ampl->InitPrev();
    copy->CopyFrom(ampl);
    copy->SetNLO(0);
    copy->SetFlag(1);
    copy->SetMS(ampl->MS());
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    copy->SetKT2(lij->Mom().Abs2());
    for (size_t i=0; i<ampl->Legs().size(); ++i)
      ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
    lij->SetStat(1|2|4);
    size_t idk(0);
    for (size_t i=0; i<copy->Legs().size(); ++i) {
      copy->Leg(i)->SetK(0);
      if (copy->Leg(i)->Id()!=idmother &&
          (copy->Leg(i)->Col().m_i==lij->Col().m_j ||
           copy->Leg(i)->Col().m_j==lij->Col().m_i))
        idk=copy->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      // Ad hoc EW partner
      size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
      if (ampl_nout==1) idk=ampl->Leg(0)->Id();
      else {
        size_t select=ampl->Legs().size();
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother ||
                 select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(copy->IdLeg(idmother));
    size_t stat1(0), stat2(0);
    d1->SetMom(RecombinedMomentum(daughters[0],photons,stat1));
    d1->SetStat(stat1);
    d1->SetFlav(daughters[0]->Flav());
    d1->SetFromDec(true);
    copy->CreateLeg(RecombinedMomentum(daughters[1],photons,stat2),
                    daughters[1]->RefFlav());
    copy->Legs().back()->SetFromDec(true);
    size_t idnew=1<<(++imax);
    copy->Legs().back()->SetId(idnew);
    copy->Legs().back()->SetStat(stat2);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  copy->IdLeg(idmother),
                                  copy->Legs().back());
    DEBUG_VAR(*copy);
    Cluster_Amplitude* tmp=copy;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
        }
      }
      DEBUG_VAR(*tmp);
    }
    std::vector<size_t> ids;
    ids.push_back(idmother);
    ids.push_back(idnew);
    while (photons.size())
      AddPhotonsClustering(copy, daughters, photons, imax, ids);
    if (daughters[0]->DecayBlob())
      AddDecayClustering(copy, daughters[0]->DecayBlob(), imax, idmother);
    if (daughters[1]->DecayBlob())
      AddDecayClustering(copy, daughters[1]->DecayBlob(), imax, idnew);
    ampl=copy;
  }
  else if (daughters.size()==3) {
    msg_Debugging()<<"1 to 3 case"<<std::endl;
    // structure m -> 0 P[->1 2]
    // propagator always combines daughters 1+2
    Cluster_Amplitude* step1=ampl->InitPrev();
    step1->CopyFrom(ampl);
    step1->SetNLO(0);
    step1->SetFlag(1);
    step1->SetMS(ampl->MS());
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    step1->SetKT2(lij->Mom().Abs2());
    if (!lij) THROW(fatal_error,"Cluster leg of id "+ToString(idmother)
                                +" not found.");
    for (size_t i=0; i<ampl->Legs().size(); ++i)
      ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
    lij->SetStat(1|2|4);
    size_t idk(0);
    for (size_t i=0; i<step1->Legs().size(); ++i) {
      step1->Leg(i)->SetK(0);
      if (step1->Leg(i)->Id()!=idmother)
	if (step1->Leg(i)->Col().m_i==lij->Col().m_j ||
	    step1->Leg(i)->Col().m_j==lij->Col().m_i) 
	  idk=step1->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      // Ad hoc EW partner
      size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
      if (ampl_nout==1) idk=ampl->Leg(0)->Id();
      else {
        size_t select=ampl->Legs().size();
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother || select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(step1->IdLeg(idmother));
    size_t stat1(0),stat2(0),stat3(0);
    d1->SetMom(RecombinedMomentum(daughters[0],photons,stat1));
    d1->SetStat(stat1);
    d1->SetFlav(daughters[0]->Flav());
    d1->SetFromDec(true);
    // todo: 1->2 qcd shower with ew fs recoil partner
    // d1->SetK(idmother);// not that simple: w->qq' has color connection in fs
    Flavour prop_flav=blob->InParticle(0)->Flav().DecayHandler()->PropFlav(blob);
    Vec4D momd2=RecombinedMomentum(daughters[1],photons,stat2);
    Vec4D momd3=RecombinedMomentum(daughters[2],photons,stat3);
    Vec4D prop_mom=momd2+momd3;
    step1->CreateLeg(prop_mom, prop_flav);
    size_t idnew1=1<<(++imax);
    step1->Legs().back()->SetId(idnew1);
    step1->Legs().back()->SetStat(0);
    step1->Legs().back()->SetFromDec(true);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  step1->IdLeg(idmother),
                                  step1->Legs().back());
    DEBUG_VAR(*step1);
    Cluster_Amplitude* tmp=step1;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew1);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew1);
        }
      }
      DEBUG_VAR(*tmp);
    }

    
    Cluster_Amplitude* step2=step1->InitPrev();
    step2->CopyFrom(step1);
    step2->SetNLO(0);
    step2->SetFlag(1);
    step2->SetMS(step1->MS());
    for (size_t i=0; i<step1->Legs().size(); ++i)
      step1->Leg(i)->SetStat(step1->Leg(i)->Stat()|1);
    step1->IdLeg(idnew1)->SetStat(1|4);
    step1->IdLeg(idnew1)->SetK(idk);
    for (size_t i=0; i<step2->Legs().size(); ++i) step2->Leg(i)->SetK(0);
    Cluster_Leg *d2(step2->IdLeg(idnew1));
    d2->SetMom(momd2);
    d2->SetStat(stat2);
    d2->SetFlav(daughters[1]->Flav());
    step2->CreateLeg(momd3, daughters[2]->Flav());
    size_t idnew2=1<<(++imax);
    step2->Legs().back()->SetId(idnew2);
    step2->Legs().back()->SetStat(stat3);
    step2->Legs().back()->SetFromDec(true);
    Cluster_Amplitude::SetColours(step1->IdLeg(idnew1),
                                  step2->IdLeg(idnew1),
                                  step2->Legs().back());
    DEBUG_VAR(*step2);
    tmp=step2;
    while (tmp->Next()) {
      tmp=tmp->Next();
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idnew1) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew2);
	}
        if (tmp->Leg(i)->K()&idnew1) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew2);
        }
      }
      DEBUG_VAR(*tmp);
    }

    std::vector<size_t> ids;
    ids.push_back(idmother);
    ids.push_back(idnew1);
    ids.push_back(idnew2);
    while (photons.size())
      AddPhotonsClustering(step2,daughters,photons,imax,ids);
    if (daughters[0]->DecayBlob())
      AddDecayClustering(step2, daughters[0]->DecayBlob(), imax, idmother);
    if (daughters[1]->DecayBlob())
      AddDecayClustering(step2, daughters[1]->DecayBlob(), imax, idnew1);
    if (daughters[2]->DecayBlob())
      AddDecayClustering(step2, daughters[2]->DecayBlob(), imax, idnew2);
    ampl=step2;
  }
  else {
    PRINT_VAR(*blob);
    THROW(fatal_error, "1 -> n not implemented yet.");
  }
}

bool DecayDefineInitialConditions(Cluster_Amplitude* ampl,
                                  Blob* initial_blob,
                                  PDF::Cluster_Definitions_Base *p_clus)
{
  DEBUG_FUNC("");
  DEBUG_VAR(*ampl);
  for (int i=0; i<initial_blob->NOutP(); ++i) {
    ampl->Leg(initial_blob->NInP()+i)->SetMom
      (initial_blob->OutParticle(i)->Momentum());
  }
  if (p_clus->ReCluster(ampl)<0) return false;
  if (ampl->NIn()==2) {
    for (Cluster_Amplitude *campl(ampl);
	 campl;campl=campl->Next()) {
      if (-campl->Leg(0)->Mom()[0]>rpa->gen.PBeam(0)[0] ||
	  -campl->Leg(1)->Mom()[0]>rpa->gen.PBeam(1)[0])
	return false;
    }
  }
  size_t imax=ampl->Legs().size()-1;
  for (int i=0; i<initial_blob->NOutP(); ++i) {
    if (initial_blob->OutParticle(i)->DecayBlob()) {
      AddDecayClustering(ampl, initial_blob->OutParticle(i)->DecayBlob(),
                         imax, 1<<(initial_blob->NInP()+i));
    }
  }
  return true;
}





// needed for hadron decay + showers:
Cluster_Amplitude* DecayClusterConfiguration(Blob *const bl)
{
  bool m_cluster=false;
  msg_Indent();
  Cluster_Amplitude* p_ampl = Cluster_Amplitude::New();
  //p_ampl->SetMS(this);
  for (int i(0);i<bl->NInP();++i) {
    Particle *p(bl->InParticle(i));
    ColorID col(p->GetFlow(2),p->GetFlow(1));
    p_ampl->CreateLeg(-p->Momentum(),p->Flav().Bar(),col,1<<i);
  }
  p_ampl->SetNIn(bl->NInP());
  for (int i(0);i<bl->NOutP();++i) {
    Particle *p(bl->OutParticle(i));
    if (p->GetFlow(1)==0 && p->GetFlow(2)==0) continue;
    ColorID col(p->GetFlow(1),p->GetFlow(2));
    p_ampl->CreateLeg(p->Momentum(),p->Flav(),col,1<<(i+p_ampl->NIn()));
  }
  while (m_cluster && p_ampl->Legs().size()>p_ampl->NIn()+2) {
    msg_Debugging()<<*p_ampl<<"\n";
    Cluster_Amplitude *ampl(p_ampl);
    p_ampl = p_ampl->InitNext();
    //p_ampl->SetMS(this);
    for (size_t i(0);i<ampl->NIn();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }
    p_ampl->SetNIn(ampl->NIn());
    Cluster_Leg *lij(NULL);
    for (size_t i(ampl->NIn());i<ampl->Legs().size()-1;++i) {
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(i+1);j<ampl->Legs().size();++j) {
        Cluster_Leg *lj(ampl->Leg(j));
        ColorID nc;
        if (li->Col().m_i==0 && li->Col().m_j==0) {
          nc=lj->Col();
        }
        else if (lj->Col().m_i==0 && lj->Col().m_j==0) {
          nc=li->Col();
        }
        else if (li->Col().m_i && li->Col().m_i==lj->Col().m_j) {
          nc.m_i=lj->Col().m_i;
          nc.m_j=li->Col().m_j;
        }
        else if (li->Col().m_j && li->Col().m_j==lj->Col().m_i) {
          nc.m_i=li->Col().m_i;
          nc.m_j=lj->Col().m_j;
        }
        if (nc.m_i>=0 && nc.m_j>=0) {
          Flavour fl(kf_photon);
          if (nc.m_i && nc.m_j) fl=Flavour(kf_gluon);
          else if (nc.m_i) fl=Flavour(kf_d);
          else if (nc.m_j) fl=Flavour(kf_d).Bar();
          p_ampl->CreateLeg(li->Mom()+lj->Mom(),fl,nc,li->Id()+lj->Id());
          lij=p_ampl->Legs().back();
          break;
        }
      }
      if (lij) break;
    }
    if (lij==NULL) THROW(fatal_error,"Internal eror");
    for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      if (cl->Id()&lij->Id()) continue;
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }
  }
  double mu2=p_ampl->Leg(0)->Mom().Abs2();
  p_ampl->SetMuF2(mu2);
  p_ampl->SetKT2(mu2);
  p_ampl->SetMuQ2(mu2);
  msg_Debugging()<<*p_ampl<<"\n";
  while (p_ampl->Prev()) {
    p_ampl=p_ampl->Prev();
    p_ampl->SetMuF2(mu2);
    p_ampl->SetKT2(mu2);
    p_ampl->SetMuQ2(mu2);
  }
  msg_Debugging()<<"}\n";
  return p_ampl;
}
