#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"

#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Single_Events/Decay_Handler_Base.H"
#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Variations.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Random.H"

#include <cassert>

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;

Perturbative_Interface::Perturbative_Interface(Matrix_Element_Handler *const meh,
                                               Hard_Decay_Handler*const dec,
                                               Shower_Handler *const psh):
  p_me(meh), p_dec(dec), p_mi(NULL), p_hd(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL)
{
  Settings& s = Settings::GetMainSettings();
  m_bbarmode = s["METS_BBAR_MODE"].SetDefault(1).Get<int>();
  m_globalkfac = s["GLOBAL_KFAC"].SetDefault(0.0).Get<double>();
  m_maxkfac = s["MENLOPS_MAX_KFAC"].SetDefault(10.0).Get<double>();
}

Perturbative_Interface::Perturbative_Interface
(MI_Handler *const mi,Shower_Handler *const psh):
  p_me(NULL), p_mi(mi), p_hd(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Decay_Handler_Base *const hdh,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(hdh), p_sc(NULL), p_shower(psh),
  p_ampl(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Soft_Collision_Handler *const sch,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(NULL), p_sc(sch), p_shower(psh),
  p_ampl(NULL) {}

Perturbative_Interface::~Perturbative_Interface() 
{
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
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
  if (!p_shower->On() ||
      (p_me && p_me->Process()->Info().m_nlomode==nlo_mode::fixedorder)) {
    m_weights = Event_Weights {};
    return Return_Value::Success;
  }
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
    p_ampl=NULL;
  }
  p_shower->CleanUp();
  msg_Indent();
  if (p_mi) {
    p_ampl=p_mi->ClusterConfiguration(blob);
    if (p_ampl==NULL) return Return_Value::Retry_Event;
    if (p_ampl->Leg(0)->Mom()[3]*p_ampl->Leg(1)->Mom()[3]>0.0) {
      msg_Tracking()<<METHOD<<"(): Invalid beams. Retry event."<<std::endl;
      return Return_Value::Retry_Event;
    }
    p_mi->SetMassMode(1);
    int stat(p_mi->ShiftMasses(p_ampl));
    if (stat<0) {
      msg_Error()<<METHOD<<"(): MI Mass shift failed. Reject event.\n"
		 <<(*blob)<<"\n";
      exit(1);
      return Return_Value::Retry_Event;
    }
    if (stat==1) {
      stat=p_mi->Shower()->GetShower()->GetClusterDefinitions()->ReCluster(p_ampl);
      if (stat!=1) {
	msg_Error()<<METHOD<<"(): MI Reclustering failed. Reject event.\n"
		   <<(*blob)<<"\n";
	exit(1);
	return Return_Value::Retry_Event;
      }
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl))
      return Return_Value::New_Event;
    return Return_Value::Success;
  }
  if (p_hd) {
    p_ampl=p_hd->ClusterConfiguration(blob);
    if (!p_shower->GetShower()->PrepareShower(p_ampl))
      return Return_Value::New_Event;
    return Return_Value::Success;
  }
  if (p_sc) {
    p_sc->SetClusterDefinitions(p_shower->GetShower()->GetClusterDefinitions());
    p_ampl=p_sc->ClusterConfiguration(blob);
    if (p_ampl==NULL) {
      msg_Error()<<METHOD<<": Soft_Collision_Handler has no amplitude.\n";
      return Return_Value::New_Event;
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl,true)) {
      msg_Error()<<METHOD<<": could not prepare shower.\n"; 
      return Return_Value::New_Event;
    }
    return Return_Value::Success;
  }
  assert(p_me != NULL);
  p_ampl=p_me->Process()->Get<Single_Process>()->Cluster
    (p_me->Process()->Integrator()->Momenta());
  if (p_ampl==NULL) {
    msg_Error()<<METHOD<<"(): Clustering failed. Reject event."<<std::endl;
    return Return_Value::New_Event;
  }
  if (p_ampl->Leg(0)->Mom()[3]*p_ampl->Leg(1)->Mom()[3]>0.0) {
    msg_Tracking()<<METHOD<<"(): Invalid beams. Reject event."<<std::endl;
    return Return_Value::New_Event;
  }
  m_weights = Event_Weights {};
  if (p_me->Process()->Info().m_ckkw&1) {
    if ((m_bbarmode&1) && p_me->HasNLO() &&
        p_me->Process()->Parent()->Info().m_fi.NLOType()==nlo_type::lo) {
      if (!LocalKFactor(p_ampl)) DEBUG_INFO("Process not found");
    }
  }
  p_me->Process()->Generator()->SetMassMode(1);
  size_t cmax(0);
  for (size_t i(0);i<p_ampl->Legs().size();++i)
    cmax=Max(cmax,(size_t)p_ampl->Leg(i)->Col().m_i);
  while (Flow::Counter()<cmax);
  p_me->Process()->Parent()->SetRBMap(p_ampl);
  if (p_dec) {
    p_dec->SetCluster(p_me->Shower()->GetShower()->GetClusterDefinitions());
    if (!p_dec->DefineInitialConditions(p_ampl, blob)) {
      msg_Tracking()<<METHOD<<"(): Decay clustering failed. Reject event."<<std::endl;
      return Return_Value::Retry_Event;
    }
  }
  while (p_ampl->Prev()) p_ampl=p_ampl->Prev();
  if (p_me->Process()->Info().m_ckkw&1) {
    blob->AddData("Sud_Weights",new Blob_Data<Event_Weights>(m_weights));
    if (p_me->EventGenerationMode()!=0) {
      const auto disc = ran->Get();
      const auto abswgt = std::abs(m_weights.Nominal());
      if (abswgt < disc) {
        return Return_Value::New_Event;
      }
      m_weights /= Min(1.0, abswgt);
    }
    Blob_Data_Base* winfo((*blob)["Weights"]);
    if (!winfo) THROW(fatal_error,"No weights information in signal blob");
    Event_Weights meweights(winfo->Get<Event_Weights>());
    blob->AddData("Weights",new Blob_Data<Event_Weights>(meweights*m_weights));
  }
  if (!p_shower->GetShower()->PrepareShower(p_ampl)) 
    return Return_Value::New_Event;
  return Return_Value::Success;
}

bool Perturbative_Interface::LocalKFactor(ATOOLS::Cluster_Amplitude* ampl)
{
  if (m_globalkfac) {
    m_weights *= m_globalkfac;
    return true;
  }
  DEBUG_FUNC(ampl->Legs().size());
  Process_Vector procs(p_me->AllProcesses());
  Process_Base::SortFlavours(ampl);
  while (ampl->Next()!=NULL) {
    ampl=ampl->Next();
    Process_Base::SortFlavours(ampl);
    if (m_bbarmode&2) {
      if (ampl->Next()) continue;
      Single_Process *proc(ampl->Proc<Single_Process>());
      if (ampl->Legs().size()-ampl->NIn()>
	  proc->Info().m_fi.NMinExternal()) break;
    }
    for (int i=procs.size()-1; i>=0; --i) {
      MCatNLO_Process* mcnloproc=dynamic_cast<MCatNLO_Process*>(procs[i]);
      if (mcnloproc) {
	Event_Weights K {mcnloproc->LocalKFactor(*ampl)};
	if (K.Nominal()==0.0 || dabs(K.Nominal())>m_maxkfac) continue;
        m_weights *= K;
        return true;
      }
    }
  }
  // no process found along ampl
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
    if (!p_hd) {
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
  PDF::Shower_Base *csh(p_shower->GetShower());
  int stat=csh->PerformShowers();
  Event_Weights weights = csh->Weights();
  m_weights *= weights;

  // NOTE: we leave shower-related weights out of the main weights, in order to
  // be able to output both ME-only and full variations
  p_hard->AddData("Shower_Weight",new Blob_Data<double>(weights.Nominal()));
  p_hard->AddData("Shower_Weights",new Blob_Data<Event_Weights>(weights));

  return stat;
}

int Perturbative_Interface::PerformDecayShowers()
{ 
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

