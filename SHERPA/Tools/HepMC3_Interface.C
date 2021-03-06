#include "SHERPA/Tools/HepMC3_Interface.H"
#ifdef USING__HEPMC3

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"


#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/Units.h"


using namespace SHERPA;
using namespace ATOOLS;
namespace HepMC = HepMC3;

EventInfo3::EventInfo3(ATOOLS::Blob * sp, const double &wgt,
                     bool namedweights,
                     bool extendedweights,
                     bool includemeonlyweights) :
  p_sp(sp),
  m_wgt(wgt),
  m_usenamedweights(namedweights),
  m_extendedweights(extendedweights),
  m_variationtypes(1, ATOOLS::Variations_Type::all),
  m_mewgt(0.), m_wgtnorm(wgt), m_ntrials(1.),
  m_pswgt(0.), m_pwgt(0.),  m_userhook(false), m_userweight(0.),
  m_mur2(0.), m_muf12(0.), m_muf22(0.),
  m_alphas(0.), m_alpha(0.), m_type(ATOOLS::nlo_type::lo),
  p_wgtinfo(NULL), p_pdfinfo(NULL), p_subevtlist(NULL)
{
  if (p_sp) {
    DEBUG_FUNC(*p_sp);
    Blob_Data_Base *db;
    ReadIn(db,"MEWeight",false);
    if (db) m_mewgt=db->Get<double>();
    m_pswgt=m_wgt/m_mewgt;
    ReadIn(db,"Weight_Norm",true);
    m_wgtnorm=db->Get<double>();
    ReadIn(db,"Trials",true);
    m_ntrials=db->Get<double>();
    ReadIn(db,"PDFInfo",false);
    if (db) {
      p_pdfinfo=&db->Get<ATOOLS::PDF_Info>();
      m_muf12=p_pdfinfo->m_muf12;
      m_muf22=p_pdfinfo->m_muf22;
    }
    ReadIn(db,"UserHook",false);
    if (db) {
      m_userhook=true;
      m_userweight=db->Get<double>();
    }
    ReadIn(db,"Renormalization_Scale",false);
    if (db) m_mur2=db->Get<double>();
    SetAlphaS();
    SetAlpha();
    if (m_extendedweights) {
      ReadIn(db,"Orders",true);
      m_orders=db->Get<std::vector<double> >();
      ReadIn(db,"MEWeightInfo",true);
      p_wgtinfo=db->Get<ME_Weight_Info*>();
    }
    ReadIn(db,"NLO_subeventlist",false);
    if (db) p_subevtlist=db->Get<NLO_subevtlist*>();
    if (p_subevtlist) m_type=p_subevtlist->Type();

    if (includemeonlyweights)
      m_variationtypes.push_back(ATOOLS::Variations_Type::main);

    ReadIn(db, "Weights", false);
    if (db) {
      m_weights = db->Get<Event_Weights>();
      if (m_weights.ContainsVariations() && !m_usenamedweights)
        THROW(fatal_error,
              "Scale and/or PDF variations cannot be written to " +
                  std::string("HepMC without using named weights. ") +
                  std::string("Try HEPMC_USE_NAMED_WEIGHTS: true"));
    }

    ReadIn(db, "Shower_Weights", false);
    m_showerweights = Event_Weights {};
    if (db) {
      m_showerweights *= db->Get<Event_Weights>();
    }
    ReadIn(db, "MC@NLO_Shower_Weights", false);
    if (db) {
      m_showerweights *= db->Get<Event_Weights>();
    }

  }
}

EventInfo3::EventInfo3(const EventInfo3 &evtinfo) :
  m_usenamedweights(evtinfo.m_usenamedweights),
  m_extendedweights(evtinfo.m_extendedweights),
  m_variationtypes(evtinfo.m_variationtypes),
  p_sp(evtinfo.p_sp),
  m_orders(evtinfo.m_orders),
  m_wgt(0.), m_mewgt(0.), m_wgtnorm(0.),
  m_ntrials(evtinfo.m_ntrials), m_pswgt(evtinfo.m_pswgt), m_pwgt(0.),
  m_mur2(0.), m_muf12(0.), m_muf22(0.),m_muq2(0.),
  m_alphas(0.), m_alpha(0.), m_type(evtinfo.m_type),
  p_wgtinfo(NULL), p_pdfinfo(evtinfo.p_pdfinfo),
  p_subevtlist(evtinfo.p_subevtlist),
  m_weights(evtinfo.m_weights),
  m_showerweights(evtinfo.m_showerweights)
{
}

void EventInfo3::ReadIn(ATOOLS::Blob_Data_Base* &db,std::string name,bool abort)
{
  db=(*p_sp)[name];
  if (abort && !db) THROW(fatal_error,name+" information missing.");
}

bool EventInfo3::WriteTo(HepMC::GenEvent &evt, const int& idx)
{
  DEBUG_FUNC("use named weights: "<<m_usenamedweights
             <<", extended weights: "<<m_extendedweights);
  std::map<std::string,double> wc;
  if (m_usenamedweights) {

 // fill standard entries to ensure backwards compatability
    wc["Weight"]=m_wgt;
    wc["MEWeight"]=m_mewgt;
    wc["WeightNormalisation"]=m_wgtnorm;
    wc["NTrials"]=m_ntrials;
    if (m_userhook) wc["UserHook"]=m_userweight;
    if (m_extendedweights) {
      wc["PSWeight"]=m_pswgt;
      // additional entries for LO/LOPS reweighting
      // x1,x2,muf2 can be found in PdfInfo; alphaS,alphaQED in their infos
      wc["MuR2"]=m_mur2;
      //Default, but why not in HepMC2 wc["MuQ2"]=m_muq2;
      wc["OQCD"]=m_orders[0];
      wc["OEW"]=m_orders[1];
      if (p_wgtinfo) {
        wc["Reweight_B"]=p_wgtinfo->m_B;
        wc["Reweight_MuR2"]=p_wgtinfo->m_mur2;
        wc["Reweight_MuF2"]=p_wgtinfo->m_muf2;
        if (p_wgtinfo->m_type&mewgttype::VI) {
          wc["Reweight_VI"]=p_wgtinfo->m_VI;
          for (size_t i=0;i<p_wgtinfo->m_wren.size();++i) {
            wc["Reweight_VI_wren_"+ToString(i)]=p_wgtinfo->m_wren[i];
          }
        }
        if (p_wgtinfo->m_type&mewgttype::KP) {
          wc["Reweight_KP"]=p_wgtinfo->m_KP;
          wc["Reweight_KP_x1p"]=p_wgtinfo->m_y1;
          wc["Reweight_KP_x2p"]=p_wgtinfo->m_y2;
          for (size_t i=0;i<p_wgtinfo->m_wfac.size();++i) {
            wc["Reweight_KP_wfac_"+ToString(i)]=p_wgtinfo->m_wfac[i];
          }
        }
        if (p_wgtinfo->m_type&mewgttype::DADS &&
            p_wgtinfo->m_dadsinfos.size()) {
          wc["Reweight_DADS_N"]=p_wgtinfo->m_dadsinfos.size();
          for (size_t i(0);i<p_wgtinfo->m_dadsinfos.size();++i) {
            wc["Reweight_DADS_"+ToString(i)+"_Weight"]
                =p_wgtinfo->m_dadsinfos[i].m_wgt;
            if (p_wgtinfo->m_dadsinfos[i].m_wgt) {
              wc["Reweight_DADS_"+ToString(i)+"_x1"]
                  =p_wgtinfo->m_dadsinfos[i].m_x1;
              wc["Reweight_DADS_"+ToString(i)+"_x2"]
                  =p_wgtinfo->m_dadsinfos[i].m_x2;
              wc["Reweight_DADS_"+ToString(i)+"_fl1"]
                  =p_wgtinfo->m_dadsinfos[i].m_fl1;
              wc["Reweight_DADS_"+ToString(i)+"_fl2"]
                  =p_wgtinfo->m_dadsinfos[i].m_fl2;
            }
          }
        }
        if (p_wgtinfo->m_type&mewgttype::METS &&
            p_wgtinfo->m_clusseqinfo.m_txfl.size()) {
          wc["Reweight_ClusterStep_N"]=p_wgtinfo->m_clusseqinfo.m_txfl.size();
          for (size_t i(0);i<p_wgtinfo->m_clusseqinfo.m_txfl.size();++i) {
            wc["Reweight_ClusterStep_"+ToString(i)+"_t"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_t;
            wc["Reweight_ClusterStep_"+ToString(i)+"_x1"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_xa;
            wc["Reweight_ClusterStep_"+ToString(i)+"_x2"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_xb;
            wc["Reweight_ClusterStep_"+ToString(i)+"_fl1"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_fla;
            wc["Reweight_ClusterStep_"+ToString(i)+"_fl2"]
                =p_wgtinfo->m_clusseqinfo.m_txfl[i].m_flb;
          }
        }
        if (p_wgtinfo->m_type&mewgttype::H) {
          wc["Reweight_RDA_N"]=p_wgtinfo->m_rdainfos.size();
          for (size_t i(0);i<p_wgtinfo->m_rdainfos.size();++i) {
            wc["Reweight_RDA_"+ToString(i)+"_Weight"]
                =p_wgtinfo->m_rdainfos[i].m_wgt;
            if (p_wgtinfo->m_rdainfos[i].m_wgt) {
              const double mur2(p_wgtinfo->m_rdainfos[i].m_mur2);
              wc["Reweight_RDA_"+ToString(i)+"_MuR2"] = mur2;
              wc["Reweight_RDA_"+ToString(i)+"_MuF12"]
                  =p_wgtinfo->m_rdainfos[i].m_muf12;
              wc["Reweight_RDA_"+ToString(i)+"_MuF22"]
                  =p_wgtinfo->m_rdainfos[i].m_muf22;
              wc["Reweight_RDA_"+ToString(i)+"_Dipole"]
                  =10000*p_wgtinfo->m_rdainfos[i].m_i
                    +100*p_wgtinfo->m_rdainfos[i].m_j
                        +p_wgtinfo->m_rdainfos[i].m_k;
              wc["Reweight_RDA_"+ToString(i)+"_AlphaS"]
                  =MODEL::s_model->ScalarFunction("alpha_S", mur2);
            }
          }
        }
        wc["Reweight_Type"]=p_wgtinfo->m_type;
      }
      if (p_subevtlist) {
        wc["Reweight_RS"]=m_pwgt;
        wc["Reweight_Type"]=64|(p_wgtinfo?p_wgtinfo->m_type:0);
      }
    }
    else {
      // if using minimal weights still dump event type if RS need correls
      if (p_subevtlist) wc["Reweight_Type"]=64;
    }
    // fill weight variations into weight container
    if (m_weights.ContainsVariations()) {
      size_t numvars = s_variations->Size();
      msg_Debugging() << "#named wgts: " << numvars << std::endl;
      for (size_t i(0); i < numvars; ++i) {
        std::string varname(s_variations->Parameters(i).m_name);
        typedef std::vector<ATOOLS::Variations_Type>::const_iterator It_type;
        for (It_type it(m_variationtypes.begin()); it != m_variationtypes.end();
             ++it) {
          const std::string typevarname((*it == ATOOLS::Variations_Type::main)
                                            ? "ME_ONLY_" + varname
                                            : varname);
          double weight {0.0};
          if (idx == -1) {
            weight = m_weights.Variation(i);
          } else {
            weight = (*p_subevtlist)[idx]->m_results.Variation(i);
          }
          if (*it == ATOOLS::Variations_Type::all) {
            weight *= m_showerweights.Variation(i);
          }
          wc[typevarname] = weight;
        }
      }
    }

  std::vector<std::string> w_names;
  std::vector<double> w_values;
  for (std::map<std::string,double>::iterator it=wc.begin(); it!=wc.end();++it)
  { w_names.push_back(it->first); w_values.push_back(it->second);}
  evt.run_info()->set_weight_names(w_names);
  evt.weights()=w_values;
  
  }
  else {
    // only offer basic event record for unnamed weights
      evt.weights().clear();
      evt.weights().push_back(m_wgt);
      evt.weights().push_back(m_mewgt);
      evt.weights().push_back(m_wgtnorm);
      evt.weights().push_back(m_ntrials);
    if (m_extendedweights) {
        evt.weights().push_back(m_pswgt);
        evt.weights().push_back(m_mur2);
        evt.weights().push_back(m_muf12);
        evt.weights().push_back(m_muf22);
        evt.weights().push_back(m_orders[0]);
        evt.weights().push_back(m_orders[1]);
    }
     evt.weights().push_back(p_subevtlist?64:0);
  }
  //evt.weights()=wc;
  if (p_pdfinfo) {
    double q(sqrt(sqrt(p_pdfinfo->m_muf12*p_pdfinfo->m_muf22)));
    HepMC::GenPdfInfoPtr pdfinfo = std::make_shared<HepMC::GenPdfInfo>();
    pdfinfo->set(p_pdfinfo->m_fl1,p_pdfinfo->m_fl2,
                           p_pdfinfo->m_x1,p_pdfinfo->m_x2,
                           q,p_pdfinfo->m_xf1,p_pdfinfo->m_xf2);
    evt.set_pdf_info(pdfinfo);
  }
   std::shared_ptr<HepMC::Attribute> a_event_scale = std::make_shared<HepMC::DoubleAttribute>(m_mur2);
   std::shared_ptr<HepMC::Attribute> a_alphas = std::make_shared<HepMC::DoubleAttribute>(m_alphas);
   std::shared_ptr<HepMC::Attribute> a_alpha = std::make_shared<HepMC::DoubleAttribute>(m_alpha);
   evt.add_attribute("alphaQCD",a_alphas);
   evt.add_attribute("alphaQED",a_alpha);
   evt.add_attribute("event_scale",a_event_scale);
  return true;
}

void EventInfo3::SetAlphaS()
{
  m_alphas=MODEL::s_model->ScalarFunction("alpha_S",m_mur2);
}

void EventInfo3::SetAlpha()
{
  m_alpha=MODEL::s_model->ScalarConstant("alpha_QED");
}

HepMC3_Interface::HepMC3_Interface() :
  m_usenamedweights(false),
  m_extendedweights(false),
  m_includemeonlyweights(false),
  m_hepmctree(false),
  p_event(NULL)
{
  Settings& s = Settings::GetMainSettings();
  m_usenamedweights =
    s["HEPMC_USE_NAMED_WEIGHTS"].SetDefault(false).Get<bool>();
  m_extendedweights =
    s["HEPMC_EXTENDED_WEIGHTS"].SetDefault(false).Get<bool>();
  m_includemeonlyweights =
    s["HEPMC_INCLUDE_ME_ONLY_VARIATIONS"].SetDefault(false).Get<bool>();
  // Switch for disconnection of 1,2,3 vertices from PS vertices
  m_hepmctree = s["HEPMC_TREE_LIKE"].SetDefault(false).Get<bool>();
}

HepMC3_Interface::~HepMC3_Interface()
{
  if (p_event) p_event->clear();
  delete p_event;
  DeleteGenSubEventList();
}


bool HepMC3_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         HepMC::GenEvent& event)
{
  const auto weight(blobs->Weight());
  event.set_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  if (!sp) sp=blobs->FindFirst(btp::Hard_Collision);
  Blob *mp(blobs->FindFirst(btp::Hard_Collision));  
  if (!mp) event.add_attribute("mpi", std::make_shared<HepMC::IntAttribute>(-1));
  
  EventInfo3 evtinfo(sp,weight,
                    m_usenamedweights,m_extendedweights,m_includemeonlyweights);
  // when subevtlist, fill hepmc-subevtlist
  if (evtinfo.SubEvtList()) return SubEvtList2ShortHepMC(evtinfo);
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  evtinfo.WriteTo(event);
  HepMC::GenVertexPtr vertex=std::make_shared<HepMC::GenVertex>();
  std::vector<HepMC::GenParticlePtr> beamparticles, inparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();
       blit!=blobs->end();++blit) {
    Blob* blob=*blit;
    if (m_ignoreblobs.count(blob->Type())) continue;
    for (int i=0;i<blob->NInP();i++) {
      if (blob->InParticle(i)->ProductionBlob()==NULL) {
        Particle* parton=blob->InParticle(i);
        auto flav = parton->Flav();
        if (flav.Kfcode() == kf_lepton) {
          if (sp->InParticle(0)->Flav().IsLepton()) {
            flav = sp->InParticle(0)->Flav();
          } else {
            flav = sp->InParticle(1)->Flav();
          }
        }
        HepMC::GenParticlePtr inpart =
            MakeGenParticle(parton->Momentum(), flav, true);
        event.add_particle(inpart);
        vertex->add_particle_in(inpart);
        //We add attributes here->
         for (int k=1;k<3;k++) {
         if (blob->InParticle(i)->GetFlow(k)>0)inpart->add_attribute("flow"+std::to_string((long long int)k),std::make_shared<HepMC::IntAttribute>(blob->InParticle(i)->GetFlow(k)));
         }
         //<-We add attributes here
         inparticles.push_back(inpart);
        // distinct because SHRIMPS has no bunches for some reason
        if (blob->Type()==btp::Beam || blob->Type()==btp::Bunch) {
          beamparticles.push_back(inpart);
        }
      }
    }
    for (int i=0;i<blob->NOutP();i++) {
      if (blob->OutParticle(i)->DecayBlob()==NULL ||
	  m_ignoreblobs.count(blob->OutParticle(i)->DecayBlob()->Type())!=0) {
        Particle* parton=blob->OutParticle(i);
        HepMC::GenParticlePtr outpart =
            MakeGenParticle(parton->Momentum(), parton->Flav(), false);
        event.add_particle(outpart);
        //We add attributes here->
        for (int k=1;k<3;k++) {
        if (blob->OutParticle(i)->GetFlow(k)>0)outpart->add_attribute("flow"+std::to_string((long long int)k),std::make_shared<HepMC::IntAttribute>(blob->OutParticle(i)->GetFlow(k)));
        }
        //<-We add attributes here 
        vertex->add_particle_out(outpart);
      }
    }
  
      if (mp==(*blit)) event.add_attribute("mpi", std::make_shared<HepMC::IntAttribute>(vertex->id()));
      if (sp==(*blit)) event.add_attribute("signal_process_vertex", std::make_shared<HepMC::IntAttribute>(vertex->id()));
  
  }
  event.add_vertex(vertex);
  vertex->add_attribute("weight0",std::make_shared<HepMC::DoubleAttribute>(1.0));
  if (beamparticles.empty() && inparticles.size()==2) {
    for (size_t j(0);j<2;++j) {
      HepMC::GenVertexPtr  beamvertex = std::make_shared<HepMC::GenVertex>();
      event.add_vertex(beamvertex);
      auto flav = (j == 0) ? rpa->gen.Beam1() : rpa->gen.Beam2();
      if (flav.Kfcode() == kf_lepton) {
        flav = sp->InParticle(j)->Flav();
      }
      HepMC::GenParticlePtr beampart =
          MakeGenParticle(rpa->gen.PBeam(j), flav, true);
      beamparticles.push_back(beampart);
      beamvertex->add_particle_in(beampart);
      beamvertex->add_particle_out(inparticles[j]);
    }
  }
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
  }

  return true;
}

bool HepMC3_Interface::SubEvtList2ShortHepMC(EventInfo3 &evtinfo)
{
  DEBUG_FUNC("subevts: "<<evtinfo.SubEvtList()->size());
  // build GenEvent for all subevts (where only the signal is available)
  // purely partonic, no beam information, may add, if needed
  for (size_t i(0);i<evtinfo.SubEvtList()->size();++i) {
    EventInfo3 subevtinfo(evtinfo);
    const NLO_subevt * sub((*evtinfo.SubEvtList())[i]);
    if (sub->m_result==0. &&
	!(sub->IsReal() && m_subeventlist.empty())) continue;
    HepMC::GenVertexPtr subvertex=std::make_shared<HepMC::GenVertex>();
    HepMC::GenEvent * subevent(new HepMC::GenEvent());
    // set the event number (could be used to identify correlated events)
    subevent->set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
    // assume that only 2->(n-2) processes, flip for Comix, flavs are correct
    HepMC::GenParticlePtr beamparticles[2]={NULL,NULL};
    for (size_t j(0);j<2;++j) {
      HepMC::GenVertexPtr  beamvertex = std::make_shared<HepMC::GenVertex>();
      subevent->add_vertex(beamvertex);
      auto flav = (j == 0) ? rpa->gen.Beam1() : rpa->gen.Beam2();
      if (flav.Kfcode() == kf_lepton) {
        flav = sub->p_fl[j];
      }
      beamparticles[j] = MakeGenParticle(rpa->gen.PBeam(j), flav, true);
      beamvertex->add_particle_in(beamparticles[j]);
      double flip(sub->p_mom[i][0]<0.);
      Vec4D momentum {flip ? -1.0 * sub->p_mom[j] : sub->p_mom[j]};
      HepMC::GenParticlePtr inpart =
          MakeGenParticle(momentum, sub->p_fl[j], true);
      subvertex->add_particle_in(inpart);
      beamvertex->add_particle_out(inpart);
//We add attributes here->
//FIXME!     for (int k=1;k<3;k++) {
//FIXME!    if (inpart->GetFlow(k)>0)outpart->add_attribute("flow"+std::to_string((long long int)k),std::make_shared<HepMC::IntAttribute>(inpart->GetFlow(k)));
//FIXME!     }
//<-We add attributes here
    }
    subevent->set_beam_particles(beamparticles[0],beamparticles[1]);
    for (size_t j(2);j<sub->m_n;++j) {
      HepMC::GenParticlePtr outpart =
          MakeGenParticle(sub->p_mom[j], sub->p_fl[j], false);
      subvertex->add_particle_out(outpart);
//We add attributes here->
//FIXME!     for (int k=1;k<3;k++) {
//FIXME!    if (outpart->GetFlow(k)>0)outpart->add_attribute("flow"+std::to_string((long long int)k),std::make_shared<HepMC::IntAttribute>(outpart->GetFlow(k)));
//FIXME!     }
//<-We add attributes here
    }
    subevent->add_vertex(subvertex);
    subvertex->add_attribute("weight0",std::make_shared<HepMC::DoubleAttribute>(1.0));
    // not enough info in subevents to set PDFInfo properly,
    // so set flavours and x1, x2 from the Signal_Process
    // reset muR, muF, alphaS, alpha
    subevtinfo.SetWeight(sub->m_result);
    subevtinfo.SetPartonicWeight(sub->m_mewgt);
    subevtinfo.SetMuR2(sub->m_mu2[stp::ren]);
    subevtinfo.SetMuF12(sub->m_mu2[stp::fac]);
    subevtinfo.SetMuF22(sub->m_mu2[stp::fac]);
    subevtinfo.SetAlphaS();
    subevtinfo.SetAlpha();
    subevtinfo.WriteTo(*subevent,i);
    m_subeventlist.push_back(subevent);
  }
  return true;
}


bool HepMC3_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs)
{
  if (blobs->empty()) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
               <<"   Empty list - nothing to translate into HepMC."<<std::endl
               <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event!=NULL) delete p_event;
  DeleteGenSubEventList();
  p_event = new HepMC::GenEvent();
  return Sherpa2ShortHepMC(blobs, *p_event);
}

std::shared_ptr<HepMC::GenParticle> HepMC3_Interface::MakeGenParticle(
    const Vec4D& mom, const Flavour& flav, bool incoming)
{
  HepMC::FourVector momentum(mom[1], mom[2], mom[3], mom[0]);
  int status = 1;
  if (incoming)
    status = (flav.StrongCharge() == 0) ? 4 : 11;
  return std::make_shared<HepMC::GenParticle>(momentum, (long int)flav, status);
}

// HS: Short-hand that takes a blob list, creates a new GenEvent and
// calls the actual Sherpa2HepMC
bool HepMC3_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
				   std::shared_ptr<HepMC::GenRunInfo> run)
{
  if (blobs->empty()) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
               <<"   Empty list - nothing to translate into HepMC."<<std::endl
               <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event) { p_event->clear(); delete p_event;}
  for (size_t i=0; i<m_subeventlist.size();++i)
    { m_subeventlist[i]->clear(); delete m_subeventlist[i];}
  m_subeventlist.clear();
  p_event = new HepMC::GenEvent(run);
  return Sherpa2HepMC(blobs, *p_event);
}

// The actual code --- calls the Blob to GenVertex code
bool HepMC3_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
                                    HepMC::GenEvent& event)
{
  const auto weight(blobs->Weight());
  DEBUG_FUNC("");
  event.set_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
  if (!m_hepmctree) event.add_attribute("cycles",std::make_shared<HepMC::IntAttribute>(1));
  // Signal Process blob --- there is only one
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  if (!sp) sp=blobs->FindFirst(btp::Hard_Collision);
  Blob *mp(blobs->FindFirst(btp::Hard_Collision));
  if (!mp) event.add_attribute("mpi", std::make_shared<HepMC::IntAttribute>(-1));
  
  
  // Meta info
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  EventInfo3 evtinfo(sp,weight,
                    m_usenamedweights,m_extendedweights,m_includemeonlyweights);
  evtinfo.WriteTo(event);
  
  m_blob2genvertex.clear();
  m_particle2genparticle.clear();
  HepMC::GenVertexPtr  vertex,psvertex;
  std::vector<HepMC::GenParticlePtr> beamparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();
       blit!=blobs->end();++blit) {
    // Call the Blob to vertex code, changes vertex pointer above
    if (Sherpa2HepMCBlobtoGenVertex(*(blit),vertex,event)) {
      event.add_vertex(vertex);
      for (size_t i(0);i<(*blit)->NInP();++i)
	if ((*blit)->InParticle(i)->ProductionBlob()==NULL) {
	  psvertex=vertex;
	  break;
	}

      if (mp==(*blit)) event.add_attribute("mpi", std::make_shared<HepMC::IntAttribute>(vertex->id()));
      if (sp==(*blit)) event.add_attribute("signal_process_vertex", std::make_shared<HepMC::IntAttribute>(vertex->id()));
      
      if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
        if ((**blit)["NLO_subeventlist"]) {
          THROW(fatal_error,"Events containing correlated subtraction events"
                +std::string(" cannot be translated into the full HepMC event")
                +std::string(" format.\n")
                +std::string("   Try 'EVENT_OUTPUT: HepMC_Short' instead."));
        }
        //event.set_signal_process_vertex(vertex);
      }
      // Find beam particles
      else if ((*blit)->Type()==ATOOLS::btp::Beam || 
	       (*blit)->Type()==ATOOLS::btp::Bunch) {
        for (auto pit: vertex->particles_in())
             {
          if (pit->production_vertex()==NULL) {
            beamparticles.push_back(pit);
          }
        }
      }
    }
  } // End Blob_List loop
  if (beamparticles.empty()) {
    HepMC::GenParticlePtr inpart[2];
    for (size_t j {0}; j < 2; ++j) {
      auto flav = (j == 0) ? rpa->gen.Beam1() : rpa->gen.Beam2();
      inpart[j] = MakeGenParticle(rpa->gen.PBeam(j), flav, true);
      psvertex->add_particle_in(inpart[j]);
    }
    event.set_beam_particles(inpart[0], inpart[1]);
  }


  // Disconnect ME, MPI and hard decay vertices from PS vertices to get a
  // tree-like record --- manipulates the final GenEvent
  // Can't use ->set_production_vertex/set_end_vertex as they are private
  // Need to use GenVertex::remove_particle(Pointer to particle)
  // But: iterator loses validity when calling GenVertex::remove_particle in the particle loop
  // Hence: fill vector with pointers and call GenVertex::remove_particle 
  if (m_hepmctree) {
    DEBUG_INFO("HEPMC_TREE_LIKE true --- straighten to "
               <<"tree enabled (disconnect 1,2,3 vertices)");
    // Iterate over all vertices to find PS vertices
    int vtx_id = -1;
    for (auto vit:  event.vertices()) {

      // Is this a PS Vertex?
      if (vit->status()==4) {
        std::vector<HepMC::GenParticlePtr > remove;
        //// Loop over outgoing particles
        for (auto  pout: vit->particles_out()) 
{
            if ( pout->end_vertex() ) {
              vtx_id = pout->end_vertex()->status(); //
              // Disconnect outgoing particle from end/decay vertex of type (1,2,3)
              if (vtx_id==1 || vtx_id==2 || vtx_id==3 )
                  remove.push_back(pout);
            }
        }
        // Loop over incoming particles
        for (auto  pin: vit->particles_in())
 {
          vtx_id = pin->production_vertex()->status();
          // Disconnect incoming particle from production vertex of type (1,2,3)  
          if (vtx_id==1 || vtx_id==2 || vtx_id==3 )
                  remove.push_back(pin);
        }
        // Iterate over Genparticle pointers to remove from current vertex (*vit)
        for (unsigned int nrem=0;nrem<remove.size();++nrem) {
            vit->remove_particle_in(remove[nrem]);
            vit->remove_particle_out(remove[nrem]);
        }
      } // Close if statement (vertex status==4)
    } // Close loop over vertices
  }

    if (beamparticles.size()==2) {
    event.add_tree( beamparticles );
  }
  
  return true;
}

// HS: this converts a Blob to GenVertex
bool HepMC3_Interface::Sherpa2HepMCBlobtoGenVertex(ATOOLS::Blob * blob, HepMC::GenVertexPtr & vertex, HepMC::GenEvent& event)
{
  if (m_ignoreblobs.count(blob->Type())) return false;
  int count = m_blob2genvertex.count(blob);
  if (count>0) {
    vertex = m_blob2genvertex[blob];
    return true;
  }
  else {
    ATOOLS::Vec4D pos = blob->Position();
    HepMC::FourVector position(pos[1],pos[2],pos[3],pos[0]);
    vertex = std::make_shared<HepMC::GenVertex>(position);
    event.add_vertex(vertex);
    vertex->add_attribute("weight0",std::make_shared<HepMC::DoubleAttribute>(1.0));
    if (blob->Type()==btp::Signal_Process)      vertex->set_status(1); // signal
    else if (blob->Type()==btp::Hard_Collision)vertex->set_status(2); // mpi
    else if (blob->Type()==btp::Hard_Decay)    vertex->set_status(3); // hard-decay
    else if (blob->Type()==btp::Shower || 
             blob->Type()==btp::QED_Radiation)  vertex->set_status(4); // PS/QED
    else if (blob->Type()==btp::Fragmentation)  vertex->set_status(5); // frag
    else if (blob->Type()==btp::Hadron_Decay)   vertex->set_status(6); // had-decay
      //{  
      //if ((*blob)["Partonic"]!=NULL) vertex->set_status(-6);
      //else vertex->set_status(6);
      //}
    else vertex->set_status(0);
  }

  bool okay = 1;
  HepMC::GenParticlePtr  _particle;
  for (int i=0;i<blob->NInP();i++) {
    if (Sherpa2HepMC(blob->InParticle(i),_particle)) {
      vertex->add_particle_in(_particle);
//We add attributes here->
    for (int k=1;k<3;k++) {
    if (blob->InParticle(i)->GetFlow(k)>0)_particle->add_attribute("flow"+std::to_string((long long int)k),std::make_shared<HepMC::IntAttribute>(blob->InParticle(i)->GetFlow(k)));
     }
//<-We add attributes here     

    }
    else okay = 0;
  }
  for (int i=0;i<blob->NOutP();i++) {
    if (Sherpa2HepMC(blob->OutParticle(i),_particle)) {
      vertex->add_particle_out(_particle);
//We add attributes here->
    for (int k=1;k<3;k++) {
    if (blob->OutParticle(i)->GetFlow(k)>0)_particle->add_attribute("flow"+std::to_string((long long int)k),std::make_shared<HepMC::IntAttribute>(blob->OutParticle(i)->GetFlow(k)));
     }
//<-We add attributes here 
    }
    else okay = 0;
  }
  m_blob2genvertex.insert(std::make_pair(blob,vertex));
  if (!okay) {
    msg_Error()<<"Error in HepMC3_Interface::Sherpa2HepMC(Blob,Vertex).\n"
               <<"    Continue event generation with new event."<<std::endl;
  }
  if (msg_LevelIsDebugging()) {
    ATOOLS::Vec4D check = blob->CheckMomentumConservation();
    double test         = ATOOLS::Vec3D(check).Abs();
    /*
    if (ATOOLS::dabs(1.-vertex->check_momentum_conservation()/test)>1.e-5 &&
        ATOOLS::dabs(test)>1.e-5) {
      msg_Error()<<"ERROR in "<<METHOD<<std::endl
                 <<"    Momentum not conserved. Continue."<<std::endl
                 <<"ERROR in Blob -> Vertex : "
                 <<vertex->check_momentum_conservation()
                 <<" <- "<<test<<" "<<check
                 <<std::endl<<(*blob)<<std::endl;
      vertex->print(msg_Error());
      msg_Error()<<"-----------------------------------------------"<<std::endl;
    }*/
  }
  return okay;
}

// HS: Sherpa Particle to HepMC::Genparticle --- fills m_particle2genparticle
// and changes the pointer reference particle ('new')
bool HepMC3_Interface::Sherpa2HepMC(ATOOLS::Particle * parton,
                                    HepMC::GenParticlePtr & particle)
{
  // HS: do nothing if parton has already been converted  
  int count = m_particle2genparticle.count(parton);
  if (count>0) {
    particle = m_particle2genparticle[parton];
    return true;
  }

  // Translate momentum vector
  ATOOLS::Vec4D mom  = parton->Momentum();
  HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);

  int status=11;
  // Assign status 1 to stable blobs (those without decays)
  // or for that Rivet specific bit, set the particle stable (1)
  // if its DecayBlob has been cut out
  if (parton->DecayBlob()==NULL ||
      m_ignoreblobs.count(parton->DecayBlob()->Type())!=0) {
    status=1;
  }
  // Non-stable particles --- what about Hard_Decay?
  else {
    if (parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Decay ||
        parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Mixing) {
      status=2;
    }
    // Set all particles going in/out of ME to status 3
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Signal_Process ||
             (parton->ProductionBlob() &&
              parton->ProductionBlob()->Type()==ATOOLS::btp::Signal_Process)) {
      status=3;
    }
    // E - gamma collider specific
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Bunch) {
      status=4;
    }
  }
  if (parton->Status()==part_status::documentation) status=20;
  particle = std::make_shared<HepMC::GenParticle>(momentum,(long int)parton->Flav(),status);
  m_particle2genparticle.insert(std::make_pair(parton,particle));
  return true;
}

void HepMC3_Interface::AddCrossSection(HepMC::GenEvent& event,
                                       const double &xs, const double &err)
{
    std::shared_ptr<HepMC::GenCrossSection> cross_section =std::make_shared<HepMC::GenCrossSection>();
    cross_section->set_cross_section(xs,err);
    event.set_cross_section(cross_section);

}

void HepMC3_Interface::DeleteGenSubEventList()
{
  for (size_t i=0; i<m_subeventlist.size();++i)
    delete m_subeventlist[i];
  m_subeventlist.clear();
}

#endif
