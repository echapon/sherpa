#include "AddOns/OpenLoops/GGH_Process_Manager.H"

#include "MODEL/Main/Model_Base.H"

#include "COMIX/Main/Process_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Exception.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"

using namespace PHASIC;
using namespace ATOOLS;

GGH_Process_Manager::GGH_Process_Manager() : p_generators(NULL)  {}

Process_Base* GGH_Process_Manager::InitializeProcess(const ATOOLS::Cluster_Amplitude& ampl, 
						     bool external, const std::vector<double>& orders){
  DEBUG_FUNC(this);
  //ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process",true);
  //My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");

  // build a process info instance
  Process_Info pi;
  pi.m_megenerator="Amegic";
  for (size_t i(0);i<ampl.NIn();++i) {
    Flavour fl(ampl.Leg(i)->Flav().Bar());
    //if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_ii.m_ps.push_back(Subprocess_Info(fl,"",""));
  }
  for (size_t i(ampl.NIn());i<ampl.Legs().size();++i) {
    Flavour fl(ampl.Leg(i)->Flav());
    //if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_fi.m_ps.push_back(Subprocess_Info(fl,"",""));
  }

  // set coupling orders correctly
  pi.m_maxcpl = pi.m_mincpl = orders;
  if(external){
    // set weirdly abused mhv-flag to get external (i.e. OpenLoops) proc
    pi.m_amegicmhv = 10;
    pi.m_loopgenerator = "OpenLoops";
    // order counting in OpenLoops is different
    for(size_t i(2); i<orders.size(); i++){
      pi.m_maxcpl[1] += pi.m_maxcpl[i];
      pi.m_mincpl[1] += pi.m_mincpl[i];
      pi.m_maxcpl[i]  = pi.m_mincpl[i] = 0.0;
    }
  }
  DEBUG_VAR(pi);
  // initialize the process
  PHASIC::Process_Base *proc= Generators()->InitializeProcess(pi,false);
  if (!proc) {
    //My_In_File::CloseDB(ATOOLS::rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
    PRINT_VAR(pi);
    THROW(fatal_error, "Could not initialize auxiliary process");
  }

  // set selector, kfactor, and scale setter
  proc->SetSelector(Selector_Key{});
  proc->SetScale(Scale_Setter_Arguments(MODEL::s_model,"VAR{sqr("+ATOOLS::ToString(rpa->gen.Ecms())+")}","Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("None"));
  
  //proc->Get<COMIX::Process_Base>()->Tests();
  //My_In_File::Close(ATOOLS::rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
  m_maps.push_back(new NLOTypeStringProcessMap_Map);
  m_procs.push_back(proc);
  proc->FillProcessMap(m_maps.back());
  return proc;
}

Process_Base* GGH_Process_Manager::GetProcess(const std::string& name, bool external){
    for(Process_Vector::const_iterator it=m_procs.begin(); it!=m_procs.end(); ++it)
    {
      // m_amegicmhv==10 signals a Single_Process_External i.e. OpenLoops proc
      if( ((*it)->Info().m_amegicmhv==10) != external )
	continue;
      NLOTypeStringProcessMap_Map::const_iterator jt = (*it)->AllProcs()->find(nlo_type::lo);
      if (jt == (*it)->AllProcs()->end()) 
	continue;
      StringProcess_Map::const_iterator kt = jt->second->find(name);
      if (kt != jt->second->end())
	return kt->second;
    }
    return NULL;
}
  
Process_Base* GGH_Process_Manager::GetProcess(const ATOOLS::Cluster_Amplitude& ampl, bool external,
					      const std::vector<double>& orders){
  std::string name = Process_Base::GenerateName(&ampl);
  Process_Base* ret = GetProcess(name, external);
  if(!ret){
    InitializeProcess(ampl, external, orders);
    ret = GetProcess(name, external);
  }
  if(!ret)
    THROW(fatal_error, "Failed to initialize process "+name);
  return ret;
}

ME_Generators* GGH_Process_Manager::Generators(){
  if(!p_generators) THROW(fatal_error, "Generators not set");
  return p_generators;
}

GGH_Process_Manager::~GGH_Process_Manager(){
  for(Map_Vector::const_iterator it=m_maps.begin(); it!=m_maps.end(); ++it){
    for(NLOTypeStringProcessMap_Map::const_iterator jt=(*it)->begin(); jt!=(*it)->end(); ++jt){
      for(StringProcess_Map::const_iterator kt=jt->second->begin(); kt!=jt->second->end(); ++kt)
	if(kt->second) delete kt->second;
      delete jt->second;
    }
    delete *it;
  }
}
