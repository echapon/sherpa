#include "AMISIC++/Perturbative/Single_Collision_Handler.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Single_Collision_Handler::Single_Collision_Handler() :
  p_processes(NULL), p_overestimator(NULL),
  m_pt2(0.), m_pt2min(0.),  
  m_pbeam1(rpa->gen.PBeam(0)), m_pbeam2(rpa->gen.PBeam(1)),
  m_S((m_pbeam1+m_pbeam2).Abs2()), m_lastpt2(m_S), 
  m_residualx1(1.), m_residualx2(1.),
  m_xt(1.), m_ymax(0.), m_y3(0.), m_y4(0.), m_x1(1.), m_x2(1.),
  m_ana(true)
{}

Single_Collision_Handler::~Single_Collision_Handler() {
  if (m_ana) FinishAnalysis();
}

void Single_Collision_Handler::Init() {
  m_pt2min = sqr((*mipars)("pt_min"));
  if (m_ana) InitAnalysis();
}

Blob * Single_Collision_Handler::NextScatter(const double & bfac) {
  // Simple logic - bfactor is provided from outside, now
  // - produce a trial kinematics (new transverse momentum smaller than the last one,
  //   supplemented with rapidities, Mandelstams ...)
  // - select a process
  // - construct its full kinematics and set the colours
  // - return a filled (at the moment 2->2 only) scattering blob
  p_overestimator->SetBFac(bfac);
  do {
    if (!SelectPT2(m_lastpt2)) return NULL;
    p_proc = p_processes->SelectProcess();
  }
  while (!p_proc || !p_proc->MakeKinematics(m_pt2,m_y3,m_y4,sqrt(m_shat)) ||
	 !p_proc->SetColours());
  return MakeBlob();
}

bool Single_Collision_Handler::SelectPT2(const double & pt2) {
  // Generate a trial kinematics
  // - produces a trial pt2 based on a fast and crude overestimator (in TrialPT2).
  //   if it falls below the minimal value, false is returned which effectively stops
  //   the generation of secondary scatters
  // - produce rapidites for the two outgoing particles (flat distribution), reconstruct
  //   Bjorken-x for the PDFs and the Mandelstams
  // - calculate the cross section summed over all parton-level processes
  // - accept or reject the kinematics with a hit-or-miss of true over overestimated
  //   cross section
  m_pt2 = pt2;
  double sigmatrue, sigmaapprox, weight;
  bool success(false);
  do {
    m_pt2  = p_overestimator->TrialPT2(m_pt2);
    m_muf2 = m_mur2 = m_pt2;
    if (m_pt2<m_pt2min) return false;
    if (!SelectRapidities() || !CalcXs() || !CalcMandelstams()) continue;
    p_processes->CalcPDFs(m_x1,m_x2,m_pt2);
    sigmatrue   = (*p_processes)(m_shat,m_that,m_uhat) * m_yvol;
    sigmaapprox = (*p_overestimator)(m_pt2);
    weight      = sigmatrue/sigmaapprox;
    if (weight>1.) {
      msg_Out()<<"      * ratio = "<<(sigmatrue/sigmaapprox)<<" "
	       <<" from sigma(true) = "<<sigmatrue<<", "
	       <<"sigma(approx) = "<<sigmaapprox<<"\n";
    }
    if (m_ana) AnalyseWeight(weight);
    if (weight > ran->Get()) success = true;
  } while (!success);
  SetLastPT2(m_pt2);
  return true;
}

Blob * Single_Collision_Handler::MakeBlob() {
  Blob * blob = new Blob();
  blob->SetType(btp::Hard_Collision);
  blob->SetStatus(blob_status::needs_showers |
		  blob_status::needs_beams |
		  blob_status::needs_reconnections |
		  blob_status::needs_hadronization);
  blob->SetId();
  blob->AddData("Weights",new Blob_Data<Event_Weights>({0, 1.0}));
  for (size_t i=0;i<2;i++) blob->AddToInParticles(p_proc->GetParticle(i));
  for (size_t i=2;i<4;i++) blob->AddToOutParticles(p_proc->GetParticle(i));
  if (m_ana) Analyse(m_pt2,blob);
  return blob;
}


bool Single_Collision_Handler::SelectRapidities() {
  // Generate two trial rapidities
  m_xt   = sqrt(4.*m_pt2/m_S);
  if (m_xt>1.) return false;
  m_ymax = log(1./m_xt*(1.+sqrt(1.-m_xt*m_xt)));
  m_yvol = sqr(2.*m_ymax);
  m_y3   = m_ymax*(2.*ran->Get()-1.);
  m_y4   = m_ymax*(2.*ran->Get()-1.);
  return true;
}

bool Single_Collision_Handler::CalcXs() {
  m_x1   = m_xt*(exp(m_y3)+exp(m_y4))/2.;
  m_x2   = m_xt*(exp(-m_y3)+exp(-m_y4))/2.;
  return (m_x1<m_residualx1 && m_x2<m_residualx2);
}

bool Single_Collision_Handler::CalcMandelstams() {
  if (m_xt*m_xt>m_x1*m_x2) return false;
  double tanhy = sqrt(1.-(m_xt*m_xt)/(m_x1*m_x2));
  m_shat = m_x1*m_x2*m_S;
  m_that = -0.5*m_shat*(1.-tanhy);
  m_uhat = -0.5*m_shat*(1.+tanhy);
  return true;
}

void Single_Collision_Handler::InitAnalysis() {
  m_histos[string("weights")]     = new Histogram(0,0.,2.,200);
  m_histos[string("weights_low")] = new Histogram(0,0.,0.1,1000);
  m_histos[string("pt")]          = new Histogram(0,0.,200.,200);
  m_histos[string("flavs")]       = new Histogram(0,-6.5,6.5,13);
}

void Single_Collision_Handler::FinishAnalysis() {
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator 
	 hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    name  = string("MPI_Analysis/")+hit->first+string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}

void Single_Collision_Handler::AnalyseWeight(const double & weight) {
  m_histos[string("weights")]->Insert(weight);
  m_histos[string("weights_low")]->Insert(weight);
}
  
void Single_Collision_Handler::Analyse(const double & pt2,Blob * blob) {
  m_histos[string("pt")]->Insert(sqrt(pt2));
  Flavour flav1 = blob->OutParticle(0)->Flav();
  Flavour flav2 = blob->OutParticle(1)->Flav();
  int fl1 = size_t(flav1.Kfcode());
  if (flav1.IsAnti()) fl1 = -fl1;
  if (fl1==21) fl1=6;
  if (fl1==22) fl1=-6;
  int fl2 = size_t(flav2.Kfcode());
  if (flav2.IsAnti()) fl2 = -fl2;
  if (fl2==21) fl2=6;
  if (fl2==22) fl2=-6;
  m_histos[string("flavs")]->Insert(fl1);
  m_histos[string("flavs")]->Insert(fl2);
}
  
void Single_Collision_Handler::Test(const double & Q2,const long int & n) {
  msg_Out()<<METHOD<<" for Q^2 = "<<Q2<<", s = "<<m_S<<".\n";
  Histogram histo(0,0.0,Q2,100);
  for (long int dryrun=0;dryrun<n;dryrun++) {
    SetLastPT2(Q2);
    bool taken(false);
    while (NextScatter() && m_pt2>0) {
      if (!taken) {
	histo.Insert(m_pt2);
	taken = true;
      }
      SetLastPT2(m_pt2);
    }
  }
  histo.Finalize();
  histo.Output("True_PT2");
  msg_Out()<<METHOD<<": finished "<<n<<" dry runs.\n";
}
