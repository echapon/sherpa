#include "RECONNECTIONS/Main/Reconnection_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Handler::Reconnection_Handler(const bool & on) :
  m_on(on), m_weights(Reconnection_Weights(this)), m_analysis(true)
{ }

Reconnection_Handler::~Reconnection_Handler() {
  if (m_analysis) {
    for (map<string,Histogram *>::iterator hit=m_histomap.begin();
	 hit!=m_histomap.end();hit++) {
      Histogram * histo = hit->second;
      string name  = string("Reconnection_Analysis/")+hit->first+string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histomap.clear();
  }
}

void Reconnection_Handler::Initialize() {
  m_weights.Initialize(); 
  if (m_analysis) {
    m_histomap[string("Reconn_MassBefore")] = new Histogram(10,0.001,100.0,200);
    m_histomap[string("Reconn_MassAfter1")] = new Histogram(10,0.001,100.0,200);
    m_histomap[string("Reconn_MassAfter2")] = new Histogram(10,0.001,100.0,200);
  }
}


void Reconnection_Handler::Reset() {
  for (size_t pos=0;pos<2;pos++) { m_cols[pos].clear(); m_parts[pos].clear(); }
  while (!m_singlets.empty()) {
    delete m_singlets.front();
    m_singlets.pop_front();
  }  
  m_weights.Reset();
}

Return_Value::code Reconnection_Handler::operator()(Blob_List *const blobs,
						    Particle_List *const parts) {
  Return_Value::code ret = Return_Value::Nothing;
  if (m_on && HarvestParticles(blobs)) {
    MakeSinglets();
    if (m_analysis) FillMassesInHistogram(m_histomap["Reconn_MassBefore"]);
    m_weights.FillTables();
    ReshuffleSinglets();
    if (m_analysis) FillMassesInHistogram(m_histomap["Reconn_MassAfter1"]);
    ReconnectSinglets();
    if (m_analysis) FillMassesInHistogram(m_histomap["Reconn_MassAfter2"]);
    AddReconnectionBlob(blobs);
    ret = Return_Value::Success; 
  }
  Reset();
  return ret;
}
  
bool Reconnection_Handler::HarvestParticles(Blob_List * blobs) {
  // Extract all coloured particles from hadronization blob(s) and fill them into
  // relevant lists.  We assume that there already exists one (or more) active
  // fragmentation blob(s), with incoming coloured particles.
  // TODO: Will have to check how this works with two hadronization blobs, for
  //       example in e+e- -> W+W- -> 4 quarks or, more tricky, pp -> WW -> 4 quarks
  Blob * blob;
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    blob = (*bit);
    if (!blob->Has(blob_status::needs_reconnections)) continue;
    blob->SetTypeSpec("Colour Reconnections");
    for (int i=0;i<blob->NInP();i++) HarvestParticleInfo(blob->InParticle(i));
    blob->UnsetStatus(blob_status::needs_reconnections |
		      blob_status::needs_hadronization);
  }
  return (!m_cols[0].empty() && !m_cols[1].empty());
}

void Reconnection_Handler::HarvestParticleInfo(ATOOLS::Particle * part) {
  // Only add strong particles with colours different from zero.
  unsigned int col[2];
  for (size_t pos=0;pos<2;pos++) col[pos] = part->GetFlow(pos+1);
  if (col[0]==0 && col[1]==0) return;
  Particle * copy = new Particle(*part);
  colpair cols = colpair(col[0],col[1]);
  // We work with colour pairs <triplet, anti-triplet>.
  // Filling maps of triplet/anti-triplet colour indices to the particles;
  // constructing the singlets will use this, the two sets m_parts[] will be
  // used for the construction of the "distances" between particles.
  for (size_t pos=0;pos<2;pos++) {
    if (col[pos]!=0) {
      m_cols[pos][col[pos]] = copy;
      m_parts[pos].insert(copy);
    }
  }
  // Make sure we know the decay blob of the original particle - we will add the
  // copy as ingoing particle to the same blob, but with possibly reshuffled colours;
  // the originl particles will end in a reconnection blob, and their copies will
  // form its outgoing particles.
  copy->SetDecayBlob(part->DecayBlob());
  copy->SetProductionBlob(NULL);
}

void Reconnection_Handler::MakeSinglets() {
  // Start singlets with triplet particles only (quarks or anti-diquarks, with GetFlow(2)==0,
  // in FindStart), and look for fitting anti-triplet colours (in FindNext), until we hit an
  // anti-triplet particle and the singlet is finished.  Once triplet particles are exhausted,
  // repeat with gluons (again, in FindStart), and make sure we book the starting gluon only
  // once.
  Particle * part, * start = FindStart();
  size_t col = start->GetFlow(1);
  Part_List * sing = new Part_List;
  m_singlets.push_back(sing);
  sing->push_back(start);
  while (start) {
    part = FindNext(col);
    if (part!=start) sing->push_back(part);
    col  = part->GetFlow(1);
    if (col==0 || part==start) {
      start = FindStart();
      if (start) {
	col  = start->GetFlow(1);
	sing = new Part_List;
	m_singlets.push_back(sing);
	sing->push_back(start);
      }
    }
    else m_cols[0].erase(col);
  }
}

Particle * Reconnection_Handler::FindStart() {
  Particle * start(NULL);
  // Look for triplet-only particles.
  for (std::map<unsigned int, ATOOLS::Particle *>::iterator cpit = m_cols[0].begin();
       cpit!=m_cols[0].end();cpit++) {
    if (cpit->second->GetFlow(2)==0) {
      start = cpit->second;
      break;
    }
  }
  // Didn't find a triplet particle, look for any particle with triplet colour
  if (start==NULL && m_cols[0].size()>0) start = m_cols[0].begin()->second;
  if (start!=0) {
    // Make sure the "start" parton is taken out of potential singlet starters.
    m_cols[0].erase(start->GetFlow(1));
    return start;
  }
  // Nothing found, cannot start other singlets.
  return NULL;
}

Particle * Reconnection_Handler::FindNext(const size_t & col) {
  std::map<unsigned int, ATOOLS::Particle *>::iterator cpit = m_cols[1].find(col);
  if (cpit==m_cols[1].end()) {
    msg_Error()<<"Error in "<<METHOD<<" did not find particle with colour[1] = "<<col<<"\n";
    exit(1);
  }
  Particle * part = cpit->second;
  m_cols[1].erase(cpit);
  return part;
}

void Reconnection_Handler::ReshuffleSinglets() {
  // Go through the singlets.  Reshuffling makes sense only if you have more than 4
  // particles in the singlet.  Keep in mind: no splitting of singlets into two as
  // result of reshuffling.  
  for (list<Part_List *>::iterator sit=m_singlets.begin();sit!=m_singlets.end();sit++) {
    if ((*sit)->size()<4) continue;
    bool hit=true;
    while(hit) hit = ReshuffleSinglet((*sit));
  }
}

bool Reconnection_Handler::ReshuffleSinglet(Part_List * singlet) {
  Part_Iterator pit[6], stopit=singlet->end(); ((stopit--)--)--;
  // Logic: check if one particle should be reinserted at another position.
  // Assume particle 4 should be inserted between 0 and 1.  Then we need to compare
  // <01> * <34> * <45> with <04> * <41> * <35>
  // Similarly, if we insert 1 between 4 and 5, we need to compare
  // <01> * <12> * <45> with <02> * <41> * <15>
  pit[0] = singlet->begin(); 
  for (size_t i=1;i<3;i++) { pit[i] = pit[i-1]; pit[i]++; }  
  do {
    if (pit[2]==singlet->end()) return false;
    pit[3] = pit[1];
    for (size_t i=4;i<6;i++) {
      pit[i] = pit[i-1]; pit[i]++;
      if (pit[i]==singlet->end() || pit[i-1]==singlet->end()) break;
    }  
    while (pit[4]!=singlet->end() && pit[5]!=singlet->end()) {
      double dist01345 = (m_weights(*pit[0],*pit[1]) *
			  m_weights(*pit[3],*pit[4]) *
			  m_weights(*pit[4],*pit[5]));
      double dist04135 = (m_weights(*pit[0],*pit[4]) *
			  m_weights(*pit[4],*pit[1]) *
			  m_weights(*pit[3],*pit[5]));
      if (dist04135!=0 && ran->Get() * dist04135 < dist01345) {
	m_weights.SetWeight(*pit[0],*pit[1],0.);
	m_weights.SetWeight(*pit[3],*pit[4],0.);
	m_weights.SetWeight(*pit[4],*pit[5],0.);
	(*pit[4])->SetFlow(2,(*pit[0])->GetFlow(1));
	(*pit[1])->SetFlow(2,(*pit[4])->GetFlow(1));
	(*pit[5])->SetFlow(2,(*pit[3])->GetFlow(1));
	singlet->insert(pit[1],*pit[4]);
	singlet->erase(pit[4]);
	return true;
      }
      double dist01245 = (m_weights(*pit[0],*pit[1]) *
			  m_weights(*pit[1],*pit[2]) *
			  m_weights(*pit[4],*pit[5]));
      double dist02415 = (m_weights(*pit[0],*pit[2]) *
			  m_weights(*pit[4],*pit[1]) *
			  m_weights(*pit[1],*pit[5]));
      if (dist02415!=0 && ran->Get() * dist02415 < dist01245) {
	m_weights.SetWeight(*pit[0],*pit[1],0.);
	m_weights.SetWeight(*pit[1],*pit[2],0.);
	m_weights.SetWeight(*pit[4],*pit[5],0.);
	(*pit[2])->SetFlow(2,(*pit[0])->GetFlow(1));
	(*pit[1])->SetFlow(2,(*pit[4])->GetFlow(1));
	(*pit[5])->SetFlow(2,(*pit[1])->GetFlow(1));
	singlet->insert(pit[5],*pit[1]);
	singlet->erase(pit[1]);
	return true;
      }
      pit[3]++;
      for (size_t i=4;i<6;i++) {pit[i] = pit[i-1]; pit[i]++; }  
    }
    for (size_t i=0;i<3;i++) pit[i]++;
  } while (pit[2]!=singlet->end());
  return false;
}

void Reconnection_Handler::ReconnectSinglets() {
  list<Part_List *>::iterator sit1, sit2;
  Part_Iterator pit11, pit12, pit21, pit22;
  double dist11, dist12, dist21, dist22;
  bool hit;
  if (m_singlets.size()<2) return;
  do {
    hit = false;
    sit1=m_singlets.begin();
    sit2=m_singlets.begin();sit2++;
    while (sit2!=m_singlets.end() && !hit) {
      // Check for pairs of consecutive partons in both singlets
      // (pit11 & pit12 in singlet 1 and pit21 & pit22 in singlet 2)
      // if they should be reconnected cross-wise, i.e., if singlet 1 contains of
      // parton up to and including pit11 and continues with pit22 to the end of the
      // original singlet 2 and vice versa.  The decision is based on distances between
      // pit11/pit12 and pit21/pit22 vs. pit11/pit22 and pit21/pit12.
      pit11  = (*sit1)->begin(); pit12=pit11; pit12++;
      dist11 = m_weights((*pit11),(*pit12));
      while (pit12!=(*sit1)->end() && !hit) {
	pit21=(*sit2)->begin(); pit22=pit21; pit22++;
	while (pit22!=(*sit2)->end() && !hit) {
	  //msg_Out()<<METHOD<<" tests to shuffle particles: "
	  //	   <<(*pit11)->Number()<<" & "<<(*pit12)->Number()<<" vs "
	  //	   <<(*pit21)->Number()<<" & "<<(*pit22)->Number()<<"\n";
	  dist12 = m_weights((*pit11),(*pit22));
	  dist21 = m_weights((*pit21),(*pit12));
	  dist22 = m_weights((*pit21),(*pit22));
	  if (ran->Get() * (dist21*dist12) < (dist11*dist22)) {
	    //msg_Out()<<"   ("<<dist11<<" * "<<dist22<<" = "<<(dist11*dist22)<<") "
	    //	     <<"vs. ("<<dist21<<" * "<<dist12<<" = "<<(dist12*dist21)<<" * "<<rand<<").\n";
	    hit = true;
	    SpliceSinglets((*sit1),(*sit2),pit12,pit22);
	    AftermathOfSlicing((*pit11),(*pit12),(*pit21),(*pit22));
	    break;
	  }
	  if (!hit) { pit21++; pit22++; }
	}
	if (!hit) { pit11++; pit12++; }
      }
      sit2++;
    }
    sit1++;
  } while (hit);
}

void Reconnection_Handler::SpliceSinglets(Part_List * sing1,Part_List * sing2,
					  Part_Iterator & pit1,Part_Iterator & pit2) {
  Part_List help;
  help.splice(help.begin(),*sing1,pit1,sing1->end());
  sing1->splice(sing1->end(),*sing2,pit2,sing2->end());
  sing2->splice(sing2->end(),help);
  // msg_Out()<<"--------- Singlet 2 with "<<sing2->size()<<" particles.\n";
  // for (Part_Iterator pit=sing2->begin();pit!=sing2->end();pit++) {
  //   msg_Out()<<"  "<<(*pit)->Number()<<" ["<<(*pit)->GetFlow(1)<<", "<<(*pit)->GetFlow(2)<<"]\n";
  // }  
}

void Reconnection_Handler::AftermathOfSlicing(Particle * part11,Particle * part12,
					      Particle * part21,Particle * part22) {
  // After a slicing has happened, we set the weights to zero to ensure that the slicing is not
  // unmade in subsequent sweeps.
  m_weights.SetWeight(part11,part12,0.);
  m_weights.SetWeight(part11,part22,0.);
  m_weights.SetWeight(part21,part22,0.);
  m_weights.SetWeight(part21,part12,0.);
  part22->SetFlow(2,part11->GetFlow(1));
  part12->SetFlow(2,part21->GetFlow(1));
}

void Reconnection_Handler::AddReconnectionBlob(Blob_List *const blobs) {
  Blob * blob = new Blob();
  blob->AddStatus(blob_status::needs_hadronization);
  blob->SetType(btp::Fragmentation);
  blob->SetId();
  Particle * part;
  while (!m_singlets.empty()) {
    for (Part_Iterator pit=m_singlets.front()->begin();
	 pit!=m_singlets.front()->end();pit++) {
      part = (*pit);
      part->DecayBlob()->AddToOutParticles(part);
      part->SetDecayBlob(NULL);
      blob->AddToInParticles(part);
    }
    delete m_singlets.front(); m_singlets.pop_front();
  }
  blobs->push_back(blob);
}

void Reconnection_Handler::FillMassesInHistogram(Histogram * histo) {
  Part_Iterator pit1, pit2;
  for (list<Part_List *>::iterator sit=m_singlets.begin();sit!=m_singlets.end();sit++) {
    pit2 = pit1 = (*sit)->begin(); pit2++;
    while (pit2!=(*sit)->end()) {
      double mij = sqrt((((*pit1)->Flav().IsGluon()?0.5:1.)*(*pit1)->Momentum()+
			 ((*pit2)->Flav().IsGluon()?0.5:1.)*(*pit2)->Momentum()).Abs2());
      histo->Insert(mij,1.);
      pit2++; pit1++;
    }
  }
}
