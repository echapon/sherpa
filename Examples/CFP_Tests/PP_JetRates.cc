// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {

  class PP_JetRates : public Analysis {
  private:

    double m_r;

    std::vector<Histo1DPtr> _h_log10_d;

  public:

    PP_JetRates(): Analysis("PP_JetRates") {}

    void init()
    {
      m_r=0.4;
      addProjection(FinalState(), "FS");
      for (size_t i=0;i<4;++i) {
	string dname="log10_d_"+to_str(i)+to_str(i+1);
	_h_log10_d.push_back(bookHisto1D(dname,100,0.0,log10(0.5*sqrtS()/GeV)));
      }
    }

    void analyze(const Event& event)
    {
      const double weight=event.weight();
      const FinalState &fs = applyProjection<FinalState>(event, "FS");
      std::vector<FourMomentum> jets;
      for (size_t i(0);i<fs.size();++i) {
	long int kfc(abs(fs.particles()[i].pdgId()));
	if ((kfc>=11 && kfc<=16) || (kfc>=22 && kfc<=25)) continue;
	MSG_DEBUG(i<<": "<<fs.particles()[i].pdgId()<<" "<<fs.particles()[i].momentum());
	jets.push_back(fs.particles()[i].momentum());
      }
      std::vector<double> dij2 = Cluster(jets);
      for (size_t i=0;i<min(_h_log10_d.size(),dij2.size());++i)
	if (dij2[dij2.size()-1-i])
	  _h_log10_d[i]->fill(0.5*log10(dij2[dij2.size()-1-i]), weight);
    }

    void finalize()
    {
      const double xsec_unitw=crossSection()/picobarn/sumOfWeights();
      for (size_t i=0;i<_h_log10_d.size();++i) scale(_h_log10_d[i],xsec_unitw);
    }

    double R2(const FourMomentum &p, const FourMomentum &q) const
    {
      double dy(p.rapidity()-q.rapidity());
      double dphi(p.phi()-q.phi());
      return (cosh(dy)-cos(dphi))/sqr(m_r);
    }

    double Q2i(const FourMomentum &p) const
    {
      return p.pT2();
    }

    double Q2ij(const FourMomentum &p,const FourMomentum &q) const
    {
      return min(p.pT2(),q.pT2())*R2(p,q);
    }

    std::vector<double> Cluster(std::vector<FourMomentum> &p)
    {
      std::vector<double> kt2;
      if (p.size()==1) kt2.push_back(Q2i(p[0]));
      if (p.size()<=1) return kt2;
      int ii=0, jj=0, n=p.size();
      std::vector<int> imap(p.size());
      for (size_t i(0);i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > kt2ij(n,std::vector<double>(n));
      double dmin=std::numeric_limits<double>::max();
      for (int i=0;i<n;++i) {
	double di=kt2ij[i][i]=Q2i(p[i]);
	if (di<dmin) { dmin=di; ii=i; jj=i; }
	for (int j=0;j<i;++j) {
	  double dij=kt2ij[i][j]=Q2ij(p[i],p[j]);
	  if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	}
      }
      while (n>0) {
	MSG_DEBUG("Q_{"<<n<<"->"<<n-1<<"} = "<<sqrt(dmin)
		  <<" <- "<<p[imap[jj]]<<" "<<p[imap[ii]]);
	if (ii!=jj) p[imap[jj]]+=p[imap[ii]];
	kt2.push_back(dmin); --n;
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	int jjx=imap[jj];
	kt2ij[jjx][jjx]=Q2i(p[jjx]);
	for (int j=0;j<jj;++j) kt2ij[jjx][imap[j]]=Q2ij(p[jjx],p[imap[j]]);
	for (int i=jj+1;i<n;++i) kt2ij[imap[i]][jjx]=Q2ij(p[jjx],p[imap[i]]);
	ii=jj=0; dmin=kt2ij[imap[0]][imap[0]];
	for (int i=0;i<n;++i) {
	  int ix=imap[i]; double di=kt2ij[ix][ix];
	  if (di<dmin) { dmin=di; ii=jj=i; }
	  for (int j=0;j<i;++j) {
	    int jx=imap[j]; double dij=kt2ij[ix][jx];
	    if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	  }
	}
      }
      return kt2;
    }

  };

  DECLARE_RIVET_PLUGIN(PP_JetRates);

}
