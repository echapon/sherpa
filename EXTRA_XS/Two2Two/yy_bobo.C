#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class yy_bobo : public ME2_Base {
  private:
    double m_mass;
    int m_charge;

    void doSelfTest();
  public:

    yy_bobo(const External_ME_Args& args);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };

  yy_bobo::yy_bobo(const External_ME_Args& args)
    : ME2_Base(args)
  {
    m_sintt=1;
    m_oew=0;
    m_oqcd=0;
    const Flavour_Vector fl = args.Flavours();
    m_mass = fl[2].HadMass();
    m_charge = fl[2].Charge();
    if (m_charge == 0) THROW(fatal_error, "You should change a known particle's mass to zero!");

    // doSelfTest(); // TODO trigger calculation based on message level
  }

  double yy_bobo::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // This is a dummy matrix element to simulate the production of charged
    // boson pairs in EPA

    double p2((momenta[0]+momenta[1]).Abs2());
    double m2 = m_mass * m_mass;

    if ((4*m2) > p2) return 0.;
    // check whether two particle production threshold is reached

    // bit more elegant calculation
    double xs = 0;
    xs-=4*m2/p2*(4-8*m2/p2);
    xs*=log(sqrt(p2)/(2*m_mass) + sqrt(p2/(4*m2) - 1));
    xs+=sqrt(1-4*m2/p2)*(2+8*m2/p2);
    xs *= 2*M_PI*0.0072992701*0.0072992701*pow(m_charge, 4.)/p2;
    //std::cout << "xs=" << xs << "\t" << "p2=" << p2 <<std::endl;
    return xs;
    /*
    // Calculation according to the paper
    double xs_par, xs_perp;
    xs_par =sqrt(1 - 4*m2/p2)*(1+6*m2/p2);
    xs_par-=4*m2/p2*(2 - 6*m2/p2)*log(sqrt(p2)/(2*m_mass) + sqrt(p2/(4*m2 )- 1));
    xs_par*=2*M_PI*0.0072992701*0.0072992701*m_charge*m_charge*m_charge*m_charge/p2;
    xs_perp =sqrt(1 - 4*m2/p2)*(1+2*m2/p2);
    xs_perp-=4*m2/p2*(2 - 2*m2/p2)*log(sqrt(p2)/(2*m_mass) + sqrt(p2/(4*m2 )- 1));
    xs_perp*=2*M_PI*0.0072992701*0.0072992701*m_charge*m_charge*m_charge*m_charge/p2;
    return xs_par + xs_perp;
    */
  }
}

DECLARE_TREEME2_GETTER(yy_bobo,"yy_bobo")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,yy_bobo>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_photon) && fl[1]==Flavour(kf_photon) &&
      (fl[2].Kfcode()==kf_pi_plus || fl[2].Kfcode()==kf_K_plus ||
       fl[2].Kfcode()==kf_D_plus || fl[2].Kfcode()==kf_D_s_plus) &&
      fl[3]==fl[2].Bar())
  {
    return new yy_bobo(args);
  }
  return NULL;
}

void yy_bobo::doSelfTest() {
  // TODO: base selftest on operator()
  std::cout << std::endl;
  std::cout << "yy_bobo - self test" << std::endl;
  std::cout << "mass=" << m_mass << std::endl;
  std::cout << "charge=" << m_charge << std::endl;
  std::cout << "p2\t xs\n\n";
  double s = 4*m_mass*m_mass;
  while(s<100) {
    double p2=s;
    s*=1.1;
    double m2 = m_mass * m_mass;
    // check whether two particle production threshold is reached
    /*
    double xs = 0;
    xs-=4*m2/p2*(4-8*m2/p2);
    xs*=log(sqrt(p2)/(2*m_mass) + sqrt(p2/(4*m2) - 1));
    xs+=sqrt(1-4*m2/p2)*(2+8*m2/p2);
    xs *= 2*M_PI*0.0072992701*0.0072992701*pow(m_charge, 4.)/p2;
    std::cout << p2 << "\t" << xs <<std::endl;
    */
    double xs_par, xs_perp;
    xs_par =sqrt(1 - 4*m2/p2)*(1+6*m2/p2);
    xs_par-=4*m2/p2*(2 - 6*m2/p2)*log(sqrt(p2)/(2*m_mass) + sqrt(p2/(4*m2 )- 1));
    xs_par*=2*M_PI*0.0072992701*0.0072992701*m_charge*m_charge*m_charge*m_charge/p2;
    xs_perp =sqrt(1 - 4*m2/p2)*(1+2*m2/p2);
    xs_perp-=4*m2/p2*(2 - 2*m2/p2)*log(sqrt(p2)/(2*m_mass) + sqrt(p2/(4*m2 )- 1));
    xs_perp*=2*M_PI*0.0072992701*0.0072992701*m_charge*m_charge*m_charge*m_charge/p2;

    std::cout << s << "\t" << xs_par + xs_perp << std::endl;
  }
  std::cout << "------------------------------------------------------------------!!\n";

}
