#ifndef ATOOLS_Phys_Variations_H
#define ATOOLS_Phys_Variations_H

#include <ostream>
#include <map>
#include <string>
#include <vector>

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Enum_Flags.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#define ENABLE_REWEIGHTING_FACTORS_HISTOGRAMS 0

namespace PDF { class PDF_Base; }
namespace MODEL { class Running_AlphaS; }
namespace BEAM { class Beam_Spectra_Handler; }


namespace ATOOLS {

  class Variations;
  class Variation_Parameters;

  extern Variations* s_variations;

  struct ScaleFactorExpansions {
    enum code {
      None            = 0,
      MuF             = 1,
      MuR             = 1<<1,
      SevenPoint      = 1<<2,  // vary independently, but exclude (up,down) and (down,up) variations
      QCUT            = 1<<3
    };
  };
  DEFINE_ENUM_FLAG_OPERATORS(ScaleFactorExpansions::code)

  //! Initialise and hold all information that is needed for varying input parameters
  class Variations {
  public:
    static bool NeedsLHAPDF6Interface();
    static void CheckConsistencyWithBeamSpectra(BEAM::Beam_Spectra_Handler *);

    Variations();
    ~Variations();
    typedef std::vector<Variation_Parameters *> Parameters_Vector;
    const Parameters_Vector * GetParametersVector() const { return &m_parameters_vector; }
    // TODO: return a const object
    Variation_Parameters& Parameters(size_t i) { return *m_parameters_vector[i]; }
    std::string
        GetVariationNameAt(Variations::Parameters_Vector::size_type) const;
    size_t Size() { return m_parameters_vector.size(); }
    template <typename Handler>
    void ForEach(Handler h);

    // register and print warnings
    void IncrementOrInitialiseWarningCounter(const std::string name) { m_warnings[name]++; }
    void PrintStatistics(std::ostream &str);

  private:
    void ReadDefaults();
#if defined USING__LHAPDF && defined USING__LHAPDF6
    void LoadLHAPDFInterfaceIfNecessary();
#endif
    void InitialiseParametersVector();
    void AddParameters(Scoped_Settings&);
    struct PDFs_And_AlphaS {
      PDFs_And_AlphaS();
      PDFs_And_AlphaS(double alphasmz);
      PDFs_And_AlphaS(std::string pdfname, int pdfmember);
      std::vector<PDF::PDF_Base *> m_pdfs;
      MODEL::Running_AlphaS *p_alphas;
      bool m_shoulddeletepdf {false};
      bool m_shoulddeletealphas {false};
    };
    void AddParameterExpandingScaleFactors(std::vector<std::string> scalestringparams,
                                           ScaleFactorExpansions::code,
                                           std::vector<Variations::PDFs_And_AlphaS>);
    void AddParameters(double, double, double, std::vector<PDFs_And_AlphaS>::const_iterator, bool deletepdfs, bool deletealphas);
    std::vector<PDFs_And_AlphaS> PDFsAndAlphaSVector(std::string pdfstringparam, bool expandpdf);

    Parameters_Vector m_parameters_vector;

    //! (name -> counter value) map used to track warnings
    std::map<std::string, unsigned long> m_warnings;

    //! whether the central value should be included when expanding VARIATIONS arguments
    bool m_includecentralvaluevariation;
    //! whether the muR variation factor should be applied to the shower AlphaS scale arguments
    bool m_reweightsplittingalphasscales;
    //! whether the muF variation factor should be applied to the shower PDF scale arguments
    bool m_reweightsplittingpdfsscales;
  };

  template <typename Handler> void Variations::ForEach(Handler h)
  {
    const auto size = Size();
    for (int i {0}; i < size; ++i) {
      h(i, Parameters(i));
    }
  }

  std::ostream& operator<<(std::ostream&, const Variations&);

#if ENABLE_REWEIGHTING_FACTORS_HISTOGRAMS
  class ReweightingFactorHistogram {

  public:

    ReweightingFactorHistogram();

    void Fill(std::string, double);
    void Write(std::string filenameaffix);

  private:

    typedef std::string EntryKey;
    typedef std::vector<long> EntryNumbers;
    typedef std::pair<double, EntryNumbers> Bin;

    std::vector<EntryKey> m_keys;
    std::vector<Bin> m_bins;
  };
#endif


  /*!
   * Hold a set of input parameters and factors for a single input parameter variation
   *
   * Used as an argument to a reweighting function.
   */
  struct Variation_Parameters {
    Variation_Parameters(double muR2fac, double muF2fac,
                         double showermuR2fac, double showermuF2fac,
			 double Qcutfac,
                         PDF::PDF_Base *pdf1,
                         PDF::PDF_Base *pdf2,
                         MODEL::Running_AlphaS *alphas,
                         bool deletepdfs,
                         bool deletealphas):
      m_muR2fac(muR2fac), m_muF2fac(muF2fac),
      m_showermuR2fac(showermuR2fac), m_showermuF2fac(showermuF2fac), m_Qcutfac(Qcutfac),
      p_pdf1(pdf1), p_pdf2(pdf2), p_alphas(alphas),
      m_deletepdfs(deletepdfs),
      m_deletealphas(deletealphas),
      m_name(GenerateName())
    {};
    ~Variation_Parameters();

    void IncrementOrInitialiseWarningCounter(const std::string name) { m_warnings[name]++; }

#if ENABLE_REWEIGHTING_FACTORS_HISTOGRAMS
    void FillReweightingFactorsHisto(std::string name, double value) { m_rewfachisto.Fill(name, value); }
#endif

    bool IsTrivial() const;

    const double m_muR2fac, m_muF2fac;
    const double m_showermuR2fac, m_showermuF2fac, m_Qcutfac;
    //! Pointers to the beam PDFs, which can be NULL for (semi-)leptonic events
    PDF::PDF_Base * const p_pdf1;
    PDF::PDF_Base * const p_pdf2;
    MODEL::Running_AlphaS * const p_alphas;
    //! Set whether the pointers to the PDFs and the AlphaS should be deleted after usage
    const bool m_deletepdfs;
    const bool m_deletealphas;
    const std::string m_name;

  private:
    //! Return a name based on the scale and PDF variations
    std::string GenerateName() const;
    //! Return part of a name based on a tag name and the associated value
    template <typename U>
    std::string GenerateNamePart(std::string tag, U value) const;
    //! (name -> counter value) map used to track warnings
    std::map<std::string, unsigned long> m_warnings;
    friend class Variations;
#if ENABLE_REWEIGHTING_FACTORS_HISTOGRAMS
    ReweightingFactorHistogram m_rewfachisto;
#endif
  };


  //! Convenience wrapper for storing and passing subevent weights corresponding to one variation
  class Subevent_Weights_Vector : public std::vector<double> {
  public:
    // inherit useful vector constructors
    //! Construct an empty weight vector
    Subevent_Weights_Vector();
    //! Construct a weight vector of size times 1.0 (i.e. we do not use the standard default 0.0!)
    Subevent_Weights_Vector(size_type count, const double& value = 1.0);
    //! Multiply all weights with a common factor
    Subevent_Weights_Vector & operator*=(const double&);
    /*!
     * Multiply all subevent weights one by one with those of another vec
     *
     * The other subevent weights must either be of equal size or they have to
     * contain only one entry. In the latter case, this single entry is used as
     * a common factor for all subevent weights.
     */
    Subevent_Weights_Vector & operator*=(const Subevent_Weights_Vector&);
    Subevent_Weights_Vector & operator+=(const Subevent_Weights_Vector &);
  };

  std::ostream& operator<<(std::ostream&, const Subevent_Weights_Vector&);


  //! Types to denote different variation sources
  enum class Variations_Type {
    all,     //!< use for weight retrieval, to get the product of the below
    main,    //!< use for everything that is not one of the below
    sudakov  //!< use for weight factors from the reweighting of the Sudakov
             //!< Veto Algorithm, i.e. the reweighting of parton shower
             //!< emissions
  };

  std::ostream& operator<<(std::ostream&, const Variations_Type&);

};

#endif // #ifndef ATOOLS_Phys_Variations_H
