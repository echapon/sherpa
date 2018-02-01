#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#else
#include <fstream>
#endif
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHERPA;
using namespace ATOOLS;

namespace SHERPA {

  class Output_Weights: public Output_Base {
  private:

    std::string m_basename, m_ext;

#ifdef USING__GZIP
    ATOOLS::ogzstream m_outstream;
#else
    std::ofstream     m_outstream;
#endif

  public:

    Output_Weights(const Output_Arguments &args):
      Output_Base("Weights")
    {
      MyStrStream basename;
      basename << args.m_outpath << "/" << args.m_outfile;
      m_ext=".wts";
#ifdef USING__GZIP
      m_ext += ".gz";
#endif
#ifdef USING__MPI
      if (MPI::COMM_WORLD.Get_size()>1) {
        basename << "_" << MPI::COMM_WORLD.Get_rank();
      }
#endif
      m_basename = basename.str();
      m_outstream.open((m_basename+m_ext).c_str());
      if (!m_outstream.good())
	THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
      int precision = args.p_reader->GetValue<int>("OUTPUT_PRECISION",12);
      m_outstream.precision(precision);
    }

    ~Output_Weights()
    {
      m_outstream.close();
    }

    void Output(Blob_List* blobs)
    {
      Blob_Data_Base *sd((*blobs->FindFirst(btp::Signal_Process))["NLO_subeventlist"]);
      NLO_subevtlist *subs=sd?sd->Get<NLO_subevtlist*>():NULL;
      Blob_Data_Base *data((*blobs->FindFirst(btp::Signal_Process))["Variation_Weights"]);
      if (data==NULL) THROW(fatal_error,"Variation weights not found.");
      ATOOLS::Variation_Weights *variationweights=&data->Get<Variation_Weights>();
      if (variationweights==NULL) THROW(fatal_error,"Variation weights not found.");
      size_t numvars = variationweights->GetNumberOfVariations();
      if (subs==NULL) {
	m_outstream<<rpa->gen.NumberOfGeneratedEvents()<<" ";
	for (size_t i(0);i<numvars;++i)
	  m_outstream<<variationweights->GetVariationNameAt(i)<<" "
		     <<variationweights->GetVariationWeightAt(i)<<" ";
      }
      else
	for (size_t j(0);j<subs->size();++j) {
	  m_outstream<<rpa->gen.NumberOfGeneratedEvents()<<" ";
	  for (size_t i(0);i<numvars;++i)
	    m_outstream
              <<variationweights->GetVariationNameAt(i)
              <<" "
              <<variationweights->GetVariationWeightAt(i,Variations_Type::all,j)
              <<" ";
	  m_outstream<<"\n";
	}
      m_outstream<<"\n";
    }

    void ChangeFile()
    {
      m_outstream.close();
      std::string newname(m_basename+m_ext);
      for (size_t i(0);FileExists(newname);
	   newname=m_basename+"."+ToString(++i)+m_ext);
      m_outstream.open(newname.c_str());
      if (!m_outstream.good())
	THROW(fatal_error, "Could not open event file "+newname+".");
    }

  };

}

DECLARE_GETTER(Output_Weights,"Weights",Output_Base,Output_Arguments);
Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_Weights>::
operator()(const Output_Arguments &args) const
{ return new Output_Weights(args); }

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_Weights>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"weight output"; }
