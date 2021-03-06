.. _External PDF:

************
External PDF
************


To use an external PDF (not included in LHAPDF) in Sherpa, you need to
provide an interface to your PDF in an external dynamic library. This
library is then loaded at runtime and it is possible within Sherpa to
access all PDFs included.

The simplest C++ code to implement your interface looks as follows

.. code-block:: c++

   #include "PDF/Main/PDF_Base.H"

   using namespace PDF;

   class Example_PDF: public PDF_Base {
   public:
     void Calculate(double x,double Q2)
     {
       // calculate values x f_a(x,Q2) for all a
     }
     double GetXPDF(const ATOOLS::Flavour a)
     {
       // return x f_a(x,Q2)
     }
     virtual PDF_Base *GetCopy()
     {
       return new Example_PDF();
     }
   };// end of class Example_PDF

   // this makes Example_PDF loadable in Sherpa
   DECLARE_PDF_GETTER(Example_PDF_Getter);
   PDF_Base *Example_PDF_Getter::operator()(const Parameter_Type &args) const
   { return new Example_PDF(); }
   // this eventually prints a help message
   void Example_PDF_Getter::PrintInfo
   (std::ostream &str,const size_t width) const
   { str<<"example PDF"; }
   // this lets Sherpa initialize and unload the library
   Example_PDF_Getter *p_get=NULL;
   extern "C" void InitPDFLib()
   { p_get = new Example_PDF_Getter("ExamplePDF"); }
   extern "C" void ExitPDFLib() { delete p_get; }

If the code is compiled into a library called libExamplePDFSherpa.so,
then this library is loaded dynamically in Sherpa using ``PDF_LIBRARY:
ExamplePDFSherpa`` either on the command line, or in the
:file:`Sherpa.yaml`. If the library is bound at compile time, like
e.g. in cmt, you may skip this step.  It is now possible to list all
accessible PDF sets by specifying :option:`SHOW_PDF_SETS: 1` on the
command line.

Finally Sherpa is instructed to retrieve the external PDF by
specifying :option:`PDF_SET: ExamplePDF` on the command line or in the
:file:`Sherpa.yaml`.
