
GITTAG = AddOns/Delphes
include ../../git.make

DELPHES_SOURCES   =	\
		Git_Info.C \
		Output_Delphes.C

if DELPHES_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaDelphesOutput.la

libSherpaDelphesOutput_la_SOURCES     =	$(DELPHES_SOURCES)
libSherpaDelphesOutput_la_LIBADD      =	@CONDITIONAL_DELPHESLIBS@
libSherpaDelphesOutput_la_CPPFLAGS    =	$(AM_CPPFLAGS) \
	@CONDITIONAL_DELPHESINCS@ @CONDITIONAL_ROOTINCS@ @CONDITIONAL_HEPMC2INCS@
else
DELPHES_EXTRA_DIST    = $(DELPHES_SOURCES)
endif

EXTRA_DIST = $(DELPHES_EXTRA_DIST)
