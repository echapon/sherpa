#!/usr/bin/python2

# sys provides argv
import sys
sys.path.append('/home/silvan/work/bsm/lib/python2.7/site-packages/SHERPA-MC')

import Sherpa
import subprocess


def test_SherPY(_argv):
    
    #Args=ArgHandler.ArgHandler()
    Generator=Sherpa.Sherpa()

    try:
        # The C++ method "InitializeTheRun(int, char**)" needs a C++ array of pointers
        # to character variables (these are the command line options that the "main" program would normally provide).
        # sys provides a list with the arguments that were provided when executing the main script.
        # The Sherpa wrapper provides a typemap, that allows to directly pass the python argument list to the 
        # method Sherpa::InitializeTheRun

#        Generator.InitializeTheRun(len(_argv),_argv)
        Generator.InitializeTheRun(1,[])
        
        #File=Sherpa.EventFile.EventFile('TEST_FILE')
    
#        for i in range(1,File.NumEvents+1):
#            Event=File.GetEvent(i)
#            Amp=Sherpa.Cluster_Amplitude.Cluster_Amplitude.New()
#            Amp.SetNIn(len(Event.HardProcess.InParticles))
            # For ingoing particles, flavour and momentum need to be reversed
            # Vec4 is a C++ template class of which each instance needs to be wrapped and renamed explicitly.
            # Vec4D is a "double" instance of the Vec4 template (Vec4<double>). Further instances
            # can easily be added (see Vec4.i)
#            for Part in Event.HardProcess.InParticles:
#                Amp.CreateLegFromPyVec4D(-Part.FourMom, Part.Flavour.Bar())
#            for Part in Event.HardProcess.OutParticles:
#                Amp.CreateLegFromPyVec4D(Part.FourMom, Part.Flavour)
#            print "-------------------------------------------------------------------------------------------"
#            print "Hard process name generated from event file data   ", Sherpa.Process_Base.Process_Base.GenerateName(Amp)
#            print "Hard process generated from event file data        ", Event.HardProcess
#            print "Weigths read-in from event file                     ", Event.Weights
            # NOTE: In case the process is not found among the initialized processes, 
            # the following methods can throw. See 'Matrix_Element_Handler.i for details.
            # In contrast to the Initialization method Sherpa::InitializeTheRun, 
            # python exceptions are thrown (instead of Sherpa's builtin exception types)
            # because the methods below are python specific extensions of the C++ class
            # Matrix_Element_Hanlder
#            print "Matrix element computed from hard process          ", Generator.GetInitHandler().GetMatrixElementHandler().GetWeight(Amp)
            #print "Matrix element computed from hard process          ", Generator.GetInitHandler().GetMatrixElementHandler().GetProcess(Amp).Differential(Amp)
#            print "-------------------------------------------------------------------------------------------\n\n"
        
    # Sherpa might throw Exceptions that need to be taken care of. Therefore, 
    # the Python module "Exception" includes a wrapper for the Exception class
    # declared in "ATOOLS/Org/Exception.H."
    except Sherpa.Exception as exc:
        # "Message.cvar" stores all global C++ objects defined in "ATOOLS/Org/Message.*".
        # "Message.cvar.msg.Error()" returns a pointer to an ostream object
        # The Methos "Stream_Exception(std::ostream, ATOOLS::Exception)" is a wrapper
        # for the operator "std::ostream &operator<<(std::ostream &str,const ATOOLS::Exception &exception)"
        # declared in "ATOOLS/Org/Exception.H". Accordingly, this wrapper is part of the python module "Exception".
        Sherpa.Stream_Exception(Sherpa.cvar.msg.Error(), exc)
        # In order to get the output before further python commands are executed, the C++ ostream object
        # needs to be flushed. This can be done using the "Flush()" method or the "Endl()" mehtod.
        # "Endl()" ends the line AND flushes the ostream (internally via "<<std::endl;")
        Sherpa.cvar.msg.Error().Endl()
        
        # SHERPA's EXEPTION TYPES:
        #  normal_exit         = 1,
        #  unknown_option      = 2,
        #  inconsistent_option = 3,
        #  not_implemented     = 4,
        #  critical_error      = 5,
        #  fatal_error         = 6,
        #  missing_input       = 7,
        #  missing_module      = 8,
        #  unknown             = 0 
        # Exception type 1 is a normal exit, most likely AMEGIC finished writing library code
        #    if exc.Type()==1:
        #        subprocess.call("./makelibs")
                
        # Note: ostream input operator for Cluster_Amplitude is also wrapped and can be used
        # for example with the global Error ostream:
        # Cluster_Amplitude.Stream_Cluster_Amplitude(Message.cvar.msg.Error(), a)
        

test_SherPY(sys.argv)
