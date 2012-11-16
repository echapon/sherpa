#!/usr/bin/python

import pygtk
import gtk



import GUI/ParameterBox.py      as parameterbox
import "GUI/Group.py"             as group
import "GUI/GenSection_GUI.py"    as gen_gui
import "GUI/BeamSection_GUI.py"   as beams_gui
import "GUI/PDFSection_GUI.py"    as pdfs_gui
import "GUI/ProcSection_GUI.py"   as procs_gui
import "GUI/MESection_GUI.py"     as mes_gui
import "GUI/ShowerSection_GUI.py" as shower_gui
import "GUI/SelSection_GUI.py"    as sels_gui
import "GUI/SoftSection_GUI.py"   as soft_gui
import "GUI/AnaSection_GUI.py"    as ana_gui

class runcard_gui_gtk():
    def __init__(self):
        self.initialize()
        self.visualSettings()
        self.constructNotebook()


    def initialize(self):
        self.parameters = parameterbox.ParameterBox()
        self.parameters.initialiseDefaults()
        self.gens  = gen_gui.generalsettings_gui(self.parameters,
                                                 "Run","Run",0)
        self.beams = beams_gui.beamsection_gui(self.parameters,
                                               "Beams","Beams",1)
        self.pdfs  = pdfs_gui.PDFsection_gui(self.parameters,
                                             "PDFs","PDFs",2)
        self.procs = procs_gui.Procsection_gui(self.parameters,
                                               "Processes","Processes",3)
        self.mes   = mes_gui.MEsection_gui(self.parameters,
                                           "ME Generators",
                                           "ME Generators",4)
        self.sels  = sels_gui.Selsection_gui(self.parameters,
                                             "Selectors",
                                             "Selectors",5)
        self.show  = shower_gui.Showersection_gui(self.parameters,
                                                  "Shower",
                                                  "Shower",6)
        self.soft  = soft_gui.Softsection_gui(self.parameters,
                                              "Nonperturbative",
                                              "Nonperturbative",7)
        self.ana   = ana_gui.Anasection_gui(self.parameters,
                                            "Analysis",
                                            "Analysis",8)

        self.guiparts = [self.gens,self.beams,self.pdfs,
                         self.procs,self.mes,self.sels,self.show,
                         self.soft,self.ana]


    def visualSettings(self):
        self.borderwidth = 20
        self.framexsize  = 700
        self.frameysize  = 600
        self.frameborder = 10

    def constructNotebook(self):
        self.window = gtk.Window()
        self.window.connect("delete-event", gtk.main_quit)
        self.window.set_border_width(10)
        table    = gtk.Table(3,6,False)
        self.window.add(table)
        notebook = gtk.Notebook() 
        notebook.set_tab_pos(gtk.POS_TOP)
        notebook.set_tab_pos((notebook.get_tab_pos()+2) %4)
        table.attach(notebook,0,6,0,1)
        notebook.show()
        sherpa      = group.groupsection()
        bufferframe = sherpa.getFrame(self.frameborder,
                                     self.framexsize,self.frameysize,
                                     self.update)
        bufferframe.show()
        bufferlabel = sherpa.getLabel()
        notebook.append_page(bufferframe,bufferlabel)

        for guipart in self.guiparts:
            bufferframe    = guipart.getFrame(self.frameborder,
                                              self.framexsize,self.frameysize,
                                              self.update)
            bufferframe.show()
            bufferlabel = gtk.Label(guipart.getName())
            notebook.append_page(bufferframe,bufferlabel)

        self.update("All")
        self.window.show_all()

        writeSbutton = gtk.Button(label="Write to Sherpa file")
        writeSbutton.connect("clicked", self.writebutton_clicked,"Sherpa")
        #writeHbutton = gtk.Button(label="Write to Herwig file")
        #writeHbutton.connect("clicked", self.writebutton_clicked,"Herwig")
        exitbutton = gtk.Button(label="Exit GUI without writing")
        exitbutton.connect("clicked", self.exitbutton_clicked)
        table.attach(writeSbutton, 1,2,1,2)
        #table.attach(writeHbutton, 2,3,1,2)
        table.attach(exitbutton, 2,3,1,2)
        writeSbutton.show()
        #writeHbutton.show()
        exitbutton.show()

    def update(self,mode):
        print 'In update with mode = %s' %(mode)
        if mode=="Beams" or mode=="All":
            collider = self.parameters.getBeamBox().getCollider()
            self.parameters.getPDFBox().initialiseDefaults(collider)
            self.pdfs.updateOptions()
            self.parameters.getProcBox().initialiseDefaults(collider)
            self.procs.updateOptions()
            self.parameters.getSoftBox().initialiseDefaults(collider)
            self.soft.updateOptions()
        if mode=="PDFs" or mode=="All":
            pdf = []
            for bid in range(0,2):
                pdf.append(self.parameters.getPDFBox().getType(bid))
                pdf.append(self.parameters.getPDFBox().getPDF(bid))
            # need to find a solution to select only valid tunes 
            # dependent on the pdfs
            # self.soft.updateModes(pdf)
        if mode=="Proc" or mode=="All":
            self.procs.updateOptions()
            self.mes.updateOptions()
            self.show.updateOptions(self.parameters.getProcBox().isNLO())
            self.sels.updateOptions(self.parameters.getBeamBox().getCollider(),
                                    self.parameters.getProcBox().getProcess(),
                                    self.parameters.getProcBox().getMinJets())


    def writebutton_clicked(self, widget, genstring):
        for guipart in self.guiparts:
            guipart.extractParameters()
        self.parameters.write(genstring)

    def exitbutton_clicked(self, widget):
        gtk.main_quit()

    

if __name__ == "__main__":
    window = runcard_gui_gtk()
    gtk.main()
