"""A ListCtrl (table) for calculating species masses from peak m/z's,
as well as storing species information such as charge and name.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx
from EditableListCtrl import EditableListCtrl
from lib import utils
from msClasses import Species
import AtroposGuiSettings
import gui.guiFunctions as gf

class ListCtrlAssigningSpecies(wx.Panel):

    def __init__(self,parent,yourself):
        wx.Panel.__init__(self,parent,-1)
        self.settings = AtroposGuiSettings.AtroposGuiSettings()

        self.listCtrl = EditableListCtrl(self)
        self.listCtrl.setColumnToIgnore(0)

        self.listCtrl.InsertColumn(0,"Option",width=80)
        self.listCtrl.InsertColumn(1,"Species",width=100)

        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0)
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)

        self.listCtrl.InsertStringItem(0, "" )
        self.listCtrl.SetStringItem(0, 0, "Name")

        self.listCtrl.InsertStringItem(1, "" )
        self.listCtrl.SetStringItem(1, 0, "Peak IDs")

        self.listCtrl.InsertStringItem(2, "" )
        self.listCtrl.SetStringItem(2, 0, "Peak m/zs")

        self.listCtrl.InsertStringItem(3, "" )
        self.listCtrl.SetStringItem(3, 0, "Mass")

        self.listCtrl.InsertStringItem(4, "" )
        self.listCtrl.SetStringItem(4, 0, "MassError")

        self.listCtrl.InsertStringItem(5, "" )
        self.listCtrl.SetStringItem(5, 0, "Charges of IDs")

        self.listCtrl.InsertStringItem(6, "" )
        self.listCtrl.SetStringItem(6, 0, "ToSimulate")
        self.listCtrl.SetStringItem(6, 1, "True")

        self.listCtrl.InsertStringItem(7, "" )
        self.listCtrl.SetStringItem(7, 0, "Peak FWHM")
        self.listCtrl.SetStringItem(7, 1, "10")

        self.listCtrl.InsertStringItem(8,"")
        self.listCtrl.SetStringItem(8, 0, "Charges to Simulate")


        self.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.OnEndLabelEdit, self.listCtrl)

        self.columnSpeciesNames = {}
        self.zsToSimulate = {}

        self.gui = yourself

    def setSettings(self,settings):
        """Set the AtroposGuiSettings object, so its functions and
        attributes can be accessed more readily.

        :parameter settings: AtroposGuiSettings() object
        """
        self.settings = settings

    def getCurrentColumn(self):
        # TODO(gns) - This hasn't been implemented yet
        # potentially a way to do this, though I don't know if it works or should be removed
        
        # return self.settings.getSpeciesCoumn()
        return 1

    def clearPeakIdFields(self):
        """Set the values for all cells in the peak ids row to ''."""
        cols = self.listCtrl.GetColumnCount()
        for column in xrange(cols):
            if column:
                self.listCtrl.SetStringItem(1,column,'')

    def getSpeciesToSimulate(self):
        """From the 'ToSimulate' column, get all of the true responses.

        :returns: List of names for the species to be simulated
        """
        cols = self.listCtrl.GetColumnCount()
        speciesToSimulate = []
        for i in range(cols-1):
            column = i+1
            toSim = self.listCtrl.GetItem(6,column).GetText()
            if utils.isBinaryResponse(toSim):
                toSim = utils.getBinaryReponse(toSim)
            else:
                toSim = False
            if toSim:
                speciesToSimulate.append(self.columnSpeciesNames[column])
        return speciesToSimulate

    def OnEndLabelEdit(self,event):
        """Context specific function for deciding what to do when a cell
        is edited depending on which row it belongs to.
        """
        row = event.m_itemIndex
        column = event.GetColumn()
        #===========================================================================
        # Species Name
        #===========================================================================
        if row == 0:
            newName = event.GetText()
            if not column in self.columnSpeciesNames:
                # name not set before
                self.columnSpeciesNames[column] = newName
            else:
                # name already set and changing it
                spD = self.settings.msPanel.ms.species
                spDsimulated = self.settings.msPanel.ms.simulatedSpecies

                if not newName in spD:
                    oldName = self.columnSpeciesNames[column]
                    # TODO (gns) - this isn't doing enough
                    # The problem must be coming from elsewhere (trying to access
                    # species 'Default09' for example)
                    # The problem occured when I reloaded the figure, so look into what
                    # happens when that button is pushed
                    spD[newName] = spD.pop(oldName)
                    spD[newName].name = newName

                    try:
                        spDsimulated[newName] = spDsimulated.pop(oldName)
                        spDsimulated[newName].name = newName
                    except:
                        spDsimulated[newName] = Species(newName)


                    self.columnSpeciesNames[column] = newName
                else:
                    print 'Provided species name is already in use'
        #===========================================================================
        # Peak finding
        #===========================================================================
        elif row == 1:
            # get text and set m/zs
            text = event.GetText()
            # Add default name if none given
            item = self.listCtrl.GetItem(0,column)
            name = item.GetText()
            if str(name) == '':
                name = 'Default%02d' %column
                self.listCtrl.SetStringItem(0,column,name)
                self.columnSpeciesNames[column] = name

            self.setPeakMzsAndSpecies(text,name, column)
        #===========================================================================
        # Entering mass
        #===========================================================================
        elif row == 3:
            '''
            use the mass to display theoretical charge states
            only works on unassigned Species # TODO(gns) - why?
            '''
            charges = self.listCtrlAssigningSpecies.listCtrl.GetItem(5,column).GetText()
            if charges == '': # TODO(gns) - why do the charges need to be blank?
                # TODO(gns) - warning dialog to say "Enter charge states before simulating"
                self.setManualMass(event.GetText())
            else:
                # Clear all the peak information
                # open warning box saying "Enter charge states before simulating"
                # Make sure that adding the charges recalculates the peak m/zs
                # and sets them properly
                event.Veto()
        #===========================================================================
        # Changing charges
        #===========================================================================
        elif row == 5:
            event.Veto()
        #===========================================================================
        # Peak FWHM
        #===========================================================================
        elif row == 7:
            # get current species name
            item = self.listCtrl.GetItem(0,column)
            name = item.GetText()
            peakFwhm = event.GetText()
            print peakFwhm, 'peakfwhm'
            if utils.isNumber(peakFwhm):
                peakFwhm = float(peakFwhm)
                if name in self.settings.msPanel.ms.species.keys():
                    self.settings.msPanel.ms.species[name].peakFwhm = peakFwhm
            else:
                self.listCtrl.SetStringItem(row,column, "10")
                event.Veto()
                print 'Peak FWHM value not a number\nReverting to default (10)'

        #===========================================================================
        # Charges to simulate
        #===========================================================================
        elif row == 8:
            # get current species name
            item = self.listCtrl.GetItem(0,column)
            name = item.GetText()

            zstext = event.GetText()
            zs = utils.getHyphenCommaList(zstext)
            if zs:
                #self.zsToSimulate[column] = zs
                if name != '':
                    self.settings.msPanel.ms.species[name].charges = zs
            else:
                print 'Invalid list of charges: %s' %zstext
                event.Veto()


    def setManualMass(self,massText):
        """For use with plotting theoretical charge state m/z values.

        :parameter massText: Mass of species as text or number
        """
        try:
            self.gui.msPanel.temporaryTheoMzs(massText)
        except:
            event.Veto()
            gf.warningDialog('Only enter numbers here!')

    def addSpeciesColumn(self):
        """Add an additional species to the ListCtrl."""
        cols = self.listCtrl.GetColumnCount()
        self.listCtrl.InsertColumn(cols,"Species %d" %cols,width=100)

        # preset values for Simulate and PeakFwhm
        columnI = self.listCtrl.GetColumnCount() -1
        self.listCtrl.SetStringItem(6,columnI,'True')
        self.listCtrl.SetStringItem(7,columnI,'10')


    def setPeakMzsAndSpecies(self,peakIdsComplexString,name,column=-1):
        """Set the peak m/zs corresponding to the peak ids entered in the
        ListCtrl. Also creates the Species() object and sets that.

        :parameter peakIdsComplexString: String containing peak ids can include commas and hyphens e.g. '4-7,9'.
        :parameter name: Species name
        :parameter column: The column to write to. Use -1 to use currently selected column
        """
        # get column to use
        if column == -1:
            self.getCurrentColumn()

        # get and set mz values
        peakIds = utils.getHyphenCommaList(peakIdsComplexString)
        gPeaks = self.settings.msPanel.ms.getgPeaksFromIds(peakIds)
        mzs =  [ "%.2f" %gPeaks[k][0] for k in gPeaks.keys()]
        mzString = ''
        for i,v in enumerate(mzs):
            if not i:
                mzString = v
            else:
                mzString += ', '+v
        self.listCtrl.SetStringItem(2,column,mzString)

        # get Peak FWHM value
        peakFwhm = self.listCtrl.GetItem(7,column).GetText()
        if utils.isNumber(peakFwhm):
            peakFwhm = float(peakFwhm)
        else:
            peakFwhm = 10
            print 'Peak Fwhm not a number\nUsing 10 ...'

        # set species
        spOb = Species(name)
        spOb.setSpecies(gPeaks, peakFwhm)
        self.settings.msPanel.ms.addSpecies(spOb,allowReplace=1)
        mass,error,zs = spOb.calculateMassAndCharges([float(mz) for mz in mzs])

        self.listCtrl.SetStringItem(3,column,"{:,}".format(round(mass,2)))
        self.listCtrl.SetStringItem(4,column,"{:,}".format(round(error,2)))

        # display zs
        if self.zsEmpty(column):
            self.setZs(zs,column)

        # put default calculated zs into ListCtrl
        self.listCtrl.SetStringItem(8,column,str(spOb.charges)[1:-1])

        # record current column
        self.settings.setSpeciesColumn(column)


    def setZs(self,zs,column=-1):
        """Set the charge values in the corresponding ListCtrl row.

        :parameter zs: List of charges
        :parameter column: The column to write to. Use -1 to use currently selected column
        """
        # get column to use
        if column == -1:
            column = self.getCurrentColumn()

        self.listCtrl.SetStringItem(5,column,str(zs)[1:-1])

    def zsEmpty(self,column=-1):
        """Check if the charges cell for this row is empty.

        :parameter column: The column to check. Use -1 to use currently selected column
        """
        # get column to use
        if column == -1:
            column = self.getCurrentColumn()
        currentzs = self.listCtrl.GetItem(5,column).GetText()
        if currentzs == '':
            return True
        else:
            return False
