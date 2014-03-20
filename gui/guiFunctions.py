import wx,os,re
from imClasses import Im


def checkIfAmphiProject(path):
    amphiProject = True
    if not path.strip("/")[-2:] == '.a':
        amphiProject = False
    filesNeeded = ['MassMobility.amphi','MassMobilityXaxis.amphi','MassMobilityYaxis.amphi']
    for f in filesNeeded:
         if not f in os.listdir(path):
            amphiProject = False
    if not amphiProject:
        print 'Not an appropriate amphitrite file: %s' %path

    return amphiProject


def openCalibrationFile(yourself):
    dlg = wx.FileDialog(yourself,message='Open Calibration File',
                        wildcard="Calibration file (*.calibration)|*.calibration",
                        style=wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
        path = dlg.GetPath()
    else:
        path = False

    return path

def openAtroposFile(yourself):
    dlg = wx.FileDialog(yourself,message='Open Atropos File',
                        wildcard="Calibration file (*.afit)|*.afit",
                        style=wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
        path = dlg.GetPath()
    else:
        path = False

    return path

def checkIfNumberTextCtrl(textCtrl):
    val = textCtrl.GetValue()
    print val
    if val == "": # stops it bugging out when blank
        val = 0.0
    try:
        val = float(val)
        return val
    except:
        return 'string error (avoid "0" being False problems)'

def checkAndSetIfNumberTextCtrl(textCtrl,originalValue):
    val = checkIfNumberTextCtrl(textCtrl)
    if type(val).__name__ != 'str':
        return float(val)
    else:
        message = 'Please only enter numbers!'
        warningDialog(message)
        textCtrl.SetValue(str(originalValue))    

def warningDialog(message):    
    dial = wx.MessageDialog(None,message, 'Warning',
                            wx.OK | wx.ICON_EXCLAMATION)
    dial.ShowModal()

