#====================================================================
#   Code Performs Analysis on Data Obtained by simulation   
#    and reconstruction for RPV bencmark 1 final states 
#   and dark photon model
#
#   Created by Nicolas Stylianou and Danny Galbinski
#
#===================================================================

import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
from array import array
import shipRoot_conf
import numpy as np
import shipDet_conf
import shipVeto
shipRoot_conf.configure()

import rpvsusy
from ROOT import TCanvas, TH1D, TF1, TMultiGraph,TGraphErrors, THStack, TFile
from ROOT import kBlack, kBlue, kRed, kGreen, kGray, kMagenta
from ROOT import gROOT, gPad, gStyle
from inspect import currentframe
import datetime


debug = False
loop=True 
PDG = ROOT.TDatabasePDG.Instance()
inputFile  = None
geoFile    = None
dy         = None
nEvents    = 9999999
fiducialCut = True
chi2Cut  = 4
measCut = 25
ecalCut = 0.150 #GeV
docaCut = 2. #cm
ipCut = 250 #cm

#e = 2.718281828459   # Euler's number
currentDate = datetime.datetime.now().strftime('%y_%m_%d_%H%M')
polyFit1 = TF1('polyFit1','pol3')
polyFit2 = TF1('polyFit2','pol3')

##############################################
#####  SETS INPUT OPTIONS AND ARGUMENTS  #####
try:
        opts, args = getopt.getopt(sys.argv[1:], 'n:f:g:Y', ['nEvents=','geoFile='])#command line options
except getopt.GetoptError:
        # print help information and exit:
        print (' enter file name')
        sys.exit()
for o, a in opts:
        if o in ('-f',):
            inputFile = a           #does this
        if o in ('-g', '--geoFile',):
            geoFile = a             #does this
        if o in ('-Y',):
            dy = float(a)           #doesnt do this
        if o in ('-n', '--nEvents=',):
            nEvents = int(a)        #doesnt do this

###############################################
########## SETS TREE FROM INPUT FILE ##########
if not inputFile.find(',')<0 :
    sTree = ROOT.TChain('cbmsim') 
    for x in inputFile.split(','):
        if x[0:4] == '/eos':         #doesnt do this
            sTree.AddFile('root://eoslhcb.cern.ch/'+x)
        else: sTree.AddFile(x)       #doesnt do this
elif inputFile[0:4] == '/eos':       #doesnt do this
    eospath = 'root://eoslhcb.cern.ch/'+inputFile
    f = ROOT.TFile.Open(eospath)
    sTree = f.cbmsim
else:                                #does this
    f = ROOT.TFile(inputFile)
    sTree = f.cbmsim

#######################################
########## LOAD ECAL GEOMETRY #########
if not geoFile:
 geoFile = inputFile.replace('ship.','geofile_full.').replace('_rec.','.')  #doesnt do this
if geoFile[0:4] == '/eos':                                                  #doesnt do this
  eospath = 'root://eoslhcb.cern.ch/'+geoFile
  fgeo = ROOT.TFile.Open(eospath)
else:                                                                       #does this
    fgeo = ROOT.TFile(geoFile)
sGeo = fgeo.FAIRGeom

if not fgeo.FindKey('ShipGeo'):                 #doesnt do this
 # old geofile, missing Shipgeo dictionary
 if sGeo.GetVolume('EcalModule3') :  
     ecalGeoFile = 'ecal_ellipse6x12m2.geo'     #doesnt do this
 else: 
     ecalGeoFile = 'ecal_ellipse5x10m2.geo'     #doesnt do this
 print ('found ecal geo for ',ecalGeoFile)
 # re-create geometry and mag. field
 if not dy:
  # try to extract from input file name
  tmp = inputFile.split('.')                    #doesnt do this
  try:
    dy = float( tmp[1]+'.'+tmp[2] )             #doesnt do this
  except:
    dy = 10.                                    #doesnt do this
 ShipGeo = ConfigRegistry.loadpy('$FAIRSHIP/geometry/geometry_config.py', Yheight = dy, EcalGeoFile = ecalGeoFile )
else:                                           #does this
 # new geofile, load Shipgeo dictionary written by run_simScript.py
  upkl    = Unpickler(fgeo)
  ShipGeo = upkl.load('ShipGeo') #load is def in def Unpickler
  ecalGeoFile = ShipGeo.ecal.File
  dy = ShipGeo.Yheight/u.m

##########################################
#########  BEGIN CREATE GEOMETRY #########
run = ROOT.FairRunSim()
modules = shipDet_conf.configure(run,ShipGeo)

gMan  = ROOT.gGeoManager
geoMat =  ROOT.genfit.TGeoMaterialInterface()
ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
volDict = {}
i=0

for x in ROOT.gGeoManager.GetListOfVolumes():
 volDict[i]=x.GetName()
 i+=1
 
bfield = ROOT.genfit.BellField(ShipGeo.Bfield.max ,ShipGeo.Bfield.z,2, ShipGeo.Yheight/2.)
fM = ROOT.genfit.FieldManager.getInstance()
fM.init(bfield)

# prepare veto decisions
veto = shipVeto.Task(sTree)
vetoDets={}

# fiducial cuts
vetoStation = ROOT.gGeoManager.GetTopVolume().GetNode('Veto_5')
vetoStation_zDown = vetoStation.GetMatrix().GetTranslation()[2]+vetoStation.GetVolume().GetShape().GetDZ()
T1Station = ROOT.gGeoManager.GetTopVolume().GetNode('Tr1_1')
T1Station_zUp = T1Station.GetMatrix().GetTranslation()[2]-T1Station.GetVolume().GetShape().GetDZ()

# calculate z front face of ecal, needed later?
top = ROOT.gGeoManager.GetTopVolume()
z_ecal      = top.GetNode('Ecal_1').GetMatrix().GetTranslation()[2]
z_ecalFront = z_ecal + top.GetNode('Ecal_1').GetVolume().GetShape().GetDZ()
z_ecalBack  = z_ecal + top.GetNode('Ecal_1').GetVolume().GetShape().GetDZ()

# initialize ecalStructure
caloTasks = []
sTree.GetEvent(0)
ecalGeo = ecalGeoFile+'z'+str(ShipGeo.ecal.z)+'.geo'
if not ecalGeo in os.listdir(os.environ['FAIRSHIP']+'/geometry'): shipDet_conf.makeEcalGeoFile(ShipGeo.ecal.z,ShipGeo.ecal.File)
ecalFiller = ROOT.ecalStructureFiller('ecalFiller', 0,ecalGeo)
ecalFiller.SetUseMCPoints(ROOT.kTRUE)
ecalFiller.StoreTrackInformation()
ecalStructure = ecalFiller.InitPython(sTree.EcalPointLite)
caloTasks.append(ecalFiller)
if sTree.GetBranch('EcalReconstructed'):
 calReco = False
 sTree.GetEvent(0)
 ecalReconstructed = sTree.EcalReconstructed
else:
 calReco = True
 print ('setup calo reconstruction of ecalReconstructed objects')
# Calorimeter reconstruction
 #GeV -> ADC conversion
 ecalDigi=ROOT.ecalDigi('ecalDigi',0)
 ecalPrepare=ROOT.ecalPrepare('ecalPrepare',0)
 ecalStructure     = ecalFiller.InitPython(sTree.EcalPointLite)
 ecalDigi.InitPython(ecalStructure)
 caloTasks.append(ecalDigi)
 ecalPrepare.InitPython(ecalStructure)
 caloTasks.append(ecalPrepare)
 # Cluster calibration
 ecalClusterCalib=ROOT.ecalClusterCalibration('ecalClusterCalibration', 0)
 #4x4 cm cells
 ecalCl3PhS=ROOT.TFormula('ecalCl3PhS', '[0]+x*([1]+x*([2]+x*[3]))')
 ecalCl3PhS.SetParameters(6.77797e-04, 5.75385e+00, 3.42690e-03, -1.16383e-04)
 ecalClusterCalib.SetStraightCalibration(3, ecalCl3PhS)
 ecalCl3Ph=ROOT.TFormula('ecalCl3Ph', '[0]+x*([1]+x*([2]+x*[3]))+[4]*x*y+[5]*x*y*y')
 ecalCl3Ph.SetParameters(0.000750975, 5.7552, 0.00282783, -8.0025e-05, -0.000823651, 0.000111561)
 ecalClusterCalib.SetCalibration(3, ecalCl3Ph)
#6x6 cm cells
 ecalCl2PhS=ROOT.TFormula('ecalCl2PhS', '[0]+x*([1]+x*([2]+x*[3]))')
 ecalCl2PhS.SetParameters(8.14724e-04, 5.67428e+00, 3.39030e-03, -1.28388e-04)
 ecalClusterCalib.SetStraightCalibration(2, ecalCl2PhS)
 ecalCl2Ph=ROOT.TFormula('ecalCl2Ph', '[0]+x*([1]+x*([2]+x*[3]))+[4]*x*y+[5]*x*y*y')
 ecalCl2Ph.SetParameters(0.000948095, 5.67471, 0.00339177, -0.000122629, -0.000169109, 8.33448e-06)
 ecalClusterCalib.SetCalibration(2, ecalCl2Ph)
 caloTasks.append(ecalClusterCalib)
 ecalReco=ROOT.ecalReco('ecalReco',0)
 caloTasks.append(ecalReco)
# Match reco to MC
 ecalMatch=ROOT.ecalMatch('ecalMatch',0)
 caloTasks.append(ecalMatch)
 ecalCalib         = ecalClusterCalib.InitPython()
 ecalReconstructed = ecalReco.InitPython(sTree.EcalClusters, ecalStructure, ecalCalib)
 ecalMatch.InitPython(ecalStructure, ecalReconstructed, sTree.MCTrack)


import TrackExtrapolateTool

#######################################
##########  HIST FUNCTIONS  ###########
h = {}
graph = {}
def create_Hists(partList):
    dictionList={partList[1]:[kBlue+3,kBlue-2], partList[2]:[kRed+3,kRed-2]}
    if not partList[3] == None:
        dictionList.update({partList[3]:[kGreen+4,kGreen+3]})
    edgesarray=[]
    edgesarray.append(0)
    for binNumber in range(0,40):
        edgesarray.append(edgesarray[binNumber]+ 0.015)
    for binNumber in range(40,86):
        edgesarray.append(edgesarray[binNumber]+0.045)
    
    ###############################
    ####  Daughter Histograms  ####
    for partName,val in dictionList.items():
        #print(partName)
        h[partName + 'StrawTime'] = TH1D(partName + 'StrawTime',partName + ' ; ' + partName + ' Straw Time [ns] ; No. of Particles',300,321.7,323.5)                       # straw time
        h[partName + 'StrawTime'].SetLineColor(val[0])
        h[partName + 'StrawTime'].SetFillColor(val[1])
        h[partName + 'EcalTime'] = TH1D(partName + 'EcalTime',partName + ' ; ' + partName + ' Ecal Time [ns] ; No. of Particles',300,359.7,361.2)                          # ecal time
        h[partName + 'EcalTime'].SetLineColor(val[0])
        h[partName + 'EcalTime'].SetFillColor(val[1])
        h[partName + 'DirDeltaTimeSmeared'] = TH1D(partName + 'DirDeltaTimeSmeared',partName + ' ; ' + partName + ' Straw-ECAL Smeared Time of Flight (directly) [ns] ; No. of Particles',300,37.8,38.4)# smeared time of flight
        h[partName + 'DirDeltaTimeSmeared'].SetLineColor(val[0])
        h[partName + 'DirDeltaTimeSmeared'].SetFillColor(val[1])
        h[partName + 'DirDeltaTime'] = TH1D(partName + 'DirDeltaTime',partName + ' ; ' + partName + ' Straw-ECAL Time of Flight (directly) [ns] ; No. of Particles',300,37.8,38.4)# time of flight
        h[partName + 'DirDeltaTime'].SetLineColor(val[0])
        h[partName + 'DirDeltaTime'].SetFillColor(val[1])
        h[partName + 'StrawHits'] = TH1D(partName + 'StrawHits',partName + ' ; ' + partName + ' No. of hits in straw tubes ; Position [cm]',300,25,50)                     #number of straw hits
        h[partName + 'StrawHits'].SetLineColor(val[0])
        h[partName + 'StrawHits'].SetFillColor(val[1])
        h[partName + 'StrawHitsMom'] = TH1D(partName + 'StrawHitsMom',partName + ' ; ' + partName + ' z-momentum through straw tubes (for particular event) [GeV/c] ; No. of Hits',500,0,45)     #momenta of straw hits
        h[partName + 'StrawHitsMom'].SetLineColor(val[0])
        h[partName + 'StrawHitsMom'].SetFillColor(val[1])
        h[partName + 'FlightLen'] = TH1D(partName + 'FlightLen',partName + ' ; ' + partName + ' Straw-ECAL Straight Flight Lenght [cm] ; No. of Particles',300,11.375,11.42)# flight Length
        h[partName + 'FlightLen'].SetLineColor(val[0])
        h[partName + 'FlightLen'].SetFillColor(val[1])
        h[partName + 'FlightLenImproved'] = TH1D(partName + 'FlightLenImproved',partName + ' ; ' + partName + ' Straw-ECAL Curved Flight Lenght [cm] ;  No. of Particles',300,11.375,11.42)      # corrected flight Length        
        h[partName + 'FlightLenImproved'].SetLineColor(val[0])
        h[partName + 'FlightLenImproved'].SetFillColor(val[1])
        h[partName + 'FlightLenDelta'] = TH1D(partName + 'FlightLenDelta',partName + ' ; ' + partName + ' Difference between straight path and better approximation [cm] ; No. of Particles',300,0,0.001)# delta flight Length        
        h[partName + 'FlightLenDelta'].SetLineColor(val[0])
        h[partName + 'FlightLenDelta'].SetFillColor(val[1])
        h[partName + 'SpeedSmeared'] = TH1D(partName + 'SpeedSmeared',partName + ' ; ' + partName + ' Smeared Beta value ; No. of Particles',300,0.997,1.001)               # smeared speed
        h[partName + 'SpeedSmeared'].SetLineColor(val[0])
        h[partName + 'SpeedSmeared'].SetFillColor(val[1])
        h[partName + 'Speed'] = TH1D(partName + 'Speed',partName + ' ; ' + partName + ' Beta value ; No. of Particles',300,0.997,1.001)                                     # speed
        h[partName + 'Speed'].SetLineColor(val[0])
        h[partName + 'Speed'].SetFillColor(val[1])
        h[partName + 'StrawMom'] = TH1D(partName + 'StrawMom',partName + ' ; ' + partName + ' Straw Momentum [GeV/c] ; No. of Particles',300,-0.05,120.)                    # straw momentum
        h[partName + 'StrawMom'].SetLineColor(val[0])
        h[partName + 'StrawMom'].SetFillColor(val[1])
        h[partName + 'EcalMom'] = TH1D(partName + 'EcalMom',partName + ' ; ' + partName + ' Ecal Momentum [GeV/c]; No. of Particles',300,-0.05,120.)                        # ecal  momentum
        h[partName + 'EcalMom'].SetLineColor(val[0])
        h[partName + 'EcalMom'].SetFillColor(val[1])
        h[partName + 'DeltaMom'] = TH1D(partName + 'DeltaMom',partName + ' ; ' + partName + ' Straw-Ecal Momentum [GeV/c]; No. of Particles',300,0.02,0.13)                 # delta momentum
        h[partName + 'DeltaMom'].SetLineColor(val[0])
        h[partName + 'DeltaMom'].SetFillColor(val[1])
        h[partName + 'RecoMom'] = TH1D(partName + 'RecoMom',partName + ' ; ' + partName + ' Reco Momentum [GeV/c] ; No. of Particles',300,-0.05,120.)                       # reco  momentum
        h[partName + 'RecoMom'].SetLineColor(val[0])
        h[partName + 'RecoMom'].SetFillColor(val[1])
        h[partName + 'TrueMom'] = TH1D(partName + 'TrueMom',partName + ' ; ' + partName + ' True Momentum [GeV/c] ; No. of Particles',300,-0.05,120.)                       # true  momentum
        h[partName + 'TrueMom'].SetLineColor(val[0])
        h[partName + 'TrueMom'].SetFillColor(val[1])
        h[partName + 'RecoMass'] = TH1D(partName + 'RecoMass',partName + ' ; ' + partName + ' Reco Mass [GeV/c2]; No. of Particles',300,0.,0.6)                             # reco  mass
        h[partName + 'RecoMass'].SetLineColor(val[0])
        h[partName + 'RecoMass'].SetFillColor(val[1])
        h[partName + 'TrueMass'] = TH1D(partName + 'TrueMass',partName + ' ; ' + partName + ' True Mass [GeV/c2] ; No. of Particles',300,0.,0.6)                            # true  mass
        h[partName + 'TrueMass'].SetLineColor(val[0])
        h[partName + 'TrueMass'].SetFillColor(val[1])
        h[partName + 'MassT'] = TH1D(partName + 'MassT',partName + ' ; ' + partName + ' Time Deduced Mass / [GeV/c2]; No. of Particles ',150,0.,3.)                         # time deduced mass
        h[partName + 'MassT'].SetLineColor(val[0])
        h[partName + 'MassT'].SetFillColor(val[1])
        h[partName + 'SmearedMass'] = TH1D(partName + 'SmearedMass',partName + ' ; ' + partName + ' Smeared Mass [GeV/c2]; No. of Particles',85,array('d',edgesarray))      # smrd  mass
        h[partName + 'SmearedMass'].SetLineColor(val[0])
        h[partName + 'SmearedMass'].SetFillColor(val[1])
        h[partName + 'ProbMeasr'] = TH1D(partName + 'ProbMeasr',partName + ' ; Mass [GeV/c2] ; Prob(particle = ' + partName + ')',85,array('d',edgesarray))                 # ID Prob
                              
    h['TotalSmearedMass'] = TH1D('TotalSmearedMass','Smeared Mass ; Smeared Mass [GeV/c2] ; No. of Particles',85,array('d',edgesarray))                             # Total mass
    h['TotalSmearedMass'].SetLineColor(1)
    h['TotalSmearedMass'].SetFillColor(kGray+2)
    h['StrawTime'] = THStack('StackStrawTime','Gaussian Straw t measurement ; Time [ns] ; No. of Particles')                                                        # straw time
    h['EcalTime'] = THStack('StackEcalTime','Gaussian Ecal t measurement ; Time [ns] ; No. of Particles')                                                           # ecal time
    h['DirDeltaTimeSmeared'] = THStack('StackDirDeltaTimeSmeared','Straw-ECAL Smeared Time of Flight (directly) ; Time of Flight [ns] ; No. of Particles')          # smeared time of flight
    h['DirDeltaTime'] = THStack('StackDirDeltaTime','Straw-ECAL Time of Flight (directly) ; Time of Flight [ns] ; No. of Particles')                                # time of flight
    h['StrawHits'] = THStack('StackStrawHits','No. of hits in straw tubes ; Position [cm] ; No. of Particles')                                                      # number of straw hits
    h['StrawHitsMom'] = THStack('StackStrawHitsMom','z-momentum through straw tubes (for particular event) ; Momentum [GeV/c] ; No. of Particles')                  # momenta of straw hits
    h['FlightLen'] = THStack('StackFlightLen','Straw-ECAL Straight Flight Lenght ; Flight Length [cm] ; No. of Particles')                                          # flight Length
    h['FlightLenImproved'] = THStack('StackFlightLenImproved','Straw-ECAL Curved Flight Lenght ; Flight Length [cm] ; No. of Particles')                            #corrected flight Length
    h['FlightLenDelta'] = THStack('StackFlightLenDelta','Difference between straight path and better approximation ; Delta Flight Length [cm] ; No. of Particles')  # delta flight Length
    h['SpeedSmeared'] = THStack('StackSpeedSmeared','Smeared Beta value ; Smeared Beta value ; No. of Particles')                                                   # smeared speed
    h['Speed'] = THStack('StackSpeed','Beta value ; Beta value ; No. of Particles')                                                                                 # speed
    h['StrawMom'] = THStack('StackStrawMom','Straw Momentum ; Straw Momentum [GeV/c] ; No. of Particles')                                                           # straw momentum
    h['EcalMom'] = THStack('StackEcalMom','Ecal Momentum ; Ecal Momentum [GeV/c] ; No. of Particles')                                                               # ecal  momentum
    h['DeltaMom'] = THStack('StackDeltaMom','Straw-Ecal Momentum ; Straw-Ecal Momentum [GeV/c] ; No. of Particles')                                                 # delta momentum
    h['RecoMom'] = THStack('StackRecoMom','Reco Momentum ; Reco Momentum [GeV/c] ; No. of Particles')                                                               # reco  momentum
    h['TrueMom'] = THStack('StackTrueMom','True Momentum ; True Momentum [GeV/c] ; No. of Particles')                                                               # true  momentum
    h['RecoMass'] = THStack('StackRecoMass','Reco Mass ; Reco Mass [GeV/c2] ; No. of Particles')                                                                    # reco  mass
    h['TrueMass'] = THStack('StackTrueMass','True Mass ; True Mass [GeV/c2] ; No. of Particles')                                                                    # true  mass
    h['SmearedMass'] = THStack('StackSmearedMass','Smeared Mass ; Smeared Mass [GeV/c2] ; No. of Particles')                                                        # smrd  mass
    h['ProbMeasr'] = THStack('StackProbMeasr','Probs identifying Particle ; Mass [GeV/c2] ; Prob(particle ID correct )')                                            # ID Prob


    #################################
    ####  Neutralino Histograms  ####
    h[partList[0] + 'TrueMass'] = TH1D(partList[0] + 'TrueMass','Monte Carlo Mass ; Invariant mass [GeV/c2] ; No. of Particles',300,0.99,1.01)                            # true mass
    h[partList[0] + 'TrueMass'].SetLineColor(1)
    h[partList[0] + 'TrueMass'].SetFillColor(kGray+2)
    h[partList[0] + 'RecoMass'] = TH1D(partList[0] + 'RecoMass','Reconstructed Mass ; Invariant mass [GeV/c2] ; No. of Particles',300,0.97,1.03)                          # reco mass
    h[partList[0] + 'RecoMass'].SetLineColor(1)
    h[partList[0] + 'RecoMass'].SetFillColor(kGray+2)
    h[partList[0] + 'TrueMom'] = TH1D(partList[0] + 'TrueMom','True (red) & Reco. (blue) Momentum ; Momentum [GeV/c] ; No. of Particles',100,0.,180.)                     # true momentum 
    h[partList[0] + 'TrueMom'].SetLineColor(1)
    h[partList[0] + 'TrueMom'].SetFillColor(kGray+2)
    h[partList[0] + 'RecoMom'] = TH1D(partList[0] + 'RecoMom','Reconstructed Momentum ; Momentum [GeV/c] ; No. of Particles',300,0.,180.)                                 # reco momentum
    h[partList[0] + 'RecoMom'].SetLineColor(1)
    h[partList[0] + 'RecoMom'].SetFillColor(kGray+2)
    h[partList[0] + 'DeltaMom'] = TH1D(partList[0] + 'DeltaMom','True/Reco Momentum Difference ; Momentum Difference [GeV/c] ; No. of Particles',300,-3.,3)               # true-reco momentum difference
    h[partList[0] + 'DeltaMom'].SetLineColor(1)
    h[partList[0] + 'DeltaMom'].SetFillColor(kGray+2)
    h[partList[0] + 'Beta'] = TH1D(partList[0] + 'Beta','Reconstructed Neutralino Beta in Z-direction; Beta; No. of Particles',100,0.994,1)
    h[partList[0] + 'Beta'].SetLineColor(1)
    h[partList[0] + 'Beta'].SetFillColor(kGray+2)
    h[partList[0] + 'Gamma'] = TH1D(partList[0] + 'Gamma','Reconstructed Neutralino Gamma in Z-direction; Gamma; No. of Particles',100,0,200)
    h[partList[0] + 'Gamma'].SetLineColor(1)
    h[partList[0] + 'Gamma'].SetFillColor(kGray+2)
    h[partList[0] + 'Theta'] = TH1D(partList[0] + 'Theta','Angle between neutralino momentum and beam line; Theta / [mrad]; No. of Particles',100,0,50)
    h[partList[0] + 'Theta'].SetLineColor(1)
    h[partList[0] + 'Theta'].SetFillColor(kGray+2)

    #######################
    ####  Veto Checks  ####
    h['IP_target'] = TH1D('IP_target','Impact parameter to target; Impact Parameter [cm]; Frequency',300,0,10)
    h['IP_target'].SetLineColor(kMagenta+4)
    h['IP_target'].SetFillColor(kMagenta-1)
    h['ecalE'] = TH1D('ecalE','Energy deposited in ECAL ; Energy [GeV/c2] ; Frequency',300,0,100)
    h['ecalE'].SetLineColor(kMagenta+4)
    h['ecalE'].SetFillColor(kMagenta-1)
    h['doca'] = TH1D('doca','Distance of closest approach between muon and kaon tracks ; Distance [cm] ; Frequency',300,0,3)
    h['doca'].SetLineColor(kMagenta+4)
    h['doca'].SetFillColor(kMagenta+-1)
    h['nmeas'] = TH1D('nmeas','No. of measurements in fitted tracks (ndf) ; ndf ; No. of tracks',300,0,50)
    h['nmeas'].SetLineColor(kMagenta+4)
    h['nmeas'].SetFillColor(kMagenta-1)
    h['Chi2'] = TH1D('Chi2','Fitted Tracks Chi Squared ; Reduced Chi Squared ; Frequency',300,0,3)
    h['Chi2'].SetLineColor(kMagenta+4)
    h['Chi2'].SetFillColor(kMagenta-1)
    h['recovertex'] = TH1D('recovertex','Reconstructed neutralino decay vertex z-coordinate ; Z  [cm] ; Frequency',100,-4000,4000)
    h['recovertex'].SetLineColor(kMagenta+4)
    h['recovertex'].SetFillColor(kMagenta-1)

    #ut.bookHist(h,'photonE','Photon Energy Distribution in ECAL',150,0,100)


    #ut.bookHist(h,partList[0] + '_no_iter','Reconstructed Mass (without track iterations)',500,0.,2.)   # reco mass(without track itrns)
    #ut.bookHist(h,'normdistr','Gaussian Distribution',500,-0.05,0.05)                               #
    #ut.bookHist(h,'smearedmass1','Time Smeared Neutralino Mass',500,0.,2.)
    #ut.bookHist(h,'smearedmass2','Time Smeared Neutralino Mass',500,0.,2.)
    #ut.bookHist(h,'smearedP1','Time Smeared Neutralino Momentum P1(red) P2(blue)',500,0.,300.)
    #ut.bookHist(h,'smearedP2','Time Smeared Neutralino Momentum',500,0.,300.)

    print('Created Histograms')
    print(partList)
    print('\n')
def makePlots2(partList):
    
    key='DAUGHTERS'

    title='Time and velocity plots'
    h[key + '_TV']=TCanvas(key + '_TV',title,1300,800)
    h[key + '_TV'].Divide(3,2)
    cv = h[key + '_TV'].cd(1)
    h['StrawTime'].Add(h[partList[1] + 'StrawTime'])
    h['StrawTime'].Add(h[partList[2] + 'StrawTime'])
    if not partList[3]==None:
        h['StrawTime'].Add(h[partList[3] + 'StrawTime'])
    h['StrawTime'].Draw('nostack')

    cv = h[key + '_TV'].cd(2)
    h['EcalTime'].Add(h[partList[1] + 'EcalTime'])
    h['EcalTime'].Add(h[partList[2] + 'EcalTime'])
    if not partList[3]==None:
        h['EcalTime'].Add(h[partList[3] + 'EcalTime'])
    h['EcalTime'].Draw('nostack')
    

    cv = h[key + '_TV'].cd(3)
    h['DirDeltaTimeSmeared'].Add(h[partList[1] + 'DirDeltaTimeSmeared'])
    h['DirDeltaTimeSmeared'].Add(h[partList[2] + 'DirDeltaTimeSmeared'])
    if not partList[3]==None:
        h['DirDeltaTimeSmeared'].Add(h[partList[3] + 'DirDeltaTimeSmeared'])
    h['DirDeltaTimeSmeared'].Draw('nostack')

    cv = h[key + '_TV'].cd(4)
    h['FlightLen'].Add(h[partList[1] + 'FlightLen'])
    h['FlightLen'].Add(h[partList[2] + 'FlightLen'])
    if not partList[3]==None:
        h['FlightLen'].Add(h[partList[3] + 'FlightLen'])
    h['FlightLen'].Draw('nostack')

    cv = h[key + '_TV'].cd(5)
    h['Speed'].Add(h[partList[1] + 'Speed'])
    h['Speed'].Add(h[partList[2] + 'Speed'])
    if not partList[3]==None:
        h['Speed'].Add(h[partList[3] + 'Speed'])
    h['Speed'].Draw('nostack')

    #h[key + '_TV'].Print('DaughterTVProp'+ currentDate + '.png')

    title='Momenta and mass plots'
    h[key + '_MOM']=TCanvas(key + '_MOM',title,1300,800)
    h[key + '_MOM'].Divide(3,2)
    cv = h[key + '_MOM'].cd(1)
    h['StrawMom'].Add(h[partList[1] + 'StrawMom'])
    h['StrawMom'].Add(h[partList[2] + 'StrawMom'])
    if not partList[3]==None:
        h['StrawMom'].Add(h[partList[3] + 'StrawMom'])
    h['StrawMom'].Draw('nostack')


    cv = h[key + '_MOM'].cd(2)
    h['EcalMom'].Add(h[partList[1] + 'EcalMom'])
    h['EcalMom'].Add(h[partList[2] + 'EcalMom'])
    if not partList[3]==None:
        h['EcalMom'].Add(h[partList[3] + 'EcalMom'])
    h['EcalMom'].Draw('nostack')

    cv = h[key + '_MOM'].cd(3)
    h['RecoMom'].Add(h[partList[1] + 'RecoMom'])
    h['RecoMom'].Add(h[partList[2] + 'RecoMom'])
    if not partList[3]==None:
        h['RecoMom'].Add(h[partList[3] + 'RecoMom'])
    h['RecoMom'].Draw('nostack')
    
    cv = h[key + '_MOM'].cd(4)
    h['DeltaMom'].Add(h[partList[1] + 'DeltaMom'])
    h['DeltaMom'].Add(h[partList[2] + 'DeltaMom'])
    if not partList[3]==None:
        h['DeltaMom'].Add(h[partList[3] + 'DeltaMom'])
    h['DeltaMom'].Draw('nostack')

    cv = h[key + '_MOM'].cd(5)
    h['RecoMass'].Add(h[partList[1] + 'RecoMass'])
    h['RecoMass'].Add(h[partList[2] + 'RecoMass'])
    if not partList[3]==None:
        h['RecoMass'].Add(h[partList[3] + 'RecoMass'])
    h['RecoMass'].Draw('nostack')

    cv = h[key + '_MOM'].cd(6)
    h[partList[1] + 'SmearedMass'].Fit('landau')
    h[partList[1] + 'SmearedMass'].GetFunction('landau').SetLineColor(kBlack)
    h[partList[2] + 'SmearedMass'].Fit('landau')
    h[partList[2] + 'SmearedMass'].GetFunction('landau').SetLineColor(kBlack)
    if not partList[3]==None:
        h[partList[3] + 'SmearedMass'].Fit('landau')
        h[partList[3] + 'SmearedMass'].GetFunction('landau').SetLineColor(kBlack)

    h['SmearedMass'].Add(h[partList[1] + 'SmearedMass'])
    h['SmearedMass'].Add(h[partList[2] + 'SmearedMass'])
    if not partList[3]==None:
        h['SmearedMass'].Add(h[partList[3] + 'SmearedMass'])
    h['SmearedMass'].Draw('nostack')

    #h['DAUGHTERS_MOM'].Print('DaughterPProp'+ currentDate + '.png')

    if partList[3]==None:
        partString=''
    else:
        partString=' or pion'
    title='Probability Plots'
    h[key + '_PROB'] = TCanvas(key + '_PROB',title,1300,800)
    h[key + '_PROB'].Divide(3,2)
    cv = h[key + '_PROB'].cd(1)
    h[partList[1] + 'ProbMeasr'].SetMarkerColor(38)
    polyFit1.SetLineColor(4)
    h[partList[1] + 'ProbMeasr'].Fit('polyFit1')
    h[partList[1] + 'ProbMeasr'].Draw('E2')
    h[partList[1] + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
    h[partList[1] + 'ProbMeasr'].SetYTitle('Prob( particle = muon )')
    h[partList[1] + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    cv = h[key + '_PROB'].cd(2)
    h[partList[2] + 'ProbMeasr'].SetMarkerColor(46)
    polyFit2.SetLineColor(2)
    h[partList[2] + 'ProbMeasr'].Fit('polyFit2')
    h[partList[2] + 'ProbMeasr'].Draw('E2')
    h[partList[2] + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
    h[partList[2] + 'ProbMeasr'].SetYTitle('Prob( particle = kaon )')
    h[partList[2] + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    if not partList[3]==None:
        cv = h[key + '_PROB'].cd(3)
        h[partList[3] + 'ProbMeasr'].SetMarkerColor(8)
        polyFit2.SetLineColor(3)
        h[partList[3] + 'ProbMeasr'].Fit('polyFit2')
        h[partList[3] + 'ProbMeasr'].Draw('E2')
        h[partList[3] + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
        h[partList[3] + 'ProbMeasr'].SetYTitle('Prob( particle = pion )')
        h[partList[3] + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    cv = h[key + '_PROB'].cd(4)
    graph['partIDProb'] = TMultiGraph()
    graph['partIDProb'].SetName('Prob(correct ID paticle)')
    x1, y1 = array( 'd' ), array( 'd' )
    ex1, ey1 = array( 'd' ), array( 'd' )
    x2, y2 = array( 'd' ), array( 'd' )
    ex2, ey2 = array( 'd' ), array( 'd' )
    if not partList[3]==None:
        x3, y3 = array( 'd' ), array( 'd' )
        ex3, ey3 = array( 'd' ), array( 'd' )
    i=0
    n=0
    numBins = h[partList[1] + 'ProbMeasr'].GetNbinsX()
    for i in range(numBins):
        x1.append(h[partList[1] + 'ProbMeasr'].GetBinCenter(i))
        ex1.append(0)
        y1.append(h[partList[1] + 'ProbMeasr'].GetBinContent(i))
        ey1.append(ROOT.TMath.Sqrt((h[partList[1] + 'ProbMeasr'].GetBinContent(i))/10))
        x2.append(h[partList[2] + 'ProbMeasr'].GetBinCenter(i))
        ex2.append(0)
        y2.append(h[partList[2] + 'ProbMeasr'].GetBinContent(i))
        ey2.append(ROOT.TMath.Sqrt((h[partList[1] + 'ProbMeasr'].GetBinContent(i))/10))
        if not partList[3]==None:
            x3.append(h[partList[3] + 'ProbMeasr'].GetBinCenter(i))
            ex3.append(0)
            y3.append(h[partList[3] + 'ProbMeasr'].GetBinContent(i))
            ey3.append(ROOT.TMath.Sqrt((h[partList[1] + 'ProbMeasr'].GetBinContent(i))/10))
        n=n+1
    graph['1'] = TGraphErrors( n, x1, y1, ex1, ey1 )
    graph['1'].SetName('Prob(ID = ' + partList[1] + ')')
    graph['1'].SetTitle('Prob(ID = ' + partList[1] + ')')
    graph['1'].GetYaxis().SetTitle( 'Prob(particle = muon)' )
    graph['1'].SetLineColor( 4 )
    graph['1'].SetLineWidth( 1 )
    graph['1'].SetMarkerColor( 4 )
    graph['1'].SetMarkerStyle( 20 )
    graph['1'].SetMarkerSize(0.5)   
    graph['2'] = TGraphErrors( n, x2, y2, ex2, ey2 )
    graph['2'].SetName('Prob(ID = ' + partList[2] + ')')
    graph['2'].SetTitle('Prob(ID = ' + partList[2] + ')')
    graph['2'].GetYaxis().SetTitle( 'Prob(particle = kaon)' )
    graph['2'].SetLineColor( 2 )
    graph['2'].SetLineWidth( 1 )
    graph['2'].SetMarkerColor( 2 )
    graph['2'].SetMarkerStyle( 20 )
    graph['2'].SetMarkerSize(0.5)
    graph['partIDProb'].Add(graph['1'], 'PC')
    graph['partIDProb'].Add(graph['2'], 'PC')
    if not partList[3]==None:
        graph['3'] = TGraphErrors( n, x3, y3, ex3, ey3 )
        graph['3'].SetName('Prob(ID = ' + partList[3] + ')')
        graph['3'].SetTitle('Prob(ID = ' + partList[3] + ')') 
        graph['3'].GetYaxis().SetTitle( 'Prob(particle = pion)' )
        graph['3'].SetLineColor( 3 )
        graph['3'].SetLineWidth( 1 )
        graph['3'].SetMarkerColor( 3 )
        graph['3'].SetMarkerStyle( 20 )
        graph['3'].SetMarkerSize(0.5)
        graph['partIDProb'].Add(graph['3'], 'PC')
    graph['partIDProb'].Draw('A pfc plc')#P PLC PFCPLC PFC
    for ckey in graph:
        graph[ckey].GetXaxis().SetTitle( 'Mass [GeV/c2]' )
    graph['partIDProb'].GetYaxis().SetTitle( 'Prob(particle=(kaon or muon' + partString + '))' )
    graph['partIDProb'].GetYaxis().SetTitleOffset(1.5)
    graph['partIDProb'].GetXaxis().SetRangeUser(0,1.5)
    gPad.BuildLegend()
    #h[key + '_PROB'].Print('DaughterProb'+ currentDate + '.png')

def createHists_MuKa_exc():
    # Canvas 1
    ut.bookHist(h,'Kaon_recomass','Reconstructed Neutral Kaon Mass',200,0.4,0.6)
    ut.bookHist(h,'RPV_recomass','Reconstructed Neutralino Mass',300,0.5,1.5)
    ut.bookHist(h,'Kaon_recomom','Reconstructed',80,0,100)
    ut.bookHist(h,'Kaon_truemom','True',80,0,100)
    ut.bookHist(h,'Piplus_recomom','Reconstructed Pi+ Momentum',80,0,40)
    ut.bookHist(h,'Piminus_recomom','Reconstructed Pi- Momentum',80,0,40)
    ut.bookHist(h,'RPV_truemom','True',100,0,300)
    ut.bookHist(h,'RPV_recomom','Reconstructed',100,0,300)

    # Canvas 2
    ut.bookHist(h,'IP_target','Impact parameter to target',120,0,10)
    ut.bookHist(h,'ecalE','Energy deposited in ECAL',150,0,100)
    ut.bookHist(h,'doca','Distance of closest approach between tracks',150,0,3)
    ut.bookHist(h,'nmeas','No. of measurements in fitted tracks (ndf)',50,0,50)
    ut.bookHist(h,'Chi2','Fitted Tracks Reduced Chi Squared',150,0,3)
    ut.bookHist(h,'recovertex','Reconstructed neutralino decay vertex z-coordinate',100,-4000,4000)
def makePlots_MuKa_exc():
    ut.bookCanvas(h,key='Exc_RPV_N',title='Results 1',nx=1500,ny=800,cx=3,cy=2)
    cv = h['Exc_RPV_N'].cd(1)
    h['RPV_recomass'].SetXTitle('Invariant mass / [GeV/c2]')
    h['RPV_recomass'].SetYTitle('No. of particles')
    h['RPV_recomass'].SetLineColor(1)
    h['RPV_recomass'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_RPV_N'].cd(2)
    h['Kaon_recomass'].SetXTitle('Invariant mass / [GeV/c2]')
    h['Kaon_recomass'].SetYTitle('No. of particles')
    h['Kaon_recomass'].SetLineColor(1)
    h['Kaon_recomass'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_RPV_N'].cd(3)
    h['RPV_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['RPV_truemom'].SetYTitle('No. of particles')
    h['RPV_truemom'].SetLineColor(2)
    h['ths1'] = ROOT.THStack('RPVmom','True & Reconstructed Neutralino Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths1'].Add(h['RPV_truemom'])
    h['ths1'].Add(h['RPV_recomom'])
    h['ths1'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_RPV_N'].cd(4)
    h['Kaon_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['Kaon_truemom'].SetYTitle('No. of particles')
    h['Kaon_truemom'].SetLineColor(2)
    h['ths2'] = ROOT.THStack('Kaonmom','True & Reconstructed K0 Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths2'].Add(h['Kaon_truemom'])
    h['ths2'].Add(h['Kaon_recomom'])
    h['ths2'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_RPV_N'].cd(5)
    h['Piplus_recomom'].SetXTitle('Momentum / [GeV/c]')
    h['Piplus_recomom'].SetYTitle('No. of particles')
    h['Piplus_recomom'].SetLineColor(3)
    h['Piplus_recomom'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_RPV_N'].cd(6)
    h['Piminus_recomom'].SetXTitle('Momentum / [GeV/c]')
    h['Piminus_recomom'].SetYTitle('No. of particles')
    h['Piminus_recomom'].SetLineColor(3)
    h['Piminus_recomom'].Draw()
    h['Exc_RPV_N'].Print('Exc_RPV_N.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='Exc_Vetos',title='Results 2',nx=1500,ny=800,cx=3,cy=2)
    cv = h['Exc_Vetos'].cd(1)
    h['IP_target'].SetXTitle('Impact Parameter / [cm]')
    h['IP_target'].SetYTitle('Frequency')
    h['IP_target'].SetLineColor(1)
    h['IP_target'].SetFillColor(17)
    h['IP_target'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_Vetos'].cd(2)
    h['ecalE'].SetXTitle('Energy / [GeV/c2]')
    h['ecalE'].SetYTitle('Frequency')
    h['ecalE'].SetLineColor(1)
    h['ecalE'].SetFillColor(17)
    h['ecalE'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_Vetos'].cd(3)
    h['doca'].SetXTitle('Distance / [cm]')
    h['doca'].SetYTitle('Frequency')
    h['doca'].SetLineColor(1)
    h['doca'].SetFillColor(17)
    h['doca'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_Vetos'].cd(4)
    h['nmeas'].SetXTitle('ndf')
    h['nmeas'].SetYTitle('No. of tracks')
    h['nmeas'].SetLineColor(1)
    h['nmeas'].SetFillColor(17)
    h['nmeas'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_Vetos'].cd(5)
    h['Chi2'].SetXTitle('Reduced Chi Squared')
    h['Chi2'].SetYTitle('Frequency')
    h['Chi2'].SetLineColor(1)
    h['Chi2'].SetFillColor(17)
    h['Chi2'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Exc_Vetos'].cd(6)
    h['recovertex'].SetXTitle('Z / [cm]')
    h['recovertex'].SetYTitle('Frequency')
    h['recovertex'].SetLineColor(1)
    h['recovertex'].SetFillColor(17)
    h['recovertex'].Draw()
    h['Exc_Vetos'].Print('Exc_Vetos.png')

def createHists_DarkPhot():
    # Canvas 1
    ut.bookHist(h,'DP_recomass','Reconstructed Dark Photon Mass',500,0,0.4)   # reconstructed mass
    ut.bookHist(h,'DP_truemom','True Dark Photon Momentum',100,0.,300.)   # true momentum distribution
    ut.bookHist(h,'DP_recomom','Reconstructed Dark Photon Momentum',100,0.,300)   # reconstructed momentum distribution
    ut.bookHist(h,'DP_mom_diff','True/Reco Dark Photon Momentum Difference',100,-3.,3)   # true/reco momentum difference
    ut.bookHist(h,'DP_beta','Reconstructed Dark Photon Beta in Z-direction',100,0.994,1)
    ut.bookHist(h,'DP_gamma','Reconstructed Dark Photon Gamma in Z-direction',100,0,200)
    ut.bookHist(h,'DP_theta','Angle between Dark Photon Momentum and Beam Line',100,0,50)

    # Canvas 2
    ut.bookHist(h,'eplus_truemom','True e+ Momentum',100,0.,140.)   # RPV muon daughter reco momentum
    ut.bookHist(h,'eminus_truemom','True e- Momentum',100,0.,140.)   # RPV pion daughter reco momentum
    ut.bookHist(h,'eplus_recomom','Reconstructed e+ Momentum',100,0.,140.)   # RPV muon daughter true momentum
    ut.bookHist(h,'eminus_recomom','Reconstructed e- Momentum',100,0.,140.)   # RPV pion daughter true momentum

    # Canvas 3
    ut.bookHist(h,'IP_target','Impact parameter to target',120,0,10)
    ut.bookHist(h,'ecalE','Energy deposited in ECAL',150,0,100)
    ut.bookHist(h,'doca','Distance of closest approach between tracks',150,0,3)
    ut.bookHist(h,'nmeas','No. of measurements in fitted tracks (ndf)',50,0,50)
    ut.bookHist(h,'Chi2','Fitted Tracks Reduced Chi Squared',150,0,3)
    ut.bookHist(h,'recovertex','Reconstructed neutralino decay vertex z-coordinate',100,-4000,4000)
def makePlots_DarkPhot():
    ut.bookCanvas(h,key='DP',title='Results 1',nx=1500,ny=800,cx=3,cy=2)
    cv = h['DP'].cd(1)
    h['DP_recomass'].SetXTitle('Invariant mass / [GeV/c2]')
    h['DP_recomass'].SetYTitle('No. of Particles')
    h['DP_recomass'].SetLineColor(1)
    h['DP_recomass'].Draw()
    print('\nDark Photon mass Gaussian fit:\n')
    fitSingleGauss('DP_recomass',0.1,0.3)
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP'].cd(2)
    h['DP_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['DP_truemom'].SetYTitle('No. of Particles')
    h['DP_truemom'].SetLineColor(2)
    h['ths1'] = ROOT.THStack('DPmom','True & Reconstructed Dark Photon Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths1'].Add(h['DP_truemom'])
    h['ths1'].Add(h['DP_recomom'])
    h['ths1'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP'].cd(3)
    h['DP_mom_diff'].SetXTitle('Momentum Difference / [GeV/c]')
    h['DP_mom_diff'].SetYTitle('No. of Particles')
    h['DP_mom_diff'].SetLineColor(1)    
    h['DP_mom_diff'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP'].cd(4)
    h['DP_beta'].SetXTitle('Beta')
    h['DP_beta'].SetYTitle('No. of Particles')
    h['DP_beta'].SetLineColor(1)
    h['DP_beta'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP'].cd(5)
    h['DP_gamma'].SetXTitle('Gamma')
    h['DP_gamma'].SetYTitle('No. of Particles')
    h['DP_gamma'].SetLineColor(1)
    h['DP_gamma'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP'].cd(6)
    h['DP_theta'].SetXTitle('Theta / [mrad]')
    h['DP_theta'].SetYTitle('No. of Particles')
    h['DP_theta'].SetLineColor(1)
    h['DP_theta'].Draw()
    h['DP'].Print('DarkPhot.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='ee',title='Results 2',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['ee'].cd(1)
    h['eplus_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['eplus_truemom'].SetYTitle('No. of Particles')
    h['eplus_truemom'].SetLineColor(2)
    h['ths2'] = ROOT.THStack('Kamom','e+ True & Reconstructed Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths2'].Add(h['eplus_truemom'])
    h['ths2'].Add(h['eplus_recomom'])
    h['ths2'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['ee'].cd(2)
    h['eminus_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['eminus_truemom'].SetYTitle('No. of Particles')
    h['eminus_truemom'].SetLineColor(2)
    h['ths3'] = ROOT.THStack('Mumom','e- True & Reconstructed Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths3'].Add(h['eminus_truemom'])
    h['ths3'].Add(h['eminus_recomom'])
    h['ths3'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    h['ee'].Print('Electrons.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='DP_Vetos',title='Results 3',nx=1500,ny=800,cx=3,cy=2)
    cv = h['DP_Vetos'].cd(1)
    h['IP_target'].SetXTitle('Impact Parameter / [cm]')
    h['IP_target'].SetYTitle('Frequency')
    h['IP_target'].SetLineColor(1)
    h['IP_target'].SetFillColor(17)
    h['IP_target'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP_Vetos'].cd(2)
    h['ecalE'].SetXTitle('Energy / [GeV/c2]')
    h['ecalE'].SetYTitle('Frequency')
    h['ecalE'].SetLineColor(1)
    h['ecalE'].SetFillColor(17)
    h['ecalE'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP_Vetos'].cd(3)
    h['doca'].SetXTitle('Distance / [cm]')
    h['doca'].SetYTitle('Frequency')
    h['doca'].SetLineColor(1)
    h['doca'].SetFillColor(17)
    h['doca'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP_Vetos'].cd(4)
    h['nmeas'].SetXTitle('ndf')
    h['nmeas'].SetYTitle('No. of tracks')
    h['nmeas'].SetLineColor(1)
    h['nmeas'].SetFillColor(17)
    h['nmeas'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP_Vetos'].cd(5)
    h['Chi2'].SetXTitle('Reduced Chi Squared')
    h['Chi2'].SetYTitle('Frequency')
    h['Chi2'].SetLineColor(1)
    h['Chi2'].SetFillColor(17)
    h['Chi2'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['DP_Vetos'].cd(6)
    h['recovertex'].SetXTitle('Z / [cm]')
    h['recovertex'].SetYTitle('Frequency')
    h['recovertex'].SetLineColor(1)
    h['recovertex'].SetFillColor(17)
    h['recovertex'].Draw()
    h['DP_Vetos'].Print('DP_Vetos.png')

######################################
#########  CHECKS FUNCTIONS  #########
def dist2InnerWall(X,Y,Z):
  dist = 0
 # return distance to inner wall perpendicular to z-axis, if outside decayVolume return 0.
  node = sGeo.FindNode(X,Y,Z)
  if ShipGeo.tankDesign < 5:
     if not 'cave' in node.GetName(): return dist  # TP 
  else:
     if not 'decayVol' in node.GetName(): return dist
  start = array('d',[X,Y,Z])
  nsteps = 8
  dalpha = 2*ROOT.TMath.Pi()/nsteps
  rsq = X**2+Y**2
  minDistance = 100 *u.m
  for n in range(nsteps):
    alpha = n * dalpha
    sdir  = array('d',[ROOT.TMath.Sin(alpha),ROOT.TMath.Cos(alpha),0.])
    node = sGeo.InitTrack(start, sdir)
    nxt = sGeo.FindNextBoundary()
    if ShipGeo.tankDesign < 5 and nxt.GetName().find('I')<0: return 0    
    distance = sGeo.GetStep()
    if distance < minDistance  : minDistance = distance
  return minDistance

def isInFiducial(X,Y,Z):
   if Z > ShipGeo.TrackStation1.z : return False
   if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
   # typical x,y Vx resolution for exclusive Neutralino decays 0.3cm,0.15cm (gaussian width)
   if dist2InnerWall(X,Y,Z)<5*u.cm: return False
   return True 

def checkFiducialVolume(sTree,trackkey,dy):
    # extrapolate track to middle of magnet and check if in decay volume
    inside = True
    if not fiducialCut: return True
    fT = sTree.FitTracks[trackkey]
    rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,ShipGeo.Bfield.z)
    if not rc: return False
    if not dist2InnerWall(pos.X(),pos.Y(),pos.Z())>0: return False
    return inside

def myVertex(t1,t2,PosDir):
 # closest distance between two tracks
    # d = |pq . u x v|/|u x v|
   a = ROOT.TVector3(PosDir[t1][0](0) ,PosDir[t1][0](1), PosDir[t1][0](2))
   u = ROOT.TVector3(PosDir[t1][1](0),PosDir[t1][1](1),PosDir[t1][1](2))
   c = ROOT.TVector3(PosDir[t2][0](0) ,PosDir[t2][0](1), PosDir[t2][0](2))
   v = ROOT.TVector3(PosDir[t2][1](0),PosDir[t2][1](1),PosDir[t2][1](2))
   pq = a-c
   uCrossv = u.Cross(v)
   dist  = pq.Dot(uCrossv)/(uCrossv.Mag()+1E-8)
   # u.a - u.c + s*|u|**2 - u.v*t    = 0
   # v.a - v.c + s*v.u    - t*|v|**2 = 0
   E = u.Dot(a) - u.Dot(c) 
   F = v.Dot(a) - v.Dot(c) 
   A,B = u.Mag2(), -u.Dot(v) 
   C,D = u.Dot(v), -v.Mag2()
   t = -(C*E-A*F)/(B*C-A*D)
   X = c.x()+v.x()*t
   Y = c.y()+v.y()*t
   Z = c.z()+v.z()*t
   return X,Y,Z,abs(dist)

def ImpactParameter(point,tPos,tMom):
  t = 0
  if hasattr(tMom,'P'): P = tMom.P()
  else:                 P = tMom.Mag()
  for i in range(3):   t += tMom(i)/P*(point(i)-tPos(i)) 
  dist = 0
  for i in range(3):   dist += (point(i)-tPos(i)-t*tMom(i)/P)**2
  dist = ROOT.TMath.Sqrt(dist)
  return dist

def fitSingleGauss(x,ba=None,be=None):
    name    = 'myGauss_'+x 
    myGauss = h[x].GetListOfFunctions().FindObject(name)
    if not myGauss:
       if not ba : ba = h[x].GetBinCenter(1) 
       if not be : be = h[x].GetBinCenter(h[x].GetNbinsX()) 
       bw    = h[x].GetBinWidth(1) 
       mean  = h[x].GetMean()
       sigma = h[x].GetRMS()
       norm  = h[x].GetEntries()*0.3
       myGauss = ROOT.TF1(name,'[0]*'+str(bw)+'/([2]*sqrt(2*pi))*exp(-0.5*((x-[1])/[2])**2)+[3]',4)
       myGauss.SetParameter(0,norm)
       myGauss.SetParameter(1,mean)
       myGauss.SetParameter(2,sigma)
       myGauss.SetParameter(3,1.)
       myGauss.SetParName(0,'Signal')
       myGauss.SetParName(1,'Mean')
       myGauss.SetParName(2,'Sigma')
       myGauss.SetParName(3,'bckgr')
    h[x].Fit(myGauss,'','',ba,be) 

def RedoVertexing(t1,t2):    
     PosDir = {} 
     for tr in [t1,t2]:
         xx  = sTree.FitTracks[tr].getFittedState()     #sets track info to xx
         PosDir[tr] = [xx.getPos(),xx.getDir()]         #position and direction of track
     xv,yv,zv,doca = myVertex(t1,t2,PosDir)
     # as we have learned, need iterative procedure
     dz = 99999.
     reps,states,newPosDir = {},{},{}
     parallelToZ = ROOT.TVector3(0., 0., 1.)
     rc = True 
     step = 0
     while dz > 0.1:
         zBefore = zv
         newPos = ROOT.TVector3(xv,yv,zv)
         # make a new rep for track 1,2
         for tr in [t1,t2]:
             xx = sTree.FitTracks[tr].getFittedState()
             reps[tr]   = ROOT.genfit.RKTrackRep(xx.getPDG())
             states[tr] = ROOT.genfit.StateOnPlane(reps[tr])
             reps[tr].setPosMom(states[tr],xx.getPos(),xx.getMom())
             try:
                 reps[tr].extrapolateToPoint(states[tr], newPos, False)
             except:
                 #print ('SHiPAna: extrapolation did not work (@RedoVertexing).')
                 rc = False
                 break
             newPosDir[tr] = [reps[tr].getPos(states[tr]),reps[tr].getDir(states[tr])]
         if not rc: break
         xv,yv,zv,doca = myVertex(t1,t2,newPosDir)
         dz = abs(zBefore-zv)
         step+=1
         if step > 10:
             #print ('Abort iteration, too many steps, pos=',xv,yv,zv,' doca=',doca,'z before and dz',zBefore,dz)
             rc = False
             break 
     if not rc: return -1,-1,-1,xv,yv,zv,doca # extrapolation failed, makes no sense to continue
     LV={}
     for tr in [t1,t2]: # from here on we have reproduced (see inv_mass()) 
         mom = reps[tr].getMom(states[tr])
         pid = abs(states[tr].getPDG()) 
         if pid == 2212: pid = 211 #why
         mass = PDG.GetParticle(pid).Mass()
         E = ROOT.TMath.Sqrt( mass*mass + mom.Mag2() )
         LV[tr] = ROOT.TLorentzVector()
         LV[tr].SetPxPyPzE(mom.x(),mom.y(),mom.z(),E)
     NeutralinoLV = LV[t1]+LV[t2]

     return NeutralinoLV,LV[t1],LV[t2],xv,yv,zv,doca

def SingleTrack_4Mom(t1):
    fittedtrack  = sTree.FitTracks[t1].getFittedState()
    Partkey = sTree.fitTrack2MC[t1]   # matches track to MC particle key
    true_part = sTree.MCTrack[Partkey]   # gives MC particle data 
         
    M = true_part.GetMass()   # mass 
    P = fittedtrack.getMomMag()   # momentum in x
    Px = fittedtrack.getMom().x()   # momentum in y
    Py = fittedtrack.getMom().y()   # momentum in z
    Pz = fittedtrack.getMom().z()   # momentum magnitude
    E = ROOT.TMath.Sqrt((M**2) + (P**2))   # energy
                                                            
    LV = ROOT.TLorentzVector()   # declares variable as TLorentzVector class
    LV.SetPxPyPzE(Px,Py,Pz,E)   # inputs 4-vector elements

    return LV

######################################
##########  TOOLS FUNCTIONS ##########
def ecalMinIon(partkey):
    ecalE_tot = 0
    if sTree.GetBranch('EcalPoint'):
        ecal_Etot = 0
        for hits in sTree.EcalPoint:
            ecal_TrackID = hits.GetTrackID()
            if ecal_TrackID == partkey:
                ecalE = hits.GetEnergyLoss()
                ecalE_tot += ecalE
    return ecalE_tot

def muonstationHits(partkey):
    if sTree.GetBranch('muonPoint'):
        true_part = sTree.MCTrack[partkey]
        if abs(true_part.GetPdgCode()) == 13:
            muonstation_z = []
            for hits in sTree.muonPoint:
                muonTrackID = hits.GetTrackID()
                if muonTrackID == partkey:
                    muonstation_z.append(hits.GetZ())
            if len(muonstation_z) > 1 and muonstation_z[0] == 4017.0: return True
            else: return False
        else: return False
    else: return False

def time_res(partMom,partName,partkey,eventN,succEventM):
    straw_time = None
    ecal_time = None
    deltaT = None
    true_part = sTree.MCTrack[partkey]   # finds MC particle
    pdg = true_part.GetPdgCode()   # identifies particles from PDG code
    c = 2.99792458*(10**8)  #speed of light

    if sTree.GetBranch('strawtubesPoint'):
        x_array = []
        y_array = []
        z_array = []
        r_array = []
        t_array = []
        px_array = []
        py_array = []
        pz_array = []
        straw_time = 0
        for k,hits in enumerate(sTree.strawtubesPoint):
            straw_TrackID = hits.GetTrackID()
            if straw_TrackID == partkey:
                x_array.append(hits.GetX())
                y_array.append(hits.GetY())
                z_array.append(hits.GetZ())
                px_array.append(hits.GetPx())
                py_array.append(hits.GetPy())
                pz_array.append(hits.GetPz())
                t_array.append(hits.GetTime())

        N = len(z_array)
        R1 = 0
        for j in range(0,N-1):
            s = ROOT.TMath.Sqrt(((0.01*x_array[j] - 0.01*x_array[j+1])**2) + ((0.01*y_array[j] - 0.01*y_array[j+1])**2) + ((0.01*z_array[j] - 0.01*z_array[j+1])**2))
            r_array.append(s)
            R1 = sum(r_array)   # total distance travelled in the straw tubes
        
        min_z_index = z_array.index(min(z_array))
        inmostStrawZ = 0.01*z_array[min_z_index]
        inmostStrawX = 0.01*x_array[min_z_index]
        inmostStrawY = 0.01*y_array[min_z_index]
        straw_time = t_array[min_z_index]
        if straw_time != None:
            smearStrawTime = np.random.normal(loc=straw_time,scale=0.01,size=None) # current width of 10 ps 

        strawPx = px_array[min_z_index]
        strawPy = py_array[min_z_index]
        strawPz = pz_array[min_z_index]
        strawP = ROOT.TMath.Sqrt((strawPx**2) + (strawPy**2) + (strawPz**2)) 

        num_hits = len(z_array)   # number of elements in the list
        if abs(pdg) == 13:   # muon
            h[partName + 'StrawHits'].Fill(num_hits)
            for hit in pz_array:
                if eventN == succEventM:
                    h[partName + 'StrawHitsMom'].Fill(hit)   # muon z-momentum through straw tubes for particular event
        if abs(pdg) == 321:   # kaon
            h[partName + 'StrawHits'].Fill(num_hits)
            for hit in pz_array:
                if eventN == succEventM:
                    h[partName + 'StrawHitsMom'].Fill(hit)   # kaon z-momentum through straw tubes for particular event

        if sTree.GetBranch('EcalPoint'):
            ecal_time = 0
            if not straw_time <= 0:
                for k,hits in enumerate(sTree.EcalPoint):
                    ecal_TrackID = hits.GetTrackID()
                    if ecal_TrackID == partkey:
                        ecal_x = 0.01*hits.GetX()
                        ecal_y = 0.01*hits.GetY()
                        ecal_z = 0.01*hits.GetZ()
                        ecal_time = hits.GetTime()
                        ecal_px = hits.GetPx()
                        ecal_py = hits.GetPy()
                        ecal_pz = hits.GetPz()
                        ecalP = ROOT.TMath.Sqrt((ecal_px**2) + (ecal_py**2) + (ecal_pz**2))

                        if ecal_time != None:
                            smearEcalTime = np.random.normal(loc=ecal_time,scale=0.01,size=None) # current width of 10 ps
                            
                            if not ecal_time <= straw_time:
                                smearDeltaT = abs(smearStrawTime - smearEcalTime)
                                deltaT = abs(straw_time - ecal_time)
                                deltaP = strawP - ecalP
                                r = ROOT.TMath.Sqrt(((ecal_x - inmostStrawX)**2) + ((ecal_y - inmostStrawY)**2) + ((ecal_z - inmostStrawZ)**2))
                                max_z_index = z_array.index(max(z_array))
                                laststraw_x = 0.01*x_array[max_z_index]
                                laststraw_y = 0.01*y_array[max_z_index]
                                laststraw_z = 0.01*z_array[max_z_index]
                                R2 = ROOT.TMath.Sqrt(((ecal_x - laststraw_x)**2) + ((ecal_y - laststraw_y)**2) + ((ecal_z - laststraw_z)**2))
                                R = R1+R2                           # better approximation of distance
                                rdiff = abs(R-r)
                                smearV = ((r/smearDeltaT)*(10**9) ) # velocity in flight
                                v = ((R/deltaT)*(10**9) )           # units of nanoseconds

    if deltaT != None:
        h[partName + 'StrawTime'].Fill(smearStrawTime)
        h[partName + 'EcalTime'].Fill(smearEcalTime)
        h[partName + 'DirDeltaTimeSmeared'].Fill(smearDeltaT)
        h[partName + 'DirDeltaTime'].Fill(deltaT)
        h[partName + 'FlightLen'].Fill(r)
        h[partName + 'FlightLenImproved'].Fill(R)
        h[partName + 'FlightLenDelta'].Fill(rdiff)
        smearBeta = smearV/c
        h[partName + 'SpeedSmeared'].Fill(smearBeta)
        beta = v/c
        h[partName + 'Speed'].Fill(beta)
        h[partName + 'StrawMom'].Fill(strawP)
        h[partName + 'EcalMom'].Fill(ecalP)
        Dp = strawP-ecalP
        h[partName + 'DeltaMom'].Fill(Dp)
        if beta < 1:
            partMass = partMom*(ROOT.TMath.Sqrt(1-(smearBeta**2)))/smearBeta
            h[partName + 'MassT'].Fill(partMass)
        if smearBeta < 1:
            smearedM = partMom*(ROOT.TMath.Sqrt(1-(smearBeta**2)))/smearBeta
            h[partName + 'SmearedMass'].Fill(smearedM)
            h['TotalSmearedMass'].Fill(smearedM)  

def track_checks(index,veto,fill):
    check = 0
    
    partkey = sTree.fitTrack2MC[index]  
    true_part = sTree.MCTrack[partkey]
    reco_part = sTree.FitTracks[index]

    fit_status = reco_part.getFitStatus()
    fit_nmeas = fit_status.getNdf()
    fit_rchi2 = fit_status.getChi2()                      
    fit_chi2 = (fit_rchi2/fit_nmeas)

    #if fit_status.isFitConverged():
       
    if fill == 1: h['Chi2'].Fill(fit_chi2)
    if not fit_chi2 < chi2Cut:
        #print('Chi squared value too high')
        if check == 0: veto[1] += 1
        check = -1

    if fill == 1: h['nmeas'].Fill(fit_nmeas)
    if not fit_nmeas > measCut:
        #print('Too few measurements')
        if check == 0: veto[2] += 1
        check = -1    
        
    if not checkFiducialVolume(sTree,index,dy): 
        #print('Track outside fiducial volume')
        if check == 0: veto[4] += 1
        check = -1

    if sTree.GetBranch('EcalPoint'):
        ecal_Etot = 0
        for hits in sTree.EcalPoint:
            ecal_TrackID = hits.GetTrackID()
            if ecal_TrackID == partkey:
                ecalE = hits.GetEnergyLoss()
                ecal_Etot += ecalE
        if fill == 1: h['ecalE'].Fill(ecal_Etot)
        if not ecal_Etot > ecalCut:
            #print('Not enough energy deposited in the ECAL')
            if check == 0: veto[5] += 1
            check = -1

    if sTree.GetBranch('muonPoint'):
        if abs(true_part.GetPdgCode()) == 13:
            muonstation_z = []
            for hits in sTree.muonPoint:
                muonTrackID = hits.GetTrackID()
                if muonTrackID == partkey:
                    muonstation_z.append(hits.GetZ())
            if not len(muonstation_z) > 2:
                #print('Did not leave hits in both the 1st and 2nd muon stations')
                if check == 0: veto[6] += 1
                check = -1
            else:
                if not muonstation_z[0] == 4017.0: #first muon station position (cm)
                    if check == 0: veto[6] += 1
                    check =-1

    return check,veto

def createRatio(h1, h2, histname):
    h3 = h1.Clone(histname)
    h3.SetMarkerStyle(20)
    h3.SetMarkerSize(0.7)
    h3.SetTitle('')
    h3.SetMinimum(0.8)
    h3.SetMaximum(1.35)
    # Set up plot for markers and errors
    h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h1,h2,1,1,'B')
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetRangeUser(-0.1,1.2)
    y.SetTitleOffset(1.)
    ## Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetRangeUser(0,1.5)
    return h3

def getRPVBranchRatio(stEntry):
    rpvsusy_instance = rpvsusy.RPVSUSY(1.,[0.2,0.03],1e3,1,True)
    bRatio = rpvsusy_instance.findDecayBranchingRatio(stEntry)
    return bRatio

def signalAcceptance(bRatio,recEntry,simEntry):
    accp=bRatio*(float(recEntry)/simEntry)
    return accp

def print_menu(): 
    print ('\n \n' + 30 * '-' + 'MENU' + 30 * '-')
    print ('1. RPV SUSY Benchmark1 --> K+/-  mu+/- visible final state')
    print ('2. RPV SUSY Benchmark1 --> K*+/- mu+/- visible final state')
    print ('3. Dark Photon         --> l+/-  l+/-  visible final state')
    print ('4. Exit')
    print (64 * '-'+ '\n') 

nEvents = min(sTree.GetEntries(),nEvents)

########################################
#########  ANALYSIS FUNCTIONS  #########
def finStateMuKa():
    if sTree.GetBranch('FitTracks'):
        particleList=['Neutralino','Mu','K+/-',None]
        create_Hists(particleList)
        successful_events = []
        muveto=10*[0]
        k2mu_MotherHP = 0
        NeutralinoMuEvents=0
        finStateEvents=0
        HP2ka_veto = 10*[0]
        veto = 10*[0]
        acceptance = 10*[0.]
        
        for n in range(nEvents):    # loops over events
            rc = sTree.GetEntry(n)    # load tree entry
            event = True

            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles                                   
                fitstatus1 = reco_part.getFitStatus()
                if not fitstatus1.isFitConverged(): continue
                muPartkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_Mu = sTree.MCTrack[muPartkey]               # gives MC particle data
                wgMu = true_Mu.GetWeight()
                if not wgMu>0.: wgMu=1.

                if abs(true_Mu.GetPdgCode()) == 13:        # checks particle is muon
                    MuMotherKey = true_Mu.GetMotherId()             # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[MuMotherKey]             # obtains mother particle data

                    #check, muveto = track_checks(index,muveto,0)   # this reduces the entries in the table
                    #if not check == 0:   # performs various checks (i.e. vertex position, fiducial volume,...)
                    #    continue

                    if true_mother.GetPdgCode() == 321:     #checks mother is kaon
                            muGrannyKey = true_mother.GetMotherId()
                            true_gran = sTree.MCTrack[muGrannyKey]
                            if true_gran.GetPdgCode() == 9900015:
                                k2mu_MotherHP+=1 # kaon has decayed to a muon in flight
                                check,HP2ka_veto = track_checks(index,HP2ka_veto,0)
                                if check == -1: HP2ka_veto[0] += 1 

                    if true_mother.GetPdgCode() == 9900015:              # checks mother is Neutralino
                        NeutralinoMuEvents+=1
                        
                        for index2,reco_part2 in enumerate(sTree.FitTracks):  # loops over index and data of track particles
                            fit_status2 = reco_part2.getFitStatus()
                            p2Partkey = sTree.fitTrack2MC[index2]                 # matches track to MC particle key
                            true_part2 = sTree.MCTrack[p2Partkey]                 # gives MC particle data
                            wgKaon = true_part2.GetWeight()
                            if not wgKaon>0.: wgKaon=1.

                            if abs(true_part2.GetPdgCode()) == 321:               # checks particle is kaon
                                part2MotherKey = true_part2.GetMotherId()            # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[part2MotherKey]          # obtains mother particle data
                                    
                                if (part2MotherKey==MuMotherKey and true_mother.GetPdgCode() == 9900015):#and daught2=='K+/-') or (true_mother.GetPdgCode() == daught2_PDG and daught2!='K+/-'):                 # check if keys are the same
                                    finStateEvents+=1

                                    ####################
                                    #####  CHECKS  #####                                        
                                    if fit_status2.isFitConverged():
                                        veto[0] += 1
                                    else: continue

                                    Neutralino_LVec,Mu_LVec,part2_LVec,X,Y,Z,doca = RedoVertexing(index,index2) # uses RedoVertexing to iterate track fitting
                                    if not Neutralino_LVec == -1: veto[9] += 1
                                    else: 
                                        print('RedoVertexing extrapolation failed (event ' + str(n) + ')')
                                        continue

                                    nmeas_muon = fitstatus1.getNdf()
                                    chi2_muon = fitstatus1.getChi2()
                                    rchi2_muon = (chi2_muon/nmeas_muon)
                                    nmeas_kaon = fit_status2.getNdf()
                                    chi2_kaon = fit_status2.getChi2()
                                    rchi2_kaon = (chi2_kaon/nmeas_kaon)
                                    h['Chi2'].Fill(rchi2_muon)
                                    h['Chi2'].Fill(rchi2_kaon)
                                    if rchi2_muon < chi2Cut and rchi2_kaon < chi2Cut:
                                        veto[1] += 1
                                    else: event = False

                                    h['nmeas'].Fill(nmeas_muon)
                                    h['nmeas'].Fill(nmeas_kaon)
                                    if event == True:
                                        if nmeas_muon > measCut and nmeas_kaon > measCut:
                                            veto[2] += 1
                                        else: event = False

                                    h['recovertex'].Fill(Z)
                                    if event == True:
                                        if isInFiducial(X,Y,Z):
                                            veto[3] += 1
                                        else: event = False

                                    if event == True:
                                        if checkFiducialVolume(sTree,index,dy) and checkFiducialVolume(sTree,index2,dy): 
                                            veto[4] += 1
                                        else: event = False

                                    ecalE_muon = ecalMinIon(muPartkey)
                                    ecalE_kaon = ecalMinIon(p2Partkey)
                                    h['ecalE'].Fill(ecalE_muon)
                                    h['ecalE'].Fill(ecalE_kaon)
                                    if event == True:
                                        if ecalE_muon > ecalCut and ecalE_kaon > ecalCut:
                                            veto[5] += 1
                                        else: event = False

                                    if event == True:
                                        if muonstationHits(muPartkey):
                                            veto[6] += 1
                                        else: event = False
                                    
                                    h['doca'].Fill(doca)
                                    if event == True:
                                        if doca < docaCut: 
                                            veto[7] +=1
                                        else: event = False

                                        
                                    NProduction_X = true_mother.GetStartX()
                                    NProduction_Y = true_mother.GetStartY()
                                    NProduction_Z = true_mother.GetStartZ()
                                    NProduction_T = true_mother.GetStartT()
                                    Neutralino_Pos = ROOT.TLorentzVector()
                                    Neutralino_Pos.SetXYZT(NProduction_X,NProduction_Y,NProduction_Z,NProduction_T)
                                    tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
                                    ip = ImpactParameter(tr,Neutralino_Pos,Neutralino_LVec)
                                    h['IP_target'].Fill(ip)
                                    if event == True:
                                        if ip < ipCut:
                                            veto[8] += 1
                                        else: event = False

                                    if event == False: continue


                                    #######################
                                    ### Reco Properties ###
                                    Neutralino_mass = Neutralino_LVec.M()
                                    Neutralino_recoMom = Neutralino_LVec.P()
                                    p2M = part2_LVec.M()
                                    p2P = part2_LVec.P()
                                    p2E = part2_LVec.E()
                                    muM = Mu_LVec.M()
                                    muP = Mu_LVec.P()
                                    muE = Mu_LVec.E()
                                    #######################
                                    ### True Properties ###
                                    p2MotherTrueMass = true_mother.GetMass()
                                    p2MotherTrueMom = true_mother.GetP()             
                                    part2TrueMom = true_part2.GetP()
                                    part2TrueMass = true_part2.GetMass()
                                    MuTrueMom = true_Mu.GetP()
                                    MuTrueMass = true_Mu.GetMass()

                                    mom_diff = p2MotherTrueMom - Neutralino_recoMom

                                    part2='K+/-'
                                    h['Mu' + 'TrueMom'].Fill(MuTrueMom)
                                    h['Mu' + 'TrueMass'].Fill(MuTrueMass)
                                    h['K+/-' + 'TrueMom'].Fill(part2TrueMom)
                                    h['K+/-' + 'TrueMass'].Fill(part2TrueMass)
                                    h['Neutralino' + 'TrueMass'].Fill(p2MotherTrueMass)            
                                    h['Neutralino' + 'TrueMom'].Fill(p2MotherTrueMom)
                                    h['Neutralino' + 'RecoMass'].Fill(Neutralino_mass)                        
                                    h['Neutralino' + 'RecoMom'].Fill(Neutralino_recoMom)                            
                                    h['Neutralino' + 'DeltaMom'].Fill(mom_diff)
                                    h['K+/-' + 'RecoMom'].Fill(p2P)
                                    h['Mu' + 'RecoMom'].Fill(muP)
                                    h['K+/-' + 'RecoMass'].Fill(p2M)
                                    h['Mu' + 'RecoMass'].Fill(muM)

                                    Neutralino_Zmom = Neutralino_LVec.Pz()
                                    Neutralino_Zbeta = 1/ROOT.TMath.Sqrt(1 + ((Neutralino_mass/Neutralino_Zmom)**2))
                                    Neutralino_Zgamma = 1/(ROOT.TMath.Sqrt(1 - (Neutralino_Zbeta**2)))
                                    cos_theta = Neutralino_Zmom/Neutralino_recoMom
                                    theta = 1000*(ROOT.TMath.ACos(cos_theta))   # angle between beam line and neutralino momentum (mrad)

                                    h['Neutralino' + 'Beta'].Fill(Neutralino_Zbeta)
                                    h['Neutralino' + 'Gamma'].Fill(Neutralino_Zgamma)
                                    h['Neutralino' + 'Theta'].Fill(theta)

                                    successful_events.append(n)   # adds entries to the list
                                    m = successful_events[0]   # arbitrarily picks the first one as an example

                                    time_res(muP,'Mu',muPartkey,n,m)
                                    time_res(p2P,'K+/-',p2Partkey, n, m)

    brRatio = getRPVBranchRatio('N -> K+ mu-')
    print(brRatio)
    #accp = signalAcceptance('N -> K+ mu-')
    for eRemn in range(9):
        if eRemn == 0:
            acceptance[eRemn] = signalAcceptance(brRatio, veto[eRemn],n+1)
        else:
            acceptance[eRemn] = signalAcceptance(brRatio, veto[eRemn],n+1)
            print(acceptance[eRemn])

                                        
    print(2*' ' + 80*'_')
    accepted = len(successful_events)
    rejected = veto[0] - accepted
    print('\n\t' + str(n+1) +  ' total number of events')
    print('\t' + str(NeutralinoMuEvents) + '  ' + 'Mu' + ' number of events (with ' + 'Neutralino' + ' mother)')
    print('\t' + str(finStateEvents) + ' number of events produced final state ')
    print('\n\t' + str(veto[0]) + ' events reconstructed for this decay mode')
    print('\t' + str(veto[9]) + ' events with successful RedoVertexing extrapolations')
    print('\t' + str(accepted) + ' events not rejected('+ str(rejected) + ' events rejected):')

    print('\n\t| Selection                       | Events remaining | Acceptance        | Selection Efficiency|')
    print('\t|---------------------------------|------------------|-------------------|---------------------|')
    print('\t| Events reconstructed            | ' + str(veto[0]) + '             | ' + str(acceptance[0]) + '   |            |')
    print('\t| Reduced chi squared < ' + str(chi2Cut) + '         | ' + str(veto[1]) + '             | ' + str(acceptance[1]) + '     |            |')
    print('\t| No. of track measurements > ' + str(measCut) + '  | ' + str(veto[2]) + '             | ' + str(acceptance[2]) + '     |            |')
    print('\t| Decay vertex in fiducial volume | ' + str(veto[3]) + '             | ' + str(acceptance[3]) +  '     |            |')
    print('\t| Both tracks in fiducial volume  | ' + str(veto[4]) + '             | ' + str(acceptance[4]) + '     |            |')
    print('\t| Each track > ' + str(ecalCut) + ' GeV in ECAL   | ' + str(veto[5]) + '             | ' + str(acceptance[5]) +  '     |            |')
    print('\t| Muon hits in 1st & 2nd stations | ' + str(veto[6]) + '             | ' + str(acceptance[6]) + '     |            |')
    print('\t| DOCA < ' + str(docaCut) + ' cm                   | ' + str(veto[7]) + '             | ' + str(acceptance[7]) + '     |            |')
    print('\t| IP to target < ' + str(ipCut) + ' cm           | ' + str(veto[8]) + '             | ' + str(acceptance[8]) + '     |            |')
    print('\n\t' + str(k2mu_MotherHP) + ' kaons (mother HP) decayed to muons before detection (' + str(k2mu_MotherHP-HP2ka_veto[0]) + ' after track checks)')
    print(2*' ' + 80*'_'+ '\n')

    h['Mu' + 'ProbMeasr'] = createRatio(h['Mu' + 'SmearedMass'],h['TotalSmearedMass'],'Mu' + 'ProbMeasr')
    h['K+/-' + 'ProbMeasr'] = createRatio(h['K+/-' + 'SmearedMass'],h['TotalSmearedMass'],'K+/-' + 'ProbMeasr')

    makePlots2(particleList)

def finStateMuKa_exc():
    if sTree.GetBranch('FitTracks'):
        print('\nRunning analysis for final state K*+- Mu-+ :\n')
        print('Decay A: N --> K*+- mu-+ --> K+- pi0 mu-+')
        print('Decay B: N --> K*+- mu-+ --> K0 pi+- mu-+\n')
        createHists_MuKa_exc()   # calls function to create histograms
        successful_events = []   # creates list of event numbers of desired decays
        A_veto = 10*[0]   # veto counter for decay A
        B_veto = 10*[0]   # veto counter for decay B
        pi_veto = 10*[0]   # creates list of veto counts for muons which decayed from pions from K* from neutrinalinos
        pi_decaycheck = 0   # variable for counting when pions from K* from neutralinos decay to muons before detection

        decayAcount_kaon = 0
        decayAcount_pion0 = 0
        decayBcount_pion = 0
        decayBcount_kaon0 = 0 

        for n in range(nEvents):   # loops over events
            rc = sTree.GetEntry(n)   # loads tree entry
            event = True
            motherEcheck = []

            #-----------------------------------------------TRACK-LOOPS------------------------------------------------

            for index,reco_muon in enumerate(sTree.FitTracks):   # loops over index and data of track particles                                   
                muPartkey = sTree.fitTrack2MC[index]   # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]   # gives MC particle data

                if abs(true_muon.GetPdgCode()) == 13:   # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()   # stores a number index of MC track of mother
                    true_motherN = sTree.MCTrack[muonMotherkey]   # obtains mother particle data
                    if true_motherN.GetPdgCode() == 9900015:

                        fitstatus_muon = reco_muon.getFitStatus()
                        if not fitstatus_muon.isFitConverged(): continue
                        
                        Muon_4Mom = SingleTrack_4Mom(index)

                        for index2,reco_part in enumerate(sTree.FitTracks):   # loops over index and data of track 
                            Partkey = sTree.fitTrack2MC[index2]   # matches track to MC particle key
                            true_part = sTree.MCTrack[Partkey]   # gives MC particle data

                            #===============================================DECAY=A========================================================

                            if abs(true_part.GetPdgCode()) == 321:   # checks particle is CHARGED KAON
                                true_kaon = true_part
                                reco_kaon = reco_part
                                kaPartkey = Partkey
                                kaonMotherkey = true_kaon.GetMotherId()
                                true_kaonEx = sTree.MCTrack[kaonMotherkey]
                                if abs(true_kaonEx.GetPdgCode()) == 323:   # checks mother is excited charged kaon
                                    kaonExMotherkey = true_kaonEx.GetMotherId()           
                                    true_motherN = sTree.MCTrack[kaonExMotherkey]
                                    if kaonExMotherkey == muonMotherkey and true_motherN.GetPdgCode() == 9900015:   # checks if mother keys are the same

                                        fitstatus_kaon = reco_kaon.getFitStatus()
                                        if not fitstatus_kaon.isFitConverged(): continue
                                        
                                        check,A_veto = track_checks(index,A_veto,0)   # performs various checks
                                        if check == -1: continue
                                        
                                        decayAcount_kaon += 1

                            #--------------------------------------------------------------------------------------------------------------

                                        for hit in sTree.EcalPoint:
                                            if hit.GetPdgCode() == 22:   # looking for photon hits in the ECAL
                                                TrackID = hit.GetTrackID()
                                                
                                                if TrackID == 7 or TrackID == 8:   # both correspond to the photons we are looking for (somehow)
                                                    photon = sTree.MCTrack[TrackID]
                                                    pion0 = sTree.MCTrack[photon.GetMotherId()]
                                                    if pion0.GetPdgCode() == 111:   # photon mother is neutral pion
                                                        pion0Motherkey = pion0.GetMotherId()
                                                        kaonEx = sTree.MCTrack[pion0Motherkey]
                                                        if kaonEx.GetPdgCode() == 323:   # neutral pion mother is excited kaon
                                                            motherN = sTree.MCTrack[kaonEx.GetMotherId()]
                                                            if pion0Motherkey == kaonMotherkey and motherN.GetPdgCode() == 9900015:
                                                                motherEcheck.append(pion0.GetEnergy())

                                        if len(motherEcheck) > 0: decayAcount_pion0 += 1

                            #===============================================DECAY=B========================================================

                            if abs(true_part.GetPdgCode()) == 13:
                                true_mother = sTree.MCTrack[true_part.GetMotherId()]
                                if abs(true_mother.GetPdgCode()) == 211:
                                    motherKaExc = sTree.MCTrack[true_mother.GetMotherId()]
                                    if abs(motherKaExc.GetPdgCode()) == 323:
                                        fitstatus = reco_part.getFitStatus()
                                        if fitstatus.isFitConverged():
                                            motherN = sTree.MCTrack[motherKaExc.GetMotherId()]
                                            if abs(motherN.GetPdgCode()) == 9900015:
                                                pi_decaycheck += 1   # pion has decayed to muon in flight
                                                check,pi_veto = track_checks(index2,pi_veto,0)
                                                if check == -1: pi_veto[0] += 1
                            
                            if abs(true_part.GetPdgCode()) == 211:   # checks particle is CHARGED PION
                                true_pion = true_part
                                reco_pion = reco_part
                                piPartkey = Partkey
                                pionMotherkey = true_pion.GetMotherId()
                                true_kaonEx = sTree.MCTrack[pionMotherkey]
                                if abs(true_kaonEx.GetPdgCode()) == 323:   # checks mother is excited charged kaon
                                    kaonExMotherkey = true_kaonEx.GetMotherId()
                                    true_motherN = sTree.MCTrack[kaonExMotherkey]
                                    if kaonExMotherkey == muonMotherkey and true_motherN.GetPdgCode() == 9900015:

                                        fitstatus_pion = reco_pion.getFitStatus()
                                        if not fitstatus_pion.isFitConverged(): continue

                                        check,pi_veto = track_checks(index2,pi_veto,0)   # just using pi_veto as we don't need it for counting
                                        if check == 0: decayBcount_pion += 1   # doing track_checks here as well to make this count correct
                                        
                                        Pion_4Mom = SingleTrack_4Mom(index2)
                                            
                            #--------------------------------------------------------------------------------------------------------------

                                        for index3,reco_piplus in enumerate(sTree.FitTracks):   # loops over index and data of track particles                                   
                                            piplusPartkey = sTree.fitTrack2MC[index3]   # matches track to MC particle key
                                            true_piplus = sTree.MCTrack[piplusPartkey]   # gives MC particle data
                            
                                            if true_piplus.GetPdgCode() == 211:   # checks particle is pion
                                                piplusMotherkey = true_piplus.GetMotherId()   # stores a number index of MC track of mother
                                                true_kaon = sTree.MCTrack[piplusMotherkey]
                                                if abs(true_kaon.GetPdgCode()) == 310 or abs(true_kaon.GetPdgCode()) == 130:   # checks mother is NEUTRAL KAON
                                                    true_kaon2 = sTree.MCTrack[true_kaon.GetMotherId()]
                                                    if true_kaon2.GetPdgCode() == 310 or abs(true_kaon2.GetPdgCode()) == 130:                                     
                                                        true_kaonEx = sTree.MCTrack[true_kaon2.GetMotherId()]
                                                        if abs(true_kaonEx.GetPdgCode()) == 323:   # checks mother is charged excited kaon   
                                                            true_motherN = sTree.MCTrack[true_kaonEx.GetMotherId()]
                                                            if true_motherN.GetPdgCode() == 9900015:

                                                                for index4,reco_piminus in enumerate(sTree.FitTracks):   # loops over index and data of track 
                                                                    piminusPartkey = sTree.fitTrack2MC[index4]   # matches track to MC particle key
                                                                    true_piminus = sTree.MCTrack[piminusPartkey]   # gives MC particle data
                                                    
                                                                    if true_piminus.GetPdgCode() == -211:   # checks particle is oppositely charged pion
                                                                        piminusMotherkey = true_piminus.GetMotherId()   # stores a number index of MC track of mother
                                                        
                                                                        if piplusMotherkey == piminusMotherkey:

                                                                            fitstatus_piplus = reco_piplus.getFitStatus() 
                                                                            fitstatus_piminus = reco_piminus.getFitStatus()
                                                                            if fitstatus_piplus.isFitConverged() and fitstatus_piminus.isFitConverged():
                                                                                B_veto[0] += 1
                                                                            else:
                                                                                #print('At least one of the track fits did not converge')
                                                                                continue

                                                                            check,B_veto = track_checks(index,B_veto,1)   # performs various track checks
                                                                            if check == -1: event = False
                                                                            if event == True: 
                                                                                check,B_veto = track_checks(index2,B_veto,1)
                                                                                if check == -1: event = False
                                                                            if event == True:
                                                                                check,B_veto = track_checks(index3,B_veto,1)
                                                                                if check == -1: event = False
                                                                            if event == True:
                                                                                check,B_veto = track_checks(index4,B_veto,1)
                                                                                if check == -1: event = False

                                                                            #---------------------------------------------PARTICLE-DATA-----------------------------------------------------

                                                                            Kaon_4Mom,Pionplus_4Mom,Pionminus_4Mom,Decay_X,Decay_Y,Decay_Z,doca = RedoVertexing(index3,index4)   # uses RedoVertexing to iterate track fitting
                                                                            if Kaon_4Mom == -1: 
                                                                                print('RedoVertexing extrapolation failed (event ' + str(n) + ')')
                                                                                B_veto[9] += 1
                                                                                continue

                                                                            KaonEx_4Mom = Kaon_4Mom + Pion_4Mom
                                                                            RPV_4Mom = KaonEx_4Mom + Muon_4Mom

                                                                            h['recovertex'].Fill(Decay_Z)
                                                                            if not isInFiducial(Decay_X,Decay_Y,Decay_Z):
                                                                                #print('Neutralino decayed outside fiducial volume')
                                                                                B_veto[6] += 1
                                                                                event = False

                                                                            h['doca'].Fill(doca)
                                                                            if not doca < docaCut: 
                                                                                #print('Distance of closest approach too large')
                                                                                B_veto[7] +=1
                                                                                event = False

                                                                            NProduction_X = true_motherN.GetStartX()
                                                                            NProduction_Y = true_motherN.GetStartY()
                                                                            NProduction_Z = true_motherN.GetStartZ()
                                                                            NProduction_T = true_motherN.GetStartT()
                                                                            RPV_Pos = ROOT.TLorentzVector()
                                                                            RPV_Pos.SetXYZT(NProduction_X,NProduction_Y,NProduction_Z,NProduction_T)
                                                                            tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
                                                                            ip = ImpactParameter(tr,RPV_Pos,RPV_4Mom)   # gives the same result as line 706 in ShipAna.py (i.e. using sTree.Particles)
                                    
                                                                            #Decay_T = true_muon.GetStartT()
                                                                            #RPV_Pos = ROOT.TLorentzVector()
                                                                            #RPV_Pos.SetXYZT(Decay_X,Decay_Y,Decay_Z,Decay_T)
                                                                            #tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
                                                                            #ip = ImpactParameter(tr,RPV_Pos,RPV_4Mom)

                                                                            h['IP_target'].Fill(ip)
                                                                            if not ip < ipCut:
                                                                                #print('Neutralino impact parameter to target too large')
                                                                                B_veto[8] += 1
                                                                                event = False

                                                                            if event == False: continue

                                                                            decayBcount_kaon0 += 1

                                                                            Kaon_recomass = Kaon_4Mom.M()   # reconstructed neutral kaon mass
                                                                            Kaon_recomom = Kaon_4Mom.P()   # reconstructed neutral kaon momentum
                                                                            Kaon_truemom = true_kaon.GetP()
                                                                            piplus_recomom = Pionplus_4Mom.P()   # reconstructed pi+ momentum
                                                                            piminus_recomom = Pionminus_4Mom.P()   # reconstructed pi- momentum
                                                                            h['Kaon_recomass'].Fill(Kaon_recomass)           
                                                                            h['Kaon_recomom'].Fill(Kaon_recomom)
                                                                            h['Kaon_truemom'].Fill(Kaon_truemom)
                                                                            h['Piplus_recomom'].Fill(piplus_recomom)
                                                                            h['Piminus_recomom'].Fill(piminus_recomom)

                                                                            RPV_recomass = RPV_4Mom.M()
                                                                            RPV_recomom = RPV_4Mom.P()
                                                                            RPV_truemom = true_motherN.GetP()
                                                                            h['RPV_recomass'].Fill(RPV_recomass)
                                                                            h['RPV_recomom'].Fill(RPV_recomom)
                                                                            h['RPV_truemom'].Fill(RPV_truemom)

                                                                            successful_events.append(n)   # adds entries to the list

        #----------------------------------------------------------------VETO-COUNTS------------------------------------------------------------------                                 

        accepted = len(successful_events)
        rejected = B_veto[0] - accepted
        print('\n\t' + str(B_veto[0]) + ' charged final states reconstructed for this decay mode')
        print('\t' + str(accepted) + ' not rejected')
        print('\t' + str(rejected) + ' rejected:')
        print('\t\t' + str(B_veto[0] - B_veto[6]) + ' events with K0 decay vertex inside fiducial volume (' + str(B_veto[6]) + ' event vetos)')
        print('\t\t' + str(B_veto[0] - B_veto[1]) + ' events with all tracks within fiducial volume (' + str(B_veto[1]) + ' event vetos)')
        print('\t\t' + str(B_veto[0] - B_veto[2]) + ' events with all tracks no. of measurements > ' + str(measCut) + ' (' + str(B_veto[2]) + ' event vetos)')
        print('\t\t' + str(B_veto[0] - B_veto[3]) + ' events with all tracks chi squared < ' + str(chi2Cut) + ' (' + str(B_veto[3]) + ' event vetos)')
        print('\t\t' + str(B_veto[0] - B_veto[4]) + ' events where all tracks leave > ' + str(ecalCut) + ' GeV in ECAL (' + str(B_veto[4]) + ' event vetos)')
        print('\t\t' + str(B_veto[0] - B_veto[5]) + ' events with hits in both the 1st and 2nd muon stations (' + str(B_veto[5]) + ' event vetos)')
        print('\t\t' + str(B_veto[0] - B_veto[7]) + ' events with DOCA between tracks < ' + str(docaCut) + ' cm (' + str(B_veto[7]) + ' event vetos)')
        print('\t\t' + str(B_veto[0] - B_veto[8]) + ' events with IP to target < ' + str(ipCut) + ' cm (' + str(B_veto[8]) + ' event vetos)')

        print('\t' + str(B_veto[9]) + ' failed RedoVertexing extrapolations')
        print('\t' + str(pi_decaycheck) + ' pions decayed to muons before detection (' + str(pi_decaycheck - pi_veto[0]) + ' after track checks)')

        print('\n\t------------------------------------')

        print('\n\t' + str(decayAcount_kaon) + ' charged kaons detected (decay A)')   # after checks and track converges
        print('\t' + str(decayAcount_pion0) + ' neutral pions detected (decay A)')
        print('\n\t' + str(decayBcount_pion) + ' charged pions detected (decay B)')
        print('\t' + str(decayBcount_kaon0) + ' neutral kaons detected (decay B)\n')

        makePlots_MuKa_exc()

def finStateDarkPhot():
    if sTree.GetBranch('FitTracks'):
        print('\nRunning analysis for dark photon :\n')
        createHists_DarkPhot()
        successful_events = []
        veto = 10*[0]   # creates list of veto counts for each possible veto cause
        for n in range(nEvents):   # loops over events
            rc = sTree.GetEntry(n)   # loads tree entry
            event = True

            #-----------------------------------------------TRACK-LOOPS------------------------------------------------

            for index,reco_part in enumerate(sTree.FitTracks):
                eminusPartkey = sTree.fitTrack2MC[index]   # matches track to MC particle key
                true_eminus = sTree.MCTrack[eminusPartkey]   # gives MC particle data              

                if true_eminus.GetPdgCode() == 11:   # checks particle is electron
                    eminusMotherkey = true_eminus.GetMotherId()   # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[eminusMotherkey]   # obtains mother particle data
                    fitstatus_eminus = reco_part.getFitStatus() 
                    if not fitstatus_eminus.isFitConverged():
                        continue

                    if true_mother.GetPdgCode() == 9900015:   # checks mother is RPV 
                        for index2,reco_part2 in enumerate(sTree.FitTracks):   # loops over index and data of track particles
                            eplusPartkey = sTree.fitTrack2MC[index2]   # matches track to MC particle key
                            true_eplus = sTree.MCTrack[eplusPartkey]   # gives MC particle data
                            if true_eplus.GetPdgCode() == -11:   # checks particle is positron
                                eplusMotherkey = true_eplus.GetMotherId()   # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[eplusMotherkey]   # obtains mother particle data

                                if eplusMotherkey == eminusMotherkey and true_mother.GetPdgCode() == 9900015:   # check if mother keys are the same
                                    fitstatus_eplus = reco_part2.getFitStatus()
                                    if fitstatus_eplus.isFitConverged():
                                        veto[0] += 1
                                    else:
                                        #print('At least one of the track fits did not converge')
                                        continue
                                    
                                    #-------------------------------------------------TRACK-CHECKS-------------------------------------------------------

                                    DP_4Mom,eminus_4Mom,eplus_4Mom,NDecay_X,NDecay_Y,NDecay_Z,doca = RedoVertexing(index,index2)   # uses RedoVertexing to iterate track fitting
                                    if DP_4Mom == -1: 
                                        print('RedoVertexing extrapolation failed (event ' + str(n) + ')')
                                        veto[9] += 1
                                        continue

                                    nmeas_eplus = fitstatus_eplus.getNdf()
                                    chi2_eplus = fitstatus_eplus.getChi2()
                                    rchi2_eplus = (chi2_eplus/nmeas_eplus)
                                    nmeas_eminus = fitstatus_eminus.getNdf()
                                    chi2_eminus = fitstatus_eminus.getChi2()
                                    rchi2_eminus = (chi2_eminus/nmeas_eminus)
                                    h['Chi2'].Fill(rchi2_eplus)
                                    h['Chi2'].Fill(rchi2_eminus)
                                    if rchi2_eplus < chi2Cut and rchi2_eminus < chi2Cut:
                                        veto[3] += 1
                                    else: event = False

                                    h['nmeas'].Fill(nmeas_eplus)
                                    h['nmeas'].Fill(nmeas_eminus)
                                    if event == True:
                                        if nmeas_eplus > measCut and nmeas_eminus > measCut:
                                            veto[2] += 1
                                        else: event = False

                                    h['recovertex'].Fill(NDecay_Z)
                                    if event == True:
                                        if isInFiducial(NDecay_X,NDecay_Y,NDecay_Z):
                                            veto[6] += 1
                                        else: event = False

                                    if event == True:
                                        if checkFiducialVolume(sTree,index,dy) and checkFiducialVolume(sTree,index2,dy): 
                                            veto[1] += 1
                                        else: event = False

                                    ecalE_eplus = ecalMinIon(eplusPartkey)
                                    ecalE_eminus = ecalMinIon(eminusPartkey)
                                    h['ecalE'].Fill(ecalE_eplus)
                                    h['ecalE'].Fill(ecalE_eminus)
                                    if event == True:
                                        if ecalE_eplus > ecalCut and ecalE_eminus > ecalCut:
                                            veto[4] += 1
                                        else: event = False

                                    #if event == True:
                                    #    if muonstationHits(muPartkey):
                                    #        veto[5] += 1
                                    #    else: event = False
                                                    
                                    h['doca'].Fill(doca)
                                    if event == True:
                                        if doca < docaCut: 
                                            veto[7] +=1
                                        else: event = False

                                    NProduction_X = true_mother.GetStartX()   # vertex coordinates of neutralino production
                                    NProduction_Y = true_mother.GetStartY()
                                    NProduction_Z = true_mother.GetStartZ()
                                    NProduction_T = true_mother.GetStartT()
                                    DP_Pos = ROOT.TLorentzVector()
                                    DP_Pos.SetXYZT(NProduction_X,NProduction_Y,NProduction_Z,NProduction_T)
                                    tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
                                    ip = ImpactParameter(tr,DP_Pos,DP_4Mom)   # gives the same result as line 706 in ShipAna.py (i.e. using sTree.Particles)
                                    h['IP_target'].Fill(ip)
                                    if event == True:
                                        if ip < ipCut:
                                            veto[8] += 1
                                        else: event = False

                                    if event == False: continue

                                    #-------------------------------------------------PARTICLE-DATA------------------------------------------------------

                                    DP_recomass = DP_4Mom.M()   # reconstructed dark photon mass
                                    DP_truemom = true_mother.GetP()   # true dark photon momentum
                                    DP_recomom = DP_4Mom.P()   # reconstructed dark photon momentum
                                    DP_momdiff = DP_truemom - DP_recomom   # dark photon true/reco momentum difference
                                    true_eplusP = true_eplus.GetP()   # true positron momentum
                                    reco_eplusP = eplus_4Mom.P()   # reconstructed positron momentum
                                    true_eminusP = true_eminus.GetP()   # true electron momentum
                                    reco_eminusP = eminus_4Mom.P()   # reconstructed electron momentum
                                                      
                                    h['DP_recomass'].Fill(DP_recomass)   # fills histograms
                                    h['DP_truemom'].Fill(DP_truemom)             
                                    h['DP_recomom'].Fill(DP_recomom)                            
                                    h['DP_mom_diff'].Fill(DP_momdiff)
                                    h['eplus_recomom'].Fill(reco_eplusP)
                                    h['eminus_recomom'].Fill(reco_eminusP)
                                    h['eplus_truemom'].Fill(true_eplusP)
                                    h['eminus_truemom'].Fill(true_eminusP)

                                    DP_Zmom = DP_4Mom.Pz()
                                    DP_Zbeta = 1/ROOT.TMath.Sqrt(1 + ((DP_recomass/DP_Zmom)**2))
                                    DP_Zgamma = 1/(ROOT.TMath.Sqrt(1 - (DP_Zbeta**2)))
                                    cos_theta = DP_Zmom/DP_recomom
                                    theta = 1000*(ROOT.TMath.ACos(cos_theta))   # angle between beam line and neutralino momentum (mrad)

                                    h['DP_beta'].Fill(DP_Zbeta)
                                    h['DP_gamma'].Fill(DP_Zgamma)
                                    h['DP_theta'].Fill(theta)

                                    successful_events.append(n)   # adds entries to the list
                                    
        #----------------------------------------------------------------VETO-COUNTS------------------------------------------------------------------

        accepted = len(successful_events)
        print('\n\t' + str(veto[0]) + ' events reconstructed for this decay mode')
        print('\t' + str(veto[0] - veto[9]) + ' events with successful RedoVertexing extrapolations')
        print('\t' + str(accepted) + ' events not rejected:')

        print('\n\t| Selection                       | Events remaining | Efficiency | Acceptance |')
        print('\t|---------------------------------|------------------|------------|------------|')
        print('\t| Events reconstructed            | ' + str(veto[0]) + '             | ' + '           |            |')
        print('\t| Reduced chi squared < ' + str(chi2Cut) + '         | ' + str(veto[3]) + '             | ' + str(round(100*veto[3]/float(veto[0]),2)) + ' %    |            |')
        print('\t| No. of track measurements > ' + str(measCut) + '  | ' + str(veto[2]) + '             | ' + str(round(100*veto[2]/float(veto[3]),2)) + ' %    |            |')
        print('\t| Decay vertex in fiducial volume | ' + str(veto[6]) + '             | ' + str(round(100*veto[6]/float(veto[2]),2)) + ' %    |            |')
        print('\t| Both tracks in fiducial volume  | ' + str(veto[1]) + '             | ' + str(round(100*veto[1]/float(veto[6]),2)) + ' %    |            |')
        print('\t| Each track > ' + str(ecalCut) + ' GeV in ECAL   | ' + str(veto[4]) + '             | ' + str(round(100*veto[4]/float(veto[1]),2)) + ' %    |            |')
        #print('\t| Muon hits in 1st & 2nd stations | ' + str(veto[5]) + '             | ' + str(round(100*veto[5]/float(veto[4]),2)) + ' %    |            |')
        print('\t| DOCA < ' + str(docaCut) + ' cm                   | ' + str(veto[7]) + '             | ' + str(round(100*veto[7]/float(veto[4]),2)) + ' %    |            |')
        print('\t| IP to target < ' + str(ipCut) + ' cm           | ' + str(veto[8]) + '             | ' + str(round(100*veto[8]/float(veto[7]),2)) + ' %    |            |')

        makePlots_DarkPhot()

########################################
#########  MAIN SECTION KINDA  #########    
while loop:        
    print_menu()
    choice = input('Enter your choice [1-4]: ')
     
    if choice==1:     
        print ('\nRPV SUSY Benchmark1 --> K+/- mu+/- visible final state selected.')
        finStateMuKa()
        #-------------------------------------------OUTPUT--------------------------------------------
        hfile = inputFile.split(',')[0].replace('_rec','_RPVeditana')#create outputFile
        if hfile[0:4] == '/eos' or not inputFile.find(',')<0:
            # do not write to eos, write to local directory 
            tmp = hfile.split('/')
            hfile = tmp[len(tmp)-1]
        ROOT.gROOT.cd()
        f = TFile(hfile,'RECREATE')
        for akey in h:
            cln = h[akey].Class().GetName()
            if not cln.find('TH')<0 or not cln.find('TP')<0:   
                h[akey].Write()
        for bkey in graph:
            cln = graph[bkey].Class().GetName()
            if not cln.find('TG')<0 or not cln.find('TM')<0:   
                f.WriteTObject(graph[bkey])
        f.Close()
        loop=False

    elif choice==2:
        print ('\nRPV SUSY Benchmark1 --> K*+/- mu+/- visible final state selected.') 
        finStateMuKa_exc()
        if hfile[0:4] == '/eos' or not inputFile.find(',')<0:
        # do not write to eos, write to local directory 
          tmp = hfile.split('/')
          hfile = tmp[len(tmp)-1] 
        ROOT.gROOT.cd()
        ut.writeHists(h,hfile)
        loop=False

    elif choice==3:
        print ('Dark Photon Model selected')
        finStateDarkPhot()
        #-------------------------------------------OUTPUT--------------------------------------------
        hfile = inputFile.split(',')[0].replace('_rec','_DarkPhot')  # Outputs histograms and ROOT file
        if hfile[0:4] == '/eos' or not inputFile.find(',')<0:
        # do not write to eos, write to local directory 
          tmp = hfile.split('/')
          hfile = tmp[len(tmp)-1] 
        ROOT.gROOT.cd()
        ut.writeHists(h,hfile)   
        loop=False

    elif choice==4:
        print ('\nEXIT was selected.')
        loop=False 

    else:
        print('\nNo such option. Try again.')




