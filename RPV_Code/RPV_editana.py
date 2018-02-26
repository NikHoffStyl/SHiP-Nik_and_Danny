###########################################
#         RPV_EDITA Test Code             #
###########################################
#Code Performs Analysis on Data Obtained by simulation and reconstruction
#Author: Nikos Stylianou

#from __future__ import print_function
import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
from array import array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH1D, TPad, TGraph, TF1, TMultiGraph,TGraphErrors, THStack, TFile
from ROOT import kBlack, kBlue, kRed, kFALSE, kSolar
from ROOT import gROOT, gPad, gStyle
from array import array
import shipRoot_conf
import shipDet_conf
import shipVeto
from inspect import currentframe
import numpy as np
import datetime
shipRoot_conf.configure()

debug = False
loop=True 
chi2CutOff  = 4.
PDG = ROOT.TDatabasePDG.Instance()
fiducialCut = True
chi2Cut  = 4
measCut = 25
ecalCut = 0.150
measCutFK = 25
measCutPR = 22
docaCut = 2.
ipCut = 250

c = 2.99792458*(10**8)
e = 2.718281828459   # Euler's number
currentDate = datetime.datetime.now().strftime("%y_%m_%d_%H%M")
polyFit1 = TF1("polyFit1","pol3")
polyFit2 = TF1("polyFit2","pol3")

def inputOptsArgs():
    inputFile  = None
    geoFile    = None
    dy         = None
    nEvents    = 9999999
    try:
            opts, args = getopt.getopt(sys.argv[1:], "n:f:g:Y", ["nEvents=","geoFile="])#command line options
    except getopt.GetoptError:
            # print help information and exit:
            print (' enter file name')
            sys.exit()
    for o, a in opts:
            if o in ("-f",):
                inputFile = a           #does this
            if o in ("-g", "--geoFile",):
                geoFile = a             #does this
            if o in ("-Y",):
                dy = float(a)           #doesnt do this
            if o in ("-n", "--nEvents=",):
                nEvents = int(a)        #doesnt do this
    return a, dy, geoFile, inputFile, nEvents, o
a, dy, geoFile, inputFile, nEvents, o = inputOptsArgs()

###############################################
########## SETS TREE FROM INPUT FILE ##########
if not inputFile.find(',')<0 :
    sTree = ROOT.TChain("cbmsim") 
    for x in inputFile.split(','):
        if x[0:4] == "/eos":         #doesnt do this
            sTree.AddFile("root://eoslhcb.cern.ch/"+x)
        else: sTree.AddFile(x)       #doesnt do this
elif inputFile[0:4] == "/eos":       #doesnt do this
    eospath = "root://eoslhcb.cern.ch/"+inputFile
    f = ROOT.TFile.Open(eospath)
    sTree = f.cbmsim
else:                                #does this
    f = ROOT.TFile(inputFile)
    sTree = f.cbmsim

#######################################
########## LOAD ECAL GEOMETRY #########
if not geoFile:
 geoFile = inputFile.replace('ship.','geofile_full.').replace('_rec.','.')  #doesnt do this
if geoFile[0:4] == "/eos":                                                  #doesnt do this
  eospath = "root://eoslhcb.cern.ch/"+geoFile
  fgeo = ROOT.TFile.Open(eospath)
else:                                                                       #does this
    fgeo = ROOT.TFile(geoFile)
sGeo = fgeo.FAIRGeom

if not fgeo.FindKey('ShipGeo'):                 #doesnt do this
 # old geofile, missing Shipgeo dictionary
 if sGeo.GetVolume('EcalModule3') :  
     ecalGeoFile = "ecal_ellipse6x12m2.geo"     #doesnt do this
 else: 
     ecalGeoFile = "ecal_ellipse5x10m2.geo"     #doesnt do this
 print ('found ecal geo for ',ecalGeoFile)
 # re-create geometry and mag. field
 if not dy:
  # try to extract from input file name
  tmp = inputFile.split('.')                    #doesnt do this
  try:
    dy = float( tmp[1]+'.'+tmp[2] )             #doesnt do this
  except:
    dy = 10.                                    #doesnt do this
 ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py", Yheight = dy, EcalGeoFile = ecalGeoFile )
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
ecalGeo = ecalGeoFile+'z'+str(ShipGeo.ecal.z)+".geo"
if not ecalGeo in os.listdir(os.environ["FAIRSHIP"]+"/geometry"): shipDet_conf.makeEcalGeoFile(ShipGeo.ecal.z,ShipGeo.ecal.File)
ecalFiller = ROOT.ecalStructureFiller("ecalFiller", 0,ecalGeo)
ecalFiller.SetUseMCPoints(ROOT.kTRUE)
ecalFiller.StoreTrackInformation()
ecalStructure = ecalFiller.InitPython(sTree.EcalPointLite)
caloTasks.append(ecalFiller)
if sTree.GetBranch("EcalReconstructed"):
 calReco = False
 sTree.GetEvent(0)
 ecalReconstructed = sTree.EcalReconstructed
else:
 calReco = True
 print ("setup calo reconstruction of ecalReconstructed objects")
# Calorimeter reconstruction
 #GeV -> ADC conversion
 ecalDigi=ROOT.ecalDigi("ecalDigi",0)
 ecalPrepare=ROOT.ecalPrepare("ecalPrepare",0)
 ecalStructure     = ecalFiller.InitPython(sTree.EcalPointLite)
 ecalDigi.InitPython(ecalStructure)
 caloTasks.append(ecalDigi)
 ecalPrepare.InitPython(ecalStructure)
 caloTasks.append(ecalPrepare)
 # Cluster calibration
 ecalClusterCalib=ROOT.ecalClusterCalibration("ecalClusterCalibration", 0)
 #4x4 cm cells
 ecalCl3PhS=ROOT.TFormula("ecalCl3PhS", "[0]+x*([1]+x*([2]+x*[3]))")
 ecalCl3PhS.SetParameters(6.77797e-04, 5.75385e+00, 3.42690e-03, -1.16383e-04)
 ecalClusterCalib.SetStraightCalibration(3, ecalCl3PhS)
 ecalCl3Ph=ROOT.TFormula("ecalCl3Ph", "[0]+x*([1]+x*([2]+x*[3]))+[4]*x*y+[5]*x*y*y")
 ecalCl3Ph.SetParameters(0.000750975, 5.7552, 0.00282783, -8.0025e-05, -0.000823651, 0.000111561)
 ecalClusterCalib.SetCalibration(3, ecalCl3Ph)
#6x6 cm cells
 ecalCl2PhS=ROOT.TFormula("ecalCl2PhS", "[0]+x*([1]+x*([2]+x*[3]))")
 ecalCl2PhS.SetParameters(8.14724e-04, 5.67428e+00, 3.39030e-03, -1.28388e-04)
 ecalClusterCalib.SetStraightCalibration(2, ecalCl2PhS)
 ecalCl2Ph=ROOT.TFormula("ecalCl2Ph", "[0]+x*([1]+x*([2]+x*[3]))+[4]*x*y+[5]*x*y*y")
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

#########  END CREATE GEOMETRY #########
########################################

import TrackExtrapolateTool

########################################
############  DEFINITIONS  #############
h = {}
#stack={}
def create_Hists(HiddPart,part1,part2, part3):
    dictionList={part1:4, part2:2}
    partList = [part1,part2]
    if not part3 == None:
        partList.append(part3)
        dictionList.update({part3:3})
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
        ut.bookHist(h,partName + 'StrawTime',partName + ' ; ' + partName + ' Straw Time [ns] ; No. of Particles',300,321.7,323.5)                                          # straw time
        h[partName + 'StrawTime'].SetLineColor(val)
        ut.bookHist(h,partName + 'EcalTime',partName + ' ; ' + partName + ' Ecal Time [ns] ; No. of Particles',300,359.7,361.2)                                            # ecal time
        h[partName + 'EcalTime'].SetLineColor(val)
        ut.bookHist(h,partName + 'DirDeltaTimeSmeared',partName + ' ; ' + partName + ' Straw-ECAL Smeared Time of Flight (directly) [ns] ; No. of Particles',300,37.8,38.4)# smeared time of flight
        h[partName + 'DirDeltaTimeSmeared'].SetLineColor(val)
        ut.bookHist(h,partName + 'DirDeltaTime',partName + ' ; ' + partName + ' Straw-ECAL Time of Flight (directly) [ns] ; No. of Particles',300,37.8,38.4)               # time of flight
        h[partName + 'DirDeltaTime'].SetLineColor(val)
        ut.bookHist(h,partName + 'StrawHits',partName + ' ; ' + partName + ' No. of hits in straw tubes ; Position [cm]',300,25,50)                                        #number of straw hits
        h[partName + 'StrawHits'].SetLineColor(val)
        ut.bookHist(h,partName + 'StrawHitsMom',partName + ' ; ' + partName + ' z-momentum through straw tubes (for particular event) [GeV/c] ; No. of Hits',500,0,45)     #momenta of straw hits
        h[partName + 'StrawHitsMom'].SetLineColor(val)
        ut.bookHist(h,partName + 'FlightLen',partName + ' ; ' + partName + ' Straw-ECAL Straight Flight Lenght [cm] ; No. of Particles',300,11.375,11.42)                  # flight Length
        h[partName + 'FlightLen'].SetLineColor(val)
        ut.bookHist(h,partName + 'FlightLenImproved',partName + ' ; ' + partName + ' Straw-ECAL Curved Flight Lenght [cm] ;  No. of Particles',300,11.375,11.42)            # corrected flight Length        
        h[partName + 'FlightLenImproved'].SetLineColor(val)
        ut.bookHist(h,partName + 'FlightLenDelta',partName + ' ; ' + partName + ' Difference between straight path and better approximation [cm] ; No. of Particles',300,0,0.001)# delta flight Length        
        h[partName + 'FlightLenDelta'].SetLineColor(val)
        ut.bookHist(h,partName + 'SpeedSmeared',partName + ' ; ' + partName + ' Smeared Beta value ; No. of Particles',300,0.997,1.001)                                     # smeared speed
        h[partName + 'SpeedSmeared'].SetLineColor(val)
        ut.bookHist(h,partName + 'Speed',partName + ' ; ' + partName + ' Beta value ; No. of Particles',300,0.997,1.001)                                                    # speed
        h[partName + 'Speed'].SetLineColor(val)
        ut.bookHist(h,partName + 'StrawMom',partName + ' ; ' + partName + ' Straw Momentum [GeV/c] ; No. of Particles',300,-0.05,120.)                                      # straw momentum
        h[partName + 'StrawMom'].SetLineColor(val)
        ut.bookHist(h,partName + 'EcalMom',partName + ' ; ' + partName + ' Ecal Momentum [GeV/c]; No. of Particles',300,-0.05,120.)                                         # ecal  momentum
        h[partName + 'EcalMom'].SetLineColor(val)
        ut.bookHist(h,partName + 'DeltaMom',partName + ' ; ' + partName + ' Straw-Ecal Momentum [GeV/c]; No. of Particles',300,0.02,0.13)                                   # delta momentum
        h[partName + 'DeltaMom'].SetLineColor(val)
        ut.bookHist(h,partName + 'RecoMom',partName + ' ; ' + partName + ' Reco Momentum [GeV/c] ; No. of Particles',300,-0.05,120.)                                        # reco  momentum
        h[partName + 'RecoMom'].SetLineColor(val)
        ut.bookHist(h,partName + 'TrueMom',partName + ' ; ' + partName + ' True Momentum [GeV/c] ; No. of Particles',300,-0.05,120.)                                        # true  momentum
        h[partName + 'TrueMom'].SetLineColor(val)
        ut.bookHist(h,partName + 'RecoMass',partName + ' ; ' + partName + ' Reco Mass [GeV/c2]; No. of Particles',300,0.,0.6)                                               # reco  mass
        h[partName + 'RecoMass'].SetLineColor(val)
        ut.bookHist(h,partName + 'TrueMass',partName + ' ; ' + partName + ' True Mass [GeV/c2] ; No. of Particles',300,0.,0.6)                                              # true  mass
        h[partName + 'TrueMass'].SetLineColor(val)
        h[partName + 'SmearedMass'] = TH1D(partName + 'SmearedMass',partName + ' ; ' + partName + ' Smeared Mass [GeV/c2]; No. of Particles',85,array('d',edgesarray))      # smrd  mass
        h[partName + 'SmearedMass'].SetLineColor(val)
        h[partName + 'ProbMeasr'] = TH1D(partName + 'ProbMeasr',partName + ' ; Mass [GeV/c2] ; Prob(particle = ' + partName + ')',85,array('d',edgesarray))                     # ID Prob
        h[partName + 'ProbMeasr'].SetLineColor(val)
                              
    h['TotalSmearedMass'] = TH1D('TotalSmearedMass','Smeared Mass ; Smeared Mass [GeV/c2] ; No. of Particles',85,array('d',edgesarray))                             # Total mass
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
    ut.bookHist(h,HiddPart + 'TrueMass','Monte Carlo Mass ; Invariant mass [GeV/c2] ; No. of Particles',300,0.99,1.01)                                              # true mass
    h[HiddPart + 'TrueMass'].SetLineColor(1)
    ut.bookHist(h,HiddPart + 'RecoMass','Reconstructed Mass ; Invariant mass [GeV/c2] ; No. of Particles',300,0.97,1.03)                                            # reco mass
    h[HiddPart + 'RecoMass'].SetLineColor(1)
    ut.bookHist(h,HiddPart + 'TrueMom','True (red) & Reco. (blue) Momentum ; Momentum [GeV/c] ; No. of Particles',100,0.,180.)                                      # true momentum 
    h[HiddPart + 'TrueMom'].SetLineColor(1)
    ut.bookHist(h,HiddPart + 'RecoMom','Reconstructed Momentum ; Momentum [GeV/c] ; No. of Particles',300,0.,180.)                                                  # reco momentum
    h[HiddPart + 'RecoMom'].SetLineColor(1)
    ut.bookHist(h,HiddPart + 'DeltaMom','True/Reco Momentum Difference ; Momentum Difference [GeV/c] ; No. of Particles',300,-3.,3)                                 # true-reco momentum difference
    h[HiddPart + 'DeltaMom'].SetLineColor(1)

    #######################
    ####  Veto Checks  ####
    ut.bookHist(h,'IP_target','Impact parameter to target; Impact Parameter [cm]; Frequency',300,0,10)
    ut.bookHist(h,'ecalE','Energy deposited in ECAL ; Energy [GeV/c2] ; Frequency',300,0,100)
    ut.bookHist(h,'doca','Distance of closest approach between muon and kaon tracks ; Distance [cm] ; Frequency',300,0,3)
    ut.bookHist(h,'nmeas','No. of measurements in fitted tracks (ndf) ; ndf ; No. of tracks',300,0,50)
    ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared ; Reduced Chi Squared ; Frequency',300,0,3)
    ut.bookHist(h,'recovertex','Reconstructed neutralino decay vertex z-coordinate ; Z  [cm] ; Frequency',100,-4000,4000)

    #ut.bookHist(h,HiddPart + '_no_iter','Reconstructed Mass (without track iterations)',500,0.,2.)   # reco mass(without track itrns)
    #ut.bookHist(h,'normdistr','Gaussian Distribution',500,-0.05,0.05)                               #
    #ut.bookHist(h,'smearedmass1','Time Smeared Neutralino Mass',500,0.,2.)
    #ut.bookHist(h,'smearedmass2','Time Smeared Neutralino Mass',500,0.,2.)
    #ut.bookHist(h,'smearedP1','Time Smeared Neutralino Momentum P1(red) P2(blue)',500,0.,300.)
    #ut.bookHist(h,'smearedP2','Time Smeared Neutralino Momentum',500,0.,300.)

    print("Created Histograms")
    print(partList)
    print('\n')

def dist2InnerWall(X,Y,Z):
  dist = 0
 # return distance to inner wall perpendicular to z-axis, if outside decayVolume return 0.
  node = sGeo.FindNode(X,Y,Z)
  if ShipGeo.tankDesign < 5:
     if not 'cave' in node.GetName():
         return dist  # TP 
  else:
     if not 'decayVol' in node.GetName():
        return dist
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
    if ShipGeo.tankDesign < 5 and nxt.GetName().find('I')<0:
        return 0    
    distance = sGeo.GetStep()
    if distance < minDistance  :
       minDistance = distance
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
                 print ('SHiPAna: extrapolation did not work (@RedoVertexing).')
                 rc = False
                 break
             newPosDir[tr] = [reps[tr].getPos(states[tr]),reps[tr].getDir(states[tr])]
         if not rc: break
         xv,yv,zv,doca = myVertex(t1,t2,newPosDir)
         dz = abs(zBefore-zv)
         step+=1
         if step > 10:
             print ('Abort iteration, too many steps, pos=',xv,yv,zv,' doca=',doca,'z before and dz',zBefore,dz)
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
     #return xv,yv,zv,doca,NeutralinoMom
     return NeutralinoLV,LV[t1],LV[t2],xv,yv,zv,doca

def time_res(partMom,partName,partkey,pdg,eventN,succEventM):
    smearStrawTime = None
    smearEcalTime = None
    smearDeltaT = None
    deltaT = None
    r = None
    R = None
    smearV = None
    v = None
    strawP = None
    ecalP = None
    deltaP = None
    if sTree.GetBranch("strawtubesPoint"):
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
            #smearStrawTime=straw_time 

        strawPx = px_array[min_z_index]
        strawPy = py_array[min_z_index]
        strawPz = pz_array[min_z_index]
        strawP = ROOT.TMath.Sqrt((strawPx**2) + (strawPy**2) + (strawPz**2)) 

        num_hits = len(z_array)   # number of elements in the list
        if pdg==13:   # muon
            h[partName + 'StrawHits'].Fill(num_hits)
            for hit in pz_array:
                if eventN == succEventM:
                    h[partName + 'StrawHitsMom'].Fill(hit)   # muon z-momentum through straw tubes for particular event
        if pdg==321:   # kaon
            h[partName + 'StrawHits'].Fill(num_hits)
            for hit in pz_array:
                if eventN == succEventM:
                    h[partName + 'StrawHitsMom'].Fill(hit)   # kaon z-momentum through straw tubes for particular event

        if sTree.GetBranch("EcalPoint"):
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
                            #smearEcalTime=ecal_time
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

    if smearStrawTime != None:
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
        if smearBeta < 1:
            #TRYING SOMETHING ELSE
            smearedM = partMom*(ROOT.TMath.Sqrt(1-(smearBeta**2)))/smearBeta
            h[partName + 'SmearedMass'].Fill(smearedM)
            h['TotalSmearedMass'].Fill(smearedM)   

def createRatio(h1, h2, histname):
    h3 = h1.Clone(histname)
    h3.SetMarkerStyle(20)
    h3.SetMarkerSize(0.7)
    h3.SetTitle("")
    h3.SetMinimum(0.8)
    h3.SetMaximum(1.35)
    # Set up plot for markers and errors
    h3.Sumw2()
    h3.SetStats(0)
    #h3.Divide(h2)
    h3.Divide(h1,h2,1,1,"B")
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetRangeUser(-0.1,1.2)
    y.SetTitleOffset(1.)
    ## Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetRangeUser(0,1.5)
    return h3

def track_checks(index,partkey,true_part,reco_part,veto,fill):
    check = 0

    if not checkFiducialVolume(sTree,index,dy): 
        #print('Track outside fiducial volume')
        veto[1] += 1
        check = -1

    fit_status = reco_part.getFitStatus() 
    fit_nmeas = fit_status.getNdf()
    if fill == 1: h['nmeas'].Fill(fit_nmeas)
    if not fit_nmeas > measCut:
        #print('Too few measurements')
        veto[2] += 1
        check = -1

    fit_rchi2 = fit_status.getChi2()                      
    fit_chi2 = (fit_rchi2/fit_nmeas)
    if fill == 1: h['Chi2'].Fill(fit_chi2)
    if not fit_chi2 < chi2Cut:
        #print('Chi squared value too high')
        veto[3] += 1
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
            veto[4] += 1
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
                veto[5] += 1
                check = -1
            else:
                if not muonstation_z[0] == 4017.0: #first muon station position
                    veto[5] += 1
                    check = -1

    return check,veto
nEvents = min(sTree.GetEntries(),nEvents)

def finState2t1t2(HiddPart,daught1,daught2):
    if HiddenPart == 'Neutralino':
        HiddPart_PDG = 9900015
    if daught1 == 'Mu':
        daught1_PDG = 13
    if  daught2 ==  'K+/-':
        daught2_PDG = 321
    elif daught2 == 'K*+/-':
        daught2_PDG = 323
    elif daught2 == 'K*0':
        daught2_PDG = 313
    if sTree.GetBranch("FitTracks"):
        totalEvents=0
        successful_events = []   # creates list of event numbers of desired decays
        daught1Events=0
        d1veto=10*[0]
        d1EventsAfterChecks=0
        k2muEvents=0
        k2mu_MotherHP = 0
        k2mu_Motherd2 = 0
        pi2muEvents=0
        pi2mu_MotherHP = 0
        pi2mu_Motherd2 = 0
        HiddPartDaught1Events=0
        daught2toKaonEvents=0
        d2toKaonveto=10*[0]
        d2toKaonEventsAfterChecks=0
        daught2toPionEvents=0
        d2toPionveto=10*[0]
        d2toPionEventsAfterChecks=0
        HP2ka_veto = 10*[0]
        d2ka_veto = 10*[0]
        HP2pi_veto = 10*[0]
        d2pi_veto = 10*[0]
        veto = 9*[0]
        
        for n in range(nEvents):
            totalEvents+=1
            rc = sTree.GetEntry(n)                              # load tree entry
            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles                                   
                d1Partkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_daught1 = sTree.MCTrack[d1Partkey]               # gives MC particle data
                if abs(true_daught1.GetPdgCode()) == daught1_PDG:        # checks particle is muon
                    daught1Events+=1
                    daught1MotherKey = true_daught1.GetMotherId()             # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[daught1MotherKey]             # obtains mother particle data

                    check, d1veto = track_checks(index,d1Partkey,true_daught1,reco_part,d1veto,0)
                    if not check == 0:   # performs various checks (i.e. vertex position, fiducial volume,...)
                        continue
                    d1EventsAfterChecks+=1

                    if true_mother.GetPdgCode() == 321:
                        d1GrannyKey = true_mother.GetMotherId()
                        true_gran = sTree.MCTrack[d1GrannyKey]
                        k2muEvents+=1
                        if true_gran.GetPdgCode() == HiddPart_PDG:
                            k2mu_MotherHP+=1
                            check,HP2ka_veto = track_checks(index,d1Partkey,true_daught1,reco_part,HP2ka_veto,0)
                            if check == -1: HP2ka_veto[0] += 1 
                        if true_gran.GetPdgCode() == daught2_PDG:
                            k2mu_Motherd2+=1
                            check,d2ka_veto = track_checks(index,d1Partkey,true_daught1,reco_part,d2ka_veto,0)
                            if check == -1: d2ka_veto[0] += 1 


                    if true_mother.GetPdgCode() == 211:
                        d1GrannyKey = true_mother.GetMotherId()
                        true_gran = sTree.MCTrack[d1GrannyKey]
                        pi2muEvents+=1
                        if true_gran.GetPdgCode() == HiddPart_PDG:
                            pi2mu_MotherHP+=1
                            check,HP2pi_veto = track_checks(index,d1Partkey,true_daught1,reco_part,HP2pi_veto,0)
                            if check == -1: HP2pi_veto[0] += 1 
                        if true_gran.GetPdgCode() == daught2_PDG:
                            pi2mu_Motherd2+=1
                            check,d2pi_veto = track_checks(index,d1Partkey,true_daught1,reco_part,d2pi_veto,0)
                            if check == -1: d2pi_veto[0] += 1

                    if true_mother.GetPdgCode() == HiddPart_PDG:              # checks mother is hidden particle
                        HiddPartDaught1Events+=1
                        event = True
                        for index2,reco_part2 in enumerate(sTree.FitTracks):  # loops over index and data of track particles
                            p2Partkey = sTree.fitTrack2MC[index2]                 # matches track to MC particle key
                            true_part2 = sTree.MCTrack[p2Partkey]                 # gives MC particle data
                            if abs(true_part2.GetPdgCode()) == 321:               # checks particle is kaon
                                part2MotherKey = true_part2.GetMotherId()            # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[part2MotherKey]          # obtains mother particle data
                                    
                                if (part2MotherKey==daught1MotherKey and daught2=='K+/-') or (true_mother.GetPdgCode() == daught2_PDG and daught2!='K+/-'):                 # check if keys are the same
                                    daught2toKaonEvents+=1
                                    #print(reco_part2)
                                    ####################
                                    #####  CHECKS  #####
                                    fit_status1 = reco_part.getFitStatus() 
                                    fit_status2 = reco_part2.getFitStatus()
                                    if fit_status1.isFitConverged() and fit_status2.isFitConverged():
                                        veto[0] += 1
                                    else:
                                        #print('At least one of the track fits did not converge')
                                        break
                                    
                                    check,veto = track_checks(index,d1Partkey,true_daught1,reco_part,veto,1)   # performs various track checks
                                    if check == -1: event = False

                                    if event == True:   # if the first track was fine, check the other one (ensures event veto counting is correct)
                                        check,veto = track_checks(index2,p2Partkey,true_part2,reco_part2,veto,1)
                                        if check == -1: event = False
                                    
                                    HiddPart_Pos = ROOT.TLorentzVector()
                                    daught1_LVec = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                    part2_LVec = ROOT.TLorentzVector()                   # declares variable as TLorentzVector class
                                    HiddPart_LVec = ROOT.TLorentzVector()                # declares variable as TLorentzVector class
                                    HiddPart_LVec,daught1_LVec,part2_LVec,X,Y,Z,doca = RedoVertexing(index,index2) # uses RedoVertexing to iterate track fitting
                                    if HiddPart_LVec == -1:
                                        print('RedoVertexing extrapolation failed')
                                        break

                                    h['recovertex'].Fill(Z)
                                    if not isInFiducial(X,Y,Z):
                                        #print('Neutralino decayed outside fiducial volume')
                                        veto[6] += 1
                                        check = -1
                                    
                                    h['doca'].Fill(doca)
                                    if not doca < docaCut: 
                                        #print('distance of closest approach too large')
                                        veto[7] +=1
                                        event = False
                                        
                                    T = true_mother.GetStartT()
                                    HiddPart_Pos.SetXYZT(X,Y,Z,T)
                                    tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
                                    ip = ImpactParameter(tr,HiddPart_Pos,HiddPart_LVec)
                                    h['IP_target'].Fill(ip)
                                    if not ip < ipCut:
                                        #print('neutralino impact parameter to target too large')
                                        veto[8] += 1
                                        event = False

                                    if event == False: break
                                    d2toKaonEventsAfterChecks+=1

                                    #check2,d2toKaonveto = track_checks(index2,p2Partkey,true_part2,reco_part2,d2toKaonveto,0)
                                    #if not check2 == 0:   # performs various checks (i.e. vertex position, fiducial volume,...)
                                    #    continue

                                    #######################
                                    ### Reco Properties ###
                                    HiddPart_mass = HiddPart_LVec.M()
                                    HiddPart_recoMom = HiddPart_LVec.P()
                                    p2M = part2_LVec.M()
                                    p2P = part2_LVec.P()
                                    p2E = part2_LVec.E()
                                    d1M = daught1_LVec.M()
                                    d1P = daught1_LVec.P()
                                    d1E = daught1_LVec.E()
                                    #######################
                                    ### True Properties ###
                                    p2MotherTrueMass = true_mother.GetMass()
                                    p2MotherTrueMom = true_mother.GetP()             
                                    part2TrueMom = true_part2.GetP()
                                    part2TrueMass = true_part2.GetMass()
                                    daught1TrueMom = true_daught1.GetP()
                                    daught1TrueMass = true_daught1.GetMass()

                                    mom_diff = p2MotherTrueMom - HiddPart_recoMom

                                    part2='K+/-'
                                    h[daught1 + 'TrueMom'].Fill(daught1TrueMom)
                                    h[daught1 + 'TrueMass'].Fill(daught1TrueMass)
                                    h[part2 + 'TrueMom'].Fill(part2TrueMom)
                                    h[part2 + 'TrueMass'].Fill(part2TrueMass)
                                    h[HiddPart + 'TrueMass'].Fill(p2MotherTrueMass)            
                                    h[HiddPart + 'TrueMom'].Fill(p2MotherTrueMom)
                                    h[HiddPart + 'RecoMass'].Fill(HiddPart_mass)                        
                                    h[HiddPart + 'RecoMom'].Fill(HiddPart_recoMom)                
                                    #h['Chi2'].Fill(d1_chi2)       
                                    #h['Chi2'].Fill(p2_chi2)                             
                                    h[HiddPart + 'DeltaMom'].Fill(mom_diff)
                                    h[part2 + 'RecoMom'].Fill(p2P)
                                    h[daught1 + 'RecoMom'].Fill(d1P)
                                    h[part2 + 'RecoMass'].Fill(p2M)
                                    h[daught1 + 'RecoMass'].Fill(d1M)

                                    successful_events.append(n)   # adds entries to the list
                                    m = successful_events[0]   # arbitrarily picks the first one as an example

                                    time_res(d1P,daught1,d1Partkey, daught1_PDG,n,m)
                                    time_res(p2P,part2,p2Partkey, 321, n, m)
                                                
                            if abs(true_part2.GetPdgCode()) == 211:     # checks particle is kaon
                                part3MotherKey = true_part2.GetMotherId()          # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[part3MotherKey]          # obtains mother particle data
                                    
                                if true_mother.GetPdgCode() == daught2_PDG and daught2!='K+/-':                 # check if keys are the same
                                    daught2toPionEvents+=1
                                    #print(reco_part2)
                                    p3MotherTrueMass = true_mother.GetMass()               # get Neutralino/final states mother mass
                                    p3MotherTrueMom = true_mother.GetP()                   # get Neutralino/final states mother mom
                                    check3,p3_chi2,d2toPionveto = track_checks(index2,true_part2,reco_part2,d2toPionveto)
                                    if not check3 == 0:   # performs various checks (i.e. vertex position, fiducial volume,...)
                                        continue
                                    d2toPionEventsAfterChecks+=1
                                    
                                    part3_LVec = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                    HiddPart_LVec = ROOT.TLorentzVector()                # declares variable as TLorentzVector class
                                    HiddPart_LVec,daught1_LVec,part3_LVec,doca = RedoVertexing(index,index2) # uses RedoVertexing to iterate track fitting
                                    if HiddPart_LVec == -1: continue
                                    if doca > 2.: 
                                        #print('distance of closest approach too large')
                                        continue

                                    HiddPart_mass = HiddPart_LVec.M()
                                    HiddPart_recoMom = HiddPart_LVec.P()
                                    mom_diff = p3MotherTrueMom - HiddPart_recoMom

                                    p3M = part3_LVec.M()
                                    p3P = part3_LVec.P()
                                    p3E = part3_LVec.E()
                                   
                                    
                                    part3TrueMom = true_part2.GetP()
                                    part3TrueMass = true_part2.GetMass()
                                    
                                    part3='pi+/-'
                                    h[part3 + 'TrueMom'].Fill(part3TrueMom)
                                    h[part3 + 'TrueMass'].Fill(part3TrueMass)
                                    h[HiddPart + 'TrueMass'].Fill(p3MotherTrueMass)            
                                    h[HiddPart + 'TrueMom'].Fill(p3MotherTrueMom)
                                    h[HiddPart + 'RecoMass'].Fill(HiddPart_mass)                        
                                    h[HiddPart + 'RecoMom'].Fill(HiddPart_recoMom)                
                                    h['Chi2'].Fill(p3_chi2)                             
                                    h[HiddPart + 'DeltaMom'].Fill(mom_diff)
                                    h[part3 + 'RecoMom'].Fill(p3P)
                                    h[part3 + 'RecoMass'].Fill(p3M)


                                    time_res(p3P,part3, p2Partkey, 211, n,m)

        print(2*' ' + 80*'_')
        print(2*' ' + '||  Total number of events = ' + str(totalEvents))
        print(2*' ' + '||  ' + daught1 + ' number of events = ' + str(daught1Events))
        print(2*' ' + '||  ' + daught1 + ' number of events after checks = ' + str(d1EventsAfterChecks))
        print(2*' ' + '||  ==>' + daught1 + ' number of events isInFiducial = ' + str(d1veto[0]))
        print(2*' ' + '||  ==>' + daught1 + ' number of events checkFiducialVolume = ' + str(d1veto[1]))
        print(2*' ' + '||  ==>' + daught1 + ' number of events fit_status = ' + str(d1veto[2]))
        print(2*' ' + '||  ==>' + daught1 + ' number of events fit_nmeas = ' + str(d1veto[3]))
        print(2*' ' + '||  ==>' + daught1 + ' number of events fit_chi2 = ' + str(d1veto[4]))
        print(2*' ' + '||  ' + daught1 + ' number of events (with K+/- mother) = ' + str(k2muEvents))
        print(2*' ' + '||  ' + daught1 + ' number of events (with K+/- mother and gran '+ HiddPart +') = ' + str(k2mu_MotherHP))
        print(2*' ' + '||  ' + daught1 + ' number of events (with K+/- mother and gran '+ daught2 +') = ' + str(k2mu_Motherd2))
        print(2*' ' + '||  ' + daught1 + ' number of events (with pi+/- mother) = ' + str(pi2muEvents))
        print(2*' ' + '||  ' + daught1 + ' number of events (with pi+/- mother and gran '+ HiddPart +') = ' + str(pi2mu_MotherHP))
        print(2*' ' + '||  ' + daught1 + ' number of events (with pi+/- mother and gran '+ daught2 +') = ' + str(pi2mu_Motherd2))
        print(2*' ' + '||  ' + daught1 + ' number of events (with ' + HiddPart + ' mother) = ' + str(HiddPartDaught1Events))
        print(2*' ' + '||  K+/- number of events (with ' + daught2 + ' mother) = ' + str(daught2toKaonEvents))
        print(2*' ' + '||  K+/- number of events (with ' + daught2 + ' mother) after checks = ' + str(d2toKaonEventsAfterChecks))
        #print(2*' ' + '||  ==> K+/- number of events (with ' + daught2 + ' mother) isInFiducial = ' + str(d2toKaonveto[0]))
        #print(2*' ' + '||  ==> K+/- number of events (with ' + daught2 + ' mother) checkFiducialVolume = ' + str(d2toKaonveto[1]))
        #print(2*' ' + '||  ==> K+/- number of events (with ' + daught2 + ' mother) fit_status = ' + str(d2toKaonveto[2]))
        #print(2*' ' + '||  ==> K+/- number of events (with ' + daught2 + ' mother) fit_nmeas = ' + str(d2toKaonveto[3]))
        #print(2*' ' + '||  ==> K+/- number of events (with ' + daught2 + ' mother) fit_chi2 = ' + str(d2toKaonveto[4]))
        print(2*' ' + '||  pi+/- number of events (with ' + daught2 + ' mother) = ' + str(daught2toPionEvents))
        print(2*' ' + '||  pi+/- number of events (with ' + daught2 + ' mother) after checks = ' + str(d2toPionEventsAfterChecks))
        print(2*' ' + '||  ==> pi+/- number of events (with ' + daught2 + ' mother) isInFiducial = ' + str(d2toPionveto[0]))
        print(2*' ' + '||  ==> pi+/- number of events (with ' + daught2 + ' mother) checkFiducialVolume = ' + str(d2toPionveto[1]))
        print(2*' ' + '||  ==> pi+/- number of events (with ' + daught2 + ' mother) fit_status = ' + str(d2toPionveto[2]))
        print(2*' ' + '||  ==> pi+/- number of events (with ' + daught2 + ' mother) fit_nmeas = ' + str(d2toPionveto[3]))
        print(2*' ' + '||  ==> pi+/- number of events (with ' + daught2 + ' mother) fit_chi2 = ' + str(d2toPionveto[4]))
        print(2*' ' + 80*'_'+ '\n')

        h[daught1 + 'ProbMeasr'] = createRatio(h[daught1 + 'SmearedMass'],h['TotalSmearedMass'],daught1 + 'ProbMeasr')
        h[part2 + 'ProbMeasr'] = createRatio(h[part2 + 'SmearedMass'],h['TotalSmearedMass'],part2 + 'ProbMeasr')
        if not daught2=='K+/-':
            h[part3 + 'ProbMeasr'] = createRatio(h[part3 + 'SmearedMass'],h['TotalSmearedMass'],part3 + 'ProbMeasr')
graph = {}
def makePlots2(HiddPart,part1,part2,part3):
    
    key='DAUGHTERS'
    title='Time and velocity plots'
    
    ut.bookCanvas(h,key + '_TV',title,nx=1300,ny=800,cx=3,cy=2)
    cv = h[key + '_TV'].cd(1)
    h['StrawTime'].Add(h[part1 + 'StrawTime'])
    h['StrawTime'].Add(h[part2 + 'StrawTime'])
    if not part3==None:
        h['StrawTime'].Add(h[part3 + 'StrawTime'])
    h['StrawTime'].Draw("nostack")

    cv = h[key + '_TV'].cd(2)
    h['EcalTime'].Add(h[part1 + 'EcalTime'])
    h['EcalTime'].Add(h[part2 + 'EcalTime'])
    if not part3==None:
        h['EcalTime'].Add(h[part3 + 'EcalTime'])
    h['EcalTime'].Draw("nostack")
    

    cv = h[key + '_TV'].cd(3)
    h['DirDeltaTimeSmeared'].Add(h[part1 + 'DirDeltaTimeSmeared'])
    h['DirDeltaTimeSmeared'].Add(h[part2 + 'DirDeltaTimeSmeared'])
    if not part3==None:
        h['DirDeltaTimeSmeared'].Add(h[part3 + 'DirDeltaTimeSmeared'])
    h['DirDeltaTimeSmeared'].Draw("nostack")

    cv = h[key + '_TV'].cd(4)
    h['FlightLen'].Add(h[part1 + 'FlightLen'])
    h['FlightLen'].Add(h[part2 + 'FlightLen'])
    if not part3==None:
        h['FlightLen'].Add(h[part3 + 'FlightLen'])
    h['FlightLen'].Draw("nostack")

    cv = h[key + '_TV'].cd(5)
    h['Speed'].Add(h[part1 + 'Speed'])
    h['Speed'].Add(h[part2 + 'Speed'])
    if not part3==None:
        h['Speed'].Add(h[part3 + 'Speed'])
    h['Speed'].Draw("nostack")

    #h[key + '_TV'].Print('DaughterTVProp'+ currentDate + '.png')

    title='Momenta and mass plots'
    ut.bookCanvas(h,key + '_MOM', title , nx=1300, ny=800, cx=3, cy=2)
    cv = h[key + '_MOM'].cd(1)
    h['StrawMom'].Add(h[part1 + 'StrawMom'])
    h['StrawMom'].Add(h[part2 + 'StrawMom'])
    if not part3==None:
        h['StrawMom'].Add(h[part3 + 'StrawMom'])
    h['StrawMom'].Draw("nostack")


    cv = h[key + '_MOM'].cd(2)
    h['EcalMom'].Add(h[part1 + 'EcalMom'])
    h['EcalMom'].Add(h[part2 + 'EcalMom'])
    if not part3==None:
        h['EcalMom'].Add(h[part3 + 'EcalMom'])
    h['EcalMom'].Draw("nostack")

    cv = h[key + '_MOM'].cd(3)
    h['RecoMom'].Add(h[part1 + 'RecoMom'])
    h['RecoMom'].Add(h[part2 + 'RecoMom'])
    if not part3==None:
        h['RecoMom'].Add(h[part3 + 'RecoMom'])
    h['RecoMom'].Draw("nostack")
    
    cv = h[key + '_MOM'].cd(4)
    h['DeltaMom'].Add(h[part1 + 'DeltaMom'])
    h['DeltaMom'].Add(h[part2 + 'DeltaMom'])
    if not part3==None:
        h['DeltaMom'].Add(h[part3 + 'DeltaMom'])
    h['DeltaMom'].Draw("nostack")

    cv = h[key + '_MOM'].cd(5)
    h['RecoMass'].Add(h[part1 + 'RecoMass'])
    h['RecoMass'].Add(h[part2 + 'RecoMass'])
    if not part3==None:
        h['RecoMass'].Add(h[part3 + 'RecoMass'])
    h['RecoMass'].Draw("nostack")

    cv = h[key + '_MOM'].cd(6)
    h[part1 + 'SmearedMass'].Fit("landau")
    h[part1 + 'SmearedMass'].GetFunction("landau").SetLineColor(kBlack)
    h[part2 + 'SmearedMass'].Fit('landau')
    h[part2 + 'SmearedMass'].GetFunction("landau").SetLineColor(kBlack)
    if not part3==None:
        h[part3 + 'SmearedMass'].Fit('landau')
        h[part3 + 'SmearedMass'].GetFunction("landau").SetLineColor(kBlack)

    h['SmearedMass'].Add(h[part1 + 'SmearedMass'])
    h['SmearedMass'].Add(h[part2 + 'SmearedMass'])
    if not part3==None:
        h['SmearedMass'].Add(h[part3 + 'SmearedMass'])
    h['SmearedMass'].Draw("nostack")

    #h['DAUGHTERS_MOM'].Print('DaughterPProp'+ currentDate + '.png')

    if part3==None:
        partString=''
    else:
        partString=' or pion'
    title='Probability Plots'
    ut.bookCanvas(h,key + '_PROB',title,nx=1300,ny=800,cx=3,cy=2)
    cv = h[key + '_PROB'].cd(1)
    h[part1 + 'ProbMeasr'].SetMarkerColor(38)
    polyFit1.SetLineColor(4)
    h[part1 + 'ProbMeasr'].Fit('polyFit1')
    h[part1 + 'ProbMeasr'].Draw('E2')
    h[part1 + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
    h[part1 + 'ProbMeasr'].SetYTitle('Prob( particle = muon )')
    h[part1 + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    cv = h[key + '_PROB'].cd(2)
    h[part2 + 'ProbMeasr'].SetMarkerColor(46)
    polyFit2.SetLineColor(2)
    h[part2 + 'ProbMeasr'].Fit('polyFit2')
    h[part2 + 'ProbMeasr'].Draw('E2')
    h[part2 + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
    h[part2 + 'ProbMeasr'].SetYTitle('Prob( particle = kaon )')
    h[part2 + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    if not part3==None:
        cv = h[key + '_PROB'].cd(3)
        h[part3 + 'ProbMeasr'].SetMarkerColor(8)
        polyFit2.SetLineColor(3)
        h[part3 + 'ProbMeasr'].Fit('polyFit2')
        h[part3 + 'ProbMeasr'].Draw('E2')
        h[part3 + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
        h[part3 + 'ProbMeasr'].SetYTitle('Prob( particle = pion )')
        h[part3 + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    cv = h[key + '_PROB'].cd(4)
    graph['partIDProb'] = TMultiGraph()
    graph['partIDProb'].SetName('Prob(correct ID paticle)')
    x1, y1 = array( 'd' ), array( 'd' )
    ex1, ey1 = array( 'd' ), array( 'd' )
    x2, y2 = array( 'd' ), array( 'd' )
    ex2, ey2 = array( 'd' ), array( 'd' )
    if not part3==None:
        x3, y3 = array( 'd' ), array( 'd' )
        ex3, ey3 = array( 'd' ), array( 'd' )
    i=0
    n=0
    numBins = h[part1 + 'ProbMeasr'].GetNbinsX()
    for i in range(numBins):
        x1.append(h[part1 + 'ProbMeasr'].GetBinCenter(i))
        ex1.append(0)
        y1.append(h[part1 + 'ProbMeasr'].GetBinContent(i))
        ey1.append(ROOT.TMath.Sqrt((h[part1 + 'ProbMeasr'].GetBinContent(i))/10))
        x2.append(h[part2 + 'ProbMeasr'].GetBinCenter(i))
        ex2.append(0)
        y2.append(h[part2 + 'ProbMeasr'].GetBinContent(i))
        ey2.append(ROOT.TMath.Sqrt((h[part1 + 'ProbMeasr'].GetBinContent(i))/10))
        if not part3==None:
            x3.append(h[part3 + 'ProbMeasr'].GetBinCenter(i))
            ex3.append(0)
            y3.append(h[part3 + 'ProbMeasr'].GetBinContent(i))
            ey3.append(ROOT.TMath.Sqrt((h[part1 + 'ProbMeasr'].GetBinContent(i))/10))
        n=n+1
    graph["1"] = TGraphErrors( n, x1, y1, ex1, ey1 )
    graph["1"].SetName('Prob(ID = ' + part1 + ')')
    graph["1"].SetTitle('Prob(ID = ' + part1 + ')')
    graph["1"].GetYaxis().SetTitle( 'Prob(particle = muon)' )
    graph["1"].SetLineColor( 4 )
    graph["1"].SetLineWidth( 1 )
    graph["1"].SetMarkerColor( 4 )
    graph["1"].SetMarkerStyle( 20 )
    graph["1"].SetMarkerSize(0.5)   
    graph["2"] = TGraphErrors( n, x2, y2, ex2, ey2 )
    graph["2"].SetName('Prob(ID = ' + part2 + ')')
    graph["2"].SetTitle('Prob(ID = ' + part2 + ')')
    graph["2"].GetYaxis().SetTitle( 'Prob(particle = kaon)' )
    graph["2"].SetLineColor( 2 )
    graph["2"].SetLineWidth( 1 )
    graph["2"].SetMarkerColor( 2 )
    graph["2"].SetMarkerStyle( 20 )
    graph["2"].SetMarkerSize(0.5)
    graph['partIDProb'].Add(graph["1"], "PC")
    graph['partIDProb'].Add(graph["2"], "PC")
    if not part3==None:
        graph["3"] = TGraphErrors( n, x3, y3, ex3, ey3 )
        graph["3"].SetName('Prob(ID = ' + part3 + ')')
        graph["3"].SetTitle('Prob(ID = ' + part3 + ')') 
        graph["3"].GetYaxis().SetTitle( 'Prob(particle = pion)' )
        graph["3"].SetLineColor( 3 )
        graph["3"].SetLineWidth( 1 )
        graph["3"].SetMarkerColor( 3 )
        graph["3"].SetMarkerStyle( 20 )
        graph["3"].SetMarkerSize(0.5)
        graph['partIDProb'].Add(graph["3"], "PC")
    graph['partIDProb'].Draw("A pfc plc")#P PLC PFCPLC PFC
    for ckey in graph:
        graph[ckey].GetXaxis().SetTitle( 'Mass [GeV/c2]' )
    graph['partIDProb'].GetYaxis().SetTitle( 'Prob(particle=(kaon or muon' + partString + '))' )
    graph['partIDProb'].GetYaxis().SetTitleOffset(1.5)
    graph['partIDProb'].GetXaxis().SetRangeUser(0,1.5)
    gPad.BuildLegend()
    h[key + '_PROB'].Print('DaughterProb'+ currentDate + '.png')
    
def print_menu(): 
    print ('\n \n' + 30 * "-" + "MENU" + 30 * "-")
    print ("1. RPV SUSY Benchmark1 --> K+/- mu+/- final state")
    print ("2. RPV SUSY Benchmark1 --> K*+/- mu+/- final state")
    print ("3. RPV SUSY Benchmark1 --> K*0 nu_mu final state")
    print ("4. Exit")
    print (64 * "-"+ '\n')

########################################
     
while loop:        
    print_menu()
    choice = input("Enter your choice [1-4]: ")
     
    if choice==1:     
        print ("RPV SUSY Benchmark1 --> K+/- mu+/- final state selected.")
        HiddenPart = 'Neutralino'
        daughter1 = 'Mu'
        daughter2 = 'K+/-'
        create_Hists(HiddenPart, daughter1, daughter2, None)
        finState2t1t2(HiddenPart,daughter1,daughter2)
        makePlots2(HiddenPart,daughter1,daughter2,None)
        loop=False
    elif choice==2:
        print ("RPV SUSY Benchmark1 --> K*+/- mu+/- final state selected.")
        HiddenPart = 'Neutralino'
        daughter1 = 'Mu'
        daughter2 = 'K*+/-'        
        create_Hists(HiddenPart, daughter1, 'K+/-', 'pi+/-')
        finState2t1t2(HiddenPart,daughter1,daughter2)
        makePlots2(HiddenPart,daughter1,'K+/-', 'pi+/-')
        loop=False
    elif choice==3:
        print ("RPV SUSY Benchmark1 --> K*0 nu_mu final state selected")
        HiddenPart = 'Neutralino'
        daughter1 = 'Mu'
        daughter2 = 'K*0'
        create_Hists(HiddenPart, daughter1, 'K+/-','pi+/-')
        finState2t1t2(HiddenPart,daughter1,daughter2)
        makePlots2(HiddenPart,daughter1,'K+/-', 'pi+/-')
        loop=False
    elif choice==4:
        print ("EXIT was selected.")
        loop=False 
    else:
        print("No such option. Try again.")

hfile = inputFile.split(',')[0].replace('_rec','_RPVeditana')#create outputFile
if hfile[0:4] == "/eos" or not inputFile.find(',')<0:
    # do not write to eos, write to local directory 
    tmp = hfile.split('/')
    hfile = tmp[len(tmp)-1]                                   #occurs only for cern users
ROOT.gROOT.cd()
f = TFile(hfile,"RECREATE")
for akey in h:
    cln = h[akey].Class().GetName()
    if not cln.find('TH')<0 or not cln.find('TP')<0:   
        h[akey].Write()
for bkey in graph:
    cln = graph[bkey].Class().GetName()
    if not cln.find('TG')<0 or not cln.find('TM')<0:   
        f.WriteTObject(graph[bkey])
f.Close()

