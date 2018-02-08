# #########################################
#         RPV_EDITA Test Code             #
###########################################
#from __future__ import print_function
import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
from array import array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TGraph, TF1, TMultiGraph
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
measCutFK = 25
measCutPR = 22
docaCut = 2.
c = 2.99792458*(10**8)
currentDate = datetime.datetime.now().strftime("%y_%m_%d_%H%M")
polyFit1 = TF1("polyFit1","pol9")
polyFit2 = TF1("polyFit2","pol9")

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
    print('Got cbmsim')

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

# calculate z front face of ecal, needed later
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
def create_Hists(HiddPart,daught1,daught2):
    ###########################
    ####  Muon Histograms  ####
    ut.bookHist(h,daught1 + 'StrawTime','Gaussian Straw t measurement',600,321.7,324.7)
    ut.bookHist(h,daught1 + 'EcalTime','Gaussian Ecal t measurement',600,359.7,361.7)
    ut.bookHist(h,daught1 + 'DirDeltaTime','Straw-ECAL Time of Flight (directly)',600,37.5,38.9)      # time of flight
    ut.bookHist(h,daught1 + 'FlightLen','Straw-ECAL Flight Lenght',600,11.36,11.47)                   # flight Length
    ut.bookHist(h,daught1 + 'Speed','Beta value',600,0.993,1.0014)                                    # speed
    ut.bookHist(h,daught1 + 'StrawMom','Straw Momentum',600,-0.05,120.)                               # straw momentum
    ut.bookHist(h,daught1 + 'EcalMom','Ecal Momentum',600,-0.05,120.)                                 # ecal  momentum
    ut.bookHist(h,daught1 + 'DeltaMom','Straw-Ecal Momentum',600,-0.1,0.21)                           # delta momentum
    ut.bookHist(h,daught1 + 'RecoMom','Reco Momentum',600,-0.05,120.)                                 # reco  momentum
    ut.bookHist(h,daught1 + 'TrueMom','True Momentum',600,-0.05,120.)                                 # true  momentum
    ut.bookHist(h,daught1 + 'RecoMass','Reco Mass',400,0.,2.)                                         # reco  mass
    ut.bookHist(h,daught1 + 'TrueMass','True Mass',400,0.,2.)                                         # true  mass
    ut.bookHist(h,daught1 + 'SmearedMass','Smeared Mass',400,0.,2.)                                   # smrd  mass
    ut.bookHist(h,daught1 + 'ProbMeasr','Probs identifying Muon',400,0.,2.)                           # ID Prob
    ###########################
    ####  Kaon Histograms  ####
    ut.bookHist(h,daught2 + 'StrawTime','Gaussian Straw t measurement',600,321.7,324.7)    
    ut.bookHist(h,daught2 + 'EcalTime','Gaussian Ecal t measurement',600,359.7,361.7)    
    ut.bookHist(h,daught2 + 'DirDeltaTime','Straw-ECAL Time of Flight (directly)',600,37.5,38.9)      # time of flight
    ut.bookHist(h,daught2 + 'FlightLen','Straw-ECAL Flight Length',600,11.36,11.47)                   # flight Length
    ut.bookHist(h,daught2 + 'Speed','Beta value',600,0.993,1.0014)                                    # speed
    ut.bookHist(h,daught2 + 'StrawMom','Straw Momentum',600,-0.05,120.)                               # straw momentum
    ut.bookHist(h,daught2 + 'EcalMom','Ecal Momentum',600,-0.05,120.)                                 # ecal  momentum
    ut.bookHist(h,daught2 + 'DeltaMom','Straw-Ecal Momentum',600,-0.1,0.21)                           # delta momentum
    ut.bookHist(h,daught2 + 'RecoMom','Reco Momentum',600,-0.05,120.)                                 # reco  momentum
    ut.bookHist(h,daught2 + 'TrueMom','True Momentum',600,-0.05,120.)                                 # true  momentum
    ut.bookHist(h,daught2 + 'RecoMass','Reco Mass',400,0.,2.)                                         # reco  mass
    ut.bookHist(h,daught2 + 'TrueMass','True Mass',400,0.,2.)                                         # true  mass
    ut.bookHist(h,daught2 + 'SmearedMass','Smeared Mass',400,0.,2.)                                   # smrd  mass
    ut.bookHist(h,daught2 + 'ProbMeasr','Prob identifying Kaon',400,0.,2.)                            # ID Prob
    ###########################
    ut.bookHist(h,'TotalSmearedMass','Smeared Mass',400,0.,2.)                                   # Total mass
    #################################
    ####  Neutralino Histograms  ####
    ut.bookHist(h,HiddPart + 'TrueMass','Monte Carlo Mass',500,0.,2.)                                # true mass
    ut.bookHist(h,HiddPart + 'RecoMass','Reconstructed Mass',500,0.,2.)                              # reco mass
    #ut.bookHist(h,HiddPart + '_no_iter','Reconstructed Mass (without track iterations)',500,0.,2.)   # reco mass(without track itrns)
    ut.bookHist(h,HiddPart + 'TrueMom','True (red) & Reco. (blue) Momentum',100,0.,300.)             # true momentum 
    ut.bookHist(h,HiddPart + 'RecoMom','Reconstructed Momentum',100,0.,300.)                         # reco momentum
    ut.bookHist(h,HiddPart + 'DeltaMom','True/Reco Momentum Difference',100,-3.,3)                   # true-reco momentum difference

    ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared',100,0.,3.)                                      # chi squared track fitting
    #ut.bookHist(h,'normdistr','Gaussian Distribution',500,-0.05,0.05)                               #
    #ut.bookHist(h,'smearedmass1','Time Smeared Neutralino Mass',500,0.,2.)
    #ut.bookHist(h,'smearedmass2','Time Smeared Neutralino Mass',500,0.,2.)
    #ut.bookHist(h,'smearedP1','Time Smeared Neutralino Momentum P1(red) P2(blue)',500,0.,300.)
    #ut.bookHist(h,'smearedP2','Time Smeared Neutralino Momentum',500,0.,300.)

    
    print("Created Histograms")

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
     if not rc: return -1,-1,-1,doca # extrapolation failed, makes no sense to continue
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
     return NeutralinoLV,LV[t1],LV[t2],doca

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

nEvents = min(sTree.GetEntries(),nEvents)

def time_res(partkey):
    smearStrawTime = None
    smearEcalTime = None
    deltaT = None
    r = None
    v = None
    strawP = None
    ecalP = None
    if sTree.GetBranch("strawtubesPoint"):
        x_array = []
        y_array = []
        z_array = []
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
        
        min_z_index = z_array.index(min(z_array))
        straw_z = 0.01*z_array[min_z_index]
        straw_x = 0.01*x_array[min_z_index]
        straw_y = 0.01*y_array[min_z_index]
        straw_time = t_array[min_z_index]
        if straw_time != None:
            smearStrawTime = np.random.normal(loc=straw_time,scale=0.01,size=None) # current width of 10 ps
            #smearStrawTime=straw_time 

        strawPx = px_array[min_z_index]
        strawPy = py_array[min_z_index]
        strawPz = pz_array[min_z_index]
        strawP = ROOT.TMath.Sqrt((strawPx**2) + (strawPy**2) + (strawPz**2)) # straw tube momentum

                

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
                        ecalP = ROOT.TMath.Sqrt((ecal_px**2) + (ecal_py**2) + (ecal_pz**2)) # straw tube momentum

                        if ecal_time != None:
                            smearEcalTime = np.random.normal(loc=ecal_time,scale=0.01,size=None) # current width of 10 ps
                            #smearEcalTime=ecal_time

                        if not ecal_time <= straw_time:
                            deltaT=abs(smearStrawTime - smearEcalTime)
                            r = ROOT.TMath.Sqrt(((ecal_x - straw_x)**2) + ((ecal_y - straw_y)**2) + ((ecal_z - straw_z)**2))
                            v=((r/deltaT)*(10**9) )# units of nanoseconds
                            
    return smearStrawTime,smearEcalTime,deltaT,r,v,strawP,ecalP

def createRatio(h1, h2):
    h3 = h1.Clone("h3")
    h3.SetMarkerStyle(20)
    h3.SetMarkerSize(0.7)
    h3.SetTitle("")
    h3.SetMinimum(0.8)
    h3.SetMaximum(1.35)
    # Set up plot for markers and errors
    h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetRangeUser(-0.1,1.2)
    y.SetTitleOffset(1.)
    ## Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetRangeUser(0,1.5)
    return h3

def finState2t1t2(HiddPart,daught1,daught2):
    if daught1 == 'Mu' and daught2 ==  'K+/-':
        daught1_PDG = 13
        daught2_PDG = 321
        HiddPart_PDG = 9900015
    elif daught2 == 'K*+/-':
        daught2_PDG=211     #values are wrong here
    elif daught2 == 'K*0':
        daught2_PDG = 311  #values are wrong here
    if sTree.GetBranch("FitTracks"):
        d2_decaycheck = 0
        for n in range(nEvents):                            # loop over events
            rc = sTree.GetEntry(n)                              # load tree entry
            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles                                   
                d1Partkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_daught1 = sTree.MCTrack[d1Partkey]                  # gives MC particle data
                #print(reco_part)
                if abs(true_daught1.GetPdgCode()) == daught1_PDG:               # checks particle is muon
                    daught1MotherKey = true_daught1.GetMotherId()             # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[daught1MotherKey]          # obtains mother particle data
                    d1MotherTrueMass = true_mother.GetMass()         # get Neutralino/final states mother mass
                    d1MotherTrueMom = true_mother.GetP()             # get Neutralino/final states mother mom
                    if true_mother.GetPdgCode() == 321:
                        #print('Kaon has decayed to a muon')
                        d2_decaycheck+=1
                    if true_mother.GetPdgCode() == HiddPart_PDG:             # checks mother is Neutralino

                        Decay_X = true_daught1.GetStartX()
                        Decay_Y = true_daught1.GetStartY()
                        Decay_Z = true_daught1.GetStartZ()
                        
                        if not isInFiducial(Decay_X,Decay_Y,Decay_Z):
                            #print('Neutralino decayed outside fiducial volume')
                            continue
                        if not checkFiducialVolume(sTree,index,dy): 
                            #print('Track outside fiducial volume')
                            continue 
                        d1_status = reco_part.getFitStatus()             
                        if not d1_status.isFitConverged():
                            #print('Fit did not converge')
                            continue
                        d1_nmeas = d1_status.getNdf()                      
                        if not d1_nmeas > 25:
                            #print('Too few measurements')
                            continue

                        d1_rchi2 = d1_status.getChi2()                      # gets chi squared value
                        d1_chi2 = (d1_rchi2/d1_nmeas)                       # gets chi value
                        if not d1_chi2 < 4:
                            #print('Chi squared value too high')
                            continue

                        for index2,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles
                            d2Partkey = sTree.fitTrack2MC[index2]                 # matches track to MC particle key
                            true_daught2 = sTree.MCTrack[d2Partkey]                  # gives MC particle data
                            if abs(true_daught2.GetPdgCode()) == daught2_PDG:               # checks particle is kaon
                                daught2MotherKey = true_daught2.GetMotherId()             # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[daught2MotherKey]          # obtains mother particle data
                                if daught2MotherKey==daught1MotherKey:                    # check if keys are the same
                                    d2MotherTrueMass = true_mother.GetMass()        # get Neutralino/final states mother mass
                                    d2MotherTrueMom = true_mother.GetP()            # get Neutralino/final states mother mom

                                    if not checkFiducialVolume(sTree,index,dy): 
                                        #print('Decay outside fiducial volume')
                                        continue 
                                    d2_status = reco_part.getFitStatus()                
                                    if not d2_status.isFitConverged():
                                        #print('Fit did not converge')
                                        continue
                                    d2_nmeas = d2_status.getNdf() 
                                    if not d2_nmeas > 25:
                                        #print('Too few measurements')
                                        continue

                                    d2_rchi2 = d2_status.getChi2()                       # chi squared value
                                    d2_chi2 = (d2_rchi2/d2_nmeas)                         # gets chi value
                                    if not d2_chi2 < 4:
                                        #print('Chi squared value too high')
                                        continue
                                    
                                    daught1_LVec = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                    daught2_LVec = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                    HiddPart_LVec = ROOT.TLorentzVector()               # declares variable as TLorentzVector class
                                    HiddPart_LVec,daught1_LVec,daught2_LVec,doca = RedoVertexing(index,index2) # uses RedoVertexing to iterate track fitting
                                    if HiddPart_LVec == -1: continue
                                    if doca > 2.: 
                                        #print('distance of closest approach too large')
                                        continue

                                    HiddPart_mass = HiddPart_LVec.M()
                                    HiddPart_recoMom = HiddPart_LVec.P()
                                    mom_diff = d2MotherTrueMom - HiddPart_recoMom

                                    d2M = daught2_LVec.M()
                                    d2P = daught2_LVec.P()
                                    d2E = daught2_LVec.E()
                                    d1M = daught1_LVec.M()
                                    d1P = daught1_LVec.P()
                                    d1E = daught1_LVec.E()
                                    
                                    daught2TrueMom = true_daught2.GetP()
                                    daught2TrueMass = true_daught2.GetMass()
                                    daught1TrueMom = true_daught1.GetP()
                                    daught1TrueMass = true_daught1.GetMass()
                                    h[daught1 + 'TrueMom'].Fill(daught1TrueMom)
                                    h[daught1 + 'TrueMass'].Fill(daught1TrueMass)
                                    h[daught2 + 'TrueMom'].Fill(daught2TrueMom)
                                    h[daught2 + 'TrueMass'].Fill(daught2TrueMass)
                                    h[HiddPart + 'TrueMass'].Fill(d2MotherTrueMass)            
                                    h[HiddPart + 'TrueMom'].Fill(d2MotherTrueMom)
                                    h[HiddPart + 'RecoMass'].Fill(HiddPart_mass)                        
                                    h[HiddPart + 'RecoMom'].Fill(HiddPart_recoMom)                
                                    h['Chi2'].Fill(d1_chi2)       
                                    h['Chi2'].Fill(d2_chi2)                             
                                    h[HiddPart + 'DeltaMom'].Fill(mom_diff)
                                    h[daught2 + 'RecoMom'].Fill(d2P)
                                    h[daught1 + 'RecoMom'].Fill(d1P)
                                    h[daught2 + 'RecoMass'].Fill(d2M)
                                    h[daught1 + 'RecoMass'].Fill(d1M)

                                    d1_strwT,d1_ecalT,d1_Dt,d1_Len,d1_v,d1_strawP,d1_ecalP = time_res(d1Partkey) 
                                    if d1_strwT!= None and d1_ecalT!= None and d1_Dt!= None and d1_Len!= None and d1_v != None and d1_strawP!=None and d1_ecalP!=None:
                                        h[daught1 + 'StrawTime'].Fill(d1_strwT)
                                        h[daught1 + 'EcalTime'].Fill(d1_ecalT)
                                        h[daught1 + 'DirDeltaTime'].Fill(d1_Dt)
                                        h[daught1 + 'FlightLen'].Fill(d1_Len)
                                        d1_beta = d1_v/c
                                        h[daught1 + 'Speed'].Fill(d1_beta)
                                        h[daught1 + 'StrawMom'].Fill(d1_strawP)
                                        h[daught1 + 'EcalMom'].Fill(d1_ecalP)
                                        d1_Dp = d1_strawP-d1_ecalP
                                        h[daught1 + 'DeltaMom'].Fill(d1_Dp)

                                        d2_strwT,d2_ecalT,d2_Dt,d2_Len,d2_v,d2_strawP,d2_ecalP = time_res(d2Partkey)
                                        if d2_strwT!= None and d2_ecalT!= None and d2_Dt!= None and d2_Len!= None and d2_v != None and d2_strawP!=None and d2_ecalP!=None:     
                                            h[daught2 + 'StrawTime'].Fill(d2_strwT)
                                            h[daught2 + 'EcalTime'].Fill(d2_ecalT)
                                            h[daught2 + 'DirDeltaTime'].Fill(d2_Dt)
                                            h[daught2 + 'FlightLen'].Fill(d2_Len)
                                            d2_beta = d2_v/c
                                            h[daught2 + 'Speed'].Fill(d2_beta)
                                            h[daught2 + 'StrawMom'].Fill(d2_strawP)
                                            h[daught2 + 'EcalMom'].Fill(d2_ecalP)
                                            d2_Dp = d2_strawP-d2_ecalP
                                            h[daught2 + 'DeltaMom'].Fill(d2_Dp)
                                              
                                            
                                            if d2_beta < 1:
                                                #TRYING SOMETHING ELSE
                                                d2_smearedM = d2P*(ROOT.TMath.Sqrt(1-(d2_beta**2)))/d2_beta
                                                h[daught2 + 'SmearedMass'].Fill(d2_smearedM)
                                                h['TotalSmearedMass'].Fill(d2_smearedM)                                                
                                                d1_smearedM = d1P*(ROOT.TMath.Sqrt(1-(d1_beta**2)))/d1_beta
                                                h[daught1 + 'SmearedMass'].Fill(d1_smearedM)
                                                h['TotalSmearedMass'].Fill(d1_smearedM)  

        print('\n'+str(d2_decaycheck) + ' K+ --> mu decays before detection\n')
        h[daught1 + 'ProbMeasr'] = createRatio(h[daught1 + 'SmearedMass'],h['TotalSmearedMass'])
        h[daught2 + 'ProbMeasr'] = createRatio(h[daught2 + 'SmearedMass'],h['TotalSmearedMass'])

def makePlots2(HiddPart,daught1,daught2):
    
    key='DAUGHTERS'
    title='Muons are Blue, Kaons are Red and so are you'
        
    ut.bookCanvas(h,key + '_TV',title,nx=1300,ny=800,cx=3,cy=2)
    cv = h[key + '_TV'].cd(1)
    h[daught1 + 'StrawTime'].SetXTitle('Time [ns]')
    h[daught1 + 'StrawTime'].SetYTitle('No. of Particles')
    h[daught1 + 'StrawTime'].Draw()
    h[daught2 + 'StrawTime'].SetLineColor(2)
    h[daught2 + 'StrawTime'].Draw('same')

    cv = h[key + '_TV'].cd(2)
    h[daught1 + 'EcalTime'].SetXTitle('Time [ns]')
    h[daught1 + 'EcalTime'].SetYTitle('No. of Particles')
    h[daught1 + 'EcalTime'].Draw()
    h[daught2 + 'EcalTime'].SetLineColor(2)
    h[daught2 + 'EcalTime'].Draw('same')

    cv = h[key + '_TV'].cd(3)
    h[daught1 + 'DirDeltaTime'].SetXTitle('Time of Flight [ns]')
    h[daught1 + 'DirDeltaTime'].SetYTitle('No. of Particles')
    h[daught1 + 'DirDeltaTime'].Draw()
    h[daught2 + 'DirDeltaTime'].SetLineColor(2)
    h[daught2 + 'DirDeltaTime'].Draw('same')

    cv = h[key + '_TV'].cd(4)
    h[daught1 + 'FlightLen'].SetXTitle('Flight Length (cm)')
    h[daught1 + 'FlightLen'].SetYTitle('No. of Particles')
    h[daught1 + 'FlightLen'].Draw()
    h[daught2 + 'FlightLen'].SetLineColor(2)
    h[daught2 + 'FlightLen'].Draw('same')

    cv = h[key + '_TV'].cd(5)
    h[daught1 + 'Speed'].SetXTitle('beta')
    h[daught1 + 'Speed'].SetYTitle('No. of Particles')
    h[daught1 + 'Speed'].Draw()
    h[daught2 + 'Speed'].SetLineColor(2)
    h[daught2 + 'Speed'].Draw('same')

    h[key + '_TV'].Print('DaughterTVProp'+ currentDate + '.png')

    ut.bookCanvas(h,key + '_MOM', title , nx=1300, ny=800, cx=3, cy=2)
    cv = h[key + '_MOM'].cd(1)
    h[daught1 + 'StrawMom'].SetXTitle('Momentum [GeV/c]')
    h[daught1 + 'StrawMom'].SetYTitle('No. of Particles')
    h[daught1 + 'StrawMom'].Draw()
    h[daught2 + 'StrawMom'].SetLineColor(2)
    h[daught2 + 'StrawMom'].Draw('same')

    cv = h[key + '_MOM'].cd(2)
    h[daught1 + 'EcalMom'].SetXTitle('Momentum [GeV/c]')
    h[daught1 + 'EcalMom'].SetYTitle('No. of Particles')
    h[daught1 + 'EcalMom'].Draw()
    h[daught2 + 'EcalMom'].SetLineColor(2)
    h[daught2 + 'EcalMom'].Draw('same')

    cv = h[key + '_MOM'].cd(3)
    h[daught1 + 'RecoMom'].SetXTitle('Momentum [GeV/c]')
    h[daught1 + 'RecoMom'].SetYTitle('No. of Particles')
    h[daught1 + 'RecoMom'].Draw()
    h[daught2 + 'RecoMom'].SetLineColor(2)
    h[daught2 + 'RecoMom'].Draw('same')

    #cv = h[key + '_MOM'].cd(4)
    #h[daught1 + 'TrueMom'].SetXTitle('Momentum [GeV/c]')
    #h[daught1 + 'TrueMom'].SetYTitle('No. of Particles')
    #h[daught1 + 'TrueMom'].Draw()
    #h[daught2 + 'TrueMom'].SetLineColor(2)
    #h[daught2 + 'TrueMom'].Draw('same')

    cv = h[key + '_MOM'].cd(4)
    h[daught1 + 'DeltaMom'].SetXTitle('Momentum [GeV/c]')
    h[daught1 + 'DeltaMom'].SetYTitle('No. of Particles')
    h[daught1 + 'DeltaMom'].Draw()
    h[daught2 + 'DeltaMom'].SetLineColor(2)
    h[daught2 + 'DeltaMom'].Draw('same')

    cv = h[key + '_MOM'].cd(5)
    h[daught1 + 'RecoMass'].SetXTitle('Mass [GeV/c2]')
    h[daught1 + 'RecoMass'].SetYTitle('No. of Particles')
    h[daught1 + 'RecoMass'].Draw()
    h[daught2 + 'RecoMass'].SetLineColor(2)
    h[daught2 + 'RecoMass'].Draw('same')

    cv = h[key + '_MOM'].cd(6)
    h[daught1 + 'SmearedMass'].SetXTitle('Mass [GeV/c2]')
    h[daught1 + 'SmearedMass'].SetYTitle('No. of Particles')
    h[daught1 + 'SmearedMass'].Draw()
    h[daught1 + 'SmearedMass'].Fit("landau")
    h[daught1 + 'SmearedMass'].GetFunction("landau").SetLineColor(kBlack)
    h[daught2 + 'SmearedMass'].SetLineColor(2)
    h[daught2 + 'SmearedMass'].Draw('same')
    h[daught2 + 'SmearedMass'].Fit("landau")
    h[daught2 + 'SmearedMass'].GetFunction("landau").SetLineColor(kBlack)

    h['DAUGHTERS_MOM'].Print('DaughterPProp'+ currentDate + '.png')

    ut.bookCanvas(h,key + '_PROB',title,nx=1300,ny=600,cx=3,cy=1)
    cv = h[key + '_PROB'].cd(1)
    h[daught1 + 'ProbMeasr'].SetMarkerColor(38)
    polyFit1.SetLineColor(4)
    h[daught1 + 'ProbMeasr'].Fit('polyFit1')
    h[daught1 + 'ProbMeasr'].Draw('E2')
    h[daught1 + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
    h[daught1 + 'ProbMeasr'].SetYTitle('Prob(particle=(kaon or muon))')
    h[daught1 + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    cv = h[key + '_PROB'].cd(2)
    h[daught2 + 'ProbMeasr'].SetMarkerColor(46)
    polyFit2.SetLineColor(2)
    h[daught2 + 'ProbMeasr'].Fit('polyFit2')
    h[daught2 + 'ProbMeasr'].Draw('E2')
    h[daught2 + 'ProbMeasr'].SetXTitle('Mass [GeV/c2]')
    h[daught2 + 'ProbMeasr'].SetYTitle('Prob(particle=(kaon or muon))')
    h[daught2 + 'ProbMeasr'].GetYaxis().SetTitleOffset(1.5)

    cv = h[key + '_PROB'].cd(3)
    multigr = TMultiGraph()
    #gStyle.SetOptTitle(kFALSE)
    #gStyle.SetPalette(kSolar)
    #n = 300
    x1, y1 = array( 'd' ), array( 'd' )
    x2, y2 = array( 'd' ), array( 'd' )
    i=0
    n=0
    for i in range(30,240,6):
        x1.append(h[daught1 + 'ProbMeasr'].GetBinCenter(i))
        y1.append(h[daught1 + 'ProbMeasr'].GetBinContent(i))
        x2.append(h[daught2 + 'ProbMeasr'].GetBinCenter(i))
        y2.append(h[daught2 + 'ProbMeasr'].GetBinContent(i))
        n=n+1
    gr1 = TGraph( n, x1, y1 )
    gr1.SetTitle('Prob(ID = ' + daught1 + ')')
    gr2 = TGraph( n, x2, y2 )
    gr2.SetTitle('Prob(ID = ' + daught2 + ')')
    gr1.SetLineColor( 4 )
    gr1.SetLineWidth( 3 )
    gr1.SetMarkerColor( 4 )
    gr1.SetMarkerStyle( 20 )
    gr1.SetMarkerSize(0.7)   
    #gr1.GetXaxis().SetRangeUser(0,1.5)   
    gr2.SetLineColor( 2 )
    gr2.SetLineWidth( 3 )
    gr2.SetMarkerColor( 2 )
    gr2.SetMarkerStyle( 20 )
    gr2.SetMarkerSize(0.7)
    #gr2.GetXaxis().SetRangeUser(0,1.5)
    multigr.Add(gr1, "PC")
    multigr.Add(gr2, "PC")
    multigr.Draw("A pfc plc")#P PLC PFCPLC PFC
    multigr.GetXaxis().SetTitle( 'Mass [GeV/c2]' )
    multigr.GetYaxis().SetTitle( 'Prob(particle=(kaon or muon))' )
    multigr.GetYaxis().SetTitleOffset(1.5)
    #gr1.Draw("CA* PLC PFC")
    #gr2.Draw("PC  PLC PFC")
    gPad.BuildLegend()
    h[key + '_PROB'].Print('DaughterProb'+ currentDate + '.png')
    
def print_menu(): 
    print (30 * "-" + "MENU" + 30 * "-")
    print ("1. RPV SUSY Benchmark1 --> K+/- mu+/- final state")
    print ("2. RPV SUSY Benchmark1 --> K*+/- mu+/- final state")
    print ("3. RPV SUSY Benchmark1 --> K*0 nu_mu final state")
    print ("4. Exit")
    print (67 * "-")

########################################
     
while loop:        
    print_menu()
    choice = input("Enter your choice [1-4]: ")
     
    if choice==1:     
        print ("RPV SUSY Benchmark1 --> K+/- mu+/- final state selected.")
        HiddenPart = 'Neutralino'
        daughter1 = 'Mu'
        daughter2 = 'K+/-'
        create_Hists(HiddenPart, daughter1, daughter2)
        finState2t1t2(HiddenPart,daughter1,daughter2)
        makePlots2(HiddenPart,daughter1,daughter2)
        loop=False
    elif choice==2:
        print ("RPV SUSY Benchmark1 --> K*+/- mu+/- final state selected.")
        HiddenPart = 'Neutralino'
        daughter1 = 'Mu'
        daughter2 = 'K*+/-'
        finState2t1t2(HiddenPart,daughter1,daughter2)
        makePlots2(HiddenPart,daughter1,daughter2)
        loop=False
    elif choice==3:
        print ("RPV SUSY Benchmark1 --> K*0 nu_mu final state selected")
        HiddenPart = 'Neutralino'
        daughter1 = 'Mu'
        daughter2 = 'K*0'
        finState2t1t2(HiddenPart,daughter1,daughter2)
        makePlots2(HiddenPart,daughter1,daughter2)
        loop=False
    elif choice==4:
        print ("EXIT was selected.")
        loop=False 
    else:
        #raw_input("No such option. Try again.")
        print("No such option. Try again.")

hfile = inputFile.split(',')[0].replace('_rec','_RPVeditana')#create outputFile
if hfile[0:4] == "/eos" or not inputFile.find(',')<0:
# do not write to eos, write to local directory 
  tmp = hfile.split('/')
  hfile = tmp[len(tmp)-1]                                   #occurs only for cern users
ROOT.gROOT.cd()
ut.writeHists(h,hfile) 