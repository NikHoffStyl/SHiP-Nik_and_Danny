# #########################################
#         RPV_EDITA Test Code             #
###########################################
import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
from array import array
import shipRoot_conf
import shipDet_conf
import shipVeto
from inspect import currentframe
import numpy as np
import datetime
#import random
shipRoot_conf.configure()

debug = False
chi2CutOff  = 4.
PDG = ROOT.TDatabasePDG.Instance()
fiducialCut = True
measCutFK = 25
measCutPR = 22
docaCut = 2.
c = 2.99792458*(10**8)
currentDate = datetime.datetime.now().strftime("_%H%M_%d_%m_%y")


def inputOptsArgs():
    inputFile  = None
    geoFile    = None
    dy         = None
    nEvents    = 9999999
    try:
            opts, args = getopt.getopt(sys.argv[1:], "n:f:g:Y", ["nEvents=","geoFile="])#command line options
    except getopt.GetoptError:
            # print help information and exit:
            print ' enter file name'
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
 print 'found ecal geo for ',ecalGeoFile
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
 print "setup calo reconstruction of ecalReconstructed objects"
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
def create_Hists():
    ###############################
    ####  Daughter Histograms  ####
    ut.bookHist(h,'MuonStrawTime','Gaussian Straw t measurement',200,310.,340.)
    ut.bookHist(h,'KaonStrawTime','Gaussian Straw t measurement',200,310.,340.)
    ut.bookHist(h,'MuonEcalTime','Gaussian Ecal t measurement',200,340.,385.)
    ut.bookHist(h,'KaonEcalTime','Gaussian Ecal t measurement',200,340.,385.)
    ut.bookHist(h,'MuonDirDeltaTime','Straw-ECAL Time of Flight (directly)',200,37.5,39.1)   # muon time of flight
    ut.bookHist(h,'KaonDirDeltaTime','Straw-ECAL Time of Flight (directly)',200,37.5,39.1)   # Kaon time of flight
    ut.bookHist(h,'MuonFlightLen','Straw-ECAL Flight Lenght',200,11.2,12.)                   # muon flight Length
    ut.bookHist(h,'KaonFlightLen','Straw-ECAL Flight Length',200,11.2,12.)                   # Kaon flight Length
    ut.bookHist(h,'MuonSpeed','Beta value',200,0.99,1.01)                                    # muon speed
    ut.bookHist(h,'KaonSpeed','Beta value',200,0.99,1.)                                    # Kaon speed
    ut.bookHist(h,'MuonStrawMom','Straw Momentum',200,-0.05,100.)          # muon momentum
    ut.bookHist(h,'KaonStrawMom','Straw Momentum',200,-0.05,100.)          # Kaon momentum
    ut.bookHist(h,'MuonEcalMom','Ecal Momentum',200,-0.05,100.)           # muon momentum
    ut.bookHist(h,'KaonEcalMom','Ecal Momentum',200,-0.05,100.)           # Kaon momentum
    ut.bookHist(h,'KaonDeltaMom','Straw-Ecal Momentum',200,-10.,100.)       # Kaon momentum
    ut.bookHist(h,'MuonDeltaMom','Straw-Ecal Momentum',200,-10.,100.)       # Kaon momentum
    ut.bookHist(h,'MuonRecoMom','Reco Momentum',200,-0.05,100.)              # muon momentum
    ut.bookHist(h,'KaonRecoMom','Reco Momentum',200,-0.05,100.)              # Kaon momentum
    ut.bookHist(h,'MuonTrueMom','True Momentum',200,-0.05,100.)              # muon momentum
    ut.bookHist(h,'KaonTrueMom','True Momentum',200,-0.05,100.)              # Kaon momentum
    ut.bookHist(h,'MuonRecoMass','Reco Mass',200,0.,2.)                      # muon mass
    ut.bookHist(h,'KaonRecoMass','Reco Mass',200,0.,2.)                      # Kaon mass
    ut.bookHist(h,'MuonTrueMass','True Mass',200,0.,2.)                      # muon mass
    ut.bookHist(h,'KaonTrueMass','True Mass',200,0.,2.)                      # Kaon mass
    ut.bookHist(h,'MuonSmearedMass','Smeared Mass',200,0.,2.)                      # muon mass
    ut.bookHist(h,'KaonSmearedMass','Smeared Mass',200,0.,2.)                      # Kaon mass
    ################################
    ut.bookHist(h,'nalino_true','Monte Carlo Mass',500,0.,2.)                                   # true mass
    ut.bookHist(h,'nalino_reco','Reconstructed Mass',500,0.,2.)                                 # reconstructed mass
    ut.bookHist(h,'nalino_no_iter','Reconstructed Mass (without track iterations)',500,0.,2.)
    ut.bookHist(h,'nalino_mom','True (red) & Reco. (blue) Momentum',100,0.,300.)                # true momentum distribution
    ut.bookHist(h,'nalino_mom_reco','Reconstructed Momentum',100,0.,300.)                       # reconstructed momentum distribution
    ut.bookHist(h,'nalino_mom_diff','True/Reco Momentum Difference',100,-3.,3)                  # true/reco momentum difference

    ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared',100,0.,3.)                                 # chi squared track fitting
    ut.bookHist(h,'normdistr','Gaussian Distribution',500,-0.05,0.05)
    ut.bookHist(h,'smearedmass1','Time Smeared Neutralino Mass',500,0.,2.)
    ut.bookHist(h,'smearedmass2','Time Smeared Neutralino Mass',500,0.,2.)
    ut.bookHist(h,'smearedP1','Time Smeared Neutralino Momentum P1(red) P2(blue)',500,0.,300.)
    ut.bookHist(h,'smearedP2','Time Smeared Neutralino Momentum',500,0.,300.)

    
    print("Created Histograms")
create_Hists()

def dist2InnerWall(X,Y,Z):
  dist = 0
 # return distance to inner wall perpendicular to z-axis, if outside decayVolume return 0.
  node = sGeo.FindNode(X,Y,Z)
  if ShipGeo.tankDesign < 5:
     if not 'cave' in node.GetName():
         return dist  # TP 
  else:
     if not 'decayVol' in node.GetName():#does this
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
        #doesnt do this
        return 0    
    distance = sGeo.GetStep()
    if distance < minDistance  :    #does this
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
                 print 'SHiPAna: extrapolation did not work'
                 rc = False
                 break
             newPosDir[tr] = [reps[tr].getPos(states[tr]),reps[tr].getDir(states[tr])]
         if not rc: break
         xv,yv,zv,doca = myVertex(t1,t2,newPosDir)
         dz = abs(zBefore-zv)
         step+=1
         if step > 10:
             print 'abort iteration, too many steps, pos=',xv,yv,zv,' doca=',doca,'z before and dz',zBefore,dz
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
     nalinoLV = LV[t1]+LV[t2]
     #return xv,yv,zv,doca,nalinoMom
     return nalinoLV,LV[t1],LV[t2],doca

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

def makePlots():
   ut.bookCanvas(h,key='DAUGHTERS_TV',title='Muons are Blue, Kaons are Red and so are you',nx=1300,ny=1000,cx=3,cy=2)
   cv = h['DAUGHTERS_TV'].cd(1)
   h['MuonStrawTime'].SetXTitle('Time [ns]')
   h['MuonStrawTime'].SetYTitle('No. of Particles')
   h['MuonStrawTime'].Draw()
   h['KaonStrawTime'].SetLineColor(2)
   h['KaonStrawTime'].Draw('same')

   cv = h['DAUGHTERS_TV'].cd(2)
   h['MuonEcalTime'].SetXTitle('Time [ns]')
   h['MuonEcalTime'].SetYTitle('No. of Particles')
   h['MuonEcalTime'].Draw()
   h['KaonEcalTime'].SetLineColor(2)
   h['KaonEcalTime'].Draw('same')

   cv = h['DAUGHTERS_TV'].cd(3)
   h['MuonDirDeltaTime'].SetXTitle('Time of Flight [ns]')
   h['MuonDirDeltaTime'].SetYTitle('No. of Particles')
   h['MuonDirDeltaTime'].Draw()
   h['KaonDirDeltaTime'].SetLineColor(2)
   h['KaonDirDeltaTime'].Draw('same')

   cv = h['DAUGHTERS_TV'].cd(4)
   h['MuonFlightLen'].SetXTitle('Flight Length (cm)')
   h['MuonFlightLen'].SetYTitle('No. of Particles')
   h['MuonFlightLen'].Draw()
   h['KaonFlightLen'].SetLineColor(2)
   h['KaonFlightLen'].Draw('same')

   cv = h['DAUGHTERS_TV'].cd(5)
   h['MuonSpeed'].SetXTitle('beta')
   h['MuonSpeed'].SetYTitle('No. of Particles')
   h['MuonSpeed'].Draw()
   h['KaonSpeed'].SetLineColor(2)
   h['KaonSpeed'].Draw('same')

   h['DAUGHTERS_TV'].Print('DaughterTVProp'+ currentDate + '.png')

   ut.bookCanvas(h,key='DAUGHTERS_MOM',title='Muons are Blue, Kaons are Red and so are you',nx=1300,ny=1000,cx=3,cy=2)
   cv = h['DAUGHTERS_MOM'].cd(1)
   h['MuonStrawMom'].SetXTitle('Momentum [GeV/c]')
   h['MuonStrawMom'].SetYTitle('No. of Particles')
   h['MuonStrawMom'].Draw()
   h['KaonStrawMom'].SetLineColor(2)
   h['KaonStrawMom'].Draw('same')

   cv = h['DAUGHTERS_MOM'].cd(2)
   h['MuonEcalMom'].SetXTitle('Momentum [GeV/c]')
   h['MuonEcalMom'].SetYTitle('No. of Particles')
   h['MuonEcalMom'].Draw()
   h['KaonEcalMom'].SetLineColor(2)
   h['KaonEcalMom'].Draw('same')

   cv = h['DAUGHTERS_MOM'].cd(3)
   h['MuonRecoMom'].SetXTitle('Momentum [GeV/c]')
   h['MuonRecoMom'].SetYTitle('No. of Particles')
   h['MuonRecoMom'].Draw()
   h['KaonRecoMom'].SetLineColor(2)
   h['KaonRecoMom'].Draw('same')

   #cv = h['DAUGHTERS_MOM'].cd(4)
   #h['MuonTrueMom'].SetXTitle('Momentum [GeV/c]')
   #h['MuonTrueMom'].SetYTitle('No. of Particles')
   #h['MuonTrueMom'].Draw()
   #h['KaonTrueMom'].SetLineColor(2)
   #h['KaonTrueMom'].Draw('same')

   cv = h['DAUGHTERS_MOM'].cd(4)
   h['MuonDeltaMom'].SetXTitle('Momentum [GeV/c]')
   h['MuonDeltaMom'].SetYTitle('No. of Particles')
   h['MuonDeltaMom'].Draw()
   h['KaonDeltaMom'].SetLineColor(2)
   h['KaonDeltaMom'].Draw('same')

   cv = h['DAUGHTERS_MOM'].cd(5)
   h['MuonRecoMass'].SetXTitle('Mass [GeV/c2]')
   h['MuonRecoMass'].SetYTitle('No. of Particles')
   h['MuonRecoMass'].Draw()
   h['KaonRecoMass'].SetLineColor(2)
   h['KaonRecoMass'].Draw('same')

   cv = h['DAUGHTERS_MOM'].cd(6)
   h['MuonSmearedMass'].SetXTitle('Mass [GeV/c2]')
   h['MuonSmearedMass'].SetYTitle('No. of Particles')
   h['MuonSmearedMass'].Draw()
   h['KaonSmearedMass'].SetLineColor(2)
   h['KaonSmearedMass'].Draw('same')
   h['DAUGHTERS_MOM'].Print('DaughterPProp'+ currentDate + '.png')


def isInFiducial(X,Y,Z):
   if Z > ShipGeo.TrackStation1.z : return False
   if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
   # typical x,y Vx resolution for exclusive nalino decays 0.3cm,0.15cm (gaussian width)
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

def finState2MuK():
    if sTree.GetBranch("FitTracks"):
        k_decaycheck = 0
        for n in range(nEvents):                            # loop over events
            rc = sTree.GetEntry(n)                              # load tree entry
            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles                                   
                muPartkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]                  # gives MC particle data
                #print(reco_part)
                if abs(true_muon.GetPdgCode()) == 13:               # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()             # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[muonMotherkey]          # obtains mother particle data
                    muonMotherTrue_mass = true_mother.GetMass()        # get nalino/final states mother mass
                    muonMotherTrue_mom = true_mother.GetP()            # get nalino/final states mother mom
                    if true_mother.GetPdgCode() == 211:
                        #print('Kaon has decayed to a muon')
                        k_decaycheck+=1
                    if true_mother.GetPdgCode() == 9900015:             # checks mother is nalino

                        Decay_X = true_muon.GetStartX()
                        Decay_Y = true_muon.GetStartY()
                        Decay_Z = true_muon.GetStartZ()
                        
                        if not isInFiducial(Decay_X,Decay_Y,Decay_Z):
                            #print('nalino decayed outside fiducial volume')
                            continue
                        if not checkFiducialVolume(sTree,index,dy): 
                            #print('Track outside fiducial volume')
                            continue 
                        mu_status = reco_part.getFitStatus()             
                        if not mu_status.isFitConverged():
                            #print('Fit did not converge')
                            continue
                        mu_nmeas = mu_status.getNdf()                      
                        if not mu_nmeas > 25:
                            #print('Too few measurements')
                            continue

                        mu_rchi2 = mu_status.getChi2()                      # gets chi squared value
                        mu_chi2 = (mu_rchi2/mu_nmeas)                       # gets chi value
                        if not mu_chi2 < 4:
                            #print('Chi squared value too high')
                            continue

                        for index2,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles
                            kPartkey = sTree.fitTrack2MC[index2]                 # matches track to MC particle key
                            true_kaon = sTree.MCTrack[kPartkey]                  # gives MC particle data
                            if abs(true_kaon.GetPdgCode()) == 321:               # checks particle is kaon
                                kaonMotherkey = true_kaon.GetMotherId()             # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[kaonMotherkey]          # obtains mother particle data
                                if kaonMotherkey==muonMotherkey:                    # check if keys are the same
                                    kaonMotherTrue_mass = true_mother.GetMass()        # get nalino/final states mother mass
                                    kaonMotherTrue_mom = true_mother.GetP()            # get nalino/final states mother mom

                                    if not checkFiducialVolume(sTree,index,dy): 
                                        #print('Decay outside fiducial volume')
                                        continue 
                                    k_status = reco_part.getFitStatus()                
                                    if not k_status.isFitConverged():
                                        #print('Fit did not converge')
                                        continue
                                    k_nmeas = k_status.getNdf() 
                                    if not k_nmeas > 25:
                                        #print('Too few measurements')
                                        continue

                                    k_rchi2 = k_status.getChi2()                       # chi squared value
                                    k_chi2 = (k_rchi2/k_nmeas)                         # gets chi value
                                    if not k_chi2 < 4:
                                        #print('Chi squared value too high')
                                        continue
                                    
                                    Muon_LVec = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                    Kaon_LVec = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                    nalino_LVec = ROOT.TLorentzVector()               # declares variable as TLorentzVector class
                                    nalino_LVec,Muon_LVec,Kaon_LVec,doca = RedoVertexing(index,index2) # uses RedoVertexing to iterate track fitting
                                    if nalino_LVec == -1: continue
                                    if doca > 2.: 
                                        #print('distance of closest approach too large')
                                        continue

                                    nalino_mass = nalino_LVec.M()
                                    nalino_reco_mom = nalino_LVec.P()
                                    mom_diff = kaonMotherTrue_mom - nalino_reco_mom

                                    kM = Kaon_LVec.M()
                                    kP = Kaon_LVec.P()
                                    kE = Kaon_LVec.E()
                                    muM = Muon_LVec.M()
                                    muP = Muon_LVec.P()
                                    muE = Muon_LVec.E()
                                    
                                    kaonTrueMom = true_kaon.GetP()
                                    kaonTrueMass = true_kaon.GetMass()
                                    muonTrueMom = true_muon.GetP()
                                    muonTrueMass = true_muon.GetMass()
                                    h['MuonTrueMom'].Fill(muonTrueMom)
                                    h['MuonTrueMass'].Fill(muonTrueMass)
                                    h['KaonTrueMom'].Fill(kaonTrueMom)
                                    h['KaonTrueMass'].Fill(kaonTrueMass)
                                    h['nalino_true'].Fill(kaonMotherTrue_mass)            
                                    h['nalino_mom'].Fill(kaonMotherTrue_mom)
                                    h['nalino_reco'].Fill(nalino_mass)                        
                                    h['nalino_mom_reco'].Fill(nalino_reco_mom)                
                                    h['Chi2'].Fill(mu_chi2)       
                                    h['Chi2'].Fill(k_chi2)                             
                                    h['nalino_mom_diff'].Fill(mom_diff)
                                    h['KaonRecoMom'].Fill(kP)
                                    h['MuonRecoMom'].Fill(muP)
                                    h['KaonRecoMass'].Fill(kM)
                                    h['MuonRecoMass'].Fill(muM)

                                    mu_strwT,mu_ecalT,mu_Dt,mu_Len,mu_v,mu_strawP,mu_ecalP = time_res(muPartkey) 
                                    if mu_strwT!= None and mu_ecalT!= None and mu_Dt!= None and mu_Len!= None and mu_v != None and mu_strawP!=None and mu_ecalP!=None:
                                        h['MuonStrawTime'].Fill(mu_strwT)
                                        h['MuonEcalTime'].Fill(mu_ecalT)
                                        h['MuonDirDeltaTime'].Fill(mu_Dt)
                                        h['MuonFlightLen'].Fill(mu_Len)
                                        mu_beta = mu_v/c
                                        h['MuonSpeed'].Fill(mu_beta)
                                        h['MuonStrawMom'].Fill(mu_strawP)
                                        h['MuonEcalMom'].Fill(mu_ecalP)
                                        mu_Dp = mu_strawP-mu_ecalP
                                        h['MuonDeltaMom'].Fill(mu_Dp)

                                        k_strwT,k_ecalT,k_Dt,k_Len,k_v,k_strawP,k_ecalP = time_res(kPartkey)
                                        if k_strwT!= None and k_ecalT!= None and k_Dt!= None and k_Len!= None and k_v != None and k_strawP!=None and k_ecalP!=None:     
                                            h['KaonStrawTime'].Fill(k_strwT)
                                            h['KaonEcalTime'].Fill(k_ecalT)
                                            h['KaonDirDeltaTime'].Fill(k_Dt)
                                            h['KaonFlightLen'].Fill(k_Len)
                                            k_beta = k_v/c
                                            h['KaonSpeed'].Fill(k_beta)
                                            h['KaonStrawMom'].Fill(k_strawP)
                                            h['KaonEcalMom'].Fill(k_ecalP)
                                            k_Dp = k_strawP-k_ecalP
                                            h['KaonDeltaMom'].Fill(k_Dp)
                                              
                                            
                                            if k_beta < 1:
                                                #TRYING SOMETHING ELSE
                                                k_smearedM = kP*(ROOT.TMath.Sqrt(1-(k_beta**2)))/k_beta
                                                h['KaonSmearedMass'].Fill(k_smearedM)
                                                mu_smearedM = muP*(ROOT.TMath.Sqrt(1-(mu_beta**2)))/mu_beta
                                                h['MuonSmearedMass'].Fill(mu_smearedM)

        #print('\n'+str(k_decaycheck) + ' K+ --> mu decays before detection\n')

#########################################


finState2MuK()          
makePlots()
hfile = inputFile.split(',')[0].replace('_rec','_RPVeditana')#create outputFile
if hfile[0:4] == "/eos" or not inputFile.find(',')<0:
# do not write to eos, write to local directory 
  tmp = hfile.split('/')
  hfile = tmp[len(tmp)-1]                                   #occurs only for cern users
ROOT.gROOT.cd()
ut.writeHists(h,hfile) 