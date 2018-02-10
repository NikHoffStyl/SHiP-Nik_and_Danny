import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
import shipRoot_conf
import numpy as np
shipRoot_conf.configure()

debug = False
PDG = ROOT.TDatabasePDG.Instance()
inputFile  = None
geoFile    = None
dy         = None
nEvents    = 9999999
fiducialCut = True
chi2Cut  = 4
measCut = 25
docaCut = 2

#--------------------------------------------------INPUT----------------------------------------------------

try:
        opts, args = getopt.getopt(sys.argv[1:], "n:f:g:A:Y:i", ["nEvents=","geoFile="])
except getopt.GetoptError:
        # print help information and exit:
        print ' enter file name'
        sys.exit()
for o, a in opts:
        if o in ("-f",):
            inputFile = a
        if o in ("-g", "--geoFile",):
            geoFile = a
        if o in ("-Y",):
            dy = float(a)
        if o in ("-n", "--nEvents=",):
            nEvents = int(a)

if not inputFile.find(',')<0 :  
  sTree = ROOT.TChain("cbmsim")
  for x in inputFile.split(','):
   if x[0:4] == "/eos":
    sTree.AddFile("root://eoslhcb.cern.ch/"+x)
   else: sTree.AddFile(x)
elif inputFile[0:4] == "/eos":
  eospath = "root://eoslhcb.cern.ch/"+inputFile
  f = ROOT.TFile.Open(eospath)
  sTree = f.cbmsim
else:
  f = ROOT.TFile(inputFile)
  sTree = f.cbmsim

#-------------------------------------------------GEOMETRY--------------------------------------------------

# try to figure out which ecal geo to load
if not geoFile:
 geoFile = inputFile.replace('ship.','geofile_full.').replace('_rec.','.')
if geoFile[0:4] == "/eos":
  eospath = "root://eoslhcb.cern.ch/"+geoFile
  fgeo = ROOT.TFile.Open(eospath)
else:  
  fgeo = ROOT.TFile(geoFile)
sGeo = fgeo.FAIRGeom

if not fgeo.FindKey('ShipGeo'):
 # old geofile, missing Shipgeo dictionary
 if sGeo.GetVolume('EcalModule3') :  ecalGeoFile = "ecal_ellipse6x12m2.geo"
 else: ecalGeoFile = "ecal_ellipse5x10m2.geo" 
 print 'found ecal geo for ',ecalGeoFile
 # re-create geometry and mag. field
 if not dy:
  # try to extract from input file name
  tmp = inputFile.split('.')
  try:
    dy = float( tmp[1]+'.'+tmp[2] )
  except:
    dy = 10.
 ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py", Yheight = dy, EcalGeoFile = ecalGeoFile )
else: 
 # new geofile, load Shipgeo dictionary written by run_simScript.py
  upkl    = Unpickler(fgeo)
  ShipGeo = upkl.load('ShipGeo')
  ecalGeoFile = ShipGeo.ecal.File
  dy = ShipGeo.Yheight/u.m

import shipDet_conf
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
import shipVeto
veto = shipVeto.Task(sTree)
vetoDets={}

# fiducial cuts
vetoStation = ROOT.gGeoManager.GetTopVolume().GetNode('Veto_5')
vetoStation_zDown = vetoStation.GetMatrix().GetTranslation()[2]+vetoStation.GetVolume().GetShape().GetDZ()
T1Station = ROOT.gGeoManager.GetTopVolume().GetNode('Tr1_1')
T1Station_zUp = T1Station.GetMatrix().GetTranslation()[2]-T1Station.GetVolume().GetShape().GetDZ()

# initialize ecal structure
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

import TrackExtrapolateTool
from array import array
c = 2.99792458*(10**8)   # speed of light
e = 2.718281828459   # Euler's number
h = {}   # creates empty dictionary

#----------------------------------------------------FUNCTIONS------------------------------------------------------------

def dist2InnerWall(X,Y,Z):
  dist = 0
 # return distance to inner wall perpendicular to z-axis, if outside decay volume return 0.
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
   # typical x,y Vx resolution for exclusive RPV decays 0.3cm,0.15cm (gaussian width)
   if dist2InnerWall(X,Y,Z)<5*u.cm: return False
   return True

def checkFiducialVolume(sTree,tkey,dy):
# extrapolate track to middle of magnet and check if in decay volume
   inside = True
   if not fiducialCut: return True
   fT = sTree.FitTracks[tkey]
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
   Z = c.z()+v.z()*t        # X,Y,Z are the coordinates of...?
   return X,Y,Z,abs(dist) 

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
      xx  = sTree.FitTracks[tr].getFittedState()
      PosDir[tr] = [xx.getPos(),xx.getDir()] # Dir = direction?
     xv,yv,zv,doca = myVertex(t1,t2,PosDir) # myVertex returns 3 X,Y,Z and abs(dist)
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
       reps[tr]   = ROOT.genfit.RKTrackRep(xx.getPDG())     # what is rep?
       states[tr] = ROOT.genfit.StateOnPlane(reps[tr])     # what is this?
       reps[tr].setPosMom(states[tr],xx.getPos(),xx.getMom())   # and this?
       try:
        reps[tr].extrapolateToPoint(states[tr], newPos, False)
       except: 
        #print('SHiPAna: extrapolation did not work)'
        rc = False  
        break # breaks if extrapolation doesn't work

       newPosDir[tr] = [reps[tr].getPos(states[tr]),reps[tr].getDir(states[tr])]
      if not rc: break
      xv,yv,zv,doca = myVertex(t1,t2,newPosDir)
      dz = abs(zBefore-zv) # repeats until dz < 0.1 unless...
      step+=1
      if step > 10:  
         # print 'abort iteration, too many steps, pos=',xv,yv,zv,' doca=',doca,'z before and dz',zBefore,dz
         rc = False
         break # ... number of iterations exceeds 10
     if not rc: return -1,-1,-1,doca # extrapolation failed, makes no sense to continue
     LV={}
     for tr in [t1,t2]: # from here on we have reproduced (see inv_mass() function in Test_1.py)     
      mom = reps[tr].getMom(states[tr])
      pid = abs(states[tr].getPDG()) 
      if pid == 2212: pid = 211 # why?
      mass = PDG.GetParticle(pid).Mass()
      E = ROOT.TMath.Sqrt(mass*mass + mom.Mag2())
      LV[tr] = ROOT.TLorentzVector()
      LV[tr].SetPxPyPzE(mom.x(),mom.y(),mom.z(),E)
     RPVMom = LV[t1]+LV[t2]

     return RPVMom,LV[t1],LV[t2],doca

def create_Hists():
    ut.bookHist(h,'RPV_true','Monte Carlo Mass',500,0.8,1.2)   # true mass
    ut.bookHist(h,'RPV_reco','Reconstructed Mass',500,0.8,1.2)   # reconstructed mass
    ut.bookHist(h,'RPV_mom','True (red) & Reco. (blue) Momentum',100,0.,300.)   # true momentum distribution
    ut.bookHist(h,'RPV_mom_reco','Reconstructed Momentum',100,0.,300)   # reconstructed momentum distribution
    ut.bookHist(h,'RPV_mom_diff','True/Reco Momentum Difference',100,-3.,3)   # true/reco momentum difference

    ut.bookHist(h,'Muon_mom_true','Muon - True (red) & Reco. (blue) Momentum',100,0.,140.)   # RPV muon daughter reco momentum
    ut.bookHist(h,'Kaon_mom_true','Kaon - True (red) & Reco. (blue) Momentum',100,0.,140.)   # RPV pion daughter reco momentum
    ut.bookHist(h,'Muon_mom','Muon - True Momentum',100,0.,140.)   # RPV muon daughter true momentum
    ut.bookHist(h,'Kaon_mom','Kaon - True Momentum',100,0.,140.)   # RPV pion daughter true momentum
    ut.bookHist(h,'ecalstraw_mom','Straw-Ecal Momentum Difference',500,0,0.4)   # includes both kaon and muon
    ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared',100,0.,3.)   # chi squared track fitting

    ut.bookHist(h,'MuonDir','Smeared Muon Straw-ECAL Time (directly)',150,37.5,40.)   # daughter muon time of flight (Gaussian blurred)
    ut.bookHist(h,'KaonDir','Smeared Kaon Straw-ECAL Time (directly)',150,37.5,40.)   # daughter kaon time of flight (Gaussian blurred)
    ut.bookHist(h,'tmass_muon','Time Deduced Muon Mass',150,0.,3.)   # time, momentum --> velocity --> gamma (L) --> mass from p=mvL
    ut.bookHist(h,'tmass_kaon','Time Deduced Kaon(red)-Muon(blue) Mass',150,0.,3.)
    ut.bookHist(h,'tsmearmass_muon','Smeared Time Deduced Kaon(red)-Muon(blue) Mass',150,0.,3.)   # same as above but using smeared time
    ut.bookHist(h,'tsmearmass_kaon','Smeared Time Deduced Kaon(red)-Muon(blue) Mass',150,0.,3.)
    ut.bookHist(h,'Daughter_masses','True Masses of Daughter Particles',100,0.,1.)   # kaon and muon true mass

    ut.bookHist(h,'MuonDir_nosmear','True Muon Straw-ECAL Time (directly)',150,37.5,40.)   # daughter muon time of flight
    ut.bookHist(h,'KaonDir_nosmear','True Kaon Straw-ECAL Time (directly)',150,37.5,40.)   # daughter kaon time of flight
    ut.bookHist(h,'num_muon','No. of muon hits in straw tubes',25,25,50)
    ut.bookHist(h,'num_kaon','No. of kaon in straw tubes',25,25,50)
    ut.bookHist(h,'track_muon','Muon z-momentum through straw tubes (for particular event)',500,20,21)
    ut.bookHist(h,'track_kaon','Kaon z-momentum through straw tubes (for particular event)',500,44,45)
    ut.bookHist(h,'straight_path','Straight distance from first straw tube hit to ECAL',200,11,12)
    ut.bookHist(h,'better_path','Curved path through straw tubes, then straight to ECAL',200,11,12)
    ut.bookHist(h,'path_diff','Difference between straight path and better approximation',100,0,0.001)

def makePlots():
    ut.bookCanvas(h,key='Test_Mass',title='Results 1',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['Test_Mass'].cd(1)
    h['RPV_true'].SetXTitle('Invariant mass [GeV/c2]')
    h['RPV_true'].SetYTitle('No. of Particles')
    h['RPV_true'].SetLineColor(1)
    h['RPV_true'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Test_Mass'].cd(2)
    h['RPV_reco'].SetXTitle('Invariant mass [GeV/c2]')
    h['RPV_reco'].SetYTitle('No. of Particles')
    h['RPV_reco'].SetLineColor(1)
    h['RPV_reco'].Draw()
    print('\nNeutralino mass Gaussian fit:')
    fitSingleGauss('RPV_reco',0.9,1.1)
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Test_Mass'].cd(3)
    h['RPV_mom'].SetXTitle('Momentum [GeV/c]')
    h['RPV_mom'].SetYTitle('No. of Particles')
    h['RPV_mom'].SetLineColor(2)
    h['RPV_mom'].Draw()
    h['RPV_mom_reco'].Draw("same")
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Test_Mass'].cd(4)
    h['RPV_mom_diff'].SetXTitle('Momentum Difference [GeV/c]')
    h['RPV_mom_diff'].SetYTitle('Frequency')
    h['RPV_mom_diff'].Draw()
    h['RPV_mom_diff'].SetLineColor(1)
    h['Test_Mass'].Print('RPV_N.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='KaMu_Graphs',title='Results 2',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['KaMu_Graphs'].cd(1)
    h['Kaon_mom_true'].SetXTitle('Momentum [GeV/c]')
    h['Kaon_mom_true'].SetYTitle('No. of particles')
    h['Kaon_mom_true'].SetLineColor(2)
    h['Kaon_mom_true'].Draw()
    h['Kaon_mom'].Draw("same")
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['KaMu_Graphs'].cd(2)
    h['Muon_mom_true'].SetXTitle('Momentum [GeV/c]')
    h['Muon_mom_true'].SetYTitle('No. of particles')
    h['Muon_mom_true'].SetLineColor(2)
    h['Muon_mom_true'].Draw()
    h['Muon_mom'].Draw("same")
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['KaMu_Graphs'].cd(3)
    h['ecalstraw_mom'].SetXTitle('Momentum Difference [GeV/c2]')
    h['ecalstraw_mom'].SetYTitle('Frequency')
    h['ecalstraw_mom'].SetLineColor(1)
    h['ecalstraw_mom'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['KaMu_Graphs'].cd(4)
    h['Chi2'].SetXTitle('Chi Squared Value')
    h['Chi2'].SetYTitle('Frequency')
    h['Chi2'].SetLineColor(1)
    h['Chi2'].Draw()
    h['KaMu_Graphs'].Print('KaMu.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='Test_Time',title='Results 3',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['Test_Time'].cd(1)
    h['MuonDir'].SetXTitle('Time [ns]')
    h['MuonDir'].SetYTitle('Frequency')
    h['MuonDir'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Test_Time'].cd(2)
    h['KaonDir'].SetXTitle('Time [ns]')
    h['KaonDir'].SetYTitle('Frequency')
    h['KaonDir'].SetLineColor(2)
    h['KaonDir'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Test_Time'].cd(3)
    h['tmass_kaon'].SetXTitle('Mass [GeV/c2]')
    h['tmass_kaon'].SetYTitle('Frequency')
    h['tmass_kaon'].SetLineColor(2)
    h['tmass_kaon'].Draw()
    h['tmass_muon'].Draw('same')
    h['Daughter_masses'].Draw('same')
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Test_Time'].cd(4)
    h['tsmearmass_kaon'].SetXTitle('Mass [GeV/c2]')
    h['tsmearmass_kaon'].SetYTitle('Frequency')
    h['tsmearmass_kaon'].SetLineColor(2)
    h['tsmearmass_kaon'].Draw()
    print('\nLandau fits for mass (time of flight):')
    h['tsmearmass_kaon'].Fit('landau')
    h['tsmearmass_kaon'].GetFunction('landau').SetLineColor(1)
    #par0 = h['tsmearmass_kaon'].GetFunction("landau").GetParameter(0)
    #par1 = h['tsmearmass_kaon'].GetFunction("landau").GetParameter(1)
    #par2 = h['tsmearmass_kaon'].GetFunction("landau").GetParameter(2)
    h['tsmearmass_muon'].Draw('same')
    h['tsmearmass_muon'].Fit('landau') # alternatively --> .Fit("pol5") for 5th order polynomial
    h['tsmearmass_muon'].GetFunction('landau').SetLineColor(1)
    h['Daughter_masses'].SetLineColor(1)
    h['Daughter_masses'].SetLineStyle(2)
    h['Daughter_masses'].Draw('same')
    h['Test_Time'].Print('Time_Mass.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='prob',title='Results 4',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['prob'].cd(1)
    h['graph_mu'].SetTitle('Probability of particle being a muon')
    h['graph_mu'].GetXaxis().SetTitle('Mass / [GeV/c2]')
    h['graph_mu'].GetYaxis().SetTitle('Probability')
    h['graph_mu'].SetMarkerStyle(33)
    h['graph_mu'].SetMarkerColor(4)
    h['graph_mu'].SetMarkerSize(1)
    h['graph_mu'].Draw('AP')
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['prob'].cd(2)
    h['graph_ka'].SetTitle('Probability of particle being a kaon')
    h['graph_ka'].GetXaxis().SetTitle('Mass / [GeV/c2]')
    h['graph_ka'].GetYaxis().SetTitle('Probability')
    h['graph_ka'].SetMarkerStyle(33)
    h['graph_ka'].SetMarkerColor(2)
    h['graph_ka'].SetMarkerSize(1)
    h['graph_ka'].Draw('AP')
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['prob'].cd(3)
    h['combo'].SetTitle('Probability of identifying muon (blue) and kaon (red)')
    h['combo'].Draw('AP')
    h['prob'].Print('Probs.png')

def track_checks(index,true_part,reco_part,veto):
    check = 0
    Decay_X = true_part.GetStartX()
    Decay_Y = true_part.GetStartY()
    Decay_Z = true_part.GetStartZ()
    if not isInFiducial(Decay_X,Decay_Y,Decay_Z):
        #print('RPV decayed outside fiducial volume')
        veto[0] += 1
        check = -1
    if not checkFiducialVolume(sTree,index,dy): 
        #print('Track outside fiducial volume')
        veto[1] += 1
        check = -1
    fit_status = reco_part.getFitStatus()             
    if not fit_status.isFitConverged():
        #print('Fit did not converge')
        veto[2] += 1
        check = -1
    fit_nmeas = fit_status.getNdf()                      
    if not fit_nmeas > measCut:
        #print('Too few measurements')
        veto[3] += 1
        check = -1
    fit_rchi2 = fit_status.getChi2()                      
    fit_chi2 = (fit_rchi2/fit_nmeas)
    if not fit_chi2 < chi2Cut:
        #print('Chi squared value too high')
        veto[4] += 1
        check = -1

    return check,fit_chi2,veto

def time_res(partkey,pdg,n,m):
    tnosmear = -1   # declares variables
    vnosmear = -1
    tsmear = -1
    vsmear = -1
    diff = -1
    if sTree.GetBranch("strawtubesPoint"):
        x_array = []   # declares lists
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
                x_array.append(hits.GetX())   # adds data to the lists
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

        min_z_index = z_array.index(min(z_array))   # gives index of the smallest element in the list
        firststraw_z = 0.01*z_array[min_z_index]   # positions, time and momenta the first straw tube hit
        firststraw_x = 0.01*x_array[min_z_index]
        firststraw_y = 0.01*y_array[min_z_index]
        straw_time = t_array[min_z_index]
        strawPx = px_array[min_z_index]
        strawPy = py_array[min_z_index]
        strawPz = pz_array[min_z_index]
        strawP = ROOT.TMath.Sqrt((strawPx**2) + (strawPy**2) + (strawPz**2)) # straw tube momentum
            
        num_hits = len(z_array)   # number of elements in the list
        if pdg==13:   # muon
            h['num_muon'].Fill(num_hits)
            for hit in pz_array:
                if n == m:
                    h['track_muon'].Fill(hit)   # muon z-momentum through straw tubes for particular event
        if pdg==321:   # kaon
            h['num_kaon'].Fill(num_hits)
            h['num_kaon'].SetLineColor(2)
            for hit in pz_array:
                if n == m:
                    h['track_kaon'].Fill(hit)   # kaon z-momentum through straw tubes for particular event
                    h['track_kaon'].SetLineColor(2)
  
        if sTree.GetBranch("EcalPoint"):
            ecal_time = 0
            if not straw_time <= 0:
                for k,hits in enumerate(sTree.EcalPoint):
                    ecal_TrackID = hits.GetTrackID()
                    if ecal_TrackID == partkey:
                        ecal_x = 0.01*hits.GetX()   # positions, time and momenta of ECAL hit
                        ecal_y = 0.01*hits.GetY()   # stored in units of cm 
                        ecal_z = 0.01*hits.GetZ()
                        ecal_time = hits.GetTime()
                        ecalPx = hits.GetPx()
                        ecalPy = hits.GetPy()
                        ecalPz = hits.GetPz()
                        ecalP = ROOT.TMath.Sqrt((ecalPx**2) + (ecalPy**2) + (ecalPz**2))   # ECAL momentum

        if not ecal_time <= 0:
            diff = strawP - ecalP   # between 1st straw tube hit and ECAL
            r = ROOT.TMath.Sqrt(((ecal_x - firststraw_x)**2) + ((ecal_y - firststraw_y)**2) + ((ecal_z - firststraw_z)**2))
            h['straight_path'].Fill(r)
            h['straight_path'].SetLineColor(1)
            max_z_index = z_array.index(max(z_array))   # gives index of the smallest element in the list
            laststraw_x = 0.01*x_array[max_z_index]
            laststraw_y = 0.01*y_array[max_z_index]
            laststraw_z = 0.01*z_array[max_z_index]
            R2 = ROOT.TMath.Sqrt(((ecal_x - laststraw_x)**2) + ((ecal_y - laststraw_y)**2) + ((ecal_z - laststraw_z)**2))
            R = R1+R2   # better approximation of distance travelled through the straw tubes
            h['better_path'].Fill(R)
            h['better_path'].SetLineColor(1)
            rdiff = abs(R-r)
            h['path_diff'].Fill(rdiff)
            h['path_diff'].SetLineColor(1)
            sigma = 0.01   # standard deviation for Gaussian
            straw_smear = np.random.normal(loc=straw_time,scale=sigma,size=None)
            ecal_smear = np.random.normal(loc=ecal_time,scale=sigma,size=None)
            tsmear = abs(straw_smear - ecal_smear)   # smeared time of flight
            vsmear = (R/tsmear)*(10**9)   # smeared velocity of flight

            tnosmear = abs(straw_time - ecal_time)   # stored in units of nanoseconds
            vnosmear = (R/tnosmear)*(10**9)   # velocity of flight
            
    return tnosmear,vnosmear,tsmear,vsmear,diff,strawP

#----------------------------------------------------EVENT-LOOPS--------------------------------------------------------

nEvents = min(sTree.GetEntries(),nEvents)   # number of generated events

def finStateMuKa():
    if sTree.GetBranch("FitTracks"):
        print('\nRunning analysis for final state K+- Mu-+ :\n')
        create_Hists()   # creates histograms for this function
        successful_events = []   # creates list of event numbers of desired decays
        veto = 6*[0]   # creates list of veto counts for each possible veto cause
        rejected = 0   # variable for total number of events which are vetoed
        ka_decaycheck = 0   # variable for counting when kaons from neutralinos decay to muons before detection
        for n in range(nEvents):   # loops over events
            rc = sTree.GetEntry(n)   # loads tree entry

            #-----------------------------------------------TRACK-LOOPS------------------------------------------------

            for index,reco_part in enumerate(sTree.FitTracks):   # loops over index and data of track particles                                   
                muPartkey = sTree.fitTrack2MC[index]   # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]   # gives MC particle data
                if abs(true_muon.GetPdgCode()) == 13:   # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()   # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[muonMotherkey]   # obtains mother particle data

                    if true_mother.GetPdgCode() == 321:
                        motherN = sTree.MCTrack[true_mother.GetMotherId()]
                        if motherN.GetPdgCode() == 9900015:
                            ka_decaycheck += 1   # kaon has decayed to a muon before detection

                    if true_mother.GetPdgCode() == 9900015:   # checks mother is RPV
                        check,mu_chi2,veto = track_checks(index,true_muon,reco_part,veto)   # performs various checks (i.e. vertex position, fiducial volume,...)
                        if not check == 0: 
                            rejected += 1
                            continue

                        for index2,reco_part2 in enumerate(sTree.FitTracks):   # loops over index and data of track particles
                            Partkey = sTree.fitTrack2MC[index2]   # matches track to MC particle key
                            true_part = sTree.MCTrack[Partkey]   # gives MC particle data
                            if abs(true_part.GetPdgCode()) == 321:   # checks particle is kaon
                                true_kaon = true_part   # renames the particle object
                                kaPartkey = Partkey   # renames the particle code
                                kaonMotherkey = true_kaon.GetMotherId()   # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[kaonMotherkey]   # obtains mother particle data

                                if kaonMotherkey == muonMotherkey:   # check if mother keys are the same
                                    check2,ka_chi2,veto = track_checks(index2,true_kaon,reco_part2,veto)
                                    if not check2 == 0:   # performs various checks (i.e. vertex position, fiducial volume,...)
                                        rejected += 1
                                        continue

                                    #---------------------------------------------PARTICLE-DATA-----------------------------------------------------

                                    kaonMotherTrue_mass = true_mother.GetMass()   # get RPV/final states mother mass
                                    kaonMotherTrue_mom = true_mother.GetP()   # get RPV/final states mother mom

                                    RPV_Vector = ROOT.TLorentzVector()   # declares variables as TLorentzVector class
                                    Muon_Vector = ROOT.TLorentzVector()
                                    Kaon_Vector = ROOT.TLorentzVector()
                                    RPV_Vector,Muon_Vector,Kaon_Vector,doca = RedoVertexing(index,index2)   # uses RedoVertexing to iterate track fitting
                                    
                                    if RPV_Vector == -1: continue
                                    if doca > docaCut: 
                                        #print('distance of closest approach too large')
                                        veto[5] +=1
                                        rejected += 1
                                        continue

                                    RPV_mass = RPV_Vector.M()   # reconstructed RPV mass
                                    RPV_reco_mom = RPV_Vector.P()   # reconstructed RPV momentum
                                    mom_diff = kaonMotherTrue_mom - RPV_reco_mom   # RPV true/reco momentum difference
                                    true_kaP = true_kaon.GetP()   # true kaon momentum
                                    reco_kaP = Kaon_Vector.P()   # reconstructed kaon momentum
                                    true_muP = true_muon.GetP()   # true muon momentum
                                    reco_muP = Muon_Vector.P()   # reconstructed muon momentum
                                    kaM = true_kaon.GetMass()   # kaon mass
                                    muM = true_muon.GetMass()   # muon mass
                                                      
                                    h['RPV_true'].Fill(kaonMotherTrue_mass)   # fills histograms
                                    h['RPV_mom'].Fill(kaonMotherTrue_mom)
                                    h['RPV_reco'].Fill(RPV_mass)                        
                                    h['RPV_mom_reco'].Fill(RPV_reco_mom)   
                                    h['Daughter_masses'].Fill(muM)
                                    h['Daughter_masses'].Fill(kaM)
                                    h['Chi2'].Fill(mu_chi2)       
                                    h['Chi2'].Fill(ka_chi2)                             
                                    h['RPV_mom_diff'].Fill(mom_diff)
                                    h['Kaon_mom'].Fill(reco_kaP)
                                    h['Muon_mom'].Fill(reco_muP)
                                    h['Kaon_mom_true'].Fill(true_kaP)
                                    h['Muon_mom_true'].Fill(true_muP)
                                    
                                    successful_events.append(n)   # adds entries to the list
                                    m = successful_events[0]   # arbitrarily picks the first one as an example
                                    
                                    #------------------------------------TIME-RESOLUTION------------------------------------------

                                    mu_t,mu_v,mu_tsmear,mu_vsmear,mu_diff,straw_muP = time_res(muPartkey,13,n,m)        
                                    if mu_t != -1:   # and mu_t < 38.05:
                                        h['MuonDir'].Fill(mu_tsmear)   # fills histogram with smeared time
                                        h['MuonDir_nosmear'].Fill(mu_t)   # fills histogram with true time

                                        beta = mu_v/c   # equations for mass calculated from true time
                                        gamma = 1/(ROOT.TMath.Sqrt(1-(beta**2)))
                                        nosmearM = straw_muP/(beta*gamma)   # previously used reco_muP
                                        beta_smear = mu_vsmear/c   # equations for mass calculated from smeared time
                                        gamma_smear = 1/(ROOT.TMath.Sqrt(1-(beta_smear**2)))
                                        smearM = straw_muP/(beta_smear*gamma_smear)
                                        
                                        h['tmass_muon'].Fill(nosmearM)   # fills histograms with mass data
                                        h['tsmearmass_muon'].Fill(smearM)

                                        ka_t,ka_v,ka_tsmear,ka_vsmear,ka_diff,straw_kaP = time_res(kaPartkey,321,n,m)      
                                        if ka_t != -1:   # and ka_t < 38.06:
                                            h['KaonDir'].Fill(ka_tsmear)   # fills histogram with smeared time
                                            h['KaonDir_nosmear'].Fill(ka_t)   # fills histogram with true time
                                            h['KaonDir_nosmear'].SetLineColor(2)

                                            beta = ka_v/c   # equations for mass calculated from true time
                                            gamma = 1/(ROOT.TMath.Sqrt(1-(beta**2)))
                                            nosmearM = straw_kaP/(beta*gamma)   # previously used reco_kaP
                                            beta_smear = ka_vsmear/c   # equations for mass calculated from smeared time
                                            gamma_smear = 1/(ROOT.TMath.Sqrt(1-(beta_smear**2)))
                                            smearM = straw_kaP/(beta_smear*gamma_smear)

                                            h['tmass_kaon'].Fill(nosmearM)   # fills histograms with mass data
                                            h['tsmearmass_kaon'].Fill(smearM)
                                            h['ecalstraw_mom'].Fill(mu_diff)   # fills histogram for momentum difference
                                            h['ecalstraw_mom'].Fill(ka_diff)

        accepted = len(successful_events)
        total = accepted + rejected
        print('\t' + str(total) + ' events reconstructed for this decay mode')
        print('\t' + str(accepted) + ' events not rejected')
        print('\t' + str(total-accepted) + ' events rejected:')
        print('\t\t' + str(total-veto[0]/2) + ' vertices inside fiducial volume (' + str(veto[0]/2) + ' vetos)')
        print('\t\t' + str(total-veto[1]) + ' tracks within fiducial volume (' + str(veto[1]) + ' vetos)')
        print('\t\t' + str(total-veto[2]) + ' fits successfully converged (' + str(veto[2]) + ' vetos)')
        print('\t\t' + str(total-veto[3]) + ' tracks with no. measurements > ' + str(measCut) + ' (' + str(veto[3]) + ' vetos)')
        print('\t\t' + str(total-veto[4]) + ' chi squared < ' + str(chi2Cut) + ' (' + str(veto[4]) + ' vetos)')
        print('\t\t' + str(total-veto[5]/2) + ' DOCA < ' + str(docaCut) + ' (' + str(veto[5]/2) + ' vetos)')
        print('\t(NOTE: ' + str(ka_decaycheck) + ' kaons decayed to muons before detection)\n')

        #-------------------------------------PROBABILITIES-----------------------------------
        
        N = (h['tsmearmass_muon'].GetSize()) - 2   # -2 for the underflow and overflow
        xRange = h['tsmearmass_muon'].GetXaxis().GetXmax()   # maximum value of the x-axis
        binsize = xRange/float(N)   # width of histogram bins
        mass = []   # list for particle masses
        prob_mu = N*[0.0]   # lists for probabilities
        prob_ka = N*[0.0]

        for bins in range(0,N):
            mass.append((xRange*float(bins)/float(N)) + binsize)   # list has N entries, last one: mass[N-1]

        emptybins = 0
        for entry,x in enumerate(mass):
            j = h['tsmearmass_muon'].GetXaxis().FindBin(x)   # jth bin
            k = h['tsmearmass_kaon'].GetXaxis().FindBin(x)   # kth bin
            num_mu = h['tsmearmass_muon'].GetBinContent(j)   # gets number of entries in jth bin
            num_ka = h['tsmearmass_kaon'].GetBinContent(k)   # gets number of entries in kth bin
            if num_mu == 0 and num_ka == 0:
                emptybins += 1
                del mass[entry]
                del prob_mu[entry]
                del prob_ka[entry]
            else:
                prob_mu[entry] = (num_mu)/(num_mu + num_ka)   # probability of being a muon
                prob_ka[entry] = (num_ka)/(num_mu + num_ka)   # probability of being a kaon
   
        x = array('f',N*[0.0])
        y_mu = array('f',N*[0.0])
        y_ka = array('f',N*[0.0])
        for i in range(0,N-emptybins):
            x[i] = mass[i]
            y_mu[i] = prob_mu[i]
            y_ka[i] = prob_ka[i]
        
        #fit = ROOT.TF1('fit','1/(1+(e**((-x*[0])+[1])))',0,2)   # sigmoid function
        #fit.SetParameters(20,10)
        fit_mu = ROOT.TF1('fit_mu','pol10',0,1.8)
        fit_mu.SetLineColor(4)
        fit_ka = ROOT.TF1('fit_ka','pol10',0,1.8)
        fit_ka.SetLineColor(2)
         
        Graph_mu = ROOT.TGraph(N,x,y_mu)
        Graph_ka = ROOT.TGraph(N,x,y_ka)
        #Graph_mu.Fit('fit_mu')
        #Graph_ka.Fit('fit_ka')
        combo = ROOT.TMultiGraph()
        combo.Add(Graph_ka)
        combo.Add(Graph_mu)
        h['graph_mu'] = Graph_mu
        h['graph_ka'] = Graph_ka
        h['combo'] = combo

        makePlots()
        
def finStateMuKa_exc():
    if sTree.GetBranch("FitTracks"):
        print('\nRunning analysis for final state K*+- Mu-+ :\n')
        decay1count_kaon = 0   # decay 1: N --> K*+ mu- --> K+ pi0 mu-
        decay1count_pion0 = 0
        decay2count_kaon0 = 0   # decay 2: N --> K*+ mu- --> K0 pi+ mu-
        decay2count_pion = 0
        veto = 6*[0]       
        for n in range(nEvents):   # loops over events
            rc = sTree.GetEntry(n)   # loads tree entry
            daughter_record1 = []
            daughter_record2 = []

            for index,reco_part in enumerate(sTree.FitTracks):   # loops over index and data of track particles                                   
                muPartkey = sTree.fitTrack2MC[index]   # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]   # gives MC particle data
                if abs(true_muon.GetPdgCode()) == 13:   # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()   # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[muonMotherkey]   # obtains mother particle data
                    check,mu_chi2,veto = track_checks(index,true_muon,reco_part,veto)   # performs various checks (i.e. vertex position, fiducial volume,...)
                    if not check == 0: continue
                    if true_mother.GetPdgCode() == 9900015:

                        for index2,reco_part2 in enumerate(sTree.FitTracks):   # loops over index and data of track 
                            Partkey = sTree.fitTrack2MC[index2]   # matches track to MC particle key
                            true_daughter = sTree.MCTrack[Partkey]   # gives MC particle data
                        
                            #===============================================DECAY=1========================================================

                            if abs(true_daughter.GetPdgCode()) == 321:   # checks particle is CHARGED KAON
                                kaonMotherkey = true_daughter.GetMotherId()   # stores number index of MC track of mother
                                true_kaonEx = sTree.MCTrack[kaonMotherkey]   # obtains mother excited kaon data
                                if abs(true_kaonEx.GetPdgCode()) == 323:   # checks mother is excited charged kaon
                                    kaonExMotherkey = true_kaonEx.GetMotherId()           
                                    true_mother_N = sTree.MCTrack[kaonExMotherkey]
                                    if kaonExMotherkey == muonMotherkey:   # checks if mother keys are the same
                                        check2,ka_chi2,veto = track_checks(index2,true_muon,reco_part2,veto)
                                        if not check2 == 0: continue  # performs various checks (i.e. vertex position, fiducial volume,...)
                                        decay1count_kaon += 1
                            #--------------------------------------------------------------------------------------------------------------
                            if true_daughter.GetPdgCode() == 11:   # checks particle is electron
                                if len(daughter_record1) > 0: continue   # to avoid double counting
                                true_photon = sTree.MCTrack[true_daughter.GetMotherId()]
                                if true_photon.GetPdgCode() == 22:   # photon
                                    true_pion = sTree.MCTrack[true_photon.GetMotherId()]
                                    if true_pion.GetPdgCode() == 111:   # neutral pion
                                        true_kaonEx = sTree.MCTrack[true_pion.GetMotherId()]
                                        if abs(true_kaonEx.GetPdgCode()) == 323:   # excited charged kaon
                                            true_N = sTree.MCTrack[true_kaonEx.GetMotherId()]
                                            if kaonExMotherkey == muonMotherkey:   # checks if mother keys are the same
                                                check2,chi2,veto = track_checks(index2,true_muon,reco_part2,veto)
                                                if not check2 == 0: continue   # performs various checks (i.e. vertex position, fiducial volume,...)
                                                #if true_N.GetPdgCode() == 9900015:
                                                decay1count_pion0 += 1
                                                daughter_record1.append(true_daughter_daughter.GetPdgCode())

                            #===============================================DECAY=2========================================================

                            if true_daughter.GetPdgCode() == 211:   # checks particle is charged positive pion (this block needs a bit more work)
                                if len(daughter_record2) > 0: continue   # to avoid double counting
                                pionMotherkey = true_daughter.GetMotherId()   # pionMotherkey != 3 , needs to be from decay of neutral kaon
                                true_kaon = sTree.MCTrack[pionMotherkey]
                                if abs(true_kaon.GetPdgCode()) == 310 or abs(true_kaon.GetPdgCode()) == 130:   # checks mother is NEUTRAL KAON
                                    true_kaon2 = sTree.MCTrack[true_kaon.GetMotherId()]
                                    if true_kaon2.GetPdgCode() == 310 or abs(true_kaon2.GetPdgCode()) == 130:
                                        true_kaonEx = sTree.MCTrack[true_kaon2.GetMotherId()]
                                        if abs(true_kaonEx.GetPdgCode()) == 323:   # checks mother is charged excited kaon   
                                            true_mother_N = sTree.MCTrack[true_kaonEx.GetMotherId()]
                                            #print(muonMotherkey,pionMotherkey,kaon0Motherkey,kaon02Motherkey,kaonExMotherkey)
                                            if kaonExMotherkey == muonMotherkey:   # checks if mother keys are the same
                                                check2,chi2,veto = track_checks(index2,true_muon,reco_part2,veto)
                                                if not check2 == 0: continue  # performs various checks (i.e. vertex position, fiducial volume,...)
                                                decay2count_kaon0 += 1
                                                daughter_record2.append(true_daughter.GetPdgCode())
                            #--------------------------------------------------------------------------------------------------------------
                            if abs(true_daughter.GetPdgCode()) == 211:   # checks particle is CHARGED PION
                                pionMotherkey = true_daughter.GetMotherId()   # pionMotherkey == 3 is what we want here
                                true_kaonEx = sTree.MCTrack[pionMotherkey] 
                                if abs(true_kaonEx.GetPdgCode()) == 323:   # checks mother is excited charged kaon
                                    kaonExMotherkey = true_kaonEx.GetMotherId()  
                                    true_mother_N = sTree.MCTrack[kaonExMotherkey]
                                    if kaonExMotherkey == muonMotherkey:   # checks if mother keys are the same
                                        check2,chi2,veto = track_checks(index2,true_muon,reco_part2,veto)
                                        if not check2 == 0: continue   # performs various checks (i.e. vertex position, fiducial volume,...)
                                        decay2count_pion += 1

                            #==============================================================================================================

        print('\t' + str(decay1count_kaon) + ' charged kaons detected')
        print('\t' + str(decay1count_pion0) + ' neutral pions detected\n')
        print('\t' + str(decay2count_pion) + ' charged pions detected')
        print('\t' + str(decay2count_kaon0) + ' neutral kaons detected\n')

#finStateMuKa_exc()
finStateMuKa()
 
#-------------------------------------------------OUTPUT----------------------------------------------------

hfile = inputFile.split(',')[0].replace('_rec','_RPV')  # Outputs histograms and ROOT file
if hfile[0:4] == "/eos" or not inputFile.find(',')<0:
# do not write to eos, write to local directory 
  tmp = hfile.split('/')
  hfile = tmp[len(tmp)-1] 
ROOT.gROOT.cd()
ut.writeHists(h,hfile)