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
chi2CutOff  = 4.
measCutFK = 25
measCutPR = 22
docaCut = 2.

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

# ------------------------------------------------GEOMETRY--------------------------------------------------

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
from ROOT import THStack

#----------------------------------------------------HISTOGRAMS-----------------------------------------------------------

h = {}
ut.bookHist(h,'RPV_true','Monte Carlo Mass',500,0.8,1.2) # true mass
ut.bookHist(h,'RPV_reco','Reconstructed Mass',500,0.8,1.2) # reconstructed mass
ut.bookHist(h,'RPV_mom','True (red) & Reco. (blue) Momentum',100,0.,300.) # true momentum distribution
ut.bookHist(h,'RPV_mom_reco','Reconstructed Momentum',100,0.,300) # reconstructed momentum distribution
ut.bookHist(h,'RPV_mom_diff','True/Reco Momentum Difference',100,-3.,3) # true/reco momentum difference

ut.bookHist(h,'Muon_mom_true','Muon - True (red) & Reco. (blue) Momentum',100,0.,140.) # RPV muon daughter reco momentum
ut.bookHist(h,'Kaon_mom_true','Kaon - True (red) & Reco. (blue) Momentum',100,0.,140.) # RPV pion daughter reco momentum
ut.bookHist(h,'Muon_mom','Muon - True Momentum',100,0.,140.) # RPV muon daughter true momentum
ut.bookHist(h,'Kaon_mom','Kaon - True Momentum',100,0.,140.) # RPV pion daughter true momentum

ut.bookHist(h,'MuonDir','Muon Straw-ECAL Time (directly)',500,37.5,40.) # muon daughter time of flight
ut.bookHist(h,'KaonDir','Kaon Straw-ECAL Time (directly)',500,37.5,40.) # kaon daughter time of flight
ut.bookHist(h,'smearedmass_muon','Time Smeared Muon Mass',100,0.,2.)
ut.bookHist(h,'smearedmass_kaon','Time Smeared Kaon(red)-Muon(blue) Mass',100,0.,2.)
ut.bookHist(h,'ecalstraw_mom','Straw-Ecal Momentum Difference',500,0,0.4)

ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared',100,0.,3.) # chi squared track fitting
ut.bookHist(h,'num_muon','No. of muon hits in straw tubes',25,25,50)
ut.bookHist(h,'num_kaon','No. of kaon in straw tubes',25,25,50)
ut.bookHist(h,'track_muon','Muon z-momentum decrease through straw tubes',100,0,100)
ut.bookHist(h,'track_kaon','Kaon z-momentum decrease through straw tubes',100,0,100)

# ---------------------------------------------------FUNCTIONS------------------------------------------------------------

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
   Z = c.z()+v.z()*t
   return X,Y,Z,abs(dist) # X,Y,Z are the coordinates of...?

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
        #print 'SHiPAna: extrapolation did not work'
        rc = False  
        break # breaks if extrapolation doesn't work

       newPosDir[tr] = [reps[tr].getPos(states[tr]),reps[tr].getDir(states[tr])]
      if not rc: break
      xv,yv,zv,doca = myVertex(t1,t2,newPosDir)
      dz = abs(zBefore-zv) # repeats until dz < 0.1 unless...
      step+=1
      if step > 10:  
         print 'abort iteration, too many steps, pos=',xv,yv,zv,' doca=',doca,'z before and dz',zBefore,dz
         rc = False
         break # ... number of iterations exceeds 10
     if not rc: return -1,-1,-1,doca # extrapolation failed, makes no sense to continue
     LV={}
     for tr in [t1,t2]: # from here on we have reproduced (see inv_mass() function)       
      mom = reps[tr].getMom(states[tr])
      pid = abs(states[tr].getPDG()) 
      if pid == 2212: pid = 211 # why
      mass = PDG.GetParticle(pid).Mass()
      E = ROOT.TMath.Sqrt(mass*mass + mom.Mag2())
      LV[tr] = ROOT.TLorentzVector()
      LV[tr].SetPxPyPzE(mom.x(),mom.y(),mom.z(),E)
     RPVMom = LV[t1]+LV[t2]

     return RPVMom,LV[t1],LV[t2],doca

def makePlots():
   ut.bookCanvas(h,key='Test_Time',title='Fit Results',nx=1000,ny=1000,cx=2,cy=2)
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
   h['smearedmass_kaon'].SetXTitle('Mass [GeV/c2]')
   h['smearedmass_kaon'].SetYTitle('Frequency')
   h['smearedmass_kaon'].SetLineColor(2)
   h['smearedmass_kaon'].Draw()
   h['smearedmass_muon'].Draw("same")
   #----------------------------------------------------------------------------------------------------------------------
   cv = h['Test_Time'].cd(4)
   h['ecalstraw_mom'].SetXTitle('Momentum Difference [GeV/c2]')
   h['ecalstraw_mom'].SetYTitle('Frequency')
   h['ecalstraw_mom'].Draw()
   h['ecalstraw_mom'].SetLineColor(1)
   h['Test_Time'].Print('Smeared_Time.png')
   #======================================================================================================================
   ut.bookCanvas(h,key='Test_Mass',title='Fit Results 2',nx=1000,ny=1000,cx=2,cy=2)
   cv = h['Test_Mass'].cd(1)
   h['RPV_true'].SetXTitle('Invariant mass [GeV/c2]')
   h['RPV_true'].SetYTitle('No. of Particles')
   h['RPV_true'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   cv = h['Test_Mass'].cd(2)
   h['RPV_reco'].SetXTitle('Invariant mass [GeV/c2]')
   h['RPV_reco'].SetYTitle('No. of Particles')
   h['RPV_reco'].Draw()
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
   h['Test_Mass'].Print('RPV_Graphs.png')
   #======================================================================================================================
   ut.bookCanvas(h,key='KaMu_Graphs',title='Fit Results 3',nx=1000,ny=500,cx=2,cy=1)
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
   h['KaMu_Graphs'].Print('KaMu_Graphs.png')

def time_res(partkey,pdg):
    tsmear = -1
    vsmear = -1
    diff = -1
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
        
        num_hits = len(z_array)
        if pdg==13:
            h['num_muon'].Fill(num_hits)
            for hit in pz_array:
                h['track_muon'].Fill(hit)
        if pdg==321:
            h['num_kaon'].Fill(num_hits)
            for hit in pz_array:
                h['track_kaon'].Fill(hit)

        min_z_index = z_array.index(min(z_array))
        straw_z = 0.01*z_array[min_z_index]
        straw_x = 0.01*x_array[min_z_index]
        straw_y = 0.01*y_array[min_z_index]
        straw_time = t_array[min_z_index]
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
                        ecal_x = 0.01*hits.GetX()           # stored in units of cm 
                        ecal_y = 0.01*hits.GetY()
                        ecal_z = 0.01*hits.GetZ()
                        ecal_time = hits.GetTime()
                        ecalPx = hits.GetPx()
                        ecalPy = hits.GetPy()
                        ecalPz = hits.GetPz()
                        ecalP = ROOT.TMath.Sqrt((ecalPx**2) + (ecalPy**2) + (ecalPz**2)) # straw tube momentum

        if not ecal_time <= 0:
            #straw_smear = np.random.normal(loc=straw_time,scale=0.1,size=None)
            #ecal_smear = np.random.normal(loc=ecal_time,scale=0.1,size=None)
            diff = strawP - ecalP
            r = ROOT.TMath.Sqrt(((ecal_x - straw_x)**2) + ((ecal_y - straw_y)**2) + ((ecal_z - straw_z)**2))
            tsmear = abs(straw_time - ecal_time)            # stored in units of nanoseconds
            vsmear = (r/tsmear)*(10**9)
            
    return tsmear,vsmear,diff

# ---------------------------------------------------EVENT-LOOP-----------------------------------------------------------

nEvents = min(sTree.GetEntries(),nEvents)
c = 2.99792458*(10**8)

def finStateMuKa():
    if sTree.GetBranch("FitTracks"):
        k_decaycheck = 0
        for n in range(nEvents):                            # loop over events
            rc = sTree.GetEntry(n)                              # load tree entry
            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles                                   
                muPartkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]                  # gives MC particle data
                if abs(true_muon.GetPdgCode()) == 13:                   # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()             # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[muonMotherkey]          # obtains mother particle data
                    if true_mother.GetPdgCode() == 211:
                        #print('Kaon has decayed to a muon before detection')
                        k_decaycheck+=1
                    if true_mother.GetPdgCode() == 9900015:             # checks mother is RPV

                        Decay_X = true_muon.GetStartX()
                        Decay_Y = true_muon.GetStartY()
                        Decay_Z = true_muon.GetStartZ()
                        if not isInFiducial(Decay_X,Decay_Y,Decay_Z):
                            #print('RPV decayed outside fiducial volume')
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
                            kaPartkey = sTree.fitTrack2MC[index2]                  # matches track to MC particle key
                            true_kaon = sTree.MCTrack[kaPartkey]                  # gives MC particle data
                            if abs(true_kaon.GetPdgCode()) == 321:              # checks particle is kaon
                                kaonMotherkey = true_kaon.GetMotherId()             # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[kaonMotherkey]          # obtains mother particle data
                                if kaonMotherkey==muonMotherkey:                    # check if keys are the same

                                    kaonMotherTrue_mass = true_mother.GetMass()         # get RPV/final states mother mass
                                    kaonMotherTrue_mom = true_mother.GetP()             # get RPV/final states mother mom

                                    if not checkFiducialVolume(sTree,index,dy): 
                                        #print('Decay outside fiducial volume')
                                        continue 
                                    ka_status = reco_part.getFitStatus()                
                                    if not ka_status.isFitConverged():
                                        #print('Fit did not converge')
                                        continue
                                    ka_nmeas = ka_status.getNdf() 
                                    if not ka_nmeas > 25:
                                        #print('Too few measurements')
                                        continue
                                    ka_rchi2 = ka_status.getChi2()                      # chi squared value
                                    ka_chi2 = (ka_rchi2/ka_nmeas)                       # gets chi value
                                    if not ka_chi2 < 4:
                                        #print('Chi squared value too high')
                                        continue

                                    RPV_Vector = ROOT.TLorentzVector()                  # declares variables as TLorentzVector class
                                    Muon_Vector = ROOT.TLorentzVector()
                                    Kaon_Vector = ROOT.TLorentzVector()
                                    RPV_Vector,Muon_Vector,Kaon_Vector,doca = RedoVertexing(index,index2)  # uses RedoVertexing to iterate track fitting
                                    
                                    if RPV_Vector == -1: continue
                                    if doca > 2.: 
                                        #print('distance of closest approach too large')
                                        continue

                                    RPV_mass = RPV_Vector.M()                           # sets RPV mass
                                    RPV_reco_mom = RPV_Vector.P()                       # sets RPV mom
                                    mom_diff = kaonMotherTrue_mom - RPV_reco_mom

                                    true_kaP = true_kaon.GetP()
                                    reco_kaP = Kaon_Vector.P()
                                    true_muP = true_muon.GetP()
                                    reco_muP = Muon_Vector.P()
                                                      
                                    h['RPV_true'].Fill(kaonMotherTrue_mass) 
                                    h['RPV_mom'].Fill(kaonMotherTrue_mom)
                                    h['RPV_reco'].Fill(RPV_mass)                        
                                    h['RPV_mom_reco'].Fill(RPV_reco_mom)                
                                    h['Chi2'].Fill(mu_chi2)       
                                    h['Chi2'].Fill(ka_chi2)                             
                                    h['RPV_mom_diff'].Fill(mom_diff)
                                    h['Kaon_mom'].Fill(reco_kaP)
                                    h['Muon_mom'].Fill(reco_muP)
                                    h['Kaon_mom_true'].Fill(true_kaP)
                                    h['Muon_mom_true'].Fill(true_muP)

                                    mu_t,mu_v,mu_diff = time_res(muPartkey,13)        
                                    if mu_t != -1: # and mu_t < 38.08:
                                        h['MuonDir'].Fill(mu_t) 
                                        beta = mu_v/c
                                        gamma = 1/(ROOT.TMath.Sqrt(1-(beta**2)))
                                        smearedM = reco_muP/(beta*gamma)
                                        h['smearedmass_muon'].Fill(smearedM)
                                        ka_t,ka_v,ka_diff = time_res(kaPartkey,321)      
                                        if ka_t != -1: # and ka_t < 38.08:
                                            h['KaonDir'].Fill(ka_t)
                                            beta = ka_v/c
                                            gamma = 1/(ROOT.TMath.Sqrt(1-(beta**2)))
                                            smearedM = reco_kaP/(beta*gamma)
                                            h['smearedmass_kaon'].Fill(smearedM)
                                            h['ecalstraw_mom'].Fill(mu_diff)
                                            h['ecalstraw_mom'].Fill(ka_diff)

        print('\n'+str(k_decaycheck) + ' muons with K+ mothers\n')

finStateMuKa()  
makePlots()

# ---------------------------------------------------OUTPUT------------------------------------------------------------

# Outputs histograms and ROOT file
hfile = inputFile.split(',')[0].replace('_rec','_RPV')
if hfile[0:4] == "/eos" or not inputFile.find(',')<0:
# do not write to eos, write to local directory 
  tmp = hfile.split('/')
  hfile = tmp[len(tmp)-1] 
ROOT.gROOT.cd()
ut.writeHists(h,hfile)
