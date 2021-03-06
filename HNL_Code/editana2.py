# #####################################
#         EDITA Test Code             #
#######################################

# example for accessing smeared hits and fitted tracks
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
shipRoot_conf.configure()

debug = False
chi2CutOff  = 4.
PDG = ROOT.TDatabasePDG.Instance()
fiducialCut = True
measCutFK = 25
measCutPR = 22
docaCut = 2.
lineLogFile=open('EDITA_lineLog.txt','w')
particleDataFile=open('partData.txt','w')
lineLogFile.write('Log of Active Tested Line numbers: \n')
particleDataFile.write('Particle Data File: \n')

def get_linenumber():
    curframe = currentframe()
    return curframe.f_back.f_lineno
def LineActivity(ns_variable,ns_line):     
    if ns_variable==ns_line:
        lineLogFile.write(str(ns_line) + '\n')
        ns_variable+=1
def inputOptsArgs():
    inputFile  = None
    geoFile    = None
    dy         = None
    nEvents    = 9999999
    try:
            opts, args = getopt.getopt(sys.argv[1:], "n:f:g:Y", ["nEvents=","geoFile="])#command line options whose API is designed to be familiar to users of the C getopt() function, consider using argparse
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

#SETS TREE FROM INPUT FILE
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

# TRIES TO FIGURE OUT WHICH ECAL GEOMETRY TO LOAD
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

# CREATE GEOMETRY
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

import TrackExtrapolateTool

### FUNCTIONS RELEVANT TO WHAT WE NEED ###
h = {}
def create_Hists():
    ut.bookHist(h,'HNL_true','Monte Carlo Mass',500,0.,2.)                      # true mass
    ut.bookHist(h,'HNL_reco','Reconstructed Mass',500,0.,2.)                    # reconstructed mass
    ut.bookHist(h,'HNL_no_iter','Reconstructed Mass (without track iterations)',500,0.,2.)
    ut.bookHist(h,'HNL_mom','True (red) & Reco. (blue) Momentum',100,0.,300.)   # true momentum distribution
    ut.bookHist(h,'HNL_mom_reco','Reconstructed Momentum',100,0.,300)           # reconstructed momentum distribution
    ut.bookHist(h,'HNL_mom_diff','True/Reco Momentum Difference',100,-3.,3)     # true/reco momentum difference

    ut.bookHist(h,'Time','Muon Straw-ECAL Time (directly)',500,36.,40.)         # muon daughter time of flight
    ut.bookHist(h,'Time2','Pion Straw-ECAL Time (directly)',500,36.,40.)        # pion daughter time of flight
    ut.bookHist(h,'Time3','Muon Straw-ECAL Time (indirectly)',500,36.,40.)      # muon daughter time of flight
    ut.bookHist(h,'Time4','Pion Straw-ECAL Time (indirectly)',500,36.,40.)      # pion daughter time of flight

    ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared',100,0.,3.)                 # chi squared track fitting

    ut.bookHist(h,'Muon_mom','Muon (HNL Daughter) Momentum',100,0.,200.)        # HNL muon daughter momentum
    ut.bookHist(h,'Pion_mom','Pion (HNL Daughter) Momentum',100,0.,200.)        # HNL pion daughter momentum
    print("Created 13 Histograms")
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

def isInFiducial(X,Y,Z):
   if Z > ShipGeo.TrackStation1.z : return False
   if Z < ShipGeo.vetoStation.z+100.*u.cm : return False
   # typical x,y Vx resolution for exclusive HNL decays 0.3cm,0.15cm (gaussian width)
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
     if not rc: return xv,yv,zv,doca,-1 # extrapolation failed, makes no sense to continue
     LV={}
     for tr in [t1,t2]: # from here on we have reproduced (see inv_mass()) 
         mom = reps[tr].getMom(states[tr])
         pid = abs(states[tr].getPDG()) 
         if pid == 2212: pid = 211 #why
         mass = PDG.GetParticle(pid).Mass()
         E = ROOT.TMath.Sqrt( mass*mass + mom.Mag2() )
         LV[tr] = ROOT.TLorentzVector()
         LV[tr].SetPxPyPzE(mom.x(),mom.y(),mom.z(),E)
     HNLMom = LV[t1]+LV[t2]
     #return xv,yv,zv,doca,HNLMom
     return doca,HNLMom

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
   ut.bookCanvas(h,key='Test_Mass',title='Fit Results',nx=1000,ny=1000,cx=2,cy=2)
   cv = h['Test_Mass'].cd(1)
   h['HNL_true'].SetXTitle('Invariant mass [GeV/c2]')
   h['HNL_true'].SetYTitle('No. of Particles')
   h['HNL_true'].Draw()
   
   cv = h['Test_Mass'].cd(2)
   h['HNL_reco'].SetXTitle('Invariant mass [GeV/c2]')
   h['HNL_reco'].SetYTitle('No. of Particles')
   h['HNL_reco'].Draw()
   fitSingleGauss('HNL_reco',0.9,1.1)
   
   cv = h['Test_Mass'].cd(3)
   h['HNL_mom'].SetXTitle('Momentum [GeV/c]')
   h['HNL_mom'].SetYTitle('No. of Particles')
   h['HNL_mom'].SetLineColor(2)
   h['HNL_mom'].Draw()
   h['HNL_mom_reco'].Draw("same")
   
   cv = h['Test_Mass'].cd(4)
   h['HNL_mom_diff'].SetXTitle('Momentum Difference [GeV/c]')
   h['HNL_mom_diff'].SetYTitle('Frequency')
   h['HNL_mom_diff'].Draw()
   h['Test_Mass'].Print('HNLgraphs.png')
   #======================================================================================================================
   ut.bookCanvas(h,key='Time_Res',title='Fit Results 2',nx=1000,ny=1000,cx=2,cy=2)
   cv = h['Time_Res'].cd(1)
   h['Time'].SetXTitle('Time [ns]')
   h['Time'].SetYTitle('Frequency')
   h['Time'].Draw()
   
   cv = h['Time_Res'].cd(2)
   h['Time2'].SetXTitle('Time [ns]')
   h['Time2'].SetYTitle('Frequency')
   h['Time2'].Draw()
   
   cv = h['Time_Res'].cd(3)
   h['Time3'].SetXTitle('Time [ns]')
   h['Time3'].SetYTitle('Frequency')
   h['Time3'].Draw()
   
   cv = h['Time_Res'].cd(4)
   h['Time4'].SetXTitle('Time [ns]')
   h['Time4'].SetYTitle('Frequency')
   h['Time4'].Draw()
   h['Time_Res'].Print('TimeRes1.png')

nEvents = min(sTree.GetEntries(),nEvents)

def time_res(partkey,v):
    t1 = None
    t2 = None
    if sTree.GetBranch("strawtubesPoint"):
        x_array = []
        y_array = []
        z_array = []
        t_array = []
        straw_time = 0
        for k,hits in enumerate(sTree.strawtubesPoint):
            straw_TrackID = hits.GetTrackID()
            if straw_TrackID == partkey:
                x_array.append(hits.GetX())
                y_array.append(hits.GetY())
                z_array.append(hits.GetZ())
                t_array.append(hits.GetTime())
        
        min_z_index = z_array.index(min(z_array))
        straw_z = 0.01*min(z_array)
        straw_x = 0.01*x_array[min_z_index]
        straw_y = 0.01*y_array[min_z_index]
        straw_time = t_array[min_z_index]

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

        if not ecal_time <= 0:
            t1 = abs(straw_time - ecal_time)
            r = ROOT.TMath.Sqrt(((ecal_x - straw_x)**2) + ((ecal_y - straw_y)**2) + ((ecal_z - straw_z)**2))
            t2 = (r/v)*(10**9) # units of nanoseconds
                            
    return t1,t2

def finState2MuPi():
    if sTree.GetBranch("FitTracks"):
        pi_decaycheck = 0
        for n in range(nEvents):                            # loop over events
            rc = sTree.GetEntry(n)                              # load tree entry
            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles                                   
                muPartkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]                  # gives MC particle data
                if abs(true_muon.GetPdgCode()) == 13:               # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()             # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[muonMotherkey]          # obtains mother particle data
                    if true_mother.GetPdgCode() == 211:
                        print('Pion has decayed to a muon')
                        pi_decaycheck+=1
                    if true_mother.GetPdgCode() == 9900015:             # checks mother is HNL

                        Decay_X = true_muon.GetStartX()
                        Decay_Y = true_muon.GetStartY()
                        Decay_Z = true_muon.GetStartZ()
                        if not isInFiducial(Decay_X,Decay_Y,Decay_Z):
                            #print('HNL decayed outside fiducial volume')
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

                        fittedstate1 = reco_part.getFittedState()           # get reconstructed muon fitted state
                        mu_M = true_muon.GetMass()                          # mass of MC muon
                        muPx = fittedstate1.getMom().x()                    # momentum in x
                        muPy = fittedstate1.getMom().y()                    # momentum in y  
                        muPz = fittedstate1.getMom().z()                    # momentum in z
                        muP = fittedstate1.getMomMag()                      # momentum magnitude
                        muE = ROOT.TMath.Sqrt((mu_M**2) + (muP**2))         # energy

                        Muon_Vector = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                        Muon_Vector.SetPxPyPzE(muPx,muPy,muPz,muE)          # inputs 4-vector elements

                        for index2,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles
                            piPartkey = sTree.fitTrack2MC[index2]                  # matches track to MC particle key
                            true_pion = sTree.MCTrack[piPartkey]                  # gives MC particle data
                            if abs(true_pion.GetPdgCode()) == 211:              # checks particle is pion
                                pionMotherkey = true_pion.GetMotherId()             # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[pionMotherkey]          # obtains mother particle data
                                if pionMotherkey==muonMotherkey:                    # check if keys are the same
                                    pionMotherTrue_mass = true_mother.GetMass()         # get HNL/final states mother mass
                                    pionMotherTrue_mom = true_mother.GetP()             # get HNL/final states mother mom

                                    if not checkFiducialVolume(sTree,index,dy): 
                                        #print('Decay outside fiducial volume')
                                        continue 
                                    pi_status = reco_part.getFitStatus()                
                                    if not pi_status.isFitConverged():
                                        #print('Fit did not converge')
                                        continue
                                    pi_nmeas = pi_status.getNdf() 
                                    if not pi_nmeas > 25:
                                        #print('Too few measurements')
                                        continue

                                    pi_rchi2 = pi_status.getChi2()                      # chi squared value
                                    pi_chi2 = (pi_rchi2/pi_nmeas)                       # gets chi value

                                    if not pi_chi2 < 4:
                                        #print('Chi squared value too high')
                                        continue

                                    fittedstate2 = reco_part.getFittedState()           # get reconstructed pion fitted state
                                    pi_M = true_pion.GetMass()                          # mass of MC pion 
                                    piP = fittedstate2.getMomMag()                      # momentum in x
                                    piPx = fittedstate2.getMom().x()                    # momentum in y
                                    piPy = fittedstate2.getMom().y()                    # momentum in z
                                    piPz = fittedstate2.getMom().z()                    # momentum magnitude
                                    piE = ROOT.TMath.Sqrt((pi_M**2) + (piP**2))         # energy

                                    piV = (3*(10**8)*piP) / ROOT.TMath.Sqrt((pi_M**2) + (piP**2))       # pion velocity
                                    muV = (3*(10**8)*muP) / ROOT.TMath.Sqrt((mu_M**2) + (muP**2))       # muon velocity

                                    Pion_Vector = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                    Pion_Vector.SetPxPyPzE(piPx,piPy,piPz,piE)          # inputs 4-vector elements
                              
                                    HNL_Vector = Muon_Vector + Pion_Vector              # adds the 4-momenta
                                    hnlmass = HNL_Vector.M() # mass without iterating over track fitting
                                    h['HNL_no_iter'].Fill(hnlmass)
                                    
                                    HNLMom_redo = ROOT.TLorentzVector()
                                    doca,HNLMom_Redo = RedoVertexing(index,index2) # uses RedoVertexing to iterate track fitting
                                    if HNLMom_Redo == -1: continue
                                    if doca > 2.: 
                                        #print('distance of closest approach too large')
                                        continue

                                    HNL_mass = HNLMom_Redo.M()                           # sets HNL mass
                                    HNL_reco_mom = HNLMom_Redo.P()                       # sets HNL mom
                                    mom_diff = pionMotherTrue_mom - HNL_reco_mom
                                    
                                    h['HNL_true'].Fill(pionMotherTrue_mass)             # fill histograms 
                                    h['HNL_mom'].Fill(pionMotherTrue_mom)
                                    h['HNL_reco'].Fill(HNL_mass)                        
                                    h['HNL_mom_reco'].Fill(HNL_reco_mom)                
                                    h['Chi2'].Fill(mu_chi2)       
                                    h['Chi2'].Fill(pi_chi2)                             
                                    h['HNL_mom_diff'].Fill(mom_diff)
                                    h['Pion_mom'].Fill(piP)
                                    h['Muon_mom'].Fill(muP)

                                    mu_t1,mu_t2 = time_res(muPartkey,muV)        
                                    if mu_t1 != None:              
                                        h['Time'].Fill(mu_t1) 
                                        h['Time3'].Fill(mu_t2)
                                    pi_t1,pi_t2 = time_res(piPartkey,piV)                            
                                    if pi_t1 != None:     
                                        h['Time2'].Fill(pi_t1)  
                                        h['Time4'].Fill(pi_t2)

        print('\n'+str(pi_decaycheck) + ' pi --> mu decays before detection\n')


finState2MuPi()          
makePlots()
hfile = inputFile.split(',')[0].replace('_rec','_NStesting')#create outputFile
if hfile[0:4] == "/eos" or not inputFile.find(',')<0:
# do not write to eos, write to local directory 
  tmp = hfile.split('/')
  hfile = tmp[len(tmp)-1]                                   #occurs only for cern users
ROOT.gROOT.cd()
ut.writeHists(h,hfile)                                      #write histograms to outputFile
lineLogFile.close()                                         #close lineTestLogFile
particleDataFile.close()                                    #close File