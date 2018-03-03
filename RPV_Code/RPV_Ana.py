import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
import shipRoot_conf
from array import array
import numpy as np
shipRoot_conf.configure()

debug = False
PDG = ROOT.TDatabasePDG.Instance()
inputFile = None
geoFile = None
dy = None
nEvents = 9999999
fiducialCut = True
chi2Cut = 4
measCut = 25
ecalCut = 0.150
docaCut = 2
ipCut = 250

#--------------------------------------------------INPUT----------------------------------------------------

try:
        opts,args = getopt.getopt(sys.argv[1:],'n:f:g:A:Y:i',['nEvents=','geoFile='])
except getopt.GetoptError:
        # print help information and exit:
        print ' enter file name'
        sys.exit()
for o,a in opts:
        if o in ('-f',):
            inputFile = a
        if o in ('-g','--geoFile',):
            geoFile = a
        if o in ('-Y',):
            dy = float(a)
        if o in ('-n','--nEvents=',):
            nEvents = int(a)

if not inputFile.find(',')<0 :  
  sTree = ROOT.TChain('cbmsim')
  for x in inputFile.split(','):
   if x[0:4] == '/eos':
    sTree.AddFile('root://eoslhcb.cern.ch/' + x)
   else: sTree.AddFile(x)
elif inputFile[0:4] == '/eos':
  eospath = 'root://eoslhcb.cern.ch/' + inputFile
  f = ROOT.TFile.Open(eospath)
  sTree = f.cbmsim
else:
  f = ROOT.TFile(inputFile)
  sTree = f.cbmsim

#-------------------------------------------------GEOMETRY--------------------------------------------------

# try to figure out which ecal geo to load
if not geoFile:
 geoFile = inputFile.replace('ship.','geofile_full.').replace('_rec.','.')
if geoFile[0:4] == '/eos':
  eospath = 'root://eoslhcb.cern.ch/' + geoFile
  fgeo = ROOT.TFile.Open(eospath)
else:  
  fgeo = ROOT.TFile(geoFile)
sGeo = fgeo.FAIRGeom

if not fgeo.FindKey('ShipGeo'):
 # old geofile, missing Shipgeo dictionary
 if sGeo.GetVolume('EcalModule3') :  ecalGeoFile = 'ecal_ellipse6x12m2.geo'
 else: ecalGeoFile = 'ecal_ellipse5x10m2.geo' 
 print 'found ecal geo for ',ecalGeoFile
 # re-create geometry and mag. field
 if not dy:
  # try to extract from input file name
  tmp = inputFile.split('.')
  try:
    dy = float( tmp[1]+'.'+tmp[2] )
  except:
    dy = 10.
 ShipGeo = ConfigRegistry.loadpy('$FAIRSHIP/geometry/geometry_config.py', Yheight = dy, EcalGeoFile = ecalGeoFile )
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
i = 0

for x in ROOT.gGeoManager.GetListOfVolumes():
 volDict[i]=x.GetName()
 i+=1

bfield = ROOT.genfit.BellField(ShipGeo.Bfield.max,ShipGeo.Bfield.z,2,ShipGeo.Yheight/2.)
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
 print 'setup calo reconstruction of ecalReconstructed objects'
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
 ecalCl2Ph=ROOT.TFormula('ecalCl2Ph','[0]+x*([1]+x*([2]+x*[3]))+[4]*x*y+[5]*x*y*y')
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
c = 2.99792458*(10**8)   # speed of light
h = {}   # creates empty dictionary for histograms and graphs

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

def checkFiducialVolume(sTree,tr,dy):
   # extrapolate track to middle of magnet and check if in decay volume
   inside = True
   if not fiducialCut: return True
   fT = sTree.FitTracks[tr]
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
   # u.a - u.c + s*|u|**2 - u.v*t = 0
   # v.a - v.c + s*v.u - t*|v|**2 = 0
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
  for i in range(3): t += tMom(i)/P*(point(i)-tPos(i)) 
  dist = 0
  for i in range(3): dist += (point(i)-tPos(i)-t*tMom(i)/P)**2
  dist = ROOT.TMath.Sqrt(dist)
  return dist

def fitSingleGauss(x,ba=None,be=None):
    name = 'myGauss_' + x 
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
        PosDir[tr] = [xx.getPos(),xx.getDir()]
    xv,yv,zv,doca = myVertex(t1,t2,PosDir)
    # as we have learned, need iterative procedure
    dz = 99999.
    reps,states,newPosDir = {},{},{}
    parallelToZ = ROOT.TVector3(0.,0.,1.)
    rc = True 
    step = 0
    while dz > 0.1:
        zBefore = zv
        newPos = ROOT.TVector3(xv,yv,zv)
        # make a new rep for track 1,2
        for tr in [t1,t2]:     
            xx = sTree.FitTracks[tr].getFittedState()
            reps[tr] = ROOT.genfit.RKTrackRep(xx.getPDG())   # what is rep?
            states[tr] = ROOT.genfit.StateOnPlane(reps[tr])   # what is this?
            reps[tr].setPosMom(states[tr],xx.getPos(),xx.getMom())   # and this?
            try:
                reps[tr].extrapolateToPoint(states[tr],newPos,False)
            except: 
                #print('SHiPAna: extrapolation did not work)'
                rc = False  
                break
            newPosDir[tr] = [reps[tr].getPos(states[tr]),reps[tr].getDir(states[tr])]
        if not rc: break
        xv,yv,zv,doca = myVertex(t1,t2,newPosDir)
        dz = abs(zBefore - zv)   # repeats until dz < 0.1 unless...
        step += 1
        if step > 10:  
            # print 'abort iteration, too many steps, pos=',xv,yv,zv,' doca=',doca,'z before and dz',zBefore,dz
            rc = False
            break
    if not rc: return -1,-1,-1,xv,yv,zv,doca   # extrapolation failed, makes no sense to continue
    LV={}
    for tr in [t1,t2]:   # from here on we have reproduced (see inv_mass() function in Test_1.py)     
        mom = reps[tr].getMom(states[tr])
        pid = abs(states[tr].getPDG()) 
        if pid == 2212: pid = 211
        mass = PDG.GetParticle(pid).Mass()
        E = ROOT.TMath.Sqrt(mass*mass + mom.Mag2())
        LV[tr] = ROOT.TLorentzVector()
        LV[tr].SetPxPyPzE(mom.x(),mom.y(),mom.z(),E)
    RPVMom = LV[t1] + LV[t2]

    return RPVMom,LV[t1],LV[t2],xv,yv,zv,doca

def SingleTrack_4Mom(tr):
    fittedtrack  = sTree.FitTracks[tr].getFittedState()
    Partkey = sTree.fitTrack2MC[tr]   # matches track to MC particle key
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

def time_res(partkey,n,m):
    tnosmear = -1   # declares variable
    true_part = sTree.MCTrack[partkey]   # finds MC particle
    pdg = true_part.GetPdgCode()   # identifies particles from PDG code

    if sTree.GetBranch('strawtubesPoint'):
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

        N = len(z_array)   # number of hits
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
        if abs(pdg) == 13:   # muon
            h['num_muon'].Fill(num_hits)
            for hit in pz_array:
                if n == m:
                    h['track_muon'].Fill(hit)   # muon z-momentum through straw tubes for particular event
        if abs(pdg) == 321:   # kaon
            h['num_kaon'].Fill(num_hits)
            for hit in pz_array:
                if n == m:
                    h['track_kaon'].Fill(hit)   # kaon z-momentum through straw tubes for particular event
  
        if sTree.GetBranch('EcalPoint'):
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
            pdiff = strawP - ecalP   # between 1st straw tube hit and ECAL

            r = ROOT.TMath.Sqrt(((ecal_x - firststraw_x)**2) + ((ecal_y - firststraw_y)**2) + ((ecal_z - firststraw_z)**2))
            h['straight_path'].Fill(r)
            max_z_index = z_array.index(max(z_array))   # gives index of the largest element in the list
            laststraw_x = 0.01*x_array[max_z_index]
            laststraw_y = 0.01*y_array[max_z_index]
            laststraw_z = 0.01*z_array[max_z_index]
            R2 = ROOT.TMath.Sqrt(((ecal_x - laststraw_x)**2) + ((ecal_y - laststraw_y)**2) + ((ecal_z - laststraw_z)**2))
            R = R1+R2   # better approximation of distance travelled through the straw tubes
            h['better_path'].Fill(R)
            rdiff = abs(R - r)
            h['path_diff'].Fill(rdiff)
            
            sigma = 0.01   # standard deviation for Gaussian
            straw_smear = np.random.normal(loc=straw_time,scale=sigma,size=None)
            ecal_smear = np.random.normal(loc=ecal_time,scale=sigma,size=None)
            tsmear = abs(straw_smear - ecal_smear)   # smeared time of flight
            vsmear = (R/tsmear)*(10**9)   # smeared velocity of flight

            tnosmear = abs(straw_time - ecal_time)   # stored in units of nanoseconds
            vnosmear = (R/tnosmear)*(10**9)   # velocity of flight

            if not tnosmear == -1:

                beta = vnosmear/c   # equations for mass calculated from true time
                gamma = 1/(ROOT.TMath.Sqrt(1 - (beta**2)))
                nosmearM = strawP/(beta*gamma)   # previously used reco_muP
            
                beta_smear = vsmear/c   # equations for mass calculated from smeared time                            
                if beta_smear < 1:
                    gamma_smear = 1/(ROOT.TMath.Sqrt(1 - (beta_smear**2)))
                    smearM = strawP/(beta_smear*gamma_smear)
                else: smearM = -1

                if abs(pdg) == 13:
                    h['MuonDir'].Fill(tsmear)   # fills histogram with smeared time
                    h['MuonDir_nosmear'].Fill(tnosmear)   # fills histogram with true time
                    h['tmass_muon'].Fill(nosmearM)   # fills histograms with mass data
                    if not smearM == -1:
                        h['Muon_SmearedMass'].Fill(smearM)
                        h['Total_SmearedMass'].Fill(smearM)
                        h['tsmearmass_muon_samebins'].Fill(smearM)
                if abs(pdg) == 321:
                    h['KaonDir'].Fill(tsmear)   # fills histogram with smeared time
                    h['KaonDir_nosmear'].Fill(tnosmear)   # fills histogram with true time
                    h['tmass_kaon'].Fill(nosmearM)   # fills histograms with mass data
                    if not smearM == -1:
                        h['Kaon_SmearedMass'].Fill(smearM)
                        h['Total_SmearedMass'].Fill(smearM)
                        h['tsmearmass_kaon_samebins'].Fill(smearM)

                h['ecalstraw_mom'].Fill(pdiff)

def track_checks(tr,veto,hist):
    check = 0

    partkey = sTree.fitTrack2MC[tr]   # matches track to MC particle key
    true_part = sTree.MCTrack[partkey]
    reco_part = sTree.FitTracks[tr]

    if not checkFiducialVolume(sTree,tr,dy): 
        #print('Track outside fiducial volume')
        veto[1] += 1
        check = -1

    fit_status = reco_part.getFitStatus()
    if fit_status.isFitConverged():

        fit_nmeas = fit_status.getNdf()
        if hist == 1: h['nmeas'].Fill(fit_nmeas)
        if not fit_nmeas > measCut:
            #print('Too few measurements')
            veto[2] += 1
            check = -1

        fit_rchi2 = fit_status.getChi2()
        fit_chi2 = (fit_rchi2/fit_nmeas)
        if hist == 1: h['Chi2'].Fill(fit_chi2)
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
        if hist == 1: h['ecalE'].Fill(ecal_Etot)
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
                if not muonstation_z[0] == 4017.0:   # z position (cm)
                    veto[5] += 1
                    check = -1

    return check,veto

def createRatio(h1,h2,histname):
    h3 = h1.Clone(histname)
    h3.SetMinimum(0.8)
    h3.SetMaximum(1.35)
    # Set up plot for markers and errors
    h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h1,h2,1,1,'B')   # divide with binomial errors
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetRangeUser(-0.1,1.2)
    y.SetTitleOffset(1.)
    ## Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetRangeUser(0,2)
    return h3

#------------------------------

def createHists_MuKa():
    # Canvas 1
    ut.bookHist(h,'RPV_recomass','Reconstructed Neutralino Mass',500,0.8,1.2)   # reconstructed mass
    ut.bookHist(h,'RPV_truemom','True Neutralino Momentum',100,0.,300.)   # true momentum distribution
    ut.bookHist(h,'RPV_recomom','Reconstructed Neutralino Momentum',100,0.,300)   # reconstructed momentum distribution
    ut.bookHist(h,'RPV_mom_diff','True/Reco Neutralino Momentum Difference',100,-3.,3)   # true/reco momentum difference
    ut.bookHist(h,'RPV_beta','Reconstructed Neutralino Beta in Z-direction',100,0.994,1)
    ut.bookHist(h,'RPV_gamma','Reconstructed Neutralino Gamma in Z-direction',100,0,200)

    # Canvas 2
    ut.bookHist(h,'Muon_truemom','True Muon Momentum',100,0.,140.)   # RPV muon daughter reco momentum
    ut.bookHist(h,'Kaon_truemom','True Kaon Momentum',100,0.,140.)   # RPV pion daughter reco momentum
    ut.bookHist(h,'Muon_recomom','Reconstructed Muon Momentum',100,0.,140.)   # RPV muon daughter true momentum
    ut.bookHist(h,'Kaon_recomom','Reconstructed Kaon Momentum',100,0.,140.)   # RPV pion daughter true momentum
    ut.bookHist(h,'ecalstraw_mom','Straw-Ecal Momentum Difference',500,0,0.4)   # includes both kaon and muon
    ut.bookHist(h,'MuonDir_nosmear','Muon Straw-ECAL Time (directly)',150,37.5,40.)   # daughter muon time of flight
    ut.bookHist(h,'KaonDir_nosmear','Kaon Straw-ECAL Time (directly)',150,37.5,40.)   # daughter kaon time of flight

    # Canvas 3
    ut.bookHist(h,'MuonDir','Smeared Muon Straw-ECAL Time',150,37.5,40.)   # daughter muon time of flight (Gaussian blurred)
    ut.bookHist(h,'KaonDir','Smeared Kaon Straw-ECAL Time',150,37.5,40.)   # daughter kaon time of flight (Gaussian blurred)
    ut.bookHist(h,'tmass_muon','Time Deduced Muon Mass',150,0.,3.)   # time, momentum --> velocity --> gamma (L) --> mass from p=mvL
    ut.bookHist(h,'tmass_kaon','Time Deduced Kaon Mass',150,0.,3.)
    ut.bookHist(h,'tsmearmass_muon_samebins','Smeared Time Deduced Muon Mass',150,0.,3.)   # same as above, but for smeared mass
    ut.bookHist(h,'tsmearmass_kaon_samebins','Smeared Time Deduced Kaon Mass',150,0.,3.)
    ut.bookHist(h,'daughter_masses','Kaon and muon true masses',50,0,1)

    edgesarray = []
    edgesarray.append(0)
    for binNumber in range(0,40):
        edgesarray.append(edgesarray[binNumber] + 0.015)
    for binNumber in range(40,86):
        edgesarray.append(edgesarray[binNumber] + 0.045)
    h['Muon_SmearedMass'] = ROOT.TH1D('Muon_SmearedMass','Muon Smeared Mass',85,array('d',edgesarray))
    h['Kaon_SmearedMass'] = ROOT.TH1D('Kaon_SmearedMass','Kaon Smeared Mass',85,array('d',edgesarray))
    h['Total_SmearedMass'] = ROOT.TH1D('Total_SmearedMass','Smeared Mass',85,array('d',edgesarray))

    # Canvas 4
    h['Muon_ProbMeasr'] = ROOT.TH1D('Muon_ProbMeasr','Probs Identifying Muon',85,array('d',edgesarray))
    h['Kaon_ProbMeasr'] = ROOT.TH1D('Kaon_ProbMeasr','Probs Identifying Kaon',85,array('d',edgesarray))

    # Canvas 5
    ut.bookHist(h,'IP_target','Impact parameter to target',120,0,10)
    ut.bookHist(h,'ecalE','Energy deposited in ECAL',150,0,100)
    ut.bookHist(h,'doca','Distance of closest approach between tracks',150,0,3)
    ut.bookHist(h,'nmeas','No. of measurements in fitted tracks (ndf)',50,0,50)
    ut.bookHist(h,'Chi2','Fitted Tracks Reduced Chi Squared',150,0,3)
    ut.bookHist(h,'recovertex','Reconstructed neutralino decay vertex z-coordinate',100,-4000,4000)

    # Not drawn on canvas
    ut.bookHist(h,'RPV_theta','Angle between neutralino momentum and beam line',100,0,50)
    ut.bookHist(h,'num_muon','No. of muon hits in straw tubes',25,25,50)
    ut.bookHist(h,'num_kaon','No. of kaon hits in straw tubes',25,25,50)
    ut.bookHist(h,'path_diff','Difference between straight path and better approximation',100,0,0.001)
    ut.bookHist(h,'straight_path','Straight distance from first straw tube hit to ECAL',200,11,12)
    ut.bookHist(h,'better_path','Curved path through straw tubes, then straight to ECAL',200,11,12)
    ut.bookHist(h,'track_muon','Muon z-momentum through straw tubes (for particular event)',500,20,21)
    ut.bookHist(h,'track_kaon','Kaon z-momentum through straw tubes (for particular event)',500,44,45)

def makePlots_MuKa():
    ut.bookCanvas(h,key='RPV_N',title='Results 1',nx=1500,ny=800,cx=3,cy=2)
    cv = h['RPV_N'].cd(1)
    h['RPV_recomass'].SetXTitle('Invariant mass / [GeV/c2]')
    h['RPV_recomass'].SetYTitle('No. of Particles')
    h['RPV_recomass'].SetLineColor(1)
    h['RPV_recomass'].Draw()
    print('\nNeutralino mass Gaussian fit:\n')
    fitSingleGauss('RPV_recomass',0.9,1.1)
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['RPV_N'].cd(2)
    h['RPV_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['RPV_truemom'].SetYTitle('No. of Particles')
    h['RPV_truemom'].SetLineColor(2)
    h['ths1'] = ROOT.THStack('RPVmom','True & Reconstructed Neutralino Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths1'].Add(h['RPV_truemom'])
    h['ths1'].Add(h['RPV_recomom'])
    h['ths1'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['RPV_N'].cd(3)
    h['RPV_mom_diff'].SetXTitle('Momentum Difference / [GeV/c]')
    h['RPV_mom_diff'].SetYTitle('No. of Particles')
    h['RPV_mom_diff'].SetLineColor(1)    
    h['RPV_mom_diff'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['RPV_N'].cd(4)
    h['RPV_beta'].SetXTitle('Beta')
    h['RPV_beta'].SetYTitle('No. of Particles')
    h['RPV_beta'].SetLineColor(1)
    h['RPV_beta'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['RPV_N'].cd(5)
    h['RPV_gamma'].SetXTitle('Gamma')
    h['RPV_gamma'].SetYTitle('No. of Particles')
    h['RPV_gamma'].SetLineColor(1)
    h['RPV_gamma'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['RPV_N'].cd(6)
    h['RPV_theta'].SetXTitle('Theta / [mrad]')
    h['RPV_theta'].SetYTitle('No. of Particles')
    h['RPV_theta'].SetLineColor(1)
    h['RPV_theta'].Draw()
    h['RPV_N'].Print('RPV_N.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='KaMu',title='Results 2',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['KaMu'].cd(1)
    h['Kaon_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['Kaon_truemom'].SetYTitle('No. of Particles')
    h['Kaon_truemom'].SetLineColor(2)
    h['ths2'] = ROOT.THStack('Kamom','Kaon True & Reconstructed Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths2'].Add(h['Kaon_truemom'])
    h['ths2'].Add(h['Kaon_recomom'])
    h['ths2'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['KaMu'].cd(2)
    h['Muon_truemom'].SetXTitle('Momentum / [GeV/c]')
    h['Muon_truemom'].SetYTitle('No. of Particles')
    h['Muon_truemom'].SetLineColor(2)
    h['ths3'] = ROOT.THStack('Mumom','Muon True & Reconstructed Momentum ; Momentum / [GeV/c] ; No. of Particles')
    h['ths3'].Add(h['Muon_truemom'])
    h['ths3'].Add(h['Muon_recomom'])
    h['ths3'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['KaMu'].cd(3)
    h['ecalstraw_mom'].SetXTitle('Momentum Difference [GeV/c]')
    h['ecalstraw_mom'].SetYTitle('Frequency')
    h['ecalstraw_mom'].SetLineColor(1)
    h['ecalstraw_mom'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['KaMu'].cd(4)
    h['MuonDir_nosmear'].SetXTitle('Time / [ns]')
    h['MuonDir_nosmear'].SetYTitle('No. of particles')
    h['KaonDir_nosmear'].SetLineColor(2)
    h['ths4'] = ROOT.THStack('MuKaDirTime','Kaon & Muon Straw-ECAL Time ; Time / [ns] ; No. of Particles')
    h['ths4'].Add(h['MuonDir_nosmear'])
    h['ths4'].Add(h['KaonDir_nosmear'])
    h['ths4'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    h['KaMu'].Print('KaMu.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='Time_Mass',title='Results 3',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['Time_Mass'].cd(1)
    h['MuonDir'].SetXTitle('Time / [ns]')
    h['MuonDir'].SetYTitle('Frequency')
    h['MuonDir'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Time_Mass'].cd(2)
    h['KaonDir'].SetXTitle('Time / [ns]')
    h['KaonDir'].SetYTitle('Frequency')
    h['KaonDir'].SetLineColor(2)
    h['KaonDir'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Time_Mass'].cd(3)
    h['tmass_kaon'].SetXTitle('Mass / [GeV/c2]')
    h['tmass_kaon'].SetYTitle('Frequency')
    h['tmass_kaon'].SetLineColor(2)
    h['daughter_masses'].SetLineColor(1)
    h['daughter_masses'].SetLineStyle(2)
    h['daughter_masses'].Draw('same')
    h['ths5'] = ROOT.THStack('tmass','Time Deduced Kaon & Muon Mass ; Mass / [GeV/c2] ; Frequency')
    h['ths5'].Add(h['tmass_kaon'])
    h['ths5'].Add(h['tmass_muon'])
    #h['ths5'].Add(h['daughter_masses'])
    h['ths5'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Time_Mass'].cd(4)
    h['tsmearmass_kaon_samebins'].SetXTitle('Mass / [GeV/c2]')
    h['tsmearmass_kaon_samebins'].SetYTitle('Frequency')
    h['tsmearmass_kaon_samebins'].SetLineColor(2)
    h['tsmearmass_muon_samebins'].SetXTitle('Mass / [GeV/c2]')
    h['tsmearmass_muon_samebins'].SetYTitle('Frequency')
    print('\nLandau fits for mass (time of flight):\n')
    h['tsmearmass_kaon_samebins'].Fit('landau')
    h['tsmearmass_kaon_samebins'].GetFunction('landau').SetLineColor(1)
    h['tsmearmass_muon_samebins'].Fit('landau')
    h['tsmearmass_muon_samebins'].GetFunction('landau').SetLineColor(1)
    #par0 = h['Muon_SmearedMass'].GetFunction('landau').GetParameter(0)
    #par1 = h['Muon_SmearedMass'].GetFunction('landau').GetParameter(1)
    #par2 = h['Muon_SmearedMass'].GetFunction('landau').GetParameter(2)
    h['ths6'] = ROOT.THStack('smeartmass','Smeared Time Deduced Kaon & Muon Mass ; Mass / [GeV/c2] ; Frequency')
    h['ths6'].Add(h['tsmearmass_kaon_samebins'])
    h['ths6'].Add(h['tsmearmass_muon_samebins'])
    #h['ths6'].Add(h['daughter_masses'])
    h['ths6'].Draw('nostack')
    ROOT.gPad.BuildLegend()
    h['Time_Mass'].Print('Time_Mass.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='Probs',title='Results 4',nx=1000,ny=1000,cx=2,cy=2)
    cv = h['Probs'].cd(1)
    h['Muon_ProbMeasr'].SetXTitle('Mass / [GeV/c2]')
    h['Muon_ProbMeasr'].SetYTitle('Probability')
    h['Muon_ProbMeasr'].SetLineColor(4)
    h['Muon_ProbMeasr'].SetMarkerColor(4)
    h['Muon_ProbMeasr'].SetMarkerStyle(33)
    h['Muon_ProbMeasr'].SetMarkerSize(1)
    print('\nPolynomial fits for probability graphs:')
    h['Muon_ProbMeasr'].Fit('pol4')
    h['Muon_ProbMeasr'].GetFunction('pol4').SetLineColor(4)
    h['Muon_ProbMeasr'].Draw('E')
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Probs'].cd(2) 
    h['Kaon_ProbMeasr'].SetXTitle('Mass / [GeV/c2]')
    h['Kaon_ProbMeasr'].SetYTitle('Probability')
    h['Kaon_ProbMeasr'].SetLineColor(2)
    h['Kaon_ProbMeasr'].SetMarkerColor(2)
    h['Kaon_ProbMeasr'].SetMarkerStyle(33)
    h['Kaon_ProbMeasr'].SetMarkerSize(1)
    h['Kaon_ProbMeasr'].Fit('pol4')
    h['Kaon_ProbMeasr'].GetFunction('pol4').SetLineColor(2)
    h['Kaon_ProbMeasr'].Draw('E')
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Probs'].cd(3)
    h['Kaon_ProbMeasr'].Draw('E')
    h['Muon_ProbMeasr'].Draw('E same')
    h['Probs'].Print('Probs.png')
    #======================================================================================================================
    ut.bookCanvas(h,key='Vetos',title='Results 5',nx=1500,ny=800,cx=3,cy=2)
    cv = h['Vetos'].cd(1)
    h['IP_target'].SetXTitle('Impact Parameter / [cm]')
    h['IP_target'].SetYTitle('Frequency')
    h['IP_target'].SetLineColor(1)
    h['IP_target'].SetFillColor(17)
    h['IP_target'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Vetos'].cd(2)
    h['ecalE'].SetXTitle('Energy / [GeV/c2]')
    h['ecalE'].SetYTitle('Frequency')
    h['ecalE'].SetLineColor(1)
    h['ecalE'].SetFillColor(17)
    h['ecalE'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Vetos'].cd(3)
    h['doca'].SetXTitle('Distance / [cm]')
    h['doca'].SetYTitle('Frequency')
    h['doca'].SetLineColor(1)
    h['doca'].SetFillColor(17)
    h['doca'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Vetos'].cd(4)
    h['nmeas'].SetXTitle('ndf')
    h['nmeas'].SetYTitle('No. of tracks')
    h['nmeas'].SetLineColor(1)
    h['nmeas'].SetFillColor(17)
    h['nmeas'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Vetos'].cd(5)
    h['Chi2'].SetXTitle('Reduced Chi Squared')
    h['Chi2'].SetYTitle('Frequency')
    h['Chi2'].SetLineColor(1)
    h['Chi2'].SetFillColor(17)
    h['Chi2'].Draw()
    #----------------------------------------------------------------------------------------------------------------------
    cv = h['Vetos'].cd(6)
    h['recovertex'].SetXTitle('Z / [cm]')
    h['recovertex'].SetYTitle('Frequency')
    h['recovertex'].SetLineColor(1)
    h['recovertex'].SetFillColor(17)
    h['recovertex'].Draw()
    h['Vetos'].Print('Vetos.png')
    #======================================================================================================================
    h['straight_path'].SetLineColor(1)
    h['better_path'].SetLineColor(1)
    h['track_kaon'].SetLineColor(2)
   
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

#----------------------------------------------------EVENT-LOOPS--------------------------------------------------------

def print_menu(): 
    print('\n' + 30*'-' + 'MENU' + 30*'-')
    print('1. RPV SUSY Benchmark 1 --> K+- mu+- fully visible final state')
    print('2. RPV SUSY Benchmark 1 --> K*+- mu+- partially visible final state')
    print('0. Exit')
    print(64 * '-' + '\n')

nEvents = min(sTree.GetEntries(),nEvents)   # number of generated events

def finStateMuKa():
    if sTree.GetBranch('FitTracks'):
        print('\nRunning analysis for final state K+- Mu-+ :\n')
        createHists_MuKa()   # calls function to create histograms
        successful_events = []   # creates list of event numbers of desired decays
        veto = 10*[0]   # creates list of veto counts for each possible veto cause
        ka_veto = 10*[0]   # creates list of veto counts for muons which decayed from kaons from neutralinos
        ka_decaycheck = 0   # variable for counting when kaons from neutralinos decay to muons before detection
        for n in range(nEvents):   # loops over events
            rc = sTree.GetEntry(n)   # loads tree entry
            event = True

            #-----------------------------------------------TRACK-LOOPS------------------------------------------------

            for index,reco_part in enumerate(sTree.FitTracks):   # loops over index and data of particle tracks                                 
                muPartkey = sTree.fitTrack2MC[index]   # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]   # gives MC particle data                      

                if abs(true_muon.GetPdgCode()) == 13:   # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()   # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[muonMotherkey]   # obtains mother particle data

                    if true_mother.GetPdgCode() == 321:
                        fitstatus = reco_part.getFitStatus()
                        if fitstatus.isFitConverged(): 
                            motherN = sTree.MCTrack[true_mother.GetMotherId()]
                            if motherN.GetPdgCode() == 9900015:
                                ka_decaycheck += 1   # kaon has decayed to a muon in flight
                                check,ka_veto = track_checks(index,ka_veto,0)
                                if check == -1: ka_veto[0] += 1 
                                
                    if true_mother.GetPdgCode() == 9900015:   # checks mother is RPV  
                        for index2,reco_part2 in enumerate(sTree.FitTracks):   # loops over index and data of track particles
                            kaPartkey = sTree.fitTrack2MC[index2]   # matches track to MC particle key
                            true_kaon = sTree.MCTrack[kaPartkey]   # gives MC particle data
                            if abs(true_kaon.GetPdgCode()) == 321:   # checks particle is kaon
                                kaonMotherkey = true_kaon.GetMotherId()   # stores a number index of MC track of mother
                                true_mother = sTree.MCTrack[kaonMotherkey]   # obtains mother particle data

                                if kaonMotherkey == muonMotherkey and true_mother.GetPdgCode() == 9900015:   # check if mother keys are the same
                                    fitstatus_muon = reco_part.getFitStatus() 
                                    fitstatus_kaon = reco_part2.getFitStatus()
                                    if fitstatus_muon.isFitConverged() and fitstatus_kaon.isFitConverged():
                                        veto[0] += 1
                                    else:
                                        #print('At least one of the track fits did not converge')
                                        continue
                                    
                                    check,veto = track_checks(index,veto,1)   # performs various track checks
                                    if check == -1: event = False
                                    else:   # if the first track was fine, check the other one (ensures event veto counting is correct)
                                        check,veto = track_checks(index2,veto,1)
                                        if check == -1: event = False
                                    
                                    #-------------------------------------------------PARTICLE-DATA------------------------------------------------------

                                    RPV_4Mom,Muon_4Mom,Kaon_4Mom,NDecay_X,NDecay_Y,NDecay_Z,doca = RedoVertexing(index,index2)   # uses RedoVertexing to iterate track fitting
                                    if RPV_4Mom == -1: 
                                        print('RedoVertexing extrapolation failed (event ' + str(n) + ')')
                                        veto[9] += 1
                                        break

                                    h['recovertex'].Fill(NDecay_Z)
                                    if not isInFiducial(NDecay_X,NDecay_Y,NDecay_Z):
                                        #print('Neutralino decayed outside fiducial volume')
                                        veto[6] += 1
                                        event = False
                                    
                                    h['doca'].Fill(doca)
                                    if not doca < docaCut: 
                                        #print('Distance of closest approach too large')
                                        veto[7] +=1
                                        event = False

                                    NProduction_X = true_mother.GetStartX()
                                    NProduction_Y = true_mother.GetStartY()
                                    NProduction_Z = true_mother.GetStartZ()
                                    NProduction_T = true_mother.GetStartT()
                                    RPV_Pos = ROOT.TLorentzVector()
                                    RPV_Pos.SetXYZT(NProduction_X,NProduction_Y,NProduction_Z,NProduction_T)
                                    tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
                                    ip = ImpactParameter(tr,RPV_Pos,RPV_4Mom)   # gives the same result as line 706 in ShipAna.py (i.e. using sTree.Particles)

                                    #NDecay_T = true_muon.GetStartT()
                                    #RPV_Pos = ROOT.TLorentzVector()
                                    #RPV_Pos.SetXYZT(NDecay_X,NDecay_Y,NDecay_Z,NDecay_T)
                                    #tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
                                    #ip = ImpactParameter(tr,RPV_Pos,RPV_4Mom)

                                    h['IP_target'].Fill(ip)
                                    if not ip < ipCut:
                                        #print('Neutralino impact parameter to target too large')
                                        veto[8] += 1
                                        event = False

                                    if event == False: continue

                                    RPV_truemass = true_mother.GetMass()   # RPV neutralino mass
                                    RPV_truemom = true_mother.GetP()   # RPV neutralino momentum
                                    RPV_recomass = RPV_4Mom.M()   # reconstructed RPV mass
                                    RPV_recomom = RPV_4Mom.P()   # reconstructed RPV momentum
                                    RPV_momdiff = RPV_truemom - RPV_recomom   # RPV true/reco momentum difference
                                    true_kaP = true_kaon.GetP()   # true kaon momentum
                                    reco_kaP = Kaon_4Mom.P()   # reconstructed kaon momentum
                                    true_muP = true_muon.GetP()   # true muon momentum
                                    reco_muP = Muon_4Mom.P()   # reconstructed muon momentum
                                    kaM = true_kaon.GetMass()   # kaon mass
                                    muM = true_muon.GetMass()   # muon mass
                                                      
                                    h['RPV_truemom'].Fill(RPV_truemom)   # fills histograms
                                    h['RPV_recomass'].Fill(RPV_recomass)                        
                                    h['RPV_recomom'].Fill(RPV_recomom)                            
                                    h['RPV_mom_diff'].Fill(RPV_momdiff)
                                    h['Kaon_recomom'].Fill(reco_kaP)
                                    h['Muon_recomom'].Fill(reco_muP)
                                    h['Kaon_truemom'].Fill(true_kaP)
                                    h['Muon_truemom'].Fill(true_muP)
                                    h['daughter_masses'].Fill(kaM)
                                    h['daughter_masses'].Fill(muM)

                                    RPV_Zmom = RPV_4Mom.Pz()
                                    RPV_Zbeta = 1/ROOT.TMath.Sqrt(1 + ((RPV_recomass/RPV_Zmom)**2))
                                    RPV_Zgamma = 1/(ROOT.TMath.Sqrt(1 - (RPV_Zbeta**2)))
                                    cos_theta = RPV_Zmom/RPV_recomom
                                    theta = 1000*(ROOT.TMath.ACos(cos_theta))   # angle between beam line and neutralino momentum (mrad)

                                    h['RPV_beta'].Fill(RPV_Zbeta)
                                    h['RPV_gamma'].Fill(RPV_Zgamma)
                                    h['RPV_theta'].Fill(theta)
                                    
                                    successful_events.append(n)   # adds entries to the list
                                    m = successful_events[0]   # arbitrarily picks the first one as an example
                                    
                                    #------------------------------------TIME-RESOLUTION------------------------------------------

                                    time_res(muPartkey,n,m)   # calculates particle time of flight and smeared mass        
                                      
                                    time_res(kaPartkey,n,m)   # calculates particle time of flight and smeared mass      

        #----------------------------------------------------------------VETO-COUNTS------------------------------------------------------------------

        accepted = len(successful_events)
        rejected = veto[0] - accepted
        print('\n\t' + str(veto[0]) + ' events reconstructed for this decay mode')
        print('\t' + str(accepted) + ' events not rejected')
        print('\t' + str(rejected) + ' events rejected:')
        print('\t\t' + str(veto[0] - veto[6]) + ' events with neutralino decay vertex inside fiducial volume (' + str(veto[6]) + ' event vetos)')
        print('\t\t' + str(veto[0] - veto[1]) + ' events with both tracks within fiducial volume (' + str(veto[1]) + ' event vetos)')
        print('\t\t' + str(veto[0] - veto[2]) + ' events with both tracks no. of measurements > ' + str(measCut) + ' (' + str(veto[2]) + ' event vetos)')
        print('\t\t' + str(veto[0] - veto[3]) + ' events with both tracks chi squared < ' + str(chi2Cut) + ' (' + str(veto[3]) + ' event vetos)')
        print('\t\t' + str(veto[0] - veto[4]) + ' events where both tracks leave > ' + str(ecalCut) + ' GeV in ECAL (' + str(veto[4]) + ' event vetos)')
        print('\t\t' + str(veto[0] - veto[5]) + ' events with hits in both the 1st and 2nd muon stations (' + str(veto[5]) + ' event vetos)')
        print('\t\t' + str(veto[0] - veto[7]) + ' events with DOCA between tracks < ' + str(docaCut) + ' cm (' + str(veto[7]) + ' event vetos)')
        print('\t\t' + str(veto[0] - veto[8]) + ' events with IP to target < ' + str(ipCut) + ' cm (' + str(veto[8]) + ' event vetos)')
        print('\t' + str(veto[9]) + ' failed RedoVertexing extrapolations')
        print('\t' + str(ka_decaycheck) + ' kaons decayed to muons before detection (' + str(ka_decaycheck-ka_veto[0]) + ' after track checks)')

        #-------------------------------------PROBABILITIES-----------------------------------
        
        h['Muon_ProbMeasr'] = createRatio(h['Muon_SmearedMass'],h['Total_SmearedMass'],'Muon_ProbMeasr')
        h['Kaon_ProbMeasr'] = createRatio(h['Kaon_SmearedMass'],h['Total_SmearedMass'],'Kaon_ProbMeasr')
        
        del h['Muon_SmearedMass']   # deletes histograms which are no longer needed
        del h['Kaon_SmearedMass']
        del h['Total_SmearedMass']

        makePlots_MuKa()

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

                                                                            if event == True: decayBcount_kaon0 += 1

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

loop = True
while loop:
    print_menu()
    choice = input('Enter your choice [0-2]: ')
     
    if choice == 1:
        finStateMuKa()
        #-----------------------------------------OUTPUT------------------------------------------
        hfile = inputFile.split(',')[0].replace('_rec','_RPV')  # Outputs histograms and ROOT file
        if hfile[0:4] == '/eos' or not inputFile.find(',')<0:
        # do not write to eos, write to local directory 
          tmp = hfile.split('/')
          hfile = tmp[len(tmp)-1] 
        ROOT.gROOT.cd()
        ut.writeHists(h,hfile)
        loop = False

    elif choice == 2:
        finStateMuKa_exc()
        #-------------------------------------------OUTPUT--------------------------------------------
        hfile = inputFile.split(',')[0].replace('_rec','_RPV_Exc')  # Outputs histograms and ROOT file
        if hfile[0:4] == '/eos' or not inputFile.find(',')<0:
        # do not write to eos, write to local directory 
          tmp = hfile.split('/')
          hfile = tmp[len(tmp)-1] 
        ROOT.gROOT.cd()
        ut.writeHists(h,hfile)
        loop = False

    elif choice == 0:
        print('\nExit selected')
        loop = False

    else: print('\nNo analysis option selected. Try again:')