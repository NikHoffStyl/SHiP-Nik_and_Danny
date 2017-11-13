import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
import shipRoot_conf
shipRoot_conf.configure()

debug = False
chi2CutOff  = 4.
PDG = ROOT.TDatabasePDG.Instance()
inputFile  = None
geoFile    = None
dy         = None
nEvents    = 9999999
fiducialCut = True
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

#----------------------------------------------------HISTOGRAMS-----------------------------------------------------------

h = {}
ut.bookHist(h,'HNL_true','Monte Carlo Mass',500,0.,2.) # true mass
ut.bookHist(h,'HNL_reco','Reconstructed Mass',500,0.,2.) # reconstructed mass
ut.bookHist(h,'HNL_mom','True & Reconstructed Momentum Distribution',100,0.,300.) # true momentum distribution
ut.bookHist(h,'HNL_mom_reco','Reconstructed Momentum',100,0.,300) # reconstructed momentum distribution
ut.bookHist(h,'HNL_mom_diff','True/Reco Momentum Difference',100,-3.,3) # true/reco momentum difference

ut.bookHist(h,'Time','Muon - Time Between Straw Tube and ECAL',100,0.,200.) # muon daughter time of flight
ut.bookHist(h,'Time2','Pion - Time Between Straw Tube and ECAL',100,0.,200.) # pion daughter time of flight
ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared',100,0.,3.) # chi squared track fitting

ut.bookHist(h,'Muon_mom','Muon (HNL Daughter) Momentum',100,0.,200.) # HNL muon daughter momentum
ut.bookHist(h,'Pion_mom','Pion (HNL Daughter) Momentum',100,0.,200.) # HNL pion daughter momentum


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

def checkFiducialVolume(sTree,tkey,dy):
# extrapolate track to middle of magnet and check if in decay volume
   inside = True
   if not fiducialCut: return True
   fT = sTree.FitTracks[tkey]
   rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,ShipGeo.Bfield.z)
   if not rc: return False
   if not dist2InnerWall(pos.X(),pos.Y(),pos.Z())>0: return False
   return inside

def match2HNL(p):
    matched = False
    hnlKey  = []
    for t in [p.GetDaughter(0),p.GetDaughter(1)]: 
      mcp = sTree.fitTrack2MC[t]
      while mcp > -0.5:
        mo = sTree.MCTrack[mcp]
        if abs(mo.GetPdgCode()) == 9900015:
           hnlKey.append(mcp)
           break  
        mcp = mo.GetMotherId()
    if len(hnlKey) == 2: 
       if hnlKey[0]==hnlKey[1]: matched = True
    return matched

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

def time_res(partkey):
    if sTree.GetBranch("strawtubesPoint"):
        z_array = []
        t_array = []
        straw_time=0
        for k,hits in enumerate(sTree.strawtubesPoint):
            TrackID = hits.GetTrackID()
            #print(TrackID)
            if TrackID == partkey:
                z_array.append(hits.GetZ())
                t_array.append(hits.GetTime())
        minZval=min(z_array)   
        min_z = z_array.index(min(z_array))
        straw_time = t_array[min_z]
    else: return None,None,None,None,None
    if sTree.GetBranch("EcalPoint"):
            if not straw_time<=0:
                for k,hits in enumerate(sTree.EcalPoint):
                    TrackID = hits.GetTrackID()
                    if TrackID == partkey:
                        ecal_time = hits.GetTime()
                        ecal_zpos = hits.GetZ()
                        if not ecal_time <= straw_time:
                            t = abs(straw_time - ecal_time)
                            return ecal_time,ecal_zpos, minZval, straw_time, t 
                    else: return None,None,None,None,None
            else: return None,None,None,None,None
    else: return None,None,None,None,None

def time_res2(partkey):
    if sTree.GetBranch("strawtubesPoint"):
        x_array = []
        y_array = []
        z_array = []
        t_array = []
        for k,hits in enumerate(sTree.strawtubesPoint):
            TrackID = hits.GetTrackID()
            #print(TrackID)
            if TrackID == partkey:
                x_array.append(hits.GetX())
                y_array.append(hits.GetY())
                z_array.append(hits.GetZ())
                t_array.append(hits.GetTime())
        
        min_z_index = z_array.index(min(z_array))
        straw_zpos = min(z_array)
        straw_xpos = x_array[min_z_index]
        straw_ypos = y_array[min_z_index]
        straw_time = t_array[min_z_index]
        
    else: return None

    if sTree.GetBranch("EcalPoint"):
            if not straw_time<=0:
                for k,hits in enumerate(sTree.EcalPoint):
                    TrackID = hits.GetTrackID()
                    if TrackID == partkey:
                        ecal_xpos = hits.GetX()
                        ecal_ypos = hits.GetY()
                        ecal_zpos = hits.GetZ()
                        ecal_time = hits.GetTime()

                        if not ecal_time <= straw_time:
                            t = abs(straw_time - ecal_time)
                            return t

            else: return None
    else: return None

def checkVetoHit(partkey):
    if sTree.GetBranch("vetoPoint"):
        veto = False
        for hits in sTree.vetoPoint:
            ID = hits.GetTrackID()
            if ID == partkey:
                veto = True
        return veto

def makePlots():
   ut.bookCanvas(h,key='Test_Mass',title='Fit Results',nx=1000,ny=1000,cx=2,cy=2)
   cv = h['Test_Mass'].cd(1)
   h['HNL_true'].SetXTitle('Invariant mass [GeV/c2]')
   h['HNL_true'].SetYTitle('No. of Particles')
   h['HNL_true'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   cv = h['Test_Mass'].cd(2)
   h['HNL_reco'].SetXTitle('Invariant mass [GeV/c2]')
   h['HNL_reco'].SetYTitle('No. of Particles')
   h['HNL_reco'].Draw()
   fitSingleGauss('HNL_reco',0.9,1.1)
   #----------------------------------------------------------------------------------------------------------------------
   cv = h['Test_Mass'].cd(3)
   h['HNL_mom'].SetXTitle('Momentum [GeV/c]')
   h['HNL_mom'].SetYTitle('No. of Particles')
   h['HNL_mom'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   #cv = h['Test_Mass'].cd(4)
   #h['HNL_mom_reco'].SetXTitle('Momentum [GeV/c]')
   #h['HNL_mom_reco'].SetYTitle('No. of Particles')
   h['HNL_mom_reco'].Draw("same")
   cv = h['Test_Mass'].cd(4)
   h['HNL_mom_diff'].SetXTitle('Momentum Difference [GeV/c]')
   h['HNL_mom_diff'].SetYTitle('Frequency')
   h['HNL_mom_diff'].Draw()
   h['Test_Mass'].Print('HNL_Graphs.png')
   #----------------------------------------------------------------------------------------------------------------------
   ut.bookCanvas(h,key='Time_Res',title='Fit Results 2',nx=1000,ny=1000,cx=2,cy=2)
   cv = h['Time_Res'].cd(1)
   h['Time'].SetXTitle('Time [ns]')
   h['Time'].SetYTitle('Frequency')
   h['Time'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   cv = h['Time_Res'].cd(2)
   h['Time2'].SetXTitle('Time [ns]')
   h['Time2'].SetYTitle('Frequency')
   h['Time2'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   cv = h['Time_Res'].cd(3)
   h['Muon_mom'].SetXTitle('Momentum [GeV/c]')
   h['Muon_mom'].SetYTitle('No. of Particles')
   h['Muon_mom'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   cv = h['Time_Res'].cd(4)
   h['Pion_mom'].SetXTitle('Momentum [GeV/c]')
   h['Pion_mom'].SetYTitle('No. of Particles')
   h['Pion_mom'].Draw()
   h['Time_Res'].Print('Time_Res.png')

# ---------------------------------------------------EVENT-LOOP-----------------------------------------------------------

nEvents = min(sTree.GetEntries(),nEvents)
import TrackExtrapolateTool
from array import array



def finStateMuPi():
    if sTree.GetBranch("FitTracks"):
        pi_decaycheck = 0                               # counter for pions decaying to muons before detection
        fiducialcheck = 0                               # counter for HNL decays outside ficucial volume
        convergecheck = 0                               # counter for failed track fits
        for n in range(nEvents):                            # loop over events
            rc = sTree.GetEntry(n)                              # load tree entry
            keylist = []                                          # create empty list
            muVector = {}                                   # create empty dictionaries
            dicMuChi2 = {}
            mupartkey = {}
            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles
                if not checkFiducialVolume(sTree,index,dy): 
                    print('Decay outside fiducial volume')
                    fiducialcheck+=1
                    continue                                    # skips to next iteration if HNL decayed in fiducial volume 
                muPartkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]                  # gives MC particle data
                #if checkVetoHit(muPartkey):
                    #print('Hit in veto tagger')
                    #continue
                if abs(true_muon.GetPdgCode()) == 13:               # checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()             # stores a number index of MC track of mother
                    keylist.append(muonMotherkey)                     # adds to list
                    true_mother = sTree.MCTrack[muonMotherkey]          # obtains mother particle data
                    if true_mother.GetPdgCode() == 9900015:             # checks mother is HNL
                        mu_status = reco_part.getFitStatus()                # gets fit status
                        if not mu_status.isFitConverged():
                            print('Fit did not converge')
                            convergecheck+=1
                            continue
                        mu_rchi2 = mu_status.getChi2()                      # gets chi squared value
                        mu_nmeas = mu_status.getNdf()                       # gets number of measurements
                        mu_chi2 = (mu_rchi2/mu_nmeas)                       # gets chi value
                        dicMuChi2[str(muonMotherkey)] = mu_chi2

                        mupartkey[str(muonMotherkey)] = muPartkey

                        fittedstate1 = reco_part.getFittedState()           # get reconstructed muon fitted state
                        mu_M = true_muon.GetMass()                          # mass of MC muon
                        muPx = fittedstate1.getMom().x()                    # momentum in x
                        muPy = fittedstate1.getMom().y()                    # momentum in y  
                        muPz = fittedstate1.getMom().z()                    # momentum in z
                        muP = fittedstate1.getMomMag()                      # momentum magnitude
                        muE = ROOT.TMath.Sqrt((mu_M**2) + (muP**2))         # energy

                        Muon_Vector = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                        Muon_Vector.SetPxPyPzE(muPx,muPy,muPz,muE)          # inputs four-momentum elements
                        muVector[str(muonMotherkey)]=Muon_Vector

                        muonMotherTrue_mass = true_mother.GetMass()         # gets HNL mass

                    if true_mother.GetPdgCode() == 211:             # checks if mother is pion not HNL
                        print('Pion has decayed to a muon')
                        pi_decaycheck+=1

            for index,reco_part in enumerate(sTree.FitTracks):  # loops over index and data of track particles
                piPartkey = sTree.fitTrack2MC[index]                  # matches track to MC particle key
                true_pion = sTree.MCTrack[piPartkey]                  # gives MC particle data
                if abs(true_pion.GetPdgCode()) == 211:              # checks particle is pion
                    pionMotherkey = true_pion.GetMotherId()             # stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[pionMotherkey]          # obtains mother particle data
                    for muonMotherkey in keylist:                       # loops through muonMother keys
                        if pionMotherkey==muonMotherkey:                    # check if keys are the same
                            #if true_mother.GetPdgCode() == 9900015:
                                pionMotherTrue_mass = true_mother.GetMass()         # get HNL/final states mother mass
                                pionMotherTrue_mom = true_mother.GetP()             # get HNL/final states mother mom
                                pi_status = reco_part.getFitStatus()                # gets fit status
                                if not pi_status.isFitConverged():
                                    print('Fit did not converge')
                                    convergecheck+=1
                                    continue
                                pi_rchi2 = pi_status.getChi2()                      # chi squared value
                                pi_nmeas = pi_status.getNdf()                       # gets number of measurements
                                pi_chi2 = (pi_rchi2/pi_nmeas)                       # gets chi value

                                fittedstate2 = reco_part.getFittedState()           # get reconstructed pion fitted state

                                pi_M = true_pion.GetMass()                          # mass of MC pion
                                piP = fittedstate2.getMomMag()                      # momentum in x
                                piPx = fittedstate2.getMom().x()                    # momentum in y
                                piPy = fittedstate2.getMom().y()                    # momentum in z
                                piPz = fittedstate2.getMom().z()                    # momentum magnitude
                                piE = ROOT.TMath.Sqrt((pi_M**2) + (piP**2))         # energy

                                Pion_Vector = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                                Pion_Vector.SetPxPyPzE(piPx,piPy,piPz,piE)          # inputs four-momentum elements

                                h['HNL_true'].Fill(pionMotherTrue_mass)             # true HNL mass
                                h['HNL_mom'].Fill(pionMotherTrue_mom)               # true HNL momentum
                                h['Pion_mom'].Fill(piP)
                                h['Muon_mom'].Fill(muP)
                            
                                HNL_Vector = muVector[str(muonMotherkey)] + Pion_Vector # adds the 4-momenta
                                HNL_mass = HNL_Vector.M()                           # sets HNL mass
                                HNL_reco_mom = HNL_Vector.P()                       # sets HNL mom
                                mom_diff = pionMotherTrue_mom - HNL_reco_mom

                                h['HNL_reco'].Fill(HNL_mass)                        # fill histograms
                                h['HNL_mom_reco'].Fill(HNL_reco_mom)                #----||-------
                                h['Chi2'].Fill(dicMuChi2[str(muonMotherkey)])       #----||-------
                                h['Chi2'].Fill(pi_chi2)                             #----||-------
                                h['HNL_mom_diff'].Fill(mom_diff)                    #----||-------
                                
                                mu_t = time_res2(mupartkey[str(muonMotherkey)])      
                                if mu_t != None:                                  
                                    h['Time'].Fill(mu_t)                                

                                pi_t = time_res2(piPartkey)                            
                                if pi_t != None:                                    
                                    h['Time2'].Fill(pi_t)      

        print('\n'+str(pi_decaycheck) + ' pi --> mu decays before detection')
        print(str(fiducialcheck) + ' HNL decays outside fiducial volume')
        print(str(convergecheck) + ' track fits failed to converge\n')
        
finStateMuPi()                       
makePlots()
print('finished creating plots')

# ---------------------------------------------------OUTPUT------------------------------------------------------------

# Outputs histograms and ROOT file
hfile = inputFile.split(',')[0].replace('_rec','_HNL')
if hfile[0:4] == "/eos" or not inputFile.find(',')<0:
# do not write to eos, write to local directory 
  tmp = hfile.split('/')
  hfile = tmp[len(tmp)-1] 
ROOT.gROOT.cd()
ut.writeHists(h,hfile)