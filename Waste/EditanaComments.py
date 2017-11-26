############################################################
#   #           #    #       ######   #######   ######     #
#    #    #    #   #   #      #          #      #          #
#     #  # #  #   #######       #        #      ####       #
#      #     #   #       #        #      #      #          #
#      #     #  #         #   ######     #      ######     #
############################################################
#for HNL in sTree.Particles:
    #t1 = HNL.GetDaughter(0) 

#MC_id = sTree.MCTrack.PdgCode()

#if abs(sTree.MCTrack[hnlkey].GetPdgCode()) == 9900015:

#mom = reps[tr].getMom(states[tr])
#pid = abs(states[tr].getPDG()) 
#if pid == 2212: pid = 211
#mass = PDG.GetParticle(pid).Mass()

#mo = sTree.MCTrack[mcp]
#if abs(mo.GetPdgCode()) == 9900015:

#pdgcode = fT.getFittedState().getPDG()
#tmp = PDG.GetParticle(pdgcode)

#tmp = PDG.GetParticle(aP.GetPdgCode())

#for hnlkey in [1,2]: 
#if abs(sTree.MCTrack[hnlkey].GetPdgCode()) == 9900015:

#idMother = abs(sTree.MCTrack[hnlkey-1].GetPdgCode())

#for k, rec_particle in enumerate(sTree.FitTracks):
#    if rec_particle.GetPdgCode()== 13 || 211: # if its a muon or a pion
#        part_key=sTree.fitTrack2MC[k]
#        rec_particle=sTree.MCTrack[part_key] #gives particle of track
#        rec_part_motherkey=rec_particle.GetMotherId()
#        rec_part_mother=sTree.MCTrack(rec_part_motherkey)
#        if rec_part_mother.GetPdgCode()== 9900015:
            #print ("Do something")

#if sTree.GetBranch("FitTracks"):
#    print('found branch FitTracks')
#    for n in range(nEvents):
#        for k, reco_part in enumerate(sTree.FitTracks):
#            if reco_part.GetPdgCode() == 13:# or reco_part.GetPdgCode() == 211:
#                print('found particle')
#                partkey = sTree.fitTrack2MC[k]
#                reco_part = sTree.MCTrack[partkey] # gives particle of track
#                motherkey = reco_part.GetMotherId() # stores the id of the mother
#                reco_mother = sTree.MCTrack[motherkey] # retrieves mother particle using id
#                if reco_mother.GetPdgCode() == 9900015:
#                    print('found mother of particle')
#if sTree.GetBranch("FitTracks"):
#    for n in range(nEvents):
#        rc = sTree.GetEntry(n)
#        for index,reco_part in enumerate(sTree.FitTracks):
#            partkey = sTree.fitTrack2MC[index]
#            true_part = sTree.MCTrack[partkey] # gives particle of track
#            if abs(true_part.GetPdgCode()) == 13 or abs(true_part.GetPdgCode()) == 211:
#                motherkey = true_part.GetMotherId() # stores the id of the mother
#                true_mother = sTree.MCTrack[motherkey] # retrieves mother particle using id
#                if true_mother.GetPdgCode() == 9900015:
#                    muonMotherTrue_mass = true_mother.GetMass()
#                    h['HNL_true'].Fill(muonMotherTrue_mass)

#if sTree.GetBranch("MCTrack"):
#    LineActivity(get_linenumber(), get_linenumber()) #doesnt do this
#    for n in range(nEvents):
#        for mc_particle in sTree.MCTrack:
#            #print(mc_particle)
#            if mc_particle.GetPdgCode() == 9900015:
#                inv_mass = mc_particle.GetMass()
#                h['HNL_reco'].Fill(inv_mass)


#def time_res():
#    if sTree.GetBranch("EcalPoint"):
#        print('found branch EcalPoint')
#        if sTree.GetBranch("strawtubesPoint"):
#            print('found branch strawtubespoint')
#            ut.bookHist(h,'time_res','Time Resolution Test',500,0.,2.)
#            for n in range(nEvents):
#                rc = sTree.GetEntry(n)
#                k=0
#                for ahit in sTree.EcalPoint:
#                    if k==0:
#                        t1 = ahit.GetTime()
#                        ecalID = ahit.GetTrackID()
#                        LineActivity(get_linenumber()+n,get_linenumber)
#                        for ahit in sTree.strawtubesPoint:
#                            if k==0:
#                                t2 = ahit.GetTime()
#                                strawID = ahit.GetTrackID()
#                                if strawID == ecalID:
#                                    time = abs(t2-t1)
#                                    h['time_res'].Fill(time)
#                                    k+=1
#time_res()




### FUNCTIONS USED IN EVENT LOOP FUNCTION ###
def HNLKinematics():
 HNLorigin={}
 ut.bookHist(h,'HNLmomNoW','momentum unweighted',100,0.,300.)
 ut.bookHist(h,'HNLmom','momentum',100,0.,300.)
 ut.bookHist(h,'HNLPtNoW','Pt unweighted',100,0.,10.)
 ut.bookHist(h,'HNLPt','Pt',100,0.,10.)
 ut.bookHist(h,'HNLmom_recTracks','momentum',100,0.,300.)
 ut.bookHist(h,'HNLmomNoW_recTracks','momentum unweighted',100,0.,300.)
 ut.bookHist(h,'HNL_sim','Simulated Mass',500,0.,2.) # new one for the simulated (Monte Carlo) mass

 for n in range(sTree.GetEntries()): 

  rc = sTree.GetEntry(n)
  for hnlkey in [1,2]: 
   if abs(sTree.MCTrack[hnlkey].GetPdgCode()) == 9900015: 
    theHNL = sTree.MCTrack[hnlkey]
    wg = theHNL.GetWeight()
    if not wg>0.: wg=1.
    idMother = abs(sTree.MCTrack[hnlkey-1].GetPdgCode())
    if not HNLorigin.has_key(idMother): HNLorigin[idMother]=0
    HNLorigin[idMother]+=wg
    P = theHNL.GetP()
    Pt = theHNL.GetPt()
    h['HNLmom'].Fill(P,wg) 
    h['HNLmomNoW'].Fill(P) 
    h['HNLPt'].Fill(Pt,wg) 
    h['HNLPtNoW'].Fill(Pt) 
    for HNL in sTree.Particles:
     t1,t2 = HNL.GetDaughter(0),HNL.GetDaughter(1) 
     for tr in [t1,t2]:
      xx  = sTree.FitTracks[tr].getFittedState()
      Prec = xx.getMom().Mag()
      h['HNLmom_recTracks'].Fill(Prec,wg) 
      h['HNLmomNoW_recTracks'].Fill(Prec)
 theSum = 0
 for x in HNLorigin: theSum+=HNLorigin[x]   
 for x in HNLorigin: print "%4i : %5.4F relative fraction: %5.4F "%(x,HNLorigin[x],HNLorigin[x]/theSum)

def VertexError(t1,t2,PosDir,CovMat,scalFac):
# with improved Vx x,y resolution
   a,u = PosDir[t1]['position'],PosDir[t1]['direction']
   c,v = PosDir[t2]['position'],PosDir[t2]['direction']
   Vsq = v.Dot(v)
   Usq = u.Dot(u)
   UV  = u.Dot(v)
   ca  = c-a
   denom = Usq*Vsq-UV**2
   tmp2 = Vsq*u-UV*v
   Va = ca.Dot(tmp2)/denom
   tmp2 = UV*u-Usq*v
   Vb = ca.Dot(tmp2)/denom
   X = (a+c+Va*u+Vb*v) * 0.5
   l1 = a - X + u*Va  # l2 = c - X + v*Vb
   dist = 2. * ROOT.TMath.Sqrt( l1.Dot(l1) )
   T = ROOT.TMatrixD(3,12)
   for i in range(3):
     for k in range(4):
       for j in range(3): 
        KD = 0
        if i==j:    #doesnt do this
            KD = 1
        if k==0 or k==2:
       # cova and covc
         temp  = ( u[j]*Vsq - v[j]*UV )*u[i] + (u[j]*UV-v[j]*Usq)*v[i]  #doesnt do this
         sign = -1
         if k==2 :  #doesnt do this
             sign = +1
         T[i][3*k+j] = 0.5*( KD + sign*temp/denom )
        elif k==1:
       # covu
         aNAZ = denom*( ca[j]*Vsq-v.Dot(ca)*v[j] )
         aZAN = ( ca.Dot(u)*Vsq-ca.Dot(v)*UV )*2*( u[j]*Vsq-v[j]*UV )
         bNAZ = denom*( ca[j]*UV+(u.Dot(ca)*v[j]) - 2*ca.Dot(v)*u[j] )
         bZAN = ( ca.Dot(u)*UV-ca.Dot(v)*Usq )*2*( u[j]*Vsq-v[j]*UV )
         T[i][3*k+j] = 0.5*( Va*KD + u[i]/denom**2*(aNAZ-aZAN) + v[i]/denom**2*(bNAZ-bZAN) )   #doesnt do this
        elif k==3:
       # covv
         aNAZ = denom*( 2*ca.Dot(u)*v[j] - ca.Dot(v)*u[j] - ca[j]*UV )
         aZAN = ( ca.Dot(u)*Vsq-ca.Dot(v)*UV )*2*( v[j]*Usq-u[j]*UV )
         bNAZ = denom*( ca.Dot(u)*u[j]-ca[j]*Usq ) 
         bZAN = ( ca.Dot(u)*UV-ca.Dot(v)*Usq )*2*( v[j]*Usq-u[j]*UV )
         T[i][3*k+j] = 0.5*(Vb*KD + u[i]/denom**2*(aNAZ-aZAN) + v[i]/denom**2*(bNAZ-bZAN) )
         LineActivity(get_linenumber()+i+k+j,get_linenumber()) #doesnt do this
   transT = ROOT.TMatrixD(12,3)
   transT.Transpose(T)
   CovTracks = ROOT.TMatrixD(12,12)
   tlist = [t1,t2]
   for k in range(2):
     for i in range(6):
       for j in range(6): 
        xfac = 1.
        if i>2:
            xfac = scalFac[tlist[k]]  
            LineActivity(get_linenumber()+k+i+j,get_linenumber()) #doesnt do this
        if j>2:
            xfac = xfac * scalFac[tlist[k]]
            LineActivity(get_linenumber()+k+i+j,get_linenumber()) #doesnt do this
        CovTracks[i+k*6][j+k*6] = CovMat[tlist[k]][i][j] * xfac
        # if i==5 or j==5 :  CovMat[tlist[k]][i][j] = 0 # ignore error on z-direction
   tmp   = ROOT.TMatrixD(3,12)
   tmp.Mult(T,CovTracks)
   covX  = ROOT.TMatrixD(3,3)
   covX.Mult(tmp,transT)
   return X,covX,dist

def ImpactParameter(point,tPos,tMom):
  t = 0
  if hasattr(tMom,'P'): P = tMom.P()
  else:                 P = tMom.Mag()
  for i in range(3):   t += tMom(i)/P*(point(i)-tPos(i)) 
  dist = 0
  for i in range(3):   dist += (point(i)-tPos(i)-t*tMom(i)/P)**2
  dist = ROOT.TMath.Sqrt(dist)
  return dist

def checkHNLorigin(sTree):
 flag = True
 if not fiducialCut: return flag
# only makes sense for signal == HNL
 if  sTree.MCTrack.GetEntries()<3: return
 # hnlkey = 2 # pythia8 cascade events
 # hnlkey = 1 # pythia8 primary events
 for hnlkey in [1,2]: 
  if abs(sTree.MCTrack[hnlkey].GetPdgCode()) == 9900015:  #would we need to use abs(), also can test this first by outputting to say  print if -9900015
   theHNLVx = sTree.MCTrack[hnlkey+1]  #this must be giving a particle class again, but don't understand name
   X,Y,Z =  theHNLVx.GetStartX(),theHNLVx.GetStartY(),theHNLVx.GetStartZ()  #gives x,y,z of paricle with index hnlkey+1 in MCTrack
   if not isInFiducial(X,Y,Z): flag = False
 return flag 

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

def ecalCluster2MC(aClus):
 # return MC track most contributing, and its fraction of energy
  trackid    = ROOT.Long()
  energy_dep = ROOT.Double()
  mcLink = {}
  for i in range( aClus.Size() ):
    mccell = ecalStructure.GetHitCell(aClus.CellNum(i))  # Get i'th cell of the cluster.
    for n in range( mccell.TrackEnergySize()):
      mccell.GetTrackEnergySlow(n, trackid, energy_dep)
      if not abs(trackid)<sTree.MCTrack.GetEntries(): tid = -1
      else: tid = int(trackid)
      if not mcLink.has_key(tid): mcLink[tid]=0
      mcLink[tid]+=energy_dep
# find trackid most contributing
  eMax,mMax = 0,-1
  for m in mcLink:
     if mcLink[m]>eMax:
        eMax = mcLink[m]
        mMax = m
  return mMax,eMax/aClus.Energy()

def getPtruthFirst(sTree,mcPartKey):
   Ptruth,Ptruthx,Ptruthy,Ptruthz = -1.,-1.,-1.,-1.
   for ahit in sTree.strawtubesPoint:
     if ahit.GetTrackID() == mcPartKey:
        Ptruthx,Ptruthy,Ptruthz = ahit.GetPx(),ahit.GetPy(),ahit.GetPz()
        Ptruth  = ROOT.TMath.Sqrt(Ptruthx**2+Ptruthy**2+Ptruthz**2)
        break
   return Ptruth,Ptruthx,Ptruthy,Ptruthz

def access2SmearedHits():
 key = 0
 for ahit in ev.SmearedHits.GetObject():
   print ahit[0],ahit[1],ahit[2],ahit[3],ahit[4],ahit[5],ahit[6]
   # follow link to true MCHit
   mchit   = TrackingHits[key]
   mctrack =  MCTracks[mchit.GetTrackID()]
   print mchit.GetZ(),mctrack.GetP(),mctrack.GetPdgCode()
   key+=1

def myEventLoop(n):
  global ecalReconstructed
  rc = sTree.GetEntry(n)
# check if tracks are made from real pattern recognition
  measCut = measCutFK
  if sTree.GetBranch("FitTracks_PR"):
      sTree.FitTracks = sTree.FitTracks_PR
      measCut = measCutPR       #doesnt do this
  if sTree.GetBranch("fitTrack2MC_PR"):
      sTree.fitTrack2MC = sTree.fitTrack2MC_PR #doesnt do this
  if sTree.GetBranch("Particles_PR"):
      sTree.Particles   = sTree.Particles_PR  #doesnt do this
  if not checkHNLorigin(sTree):
      return   #does this
  wg = sTree.MCTrack[1].GetWeight()
  if not wg>0.:  #doesnt do this
      wg=1.

# make some ecal cluster analysis if exist
  if sTree.FindBranch("EcalClusters"):
   LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
   if calReco:
       LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
       ecalReconstructed.Delete()
   else:
       LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
       ecalReconstructed = sTree.EcalReconstructed
   for x in caloTasks:
    if x.GetName() == 'ecalFiller': 
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        x.Exec('start',sTree.EcalPointLite)
    elif x.GetName() == 'ecalMatch':
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        x.Exec('start',ecalReconstructed,sTree.MCTrack)
    else : 
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        x.Exec('start')
   for aClus in ecalReconstructed:
    mMax = aClus.MCTrack()
    if mMax <0 or mMax > sTree.MCTrack.GetEntries():
        aP = None # this should never happen, otherwise the ECAL MC matching has a bug
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
    else: 
        aP = sTree.MCTrack[mMax]
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
    if aP:    
      tmp = PDG.GetParticle(aP.GetPdgCode())
      LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
      if tmp:
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          pName = 'ecalReconstructed_'+tmp.GetName()
      else:
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          pName = 'ecalReconstructed_'+str(aP.GetPdgCode())
    else:
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        pName = 'ecalReconstructed_unknown' 
    if not h.has_key(pName): 
      ut.bookHist(h,pName,'x/y and energy for '+pName.split('_')[1],50,-3.,3.,50,-6.,6.)
      LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
    rc = h[pName].Fill(aClus.X()/u.m,aClus.Y()/u.m,aClus.RecoE()/u.GeV)
# look at distance to tracks 
    for fT in sTree.FitTracks:
     rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(fT,z_ecal)
     if rc:
      pdgcode = fT.getFittedState().getPDG()
      tmp = PDG.GetParticle(pdgcode)
      if tmp:
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          tName = 'ecalReconstructed_dist_'+tmp.GetName()
      else:
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          tName = 'ecalReconstructed_dist_'+str(aP.GetPdgCode())
      if not h.has_key(tName):
       LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
       p = tName.split('dist_')[1]
       ut.bookHist(h,tName,'Ecal cluster distance t0 '+p,100,0.,100.*u.cm)
       ut.bookHist(h,tName.replace('dist','distx'),'Ecal cluster distance to '+p+' in X ',100,-50.*u.cm,50.*u.cm)
       ut.bookHist(h,tName.replace('dist','disty'),'Ecal cluster distance to '+p+' in Y ',100,-50.*u.cm,50.*u.cm)
      dist = ROOT.TMath.Sqrt( (aClus.X()-pos.X())**2+(aClus.Y()-pos.Y())**2 )
      rc = h[tName].Fill(dist)
      rc = h[tName.replace('dist','distx')].Fill( aClus.X()-pos.X() )
      rc = h[tName.replace('dist','disty')].Fill( aClus.Y()-pos.Y() )
# compare with old method
   for aClus in sTree.EcalClusters:
     LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
     rc = h['ecalClusters'].Fill(aClus.X()/u.m,aClus.Y()/u.m,aClus.Energy()/u.GeV)
     mMax,frac = ecalCluster2MC(aClus)
# return MC track most contributing, and its fraction of energy
     if mMax>0:    
      aP = sTree.MCTrack[mMax]   
      tmp = PDG.GetParticle(aP.GetPdgCode())
      LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
      if tmp: 
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          pName = 'ecalClusters_'+tmp.GetName()
      else:
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          pName = 'ecalClusters_'+str(aP.GetPdgCode())
     else:
         LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
         pName = 'ecalClusters_unknown' 
     if not h.has_key(pName):
         LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
         ut.bookHist(h,pName,'x/y and energy for '+pName.split('_')[1],50,-3.,3.,50,-6.,6.)
     rc = h[pName].Fill(aClus.X()/u.m,aClus.Y()/u.m,aClus.Energy()/u.GeV)
     
# make some straw hit analysis
  hitlist = {}
  for ahit in sTree.strawtubesPoint:
     detID = ahit.GetDetectorID()
     top = ROOT.TVector3()
     bot = ROOT.TVector3()
     modules["Strawtubes"].StrawEndPoints(detID,bot,top)
     dw  = ahit.dist2Wire()
     if detID < 50000000 :
      LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
      if abs(top.y())==abs(bot.y()):
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          h['disty'].Fill(dw)
      if abs(top.y())>abs(bot.y()):
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          h['distu'].Fill(dw)
      if abs(top.y())<abs(bot.y()):
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          h['distv'].Fill(dw)
#
     trID = ahit.GetTrackID()
     if not trID < 0 :
         if hitlist.has_key(trID):
             LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
             hitlist[trID]+=1
         else:
             LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
             hitlist[trID]=1
  for tr in hitlist:
      LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
      h['meanhits'].Fill(hitlist[tr])
  key = -1
  fittedTracks = {}
  for atrack in sTree.FitTracks:
   key+=1
# kill tracks outside fiducial volume
   if not checkFiducialVolume(sTree,key,dy):
       LineActivity(get_linenumber()+n+key, get_linenumber()) #doesnt do this
       continue
   fitStatus   = atrack.getFitStatus()
   nmeas = fitStatus.getNdf()
   h['meas'].Fill(nmeas)
   if not fitStatus.isFitConverged() :
       LineActivity(get_linenumber()+n+key, get_linenumber()) #doesnt do this
       continue
   h['meas2'].Fill(nmeas)
   if nmeas < measCut:
       LineActivity(get_linenumber()+n+key, get_linenumber()) #doesnt do this
       continue
   fittedTracks[key] = atrack
# needs different study why fit has not converged, continue with fitted tracks
   rchi2 = fitStatus.getChi2()
   prob = ROOT.TMath.Prob(rchi2,int(nmeas))
   h['prob'].Fill(prob)
   chi2 = rchi2/nmeas
   fittedState = atrack.getFittedState()
   h['chi2'].Fill(chi2,wg)
   h['measVSchi2'].Fill(atrack.getNumPoints(),chi2)
   P = fittedState.getMomMag()
   Px,Py,Pz = fittedState.getMom().x(),fittedState.getMom().y(),fittedState.getMom().z()
   cov = fittedState.get6DCov()
   if len(sTree.fitTrack2MC)-1<key:
       LineActivity(get_linenumber()+n+key, get_linenumber()) #doesnt do this
       continue
   mcPartKey = sTree.fitTrack2MC[key]
   mcPart    = sTree.MCTrack[mcPartKey]
   if not mcPart :
       LineActivity(get_linenumber()+n+key, get_linenumber()) #doesnt do this
       continue
   Ptruth_start     = mcPart.GetP()
   Ptruthz_start    = mcPart.GetPz()
   # get p truth from first strawpoint
   Ptruth,Ptruthx,Ptruthy,Ptruthz = getPtruthFirst(sTree,mcPartKey)
   delPOverP = (Ptruth - P)/Ptruth
   h['delPOverP'].Fill(Ptruth,delPOverP)
   delPOverPz = (1./Ptruthz - 1./Pz) * Ptruthz
   h['pullPOverPx'].Fill( Ptruth,(Ptruthx-Px)/ROOT.TMath.Sqrt(cov[3][3]) )   
   h['pullPOverPy'].Fill( Ptruth,(Ptruthy-Py)/ROOT.TMath.Sqrt(cov[4][4]) )   
   h['pullPOverPz'].Fill( Ptruth,(Ptruthz-Pz)/ROOT.TMath.Sqrt(cov[5][5]) )   
   h['delPOverPz'].Fill(Ptruthz,delPOverPz)
   if chi2>chi2CutOff:
       LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
       continue
   h['delPOverP2'].Fill(Ptruth,delPOverP)
   h['delPOverP2z'].Fill(Ptruth,delPOverPz)
# try measure impact parameter
   trackDir = fittedState.getDir()
   trackPos = fittedState.getPos()
   vx = ROOT.TVector3()
   mcPart.GetStartVertex(vx)
   t = 0
   for i in range(3):
       LineActivity(get_linenumber()+n+key+i, get_linenumber()) #doesnt do this
       t += trackDir(i)*(vx(i)-trackPos(i)) 
   dist = 0
   for i in range(3):
       LineActivity(get_linenumber()+n+key+i, get_linenumber()) #doesnt do this
       dist += (vx(i)-trackPos(i)-t*trackDir(i))**2
   dist = ROOT.TMath.Sqrt(dist)
   h['IP'].Fill(dist) 
# ---
# loop over particles, 2-track combinations
  for HNL in sTree.Particles:
    t1,t2 = HNL.GetDaughter(0),HNL.GetDaughter(1) 
# kill tracks outside fiducial volume, if enabled
    if not checkFiducialVolume(sTree,t1,dy) or not checkFiducialVolume(sTree,t2,dy) :
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        continue
    checkMeasurements = True
# cut on nDOF
    for tr in [t1,t2]:
      fitStatus  = sTree.FitTracks[tr].getFitStatus()
      nmeas = fitStatus.getNdf()
      if nmeas < measCut:
          LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
          checkMeasurements = False
    if not checkMeasurements:
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        continue
# check mc matching 
    if not match2HNL(HNL):
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        continue
    HNLPos = ROOT.TLorentzVector()
    HNL.ProductionVertex(HNLPos)
    HNLMom = ROOT.TLorentzVector()
    HNL.Momentum(HNLMom)
# check if DOCA info exist
    if hasattr(HNL,"GetDoca"):
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        doca  =  HNL.GetDoca()
    elif HNL.GetMother(1)==99 :
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        doca  =  HNLPos.T()
    else:
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        # redo doca calculation
        xv,yv,zv,doca,HNLMom  = RedoVertexing(t1,t2)
        if HNLMom == -1:
            LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
            continue
 # check if decay inside decay volume
    if not isInFiducial(HNLPos.X(),HNLPos.Y(),HNLPos.Z()):
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        continue  
    h['Doca'].Fill(doca) 
    if  doca > docaCut :
        LineActivity(get_linenumber()+n, get_linenumber()) #doesnt do this
        continue
    tr = ROOT.TVector3(0,0,ShipGeo.target.z0)
    dist = ImpactParameter(tr,HNLPos,HNLMom)
    #maybe add if statement here?
    #if (HNL.GetPdgCode() == 9900015) or ((HNL.GetPdgCode() == 13) and (HNL.GetMotherId() == 9900015)) or ((HNL.GetPdgCode() == 211) and (HNL.GetMotherId() == 9900015)):
    mass = HNLMom.M()

    h['IP0'].Fill(dist)  
    h['IP0/mass'].Fill(mass,dist)
    h['HNL'].Fill(mass)
    h['HNLw'].Fill(mass,wg)
#
    vetoDets['SBT'] = veto.SBT_decision()
    vetoDets['SVT'] = veto.SVT_decision()
    vetoDets['UVT'] = veto.UVT_decision()
    vetoDets['RPC'] = veto.RPC_decision()
    vetoDets['TRA'] = veto.Track_decision()
    h['nrtracks'].Fill(vetoDets['TRA'][2])
    h['nrSVT'].Fill(vetoDets['SVT'][2])
    h['nrUVT'].Fill(vetoDets['UVT'][2])
    h['nrSBT'].Fill(vetoDets['SBT'][2])
    h['nrRPC'].Fill(vetoDets['RPC'][2])
#   HNL true
    mctrack = sTree.MCTrack[sTree.fitTrack2MC[t1]]
    h['Vzresol'].Fill( (mctrack.GetStartZ()-HNLPos.Z())/u.cm )
    h['Vxresol'].Fill( (mctrack.GetStartX()-HNLPos.X())/u.cm )
    h['Vyresol'].Fill( (mctrack.GetStartY()-HNLPos.Y())/u.cm )
    PosDir,newPosDir,CovMat,scalFac = {},{},{},{}
# opening angle at vertex
    newPos = ROOT.TVector3(HNLPos.X(),HNLPos.Y(),HNLPos.Z())
    st1,st2 = sTree.FitTracks[t1].getFittedState(),sTree.FitTracks[t2].getFittedState()
    PosDir[t1] = {'position':st1.getPos(),'direction':st1.getDir(),'momentum':st1.getMom()}
    PosDir[t2] = {'position':st2.getPos(),'direction':st2.getDir(),'momentum':st2.getMom()}
    CovMat[t1] = st1.get6DCov() 
    CovMat[t2] = st2.get6DCov() 
    rep1,rep2 = ROOT.genfit.RKTrackRep(st1.getPDG()),ROOT.genfit.RKTrackRep(st2.getPDG())  
    state1,state2 = ROOT.genfit.StateOnPlane(rep1),ROOT.genfit.StateOnPlane(rep2)
    rep1.setPosMom(state1,st1.getPos(),st1.getMom())
    rep2.setPosMom(state2,st2.getPos(),st2.getMom())
    try:
     rep1.extrapolateToPoint(state1, newPos, False)
     rep2.extrapolateToPoint(state2, newPos, False)
     mom1,mom2 = rep1.getMom(state1),rep2.getMom(state2)
    except:
     mom1,mom2 = st1.getMom(),st2.getMom()
    newPosDir[t1] = {'position':rep1.getPos(state1),'direction':rep1.getDir(state1),'momentum':mom1}
    newPosDir[t2] = {'position':rep2.getPos(state2),'direction':rep2.getDir(state2),'momentum':mom2}
    oa = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag()) 
    h['oa'].Fill(oa)
#
    covX = HNL.GetCovV()
    dist = HNL.GetDoca()
    h['Vzpull'].Fill( (mctrack.GetStartZ()-HNLPos.Z())/ROOT.TMath.Sqrt(covX[2][2]) )
    h['Vxpull'].Fill( (mctrack.GetStartX()-HNLPos.X())/ROOT.TMath.Sqrt(covX[0][0]) )
    h['Vypull'].Fill( (mctrack.GetStartY()-HNLPos.Y())/ROOT.TMath.Sqrt(covX[1][1]) )
h = {}
def create_Hists():
    #DECLARE HISTOGRAMS
    ut.bookHist(h,'delPOverP','delP / P',400,0.,200.,100,-0.5,0.5)
    ut.bookHist(h,'pullPOverPx','delPx / sigma',400,0.,200.,100,-3.,3.)
    ut.bookHist(h,'pullPOverPy','delPy / sigma',400,0.,200.,100,-3.,3.)
    ut.bookHist(h,'pullPOverPz','delPz / sigma',400,0.,200.,100,-3.,3.)
    ut.bookHist(h,'delPOverP2','delP / P chi2/nmeas<'+str(chi2CutOff),400,0.,200.,100,-0.5,0.5)
    ut.bookHist(h,'delPOverPz','delPz / Pz',400,0.,200.,100,-0.5,0.5)
    ut.bookHist(h,'delPOverP2z','delPz / Pz chi2/nmeas<'+str(chi2CutOff),400,0.,200.,100,-0.5,0.5)
    ut.bookHist(h,'chi2','chi2/nmeas after trackfit',100,0.,10.)##
    ut.bookHist(h,'prob','prob(chi2)',100,0.,1.)##
    ut.bookHist(h,'IP','Impact Parameter',100,0.,10.)
    ut.bookHist(h,'Vzresol','Vz reco - true [cm]',100,-50.,50.)
    ut.bookHist(h,'Vxresol','Vx reco - true [cm]',100,-10.,10.)
    ut.bookHist(h,'Vyresol','Vy reco - true [cm]',100,-10.,10.)
    ut.bookHist(h,'Vzpull','Vz pull',100,-5.,5.)
    ut.bookHist(h,'Vxpull','Vx pull',100,-5.,5.)
    ut.bookHist(h,'Vypull','Vy pull',100,-5.,5.)
    ut.bookHist(h,'Doca','Doca between two tracks',100,0.,10.)
    ut.bookHist(h,'IP0','Impact Parameter to target',100,0.,100.)
    ut.bookHist(h,'IP0/mass','Impact Parameter to target vs mass',100,0.,2.,100,0.,100.)
    #
    ut.bookHist(h,'HNL','Original Reconstructed Mass',500,0.,2.) # original one for the reconstructed mass
    ut.bookHist(h,'HNL_true','True Mass',500,0.,2.) # new one for the simulated (Monte Carlo) mass
    ut.bookHist(h,'HNL_reco','Reconstructed Mass',500,0.,2.) # new one for the reconstructed mass
    ut.bookHist(h,'Chi2','Fitted Tracks Chi Squared',100,0.,3.)
    ut.bookHist(h,'HNL_mom','Monte Carlo Momentum',100,0.,300.)
    ut.bookHist(h,'Time','Muon - Time Between Straw Tube and ECAL',500,36.,42.)
    ut.bookHist(h,'Time2','Pion - Time Between Straw Tube and ECAL',500,36.,42.)
    ut.bookHist(h,'HNL_mom_reco','Reconstructed Momentum',100,0.,300)
    ut.bookHist(h,'Time_man','Analytical Time',500,36.,42.)
    #
    ut.bookHist(h,'HNLw','Reconstructed Mass with weights',500,0.,2.)
    ut.bookHist(h,'meas','number of measurements',40,-0.5,39.5)##
    ut.bookHist(h,'meas2','number of measurements, fitted track',40,-0.5,39.5)##
    ut.bookHist(h,'measVSchi2','number of measurements vs chi2/meas',40,-0.5,39.5,100,0.,10.)##
    ut.bookHist(h,'distu','distance to wire',100,0.,1.)
    ut.bookHist(h,'distv','distance to wire',100,0.,1.)
    ut.bookHist(h,'disty','distance to wire',100,0.,1.)
    ut.bookHist(h,'meanhits','mean number of hits / track',50,-0.5,49.5)
    ut.bookHist(h,'ecalClusters','x/y and energy',50,-3.,3.,50,-6.,6.)
    ut.bookHist(h,'oa','cos opening angle',100,0.999,1.)
    ut.bookHist(h,'nrtracks','nr of tracks in signal selected',10,-0.5,9.5)
    ut.bookHist(h,'nrSVT','nr of hits in SVT',10,-0.5,9.5)
    ut.bookHist(h,'nrUVT','nr of hits in UVT',100,-0.5,99.5)
    ut.bookHist(h,'nrSBT','nr of hits in SBT',100,-0.5,99.5)
    ut.bookHist(h,'nrRPC','nr of hits in RPC',100,-0.5,99.5)
create_Hists()

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

   #cv = h['Test_Mass'].cd(3)
   #h['Chi2'].Draw()
   cv = h['Test_Mass'].cd(3)
   h['HNL_mom'].SetXTitle('Momentum [GeV/c]')
   h['HNL_mom'].SetYTitle('No. of Particles')
   h['HNL_mom'].Draw()

   cv = h['Test_Mass'].cd(4)
   h['HNL_mom_reco'].SetXTitle('Momentum [GeV/c]')
   h['HNL_mom_reco'].SetYTitle('No. of Particles')
   h['HNL_mom_reco'].Draw()
   h['Test_Mass'].Print('MassAndMom.png')
   #----------------------------------------------------------------------------------------------------------------------
   ut.bookCanvas(h,key='Time_Res',title='Time Plots',nx=900,ny=700,cx=1,cy=1)
   #cv = h['Time_Res'].cd(1)
   h['Time'].SetXTitle('Time [ns]')
   h['Time'].SetYTitle('Frequency')
   h['Time'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   #cv = h['Time_Res'].cd(2)
   h['Time2']. SetLineColor(2)
   h['Time2'].SetXTitle('Time [ns]')
   h['Time2'].SetYTitle('Frequency')
   h['Time2'].Draw('same')
   #h['Time2'].Draw()
   #----------------------------------------------------------------------------------------------------------------------
   #cv = h['Time_Res'].cd(3)
   h['Time_man'].SetLineColor(1)
   h['Time_man'].SetXTitle('Time [ns]')
   h['Time_man'].SetYTitle('Frequency')
   h['Time_man'].Draw('same')
   #h['Time_man'].Draw()
   h['Time_Res'].Print('TimeRes.png')
   print ('Made the plots')

for n in (nEvents):
    rc = sTree.GetEntry(n)
    if sTree.GetBranch("strawtubesPoint"):
        z_array = []
        t_array = []
        straw_time=0
        for k,hits in enumerate(sTree.strawtubesPoint):
            TrackID = hits.GetTrackID()
            if TrackID == partkey:
                z_array.append(hits.GetZ())
                t_array.append(hits.GetTime())
        straw_zpos=min(z_array)   
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
                            return ecal_time,ecal_zpos, straw_zpos, straw_time, t 
                    else: return None,None,None,None,None
            else: return None,None,None,None,None
    else: return None,None,None,None,None
HNLKinematics()  

def time_res(partkey):
    if sTree.GetBranch("strawtubesPoint"):
        z_array = []
        t_array = []
        straw_time=0
        for k,hits in enumerate(sTree.strawtubesPoint):
            TrackID = hits.GetTrackID()
            #print(TrackID)
            if TrackID == partkey:
                LineActivity(get_linenumber(),get_linenumber())
                z_array.append(hits.GetZ())
                t_array.append(hits.GetTime())
            else: LineActivity(get_linenumber(),get_linenumber())
        straw_zpos=min(z_array)   
        min_z = z_array.index(min(z_array))
        straw_time = t_array[min_z]
    else:
        LineActivity(get_linenumber(),get_linenumber())
        return None,None,None,None,None
    if sTree.GetBranch("EcalPoint"):
            if not straw_time<=0:
                for k,hits in enumerate(sTree.EcalPoint):
                    TrackID = hits.GetTrackID()
                    if TrackID == partkey:
                        LineActivity(get_linenumber(),get_linenumber())
                        ecal_time = hits.GetTime()
                        ecal_zpos = hits.GetZ()
                        if not ecal_time <= straw_time:
                            t = abs(straw_time - ecal_time)
                            return ecal_time,ecal_zpos, straw_zpos, straw_time, t
                        else:LineActivity(get_linenumber(),get_linenumber())
                    else:
                        LineActivity(get_linenumber(),get_linenumber())
                        return None,None,None,None,None
            else:
                LineActivity(get_linenumber(),get_linenumber())
                return None,None,None,None,None
    else:
        LineActivity(get_linenumber(),get_linenumber())
        return None,None,None,None,None

def time_resVrs2(partkey):
    if sTree.GetBranch("strawtubesPoint"):
        x_array = []
        y_array = []
        z_array = []
        t_array = []
        straw_time=0
        for k,hits in enumerate(sTree.strawtubesPoint):
            #print(k,hits)
            TrackID = hits.GetTrackID()
            #print(TrackID)
            if TrackID == partkey:
                LineActivity(get_linenumber(),get_linenumber())
                x_array.append(hits.GetX())
                y_array.append(hits.GetY())
                z_array.append(hits.GetZ())
                t_array.append(hits.GetTime())
            else: LineActivity(get_linenumber(),get_linenumber())
        min_z = z_array.index(min(z_array))
        straw_zpos=min(z_array)
        straw_xpos = x_array[min_z]
        straw_ypos = y_array[min_z]
        straw_time = t_array[min_z]
    else:
        LineActivity(get_linenumber(),get_linenumber())
        return None,None,None,None,None,None,None,None,None
    if sTree.GetBranch("EcalPoint"):
            if not straw_time<=0:
                for k,hits in enumerate(sTree.EcalPoint):
                    #print(k,hits)
                    TrackID = hits.GetTrackID()
                    #print(TrackID)
                    ecal_time = hits.GetTime()
                    ecal_zpos = hits.GetZ()
                    ecal_xpos = hits.GetX()
                    ecal_ypos = hits.GetY()
                    if not ecal_time <= straw_time:
                        t = abs(straw_time - ecal_time)
                        return ecal_time, ecal_xpos, ecal_ypos, ecal_zpos, straw_xpos,straw_ypos,straw_zpos, straw_time, t 
            else:
                LineActivity(get_linenumber(),get_linenumber())
                return None,None,None,None,None,None,None,None,None
    else:
        LineActivity(get_linenumber(),get_linenumber())
        return None,None,None,None,None,None,None,None,None

def finStateMuPi():
    #INVARIANT MASS OF HNL THAT DECAYED TO MU-PI
    if sTree.GetBranch("FitTracks"):
        pi_decaycheck = 0
        for n in range(nEvents):                            #loop over events
            rc = sTree.GetEntry(n)                              #load tree entry
            keylist=[]                                          #create empty list
            muVector = {}                                       #create empty dictionaries
            dicMuChi2 = {}
            mupartkey={}
            for index,reco_part in enumerate(sTree.FitTracks):  #loops over index and data of track particles
                if not checkFiducialVolume(sTree,index,dy): break   #exits loop if HNL decayed in fiducial volume 
                muPartkey = sTree.fitTrack2MC[index]                #mathches track to MC particle key
                true_muon = sTree.MCTrack[muPartkey]                #gives MC particle data
                if abs(true_muon.GetPdgCode()) == 13:               #checks particle is muon
                    muonMotherkey = true_muon.GetMotherId()             #stores a number index of MC track of mother
                    keylist.append(muonMotherkey)                       #adds to list
                    if true_mother.GetPdgCode() == 211:                 #checks mother is HNL
                        print('Pion has decayed to a muon')
                        pi_decaycheck+=1
                    true_mother = sTree.MCTrack[muonMotherkey]          #obtains mother particle data
                    if true_mother.GetPdgCode() == 9900015:             #checks mother is HNL
                        Decay_X = true_muon.GetStartX()
                        Decay_Y = true_muon.GetStartY()
                        Decay_Z = true_muon.GetStartZ()
                        if not isInFiducial(Decay_X,Decay_Y,Decay_Z):
                            #print('HNL decayed outside fiducial volume')
                            continue
                        if not checkFiducialVolume(sTree,index,dy): 
                            #print('Track outside fiducial volume')
                            continue 
                        mu_status = reco_part.getFitStatus()                #gets fit status           
                        if not mu_status.isFitConverged():
                            #print('Fit did not converge')
                            continue
                        mu_nmeas = mu_status.getNdf()                       #gets number of measurements                      
                        if not mu_nmeas > 25:
                            #print('Too few measurements')
                            continue

                        mu_rchi2 = mu_status.getChi2()                      #gets chi squared value
                        mu_chi2 = (mu_rchi2/mu_nmeas)                       #gets chi value
                        if not mu_chi2 < 4:
                            #print('Chi squared value too high')
                            dicMuChi2[str(muonMotherkey)]=mu_chi2
                            continue

                        mupartkey[str(muonMotherkey)] = muPartkey           #muon Track Id

                        fittedstate1 = reco_part.getFittedState()           #get reconstructed muon fitted state
                        mu_M = true_muon.GetMass()                          #gets mass of MC muon
                        muPx = fittedstate1.getMom().x()                    #then its momentum in x,
                        muPy = fittedstate1.getMom().y()                    #y  
                        muPz = fittedstate1.getMom().z()                    #and z
                        muP = fittedstate1.getMomMag()                      #then its momentum magnitude
                        muE = ROOT.TMath.Sqrt((mu_M**2) + (muP**2))         #then its energy

                        Muon_Vector = ROOT.TLorentzVector()                 # declares variable as TLorentzVector class
                        Muon_Vector.SetPxPyPzE(muPx,muPy,muPz,muE)          # inputs four-momentum elements
                        
                        muVector[str(muonMotherkey)]=Muon_Vector
                        muonMotherTrue_mass = true_mother.GetMass()         #gets HNL mass

            for index,reco_part in enumerate(sTree.FitTracks):  #loops over index and data of track particles
                piPartkey = sTree.fitTrack2MC[index]                  #mathches track to MC particle key
                true_pion = sTree.MCTrack[piPartkey]                  #gives MC particle data
                if abs(true_pion.GetPdgCode()) == 211:              #checks particle is pion
                    pionMotherkey = true_pion.GetMotherId()             #stores a number index of MC track of mother
                    true_mother = sTree.MCTrack[pionMotherkey]          #obtains mother particle data
                    for muonMotherkey in keylist:                       #loops through muonMother keys
                        if pionMotherkey==muonMotherkey:                    #check if keys are the same
                            if true_mother.GetPdgCode() == 9900015:             #check if mother is HNL
                                pionMotherTrue_mass = true_mother.GetMass()         #get HNL/final states mother mass
                                pionMotherTrue_mom = true_mother.GetP()             #get HNL/final states mother mom
                                h['HNL_true'].Fill(pionMotherTrue_mass)             #true HNL mass
                                h['HNL_mom'].Fill(pionMotherTrue_mom)               # true HNL momentum

                                if not checkFiducialVolume(sTree,index,dy): 
                                        #print('Decay outside fiducial volume')
                                        continue 
                                pi_status = reco_part.getFitStatus()                #gets fit status                
                                if not pi_status.isFitConverged():
                                    #print('Fit did not converge')
                                    continue
                                pi_nmeas = pi_status.getNdf()                       #gets number of measurements 
                                if not pi_nmeas > 25:
                                    #print('Too few measurements')
                                    continue

                                pi_rchi2 = pi_status.getChi2()                      #chi squared value
                                pi_chi2 = (pi_rchi2/pi_nmeas)                       #gets chi value
                                if not pi_chi2 < 4:
                                        #print('Chi squared value too high')
                                        continue

                                fittedstate2 = reco_part.getFittedState()           #get reconstructed pion fitted state
                                pi_M = true_pion.GetMass()                          #gets mass of MC muon
                                piP = fittedstate2.getMomMag()                      #then its momentum in x,
                                piPx = fittedstate2.getMom().x()                    #y
                                piPy = fittedstate2.getMom().y()                    #and z
                                piPz = fittedstate2.getMom().z()                    #then its momentum magnitude
                                #print(piPz)                                         #prints TrackID of pion
                                piE = ROOT.TMath.Sqrt((pi_M**2) + (piP**2))         #then its energy

                                Pion_Vector = ROOT.TLorentzVector()                 #declares variable as TLorentzVector class
                                Pion_Vector.SetPxPyPzE(piPx,piPy,piPz,piE)          #inputs four-momentum elements
                            
                                HNL_Vector = muVector[str(muonMotherkey)] + Pion_Vector #adds the 4-momenta
                                HNL_mass = HNL_Vector.M()                           #sets HNL mass
                                HNL_reco_mom = HNL_Vector.P()                       #sets HNL mom
                                h['HNL_reco'].Fill(HNL_mass)                        #fill histograms
                                h['HNL_mom_reco'].Fill(HNL_reco_mom)                #----||-------
                                h['Chi2'].Fill(dicMuChi2[str(muonMotherkey)])       #----||-------
                                h['Chi2'].Fill(pi_chi2)                             #----||-------
                                #print(mupartkey[str(muonMotherkey)])                #prints Track ID of muon
                                
                                muEcalT,muEcalZ, muMinStrawZ, muStrawT, mu_t = time_res(mupartkey[str(muonMotherkey)])      #
                                #particleDataFile.write('mu: \t' + str(muEcalT) + '\t' + str(muEcalZ) + '\t' + str(muMinStrawZ) + '\t' + str(muStrawT) + '\t' + str(mu_t) + '\n')
                                if mu_t != None:                                    #
                                    h['Time'].Fill(mu_t)                            #
                                speedOfLight=3*(10**8)
                                piEcalT, piEcalX, piEcalY, piEcalZ, pistraw_xpos,pistraw_ypos,pistraw_zpos, piStrawT, pi_t  = time_resVrs2(piPartkey)                          #
                                piDeltaPos=(ROOT.TMath.Sqrt(((piEcalX-pistraw_xpos)**2)+((piEcalY-pistraw_ypos)**2)+((piEcalZ-pistraw_zpos)**2)))/100
                                piPjoules=piP*1.602*(10**(-13))
                                piGammaVel=piP/(pi_M/(3*(10**8)))
                                Betta=piP/ROOT.TMath.Sqrt((pi_M**2)+(piP**2))
                                manualCalTime=(piDeltaPos*(10**9))/(Betta*speedOfLight)
                                #pideltaZz=piEcalZ-piMinStrawZ
                                #particleDataFile.write('pi: \t' + str(piEcalT) + '\t' + str(piEcalZ) + '\t' + str(piMinStrawZ) + '\t' + str(piStrawT) + '\t' + str(pi_t) + '\t' +  str(piPz) + '\n')
                                particleDataFile.write('pi: \t'+ str(piDeltaPos)+ '\t'  +  str(piP) + '\t'  +  str(piPjoules) +'\t' +  str(manualCalTime) + '\n')
                                h['Time_man'].Fill(manualCalTime)
                                if pi_t != None:                                     #
                                    h['Time2'].Fill(pi_t) 
                                    
#for n in range(nEvents):
#    myEventLoop(n)
#    sTree.FitTracks.Delete()