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