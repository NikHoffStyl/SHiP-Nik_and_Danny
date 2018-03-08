import ROOT, os, sys
import shipunit as u
import rpvsusy
from pythia8_conf_utils import *
import readDecayTable

rpvsusy_instance = rpvsusy.RPVSUSY(1.,[1,1],1e3,1,True)
bRatio = rpvsusy_instance.findBrancingRatio('N -> K+ mu-')
print('Shit Happens Ratio is:')
#print(bRatio)