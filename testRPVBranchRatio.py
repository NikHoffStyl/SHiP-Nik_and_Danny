#import ROOT, os, sys
#import shipunit as u
import rpvsusy
#from pythia8_conf_utils import *
#import readDecayTable

rpvsusy_instance = rpvsusy.RPVSUSY(1.,[0.2,0.03],1e3,1,True)
bRatio = rpvsusy_instance.findDecayBranchingRatio('N -> K+ mu-')
print('Ratio is:')
print(bRatio)