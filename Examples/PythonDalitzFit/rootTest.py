import ROOT
from ROOT import gROOT, TCanvas, TF1, TVector3, TTree

import Dalitz_ext
from Dalitz_ext import *

import numpy as np

fit = PythonFit()

print "TVector 3 from ComPWA"
vecIn = TVector3(1,2,3)
mult = vecIn.X()*vecIn.Y()*vecIn.Z()
print mult
vec = fit.testVec(vecIn)
mult = vec.X()*vec.Y()*vec.Z()
print mult

print "TTree from ComPWA"
tt = fit.testTree();

print "Writing a tree"

f = ROOT.TFile("tree.root", "recreate")
t = ROOT.TTree("name_of_tree", "tree title")
#newtree = tt.CloneTree(0)
#cobj_tree = ROOT.AsCObject( tree )
#myroot.printTree( newtree )

# create 1 dimensional float arrays (python's float datatype corresponds to c++ doubles)
# as fill variables
n = np.zeros(1, dtype=float)
u = np.zeros(1, dtype=float)

# create the branches and assign the fill-variables to them
t.Branch('normal', n, 'normal/D')
t.Branch('uniform', u, 'uniform/D')

# create some random numbers, fill them into the fill varibles and call Fill()
for i in xrange(100000):
	n[0] = ROOT.gRandom.Gaus()
	u[0] = ROOT.gRandom.Uniform()
	t.Fill()

# write the tree into the output file and close the file
tt.Write()
f.Write()
f.Close()

print "test ComPWA fit"

#fit.StartFit()
exit()
 