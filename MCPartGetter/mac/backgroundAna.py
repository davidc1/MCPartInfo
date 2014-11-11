import os, ROOT, sys
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.path import Path
from matplotlib.lines import Line2D

def col(x):
    c = ""
    if ( (x == 13) or (x == -13) ):
        c = "r-"
    if (x == 22):
        c = "k-"
    if ( (x == 11) or (x == -11) ):
        c = "b-"

    return c


# Now import ana_processor & your class. For this example, ana_base.
gSystem.Load("libLArLite_Analysis")
gSystem.Load("libLArLite_LArUtil")

geom = larutil.Geometry.GetME()
detLength = geom.DetLength()
detHalfWidth  = geom.DetHalfWidth()
detHalfHeight = geom.DetHalfHeight()


ProfYX = plt.Rectangle( (0,-detHalfHeight), 2*detHalfWidth, 2*detHalfHeight, facecolor="#aaaaaa")
ProfYZ = plt.Rectangle( (0,-detHalfHeight), detLength, 2*detHalfHeight, facecolor="#aaaaaa")
ProfXZ = plt.Rectangle( (0,0), detLength, 2*detHalfWidth, facecolor="#aaaaaa")

plt.gca().add_patch(ProfYZ)

#legend:
lineBlue = Line2D([],[],color='blue',linewidth=5)
lineRed = Line2D([],[],color='red',linewidth=5)
lineBlack = Line2D([],[],color='black',linewidth=5)
plt.legend([ProfYZ,lineBlue,lineRed,lineBlack], ["TPC Boundaries","Electrons","Muons","Photons"],loc=4,fontsize=20)

f = ROOT.TFile("tot.root")
t = f.Get("ana_tree")
m = f.Get("muon_tree")


muonEntries = m.GetEntries()
for j in range(muonEntries):

numEntries = t.GetEntries()
print numEntries
for i in range(numEntries):
    t.GetEntry(i)
    #      #if in TPC   # and either parent is charged and too far away or parent is not charged
    stri = str(t.ProcHist)
    #print stri[-8:]
    if ((t.inTPC == 1) and (t.Energy>0.01)):
        print "Found!"
#    if ((t.inTPC == 1) and (t.isCut==0) and (stri[-8:] == "compt:11")):#and ( (t.MotherDist >= 5. and t.MotherPDG != 22) or (t.MotherPDG == 22) ) ): 
        #((t.ProcHist == "primary:13_muIoni:11_eBrem:22_compt:11") or (t.ProcHist == "primary:-13_muIoni:11_eBrem:22_compt:11")) ):
        MotherSteps = t.MotherTraj.size()
        PartSteps   = t.PartTraj.size()
        AncestorSteps = t.AncestorTraj.size()
        XpointsM = np.array([])
        YpointsM = np.array([])
        ZpointsM = np.array([])
        XpointsP = np.array([])
        YpointsP = np.array([])
        ZpointsP = np.array([])
        XpointsA = np.array([])
        YpointsA = np.array([])
        ZpointsA = np.array([])
        
        XpointsP = np.append( XpointsP, t.StartX )
        YpointsP = np.append( XpointsP, t.StartY )
        ZpointsP = np.append( XpointsP, t.StartZ )
        plt.plot(ZpointsP,YpointsP,"bo")
        '''
        if (AncestorSteps > 0):
            for n in range(AncestorSteps):
                XpointsA = np.append( XpointsA, t.AncestorTraj.at(n).at(0) )
                YpointsA = np.append( YpointsA, t.AncestorTraj.at(n).at(1) )
                ZpointsA = np.append( ZpointsA, t.AncestorTraj.at(n).at(2) )
            plt.plot(ZpointsA,YpointsA,col(t.MotherPDG),linewidth=2)
        if (MotherSteps > 0):
            for n in range(MotherSteps):
                XpointsM = np.append( XpointsM, t.MotherTraj.at(n).at(0) )
                YpointsM = np.append( YpointsM, t.MotherTraj.at(n).at(1) )
                ZpointsM = np.append( ZpointsM, t.MotherTraj.at(n).at(2) )
            plt.plot(ZpointsM,YpointsM,col(t.MotherPDG),linewidth=2)
        if (PartSteps > 0):
            for m in range(PartSteps):
                XpointsP = np.append( XpointsP, t.PartTraj.at(m).at(0) )
                YpointsP = np.append( YpointsP, t.PartTraj.at(m).at(1) )
                ZpointsP = np.append( ZpointsP, t.PartTraj.at(m).at(2) )
            plt.plot(ZpointsP,YpointsP,"b-",linewidth=3)
        '''


plt.title("MCParticle Tracks - YZ Projection - All Compton e- > 10 MeV (100 Events)",fontsize=20)
#plt.xlim([-100,1100])
plt.xlim([0-50,detLength+50])
plt.ylim([-detHalfHeight-50,detHalfHeight+50])
plt.xlabel("Z Position - Aligned with Beam Direction - [cm]",fontsize=20)
plt.ylabel("Y Position - Vertical Direction - [cm]",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()
