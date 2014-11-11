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
#plt.legend([ProfYZ,lineBlue,lineRed], ["TPC Boundaries","Max Deflection Point","Muons"],loc=4,fontsize=20)
plt.legend([ProfYZ,lineBlue], ["TPC Boundaries","Max Deflection Point"],loc=4,fontsize=20)

f = ROOT.TFile("cosmic_AllTrack.root")

m = f.Get("muon_tree")


muonEntries = m.GetEntries()
for j in range(200):
    
    m.GetEntry(j)

    XpointsMu = np.array([])
    YpointsMu = np.array([])
    ZpointsMu = np.array([])


    if (m.MuonTraj and (m.inTPC==1)):
        muonSteps = m.MuonTraj.size()
        isinside = 0
        exited = 0

        
        if (muonSteps > 0):

            '''
            XpointsMu = np.append( XpointsMu, m.MuonTraj.at(0).at(0) )
            YpointsMu = np.append( YpointsMu, m.MuonTraj.at(0).at(1) )
            ZpointsMu = np.append( ZpointsMu, m.MuonTraj.at(0).at(2) )

            XpointsMu = np.append( XpointsMu, m.maxDeflectionX )
            YpointsMu = np.append( YpointsMu, m.maxDeflectionY )
            ZpointsMu = np.append( ZpointsMu, m.maxDeflectionZ )

            XpointsMu = np.append( XpointsMu, m.MuonTraj.at(muonSteps-1).at(0) )
            YpointsMu = np.append( YpointsMu, m.MuonTraj.at(muonSteps-1).at(1) )
            ZpointsMu = np.append( ZpointsMu, m.MuonTraj.at(muonSteps-1).at(2) )
            '''
            #plt.plot(m.maxDeflectionZ,m.maxDeflectionY,"bo",markersize=4);
            #plt.plot(ZpointsMu,YpointsMu,"r-",linewidth=2)
            
            stepsin = 0
            
            Xenter = 0
            Yenter = 0
            Zenter = 0
            Xexit  = 0
            Yexit  = 0
            Zexit  = 0
            
            for n in range(muonSteps):
                XpointsMu = np.append( XpointsMu, m.MuonTraj.at(n).at(0) )
                YpointsMu = np.append( YpointsMu, m.MuonTraj.at(n).at(1) )
                ZpointsMu = np.append( ZpointsMu, m.MuonTraj.at(n).at(2) )
                if ( (m.MuonTraj.at(n).at(0) > 0) and (m.MuonTraj.at(n).at(0) < 2*detHalfWidth) and
                     (m.MuonTraj.at(n).at(1) > -detHalfHeight) and (m.MuonTraj.at(n).at(1) < detHalfHeight) and
                     (m.MuonTraj.at(n).at(2) > 0) and (m.MuonTraj.at(n).at(2) < detLength) ):
                    if (isinside == 0):
                        Xenter = m.MuonTraj.at(n).at(0)
                        Yenter = m.MuonTraj.at(n).at(1)
                        Zenter = m.MuonTraj.at(n).at(2)
                    stepsin += 1
                    isinside = 1
                else:
                    if ((isinside == 1) and (exited == 0)):
                        exited = 1
                        Xexit = m.MuonTraj.at(n).at(0)
                        Yexit = m.MuonTraj.at(n).at(1)
                        Zexit = m.MuonTraj.at(n).at(2)

            if (Xexit == 0):
                Xexit = m.MuonTraj.back().at(0)
                Yexit = m.MuonTraj.back().at(1)
                Zexit = m.MuonTraj.back().at(2)

            

            if (isinside == 1):
                print "Energy: {0}   Steps: {1}   StepsinTPC: {2}".format(m.muonE,muonSteps,stepsin)
                print "Entry: [{0}, {1}, {2}]  --> Exit: [{3}, {4}, {5}]".format(Xenter,Yenter,Zenter,Xexit,Yexit,Zexit)
                print "Distance: {0}".format(np.sqrt( (Xenter-Xexit)*(Xenter-Xexit) + (Yenter-Yexit)*(Yenter-Yexit) + (Zenter-Zexit)*(Zenter-Zexit) ) )
                print
                plt.plot(ZpointsMu,YpointsMu,"r-",linewidth=2)


plt.title("Muon Tracks - YZ Projection - All Muons that Enter TPC",fontsize=20)
#plt.xlim([0-1100,2*detHalfWidth+1100])
plt.xlim([0-1100,detLength+1100])
plt.ylim([-detHalfHeight-1100,detHalfHeight+1100])
plt.xlabel("Z Position - Aligned with Beam Direction - [cm]",fontsize=20)
#plt.xlabel("X Position - Aligned with Electric Field - [cm]",fontsize=20)
plt.ylabel("Y Position - Vertical Direction - [cm]",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()
