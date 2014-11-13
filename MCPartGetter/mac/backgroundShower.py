import os, ROOT, sys
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.path import Path
from matplotlib.lines import Line2D

def col(x):
    c = "k"
    if (x == 11):
        c = "b"
    if (x == -11):
        c = "r"

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

f = ROOT.TFile("makeshowers.root")

m = f.Get("shower_tree")

colorlist = ['b','g','y','k','c','r','w']
color = 0

showerEntries = m.GetEntries()
for j in range(showerEntries):
    
    m.GetEntry(j)

    color += 1

    if ( m.ShowerTraj and (m.inTPC==1) ):
        showerTracks = m.ShowerTraj.size()
        
        for x in range(showerTracks):

            thistrack = m.ShowerTraj.at(x)

            XpointsShower = np.array([])
            YpointsShower = np.array([])
            ZpointsShower = np.array([])

            for n in range(thistrack.size()):

                XpointsShower = np.append( XpointsShower, thistrack.at(n).at(0) )
                YpointsShower = np.append( YpointsShower, thistrack.at(n).at(1) )
                ZpointsShower = np.append( ZpointsShower, thistrack.at(n).at(2) )

            plt.plot(ZpointsShower,YpointsShower,colorlist[color%7],markersize=2);


plt.title("Muon Tracks - YZ Projection - All Muons that Enter TPC",fontsize=20)
#plt.xlim([0-1100,2*detHalfWidth+1100])
plt.xlim([0-50,detLength+50])
plt.ylim([-detHalfHeight-50,detHalfHeight+50])
plt.xlabel("Z Position - Aligned with Beam Direction - [cm]",fontsize=20)
#plt.xlabel("X Position - Aligned with Electric Field - [cm]",fontsize=20)
plt.ylabel("Y Position - Vertical Direction - [cm]",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()
