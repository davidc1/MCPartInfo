import os, ROOT, sys
from ROOT import *
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D

def col(x):
    c = "k"
    if (x == 11):
        c = "b"
    if (x == -11):
        c = "r"

    return c

#find shower length from energy
#E in MeV
def findLength(E,pdg):
    a = 0
    if (pdg==11 or pdg==-11):
        a = 0.5
    else:
        a = 1
        
    t = np.log((E/1000.)/5.7)
    
    L = t - a + 0.08*18 + 9.6

    return L

#in 2 dimensions
def getTriangle(start,end):

    slope = (end[1]-start[1])/(end[0]-start[0])
    antislope = -1/slope

    theta = np.arctan(antislope)
    '''
    intercept = end[1]-antislope*end[0]

    mom = np.array([end[0]-start[0], end[1]-start[1]])
    dirmag = np.sqrt( mom[0]*mom[0] + mom[1]*mom[1] )
    
    angleWithX = np.arccos(mom[0]/dirmag)
    
    theta = np.pi/2. - angleWithX
    '''
    point1 = np.array( [end[0]+5*np.cos(theta), end[1]+5*np.sin(theta)] )
    point2 = np.array( [end[0]-5*np.cos(theta), end[1]-5*np.sin(theta)] )

    tri = np.array([ start, point1, point2 ])

    return tri 

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

#fig = plt.figure(figsize=(12,10))
#plt.ion()


#legend:
lineBlue = Line2D([],[],color='blue',linewidth=5)
lineRed = Line2D([],[],color='red',linewidth=5)
lineBlack = Line2D([],[],color='black',linewidth=5)


f = ROOT.TFile("shrBackground.root")

m = f.Get("ana_tree")

colorlist = ['b','g','y','k','c','r','w']
color = 0
event = 0

showerEntries = m.GetEntries()
for j in range(showerEntries):
    
    m.GetEntry(j)

    if ( (m.inActiveVolume == 1) and (m.E > 100) ):

        #get shower start point and momentum
        shrX = m.X
        shrY = m.Y
        shrZ = m.Z
        
        shrPx = m.Px
        shrPy = m.Py
        shrPz = m.Pz
        
        shrMom = np.sqrt(m.Px*m.Px + m.Py*m.Py + m.Pz*m.Pz)

        shrEndX = shrX + (shrPx/shrMom)*(m.E/10.)
        shrEndY = shrY + (shrPy/shrMom)*(m.E/10.)
        shrEndZ = shrZ + (shrPz/shrMom)*(m.E/10.)
        '''
        shrEndX = shrX + (shrPx/shrMom)*findLength(m.E,m.PDG)
        shrEndY = shrY + (shrPy/shrMom)*findLength(m.E,m.PDG)
        shrEndZ = shrZ + (shrPz/shrMom)*findLength(m.E,m.PDG)
        '''

        shrStart = np.array([shrZ,shrY])
        shrEnd = np.array([shrEndZ,shrEndY])

        triangle = getTriangle(shrStart,shrEnd)
        
        triZ = np.array( [ triangle[0][0], triangle[1][0], triangle[2][0] ] )
        triY = np.array( [ triangle[0][1], triangle[1][1], triangle[2][1] ] )
        #        plt.plot(triZ, triY, "b", markersize=2)
        plt.plot(shrZ, shrY, "ko")
        plt.plot(shrEndZ, shrEndY, "ko")
        plt.gca().add_patch(Polygon(triangle,closed=True,fill=True))
        
        #patch = patches.PathPatch(getTriangle(shrStart,shrEnd),facecolor='b',lw=2)
        #plt.gca().add_patch(patch)
        '''
        shrTrajX = np.array([shrX,shrEndX])
        shrTrajY = np.array([shrY,shrEndY])
        shrTrajZ = np.array([shrZ,shrEndZ])

        plt.plot(shrTrajZ,shrTrajY,'b',markersize=2)
        '''
        color += 1

        #get shower ancestor track
        ancestorTrack = m.AncestorTraj
        if (ancestorTrack.size() > 0):
            
            XpointsAnc = np.array([])
            YpointsAnc = np.array([])
            ZpointsAnc = np.array([])
            
            steps = ancestorTrack.size()

            for x in range(steps):

                XpointsAnc = np.append( XpointsAnc, ancestorTrack.at(x).at(0) )
                YpointsAnc = np.append( YpointsAnc, ancestorTrack.at(x).at(1) )
                ZpointsAnc = np.append( ZpointsAnc, ancestorTrack.at(x).at(2) )

            plt.plot(ZpointsAnc,YpointsAnc,'r',markersize=2)



plt.grid(1,which='both')
plt.title("Muon Tracks - YZ Projection - All Muons that Enter TPC",fontsize=20)
#plt.xlim([0-1100,2*detHalfWidth+1100])
plt.xlim([0-50,detLength+50])
plt.ylim([-detHalfHeight-50,detHalfHeight+50])
plt.xlabel("Z Position - Aligned with Beam Direction - [cm]",fontsize=20)
#plt.xlabel("X Position - Aligned with Electric Field - [cm]",fontsize=20)
plt.ylabel("Y Position - Vertical Direction - [cm]",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()

