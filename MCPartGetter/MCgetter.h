/**
 * \file MCgetter.h
 *
 * \ingroup Analysis
 * 
 * \brief Class def header for a class MCgetter
 *
 * @author David Caratelli
 */

/** \addtogroup Analysis

    @{*/

#ifndef MCGETTER_H
#define MCGETTER_H

#include "Analysis/ana_base.h"
#include "TreeTop.h"
#include "LArUtil/Geometry.h"
#include <vector>
#include <stdexcept>

namespace larlite {
  /**
     \class MCgetter
     User custom analysis class made by David Caratelli
   */
  class MCgetter{
  
  public:

    /// Default constructor
    MCgetter(){ _Ecut=0; _searchPDG=false;};

    // Normal constructor
    MCgetter(event_mcpart *event_part) { _event_part = event_part; }

    /// Default destructor
    virtual ~MCgetter(){};

    /// set energy cut
    void SetECut(double E) { _Ecut=E; }

    /// Reset
    void Reset(event_mcpart *event_part) { _event_part = event_part, fillParticleMap(); setTreeTops(); };

    /// Fill Particle Map
    void fillParticleMap();

    /// Set TreeTops
    void setTreeTops();

    /// Add nodes
    void AddNodes(mcpart part, TreeNode& parentNode, int ancestor);

    /// Print Daughter Information for user-selected particle
    void drawDaughters(mcpart paticle);

    /// Search Particle Map: return mcpart index in event_part object
    int searchParticleMap(const int trackID);

    /// Search Particle Map: return mcpart index in event_part object
    int searchParticleMap(std::set<int>::iterator it){ return _particleMap.find(*it)->second; }

    /// Get Mother particle. returns element in event_part
    int getMother(const mcpart particle);

    /// Find particle's ancestor: return ancestor's TrackId
    int getAncestor(const mcpart particle);

    /// get TreeTops
    std::vector<TreeTop> getTreeTops() { return _treeTops; }

    /// find process
    std::vector<int> findProcess(const std::vector<std::pair<int, std::string> > procList);

    /// find process
    void findProcess(TreeNode currentNode, // Which node is currently under examination
		     int currentProc,      // Which step in the process is being searched for
		     int branchPosition,   // branch position being looked at
		     std::vector<TreeNode> nodePath, // List of nodes succesfully matched
		     std::vector<int> branchPath,    // List of branches followed
		     const std::vector<std::pair<int,std::string> > procList,
		     std::vector<int>& outputs);

    /// get Trajectory Points
    std::vector<std::vector<double> > getTrajectoryPoints(mcpart *part);

    /// get Trajectory Points - if in TPC (buffer denotes extra distance, in all directions, to save track)
    std::vector<std::vector<double> > getTrajectoryPointsInTPC(mcpart *part, double buffer);

    /// get Trajectory Points before TPC - if in TPC (buffer denotes extra distance, in all directions, to save track)
    std::vector<std::vector<double> > getTrajectoryPointsBeforeTPC(mcpart *part, double buffer, int& intpc);

    /// get Energy Points - if in TPC (buffer denotes extra distance, in all directions, to save track)
    std::vector<double> getEnergyPointsInTPC(mcpart *part, double buffer);

    /// get Start-End Trajectory Points
    std::vector<std::vector<double> > getStartEndTrajectory(mcpart *part);

    /// Get All particles with a certain PDG code - returns the TrackId
    void getAllPDGs(std::vector<int> pdgs) { _searchPDG = true; _getPDGs = pdgs; }
    
    /// Get PDG list
    std::vector<std::vector<int> > getPDGlist() { return _PDGlist; }

    /// Get TreeNodes list
    std::vector<std::vector<TreeNode> > getTreeNodelist() { return _TreeNodes; }

    /// find MC showers
    void findMCShowers(TreeNode node);

    /// make MC shower
    void makeMCShower(TreeNode node);

    /// get shower particles
    void getAllShowerParticles(TreeNode node, std::vector<int> &thisShower);

    /// get All MCShowers
    std::vector<std::vector<int> > getAllShowers() { return _MCShowers; }

    /// Find PDG list
    std::vector<int> findPDGlist(int pdgcode);

    /// Get String with particle's full process info
    std::string getFullProcess(mcpart part);
    
    /// set verbose
    void SetVerbose(bool on) { _verbose = on; }

    /// create a list of all processes involved that can be searched

    /// function that returns true/false if particle is in list

    protected:

    // Map for particle element in list and trackID
    // MAP: <trackID, Position in event_part>
    std::map<int,int> _particleMap;

    // A pointer to the event_particle list
    event_mcpart *_event_part;
    
    /// Cut on Energy for particle info to be displayed
    double _Ecut;

    /// Vector of TreeTops: one for each "primary" particle
    std::vector<TreeTop> _treeTops;

    /// Verobse flag
    bool _verbose;

    /// Bool for whether to get particles w/ certain PDG
    bool _searchPDG;

    /// PDG code of particle to be searched - list for all codes to be searched
    std::vector<int> _getPDGs;

    /// PDGlist: vector vector of particle trackIds that match pdg list being searched. One element of outer vector per PDG given
    std::vector<std::vector<int> > _PDGlist;

    /// PDGlist: vector vector of particle TreeNodes that match pdg list being searched. One element of outer vector per PDG given
    std::vector<std::vector<TreeNode> > _TreeNodes;

    /// List of MCShowers. Each element is a list of TrackIDs
    /// The first entry is the first charged particle in the shower
    std::vector<std::vector<int> > _MCShowers;

  };
}
#endif

/** @} */ // end of doxygen group 
