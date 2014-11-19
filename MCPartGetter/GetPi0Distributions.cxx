#ifndef GETPI0DISTRIBUTIONS_CXX
#define GETPI0DISTRIBUTIONS_CXX

#include "GetPi0Distributions.h"

namespace larlite {

  bool GetPi0Distributions::initialize() {

    if (!_tree){
      _tree = new TTree("ana_tree","");
      
      _tree->Branch("_totE",&_totE,"totE/D");
      _tree->Branch("_Nshowers",&_Nshowers,"NShowers/I");
      _tree->Branch("_ShowerE1",&_ShowerE1,"ShowerE1/D");
      _tree->Branch("_ShowerE2",&_ShowerE2,"ShowerE2/D");
      _tree->Branch("_Angle",&_Angle,"Angle/D");
      _tree->Branch("_VtxDist1",&_VtxDist1,"VtxDist1/D");
      _tree->Branch("_VtxDist2",&_VtxDist2,"VtxDist2/D");
    }	

    return true;
  }
  
  bool GetPi0Distributions::analyze(storage_manager* storage) {

    clearTree();

    // get MCShowers
    auto evt_mcshower = storage->get_data<event_mcshower>("mcreco");

    if (!evt_mcshower){
      std::cout << "No showers of this type!" << std::endl;
      return false;
    }
      
  
    _Nshowers = evt_mcshower->size();

    if (evt_mcshower->size() == 2){
      mcshower shr1 = evt_mcshower->at(0);
      mcshower shr2 = evt_mcshower->at(1);

      mcstep start1 = shr1.Start();
      mcstep start2 = shr2.Start();
      
      mcstep dep1 = shr1.DetProfile();
      mcstep dep2 = shr2.DetProfile();

      // a pi0! continue
      _ShowerE1 = dep1.E();
      _ShowerE2 = dep2.E();

      _totE = _ShowerE1+_ShowerE2;
      
      _VtxDist1 = sqrt( (dep1.X() - start1.X())*(dep1.X() - start1.X()) +
			(dep1.Y() - start1.Y())*(dep1.Y() - start1.Y()) +
			(dep1.Z() - start1.Z())*(dep1.Z() - start1.Z())  );

      _VtxDist2 = sqrt( (dep2.X() - start2.X())*(dep2.X() - start2.X()) +
			(dep2.Y() - start2.Y())*(dep2.Y() - start2.Y()) +
			(dep2.Z() - start2.Z())*(dep2.Z() - start2.Z())  );

      double mom1 = sqrt( dep1.Px()*dep1.Px() +
			  dep1.Py()*dep1.Py() +
			  dep1.Pz()*dep1.Pz() );
			    
      double mom2 = sqrt( dep2.Px()*dep2.Px() +
			  dep2.Py()*dep2.Py() +
			  dep2.Pz()*dep2.Pz() );


      _Angle = (180./3.14)*acos( (dep1.Px()*dep2.Px() + dep1.Py()*dep2.Py() + dep1.Pz()*dep2.Pz())/(mom1*mom2) );

    }

    if (_tree)
      _tree->Fill();

    return true;
  }

  bool GetPi0Distributions::finalize() {

    if (_tree)
      _tree->Write();

    return true;
  }

  void GetPi0Distributions::clearTree() {

    _Nshowers = -1;
    _totE = -1;
    _ShowerE1 = -1;
    _ShowerE2 = -1;
    _Angle = -1;
    _VtxDist1 = -1;
    _VtxDist2 = -1;
    
    return;
  }

}
#endif
