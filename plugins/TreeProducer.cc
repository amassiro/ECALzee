// -*- C++ -*-
//
// Package:    EcalZee
// Class:      TreeProducer
// 
/**\class TreeProducer TreeProducer.cc EcalZee/plugins/TreeProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Massironi
//         Created:  Thu, 17 Nov 2016 10:09:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



// ECAL specific

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"


#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"



#include "FWCore/Framework/interface/ESHandle.h"



#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"   //     ----> CLHEP/Geometry/Point3D.h issue
// #include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"    ----> CLHEP/Geometry/Point3D.h issue
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
// #include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"      ----> CLHEP/Geometry/Point3D.h issue
// #include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"      ----> CLHEP/Geometry/Point3D.h issue




// #include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
// #include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"







//---- for TP

// #include "CondFormats/EcalObjects/interface/EcalTPGLutIdMap.h"
// #include "CondFormats/EcalObjects/interface/EcalTPGLutGroup.h"
// #include "CondFormats/EcalObjects/interface/EcalTPGPhysicsConst.h"
// #include "CondFormats/DataRecord/interface/EcalTPGLutIdMapRcd.h"
// #include "CondFormats/DataRecord/interface/EcalTPGLutGroupRcd.h"
// #include "CondFormats/DataRecord/interface/EcalTPGPhysicsConstRcd.h"
// 
// 
// #include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"






#include "TTree.h"




//---- for Zee part
#include "DataFormats/PatCandidates/interface/Electron.h"






//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TreeProducer(const edm::ParameterSet&);
      ~TreeProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      edm::InputTag _vtxTag;
      edm::InputTag _EleTag;
      edm::InputTag _Supercluster_EB_Tag;
      edm::InputTag _Supercluster_EE_Tag;
      
      edm::EDGetTokenT<reco::VertexCollection> _Token_vtxTag;
      edm::EDGetTokenT<std::vector<pat::Electron> > _Token_EleTag;
//       edm::EDGetTokenT<std::vector<reco::PFCandidate> > _Token_EleTag;
      edm::EDGetTokenT<reco::SuperClusterCollection> _Token_Supercluster_EB_Tag;
      edm::EDGetTokenT<reco::SuperClusterCollection> _Token_Supercluster_EE_Tag;
      
      
      
      TTree *_outTree;
      
      UInt_t _run;
      UShort_t _lumi;
      UShort_t _bx;
      UShort_t _event;      
   

      std::vector<float> _std_vector_Ele_pt;
      std::vector<float> _std_vector_Ele_eta;
      std::vector<float> _std_vector_Ele_phi;
      std::vector<float> _std_vector_Ele_r9;
      std::vector<float> _std_vector_Ele_sigmaIetaIeta;
      std::vector<float> _std_vector_Ele_sigmaIphiIphi;
      std::vector<float> _std_vector_Ele_dr03EcalRecHitSumEt;
      
      std::vector<float> _std_vector_SC_EB_raw_et;
      std::vector<float> _std_vector_SC_EB_et;
      std::vector<float> _std_vector_SC_EB_eta;
      std::vector<float> _std_vector_SC_EB_phi;

      std::vector<float> _std_vector_SC_EE_raw_et;
      std::vector<float> _std_vector_SC_EE_et;
      std::vector<float> _std_vector_SC_EE_eta;
      std::vector<float> _std_vector_SC_EE_phi;
      
      float _mll;
      
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TreeProducer::TreeProducer(const edm::ParameterSet& iConfig)

{
  
  _vtxTag = iConfig.getParameter<edm::InputTag>("vtxTag");
  _EleTag = iConfig.getParameter<edm::InputTag>("EleTag"); 
  _Supercluster_EB_Tag = iConfig.getParameter<edm::InputTag>("SuperClusterEBTag"); 
  _Supercluster_EE_Tag = iConfig.getParameter<edm::InputTag>("SuperClusterEETag"); 
  
  
  _Token_vtxTag   = consumes<reco::VertexCollection>(_vtxTag);
  _Token_EleTag   = consumes<std::vector<pat::Electron> >(_EleTag);
//   _Token_EleTag   = consumes<std::vector<reco::PFCandidate> >(_EleTag);
  _Token_Supercluster_EB_Tag = consumes<reco::SuperClusterCollection> (_Supercluster_EB_Tag);
  _Token_Supercluster_EE_Tag = consumes<reco::SuperClusterCollection> (_Supercluster_EE_Tag);
  
 
   //now do what ever initialization is needed
   usesResource("TFileService");
   
   
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   
   _outTree = fs->make<TTree>("tree","tree");
   
   _outTree->Branch("run",               &_run,             "run/i");
   _outTree->Branch("lumi",              &_lumi,            "lumi/s");
   _outTree->Branch("bx",                &_bx,              "bx/s");
   _outTree->Branch("event",             &_event,           "event/i");
 
   
   _outTree -> Branch("std_vector_Ele_pt"      , "std::vector<float>",   &_std_vector_Ele_pt);
   _outTree -> Branch("std_vector_Ele_eta"     , "std::vector<float>",   &_std_vector_Ele_eta);
   _outTree -> Branch("std_vector_Ele_phi"     , "std::vector<float>",   &_std_vector_Ele_phi);
   _outTree -> Branch("std_vector_Ele_r9"                      , "std::vector<float>",   &_std_vector_Ele_r9);
   _outTree -> Branch("std_vector_Ele_sigmaIetaIeta"           , "std::vector<float>",   &_std_vector_Ele_sigmaIetaIeta);
   _outTree -> Branch("std_vector_Ele_sigmaIphiIphi"           , "std::vector<float>",   &_std_vector_Ele_sigmaIphiIphi);
   _outTree -> Branch("std_vector_Ele_dr03EcalRecHitSumEt"     , "std::vector<float>",   &_std_vector_Ele_dr03EcalRecHitSumEt);
   
   _outTree -> Branch("std_vector_SC_EB_raw_et"  , "std::vector<float>",   &_std_vector_SC_EB_raw_et);
   _outTree -> Branch("std_vector_SC_EB_et"      , "std::vector<float>",   &_std_vector_SC_EB_et);
   _outTree -> Branch("std_vector_SC_EB_eta"     , "std::vector<float>",   &_std_vector_SC_EB_eta);
   _outTree -> Branch("std_vector_SC_EB_phi"     , "std::vector<float>",   &_std_vector_SC_EB_phi);

   _outTree -> Branch("std_vector_SC_EE_raw_et"  , "std::vector<float>",   &_std_vector_SC_EE_raw_et);
   _outTree -> Branch("std_vector_SC_EE_et"      , "std::vector<float>",   &_std_vector_SC_EE_et);
   _outTree -> Branch("std_vector_SC_EE_eta"     , "std::vector<float>",   &_std_vector_SC_EE_eta);
   _outTree -> Branch("std_vector_SC_EE_phi"     , "std::vector<float>",   &_std_vector_SC_EE_phi);
   
   _outTree -> Branch("mll",&_mll,"mll/F");
   
   
}


TreeProducer::~TreeProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  _run = iEvent.eventAuxiliary().run();
  _lumi = iEvent.eventAuxiliary().luminosityBlock();
  _bx = iEvent.eventAuxiliary().bunchCrossing();
  _event = iEvent.eventAuxiliary().event();
  
  
  
  ///---- Vertexes
  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByToken(_Token_vtxTag,vtxH);
  
  
  //---- Electrons
  edm::Handle<std::vector <pat::Electron> > electronHandle;
//   edm::Handle<std::vector <reco::PFCandidate> > electronHandle;
  iEvent.getByToken(_Token_EleTag,electronHandle);
  
  std::vector<pat::Electron> electrons = *electronHandle;
//   std::vector<reco::PFCandidate> electrons = *electronHandle;
  
  _std_vector_Ele_pt.clear();
  _std_vector_Ele_eta.clear();
  _std_vector_Ele_phi.clear();
  _std_vector_Ele_r9                  .clear();
  _std_vector_Ele_sigmaIetaIeta       .clear();
  _std_vector_Ele_sigmaIphiIphi       .clear();
  _std_vector_Ele_dr03EcalRecHitSumEt .clear();
  
  
  
  
  //---- the first two are used to build mll
  _mll = -999;
  bool notfound = true;
  reco::Candidate::LorentzVector L1;
  for ( unsigned int i=0; i<electrons.size(); ++i ){
    pat::Electron electron = electrons.at(i);
//     reco::PFCandidate electron = electrons.at(i);
    
    _std_vector_Ele_pt.push_back(electron.pt());
    _std_vector_Ele_eta.push_back(electron.eta());
    _std_vector_Ele_phi.push_back(electron.phi());
    _std_vector_Ele_r9                 .push_back(electron.r9                 ());
    _std_vector_Ele_sigmaIetaIeta      .push_back(electron.sigmaIetaIeta      ());
    _std_vector_Ele_sigmaIphiIphi      .push_back(electron.sigmaIphiIphi      ());
    _std_vector_Ele_dr03EcalRecHitSumEt.push_back(electron.dr03EcalRecHitSumEt());
    
    
    if (fabs(electron.eta()) < 1.5 && notfound) {
      L1 = electron.p4();
      notfound = false;
//       break;
    }    
  }
  
  
  reco::Candidate::LorentzVector L2;
  
  
  _std_vector_SC_EB_raw_et.clear();
  _std_vector_SC_EB_et.clear();
  _std_vector_SC_EB_eta.clear();
  _std_vector_SC_EB_phi.clear();

  _std_vector_SC_EE_raw_et.clear();
  _std_vector_SC_EE_et.clear();
  _std_vector_SC_EE_eta.clear();
  _std_vector_SC_EE_phi.clear();
  
  
  // Super Clusters
  // ... barrel
  edm::Handle<reco::SuperClusterCollection> superClusters_EB_Handle;
  iEvent.getByToken( _Token_Supercluster_EB_Tag, superClusters_EB_Handle );
  const reco::SuperClusterCollection* theBarrelSuperClusters = superClusters_EB_Handle.product () ;
  if ( ! superClusters_EB_Handle.isValid() ) {
    std::cerr << "EcalValidation::analyze --> superClusters_EB_Handle not found" << std::endl; 
  }
  else {
       
      for (reco::SuperClusterCollection::const_iterator itSC = theBarrelSuperClusters->begin(); itSC != theBarrelSuperClusters->end(); ++itSC ) {
        double scEt = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
        double scRawEt = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));
        
        _std_vector_SC_EB_raw_et.push_back(scRawEt);
        _std_vector_SC_EB_et.push_back(scEt);
        _std_vector_SC_EB_eta.push_back(itSC->position().eta());
        _std_vector_SC_EB_phi.push_back(itSC->position().phi());
      }
        
  }
  

  // ... endcap
  edm::Handle<reco::SuperClusterCollection> superClusters_EE_Handle;
  iEvent.getByToken( _Token_Supercluster_EE_Tag, superClusters_EE_Handle );
  const reco::SuperClusterCollection* theEndcapSuperClusters = superClusters_EE_Handle.product () ;
  if ( ! superClusters_EE_Handle.isValid() ) {
    std::cerr << "EcalValidation::analyze --> superClusters_EE_Handle not found" << std::endl; 
  }
  else {
    
    for (reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin(); itSC != theEndcapSuperClusters->end(); ++itSC ) {
      double scEt = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
      double scRawEt = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));
      
      _std_vector_SC_EE_raw_et.push_back(scRawEt);
      _std_vector_SC_EE_et.push_back(scEt);
      _std_vector_SC_EE_eta.push_back(itSC->position().eta());
      _std_vector_SC_EE_phi.push_back(itSC->position().phi());
      
      
      if (_mll == -999 && fabs(itSC->position().eta()) > 2.5) {
        reco::Candidate::LorentzVector L2;
        L2.SetE(itSC -> energy());
        L2.SetPx(scEt * cos (itSC->position().phi()));
        L2.SetPy(scEt * sin (itSC->position().phi()));
        L2.SetPz(itSC -> energy() * cos (  2.*atan( exp(- itSC->position().eta() ))   ));
        _mll = (L1+L2).mass();
        if (_mll < 10) _mll = -999;
      }
      
      
      
    }
    
  }
  
  
//   int nSCcleaned = 0;
//   int nSC = 0;
//   
//   for (reco::SuperClusterCollection::const_iterator itSC = theBarrelSuperClusters->begin(); 
//        itSC != theBarrelSuperClusters->end(); ++itSC ) {
//     
//     double scEt = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
//   double scRawEt = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));
  
  
  
  
  
  
  _outTree->Fill();
  
   
}




// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer);
