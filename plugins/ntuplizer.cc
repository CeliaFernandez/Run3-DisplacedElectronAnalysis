#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"


struct SC {
  std::vector<float> et, eta, phi;
  std::vector<float> seedTime;

  void clear() {
    et.clear(); eta.clear(); phi.clear();
    seedTime.clear();
  };

};


class ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ntuplizer(const edm::ParameterSet&);
      ~ntuplizer();

      edm::ConsumesCollector iC = consumesCollector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::ParameterSet parameters;

      //
      // --- Tokens and Handles
      //

      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone> > triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>  triggerPrescales_;

      edm::EDGetTokenT<edm::View<pat::Photon> > photonToken;
      edm::Handle<edm::View<pat::Photon> > photons;

      edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollectionToken_;
      edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollectionToken_;
      // edm::Handle<EcalRecHitCollection> ebRecHits_;
      // edm::Handle<EcalRecHitCollection> ecRecHits_;

      EcalRecHitCollection ebRecHits_;
      EcalRecHitCollection eeRecHits_;

      //const EcalClusterLazyTools::ESGetTokens ecalClusterToolsESGetTokens_;

      //
      // --- Variables
      //

      bool isData = false;

      // Event
      Int_t event = 0;
      Int_t lumiBlock = 0;
      Int_t run = 0;

      // Branch variables
      int nSC;
      SC SCs;

      //
      // --- Output
      //
      std::string output_filename;
      TH1F *counts;
      TFile *file_out;
      TTree *tree_out;

};

// Constructor
// ntuplizer::ntuplizer(const edm::ParameterSet& iConfig)  {
ntuplizer::ntuplizer(const edm::ParameterSet& iConfig) 
  : reducedBarrelRecHitCollectionToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection"))),
    reducedEndcapRecHitCollectionToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection"))) {
//    ecalClusterToolsESGetTokens_{iC} {

   std::cout << "Constructor" << std::endl;
   usesResource("TFileService");

   parameters = iConfig;

   counts = new TH1F("counts", "", 1, 0, 1);

   isData = parameters.getParameter<bool>("isData");
   photonToken = consumes<edm::View<pat::Photon> >  (parameters.getParameter<edm::InputTag>("photonCollection"));

}


// Destructor
ntuplizer::~ntuplizer() {
}


// beginJob (Before first event)
void ntuplizer::beginJob() {

   std::cout << "Begin Job" << std::endl;

   // Init the file and the TTree
   output_filename = parameters.getParameter<std::string>("nameOfOutput");
   file_out = new TFile(output_filename.c_str(), "RECREATE");
   tree_out = new TTree("Events", "Events");

   // Analyzer parameters
   isData = parameters.getParameter<bool>("isData");


   // TTree branches
   tree_out->Branch("event", &event, "event/I");
   tree_out->Branch("lumiBlock", &lumiBlock, "lumiBlock/I");
   tree_out->Branch("run", &run, "run/I");

   tree_out->Branch("nSC", &nSC);
   tree_out->Branch("SC_et", &SCs.et);
   tree_out->Branch("SC_eta", &SCs.eta);
   tree_out->Branch("SC_phi", &SCs.phi);
   tree_out->Branch("SC_seedTime", &SCs.seedTime);


}

// endJob (After event loop has finished)
void ntuplizer::endJob()
{

    std::cout << "End Job" << std::endl;
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();

}


// fillDescriptions
void ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// Analyze (per event)
void ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   iEvent.getByToken(photonToken, photons);
   ebRecHits_ = iEvent.get(reducedBarrelRecHitCollectionToken_);
   eeRecHits_ = iEvent.get(reducedEndcapRecHitCollectionToken_);
   //iEvent.getByToken(reducedBarrelRecHitCollectionToken_, ebRecHits_);
   //iEvent.getByToken(reducedEndcapRecHitCollectionToken_, eeRecHits_);

   // Clear all variables
   nSC = 0;
   SCs.clear();

   // Count number of events read
   counts->Fill(0);

   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();


   // -> Superclusters
   //EcalClusterLazyTools lazyTool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
   //edm::InputTag reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
   //edm::InputTag reducedEERecHitCollection(string("reducedEcalRecHitsEE"));
   //EcalClusterLazyTools lazyTools(iEvent, ecalClusterToolsESGetTokens_.get(iSetup), reducedBarrelRecHitCollectionToken_, reducedEndcapRecHitCollectionToken_);
   for (const auto& pho : *photons) {
     nSC++;
     SCs.et.push_back(pho.et());
     SCs.eta.push_back(pho.eta());
     SCs.phi.push_back(pho.phi());

     reco::SuperClusterRef clusterRef = pho.superCluster();
     reco::SuperCluster cluster = *clusterRef;
     if (cluster.size() > 0) {
       DetId id = (cluster.hitsAndFractions()[0]).first;
       EcalRecHitCollection *recHits = nullptr;
       if (id.subdetId() == EcalBarrel) {
         recHits = &ebRecHits_;
       } else if (id.subdetId() == EcalEndcap) {
         recHits = &eeRecHits_;
       } else {
         std::cout << "Invalid subdet, should never happen" << std::endl;
       }

       auto theSeedHit = recHits->find(id);
       if (theSeedHit != recHits->end())
         SCs.seedTime.push_back((*theSeedHit).time());
       else
         SCs.seedTime.push_back(-99.);
     } else {
       SCs.seedTime.push_back(-99.);
     }

     // Get the time
     //std::cout << lazyTools.SuperClusterTime(*pho.superCluster(), iEvent) << std::endl;
     /*
     DetId seed = (pho.superCluster()->seed()->hitsAndFractions())[0].first;
     bool isBarrel = seed.subdetId() == EcalBarrel;
     const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
     EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
     if (theSeedHit != rechits->end()) {
     phoSeedTime_  .push_back((*theSeedHit).time());
       phoSeedEnergy_.push_back((*theSeedHit).energy());
     } else{
       phoSeedTime_  .push_back(-99.);
       phoSeedEnergy_.push_back(-99.);
     }
     */


   }

   // -> Fill tree
   tree_out->Fill();

}

DEFINE_FWK_MODULE(ntuplizer);
