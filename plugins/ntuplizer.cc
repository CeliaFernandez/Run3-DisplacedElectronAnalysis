#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
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


struct electron {
  std::vector<float> pt, et, eta, phi;
  std::vector<float> dB, edB, ip3d;
  std::vector<float> sigmaIetaIphi, full5x5_sigmaIetaIphi;
  std::vector<float> ecalRegressionEnergy, ecalRegressionError, ecalTrackRegressionEnergy, ecalTrackRegressionError;
  std::vector<float> ecalScale, ecalSmear, ecalRegressionScale, ecalRegressionSmear;
  std::vector<bool> passConversionVeto;
  std::vector<float> trackIso, ecalIso, hcalIso; 
  std::vector<int> nSeeds, nTotalSeeds;
  std::vector<int> nBasicClusters, nBasicClustersWithSeed;
  std::vector<float> seedTime, clusterX, clusterY, clusterZ;
  std::vector<float> avgSeedTime, avgClusterX, avgClusterY, avgClusterZ;

  void clear() {
    pt.clear(); et.clear(); eta.clear(); phi.clear();
    dB.clear(); edB.clear(); ip3d.clear();
    sigmaIetaIphi.clear(); full5x5_sigmaIetaIphi.clear();
    ecalRegressionEnergy.clear(); ecalRegressionError.clear(); ecalTrackRegressionEnergy.clear(); ecalTrackRegressionError.clear();
    ecalScale.clear(); ecalSmear.clear(); ecalRegressionScale.clear(); ecalRegressionSmear.clear();
    passConversionVeto.clear();
    trackIso.clear(); ecalIso.clear(); hcalIso.clear();
    nSeeds.clear(); nTotalSeeds.clear();
    nBasicClusters.clear();
    nBasicClustersWithSeed.clear();
    seedTime.clear(); clusterX.clear(); clusterY.clear(); clusterZ.clear();
    avgSeedTime.clear(); avgClusterX.clear(); avgClusterY.clear(); avgClusterZ.clear();
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

      edm::EDGetTokenT<edm::View<pat::Electron> > lowPtElectronToken;
      edm::Handle<edm::View<pat::Electron> > lowPtElectrons;

      edm::EDGetTokenT<edm::View<pat::Electron> > electronToken;
      edm::Handle<edm::View<pat::Electron> > electrons;

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

      int nLowPtElectron;
      electron lowPtElecs;

      int nStdElectron;
      electron stdElecs;

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
   lowPtElectronToken = consumes<edm::View<pat::Electron> >  (parameters.getParameter<edm::InputTag>("lowPtElectronCollection"));
   electronToken = consumes<edm::View<pat::Electron> >  (parameters.getParameter<edm::InputTag>("electronCollection"));

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

   tree_out->Branch("nLowPtElectron", &nLowPtElectron);
   tree_out->Branch("lowPtElectron_pt", &lowPtElecs.pt);
   tree_out->Branch("lowPtElectron_et", &lowPtElecs.et);
   tree_out->Branch("lowPtElectron_eta", &lowPtElecs.eta);
   tree_out->Branch("lowPtElectron_phi", &lowPtElecs.phi);
   //tree_out->Branch("lowPtElectron_nBasicClusters", &lowPtElecs.nBasicClusters);
   //tree_out->Branch("lowPtElectron_nBasicClustersWithSeed", &lowPtElecs.nBasicClustersWithSeed);
   tree_out->Branch("lowPtElectron_seedTime", &lowPtElecs.seedTime);
   tree_out->Branch("lowPtElectron_nSeeds", &lowPtElecs.nSeeds);
   tree_out->Branch("lowPtElectron_nTotalSeeds", &lowPtElecs.nTotalSeeds);
   //tree_out->Branch("lowPtElectron_clusterX", &lowPtElecs.clusterX);
   //tree_out->Branch("lowPtElectron_clusterY", &lowPtElecs.clusterY);
   //tree_out->Branch("lowPtElectron_clusterZ", &lowPtElecs.clusterZ);
   //tree_out->Branch("lowPtElectron_avgSeedTime", &lowPtElecs.avgSeedTime);
   //tree_out->Branch("lowPtElectron_avgClusterX", &lowPtElecs.avgClusterX);
   //tree_out->Branch("lowPtElectron_avgClusterY", &lowPtElecs.avgClusterY);
   //tree_out->Branch("lowPtElectron_avgClusterZ", &lowPtElecs.avgClusterZ);
   tree_out->Branch("lowPtElectron_dB", &lowPtElecs.dB);
   tree_out->Branch("lowPtElectron_edB", &lowPtElecs.edB);
   tree_out->Branch("lowPtElectron_ip3d", &lowPtElecs.ip3d);
   tree_out->Branch("lowPtElectron_sigmaIetaIphi", &lowPtElecs.sigmaIetaIphi);
   tree_out->Branch("lowPtElectron_full5x5_sigmaIetaIphi", &lowPtElecs.full5x5_sigmaIetaIphi);
   tree_out->Branch("lowPtElectron_ecalRegressionEnergy", &lowPtElecs.ecalRegressionEnergy);
   tree_out->Branch("lowPtElectron_ecalRegressionError", &lowPtElecs.ecalRegressionError);
   tree_out->Branch("lowPtElectron_ecalTrackRegressionEnergy", &lowPtElecs.ecalTrackRegressionEnergy);
   tree_out->Branch("lowPtElectron_ecalTrackRegressionError", &lowPtElecs.ecalTrackRegressionError);
   tree_out->Branch("lowPtElectron_ecalScale", &lowPtElecs.ecalScale);
   tree_out->Branch("lowPtElectron_ecalSmear", &lowPtElecs.ecalSmear);
   tree_out->Branch("lowPtElectron_ecalRegressionScale", &lowPtElecs.ecalRegressionScale);
   tree_out->Branch("lowPtElectron_ecalRegressionSmear", &lowPtElecs.ecalRegressionSmear);
   tree_out->Branch("lowPtElectron_trackIso", &lowPtElecs.trackIso);
   tree_out->Branch("lowPtElectron_ecalIso", &lowPtElecs.ecalIso);
   tree_out->Branch("lowPtElectron_hcalIso", &lowPtElecs.hcalIso);
   tree_out->Branch("lowPtElectron_passConversionVeto", &lowPtElecs.passConversionVeto);

   tree_out->Branch("nStdElectron", &nStdElectron);
   tree_out->Branch("stdElectron_et", &stdElecs.et);
   tree_out->Branch("stdElectron_eta", &stdElecs.eta);
   tree_out->Branch("stdElectron_phi", &stdElecs.phi);
   tree_out->Branch("stdElectron_seedTime", &stdElecs.seedTime);

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
   iEvent.getByToken(lowPtElectronToken, lowPtElectrons);
   iEvent.getByToken(electronToken, electrons);
   ebRecHits_ = iEvent.get(reducedBarrelRecHitCollectionToken_);
   eeRecHits_ = iEvent.get(reducedEndcapRecHitCollectionToken_);
   //iEvent.getByToken(reducedBarrelRecHitCollectionToken_, ebRecHits_);
   //iEvent.getByToken(reducedEndcapRecHitCollectionToken_, eeRecHits_);

   // Clear all variables
   nSC = 0;
   SCs.clear();
   nLowPtElectron = 0;
   lowPtElecs.clear();
   nStdElectron = 0;
   stdElecs.clear();

   // Count number of events read
   counts->Fill(0);

   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();

   // Low Pt Electrons
   for (const auto& ele : *lowPtElectrons) {
     nLowPtElectron++;
     lowPtElecs.pt.push_back(ele.pt());
     lowPtElecs.et.push_back(ele.et());
     lowPtElecs.phi.push_back(ele.phi());
     lowPtElecs.eta.push_back(ele.eta());
     lowPtElecs.sigmaIetaIphi.push_back(ele.sigmaIetaIphi());
     lowPtElecs.full5x5_sigmaIetaIphi.push_back(ele.full5x5_sigmaIetaIphi());
     lowPtElecs.dB.push_back(ele.dB());
     lowPtElecs.edB.push_back(ele.edB());
     lowPtElecs.ip3d.push_back(ele.ip3d());
     lowPtElecs.trackIso.push_back(ele.trackIso());
     lowPtElecs.ecalIso.push_back(ele.ecalIso());
     lowPtElecs.hcalIso.push_back(ele.hcalIso());
     lowPtElecs.ecalRegressionEnergy.push_back(ele.ecalRegressionEnergy());
     lowPtElecs.ecalRegressionError.push_back(ele.ecalRegressionError());
     lowPtElecs.ecalTrackRegressionEnergy.push_back(ele.ecalTrackRegressionEnergy());
     lowPtElecs.ecalTrackRegressionError.push_back(ele.ecalTrackRegressionError());
     lowPtElecs.ecalScale.push_back(ele.ecalScale());
     lowPtElecs.ecalSmear.push_back(ele.ecalSmear());
     lowPtElecs.ecalRegressionScale.push_back(ele.ecalRegressionScale());
     lowPtElecs.ecalRegressionSmear.push_back(ele.ecalRegressionSmear());
     
     // Get the time
     float ttime = 0.;
     //float tcx = 0.;
     //float tcy = 0.;
     //float tcz = 0.;
     //int nbsc = 0;
     //int nbsc_seed = 0;
     int nseed = 0;
     reco::CaloClusterPtr seedPtr = ele.seed();
     reco::CaloCluster seed = *seedPtr;
     if (seed.size() > 0 ) {
       for (const auto &hit:seed.hitsAndFractions()) {
         DetId id = (hit).first;
         float fraction = (hit).first;
         EcalRecHitCollection *recHits = nullptr;
         if (id.subdetId() == EcalBarrel) {
           recHits = &ebRecHits_;
         } else if (id.subdetId() == EcalEndcap) {
           recHits = &eeRecHits_;
         } else {
           std::cout << "Invalid subdet, should never happen" << std::endl;
         }
         auto theSeedHit = recHits->find(id);
         float time = -99.;
         if (theSeedHit != recHits->end()) {
           time = (*theSeedHit).time();
           ttime = ttime + time;
           nseed++;
         }
       }
     }
     if (nseed>0)
       lowPtElecs.seedTime.push_back(ttime/nseed);
     else
       lowPtElecs.seedTime.push_back(-99.);
     lowPtElecs.nSeeds.push_back(nseed);
     lowPtElecs.nTotalSeeds.push_back(seed.size());
     //std::cout << "Seed size: " << seed.size() << std::endl;
     //for (unsigned int i=0; i < seed.size(); i++)
     //  std::cout << seed.printHitAndFraction(i) << std::endl;

     //std::vector<reco::CaloCluster> clusters = ele.basicClusters();
     //lowPtElecs.nBasicClusters.push_back(clusters.size());
     /*
     if (clusters.size() > 0 ) {
       for (const auto& cluster:clusters) {
         DetId id = (cluster.hitsAndFractions()[0]).first;
         float cx = cluster.position().x();
         float cy = cluster.position().y();
         float cz = cluster.position().z();
         EcalRecHitCollection *recHits = nullptr;
         if (id.subdetId() == EcalBarrel) {
           recHits = &ebRecHits_;
         } else if (id.subdetId() == EcalEndcap) {
           recHits = &eeRecHits_;
         } else {
           std::cout << "Invalid subdet, should never happen" << std::endl;
         }

         auto theSeedHit = recHits->find(id);
         float time = -99.;
         if (theSeedHit != recHits->end()) {
           //lowPtElecs.seedTime.push_back((*theSeedHit).time());
           time = (*theSeedHit).time();
           ttime = ttime + time;
           tcx = tcx + cx;
           tcy = tcy + cy;
           tcz = tcz + cz;
           nbsc_seed++;
         }
         nbsc++;
         
         // Save first seed
         if (nbsc_seed==1 and (lowPtElecs.seedTime.size() < lowPtElecs.pt.size())) {
           lowPtElecs.seedTime.push_back(time);
           lowPtElecs.clusterX.push_back(cx);
           lowPtElecs.clusterY.push_back(cy);
           lowPtElecs.clusterZ.push_back(cz);
         }
       }
     } else {
       lowPtElecs.seedTime.push_back(-99.);
     }
     */
     // Save average cluster seed props
    // lowPtElecs.avgSeedTime.push_back(ttime/nbsc_seed);
    // lowPtElecs.avgClusterX.push_back(tcx/nbsc_seed);
    // lowPtElecs.avgClusterY.push_back(tcy/nbsc_seed);
    // lowPtElecs.avgClusterZ.push_back(tcz/nbsc_seed);
   }

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
