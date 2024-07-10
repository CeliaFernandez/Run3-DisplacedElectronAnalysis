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

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

struct genPart {
  std::vector<float> pt, eta, phi, m;
  std::vector<float> vx, vy, vz, lxy;
  std::vector<int> status, index, pdgId, motherIndex, motherPdgId;

  void clear() {
    pt.clear(); eta.clear(); phi.clear(); m.clear();
    vx.clear(); vy.clear(); vz.clear(); lxy.clear();
    status.clear(); index.clear(); pdgId.clear(); motherIndex.clear(); motherPdgId.clear();
  }
};

struct PV {
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  //double x_temp; double y_temp; double z_temp; 

  void clear() {
    x.clear(); y.clear(); z.clear();
    //x_temp = 0; y_temp = 0; z_temp = 0;
  };
};


struct SC {
  std::vector<float> et, eta, phi;
  std::vector<float> seedTime;

  void clear() {
    et.clear(); eta.clear(); phi.clear();
    seedTime.clear();
  };

};


struct electron {
  std::vector<float> charge;
  std::vector<float> Lxy_SV;
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
    charge.clear(); Lxy_SV.clear();
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
      edm::EDGetTokenT<edm::View<reco::Vertex> > primaryVertices;
      edm::Handle<edm::View<reco::Vertex> > pVtxs;

      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleToken;
      edm::Handle<edm::View<reco::GenParticle> > genParticles;

      edm::EDGetTokenT<edm::View<pat::Photon> > photonToken;
      edm::Handle<edm::View<pat::Photon> > photons;

      edm::EDGetTokenT<edm::View<pat::Electron> > lowPtElectronToken;
      edm::Handle<edm::View<pat::Electron> > lowPtElectrons;

      edm::EDGetTokenT<edm::View<pat::Electron> > electronToken;
      edm::Handle<edm::View<pat::Electron> > electrons;

      const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbESToken_;

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
      int nGP;
      genPart genParts;

      int nSC;
      SC SCs;

      int nLowPtElectron;
      electron lowPtElecs;

      int nStdElectron;
      electron stdElecs;

      PV PVs;

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
  : ttbESToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    reducedBarrelRecHitCollectionToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection"))),
    reducedEndcapRecHitCollectionToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection"))) {
//    ecalClusterToolsESGetTokens_{iC} {

   std::cout << "Constructor" << std::endl;
   usesResource("TFileService");

   parameters = iConfig;

   counts = new TH1F("counts", "", 1, 0, 1);

   isData = parameters.getParameter<bool>("isData");
   genParticleToken = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genCollection"));
   photonToken = consumes<edm::View<pat::Photon> >  (parameters.getParameter<edm::InputTag>("photonCollection"));
   lowPtElectronToken = consumes<edm::View<pat::Electron> >  (parameters.getParameter<edm::InputTag>("lowPtElectronCollection"));
   electronToken = consumes<edm::View<pat::Electron> >  (parameters.getParameter<edm::InputTag>("electronCollection"));
   primaryVertices = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("primaryVertexCollection"));

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

   tree_out->Branch("genParticle_pt", &genParts.pt);
   tree_out->Branch("genParticle_eta", &genParts.eta);
   tree_out->Branch("genParticle_phi", &genParts.phi);
   tree_out->Branch("genParticle_m", &genParts.m);
   tree_out->Branch("genParticle_vx", &genParts.vx);
   tree_out->Branch("genParticle_vy", &genParts.vy);
   tree_out->Branch("genParticle_vz", &genParts.vz);
   tree_out->Branch("genParticle_lxy", &genParts.lxy);
   tree_out->Branch("genParticle_status", &genParts.status);
   tree_out->Branch("genParticle_index", &genParts.index);
   tree_out->Branch("genParticle_pdgId", &genParts.pdgId);
   tree_out->Branch("genParticle_motherIndex", &genParts.motherIndex);
   tree_out->Branch("genParticle_motherPdgId", &genParts.motherPdgId);

   tree_out->Branch("primary_Vertex_x", &PVs.x);
   tree_out->Branch("primary_Vertex_y", &PVs.y);

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
   tree_out->Branch("lowPtElectron_charge", &lowPtElecs.charge);
   tree_out->Branch("lowPtElectron_Lxy_SV", &lowPtElecs.Lxy_SV);
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
   //int type = 99;  // 0:electrons 1:muons
   //bool canFitVertex = false;
   //bool hasValidVertex = false;
   double temp_normalizedChi2;
   double temp_Lxy_SV = 0;
   double normalizedChi2 = 0;
   double vx = 0;                  // x coordinate of dilepton vertex
   double vy = 0;                  // y coordinate of dilepton vertex
   double px = 0;                  // x coordinate of primary vertex
   double py = 0;                  // y coordinate of primary vertex
   
   double Lxy_SV = -99; //Lxy of the secondary vertex


   iEvent.getByToken(genParticleToken, genParticles);
   iEvent.getByToken(photonToken, photons);
   iEvent.getByToken(lowPtElectronToken, lowPtElectrons);
   iEvent.getByToken(electronToken, electrons);
   iEvent.getByToken(primaryVertices, pVtxs);
   
   ebRecHits_ = iEvent.get(reducedBarrelRecHitCollectionToken_);
   eeRecHits_ = iEvent.get(reducedEndcapRecHitCollectionToken_);
   //iEvent.getByToken(reducedBarrelRecHitCollectionToken_, ebRecHits_);
   //iEvent.getByToken(reducedEndcapRecHitCollectionToken_, eeRecHits_);

   const TransientTrackBuilder* theB = &iSetup.getData(ttbESToken_);

   // Clear all variables
   nGP = 0;
   genParts.clear();
   nSC = 0;
   SCs.clear();
   PVs.clear();
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

   // Primary Vertex
  
   // Generated particles
   if (!isData) {
     for (unsigned int iGen=0; iGen<genParticles->size(); iGen++) {
       auto genpart = genParticles->at(iGen);
       if (abs(genpart.pdgId())!=13 && // Muon
           abs(genpart.pdgId())!=11 && // Electron
           abs(genpart.pdgId())!=443) // JPsi
         continue;
       if (!genpart.isLastCopy())
         continue;

       int motherIdx = -1, motherPdgId = 0; // Default value
       if (genpart.numberOfMothers() > 0) {
         motherIdx = genpart.motherRef().index();
         while (genParticles->at(motherIdx).pdgId() == genpart.pdgId()) {
           motherIdx = genParticles->at(motherIdx).motherRef().index();
         }
         motherPdgId = genParticles->at(motherIdx).pdgId();
       }
       genParts.pt.push_back(genpart.pt());
       genParts.eta.push_back(genpart.eta());
       genParts.phi.push_back(genpart.phi());
       genParts.m.push_back(genpart.mass());
       genParts.vx.push_back(genpart.vx());
       genParts.vy.push_back(genpart.vy());
       genParts.vz.push_back(genpart.vz());
       genParts.lxy.push_back(TMath::Hypot(genpart.vx(), genpart.vy()));
       genParts.status.push_back(genpart.status());
       genParts.index.push_back(iGen);
       genParts.pdgId.push_back(genpart.pdgId());
       genParts.motherIndex.push_back(motherIdx);
       genParts.motherPdgId.push_back(motherPdgId);
     }
   }
     
   // Low Pt Electrons
   for (const auto& ele : *lowPtElectrons) {
     std::cout << ele.pt() << std::endl;
     nLowPtElectron++;
     lowPtElecs.pt.push_back(ele.pt());
     lowPtElecs.et.push_back(ele.et());
     lowPtElecs.phi.push_back(ele.phi());
     lowPtElecs.eta.push_back(ele.eta());
     lowPtElecs.charge.push_back(ele.charge());
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

   // Pairing electrons and finding tracks
   for (int i=0; i<nLowPtElectron; i++){
    //Gets first electron
    //std::cout << "Electron" << std::endl;
    //std::cout << i << std::endl;
    Lxy_SV = -99;
    normalizedChi2 = 99;
    //hasValidVertex = false;
    for (int j=i+1; j<nLowPtElectron; j++){
      //std::cout << "Paied with Electron" << std::endl;
      //std::cout << j << std::endl;
      double charge1 = lowPtElectrons->at(i).charge();
      double charge2 = lowPtElectrons->at(j).charge();
      if (charge1*charge2>0){
        //std::cout << "Invalid Pair" << std::endl;
        continue;
      };
        
    
      // Get tracks:
      const reco::Track &tr_A = *(lowPtElectrons->at(i).gsfTrack());
      const reco::Track &tr_B = *(lowPtElectrons->at(j).gsfTrack());
      std::vector<reco::TransientTrack> vec_refitTracks;
      reco::TransientTrack isotransienttrackA = theB->build(tr_A);
      reco::TransientTrack isotransienttrackB = theB->build(tr_B);
      vec_refitTracks.push_back(isotransienttrackA); vec_refitTracks.push_back(isotransienttrackB);
    
      // Fit tracks:
      KalmanVertexFitter thefitterll;
      TransientVertex myVertex = thefitterll.vertex(vec_refitTracks);
      const reco::Vertex secV = myVertex;
      const reco::Vertex PrimV = (*pVtxs)[0];


      if (secV.isValid()) {
      
        //hasValidVertex = true;

        // Define the axis along the direction of the distance is defined:
        GlobalVector axis(0,0,0);
        axis = GlobalVector(secV.x(),secV.y(),secV.z());
        
        // Define the errors and points of the CMS centre point and the beam spot:
        /*//math::Error<3>::type e0; // dummy
        //math::Error<3>::type cov = bs.covariance3D();
        //math::XYZPoint pbs(bs.x0(), bs.y0(), bs.z0());
        math::XYZPoint p0(0.0, 0.0, 0.0);

        // Fake vertices for the CMS centre and beam spot:
        const reco::Vertex v0(p0, e0);
        const reco::Vertex vbs(pbs, cov);

        // Measurements:
        Measurement1D vMeas_PV = reco::SecondaryVertex::computeDist2d(pv,secV,axis,true);
        Measurement1D vMeas_0 = reco::SecondaryVertex::computeDist2d(v0,secV,axis,false);
        Measurement1D vMeas_BS = reco::SecondaryVertex::computeDist2d(vbs,secV,axis,true);
        

        // Distance values:
        Lxy_0 = vMeas_0.value();
        Ixy_0 = vMeas_0.significance();
        Lxy_PV = vMeas_PV.value();
        Ixy_PV = vMeas_PV.significance();
        Lxy_BS = vMeas_BS.value();
        Ixy_BS = vMeas_BS.significance();

        */// Vertex position and fit details:
        temp_normalizedChi2 = myVertex.normalisedChiSquared();         
        vx = secV.x();
        vy = secV.y();
        
        //std::cout << PrimV.x() << std::endl;
        //std::cout << PVs.x << std::endl;
        temp_Lxy_SV = std::sqrt(std::pow(vx-PrimV.x(), 2) + std::pow(vy-PrimV.y(), 2));
        if (temp_normalizedChi2 < normalizedChi2){
          Lxy_SV = temp_Lxy_SV;
          normalizedChi2 = temp_normalizedChi2;
          px = PrimV.x();
          py = PrimV.y();
        };
        /*

        // Kinematics: 
        leadingPt = (tr_A.pt()>tr_B.pt())? tr_A.pt(): tr_B.pt();
        subleadingPt = (tr_A.pt()<tr_B.pt())? tr_A.pt(): tr_B.pt();
        etaA = tr_A.eta();
        etaB = tr_B.eta();   

        // Vector angles:  
        TVector3 vec3A(tr_A.px(), tr_A.py(), tr_A.pz());
        TVector3 vec3B(tr_B.px(), tr_B.py(), tr_B.pz());
        TVector3 divec3 = vec3A + vec3B;
        TVector3 vtxvec3(secV.x() - pv.x(), secV.y() - pv.y(), secV.z() - pv.z());
        cosAlpha = TMath::Cos(vec3A.Angle(vec3B));
        dPhi = divec3.DeltaPhi(vtxvec3);
        dR = vec3A.DeltaR(vec3B);
        */
      };
    };
    lowPtElecs.Lxy_SV.push_back(Lxy_SV);
    PVs.x.push_back(px);
    PVs.y.push_back(py);
   };

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
       if (theSeedHit != recHits->end()){
         SCs.seedTime.push_back((*theSeedHit).time());
       } else{
         SCs.seedTime.push_back(-99.);
       }
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
