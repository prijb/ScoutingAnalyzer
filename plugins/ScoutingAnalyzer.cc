// -*- C++ -*-
//
// Package:    Scouting2022Analyzer/ScoutingAnalyzer
// Class:      ScoutingAnalyzer
//
/**\class ScoutingAnalyzer ScoutingAnalyzer.cc Scouting2022Analyzer/ScoutingAnalyzer/plugins/ScoutingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Prijith Babu Pradeep
//         Created:  Mon, 19 Jun 2023 00:18:32 GMT
//
//

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

// dataformats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/fillCovariance.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/Common/interface/AssociationMap.h"

// dataformats (PAT stuff?)
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// dataformats (scouting specific?)
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingHitPatternPOD.h"

// dataformats (trigger)
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// what?? (gen??)
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// tracking tools 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// CMSSW stuff
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace std;

using reco::TrackCollection;

class ScoutingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingAnalyzer(const edm::ParameterSet&);
  ~ScoutingAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // Stuff that involves inheritance from WatchRuns and WatchLuminosityBlocks
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void clearVars();

  // ----------member data ---------------------------
  // The darn tokens (Will be initialised later)
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> PVToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> SVToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> pfToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron>> electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonsToken;

  const edm::EDGetTokenT<double> rhoToken;
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;  //HLT results token

  // Mapping trigger paths
  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  // L1 stuff
  bool doL1;
  triggerExpression::Data triggerCache_;

  // Some malarkey involving triggers
  unsigned char trig;
  edm::InputTag algInputTag_;
  edm::InputTag extInputTag_;
  edm::EDGetToken algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string> l1Seeds_;
  std::vector<bool> l1Result_;


  // This thing is also an ntuplizer, so now initialise the variables to fill a tree
  // Primary vertex details
  UInt_t n_PV;
  vector<Float16_t> PV_x;
  vector<Float16_t> PV_y;
  vector<Float16_t> PV_z;
  vector<Float16_t> PV_xError;
  vector<Float16_t> PV_yError;
  vector<Float16_t> PV_zError;
  vector<Int_t> PV_trksize;
  vector<Float16_t> PV_chi2;
  vector<Int_t> PV_ndof;
  vector<Bool_t> PV_isvalidvtx;

  //Secondary vertex details
  UInt_t n_SV;
  vector<Float16_t> SV_x;
  vector<Float16_t> SV_y;
  vector<Float16_t> SV_z;
  vector<Float16_t> SV_dxy;
  vector<Float16_t> SV_dxySig;
  vector<Float16_t> SV_xError;
  vector<Float16_t> SV_yError;
  vector<Float16_t> SV_zError;
  vector<Int_t> SV_trksize;
  vector<Float16_t> SV_chi2;
  vector<Int_t> SV_ndof;
  vector<Bool_t> SV_isvalidvtx;

  vector<Int_t> SV_nMuon;
  vector<Int_t> SV_mu1idx;
  vector<Int_t> SV_mu2idx;
  vector<Float16_t> SV_mass;

  //Secondary vertex (overlap) details
  UInt_t n_SVOverlap;
  vector<Float16_t> SVOverlap_x;
  vector<Float16_t> SVOverlap_y;
  vector<Float16_t> SVOverlap_z;
  vector<Float16_t> SVOverlap_dxy;
  vector<Int_t> SVOverlap_sv1idx;
  vector<Int_t> SVOverlap_sv2idx;
  vector<Float16_t> SVOverlap_mass;

  /*
  //PF candidates
  UInt_t n_pf;
  vector<Float16_t> PFParticle_pt;
  vector<Float16_t> PFParticle_eta;
  vector<Float16_t> PFParticle_phi;
  vector<Int_t> PFParticle_pdgId;
  vector<Int_t> PFParticle_vertex;
  vector<Float16_t> PFParticle_normchi2;
  vector<Float16_t> PFParticle_dz;
  vector<Float16_t> PFParticle_dxy;
  vector<Float16_t> PFParticle_dzsig;
  vector<Float16_t> PFParticle_dxysig;
  vector<UInt_t>	PFParticle_lostInnerHits;
  vector<UInt_t>	PFParticle_quality;
  vector<Float16_t> PFParticle_trkpt;
  vector<Float16_t> PFParticle_trketa;
  vector<Float16_t> PFParticle_trkphi;
  vector<Bool_t> PFParticle_relativetrkvars;
  */

  //Electron colletion details (what is max_ele for??)
  
  const static int 	max_ele = 1000;
  UInt_t n_ele;
  vector<Float16_t> Electron_pt;
  vector<Float16_t> Electron_eta;
  vector<Float16_t> Electron_phi;
  vector<Float16_t> Electron_m;
  vector<Int_t> Electron_charge;
  vector<Float16_t> Electron_detain;
  vector<Float16_t> Electron_dphiin;
  vector<Float16_t> Electron_sigmaietaieta;
  vector<Float16_t> Electron_hoe;
  vector<Float16_t> Electron_ooemoop;
  vector<Int_t>	Electron_missinghits;
  vector<Float16_t> Electron_ecaliso;
  vector<Float16_t> Electron_hcaliso;
  vector<Float16_t> Electron_tkiso;
  vector<Float16_t> Electron_r9;
  vector<Float16_t> Electron_smin;
  vector<Float16_t> Electron_smaj;
  vector<UInt_t> Electron_seedid;
  vector<Bool_t> Electron_rechitzerosuppression;
  

  //Muon collection details (added vertex id)
  const static int 	max_mu = 1000;
  UInt_t n_mu;
  vector<Float16_t> Muon_pt;
  vector<Float16_t> Muon_eta;
  vector<Float16_t> Muon_phi;
  vector<Float16_t> Muon_m;
  vector<Float16_t> Muon_ecaliso;
  vector<Float16_t> Muon_hcaliso;
  vector<Float16_t> Muon_trkiso;
  vector<Float16_t> Muon_chi2;
  vector<Float16_t> Muon_ndof;
  vector<Float16_t> Muon_charge;
  vector<Float16_t> Muon_dxy;
  vector<Float16_t> Muon_dz;
  vector<Float16_t> Muon_dxyerror;
  vector<Float16_t> Muon_dzerror;
  vector<Float16_t> Muon_nvalidmuon_hits;
  vector<Float16_t> Muon_nvalidpixelhits;
  
  vector<Float16_t> Muon_nmatchedstations;
  vector<Float16_t> Muon_type;
  vector<Float16_t> Muon_nvalidstriphits;
  vector<Float16_t> Muon_trkqoverp;
  vector<Float16_t> Muon_trklambda;
  vector<Float16_t> Muon_trkpt;
  vector<Float16_t> Muon_trkphi;
  vector<Float16_t> Muon_trketa;
  vector<Float16_t> Muon_trkqoverperror;
  vector<Float16_t> Muon_trklambdaerror;
  vector<Float16_t> Muon_trkpterror;
  vector<Float16_t> Muon_trkphierror;
  vector<Float16_t> Muon_trketaerror;
  vector<Float16_t> Muon_trkdszerror;
  vector<Float16_t> Muon_trkdsz;

  //Vertex id
  vector<std::vector<int>> Muon_vtxIndx;

  //Hit pattern info
  vector<uint8_t> Muon_hitCount;
  vector<uint8_t> Muon_beginTrackHits;
  vector<uint8_t> Muon_endTrackHits;
  vector<uint8_t> Muon_beginInner;
  vector<uint8_t> Muon_endInner;
  vector<uint8_t> Muon_beginOuter;
  vector<uint8_t> Muon_endOuter;


  //Rho info
  UInt_t n_rhoval;
  vector<Float16_t> rho;

  //The actual tree declared. We're gonna make it in the constructor dw
  TTree* tree;

  //Run and lumisection
  Int_t run;
  Int_t event;
  Int_t lumSec;




};

//
// constructors and destructor (a lotta initialisations of all the tokens)
// Initialisation is of the form
// tokenName(consumes<type>(iConfig.getParameter<edm::InputTag>(string you pass in python config)))
ScoutingAnalyzer::ScoutingAnalyzer(const edm::ParameterSet& iConfig):
  PVToken(consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("primaryVtx"))),
  SVToken(consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("secondaryVtx"))),
  //pfToken(consumes<std::vector<Run3ScoutingParticle>>(iConfig.getParameter<edm::InputTag>("pfparticles"))),
  electronsToken(consumes<std::vector<Run3ScoutingElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonsToken(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  //Get trigger related tokens
  triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerresults")),
  triggerResultsToken(consumes<edm::TriggerResults>(triggerResultsTag)),
  //What sorcery is this?
  doL1(iConfig.existsAs<bool>("doL1")?iConfig.getParameter<bool>("doL1"):false)
  {

  //now do what ever initialization is needed
  // What does usesResource mean??
  usesResource("TFileService");

  // L1 stuff?
  if (doL1){
    algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
    extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
    algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
    l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(
    iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }

  // Access the TFileService (help...)
  edm::Service<TFileService> fs;

  // Make the tree using the file service (huh....)
  tree = fs->make<TTree>("tree", "tree");

  // These three branches contain the event number, and what run and lumisection they belong to?
  tree->Branch("lumSec", &lumSec, "lumSec/i" );
  tree->Branch("run", &run, "run/i" );
  tree->Branch("event", &event, "event/i" );

  // The interesting stuff

  // Triggers
  tree->Branch("trig", &trig, "trig/b");
  tree->Branch("l1Result", "std::vector<bool>" ,&l1Result_, 32000, 0);

  // Primary vertex info (why does only the first branch have a /i?)
  tree->Branch("n_PV", &n_PV, "n_PV/i");
  tree->Branch("PV_x", &PV_x);
  tree->Branch("PV_y", &PV_y);
  tree->Branch("PV_z", &PV_z);
  tree->Branch("PV_xError", &PV_xError);
  tree->Branch("PV_yError", &PV_yError);
  tree->Branch("PV_zError", &PV_zError);
  tree->Branch("PV_trksize", &PV_trksize);
  tree->Branch("PV_chi2", &PV_chi2);
  tree->Branch("PV_ndof", &PV_ndof);
  tree->Branch("PV_isvalidvtx", &PV_isvalidvtx);

  // Secondary vertex info 
  tree->Branch("n_SV", &n_SV, "n_SV/i");
  tree->Branch("SV_x", &SV_x);
  tree->Branch("SV_y", &SV_y);
  tree->Branch("SV_z", &SV_z);
  tree->Branch("SV_xError", &SV_xError);
  tree->Branch("SV_yError", &SV_yError);
  tree->Branch("SV_zError", &SV_zError);
  tree->Branch("SV_dxy", &SV_dxy);
  tree->Branch("SV_dxySig", &SV_dxySig);
  tree->Branch("SV_trksize", &SV_trksize);
  tree->Branch("SV_chi2", &SV_chi2);
  tree->Branch("SV_ndof", &SV_ndof);
  tree->Branch("SV_isvalidvtx", &SV_isvalidvtx);

  tree->Branch("SV_nMuon", &SV_nMuon);
  tree->Branch("SV_mu1idx", &SV_mu1idx);
  tree->Branch("SV_mu2idx", &SV_mu2idx);
  tree->Branch("SV_mass", &SV_mass);

  // Overlap secondary vertex info
  tree->Branch("n_SVOverlap", &n_SVOverlap, "n_SVOverlap/i");
  tree->Branch("SVOverlap_x", &SVOverlap_x);
  tree->Branch("SVOverlap_y", &SVOverlap_y);
  tree->Branch("SVOverlap_z", &SVOverlap_z);
  tree->Branch("SVOverlap_dxy", &SVOverlap_dxy);
  tree->Branch("SVOverlap_sv1idx", &SVOverlap_sv1idx);
  tree->Branch("SVOverlap_sv2idx", &SVOverlap_sv2idx);
  tree->Branch("SVOverlap_mass", &SVOverlap_mass);

  /*
  // Pf particles
  tree->Branch("n_pf", &n_pf, "n_pf/i");
  tree->Branch("PFParticle_pt", &PFParticle_pt);
  tree->Branch("PFParticle_eta", &PFParticle_eta);
  tree->Branch("PFParticle_phi", &PFParticle_phi);
  tree->Branch("PFParticle_pdgId", &PFParticle_pdgId);
  tree->Branch("PFParticle_vertex", &PFParticle_vertex);
  tree->Branch("PFParticle_normchi2", &PFParticle_normchi2);
  tree->Branch("PFParticle_dz", &PFParticle_dz);
  tree->Branch("PFParticle_dxy", &PFParticle_dxy);
  tree->Branch("PFParticle_dzsig", &PFParticle_dzsig);
  tree->Branch("PFParticle_dxysig", &PFParticle_dxysig);
  tree->Branch("PFParticle_lostInnerHits", &PFParticle_lostInnerHits);
  tree->Branch("PFParticle_quality", &PFParticle_quality);
  tree->Branch("PFParticle_trkpt", &PFParticle_trkpt);
  tree->Branch("PFParticle_trketa", &PFParticle_trketa);
  tree->Branch("PFParticle_trkphi", &PFParticle_trkphi);
  tree->Branch("PFParticle_relativetrkvars", &PFParticle_relativetrkvars);
  */



  // Electrons
  // Note: Some shenanigans involving electron variables
  
  tree->Branch("n_ele", &n_ele, "n_ele/i");
  tree->Branch("Electron_pt", &Electron_pt);
  tree->Branch("Electron_eta", &Electron_eta);
  tree->Branch("Electron_phi", &Electron_phi);
  tree->Branch("Electron_m", &Electron_m);
  tree->Branch("Electron_charge", &Electron_charge);
  tree->Branch("Electron_detain", &Electron_detain);
  tree->Branch("Electron_dphiin", &Electron_dphiin);
  tree->Branch("Electron_sigmaietaieta", &Electron_sigmaietaieta);
  tree->Branch("Electron_hoe", &Electron_hoe);
  tree->Branch("Electron_ooemoop", &Electron_ooemoop);
  tree->Branch("Electron_missinghits", &Electron_missinghits);
  tree->Branch("Electron_ecaliso", &Electron_ecaliso);
  tree->Branch("Electron_hcaliso", &Electron_hcaliso);
  tree->Branch("Electron_tkiso", &Electron_tkiso);
  tree->Branch("Electron_r9", &Electron_r9);
  tree->Branch("Electron_smin", &Electron_smaj);
  tree->Branch("Electron_smaj", &Electron_smin);
  tree->Branch("Electron_seedid", &Electron_seedid);
  tree->Branch("Electron_rechitzerosuppression", &Electron_rechitzerosuppression);
  
  // Muons
  tree->Branch("n_mu",&n_mu,"n_mu/i");
  tree->Branch("Muon_pt", &Muon_pt);
  tree->Branch("Muon_eta", &Muon_eta);
  tree->Branch("Muon_phi", &Muon_phi);
  tree->Branch("Muon_m", &Muon_m);
  tree->Branch("Muon_ecaliso", &Muon_ecaliso);
  tree->Branch("Muon_hcaliso", &Muon_hcaliso);
  tree->Branch("Muon_trkiso", &Muon_trkiso);
  tree->Branch("Muon_chi2", &Muon_chi2);
  tree->Branch("Muon_ndof", &Muon_ndof);
  tree->Branch("Muon_charge", &Muon_charge);
  tree->Branch("Muon_dxy", &Muon_dxy);
  tree->Branch("Muon_dz", &Muon_dz);
  tree->Branch("Muon_dxyerror", &Muon_dxyerror);
  tree->Branch("Muon_dzerror", &Muon_dzerror);
  tree->Branch("Muon_nvalidmuon_hits", &Muon_nvalidmuon_hits);
  tree->Branch("Muon_validpixelhits", &Muon_nvalidpixelhits );
  
  tree->Branch("Muon_nmatchedstations", &Muon_nmatchedstations);
  tree->Branch("Muon_type", &Muon_type);
  tree->Branch("Muon_nvalidstriphits", &Muon_nvalidstriphits);
  tree->Branch("Muon_trkqoverp", &Muon_trkqoverp);
  tree->Branch("Muon_trklambda", &Muon_trklambda);
  tree->Branch("Muon_trkpt", &Muon_trkpt);
  tree->Branch("Muon_trkphi", &Muon_trkphi);
  tree->Branch("Muon_trketa", &Muon_trketa);
  tree->Branch("Muon_trkqoverperror", &Muon_trkqoverperror);
  tree->Branch("Muon_trklambdaerror", &Muon_trklambdaerror);
  tree->Branch("Muon_trkpterror", &Muon_trkpterror);
  tree->Branch("Muon_trkphierror", &Muon_trkphierror);
  tree->Branch("Muon_trketaerror", &Muon_trketaerror);
  tree->Branch("Muon_trkdzerror", &Muon_trkdszerror);
  tree->Branch("Muon_trkdz", &Muon_trkdsz);

  tree->Branch("Muon_vtxIndx", &Muon_vtxIndx);

  tree->Branch("Muon_hitCount", &Muon_hitCount);
  tree->Branch("Muon_beginTrackHits", &Muon_beginTrackHits);
  tree->Branch("Muon_endTrackHits", &Muon_endTrackHits);
  tree->Branch("Muon_beginInner", &Muon_beginInner);
  tree->Branch("Muon_endInner", &Muon_endInner);
  tree->Branch("Muon_beginOuter", &Muon_beginOuter);
  tree->Branch("Muon_endOuter", &Muon_endOuter);

  // Rho
  tree->Branch("n_rhoval", &n_rhoval, "n_rhoval/i");
  tree->Branch("rho", &rho);

}

// Destructor!
ScoutingAnalyzer::~ScoutingAnalyzer() {
}

//
// member functions
//

// Fill stuff

// ------------ method called for each event  ------------
void ScoutingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  using namespace std;
  using namespace reco;

  //Token's already initalized, so now we use this token to fill a newly initialized handle
  //This is the FWLite stuff you're more familiar with but FWLite gets by label rather than token

  //Get trigger handles
  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken(triggerResultsToken, triggerResultsHandle);
  bool triggerResultsValid = triggerResultsHandle.isValid();

  Handle<vector<Run3ScoutingVertex>> PVHandle;
  iEvent.getByToken(PVToken, PVHandle);
  bool PVValid = PVHandle.isValid();

  Handle<vector<Run3ScoutingVertex>> SVHandle;
  iEvent.getByToken(SVToken, SVHandle);
  bool SVValid = SVHandle.isValid();

  //Handle<vector<Run3ScoutingParticle>> pfHandle;
  //iEvent.getByToken(pfToken, pfHandle);
  //bool pfValid = pfHandle.isValid();

  
  Handle<vector<Run3ScoutingElectron>> electronsHandle;
  iEvent.getByToken(electronsToken, electronsHandle);
  bool electronValid = electronsHandle.isValid();
  

  Handle<vector<Run3ScoutingMuon>> muonsHandle;
  iEvent.getByToken(muonsToken, muonsHandle);
  bool muonValid = muonsHandle.isValid();

  Handle<double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  bool rhoValid = rhoHandle.isValid();

  //No clue what this is but it directly feeds into the branch without handles
  run = iEvent.eventAuxiliary().run();
  event = iEvent.eventAuxiliary().event();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();
  
  //Now fill collections!

  //HLT. Why is trig summed over all events this way without resetting?
  trig = 0;
  //Check how triggers are fired?
  if(triggerResultsValid){
    for (size_t i = 0; i < triggerPathsVector.size(); i++){
      //??
      if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
      //If the 
      if (i == 0  && triggerResultsHandle->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   1; // DST_PixelTracking
      if (i == 1  && triggerResultsHandle->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   2; // DST_HLTMu
    }
  }
  
  bool mu_or_bit = 0;

  //L1 
  if (doL1) {
   //Get l1GtUtils from the event, event setup and the token
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);

  
  //Runs over all the L1 algorithms and checks if event passes or not
  /*
    for( int r = 0; r<280; r++){
	    string name ("empty");
	    bool algoName_ = false;
	    algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
	    cout << "getAlgNameFromBit = " << algoName_  << endl;
	    cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
	  }
  */

    //I think this fills the data with the L1 decisions?
    for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
      bool l1htbit = 0;	
      l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
      //cout<<string(l1Seeds_[iseed])<<"  "<<l1htbit<<endl;
      //cout << "L1 seed name: " << string(l1Seeds_[iseed]) << endl;
      //cout << "L1 decision: " <<  l1htbit << endl;
      //cout << "" << endl;
      l1Result_.push_back( l1htbit );
      mu_or_bit += l1htbit;
      }
  }


  //Fill PV where you count and add that to event as well
  //Check if handle is valid
  //* dereferences handle and gets collection. Then pass each element by reference?
  n_PV = 0;
  if(PVValid){
    for (auto &PV: *PVHandle){
      //Use the pointer to PV 
      auto *PV_iter = &PV;
      PV_x.push_back(PV_iter->x());
      PV_y.push_back(PV_iter->y());
      PV_z.push_back(PV_iter->z());
      PV_xError.push_back(PV_iter->xError());
      PV_yError.push_back(PV_iter->yError());
      PV_zError.push_back(PV_iter->zError());
      PV_trksize.push_back(PV_iter->tracksSize());
      PV_chi2.push_back(PV_iter->chi2());
      PV_ndof.push_back(PV_iter->ndof());
      PV_isvalidvtx.push_back(PV_iter->isValidVtx());
      n_PV++;
    }
  }

  //Fill SV
  n_SV = 0;

  if(SVValid){
    for (auto &SV: *SVHandle){
      auto *SV_iter = &SV;
      SV_x.push_back(SV_iter->x());
      SV_y.push_back(SV_iter->y());
      SV_z.push_back(SV_iter->z());
      //Calculate dxy
      Float16_t dx = SV_iter->x() - PV_x.at(0);
      Float16_t dy = SV_iter->y() - PV_y.at(0);
      Float16_t dxy = TMath::Sqrt(dx*dx + dy*dy);
      SV_dxy.push_back(dxy);
      SV_xError.push_back(SV_iter->xError());
      SV_yError.push_back(SV_iter->yError());
      SV_zError.push_back(SV_iter->zError());
      //Calculate dxysig
      Float16_t dxerr = TMath::Sqrt(std::pow(SV_iter->xError(), 2) + std::pow(PV_xError.at(0), 2));
      Float16_t dyerr = TMath::Sqrt(std::pow(SV_iter->yError(), 2) + std::pow(PV_yError.at(0), 2));
      Float16_t dxyerr = TMath::Sqrt(std::pow(dx*dxerr, 2) + std::pow(dy*dyerr, 2))/dxy;
      Float16_t dxysig = dxy/dxyerr;
      SV_dxySig.push_back(dxysig);
      SV_trksize.push_back(SV_iter->tracksSize());
      SV_chi2.push_back(SV_iter->chi2());
      SV_ndof.push_back(SV_iter->ndof());
      SV_isvalidvtx.push_back(SV_iter->isValidVtx());

      //Fill with mass
      TLorentzVector SV_Sum;
      UInt_t n_SV_mu = 0;
      //Muon index vector
      std::vector<int> muidx_vec(2,-1);


      if (muonValid){
        int muidx = 0;
        for (auto &muo: *muonsHandle){
          auto *muons_iter = &muo;
          const auto& Muon_vtxIndx_inter = muons_iter->vtxIndx();

          if(std::find(Muon_vtxIndx_inter.begin(), Muon_vtxIndx_inter.end(), n_SV) != Muon_vtxIndx_inter.end()){
            //Assign the muon index (skip if already two muons have been associated)
            if(n_SV_mu < 2) muidx_vec[n_SV_mu] = muidx;
            //Sum up to compute mass
            TLorentzVector muon_fourvec;
            muon_fourvec.SetPtEtaPhiM(muons_iter->pt(), muons_iter->eta(), muons_iter->phi(), muons_iter->m());
            SV_Sum += muon_fourvec;
            n_SV_mu += 1;
          }
          muidx += 1;
        }
      }

      SV_nMuon.push_back(n_SV_mu);
      SV_mu1idx.push_back(muidx_vec[0]);
      SV_mu2idx.push_back(muidx_vec[1]);
      SV_mass.push_back(SV_Sum.M());
      n_SV++;
    }
  }

  //Make OverlapSVs
  float highest_sv_uncert_x, highest_sv_uncert_y, highest_sv_uncert_z;
  float sv_delx, sv_dely, sv_delz;
  int svoverlap_mu1idx, svoverlap_mu2idx, svoverlap_mu3idx, svoverlap_mu4idx;

  n_SVOverlap = 0;

  for(UInt_t sv_i=0; sv_i<n_SV; sv_i++){
    for(UInt_t sv_j=(sv_i+1); sv_j<n_SV; sv_j++){
      highest_sv_uncert_x = (SV_xError.at(sv_i) > SV_xError.at(sv_j)) ? SV_xError.at(sv_i) : SV_xError.at(sv_j); 
      highest_sv_uncert_y = (SV_yError.at(sv_i) > SV_yError.at(sv_j)) ? SV_yError.at(sv_i) : SV_yError.at(sv_j); 
      highest_sv_uncert_z = (SV_zError.at(sv_i) > SV_zError.at(sv_j)) ? SV_zError.at(sv_i) : SV_zError.at(sv_j); 

      sv_delx = TMath::Abs(SV_x.at(sv_i) - SV_x.at(sv_j));
      sv_dely = TMath::Abs(SV_y.at(sv_i) - SV_y.at(sv_j));
      sv_delz = TMath::Abs(SV_z.at(sv_i) - SV_z.at(sv_j));

      //If the positions are closer than max uncertainty, they overlap
      if(((sv_delx<highest_sv_uncert_x)&&(sv_dely<highest_sv_uncert_y))&&(sv_delz<highest_sv_uncert_z)){
        SVOverlap_x.push_back(0.5*(SV_x.at(sv_i) + SV_x.at(sv_j)));
        SVOverlap_y.push_back(0.5*(SV_y.at(sv_i) + SV_y.at(sv_j)));
        SVOverlap_z.push_back(0.5*(SV_z.at(sv_i) + SV_z.at(sv_j)));
        SVOverlap_dxy.push_back(0.5*(SV_dxy.at(sv_i) + SV_dxy.at(sv_j)));
        SVOverlap_sv1idx.push_back(sv_i);
        SVOverlap_sv2idx.push_back(sv_j);

        svoverlap_mu1idx = SV_mu1idx.at(sv_i);
        svoverlap_mu2idx = SV_mu2idx.at(sv_i);
        svoverlap_mu3idx = SV_mu1idx.at(sv_j);
        svoverlap_mu4idx = SV_mu2idx.at(sv_j);

        //Mass calculation if all four ids are valid
        if((svoverlap_mu1idx==-1)||(svoverlap_mu2idx==-1)||(svoverlap_mu3idx==-1)||(svoverlap_mu4idx==-1)) SVOverlap_mass.push_back(-1.);

        else{
          TLorentzVector SVOverlap_Sum;
          auto muon1 = muonsHandle->at(svoverlap_mu1idx);
          auto muon2 = muonsHandle->at(svoverlap_mu2idx);
          auto muon3 = muonsHandle->at(svoverlap_mu3idx);
          auto muon4 = muonsHandle->at(svoverlap_mu4idx);

          TLorentzVector muon_fourvec1, muon_fourvec2, muon_fourvec3, muon_fourvec4;
          
          muon_fourvec1.SetPtEtaPhiM(muon1.pt(), muon1.eta(), muon1.phi(), muon1.m());
          muon_fourvec2.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), muon2.m());
          muon_fourvec3.SetPtEtaPhiM(muon3.pt(), muon3.eta(), muon3.phi(), muon3.m());
          muon_fourvec4.SetPtEtaPhiM(muon4.pt(), muon4.eta(), muon4.phi(), muon4.m());

          SVOverlap_Sum = muon_fourvec1 + muon_fourvec2 + muon_fourvec3 + muon_fourvec4;
          SVOverlap_mass.push_back(SVOverlap_Sum.M());
        }

        n_SVOverlap++;
      }
    } 
  }

  /*
  n_pf = 0;
  if(pfValid) {
    for (auto &pfcand : *pfHandle) {
      auto *pfcand_iter = &pfcand;
      PFParticle_pt.push_back(pfcand_iter->pt());
      PFParticle_eta.push_back(pfcand_iter->eta());
      PFParticle_phi.push_back(pfcand_iter->phi());
      PFParticle_pdgId.push_back(pfcand_iter->pdgId());
      PFParticle_vertex.push_back(pfcand_iter->vertex());
      PFParticle_normchi2.push_back(pfcand_iter->normchi2());
      PFParticle_dz.push_back(pfcand_iter->dz());
      PFParticle_dxy.push_back(pfcand_iter->dxy());
      PFParticle_dzsig.push_back(pfcand_iter->dzsig());
      PFParticle_dxysig.push_back(pfcand_iter->dxysig());
      PFParticle_lostInnerHits.push_back(pfcand_iter->lostInnerHits());
      PFParticle_quality.push_back(pfcand_iter->quality());
      PFParticle_trkpt.push_back(pfcand_iter->trk_pt());
      PFParticle_trketa.push_back(pfcand_iter->trk_eta());
      PFParticle_trkphi.push_back(pfcand_iter->trk_phi());
      PFParticle_relativetrkvars.push_back(pfcand_iter->relative_trk_vars());
      n_pf++;
    } 
  } 
  */

  
  //Fill muons
  n_mu = 0;
  if(muonValid) {
    for (auto &muo : *muonsHandle) {
      auto *muons_iter = &muo;
      Muon_pt.push_back(muons_iter->pt());
      Muon_eta.push_back(muons_iter->eta());
      Muon_phi.push_back(muons_iter->phi());	
      Muon_m.push_back(muons_iter->m());
      Muon_ecaliso.push_back(muons_iter->ecalIso());
      Muon_hcaliso.push_back(muons_iter->hcalIso());
      Muon_trkiso.push_back(muons_iter->trackIso());
      Muon_chi2.push_back(muons_iter->normalizedChi2());
      Muon_ndof.push_back(muons_iter->trk_ndof());
      Muon_charge.push_back(muons_iter->charge());
      Muon_dxy.push_back(muons_iter->trk_dxy());
      Muon_dxyerror.push_back(muons_iter->trk_dxyError());
      Muon_dz.push_back(muons_iter->trk_dz());
      Muon_dzerror.push_back(muons_iter->trk_dzError());
      Muon_nvalidmuon_hits.push_back(muons_iter->nValidRecoMuonHits());
      Muon_nvalidpixelhits.push_back(muons_iter->nValidPixelHits());
      Muon_nvalidstriphits.push_back(muons_iter->nValidStripHits());
      Muon_nmatchedstations.push_back(muons_iter->nRecoMuonMatchedStations());

      Muon_type.push_back(muons_iter->type());
      Muon_trkqoverp.push_back(muons_iter->trk_qoverp());
      Muon_trklambda.push_back(muons_iter->trk_lambda());
      Muon_trkpt.push_back(muons_iter->trk_pt());
      Muon_trkphi.push_back(muons_iter->trk_phi());
      Muon_trketa.push_back(muons_iter->trk_eta());
      Muon_trkqoverperror.push_back(muons_iter->trk_dxyError());
      Muon_trklambdaerror.push_back(muons_iter->trk_dzError());
      Muon_trkpterror.push_back(muons_iter->trk_qoverpError());
      Muon_trkphierror.push_back(muons_iter->trk_lambdaError());
      Muon_trketaerror.push_back(muons_iter->trk_phiError());
      Muon_trkdsz.push_back(muons_iter->trk_dsz());
      Muon_trkdszerror.push_back(muons_iter->trk_dszError());

      Muon_vtxIndx.push_back(muons_iter->vtxIndx());

      Muon_hitCount.push_back((muons_iter->trk_hitPattern()).hitCount);
      Muon_beginTrackHits.push_back((muons_iter->trk_hitPattern()).beginTrackHits);
      Muon_endTrackHits.push_back((muons_iter->trk_hitPattern()).endTrackHits);
      Muon_beginInner.push_back((muons_iter->trk_hitPattern()).beginInner);
      Muon_endInner.push_back((muons_iter->trk_hitPattern()).endInner);
      Muon_beginOuter.push_back((muons_iter->trk_hitPattern()).beginOuter);
      Muon_endOuter.push_back((muons_iter->trk_hitPattern()).endOuter);

      
      n_mu++;
    } 
  } 

  //Fill rho
  n_rhoval = 0;
  if(rhoValid) {
    rho.push_back(*rhoHandle);
    n_rhoval++;
  }

  //To study events that mass muon OR
  /*
  if(mu_or_bit == true){
  cout << "Event passes Muon OR" << endl;
  cout << "Number of muons in event: " << n_mu << endl;
  }
  */
  //Finally fill the tree! (Clear vars does the cleaning)
  tree->Fill();
  clearVars();

}

// ------------ clear variables

void ScoutingAnalyzer::clearVars() {

  l1Result_.clear();

  PV_x.clear();
  PV_y.clear();
  PV_z.clear();
  PV_xError.clear();
  PV_yError.clear();
  PV_zError.clear();
  PV_trksize.clear();
  PV_chi2.clear();
  PV_ndof.clear();
  PV_isvalidvtx.clear();
  
  SV_x.clear();
  SV_y.clear();
  SV_z.clear();
  SV_xError.clear();
  SV_yError.clear();
  SV_zError.clear();
  SV_dxy.clear();
  SV_dxySig.clear();
  SV_trksize.clear();
  SV_chi2.clear();
  SV_ndof.clear();
  SV_isvalidvtx.clear();
  SV_nMuon.clear();
  SV_mu1idx.clear();
  SV_mu2idx.clear();
  SV_mass.clear();


  SVOverlap_x.clear();
  SVOverlap_y.clear();
  SVOverlap_z.clear();
  SVOverlap_dxy.clear();
  SVOverlap_sv1idx.clear();
  SVOverlap_sv2idx.clear();
  SVOverlap_mass.clear();


  

  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_m.clear();
  Muon_ecaliso.clear();
  Muon_hcaliso.clear();
  Muon_trkiso.clear();
  Muon_chi2.clear();
  Muon_ndof.clear();
  Muon_charge.clear();
  Muon_dxy.clear();
  Muon_dxyerror.clear();
  Muon_dz.clear();
  Muon_dzerror.clear();
  Muon_nvalidmuon_hits.clear();
  Muon_nvalidpixelhits.clear();
  Muon_nmatchedstations.clear();
  Muon_type.clear();
  Muon_nvalidstriphits.clear();
  Muon_trkqoverp.clear();
  Muon_trklambda.clear();
  Muon_trkpt.clear();
  Muon_trkphi.clear();
  Muon_trketa.clear();
  Muon_trkqoverperror.clear();
  Muon_trklambdaerror.clear();
  Muon_trkpterror.clear();
  Muon_trkphierror.clear();
  Muon_trketaerror.clear();
  Muon_trkdszerror.clear();
  Muon_trkdsz.clear();

  Muon_vtxIndx.clear();

  Muon_hitCount.clear();
  Muon_beginTrackHits.clear();
  Muon_endTrackHits.clear();
  Muon_beginInner.clear();
  Muon_endInner.clear();
  Muon_beginOuter.clear();
  Muon_endOuter.clear();

  rho.clear();

}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingAnalyzer::endJob() {
  // please remove this method if not needed
}

void ScoutingAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  // HLT paths

  //For 2022
  triggerPathsVector.push_back("DST_Run3_PFScoutingPixelTracking_v*");
  triggerPathsVector.push_back("DST_HLTMuon_Run3_PFScoutingPixelTracking_v*");
 
  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
      std::string pathName = hltConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }
}

void ScoutingAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingAnalyzer);