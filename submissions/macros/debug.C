#include <iostream>
#include <map>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include <chrono>
#include <cmath>
#include "Math/Vector4D.h"

#include "MuPOG.C"

using namespace std;

#include "correction.h"
using correction::CorrectionSet;

#include "Utils.C"

void analyze(TString srcfile, int sample, float xsec_, int sampleevents_) {
  //TString filename = "root://cms-xrd-global.cern.ch/" + srcfile;
  TString filename = "root://xrootd-cms.infn.it/" + srcfile;
  TFile *f = TFile::Open(filename);
  if(f) cout<<"File has been opened!"<<endl;
  TFile *output= new TFile("testoutput.root","RECREATE");
  TTree* tree = (TTree*)f->Get("Events");
  TTree* runtree = (TTree*)f->Get("Runs");

  TTree *outtree = new TTree("outtree", "outtree");

  Int_t O_njets, O_sample;
  Float_t O_mvis, O_taupt, O_taueta, O_leppt, O_mjj, O_deltaRjj;
  Float_t O_tauphi, O_lepeta, O_lepphi, O_metpt, O_metphi;
  Float_t O_jet1eta, O_jet1phi, O_jet2eta, O_jet2phi;
  Float_t O_weight;
  UChar_t O_tau_genflav, O_lep_genflav;
 
  bool O_ismuon;
 
  outtree->Branch("weight",&O_weight,"weight/F");
  outtree->Branch("taupt",&O_taupt,"taupt/F");
  outtree->Branch("ismuon",&O_ismuon,"ismuon/O");
 
  outtree->Branch("taueta",  &O_taueta,  "taueta/F");
  outtree->Branch("tauphi",  &O_tauphi,  "tauphi/F");
  outtree->Branch("tau_genflav", &O_tau_genflav, "tau_genflav/b");
  outtree->Branch("lep_genflav", &O_lep_genflav, "lep_genflav/b");
 
  Int_t   O_nele;
  Float_t O_ele_pt_arr[20], O_ele_eta_arr[20];
  Int_t  O_ele_id_arr[20];
  Int_t   O_nmuon;
  Float_t O_muon_pt_arr[20], O_muon_eta_arr[20];
  Bool_t  O_muon_id_arr[20];
 
  Int_t O_eleindex, O_muindex;
 
  outtree->Branch("nele",   &O_nele, "nele/I");
  outtree->Branch("ele_pt_arr",       O_ele_pt_arr,     "ele_pt_arr[nele]/F");
  outtree->Branch("ele_eta_arr",      O_ele_eta_arr,    "ele_eta_arr[nele]/F");
  outtree->Branch("ele_id_arr",       O_ele_id_arr,     "ele_id_arr[nele]/I");
  outtree->Branch("nmuon",   &O_nmuon, "nmuon/I");
  outtree->Branch("muon_pt_arr",       O_muon_pt_arr,     "muon_pt_arr[nmuon]/F");
  outtree->Branch("muon_eta_arr",      O_muon_eta_arr,    "muon_eta_arr[nmuon]/F");
  outtree->Branch("muon_id_arr",       O_muon_id_arr,     "muon_id_arr[nmuon]/O");
 
  outtree->Branch("eleindex", &O_eleindex, "eleindex/I");
  outtree->Branch("muindex",  &O_muindex,  "muindex/I");
 
  //-------------------------- NEW: GenPart / GenJet matching for fake tau study --------------------
 
  Bool_t  O_tau_matched_genpart;   // se Tau_genPartIdx >= 0
  Int_t   O_tau_genpart_pdgid;     // pdgId della particella matchata (se esiste)
  Bool_t  O_tau_matched_genjet;    // se esiste un GenJet vicino (dR<0.4)
  Float_t O_tau_genjet_dr;         // dR tra tau ricostruito e GenJet piu' vicino
  Float_t O_tau_genjet_pt;         // pt del GenJet piu' vicino
  Float_t O_tau_genjet_eta;        // eta del GenJet piu' vicino
  Short_t O_tau_genjet_partonflav; // sapore del parton del GenJet (utile per capire se e' ISR da gluone/quark)
 
  outtree->Branch("tau_matched_genpart",   &O_tau_matched_genpart,   "tau_matched_genpart/O");
  outtree->Branch("tau_genpart_pdgid",     &O_tau_genpart_pdgid,     "tau_genpart_pdgid/I");
  outtree->Branch("tau_matched_genjet",    &O_tau_matched_genjet,    "tau_matched_genjet/O");
  outtree->Branch("tau_genjet_dr",         &O_tau_genjet_dr,         "tau_genjet_dr/F");
  outtree->Branch("tau_genjet_pt",         &O_tau_genjet_pt,         "tau_genjet_pt/F");
  outtree->Branch("tau_genjet_eta",        &O_tau_genjet_eta,        "tau_genjet_eta/F");
  outtree->Branch("tau_genjet_partonflav", &O_tau_genjet_partonflav, "tau_genjet_partonflav/S");

  //-------------------------------------------------------------------------------------------------

  tree->SetBranchStatus("*", 0);	//Turn off all the Branches and after turn on only what i need

  tree->SetBranchStatus("genWeight", 1);
  Float_t genweight_;
  tree->SetBranchAddress("genWeight",&genweight_);

  tree->SetBranchStatus("run", 1);
  UInt_t run_;
  tree->SetBranchAddress("run",&run_);

  tree->SetBranchStatus("event", 1);
  ULong64_t event_;
  tree->SetBranchAddress("event",&event_);

  tree->SetBranchStatus("luminosityBlock", 1);
  UInt_t ls_;
  tree->SetBranchAddress("luminosityBlock",&ls_);

//-------------------------- ELECTRONS ------------------------------------------------------------------------------------

  tree->SetBranchStatus("Electron_pt", 1);
  Float_t ele_pt_[128];
  tree->SetBranchAddress("Electron_pt",&ele_pt_);

  tree->SetBranchStatus("Electron_eta", 1);
  Float_t ele_eta_[128];
  tree->SetBranchAddress("Electron_eta",&ele_eta_);

  tree->SetBranchStatus("Electron_cutBased", 1);
  UChar_t ele_id_[128];
  tree->SetBranchAddress("Electron_cutBased",&ele_id_);

  tree->SetBranchStatus("Electron_charge", 1);
  Int_t ele_charge_[128];
  tree->SetBranchAddress("Electron_charge",&ele_charge_);

  tree->SetBranchStatus("Electron_mass", 1);
  Float_t ele_mass_[128];
  tree->SetBranchAddress("Electron_mass",&ele_mass_);

  tree->SetBranchStatus("Electron_phi", 1);
  Float_t ele_phi_[128];
  tree->SetBranchAddress("Electron_phi",&ele_phi_);

  tree->SetBranchStatus("Electron_dxy", 1);
  Float_t ele_dxy_[128];
  tree->SetBranchAddress("Electron_dxy",&ele_dxy_);

  tree->SetBranchStatus("Electron_dz", 1);
  Float_t ele_dz_[128];
  tree->SetBranchAddress("Electron_dz",&ele_dz_);

  tree->SetBranchStatus("Electron_convVeto", 1);
  Bool_t ele_conv_[128];
  tree->SetBranchAddress("Electron_convVeto",&ele_conv_);
  
  tree->SetBranchStatus("Electron_mvaIso_WP80", 1);
  Bool_t ele_mvaid_[128];
  tree->SetBranchAddress("Electron_mvaIso_WP80", &ele_mvaid_);

  tree->SetBranchStatus("Electron_r9", 1);
  Float_t ele_r9_[128];
  tree->SetBranchAddress("Electron_r9",&ele_r9_);

  tree->SetBranchStatus("Electron_seedGain", 1);
  UChar_t ele_gain_[128];
  tree->SetBranchAddress("Electron_seedGain",&ele_gain_);

  tree->SetBranchStatus("nElectron", 1);
  Int_t nelectrons_;
  tree->SetBranchAddress("nElectron",&nelectrons_);

//-------------------------- MUONS ----------------------------------------------------------------------------------------
  
  tree->SetBranchStatus("Muon_pt", 1);
  Float_t muon_pt_[128];
  tree->SetBranchAddress("Muon_pt",&muon_pt_);

  tree->SetBranchStatus("Muon_eta", 1);
  Float_t muon_eta_[128];
  tree->SetBranchAddress("Muon_eta",&muon_eta_);

  tree->SetBranchStatus("Muon_charge", 1);
  Int_t muon_charge_[128];
  tree->SetBranchAddress("Muon_charge",&muon_charge_);

  tree->SetBranchStatus("Muon_mass", 1);
  Float_t muon_mass_[128];
  tree->SetBranchAddress("Muon_mass",&muon_mass_);

  tree->SetBranchStatus("Muon_phi", 1);
  Float_t muon_phi_[128];
  tree->SetBranchAddress("Muon_phi",&muon_phi_);

  tree->SetBranchStatus("Muon_pfRelIso04_all", 1);
  Float_t muon_isoscore_[128];
  tree->SetBranchAddress("Muon_pfRelIso04_all",&muon_isoscore_);

  tree->SetBranchStatus("Muon_miniPFRelIso_all", 1);
  Float_t muon_iso_[128];
  tree->SetBranchAddress("Muon_miniPFRelIso_all", &muon_iso_);

  tree->SetBranchStatus("Muon_mediumId", 1);
  Bool_t muon_mediumid_[128];
  tree->SetBranchAddress("Muon_mediumId", &muon_mediumid_);
  
  tree->SetBranchStatus("Muon_dxy", 1);
  Float_t muon_dxy_[128];
  tree->SetBranchAddress("Muon_dxy", &muon_dxy_);

  tree->SetBranchStatus("Muon_dz", 1);
  Float_t muon_dz_[128];
  tree->SetBranchAddress("Muon_dz", &muon_dz_);

  tree->SetBranchStatus("Muon_nTrackerLayers", 1);
  UChar_t muon_ntracklayers_[128];
  tree->SetBranchAddress("Muon_nTrackerLayers", &muon_ntracklayers_);

  tree->SetBranchStatus("nMuon", 1);
  Int_t nmuons_;
  tree->SetBranchAddress("nMuon",&nmuons_);

  tree->SetBranchStatus("Muon_looseId", 1);
  Bool_t muon_looseid_[128];
  tree->SetBranchAddress("Muon_looseId",&muon_looseid_);

//-------------------------- TAUS -----------------------------------------------------------------------------------------
  
  tree->SetBranchStatus("Tau_pt", 1);
  Float_t tau_pt_[128];
  tree->SetBranchAddress("Tau_pt",&tau_pt_);

  tree->SetBranchStatus("Tau_eta", 1);
  Float_t tau_eta_[128];
  tree->SetBranchAddress("Tau_eta",&tau_eta_);

  tree->SetBranchStatus("Tau_mass", 1);
  Float_t tau_mass_[128];
  tree->SetBranchAddress("Tau_mass",&tau_mass_);

  tree->SetBranchStatus("Tau_phi", 1);
  Float_t tau_phi_[128];
  tree->SetBranchAddress("Tau_phi",&tau_phi_);

  tree->SetBranchStatus("Tau_decayMode", 1);
  UChar_t tau_decay_[128];
  tree->SetBranchAddress("Tau_decayMode",&tau_decay_);

  tree->SetBranchStatus("Tau_idDeepTau2018v2p5VSe", 1);
  UChar_t tauidvse_[128];
  tree->SetBranchAddress("Tau_idDeepTau2018v2p5VSe",&tauidvse_);

  tree->SetBranchStatus("Tau_idDeepTau2018v2p5VSmu", 1);
  UChar_t tauidvsmu_[128];
  tree->SetBranchAddress("Tau_idDeepTau2018v2p5VSmu",&tauidvsmu_);

  tree->SetBranchStatus("Tau_idDeepTau2018v2p5VSjet", 1);
  UChar_t tauidvsjet_[128];
  tree->SetBranchAddress("Tau_idDeepTau2018v2p5VSjet",&tauidvsjet_);

  tree->SetBranchStatus("Tau_charge", 1);
  Short_t tau_charge_[128];
  tree->SetBranchAddress("Tau_charge",&tau_charge_);

  tree->SetBranchStatus("Tau_dxy", 1);
  Float_t tau_dxy_[128];
  tree->SetBranchAddress("Tau_dxy",&tau_dxy_);

  tree->SetBranchStatus("Tau_dz", 1);
  Float_t tau_dz_[128];
  tree->SetBranchAddress("Tau_dz",&tau_dz_);

  tree->SetBranchStatus("nTau", 1);
  Int_t ntaus_;
  tree->SetBranchAddress("nTau",&ntaus_);

//-------------------------- MC TRUTH -------------------------------------------------------------------------------------

  tree->SetBranchStatus("Tau_genPartFlav", 1);
  UChar_t tau_source_[128];
  tree->SetBranchAddress("Tau_genPartFlav",&tau_source_);

  tree->SetBranchStatus("Electron_genPartFlav", 1);
  UChar_t ele_source_[128];
  tree->SetBranchAddress("Electron_genPartFlav",&ele_source_);

  tree->SetBranchStatus("Muon_genPartFlav", 1);
  UChar_t muon_source_[128];
  tree->SetBranchAddress("Muon_genPartFlav",&muon_source_);

 // -- NUOVO: Tau_genPartIdx per il matching diretto --
  tree->SetBranchStatus("Tau_genPartIdx", 1);
  Short_t tau_genpartidx_[128];
  tree->SetBranchAddress("Tau_genPartIdx", &tau_genpartidx_);
 
  // -- NUOVO: GenPart per identificare la particella matchata --
  tree->SetBranchStatus("GenPart_pdgId", 1);
  Int_t genpart_pdgid_[1024];
  tree->SetBranchAddress("GenPart_pdgId", &genpart_pdgid_);
 
  tree->SetBranchStatus("nGenPart", 1);
  Int_t ngenpart_;
  tree->SetBranchAddress("nGenPart", &ngenpart_);
 
  // -- NUOVO: GenJet per il matching spaziale (ISR check) --
  tree->SetBranchStatus("GenJet_pt", 1);
  Float_t genjet_pt_[128];
  tree->SetBranchAddress("GenJet_pt", &genjet_pt_);
 
  tree->SetBranchStatus("GenJet_eta", 1);
  Float_t genjet_eta_[128];
  tree->SetBranchAddress("GenJet_eta", &genjet_eta_);
 
  tree->SetBranchStatus("GenJet_phi", 1);
  Float_t genjet_phi_[128];
  tree->SetBranchAddress("GenJet_phi", &genjet_phi_);
 
  tree->SetBranchStatus("GenJet_partonFlavour", 1);
  Short_t genjet_partonflav_[128];
  tree->SetBranchAddress("GenJet_partonFlavour", &genjet_partonflav_);
 
  tree->SetBranchStatus("nGenJet", 1);
  Int_t ngenjet_;
  tree->SetBranchAddress("nGenJet", &ngenjet_);

//-------------------------- JETS -----------------------------------------------------------------------------------------

  tree->SetBranchStatus("Jet_pt", 1);
  Float_t jet_pt_[128];
  tree->SetBranchAddress("Jet_pt",&jet_pt_);

  tree->SetBranchStatus("Jet_eta", 1);
  Float_t jet_eta_[128];
  tree->SetBranchAddress("Jet_eta",&jet_eta_);

  tree->SetBranchStatus("Jet_phi", 1);
  Float_t jet_phi_[128];
  tree->SetBranchAddress("Jet_phi",&jet_phi_);
  
  tree->SetBranchStatus("Jet_mass", 1);
  Float_t jet_mass_[128];
  tree->SetBranchAddress("Jet_mass",&jet_mass_);

  tree->SetBranchStatus("Jet_hadronFlavour", 1);
  UChar_t jet_flav_[128];
  tree->SetBranchAddress("Jet_hadronFlavour",&jet_flav_);

  //tree->SetBranchStatus("Jet_btagDeepFlavB", 1); //Jet_btagDeepFlavB
  //Float_t jet_btag_[128];
  //tree->SetBranchAddress("Jet_btagDeepFlavB",&jet_btag_);

  tree->SetBranchStatus("Jet_btagUParTAK4B", 1); //Jet_btagUParTAK4B
  Float_t jet_btag_[128];
  tree->SetBranchAddress("Jet_btagUParTAK4B",&jet_btag_);

  tree->SetBranchStatus("Jet_rawFactor", 1);
  Float_t jet_raw_[128];
  tree->SetBranchAddress("Jet_rawFactor",&jet_raw_);

  tree->SetBranchStatus("GenJet_partonFlavour", 1);
  Short_t genjet_flav_[128];
  tree->SetBranchAddress("GenJet_partonFlavour",&genjet_flav_);

  tree->SetBranchStatus("nJet", 1);
  Int_t njets_;
  tree->SetBranchAddress("nJet",&njets_);

//-------------------------- HIGH LEVEL TRIGGER ---------------------------------------------------------------------------
  
  tree->SetBranchStatus("HLT_IsoMu24", 1);
  Bool_t mutri_;
  tree->SetBranchAddress("HLT_IsoMu24",&mutri_);

  tree->SetBranchStatus("HLT_IsoMu27", 1);
  Bool_t mutri2_;
  tree->SetBranchAddress("HLT_IsoMu27",&mutri2_);

  tree->SetBranchStatus("HLT_Ele30_WPTight_Gsf", 1);
  Bool_t eletri_;
  tree->SetBranchAddress("HLT_Ele30_WPTight_Gsf",&eletri_);

  tree->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1);
  Bool_t eletri2_;
  tree->SetBranchAddress("HLT_Ele32_WPTight_Gsf",&eletri2_);

//-------------------------- FLAGS ----------------------------------------------------------------------------------------

  tree->SetBranchStatus("Flag_goodVertices", 1);
  Bool_t flag1_;
  tree->SetBranchAddress("Flag_goodVertices",&flag1_);

  tree->SetBranchStatus("Flag_globalSuperTightHalo2016Filter", 1);
  Bool_t flag2_;
  tree->SetBranchAddress("Flag_globalSuperTightHalo2016Filter",&flag2_);

  tree->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
  Bool_t flag3_;
  tree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",&flag3_);

  tree->SetBranchStatus("Flag_BadPFMuonFilter", 1);
  Bool_t flag4_;
  tree->SetBranchAddress("Flag_BadPFMuonFilter",&flag4_);

  tree->SetBranchStatus("Flag_BadPFMuonDzFilter", 1);
  Bool_t flag5_;
  tree->SetBranchAddress("Flag_BadPFMuonDzFilter",&flag5_);

  tree->SetBranchStatus("Flag_hfNoisyHitsFilter", 1);
  Bool_t flag6_;
  tree->SetBranchAddress("Flag_hfNoisyHitsFilter",&flag6_);

  tree->SetBranchStatus("Flag_eeBadScFilter", 1);
  Bool_t flag7_;
  tree->SetBranchAddress("Flag_eeBadScFilter",&flag7_);

  tree->SetBranchStatus("Flag_ecalBadCalibFilter", 1);
  Bool_t flag8_;
  tree->SetBranchAddress("Flag_ecalBadCalibFilter",&flag8_);

//-------------------------- RANDOM STUFF ---------------------------------------------------------------------------------

  tree->SetBranchStatus("Rho_fixedGridRhoFastjetCentralCalo", 1);
  Float_t rho_calo_;
  tree->SetBranchAddress("Rho_fixedGridRhoFastjetCentralCalo",&rho_calo_);
  
  tree->SetBranchStatus("Pileup_nTrueInt", 1);
  Float_t npu2_;
  tree->SetBranchAddress("Pileup_nTrueInt",&npu2_);

  runtree->SetBranchStatus("genEventSumw", 1);
  Double_t sumgenw_;
  runtree->SetBranchAddress("genEventSumw",&sumgenw_);

  tree->SetBranchStatus("PuppiMET_pt",  1);
  Float_t met_pt_;
  tree->SetBranchAddress("PuppiMET_pt", &met_pt_);

  tree->SetBranchStatus("PuppiMET_phi",  1);
  Float_t met_phi_;
  tree->SetBranchAddress("PuppiMET_phi", &met_phi_);

//-------------------------- JET ID ---------------------------------------------------------------------------------------

  tree->SetBranchStatus("Jet_neHEF", 1);
  Float_t jet_neHEF_[128];
  tree->SetBranchAddress("Jet_neHEF", &jet_neHEF_);

  tree->SetBranchStatus("Jet_chHEF", 1);
  Float_t jet_chHEF_[128];
  tree->SetBranchAddress("Jet_chHEF", &jet_chHEF_);

  tree->SetBranchStatus("Jet_neEmEF", 1);
  Float_t jet_neEmEF_[128];
  tree->SetBranchAddress("Jet_neEmEF", &jet_neEmEF_);

  tree->SetBranchStatus("Jet_chEmEF", 1);
  Float_t jet_chEmEF_[128];
  tree->SetBranchAddress("Jet_chEmEF", &jet_chEmEF_);

  tree->SetBranchStatus("Jet_muEF", 1);
  Float_t jet_muEF_[128];
  tree->SetBranchAddress("Jet_muEF", &jet_muEF_);

  tree->SetBranchStatus("Jet_neMultiplicity", 1);
  UChar_t jet_neMultiplicity_[128];
  tree->SetBranchAddress("Jet_neMultiplicity", &jet_neMultiplicity_);

  tree->SetBranchStatus("Jet_chMultiplicity", 1);
  UChar_t jet_chMultiplicity_[128];
  tree->SetBranchAddress("Jet_chMultiplicity", &jet_chMultiplicity_);

  runtree->GetEntry(0);

  float weightscale_=1/sumgenw_;

  float lumi = 220.79; //2024+2025 luminosity (in fb^-1)

  weightscale_*=lumi*xsec_*1000;
  Long64_t numEntries = tree->GetEntries();
  weightscale_*=numEntries/(float)sampleevents_;

  TRandom3 *tr3 = new TRandom3;

  auto startTime = std::chrono::high_resolution_clock::now();

  cout<<"Started!"<<endl;

  float totalll1=0., totalll2=0.;
  int nummm=0;

  for (Long64_t i = 0; i < numEntries; ++i) {

    //if(i>5000) break;
    if(i%10000==0) cout<<"Analyzing event n. "<<i<<"/"<<numEntries<<endl;
    
    tree->GetEntry(i);

    bool excflag=0;
	bool btagflag=0;

    Bool_t METfilters= (flag1_ && flag2_ && flag3_ && flag4_ && flag5_ && flag6_ && flag7_ && flag8_);
    if(!METfilters) excflag=1;

    Float_t weight_=genweight_*weightscale_;

    weight_*=pu_SF->evaluate({npu2_,"nominal"});

    Int_t lepton=0, jets;

    Int_t ntaus=0, nbtags=0, taucharge=0, nelectrons=0, nmuons=0, lepcharge=0, njets=0;
    Float_t selectedtaupt=0., selectedleppt=0.;
	Int_t jet1index = -1, jet2index = -1;
	Bool_t ismuon = true;
	Int_t   nleps=0;

    ROOT::Math::PtEtaPhiMVector p4tau, p4lep, p4jet1, p4jet2;
    
    int tauindex=0;
    for(int j=0; j<ntaus_; j++){
      Float_t taupt = tau_pt_[j];
      Bool_t pass = TauSelector(taupt, tau_eta_[j], tauidvse_[j], tauidvsmu_[j], tauidvsjet_[j], tau_source_[j], tau_decay_[j], tau_dz_[j], weight_);
      if(pass){
		ntaus++;
        tauindex=j;
		taucharge=tau_charge_[j];
		p4tau = ROOT::Math::PtEtaPhiMVector(taupt,tau_eta_[j],tau_phi_[j],tau_mass_[j]);
      }
    }
        
    int typeevent=0; //1=mutauh, 2=eletauh
    int eleindex=200;
    int muindex=200;
      
    if(ntaus == 1){

    	for(int j=0; j<nelectrons_; j++){
			Float_t elept = ele_pt_[j];
			Bool_t pass = ElectronSelector(elept, ele_eta_[j], ele_phi_[j], ele_mvaid_[j], ele_dxy_[j], ele_dz_[j], ele_conv_[j], ele_r9_[j], ele_gain_[j], run_, weight_);
			if(pass){
	  			nelectrons++;
	  			eleindex=j;
	  			lepcharge=ele_charge_[j];
	  			p4lep = ROOT::Math::PtEtaPhiMVector(elept,ele_eta_[eleindex],ele_phi_[eleindex],ele_mass_[eleindex]);
			}
    	}
      
    	for(int j=0; j<nmuons_; j++){
			Float_t muonpt = muon_pt_[j];
			Bool_t pass = MuonSelector(muonpt,muon_eta_[j],muon_phi_[j],muon_mediumid_[j],muon_dxy_[j],muon_dz_[j],muon_isoscore_[j],muon_charge_[j],muon_ntracklayers_[j], event_, ls_, weight_);
        	if(pass){
        		nmuons++;
	  			muindex=j;
	  			lepcharge=muon_charge_[j];
	  			p4lep = ROOT::Math::PtEtaPhiMVector(muonpt,muon_eta_[muindex],muon_phi_[muindex],muon_mass_[muindex]);
        	}
    	}

		if(nelectrons_>0){
			for(int l=0; l<nelectrons_; l++){
  			    O_ele_pt_arr[l]  = ele_pt_[l];
				O_ele_eta_arr[l] = ele_eta_[l];
            	O_ele_id_arr[l]  = static_cast<int>(ele_id_[l]);
 			}
		}

		if(nmuons_>0){
			for(int l=0; l<nmuons_; l++){
  			    O_muon_pt_arr[l]  = muon_pt_[l];
				O_muon_eta_arr[l] = muon_eta_[l];
            	O_muon_id_arr[l]  = muon_looseid_[l];
 			}
		}

		//----------------------- Veto on additional Loose leptons -------------------------------------

		for(int i=0; i<nelectrons_; i++){
  			if(i==eleindex) continue;
  			if(ele_pt_[i] > 10 && abs(ele_eta_[i]) < 2.4 && static_cast<int>(ele_id_[i])>=1) excflag=1;
 		}

		for(int i=0; i<nmuons_; i++){
  			if(i==muindex) continue;
  			if(muon_pt_[i] > 10 && abs(muon_eta_[i]) < 2.4 && muon_looseid_[i]) excflag=1;
 		}

		//----------------------------------------------------------------------------------------------

    	if(nelectrons+nmuons==1 and tau_charge_[tauindex]==lepcharge){
			O_tau_genflav = tau_source_[tauindex];
    		if(nmuons==1){ 
				typeevent=1; //Mu-tauh events
				ismuon = true;
				O_lep_genflav = muon_source_[muindex];
			}else{ 
				typeevent=2; //E-tauh events
				ismuon = false;
				O_lep_genflav = ele_source_[eleindex];
			}

          //--------------- NUOVO: GenPart / GenJet matching per il tau ---------------
 
            O_tau_matched_genpart = false;
            O_tau_genpart_pdgid = 0;
            O_tau_matched_genjet = false;
        	O_tau_genjet_dr = -1;
        	O_tau_genjet_pt = -1;
        	O_tau_genjet_eta = -99;
        	O_tau_genjet_partonflav = -99;
 
        	Short_t genidx = tau_genpartidx_[tauindex];
        	if(genidx >= 0 && genidx < ngenpart_){
          		O_tau_matched_genpart = true;
          		O_tau_genpart_pdgid = genpart_pdgid_[genidx];
        	}
 
        	// Cerca il GenJet piu' vicino al tau ricostruito (utile soprattutto se unmatched a livello di GenPart)
        	float min_dr = 999;
        	int best_genjet = -1;
        	for(int g=0; g<ngenjet_ && g<128; g++){
          		ROOT::Math::PtEtaPhiMVector p4genjet(genjet_pt_[g], genjet_eta_[g], genjet_phi_[g], 0);
          		float dr = deltaR(p4tau, p4genjet);
          		if(dr < min_dr){
            		min_dr = dr;
            		best_genjet = g;
          		}
        	}
        	if(best_genjet >= 0 && min_dr < 0.4){
          		O_tau_matched_genjet = true;
          		O_tau_genjet_dr = min_dr;
          		O_tau_genjet_pt = genjet_pt_[best_genjet];
          		O_tau_genjet_eta = genjet_eta_[best_genjet];
          		O_tau_genjet_partonflav = genjet_partonflav_[best_genjet];
        	}
    	}
      
    	for(int j=0; j<njets_; j++){
			
			ROOT::Math::PtEtaPhiMVector p4jet(jet_pt_[j],jet_eta_[j],jet_phi_[j],jet_mass_[j]);
			if(deltaR(p4jet,p4tau)<0.4) { // Reject jets that overlap with the tau
	  			continue;
			}
			if(deltaR(p4jet,p4lep)<0.4) { // Reject jets that overlap with the other lepton
	  			continue;
			}
			Float_t jetpt = jet_pt_[j];
			Bool_t pass = JetSelector(jetpt,jet_eta_[j],jet_phi_[j],jet_raw_[j],rho_calo_);
			if(pass) {
				if(!JetIdTightLepVeto(jet_eta_[j], jet_neHEF_[j], jet_neEmEF_[j], jet_chEmEF_[j], jet_muEF_[j], jet_chHEF_[j], jet_neMultiplicity_[j], jet_chMultiplicity_[j])) continue;

	  			njets++;

				if(jet_btag_[j] > WP_M && TMath::Abs(jet_eta_[j]) < 2.5)	btagflag = 1;

	  			if(njets==1){
					p4jet1 = ROOT::Math::PtEtaPhiMVector(jetpt,jet_eta_[j],jet_phi_[j],jet_mass_[j]);
				}else if(njets==2){
					p4jet2 = ROOT::Math::PtEtaPhiMVector(jetpt,jet_eta_[j],jet_phi_[j],jet_mass_[j]);
				}

			}
    	}
	}
    
    bool trigpath=false; 
    if(typeevent==1) trigpath=mutri_; 
    else if(typeevent==2) trigpath=eletri_;
    
    if(trigpath and typeevent>0 and !excflag and njets>=2 and !btagflag){

      float puweight=pu_SF->evaluate({npu2_,"nominal"});

      ROOT::Math::PtEtaPhiMVector p4=p4lep+p4tau;
      float mvis_ = p4.M(); // (Visible) invariant mass of the e-tau or mu-tau system

      ROOT::Math::PtEtaPhiMVector p4jets=p4jet1+p4jet2;
      float mjj_ = p4jets.M(); // Invariant mass of the dijet system
      float deltaRjj_ = deltaR(p4jet1,p4jet2);
	
      O_weight   = weight_;
      O_taupt    = p4tau.Pt();
	  O_taueta   = p4tau.Eta();
	  O_tauphi   = p4tau.Phi();
	  O_ismuon   = ismuon;

	  O_nele     = nelectrons_;
	  O_nmuon    = nmuons_;
	  O_eleindex = eleindex;
      O_muindex  = muindex;
      outtree->Fill();

    }
 }
  
  outtree->Write();
  auto endTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

  std::cout << "Total execution time: " << duration << " milliseconds" << std::endl;

  f->Close();
  output->Close();
}


int main(int argc, char* argv[]) {

  TString srcfile=argv[1];
  int sample= std::atoi(argv[2]);
  float xsec_= std::atof(argv[3]);
  int sampleevents_= std::atoi(argv[4]);

  analyze(srcfile,sample,xsec_,sampleevents_);

  return 1;
}
