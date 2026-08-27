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
  Float_t O_jeteta[20], O_jetphi[20], O_jetpt[20];
  Float_t O_weight;
  UChar_t O_tau_genflav, O_lep_genflav;
  bool O_ismuon, O_excflag;
  Bool_t O_lep_charge_flip;
  Int_t  O_lep_gen_charge;
  Int_t  O_nGenJet;
  Float_t O_minor_deltaR[128];

  outtree->Branch("weight",&O_weight,"weight/F");
  outtree->Branch("taupt",&O_taupt,"taupt/F");
  outtree->Branch("taueta",  &O_taueta,  "taueta/F");
  outtree->Branch("tauphi",  &O_tauphi,  "tauphi/F");
  outtree->Branch("lepeta",  &O_lepeta,  "lepeta/F");
  outtree->Branch("lepphi",  &O_lepphi,  "lepphi/F");
  outtree->Branch("njets",&O_njets,"njets/I");
  outtree->Branch("sample",&O_sample,"sample/I");
  outtree->Branch("ismuon",&O_ismuon,"ismuon/O");
  outtree->Branch("tau_genflav",     &O_tau_genflav,     "tau_genflav/b");
  outtree->Branch("lep_genflav",     &O_lep_genflav,     "lep_genflav/b");
  outtree->Branch("nGenJet",         &O_nGenJet,         "nGenJet/I");
  outtree->Branch("jetpt",           O_jetpt,            "jetpt[njets]/F");
  outtree->Branch("jeteta",          O_jeteta,           "jeteta[njets]/F");
  outtree->Branch("jetphi",          O_jetphi,           "jetphi[njets]/F");

  Int_t   O_nGenLep;
  Int_t   O_genlep_pdgid[20];
  Float_t O_genlep_pt[20];
  Float_t O_genlep_eta[20];
  Float_t O_genlep_phi[20];

 outtree->Branch("nGenLep",       &O_nGenLep,      "nGenLep/I"); 
 outtree->Branch("genlep_pdgid",  O_genlep_pdgid,  "genlep_pdgid[nGenLep]/I");
 outtree->Branch("genlep_pt",     O_genlep_pt,     "genlep_pt[nGenLep]/F");
 outtree->Branch("genlep_eta",    O_genlep_eta,    "genlep_eta[nGenLep]/F");
 outtree->Branch("genlep_phi",    O_genlep_phi,    "genlep_phi[nGenLep]/F");

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

  tree->SetBranchStatus("Jet_btagUParTAK4B", 1); //Jet_btagUParTAK4B
  Float_t jet_btag_[128];
  tree->SetBranchAddress("Jet_btagUParTAK4B",&jet_btag_);

  tree->SetBranchStatus("Jet_rawFactor", 1);
  Float_t jet_raw_[128];
  tree->SetBranchAddress("Jet_rawFactor",&jet_raw_);

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

  tree->SetBranchStatus("HLT_PFJet40", 1);
  Bool_t HLT_PFJet40_;
  tree->SetBranchAddress("HLT_PFJet40",&HLT_PFJet40_);

  tree->SetBranchStatus("HLT_PFJet60", 1);
  Bool_t HLT_PFJet60_;
  tree->SetBranchAddress("HLT_PFJet60",&HLT_PFJet60_);

  tree->SetBranchStatus("HLT_PFJet80", 1);
  Bool_t HLT_PFJet80_;
  tree->SetBranchAddress("HLT_PFJet80",&HLT_PFJet80_);

  tree->SetBranchStatus("HLT_PFJet110", 1);
  Bool_t HLT_PFJet110_;
  tree->SetBranchAddress("HLT_PFJet110",&HLT_PFJet110_);

  tree->SetBranchStatus("HLT_PFJet140", 1);
  Bool_t HLT_PFJet140_;
  tree->SetBranchAddress("HLT_PFJet140",&HLT_PFJet140_);

  tree->SetBranchStatus("HLT_PFJet200", 1);
  Bool_t HLT_PFJet200_;
  tree->SetBranchAddress("HLT_PFJet200",&HLT_PFJet200_);

  tree->SetBranchStatus("HLT_PFJet260", 1);
  Bool_t HLT_PFJet260_;
  tree->SetBranchAddress("HLT_PFJet260",&HLT_PFJet260_);

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

//-------------------------- OUTPUT TTREE ---------------------------------------------------------------------------------

  Int_t O_ntaus_den, O_ntaus_num;
  Float_t O_taupt_den[20], O_taueta_den[20], O_tauphi_den[20];
  Float_t O_taupt_num[20], O_taueta_num[20], O_tauphi_num[20];
  Int_t O_tauDM_den[20],O_tauDM_num[20];

  Int_t O_neles_den, O_neles_num;
  Float_t O_elept_den[20], O_eleeta_den[20], O_elephi_den[20];
  Float_t O_elept_num[20], O_eleeta_num[20], O_elephi_num[20];

  Int_t O_nmuons_den, O_nmuons_num;
  Float_t O_muonpt_den[20], O_muoneta_den[20], O_muonphi_den[20];
  Float_t O_muonpt_num[20], O_muoneta_num[20], O_muonphi_num[20];

  outtree->Branch("ntaus_den",      &O_ntaus_den,   "ntaus_den/I");
  outtree->Branch("taupt_den",       O_taupt_den,   "taupt_den[ntaus_den]/F");
  outtree->Branch("taueta_den",      O_taueta_den,  "taueta_den[ntaus_den]/F");
  outtree->Branch("tauphi_den",      O_tauphi_den,  "tauphi_den[ntaus_den]/F");
  outtree->Branch("tauDM_den",       O_tauDM_den,   "tauDM_den[ntaus_den]/I");

  outtree->Branch("ntaus_num",      &O_ntaus_num,   "ntaus_num/I");
  outtree->Branch("taupt_num",       O_taupt_num,   "taupt_num[ntaus_num]/F");
  outtree->Branch("taueta_num",      O_taueta_num,  "taueta_num[ntaus_num]/F");
  outtree->Branch("tauphi_num",      O_tauphi_num,  "tauphi_num[ntaus_num]/F");
  outtree->Branch("tauDM_num",       O_tauDM_num,   "tauDM_num[ntaus_num]/I");

  outtree->Branch("neles_den",      &O_neles_den,   "neles_den/I");
  outtree->Branch("elept_den",       O_elept_den,   "elept_den[neles_den]/F");
  outtree->Branch("eleeta_den",      O_eleeta_den,  "eleeta_den[neles_den]/F");
  outtree->Branch("elephi_den",      O_elephi_den,  "elephi_den[neles_den]/F");

  outtree->Branch("neles_num",      &O_neles_num,   "neles_num/I");
  outtree->Branch("elept_num",       O_elept_num,   "elept_num[neles_num]/F");
  outtree->Branch("eleeta_num",      O_eleeta_num,  "eleeta_num[neles_num]/F");
  outtree->Branch("elephi_num",      O_elephi_num,  "elephi_num[neles_num]/F");

  outtree->Branch("nmuons_den",      &O_nmuons_den,   "nmuons_den/I");
  outtree->Branch("muonpt_den",       O_muonpt_den,   "muonpt_den[nmuons_den]/F");
  outtree->Branch("muoneta_den",      O_muoneta_den,  "muoneta_den[nmuons_den]/F");
  outtree->Branch("muonphi_den",      O_muonphi_den,  "muonphi_den[nmuons_den]/F");

  outtree->Branch("nmuons_num",      &O_nmuons_num,   "nmuons_num/I");
  outtree->Branch("muonpt_num",       O_muonpt_num,   "muonpt_num[nmuons_num]/F");
  outtree->Branch("muoneta_num",      O_muoneta_num,  "muoneta_num[nmuons_num]/F");
  outtree->Branch("muonphi_num",      O_muonphi_num,  "muonphi_num[nmuons_num]/F");

\\-------------------------------------------------------------------------------------------------------------------------
  runtree->GetEntry(0);

  float weightscale_=1/sumgenw_;

  Long64_t numEntries = tree->GetEntries();
  weightscale_*=numEntries/(float)sampleevents_;

  for (Long64_t i = 0; i < numEntries; ++i) {
    tree->GetEntry(i);

    bool excflag=0;

    Bool_t METfilters= (flag1_ && flag2_ && flag3_ && flag4_ && flag5_ && flag6_ && flag7_ && flag8_);
    Float_t lumi_eff = trigpath_Jet(HLT_PFJet40_, HLT_PFJet60_, HLT_PFJet80_, HLT_PFJet110_, HLT_PFJet140_, HLT_PFJet200_, HLT_PFJet260_);
    if(!METfilters || Lumi_eff == -200.) continue;

    weightscale_*=lumi_eff*xsec_*1000;

    Float_t weight_=genweight_*weightscale_;
    weight_*=pu_SF->evaluate({npu2_,"nominal"});

    vector<int> tauidx_den;
    vector<int> tauidx_num;
    vector<int> eleidx_den;
    vector<int> eleidx_num;
    vector<int> muonidx_den;
    vector<int> muonidx_num;
 
    
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
