#include <iostream>
#include <map>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1I.h"
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
  TFile *output= new TFile("testoutput.root","RECREATE");
  TTree* tree = (TTree*)f->Get("Events");
  TTree* runtree = (TTree*)f->Get("Runs");

  TTree *outtree = new TTree("outtree", "outtree");

  tree->SetBranchStatus("*", 0);	//Turn off all the Branches and after turn on only what i need

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

//-------------------------- HIGH LEVEL TRIGGER (HLT)---------------------------------------------------------------------------
  
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

  tree->SetBranchStatus("HLT_PFHT180", 1);
  Bool_t HLT_PFHT180_;
  tree->SetBranchAddress("HLT_PFHT180",&HLT_PFHT180_);

  tree->SetBranchStatus("HLT_PFHT250", 1);
  Bool_t HLT_PFHT250_;
  tree->SetBranchAddress("HLT_PFHT250",&HLT_PFHT250_);

  tree->SetBranchStatus("HLT_PFHT350", 1);
  Bool_t HLT_PFHT350_;
  tree->SetBranchAddress("HLT_PFHT350",&HLT_PFHT350_);

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

//-------------------------- MET ---------------------------------------------------------------------------------------

  tree->SetBranchStatus("PuppiMET_pt",  1);
  Float_t met_pt_;
  tree->SetBranchAddress("PuppiMET_pt", &met_pt_);

  tree->SetBranchStatus("PuppiMET_phi",  1);
  Float_t met_phi_;
  tree->SetBranchAddress("PuppiMET_phi", &met_phi_);

  runtree->GetEntry(0);

  Long64_t numEntries = tree->GetEntries();

  auto goldenMap = loadGoldenJSON(Goldenjson_2024);

//-------------------------- OUTPUT HISTOS ------------------------------------------------------------------------------

	Float_t edge_pt1[]  = {20.,25.,30.,35.,40.,45.,50.,55.,60.,70.,80.,100.,120.,140.,160.,180.,200.};
	Float_t edge_pt2[]  = {20.,40.,60.,80.,100.,120.,145.,200.};
	Float_t edge_eta1[]  = {-2.5,-2.,-1.479,-1.,0.,1.,1.479,2.,2.5};
	Float_t edge_eta2[]  = {0.,1.,1.479,2.,2.5};

	const int n_pt1  = sizeof(edge_pt1)/sizeof(edge_pt1[0])  - 1;
	const int n_pt2  = sizeof(edge_pt2)/sizeof(edge_pt2[0])  - 1;
	const int n_eta1  = sizeof(edge_eta1)/sizeof(edge_eta1[0]) - 1;
	const int n_eta2  = sizeof(edge_eta2)/sizeof(edge_eta2[0]) - 1;
	
	TH1F *h_tau_pt_den   = new TH1F("h_tau_pt_den","h_tau_pt_den",n_pt1,edge_pt1);
	TH1F *h_tau_eta_den  = new TH1F("h_tau_eta_den","h_tau_eta_den",n_eta1,edge_eta1);
	TH2F *h_tau_2d_den   = new TH2F("h_tau_2d_den","h_tau_2d_den",n_pt2,edge_pt2,n_eta2,edge_eta2);
	
	TH1F *h_tau_pt_num   = new TH1F("h_tau_pt_num","h_tau_pt_num",n_pt1,edge_pt1);
	TH1F *h_tau_eta_num  = new TH1F("h_tau_eta_num","h_tau_eta_num",n_eta1,edge_eta1);
	TH2F *h_tau_2d_num   = new TH2F("h_tau_2d_num","h_tau_2d_num",n_pt2,edge_pt2,n_eta2,edge_eta2);
	
	TH1F *h_ele_pt_den   = new TH1F("h_ele_pt_den","h_ele_pt_den",n_pt1,edge_pt1);
	TH1F *h_ele_eta_den  = new TH1F("h_ele_eta_den","h_ele_eta_den",n_eta1,edge_eta1);
	TH2F *h_ele_2d_den   = new TH2F("h_ele_2d_den","h_ele_2d_den",n_pt2,edge_pt2,n_eta2,edge_eta2);
	
	TH1F *h_ele_pt_num   = new TH1F("h_ele_pt_num","h_ele_pt_num",n_pt1,edge_pt1);
	TH1F *h_ele_eta_num  = new TH1F("h_ele_eta_num","h_ele_eta_num",n_eta1,edge_eta1);
	TH2F *h_ele_2d_num   = new TH2F("h_ele_2d_num","h_ele_2d_num",n_pt2,edge_pt2,n_eta2,edge_eta2);
	
	TH1F *h_muon_pt_den  = new TH1F("h_muon_pt_den","h_muon_pt_den",n_pt1,edge_pt1);
	TH1F *h_muon_eta_den = new TH1F("h_muon_eta_den","h_muon_eta_den",n_eta1,edge_eta1);
	TH2F *h_muon_2d_den  = new TH2F("h_muon_2d_den","h_muon_2d_den",n_pt2,edge_pt2,n_eta2,edge_eta2);
	
	TH1F *h_muon_pt_num  = new TH1F("h_muon_pt_num","h_muon_pt_num",n_pt1,edge_pt1);
	TH1F *h_muon_eta_num = new TH1F("h_muon_eta_num","h_muon_eta_num",n_eta1,edge_eta1);
	TH2F *h_muon_2d_num  = new TH2F("h_muon_2d_num","h_muon_2d_num",n_pt2,edge_pt2,n_eta2,edge_eta2);

//-------------------------------------------------------------------------------------------------------------------------

  for (Long64_t i = 0; i < numEntries; ++i) {
    tree->GetEntry(i);

    Bool_t golden_event = is_valid_event(goldenMap, run_, ls_);
    if(!golden_event) continue;
    
    Bool_t excflag = 0;

    Bool_t METfilters = (flag1_ && flag2_ && flag3_ && flag4_ && flag5_ && flag6_ && flag7_ && flag8_);
    Bool_t trigpath = (HLT_PFJet40_ || HLT_PFJet60_ || HLT_PFJet80_ || HLT_PFJet110_ || HLT_PFJet140_ || HLT_PFJet200_ || HLT_PFJet260_ || HLT_PFHT180_ || HLT_PFHT250_ || HLT_PFHT350_);
    if(!METfilters || !trigpath) continue;

	Int_t ntaus=0, neles=0, nmuons=0;

    //-------------------------- TAU ----------------------------------------------------------------
    for(int j=0; j<ntaus_; j++){
	  if(ntaus == 1) break;

      int vse   = static_cast<int>(tauidvse_[j]);
      int vsmu  = static_cast<int>(tauidvsmu_[j]);
      int vsjet = static_cast<int>(tauidvsjet_[j]);
      int DM    = static_cast<int>(tau_decay_[j]);
      if(DM==2 || DM==5 || DM==6) continue;

      if(vse>=6 && vsmu>=4 && tau_pt_[j]>20 && abs(tau_eta_[j])<2.3 && abs(tau_dz_[j])<0.2){
		
		Float_t mT = m_T(tau_pt_[j], met_phi_, tau_phi_[j], met_phi_);
		if(mT > 50) continue;

        if(vsjet>=4){
          h_tau_pt_den->Fill(tau_pt_[j]);
          h_tau_eta_den->Fill(tau_eta_[j]);
          h_tau_2d_den->Fill(tau_pt_[j], abs(tau_eta_[j]));
		  		ntaus++;
        }
        if(vsjet>=5){
          h_tau_pt_num->Fill(tau_pt_[j]);
          h_tau_eta_num->Fill(tau_eta_[j]);
          h_tau_2d_num->Fill(tau_pt_[j], abs(tau_eta_[j]));
        }
      }
    }
 
    //-------------------------- ELECTRON -----------------------------------------------------------
    for(int j=0; j<nelectrons_; j++){
		if(neles == 1) break;

		if(ele_pt_[j]>30 && abs(ele_eta_[j])<2.5 && abs(ele_dxy_[j])<0.1 && abs(ele_dz_[j])<0.2 && ele_conv_[j]){

		Float_t mT = m_T(ele_pt_[j], met_phi_, ele_phi_[j], met_phi_);
	  	if(mT > 50) continue;

        if(ele_mvaid90_[j]){
          h_ele_pt_den->Fill(ele_pt_[j]);
          h_ele_eta_den->Fill(ele_eta_[j]);
          h_ele_2d_den->Fill(ele_pt_[j], abs(ele_eta_[j]));
		  		neles++;
        }
        if(ele_mvaid80_[j]){
          h_ele_pt_num->Fill(ele_pt_[j]);
          h_ele_eta_num->Fill(ele_eta_[j]);
          h_ele_2d_num->Fill(ele_pt_[j], abs(ele_eta_[j]));
        }
      }
    }
 
    //-------------------------- MUON ---------------------------------------------------------------
    for(int j=0; j<nmuons_; j++){
		if(nmuons == 1) break;

 		if(muon_pt_[j]>30 && abs(muon_eta_[j]) < 2.4 && abs(muon_dxy_[j])<0.1 && abs(muon_dz_[j])<0.2){

		Float_t mT = m_T(muon_pt_[j], met_phi_, muon_phi_[j], met_phi_);
	  	if(mT > 50) continue;

        if(muon_looseid_[j] && muon_isoscore_[j]<0.4){
          h_muon_pt_den->Fill(muon_pt_[j]);
          h_muon_eta_den->Fill(muon_eta_[j]);
          h_muon_2d_den->Fill(muon_pt_[j], abs(muon_eta_[j]));
		  nmuons++;
        }
        if(muon_mediumid_[j] && muon_isoscore_[j]<0.15){
          h_muon_pt_num->Fill(muon_pt_[j]);
          h_muon_eta_num->Fill(muon_eta_[j]);
          h_muon_2d_num->Fill(muon_pt_[j], abs(muon_eta_[j]));
        }
      }
    }
	}
 
  h_tau_pt_den->Write();   h_tau_eta_den->Write();   h_tau_2d_den->Write();
  h_tau_pt_num->Write();   h_tau_eta_num->Write();   h_tau_2d_num->Write();
  h_ele_pt_den->Write();   h_ele_eta_den->Write();   h_ele_2d_den->Write();
  h_ele_pt_num->Write();   h_ele_eta_num->Write();   h_ele_2d_num->Write();
  h_muon_pt_den->Write();  h_muon_eta_den->Write();  h_muon_2d_den->Write();
  h_muon_pt_num->Write();  h_muon_eta_num->Write();  h_muon_2d_num->Write();

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
