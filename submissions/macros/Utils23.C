#include <iostream>
#include <string>
#include "TString.h"
#include <vector>
#include <array>
#include <map>
#include <fstream>
#include <algorithm>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

string JEC_json = "jsons_archive_2023/jet_jerc_2023BPix.json.gz";
auto JEC_c_set = CorrectionSet::from_file(JEC_json);

string MCera="Summer23BPixPrompt23";
string MCeraJER="Summer23BPixPrompt23_RunD";
string tauera="2022_postEE";

auto JEC_L2MC = JEC_c_set->at(MCera+"_V1_MC_L2Relative_AK4PFPuppi");

auto JER = JEC_c_set->at(MCeraJER+"_JRV1_MC_PtResolution_AK4PFPuppi");
auto JER_SF = JEC_c_set->at(MCeraJER+"_JRV1_MC_ScaleFactor_AK4PFPuppi");

std::array<decltype(JEC_c_set->begin()->second), 27> jec_syst = {{
    JEC_c_set->at(MCera+"_V1_MC_AbsoluteStat_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_AbsoluteScale_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_AbsoluteMPFBias_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_Fragmentation_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_SinglePionECAL_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_SinglePionHCAL_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_FlavorQCD_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_TimePtEta_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeJEREC1_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeJEREC2_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeJERHF_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativePtBB_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativePtEC1_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativePtEC2_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativePtHF_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeBal_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeSample_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeFSR_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeStatFSR_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeStatEC_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_RelativeStatHF_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_PileUpDataMC_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_PileUpPtRef_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_PileUpPtBB_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_PileUpPtEC1_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_PileUpPtEC2_AK4PFPuppi"),
    JEC_c_set->at(MCera+"_V1_MC_PileUpPtHF_AK4PFPuppi")
  }};

string JET_vetomap_2023 = "jsons/jetvetomaps_2023BPix.json.gz";
auto JETveto_c_set_2023 = CorrectionSet::from_file(JET_vetomap_2023);
auto JET_veto_2023      = JETveto_c_set_2023->at("Summer23BPixPrompt23_RunD_V1");

string JET_vetomap_2024 = "jsons/jetvetomaps_2024.json.gz";
auto JETveto_c_set_2024 = CorrectionSet::from_file(JET_vetomap_2024);
auto JET_veto_2024      = JETveto_c_set_2024->at("Summer24Prompt24_RunBCDEFGHI_V1");

string JET_vetomap_2025 = "jsons/jetvetomaps_2025.json.gz";
auto JETveto_c_set_2025 = CorrectionSet::from_file(JET_vetomap_2025);
auto JET_veto_2025      = JETveto_c_set_2025->at("Summer24Prompt25_RunCDEFG_V1");

string pu_file= "jsons_archive_2023/puWeights_2023BPix.json.gz";
string pucorr="Collisions2023_369803_370790_eraD_GoldenJson";
auto pu_c_set = CorrectionSet::from_file(pu_file);

auto pu_SF = pu_c_set->at(pucorr);

string bjet_file= "jsons_archive_2023/btagging_2023BPix.json.gz";
auto bjet_c_set = CorrectionSet::from_file(bjet_file);
auto bjet_SF = bjet_c_set->at("deepJet_comb"); //deepJet                                                                                                                                                    
auto lightjet_SF = bjet_c_set->at("deepJet_light");
auto bjet_shape = bjet_c_set->at("deepJet_shape");

auto btag_thr_getter=bjet_c_set->at("deepJet_wp_values");
float btag_thr_=btag_thr_getter->evaluate({"L"});

string ele_idfile = "jsons_archive_2023/electron_2023BPix.json.gz";
auto ele_c_set = CorrectionSet::from_file(ele_idfile);
auto ele_SF= ele_c_set->at("Electron-ID-SF");

string ele_ssfile = "jsons_archive_2023/electronSS.json.gz";
string ele_HLTfile = "jsons_archive_2023/electronHlt_2023BPix.json.gz";
string scalecset = "2022Re-recoE+PromptFG_ScaleJSON";
string smearingcset = "2022Re-recoE+PromptFG_SmearingJSON";
auto ele_HLT_c_set = CorrectionSet::from_file(ele_HLTfile);
auto ele_ss_c_set = CorrectionSet::from_file(ele_ssfile);
auto ele_HLT= ele_HLT_c_set->at("Electron-HLT-SF");
auto ele_scale= ele_ss_c_set->at(scalecset);
auto ele_smearing= ele_ss_c_set->at(smearingcset);

string muon_idfile = "jsons_archive_2023/muon_2023BPix.json.gz";
string muon_HLTfile = "jsons_archive_2023/MuTri_2023BPix.json.gz";
string muon_ssfile = "jsons_archive_2023/muon_scalesmearing_Summer23BPix.json.gz";
auto muon_c_set = CorrectionSet::from_file(muon_idfile);
auto muon_HLT_c_set = CorrectionSet::from_file(muon_HLTfile);
auto muon_ss_c_set = CorrectionSet::from_file(muon_ssfile);
auto muon_SF1= muon_c_set->at("NUM_MediumID_DEN_TrackerMuons");
auto muon_SF2= muon_c_set->at("NUM_TightPFIso_DEN_MediumID");
auto muon_HLT_SF= muon_HLT_c_set->at("NUM_IsoMu24_DEN_customloosesel");
auto muon_adata = muon_ss_c_set->at("a_data");
auto muon_amc = muon_ss_c_set->at("a_mc");
auto muon_Mdata = muon_ss_c_set->at("m_data");
auto muon_Mmc = muon_ss_c_set->at("m_mc");
auto muon_cbparams = muon_ss_c_set->at("cb_params");
auto muon_kdata = muon_ss_c_set->at("k_data");
auto muon_kmc = muon_ss_c_set->at("k_mc");
auto muon_polyparams = muon_ss_c_set->at("poly_params");

string tau_idfile = "jsons_archive_2023/tau_DeepTau2018v2p5_2023_postBPix_FIXED.json.gz";
string tau_idfile2 = "jsons_archive_2023/test_tau_pt-dm_2023_postBPix.json.gz";
string tes_file = "jsons_archive_2023/test_tes_tau_pt-dm_2023_postBPix.json.gz";
auto tau_c_set = CorrectionSet::from_file(tau_idfile);
auto tau_c_set2 = CorrectionSet::from_file(tau_idfile2);
auto tes_c_set = CorrectionSet::from_file(tes_file);
auto tau_SFvsjet= tau_c_set->at("DeepTau2018v2p5VSjet");
auto tau_SFvsjet2= tau_c_set2->at("test_DeepTau2018v2p5VSjet_pt-dm");
auto tau_TES= tes_c_set->at("test_DeepTau2018v2p5VSjet_tes_pt-dm");
auto tau_SFvse= tau_c_set->at("DeepTau2018v2p5VSe");
auto tau_SFvsmu= tau_c_set->at("DeepTau2018v2p5VSmu");
auto tau_energyscale= tau_c_set->at("tau_energy_scale");

string DY_ptfile = "jsons_archive_2023/DY_pTll_weights_2023postBPix.json.gz";
auto DY_c_set = CorrectionSet::from_file(DY_ptfile);
auto DY_SF= DY_c_set->at("DY_pTll_reweighting");


//-----------------------------------------------------------------------------------------------------------------------------

struct YearConfig{
    string goldenjson;
    float lumi_recorded;
    decltype(JET_veto_2024) jetvetomap;
    array<float, 10> lumi_eff;
};

// lumi_eff array = {PFjet40,PFJet60,PFjet80,PFJet110,PFjet140,PFJet200,PFjet260,PFHT180,PFHT250,PFHT350}
map<string, YearConfig> YearConfig_map = {
    {"2024", {"jsons/Cert_Collisions2024_378981_386951_Golden.json",
              109.95,
              JET_veto_2024,
              {0.00022,0.00166,0.00640,0.02428,0.07285,0.31221,0.85371,0.00959,0.02665,0.42686}}
    },
    {"2025", {"jsons/Cert_Collisions2025_391658_398903_Golden.json",
              110.63,
              JET_veto_2025,
              {0.00023,0.00152,0.00456,0.02050,0.06149,0.25740,0.86425,0.00851,0.02459,0.43213}}
    }
};

//Practical functions for some computations

Float_t deltaPhi(const Float_t& phi1, const Float_t &phi2) {
  Float_t result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

void JetPhi(Float_t &phi){  // Safeguard against cases where the phi is not in the range of the vetomap corrections
  if(abs(phi)>3.141592653589793){ 
    if(std::signbit(phi)) phi+=6.2831853;
    else phi-=6.2831853;
  }
}

double deltaR(const ROOT::Math::PtEtaPhiMVector &a, const ROOT::Math::PtEtaPhiMVector &b){
  double dPhi = deltaPhi(a.Phi(), b.Phi());
  return sqrt( dPhi*dPhi + (a.Eta() - b.Eta())*(a.Eta() - b.Eta()) );
}

double get_rndm(double mean, double sigma, double n, double alpha, double phi, int evtNumber, int lumiNumber) {
  // instantiate CB and get random number following the CB
  CrystalBall cb(mean, sigma, alpha, n);
  int64_t phi_seed = static_cast<int64_t>((phi / M_PI) * ((1LL << 31) - 1)) & 0xFFF;
  SeedSequence seq{static_cast<uint32_t>(evtNumber), static_cast<uint32_t>(lumiNumber), static_cast<uint32_t>(phi_seed)};
  uint32_t seed;
  seq.generate(&seed, &seed + 1);
  TRandom3 rnd(seed);
  double rndm = rnd.Rndm();
  return cb.invcdf(rndm);
}

// B-tagging efficiencies (needed to apply SFs correctly)

Float_t getBTagEff(Float_t pt, Float_t eta, int flav){

  const Double_t ptlimits[]={20.0,30.0,50.0,70.0,100.0,140.0,200.0,300.0,600.0,1000.0};
  const Double_t etalimits[]={0.0,0.9,1.5,2.1,2.4};
  
  float mapb[9][4]={{0.924,0.915,0.911,0.908},{0.930,0.919,0.913,0.903},{0.937,0.927,0.919,0.909},{0.942,0.931,0.924,0.914},{0.946,0.936,0.929,0.919},{0.948,0.939,0.932,0.922},{0.949,0.939,0.932,0.924},{0.946,0.935,0.933,0.924},{0.948,0.941,0.947,0.980}};;
  float mapc[9][4]={{0.593,0.589,0.606,0.638},{0.551,0.545,0.553,0.570},{0.523,0.523,0.526,0.529},{0.504,0.507,0.510,0.509},{0.492,0.499,0.504,0.505},{0.499,0.504,0.510,0.508},{0.522,0.521,0.533,0.543},{0.566,0.576,0.602,0.631},{0.639,0.655,0.714,0.817}};;
  float maplight[9][4]={{0.184,0.220,0.286,0.384},{0.109,0.138,0.184,0.253},{0.077,0.101,0.136,0.186},{0.062,0.082,0.112,0.157},{0.056,0.073,0.104,0.151},{0.058,0.078,0.114,0.162},{0.074,0.102,0.148,0.214},{0.121,0.171,0.248,0.350},{0.232,0.339,0.471,0.625}};;

  for(int i=0; i<9;i++){
    if(pt<ptlimits[i+1]){
      for(int j=0; j<4;j++){
        if(eta<etalimits[j+1]){
          if(flav==5) return mapb[i][j];
          else if(flav==4) return mapc[i][j];
          else if(flav==0) return maplight[i][j];
          break;}
      }
      break;
    }

  }
  return 0.;
}

//Object selectors

Bool_t TauSelector(Float_t &pt, Float_t eta, UChar_t vse_, UChar_t vsmu_, UChar_t vsjet_, UChar_t source_, UChar_t DM_, Float_t dz, Float_t &weight){
  int DM = static_cast<int>(DM_);
  if(DM==2 or DM==5 or DM==6) return false;
  int vse = static_cast<int>(vse_);
  int vsmu = static_cast<int>(vsmu_);
  int vsjet = static_cast<int>(vsjet_);

  if(pt>20 and abs(eta)<2.3 and vse>=6 and vsmu>=4 and vsjet>=7 and abs(dz)<0.2){
    int source = static_cast<int>(source_);

    //Scale factor for genuine taus
    weight*=tau_SFvsjet2->evaluate({pt,DM,source,"VTight","Tight","default","dm"});

    //Scale factors for misidentified taus
    if(source==2 || source==4) {
      weight*=tau_SFvsmu->evaluate({eta,source,"Tight","nom"});
    }
    else if(source==1 || source==3) {
      weight*=tau_SFvse->evaluate({eta,DM,source,"Tight","nom"});
    }
    //Energy scale correction
    else if(source==5){
      pt*=tau_TES->evaluate({pt,DM,source,"VTight","Tight","default","dm"}); 
    }
    return true;
  }
  else return false;
}

Bool_t ElectronSelector(Float_t &pt, Float_t eta, Float_t phi, Bool_t id, Float_t dxy, Float_t dz, Bool_t convveto, Float_t r9, UChar_t gain, UInt_t run, Float_t &weight){
  if(pt>30 and abs(eta)<2.5 and id and abs(dxy)<0.1 and abs(dz)<0.2 and convveto){
    weight*=ele_HLT->evaluate({"2023PromptD","sf","HLT_SF_Ele30_MVAiso80ID",eta,pt}); //Trigger SF
    weight*=ele_SF->evaluate({"2023PromptD","sf","wp80iso",eta,pt,phi}); //ID SF

    if(pt<75) {
      weight*=ele_SF->evaluate({"2023PromptD","sf","Reco20to75",eta,pt,phi}); //reco SF    }
    else {
      weight*=ele_SF->evaluate({"2023PromptD","sf","RecoAbove75",eta,pt,phi});
    }

    pt*=ele_scale->evaluate({"total_correction",static_cast<int>(gain),static_cast<double>(run),eta,r9,pt}); //Momentum scale correction
    float sig_smear=ele_smearing->evaluate({"rho",eta,r9});                                                                                                             
    float ran=gRandom->Gaus(1.,sig_smear);                                                                                                                                          
    pt*=ran; //Momentum smearing correction
    return true;
  }
  else return false;  
}

Bool_t MuonSelector(Float_t &pt, Float_t eta, Float_t phi, Bool_t id, Float_t dxy, Float_t dz, Float_t isoscore, Int_t charge, UChar_t tracklayers_char, ULong64_t event, UInt_t ls, Float_t &weight){
  if(pt > 30 && abs(eta) < 2.4 && id && abs(dxy)<0.1 && abs(dz)<0.2 && isoscore<0.15){
    weight *= muon_SF1->evaluate({abs(eta),pt,"nominal"}); // ID SF
    weight *= muon_SF2->evaluate({abs(eta),pt,"nominal"}); // ISO SF
    weight *= muon_HLT_SF->evaluate({abs(eta),pt,"nominal"}); // TRIG SF
    float aMC = muon_amc->evaluate({eta,phi,"nom"});
    float MMC = muon_Mmc->evaluate({eta,phi,"nom"});
    pt = 1/((MMC/pt)+aMC*charge); // Momentum scale correction 
    float ntracklayers = static_cast<float>(tracklayers_char);
    float mean = muon_cbparams->evaluate({abs(eta),ntracklayers,0});
    float sigma = muon_cbparams->evaluate({abs(eta),ntracklayers,1});
    float n = muon_cbparams->evaluate({abs(eta),ntracklayers,2});
    float alpha = muon_cbparams->evaluate({abs(eta),ntracklayers,3});
    float kDATA = muon_kdata->evaluate({abs(eta),"nom"});
    float kMC = muon_kmc->evaluate({abs(eta),"nom"});
    float polyparam0 = muon_polyparams->evaluate({abs(eta),ntracklayers,0});
    float polyparam1 = muon_polyparams->evaluate({abs(eta),ntracklayers,1});
    float polyparam2 = muon_polyparams->evaluate({abs(eta),ntracklayers,2});
    float std = polyparam0+pt*polyparam1+pt*pt*polyparam2;
    float kfactor=0;
    if(kDATA>kMC) kfactor=sqrt(kDATA*kDATA-kMC*kMC);
    float rndm = get_rndm(mean, sigma, n, alpha, phi, static_cast<int>(event), ls);
    pt*=(1+kfactor*std*rndm); // Momentum smearing correction
    return true;
  }
  else return false;
}

Bool_t JetSelector(Float_t &pt, Float_t eta, Float_t &phi, Float_t rawfactor, Float_t rhocalo, Float_t neHEF, Float_t neEmEF, Float_t chEmEF, Float_t muEF, Float_t chHEF, Int_t neMultiplicity, Int_t chMultiplicity){
  if(abs(eta)<5.1){
    JetPhi(phi);

    pt*=(1-rawfactor);
    pt*= JEC_L2MC->evaluate({eta,phi,pt}); // Jet pt scale correction
    float ran=gRandom->Gaus(0,JER->evaluate({eta,pt,rhocalo}));
    float JERSF = JER_SF->evaluate({eta,pt,"nom"});
    pt*=1+ran*sqrt(JERSF*JERSF-1);

    bool jetid = jet_id_tight_lepveto->evaluate({abs(eta), chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, chMultiplicity+neMultiplicity});

    if(!jetid) return false;
    if(pt<=50 && abs(eta)>2.5 && abs(eta)<3) return false;    // spikes in that abs(eta) range (for 2025 it should be resolved)    
    
    if(pt>30) return true;
    else return false;
  }
  else return false;
}

map<UInt_t, vector<pair<UInt_t,UInt_t>>> loadGoldenJSON(string filename){
    ifstream f(filename); 
    json data = json::parse(f);
    map<UInt_t, vector<pair<UInt_t,UInt_t>>> goldenMap;

    for (auto& [key, ranges] : data.items()) {
        UInt_t run = stoi(key);
        for (auto& ls : ranges) goldenMap[run].push_back({ls[0], ls[1]});
    }

    return goldenMap;
}

Bool_t is_valid_event(const std::map<UInt_t, std::vector<std::pair<UInt_t, UInt_t>>>& goldenMap, UInt_t run, UInt_t lumi) {
    auto it = goldenMap.find(run);
    if (it == goldenMap.end()) return false; // La run non esiste

    for (const auto& [in, fin] : it->second) {
        if (lumi >= in && lumi <= fin) {
            return true; // L'evento cade in un intervallo valido
        }
    }
    return false;
}

string Run(TString file) {
    string year;
    
    if (file.Contains("Run2022") || file.Contains("RunIII2022"))      year = "2022";
    else if (file.Contains("Run2023") || file.Contains("RunIII2023")) year = "2023";
    else if (file.Contains("Run2024") || file.Contains("RunIII2024")) year = "2024";
    else if (file.Contains("Run2025") || file.Contains("RunIII2025")) year = "2025";
    else if (file.Contains("Run2026") || file.Contains("RunIII2026")) year = "2026";
    else year = "Run year wasn't found!";

    return year;
}

Bool_t isMC(TString file){
  if(file.Contains("Summer2024")) return true;
  else return false;
}