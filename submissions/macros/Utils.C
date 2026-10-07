#include <iostream>
#include <string>
#include "TString.h"
#include <vector>
#include <array>
#include <map>
#include <fstream>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <cmath>
#include <variant>
using json = nlohmann::json;

// Scale Factors for the analysis, where the year is not specified is for 2024 

//------------- JET VETO MAPS ---------------------------------------------------
string JET_vetomap_2023 = "jsons/jetvetomaps_2023BPix.json.gz";
auto JETveto_c_set_2023 = CorrectionSet::from_file(JET_vetomap_2023);
auto JET_veto_2023      = JETveto_c_set_2023->at("Summer23BPixPrompt23_RunD_V1");

string JET_vetomap_2024 = "jsons/jetvetomaps_2024.json.gz";
auto JETveto_c_set_2024 = CorrectionSet::from_file(JET_vetomap_2024);
auto JET_veto_2024      = JETveto_c_set_2024->at("Summer24Prompt24_RunBCDEFGHI_V1");

string JET_vetomap_2025 = "jsons/jetvetomaps_2025.json.gz";
auto JETveto_c_set_2025 = CorrectionSet::from_file(JET_vetomap_2025);
auto JET_veto_2025      = JETveto_c_set_2025->at("Summer24Prompt25_RunCDEFG_V1");


//------------- Jet SFs ---------------------------------------------------------
string JEC_file           = "jsons/jet_jerc.json.gz";
auto JEC_c_set            = CorrectionSet::from_file(JEC_file);
auto JEC_L2MC             = JEC_c_set->at("Summer24Prompt24_V5_MC_L2Relative_AK4PFPuppi");
auto JER                  = JEC_c_set->at("Summer24Prompt24_JRV2_MC_PtResolution_AK4PFPuppi");
auto JER_SF               = JEC_c_set->at("Summer24Prompt24_JRV2_MC_ScaleFactor_AK4PFPuppi");

string jet_idfile           = "jsons/jetid.json.gz";
auto jet_id_c_set         = CorrectionSet::from_file(jet_idfile);
auto jet_id_tight_lepveto = jet_id_c_set->at("AK4PUPPI_TightLeptonVeto");
auto jet_id_tight         = jet_id_c_set->at("AK4PUPPI_Tight");


//------------- Pile-Up SFs -----------------------------------------------------
string pu_file = "jsons/puWeights_BCDEFGHI.json.gz";
auto pu_c_set  = CorrectionSet::from_file(pu_file);
auto pu_SF     = pu_c_set->at("Collisions24_BCDEFGHI_goldenJSON");


//------------- B-Tagging SFs ---------------------------------------------------
//Non dovrebbero servirmi a nulla poichè veto ogni evento con anche solo un b-jet
string bjet_file     = "jsons/btagging.json.gz";
auto bjet_c_set      = CorrectionSet::from_file(bjet_file);                                                                                                                                                    
auto bcjet_SF        = bjet_c_set->at("UParTAK4_comb");                                                                                                                                                    
auto lightjet_SF     = bjet_c_set->at("UParTAK4_light");      

auto btag_thr_getter = bjet_c_set->at("UParTAK4_wp_values"); 
float btag_thr_      = btag_thr_getter->evaluate({"M"});


//------------- Electron SFs ----------------------------------------------------
string ele_idfile  = "jsons/electron.json.gz";
auto ele_c_set     = CorrectionSet::from_file(ele_idfile);
auto ele_SF        = ele_c_set->at("Electron-ID-SF");

string ele_HLTfile = "jsons/electronHlt.json.gz";
auto ele_HLT_c_set = CorrectionSet::from_file(ele_HLTfile);
auto ele_HLT       = ele_HLT_c_set->at("Electron-HLT-SF");

string ele_ssfile  = "jsons/electronSS_EtDependent.json.gz";
auto ele_ss_c_set  = CorrectionSet::from_file(ele_ssfile);
auto ele_scale     = ele_ss_c_set->at("EGMScale_ElePT_2024");     //???
auto ele_smearing  = ele_ss_c_set->at("EGMSmearAndSyst_ElePT_2024");  //???


//------------- Muon SFs --------------------------------------------------------
string muon_idfile   = "jsons/muon_Z.json.gz";
auto muon_c_set      = CorrectionSet::from_file(muon_idfile);
auto muon_SF1        = muon_c_set->at("NUM_MediumID_DEN_TrackerMuons");
auto muon_SF2        = muon_c_set->at("NUM_TightPFIso_DEN_MediumID");  //??? è effettivamente Tight ???
auto muon_HLT_SF     = muon_c_set->at("NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoTight");
auto muon_SF1den        = muon_c_set->at("NUM_LooseID_DEN_TrackerMuons");
auto muon_SF2den        = muon_c_set->at("NUM_LoosePFIso_DEN_LooseID");

string muon_ssfile   = "jsons/muon_scalesmearing.json.gz";
auto muon_ss_c_set   = CorrectionSet::from_file(muon_ssfile);
auto muon_adata      = muon_ss_c_set->at("a_data");
auto muon_amc        = muon_ss_c_set->at("a_mc");
auto muon_Mdata      = muon_ss_c_set->at("m_data");
auto muon_Mmc        = muon_ss_c_set->at("m_mc");
auto muon_cbparams   = muon_ss_c_set->at("cb_params");
auto muon_kdata      = muon_ss_c_set->at("k_data");
auto muon_kmc        = muon_ss_c_set->at("k_mc");
auto muon_polyparams = muon_ss_c_set->at("poly_params");


//------------- Tau SFs ---------------------------------------------------------
string tau_idfile    = "jsons/tau.json.gz";
auto tau_c_set       = CorrectionSet::from_file(tau_idfile);
auto tau_SFvsjet     = tau_c_set->at("DeepTau2018v2p5VSjet");
auto tau_SFvse       = tau_c_set->at("DeepTau2018v2p5VSe");
auto tau_SFvsmu      = tau_c_set->at("DeepTau2018v2p5VSmu");
auto tau_energyscale = tau_c_set->at("tau_energy_scale");


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


//Object selectors

Bool_t TauSelector(Float_t &pt, Float_t eta, UChar_t vse_, UChar_t vsmu_, UChar_t vsjet_, UChar_t source_, UChar_t DM_, Float_t dz, Float_t &weight){
  int DM = static_cast<int>(DM_);
  if(DM==2 or DM==5 or DM==6) return false;
  int vse = static_cast<int>(vse_);
  int vsmu = static_cast<int>(vsmu_);
  int vsjet = static_cast<int>(vsjet_);

  if(pt>20 and abs(eta)<2.3 and vse>=6 and vsmu>=4 and vsjet>=5 and abs(dz)<0.2){
    int source = static_cast<int>(source_);

    //Scale factor for genuine taus
    weight *= tau_SFvsjet->evaluate({pt,DM,source,"Medium","Tight","nom","dm"});

    //Scale factors for misidentified taus
    if(source==2 || source==4) {
      weight *= tau_SFvsmu->evaluate({abs(eta),source,"Tight","Tight","Medium","nom"});
    }
    else if(source==1 || source==3) {
      weight *= tau_SFvse->evaluate({abs(eta),DM,source,"Tight","nom"});
    }
    //Energy scale correction
    else if(source==5){
      pt *= tau_energyscale->evaluate({pt,abs(eta),DM,source,"DeepTau2018v2p5","Medium","Tight","nom"});  
    }
    return true;
  }
  else return false;
}

Bool_t ElectronSelector(Float_t &pt, Float_t eta, Float_t phi, Bool_t id, Float_t dxy, Float_t dz, Bool_t convveto, Float_t r9, UChar_t gain, UInt_t run, Float_t &weight){
  if(pt>30 and abs(eta)<2.5 and id and abs(dxy)<0.1 and abs(dz)<0.2 and convveto){
    weight *= ele_HLT->evaluate({"2024Prompt","sf","HLT_SF_Ele30_MVAiso80ID",eta,pt}); //Trigger SF
    weight *= ele_SF->evaluate({"2024Prompt","sf","wp80iso",eta,pt}); //ID SF

    if(pt<75) {
      weight *= ele_SF->evaluate({"2024Prompt","sf","Reco20to75",eta,pt}); //reco SF
    }
    else {
      weight *= ele_SF->evaluate({"2024Prompt","sf","RecoAbove75",eta,pt});
    }

    pt *= ele_scale->evaluate({pt,r9,eta}); //Momentum scale correction
    float sig_smear = ele_smearing->evaluate({"smear",pt,r9,eta});                                                                                                             
    float ran = gRandom->Gaus(1.,sig_smear);                                                                                                                                          
    pt *= ran; //Momentum smearing correction
    return true;
  }
  else return false;  
}

Bool_t MuonSelector(Float_t &pt, Float_t eta, Float_t phi, Bool_t id, Float_t dxy, Float_t dz, Float_t isoscore, Int_t charge, UChar_t tracklayers_char, ULong64_t event, UInt_t ls, Float_t &weight){
  if(pt > 30 && abs(eta) < 2.4 && id && abs(dxy)<0.1 && abs(dz)<0.2 && isoscore<0.15){
    weight *= muon_SF1->evaluate({eta,pt,"nominal"}); // ID SF
    weight *= muon_SF2->evaluate({eta,pt,"nominal"}); // ISO SF
    weight *= muon_HLT_SF->evaluate({eta,pt,"nominal"}); // TRIG SF

    float aMC          = muon_amc->evaluate({eta,phi,"nom"});
    float MMC          = muon_Mmc->evaluate({eta,phi,"nom"});
    pt = 1/((MMC/pt)+aMC*charge); // Momentum scale correction 

    float ntracklayers = static_cast<float>(tracklayers_char);
    float mean         = muon_cbparams->evaluate({abs(eta),ntracklayers,0});
    float sigma        = muon_cbparams->evaluate({abs(eta),ntracklayers,1});
    float n            = muon_cbparams->evaluate({abs(eta),ntracklayers,2});
    float alpha        = muon_cbparams->evaluate({abs(eta),ntracklayers,3});
    float kDATA        = muon_kdata->evaluate({abs(eta),"nom"});
    float kMC          = muon_kmc->evaluate({abs(eta),"nom"});
    float polyparam0   = muon_polyparams->evaluate({abs(eta),ntracklayers,0});
    float polyparam1   = muon_polyparams->evaluate({abs(eta),ntracklayers,1});
    float polyparam2   = muon_polyparams->evaluate({abs(eta),ntracklayers,2});
    float std          = polyparam0+pt*polyparam1+pt*pt*polyparam2;
    float kfactor      = 0.;
    if(kDATA>kMC) kfactor = sqrt(kDATA*kDATA-kMC*kMC);
    float rndm         = get_rndm(mean, sigma, n, alpha, phi, static_cast<int>(event), ls);
    pt *= (1+kfactor*std*rndm); // Momentum smearing correction

    return true;
  }
  else return false;
}

Bool_t JetSelector(Float_t &pt, Float_t eta, Float_t &phi, Float_t rawfactor, Float_t rhocalo, Float_t neHEF, Float_t neEmEF, Float_t chEmEF, Float_t muEF, Float_t chHEF, Int_t neMultiplicity, Int_t chMultiplicity){
  if(abs(eta)<5.1){
    JetPhi(phi);

    pt *= (1-rawfactor);
    pt *= JEC_L2MC->evaluate({eta,phi,pt}); // Jet pt scale correction

    float ran = gRandom->Gaus(0,JER->evaluate({eta,pt,rhocalo}));
    float JERSF = JER_SF->evaluate({eta,pt});
    pt *= 1+ran*sqrt(JERSF*JERSF-1);

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

// Prescaled luminosity for the hadronic HLTs
Float_t trigpath_Jet(string year, const Bool_t HLT_PFJet40_, const Bool_t HLT_PFJet60_, const Bool_t HLT_PFJet80_, const Bool_t HLT_PFJet110_, const Bool_t HLT_PFJet140_, const Bool_t HLT_PFJet200_, const Bool_t HLT_PFJet260_, const Bool_t HLT_PFHT180_, const Bool_t HLT_PFHT250_, const Bool_t HLT_PFHT350_){
    vector<float> lumi = {-200.};
    if(HLT_PFJet40_)  lumi.push_back(YearConfig_map[year].lumi_eff[0]);
    if(HLT_PFJet60_)  lumi.push_back(YearConfig_map[year].lumi_eff[1]);
    if(HLT_PFJet80_)  lumi.push_back(YearConfig_map[year].lumi_eff[2]);
    if(HLT_PFJet110_) lumi.push_back(YearConfig_map[year].lumi_eff[3]);
    if(HLT_PFJet140_) lumi.push_back(YearConfig_map[year].lumi_eff[4]);
    if(HLT_PFJet200_) lumi.push_back(YearConfig_map[year].lumi_eff[5]);
    if(HLT_PFJet260_) lumi.push_back(YearConfig_map[year].lumi_eff[6]);
    if(HLT_PFHT180_)  lumi.push_back(YearConfig_map[year].lumi_eff[7]);
    if(HLT_PFHT250_)  lumi.push_back(YearConfig_map[year].lumi_eff[8]);
    if(HLT_PFHT350_)  lumi.push_back(YearConfig_map[year].lumi_eff[9]);

    Float_t max = *max_element(lumi.begin(), lumi.end());

    return max;
}

Bool_t TauSelector_FR(Float_t &pt, Float_t eta, UChar_t vse_, UChar_t vsmu_, UChar_t source_, UChar_t DM_, Float_t dz, Float_t &weight_den, Float_t &weight_num){
  int DM = static_cast<int>(DM_);
  if(DM==2 or DM==5 or DM==6) return false;
  int vse = static_cast<int>(vse_);
  int vsmu = static_cast<int>(vsmu_);
  float energy_SF = 0.;

  if(!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(dz)){
    cout << "TauSelector_FR: ingresso non finito pt=" << pt << " eta=" << eta << " dz=" << dz << endl;
  }
  
  if(pt>20 and abs(eta)<2.3 and vse>=6 and vsmu>=4 and abs(dz)<0.2){
    int source = static_cast<int>(source_);

    //Scale factor for genuine taus
    weight_den *= tau_SFvsjet->evaluate({pt,DM,source,"Loose","Tight","nom","dm"});
    weight_num *= tau_SFvsjet->evaluate({pt,DM,source,"Medium","Tight","nom","dm"});

    //Scale factors for misidentified taus
    if(source==2 || source==4) {
      weight_den *= tau_SFvsmu->evaluate({abs(eta),source,"Tight","Tight","Loose","nom"});
      weight_num *= tau_SFvsmu->evaluate({abs(eta),source,"Tight","Tight","Medium","nom"});
    }
    else if(source==1 || source==3) {
      weight_den *= tau_SFvse->evaluate({abs(eta),DM,source,"Tight","nom"});
      weight_num *= tau_SFvse->evaluate({abs(eta),DM,source,"Tight","nom"});
    }
    //Energy scale correction
    else if(source==5){
      energy_SF = tau_energyscale->evaluate({pt,abs(eta),DM,source,"DeepTau2018v2p5","Loose","Tight","nom"});
      pt *= energy_SF;
    }

    if(!std::isfinite(pt)){
    cout << "TauSelector_FR: uscita non finito pt =" << pt << endl;
    cout << "Energy SF = " << energy_SF << endl;      
    }

    return true;
  }
  else return false;
}

Bool_t ElectronSelector_FR(Float_t &pt, Float_t eta, Float_t phi, Float_t dxy, Float_t dz, Bool_t convveto, Float_t r9, UChar_t gain, UInt_t run, Float_t &weight_den, Float_t &weight_num){  
  if(pt>30 and abs(eta)<2.5 and abs(dxy)<0.1 and abs(dz)<0.2 and convveto){
    weight_den *= ele_SF->evaluate({"2024Prompt","sf","wp90iso",eta,pt}); //ID SF
    weight_num *= ele_SF->evaluate({"2024Prompt","sf","wp80iso",eta,pt}); //ID SF

    if(!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(dz)){
      cout << "TauSelector_FR: ingresso non finito pt=" << pt << " eta=" << eta << " dz=" << dz << endl;
    }

    if(pt<75) {
      weight_den *= ele_SF->evaluate({"2024Prompt","sf","Reco20to75",eta,pt}); //reco SF
      weight_num *= ele_SF->evaluate({"2024Prompt","sf","Reco20to75",eta,pt}); //reco SF
    }
    else {
      weight_den *= ele_SF->evaluate({"2024Prompt","sf","RecoAbove75",eta,pt});
      weight_num *= ele_SF->evaluate({"2024Prompt","sf","RecoAbove75",eta,pt});
    }

    pt *= ele_scale->evaluate({pt,r9,eta}); //Momentum scale correction
    float sig_smear = ele_smearing->evaluate({"smear",pt,r9,eta});                                                                                                             
    float ran = gRandom->Gaus(1.,sig_smear);                                                                                                                                          
    pt *= ran; //Momentum smearing correction

    if(!std::isfinite(pt)){
      cout << "TauSelector_FR: uscita non finito pt =" << pt << endl;
      cout << "Energy SF = " << ran << endl;      
    }

    return true;
  }
  else return false;   
}

Bool_t MuonSelector_FR(Float_t &pt, Float_t eta, Float_t phi, Float_t dxy, Float_t dz, Int_t charge, UChar_t tracklayers_char, ULong64_t event, UInt_t ls, Float_t &weight_den, Float_t &weight_num){
  if(pt > 30 && abs(eta) < 2.4 && abs(dxy)<0.1 && abs(dz)<0.2){
    weight_den *= muon_SF1den->evaluate({eta,pt,"nominal"}); // ID SF
    weight_den *= muon_SF2den->evaluate({eta,pt,"nominal"}); // ISO SF
    weight_num *= muon_SF1->evaluate({eta,pt,"nominal"}); // ID SF
    weight_num *= muon_SF2->evaluate({eta,pt,"nominal"}); // ISO SF

    if(!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(dz)){
      cout << "TauSelector_FR: ingresso non finito pt=" << pt << " eta=" << eta << " dz=" << dz << endl;
    }

    float aMC          = muon_amc->evaluate({eta,phi,"nom"});
    float MMC          = muon_Mmc->evaluate({eta,phi,"nom"});
    pt = 1/((MMC/pt)+aMC*charge); // Momentum scale correction 

    float ntracklayers = static_cast<float>(tracklayers_char);
    float mean         = muon_cbparams->evaluate({abs(eta),ntracklayers,0});
    float sigma        = muon_cbparams->evaluate({abs(eta),ntracklayers,1});
    float n            = muon_cbparams->evaluate({abs(eta),ntracklayers,2});
    float alpha        = muon_cbparams->evaluate({abs(eta),ntracklayers,3});
    float kDATA        = muon_kdata->evaluate({abs(eta),"nom"});
    float kMC          = muon_kmc->evaluate({abs(eta),"nom"});
    float polyparam0   = muon_polyparams->evaluate({abs(eta),ntracklayers,0});
    float polyparam1   = muon_polyparams->evaluate({abs(eta),ntracklayers,1});
    float polyparam2   = muon_polyparams->evaluate({abs(eta),ntracklayers,2});
    float std          = polyparam0+pt*polyparam1+pt*pt*polyparam2;
    float kfactor      = 0.;
    if(kDATA>kMC) kfactor = sqrt(kDATA*kDATA-kMC*kMC);
    float rndm         = get_rndm(mean, sigma, n, alpha, phi, static_cast<int>(event), ls);
    pt *= (1+kfactor*std*rndm); // Momentum smearing correction

    if(!std::isfinite(pt)){
      cout << "TauSelector_FR: uscita non finito pt =" << pt << endl;
      cout << "Energy SF = " << (1+kfactor*std*rndm) << endl;      
    }

    return true;
  }
  else return false;
}

Float_t m_T(Float_t leppt, Float_t metpt, Float_t lepphi, Float_t metphi){
    Float_t deltaphi = deltaPhi(lepphi,metphi);

    return TMath::Sqrt(2*leppt*metpt*(1 - TMath::Cos(deltaphi)));
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

Float_t getBTagEff(Float_t pt, Float_t eta, int flav){

  const Double_t ptlimits[]={20.0,30.0,50.0,70.0,100.0,140.0,200.0,300.0,600.0,1000.0};
  const Double_t etalimits[]={0.0,0.9,1.5,2.1,2.4};
  
  float mapb[9][4]={{0.783,0.765,0.742,0.710},{0.837,0.817,0.798,0.771},{0.869,0.849,0.834,0.811},{0.884,0.865,0.853,0.830},{0.891,0.874,0.863,0.844},{0.896,0.881,0.868,0.848},{0.893,0.878,0.861,0.837},{0.886,0.867,0.843,0.813},{0.878,0.853,0.824,0.766}};
  float mapc[9][4]={{0.119,0.136,0.144,0.154},{0.085,0.101,0.111,0.125},{0.072,0.085,0.095,0.107},{0.067,0.079,0.087,0.099},{0.064,0.076,0.083,0.098},{0.065,0.078,0.086,0.103},{0.075,0.090,0.098,0.117},{0.106,0.127,0.137,0.154},{0.173,0.194,0.202,0.189}};
  float maplight[9][4]={{0.013,0.016,0.021,0.028},{0.008,0.010,0.014,0.019},{0.006,0.008,0.011,0.014},{0.005,0.007,0.009,0.012},{0.005,0.006,0.009,0.012},{0.005,0.006,0.009,0.012},{0.005,0.007,0.010,0.015},{0.008,0.013,0.017,0.023},{0.019,0.028,0.033,0.038}};

  for(int i=0; i<9;i++){
    if(pt < ptlimits[i+1]){
      for(int j=0; j<4;j++){
        if(eta < etalimits[j+1]){
          if(flav==5) return mapb[i][j];
          else if(flav==4) return mapc[i][j];
          else if(flav==0) return maplight[i][j];
          break;
        }
      }
      break;
    }

  }
  return 0.;
}

void bvetoSelector(Float_t jet_btag, UChar_t flavour_, Float_t pt, Float_t eta, Float_t &weight, Bool_t &btagflag){
  Int_t flavour = static_cast<int>(flavour_);
  Float_t eff = 0.;
  Float_t SF = 0.;

  eff = getBTagEff(pt, abs(eta), flavour);
  if(eff == 0.){
    cout<< "Error: Unknown jet flavour " << flavour << "|" << " or pt out of bounds " << pt << endl;
    btagflag = true;
    return;
  }

  if(jet_btag >= btag_thr_){
    btagflag = true;
    return;
  }else{
    if(flavour == 5 || flavour == 4){
      SF = bcjet_SF->evaluate({"central","M",flavour,abs(eta),pt});
      weight *= (1 - SF*eff)/(1 - eff);
    }else if(flavour == 0){
      SF = lightjet_SF->evaluate({"central","M",flavour,abs(eta),pt});
      weight *= (1 - SF*eff)/(1 - eff);
    }
    return;
  }
}