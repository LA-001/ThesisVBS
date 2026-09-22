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

auto jet_idfile           = "jsons/jetid.json.gz";
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
auto bjet_SF         = bjet_c_set->at("UParTAK4_comb");                                                                                                                                                    
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

  if(pt>20 and abs(eta)<2.3 and vse>=6 and vsmu>=4 and vsjet>=5 and abs(dz)<0.2){
    int source = static_cast<int>(source_);

    //Scale factor for genuine taus
    if(pt <= 140.) weight *= tau_SFvsjet->evaluate({pt,DM,source,"Medium","Tight","nom","dm"});
    else if(pt > 140.) weight *= tau_SFvsjet->evaluate({pt,DM,source,"Medium","Tight","nom","pt"});

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

Bool_t JetSelector(string year, Float_t &pt, Float_t eta, Float_t &phi, Float_t rawfactor, Float_t rhocalo,
                  Float_t neHEF, Float_t neEmEF, Float_t chEmEF, Float_t muEF, Float_t chHEF, Int_t neMultiplicity, Int_t chMultiplicity){
  if(abs(eta)<5.1){
    if(abs(phi)>3.141592653589793){ // Safeguard against cases where the phi is not in the range of the vetomap corrections
      if(std::signbit(phi)) phi+=6.2831853;
      else phi-=6.2831853;
    }
    bool jetveto = YearConfig_map[year].jetvetomap->evaluate({"jetvetomap",eta,phi});
    if(jetveto) return false;

    pt *= (1-rawfactor);
    pt *= JEC_L2MC->evaluate({eta,phi,pt}); // Jet pt scale correction

    float ran = gRandom->Gaus(0,JER->evaluate({eta,pt,rhocalo}));
    float JERSF = JER_SF->evaluate({eta,pt});
    pt *= 1+ran*sqrt(JERSF*JERSF-1);

    if(pt<=50 && abs(eta)>2.5 && abs(eta)<3) return false;    // spikes in that abs(eta) range (for 2025 it should be resolved)

    bool jetid = jet_id_tight_lepveto->evaluate({abs(eta), chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, chMultiplicity+neMultiplicity});
    if(!jetid) return false;

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

/*    // Ancora da guardare questi
Bool_t TauSelector_prompt(Float_t &pt, Float_t eta, UChar_t vse_, UChar_t vsmu_, UChar_t source_, UChar_t DM_, Float_t dz, Float_t &weight){
  
  int DM = static_cast<int>(DM_);
  if(DM==2 or DM==5 or DM==6) return false;
  int vse = static_cast<int>(vse_);
  int vsmu = static_cast<int>(vsmu_);
  
  if(pt>20 and abs(eta)<2.3 and vse>=6 and vsmu>=4 and abs(dz)<0.2){
    int source = static_cast<int>(source_);
    //Scale factor for genuine taus
    weight*=tau_SFvsjet3->evaluate({pt,DM,source,"Loose","Tight","nom","dm"});            //wp Loose per lo studio dei leptoni fake, poi rimettere Medium
    //Scale factors for misidentified taus
    if(source==2 || source==4) {
      weight*=tau_SFvsmu3->evaluate({abs(eta),source,"Tight","Tight","Loose","nom"});
    }
    else if(source==1 || source==3) {
      weight*=tau_SFvse3->evaluate({abs(eta),DM,source,"Tight","nom"});
    }
    //Energy scale correction
    else if(source==5){
      pt*=tau_energyscale3->evaluate({pt,abs(eta),DM,source,"DeepTau2018v2p5","Loose","Tight","nom"});
    }
    return true;
  }
  else return false;
}

Bool_t ElectronSelector_prompt(Float_t &pt, Float_t eta, Float_t phi, Float_t dxy, Float_t dz, Bool_t convveto, Float_t r9, UChar_t gain, UInt_t run, Float_t &weight){
  if(pt>30 and abs(eta)<2.5 and abs(dxy)<0.1 and abs(dz)<0.2 and convveto){
    if(pt<75) {
      weight*=ele_SF->evaluate({"2023PromptD","sf","Reco20to75",eta,pt,phi}); //reco SF
    }
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

void Electron_weightSF(Float_t &weight, Float_t eta, Float_t pt, Float_t phi, Bool_t iswp80){
    if(iswp80){
        weight *= ele_SF->evaluate({"2023PromptD","sf","wp80iso",eta,pt,phi});
    }else{
        weight *= ele_SF->evaluate({"2023PromptD","sf","wp90iso",eta,pt,phi});
    }
}

Bool_t MuonSelector_prompt(Float_t &pt, Float_t eta, Float_t phi, Float_t dxy, Float_t dz, Int_t charge, UChar_t tracklayers_char, ULong64_t event, UInt_t ls, Float_t &weight){
  if(pt > 30 && abs(eta) < 2.4 && abs(dxy)<0.1 && abs(dz)<0.2){
    weight *= muon_SF1->evaluate({abs(eta),pt,"nominal"}); // ID SF
    weight *= muon_SF2->evaluate({abs(eta),pt,"nominal"}); // ISO SF
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
*/

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
