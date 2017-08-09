#include "../include/PUWeight.h"

#include "TFile.h"

#include <iostream>

//==============================================================================================
// Get MC pile-up scenario from string representation
PUWeight::Scenario PUWeight::toScenario(const std::string& str) {
  PUWeight::Scenario sc = Summer16_25ns;
  if( str == "PUS25ns" ) sc = Summer16_25ns;
  else {
    std::cerr << "\n\nERROR unknown scenario '" << str << "'" << std::endl;
    throw std::exception();
  }

  return sc;
}

//==============================================================================================
// MC pile-up scenario to string representation
std::string PUWeight::toString(const PUWeight::Scenario sc) {
  std::string str;
  if( sc == Summer16_25ns ) str = "PUS25ns";
  else {
    std::cerr << "\n\nERROR unknown scenario '" << sc << "'" << std::endl;
    throw std::exception();
  }

  return str;
}

//==============================================================================================
// Constructor. Initializes default behaviour to return PU weight of 1
PUWeight::PUWeight()
  : isInit_(false), nPUMax_(0) {}

//==============================================================================================
// Initialise weights for a given MC pile-up scenario. Can only be
// called once.
void PUWeight::initPUWeights(const std::string& dataRootFileName, const std::string& dataRootHistName, const std::string& mcScenario) {

  // if( isInit_ ) {
  //   std::cerr << "\n\nERROR in PUWeight: weights already initialised" << std::endl;
  //   throw std::exception();
  // }

  // Get data distribution from file
  TFile file(dataRootFileName.c_str(), "READ");
  TH1* h = NULL;
  file.GetObject(dataRootHistName.c_str(),h);
  if( h == NULL ) {
    std::cerr << "\n\nERROR in PUWeight: Histogram " << dataRootHistName << " does not exist in file '" << dataRootFileName << "'\n.";
    throw std::exception();
  }
  h->SetDirectory(0);
  file.Close();
  
  PUWeight::Scenario sc = toScenario(mcScenario);

  // Computing weights
  puWeigths_ = generateWeights(sc,h);
  nPUMax_ = puWeigths_.size();

  // Clean up
  delete h;

  isInit_ = true;
}

//==============================================================================================
// Get weight factor dependent on number of added PU interactions
double PUWeight::getPUWeight(const int nPU) const {

  double w = 1.;
  if( isInit_ ) {
    if( nPU >= nPUMax_ ) {
      //std::cerr << "WARNING: Number of PU vertices = " << nPU << " out of histogram binning." << std::endl;
      // In case nPU is out-of data-profile binning,
      // use weight from last bin
      w = puWeigths_.back();
    } else {
      w = puWeigths_.at(nPU);
    }
  }

  return w;
}

//==============================================================================================
// Generate weights for given data PU distribution
// Scenarios from: https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
// Code adapted from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
std::vector<double> PUWeight::generateWeights(const PUWeight::Scenario sc, const TH1* data_npu_estimated) const {

  // Store probabilites for each pu bin
  unsigned int nPUMax = 0;
  double *npuProbs = 0;

  nPUMax = 70;
  double npuSummer16_25ns[70] = {
    0.0,                                                                                                                                                                               
    0.0,                                                                                                                                                                              
    6.71208510924e-05,                                                                                                                                                               
    0.000134241702185,                                                                                                                                                              
    3.35604255462e-05,                                                                                                                                                               
    0.000201362553277,                                                                                                                                                               
    0.000402725106554,                                                                                                                                                               
    0.000671208510924,                                                                                                                                                            
    0.00127529617076,                                                                                                                                                               
    0.00214786723496,                                                                                                                                                               
    0.00543678893848,                                                                                                                                                               
    0.0087928314931,                                                                                                                                                               
    0.0143974225593,                                                                                                                                                               
    0.0209081451153,                                                                                                                                                               
    0.0273517468201,                                                                                                                                                               
    0.035574051079,                                                                                                                                                               
    0.0383931268248,                                                                                                                                                               
    0.044299761721,                                                                                                                                                               
    0.0498707923616,                                                                                                                                                              
    0.0516830553411,                                                                                                                                                               
    0.0561801523643,                                                                                                                                                               
    0.0559452293855,              
    0.0569184817263,              
    0.05560962513,                
    0.0527905493842,              
    0.045776420445,               
    0.0451723327852,              
    0.0383595663993,              
    0.0338289089506,              
    0.0327885357586,              
    0.0266469778837,              
    0.0227539685203,              
    0.0199013323489,              
    0.0185253549015,              
    0.0158405208578,              
    0.014867268517,               
    0.0128200825586,              
    0.0121488740477,              
    0.0107057757492,              
    0.00943047957848,             
    0.00956472128067,             
    0.00721549149243,             
    0.00714837064134,             
    0.00667852468369,             
    0.00550390978958,             
    0.00496694298084,             
    0.00432929489546,             
    0.00332248212907,             
    0.0030539987247,              
    0.00181226297949,             
    0.00228210893714,             
    0.00137597744739,             
    0.00137597744739,             
    0.00127529617076,             
    0.000570527234285,            
    0.000134241702185,            
    0.000234922978823,            
    0.00026848340437,             
    6.71208510924e-05,            
    0.0,                          
    6.71208510924e-05,            
    3.35604255462e-05,            
    0.0,                          
    3.35604255462e-05,            
    0.0,                          
    0.0,                          
    0.0,                          
    3.35604255462e-05,            
    0.0,                          
    0.0                           
  };
  npuProbs = npuSummer16_25ns;

  //}


  // Check that binning of data-profile matches MC scenario
  if( nPUMax != static_cast<unsigned int>(data_npu_estimated->GetNbinsX()) ) {
    std::cerr << "\n\nERROR number of bins (" << data_npu_estimated->GetNbinsX() << ") in data PU-profile does not match number of bins (" << nPUMax << ") in MC scenario " << toString(sc) << std::endl;
    throw std::exception();
  }

  std::vector<double> result(nPUMax,0.);
  double s = 0.;
  for(unsigned int npu = 0; npu < nPUMax; ++npu) {
    const double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
    result[npu] = npu_estimated / npuProbs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(unsigned int npu = 0; npu < nPUMax; ++npu) {
    result[npu] /= s;
  }

  return result;
}
