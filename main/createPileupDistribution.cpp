#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "interface/PUReweighting.h"
#include "interface/OutTree.h"
#include "interface/StdHistoSet.h"
#include "src/RoccoR.cc"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TRandom3.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> createPileupDistribution.cpp::usage:   " << argv[0] << " configFileName   [default=0/debug=1]" << std::endl;
    return -1;
  }
  
  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  //--- open files and get trees
  std::string inputFileListPU = opts.GetOpt<std::string>("Input.inputFileListPU");
  
  //--- define histograms
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  TFile* outFile = TFile::Open(Form("%s.root",outputFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  
  //----------------------
  //--- pileup reweighting
  TH1F* pileup_mc = NULL;
  std::string fileName;
  std::ifstream listPU(inputFileListPU.c_str(),std::ios::in);
  while(1)
  {
    getline(listPU,fileName,'\n');
    if( !listPU.good() ) break;
    
    std::cout << fileName << std::endl;
    TFile* tempInFile = TFile::Open((fileName).c_str(),"READ");
    TH1F* tempHisto = (TH1F*)( tempInFile->Get("DumpPU/pileup") );
    
    if( pileup_mc == NULL )
    {
      outFile -> cd();
      pileup_mc = new TH1F(Form("pileup_mc"),"",tempHisto->GetNbinsX(),tempHisto->GetBinLowEdge(1),tempHisto->GetBinLowEdge(tempHisto->GetNbinsX()+tempHisto->GetBinWidth(tempHisto->GetNbinsX())));
    }
    pileup_mc -> Add(tempHisto);
    
    tempInFile -> Close();
  }
  
  outFile -> cd();
  pileup_mc -> Write();
  outFile -> Close();
  
  return 0;
}
