/** Macro to run Atmospheric Systematics at Low and High scales. 
 *  Runs for data and MC to produce final multiplicty plots. 
 *  @authors Jyotirmai Singh, Max Smiley  */

#include <string> 
#include <iostream>
#include <vector>
#include <TH1.h>
#ifdef __MAKECINT__
#pragma link C++ class std::vector<TH1*>;
#endif
void runAll(std::string dataName, std::string mcName) {
  gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine(".x sysStyle.C");
  gROOT->SetStyle("sysStyle");
  gROOT->ProcessLine(".L Systematic.cpp");
  gROOT->ProcessLine(".L HighScaleSystematics.cpp");
  gROOT->ProcessLine(".L LowScaleSystematics.cpp");
  gROOT->ProcessLine(".L HistogramMaker.cpp");
  gROOT->ProcessLine(".L HistogramMerger.cpp");
  gROOT->ProcessLine(".L HistogramOverlayer.cpp");
  std::string dataCommand = ".x loop.cc(true, \"" + dataName + "\")";
  std::cout << dataCommand << std::endl;
  gROOT->ProcessLine(dataCommand.c_str());
  gRandom->SetSeed(42);
  std::string mcCommand = ".x loop.cc(false, \"" +  mcName + "\")";
  gROOT->ProcessLine(mcCommand.c_str());
  gROOT->ProcessLine(".x superimposeDataMC.cpp(\"Systematics.root\", \"SystematicsData.root\")");
}
