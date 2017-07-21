/** @class HistogramOverlayer
 *  Takes a data file and MC file (or any two files) and
 *  overlays the 1D or profile histograms within, 
 *  and writes out any 2D histograms as before
 *  @author Max Smiley */

#include <cmath>
#include "HistogramOverlayer.h"
#include "HistogramMaker.h"

HistogramOverlayer::HistogramOverlayer(TFile* dataFile, TFile* mcFile, TFile* outFile)
{
    _dataFile = dataFile;
    _mcFile = mcFile;
    _outFile = outFile;

    std::string outfilename = std::string(_outFile->GetName());
    size_t index = outfilename.find(".");
    _pdfFilename = outfilename.substr(0, index) + ".pdf";

}

void renameHist(TH1* hist, std::string addon) {
    hist->SetName((string(hist->GetName())+addon).c_str());
    hist->SetTitle((string(hist->GetTitle()) + string(" ") + addon).c_str());
}

void logCheck(TCanvas* canv, TH1* histData, TH1* histMC) {
    if ((string(histData->GetName())).find("nfollowers") != std::string::npos &&
        (string(histData->GetName())).find("eeffenergy") != std::string::npos) {
      canv->SetLogx();
    }
    else {
      canv->SetLogx(0);
    }
}

void setDrawStyles(TCanvas* canv, TH1F* histData, TH1F* histMC) {
    std::string styleMC;
    std::string styleData = "e1";
    histMC->SetMarkerColor(kRed);
    histMC->SetLineColor(kRed);

    if ((string(histData->GetName())).find("nfollowers") != std::string::npos &&
        (string(histData->GetName())).find("eeffenergy") != std::string::npos) {
      histMC->SetFillColor(kRed);
      histMC->SetFillStyle(3001);
      styleMC = "e2";
    }
    else {
      histData->Scale(1./histData->Integral(), "width");
      histMC->Scale(1./histMC->Integral(), "width");
      styleMC = "hist";
    }

    if (histMC->GetMaximum() >= histData->GetMaximum()) {
      boxMover(canv, histMC, histData, styleMC, styleData);
    }
    else {
      boxMover(canv, histData, histMC, styleData, styleMC);
    }
}

void boxMover(TCanvas* canv, TH1F* hist1, TH1F* hist2, std::string style1, std::string style2) {
        TPaveStats *st1;
        TPaveStats *st2;
        canv->cd();

        hist1->Draw(style1.c_str());
        canv->Modified();
        canv->Update();
        st1 = (TPaveStats*) canv->GetPrimitive("stats");
        st1->SetName((string(st1->GetName())+string(hist1->GetName())).c_str());
        hist2->Draw((style2 +string("sames")).c_str());
        canv->Modified();
        canv->Update();
        st2 = (TPaveStats*) canv->GetPrimitive("stats");
        st2->SetName((string(st2->GetName())+string(hist2->GetName())).c_str());

        Double_t newY2 = st1->GetY1NDC();
        Double_t newY1 = newY2 - (st2->GetY2NDC() - st2->GetY1NDC());
        st2->SetY1NDC(newY1);
        st2->SetY2NDC(newY2);
        canv->Modified();
        canv->Update();

        TLegend* leg = canv->BuildLegend(st2->GetX1NDC(), st2->GetY1NDC() - (st2->GetY2NDC() - st2->GetY1NDC()),st1->GetX2NDC(), st2->GetY1NDC());
        canv->Modified();
        canv->Update();
}

void HistogramOverlayer::printCanvas(TCanvas* canv, int index, int maxsize, bool first, bool oftwo) {
    _outFile->WriteTObject(canv);
    if (index == 0) {
      if (first) canv->Print((_pdfFilename+string("(")).c_str(), "pdf");
      else canv->Print(_pdfFilename.c_str(), "pdf");
    }
    else if (index == maxsize - 1) {
      if (first && oftwo) canv->Print(_pdfFilename.c_str(), "pdf");
      else canv->Print((_pdfFilename+string(")")).c_str(), "pdf");
    }
    else {
        canv->Print(_pdfFilename.c_str(), "pdf");
    }
    canv->Clear();
}

/** Main method, makes all histograms and draws them on
 *  canvases. Saves canvases to _outFile for later access. Applies */
void HistogramOverlayer::makeHistograms()
{
    HistogramMaker* histsData = new HistogramMaker(_dataFile);
    HistogramMaker* histsMC = new HistogramMaker(_mcFile);

    TCanvas* c1 = new TCanvas();
    TCanvas* c2 = new TCanvas();

    if (histsData->hists.size() != histsMC->hists.size()) return;
    size_t numHists = histsData->hists.size();

    for (unsigned i =0; i < numHists; i++) {
      if (histsData->hists.at(i)->IsA() == TH1F::Class() && histsMC->hists.at(i)->IsA() == TH1F::Class()) {
        TH1F* hData = (TH1F*) histsData->hists.at(i);
        TH1F* hMC = (TH1F*) histsMC->hists.at(i);
        hMC->Sumw2();
        hData->Sumw2();

        logCheck(c1, hData, hMC);
        c1->SetName(hData->GetName());
        renameHist(hData, "Data");
        renameHist(hMC, "MC");

        c1->cd();
        setDrawStyles(c1, hData, hMC);
        this->printCanvas(c1, i, numHists, true, false);
      }

      else if (histsData->hists.at(i)->IsA() == TH2F::Class() && histsMC->hists.at(i)->IsA() == TH2F::Class()) {
        TH2F* h2Data = (TH2F*) histsData->hists.at(i);
        TH2F* h2MC = (TH2F*) histsMC->hists.at(i);
        c1->SetName(h2Data->GetName());
        c2->SetName(h2MC->GetName());
        renameHist(h2Data, "Data");
        renameHist(h2MC, "MC");
        logCheck(c1, h2Data, h2MC);
        logCheck(c2, h2Data, h2MC);

        c1->cd();
        h2Data->Draw("colz");
        c2->cd();
        h2MC->Draw("colz");
        this->printCanvas(c1, i, numHists, true, true);
        this->printCanvas(c2, i, numHists, false, true);
      }

      else if (i == numHists - 1) {
         c1->Print((_pdfFilename+string("]")).c_str(), "pdf");
      }
    }

    _outFile->Close();
}
