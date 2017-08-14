/** @class HistogramMerger
 *  Takes a nominal, lower, and upper file and
 *  merges the histograms as appropriate, changing errors
 *  where needed.
 *  @authors Jyotirmai Singh, Max Smiley */

#include <cmath>
#include "HistogramMerger.h"
#include "HistogramMaker.h"

HistogramMerger::HistogramMerger(TFile* nomFile, TFile* upFile, TFile* lowFile, TFile* outFile, bool useData)
{
    _nomFile = nomFile;
    _upFile = upFile;
    _lowFile = lowFile;
    _outFile = outFile;
    _data = useData;
}

/** Applies bin error correction to a nominal histogram, changing its bin errors
 *  to abs(upperError - lowerError). If this quantity is 0, we don't change the
 *  bin error since the SetBinError function behaves strangely when the error is
 *  0 and the error bars are not drawn properly.
 *  @param nominal The nominal histogram whose errors are changed.
 *  @param upper The histogram resulting from applying upper systematics.
 *  @param lower The histogram resulting from applying lower systematics. */
void HistogramMerger::applyBinErrorCorrection(TH1F* nominal, TH1F* upper, TH1F* lower)
{
    int nBins = nominal->GetSize();
    for (int i = 0; i < nBins; i++) {
        double upperErr = upper->GetBinContent(i);
        double lowerErr = lower->GetBinContent(i);
        double error = fabs(upperErr - lowerErr);
	if (error == 0.0) {
	    continue;
	}
        nominal->SetBinError(i, error);
    }
}

/** Overlays nominal with statistical and systematic errors.
 *  @param nominal The nominal histogram to be overlaid.
 *  @param upper The histogram resulting from applying upper systematics.
 *  @param lower The histogram resulting from applying lower systematics.
 *  @param canv The canvas on which the overlaid nominal plots are drawn. */
void HistogramMerger::mergeHistograms(TH1F* nominal, TH1F* upper, TH1F* lower,
                                      TCanvas* canv)
{
    TCanvas* canvas = canv;
    TH1F* clonedHisto = nominal->Clone();
    clonedHisto->Draw("e1");
    this->applyBinErrorCorrection(nominal, upper, lower);
    nominal->SetFillColor(kRed);
    nominal->SetFillStyle(3001);
    nominal->Draw("e2 same");
    canv->Update();
}

/** Overlays nominal with statistical and systematic errors. Specialised
 *  for mean histograms where a division by a norm histogram is required.
 *  Performs division by norm first and then calls mergeHistograms.
 *  @param nominal The nominal histogram to be overlaid.
 *  @param upper The histogram resulting from applying upper systematics.
 *  @param lower The histogram resulting from applying lower systematics.
 *  @param nominalNorm What nominal is divided by.
 *  @param upperNorm What upper is divided by.
 *  @param lowerNorm What lower is divided by.
 *  @param canv The canvas on which the overlaid nominal plots are drawn. */
void HistogramMerger::mergeNormDivisionHistograms(TH1F* nominal, TH1F* upper, TH1F* lower,
                                                  TH1F* nominalNorm, TH1F* upperNorm, TH1F* lowerNorm,
                                                  TCanvas* canv)
{
    nominal->Sumw2();
    nominalNorm->Sumw2();
    upper->Sumw2();
    upperNorm->Sumw2();
    lower->Sumw2();
    lowerNorm->Sumw2();

    nominal->Divide(nominalNorm);
    upper->Divide(upperNorm);
    lower->Divide(lowerNorm);

    this->mergeHistograms(nominal, upper, lower, canv);
}

void HistogramMerger::plotPromptHistogram(TH1F* nominal, TH1F* upper, TH1F* lower, TCanvas* canv)
{

    if (_data) {
        TCanvas* c = canv;
        nominal->Draw("e1");
    } else {
        this->mergeHistograms(nominal, upper, lower, canv);
    }
}

void HistogramMerger::plotMeanFollowerHistogram(TH1F *nominal, TH1F *upper, TH1F *lower, TH1F *nominalNorm,
                                                TH1F *upperNorm, TH1F *lowerNorm, TCanvas *canv)
{
    if (_data) {
        TCanvas* c = canv;
        nominal->Sumw2();
        nominalNorm->Sumw2();
        nominal->Divide(nominalNorm);
        nominal->Draw("e1");
    } else {
        this->mergeNormDivisionHistograms(nominal, upper, lower,
                                          nominalNorm, upperNorm, lowerNorm,
                                          canv);
    }
}


/** Main method, makes all histograms and draws them on
 *  canvases. Saves canvases to _outFile for later access. Applies
 *  merging of histograms when necessary. */
void HistogramMerger::makeHistograms()
{
    HistogramMaker* histsNominal = new HistogramMaker(_nomFile);
    HistogramMaker* histsUpper = new HistogramMaker(_upFile);
    HistogramMaker* histsLower = new HistogramMaker(_lowFile);

    TCanvas* cfitpos = new TCanvas("cfitpos", "Fitted Position for Prompt", 1600, 800);
    cfitpos->Divide(3,2);
    cfitpos->cd(1);
    histsNominal->hfitpos[0]->Draw();
    cfitpos->cd(2);
    histsNominal->hfitpos[1]->Draw();
    cfitpos->cd(3);
    histsNominal->hfitpos[2]->Draw();
    cfitpos->cd(4);
    histsNominal->hfitposdiff[0]->Draw();
    cfitpos->cd(5);
    histsNominal->hfitposdiff[1]->Draw();
    cfitpos->cd(6);
    histsNominal->hfitposdiff[2]->Draw();

    TCanvas* cseedpos = new TCanvas("cseedpos", "Seed Position for Prompt", 1600, 800);
    cseedpos->Divide(3,2);
    cseedpos->cd(1);
    histsNominal->hseedpos[0]->Draw();
    cseedpos->cd(2);
    histsNominal->hseedpos[1]->Draw();
    cseedpos->cd(3);
    histsNominal->hseedpos[2]->Draw();
    cseedpos->cd(4);
    histsNominal->hseedposdiff[0]->Draw();
    cseedpos->cd(5);
    histsNominal->hseedposdiff[1]->Draw();
    cseedpos->cd(6);
    histsNominal->hseedposdiff[2]->Draw(); 


    TCanvas* cprompt = new TCanvas("cprompt","Prompt Distributions", 1600, 800);
    cprompt->Divide(3,2);
    cprompt->cd(1);
    histsNominal->hprompt_nhits->Draw("e1");
    cprompt->cd(2);
    histsNominal->hprompt_nrings->Draw("e1");
    cprompt->cd(3);
    histsNominal->hprompt_pid->Draw("e1");
    cprompt->cd(4);
    plotPromptHistogram(histsNominal->hprompt_eeffenergy,
                        histsUpper->hprompt_eeffenergy,
                        histsLower->hprompt_eeffenergy,
                        cprompt);

    cprompt->cd(5);
    plotPromptHistogram(histsNominal->hprompt_ueffenergy,
                        histsUpper->hprompt_ueffenergy,
                        histsLower->hprompt_ueffenergy,
                        cprompt);
    cprompt->cd(6);
    plotPromptHistogram(histsNominal->hprompt_meffenergy,
                        histsUpper->hprompt_meffenergy,
                        histsLower->hprompt_meffenergy,
                        cprompt);
    cprompt->Update();

    TCanvas* cmichele = new TCanvas("cmichele","Michel-e Distributions",900,600);
    cmichele->Divide(3,2);
    cmichele->cd(1);
    histsNominal->hmichele_nhits->Draw("e1");
    cmichele->cd(2);
    histsNominal->hmichele_deltat->Draw("e1");
    cmichele->cd(4);
    histsNominal->hmichele_energy0->Draw("e1");
    cmichele->cd(5);
    histsNominal->hmichele_energy1->Draw("e1");
    cmichele->cd(6);
    histsNominal->hmichele_energy2->Draw("e1");


    TCanvas* cfollowers = new TCanvas("cfollowers","Neutron Followers Distributions",1200,600);
    cfollowers->Divide(4,2);
    cfollowers->cd(1);
    histsNominal->hfollowers_nhits->Draw("e1");
    cfollowers->cd(2);
    histsNominal->hfollowers_energy0->Draw("e1");
    cfollowers->cd(3);
    histsNominal->hfollowers_deltat->Draw("e1");
    cfollowers->cd(4);
    histsNominal->hfollowers_dist->Draw("e1");
    cfollowers->cd(5);
    histsNominal->hfollowers_energy0->Draw("e1");
    cfollowers->cd(6);
    histsNominal->hfollowers_energy1->Draw("e1");
    cfollowers->cd(7);
    histsNominal->hfollowers_energy2->Draw("e1");

    TCanvas* cnfollow = new TCanvas("cnfollow","Number of Followers",900,300);
    cnfollow->Divide(3,1);
    cnfollow->cd(1);
    histsNominal->nfollowers_tot->Draw("e1");
    cnfollow->cd(2);
    histsNominal->nfollowers_sring->Draw("e1");
    cnfollow->cd(3);
    histsNominal->nfollowers_mring->Draw("e1");

    TCanvas* cnhit_nofollow = new TCanvas("cnhit_nofollow","NHit for Prompt Events with NO Followers",900,300);
    cnhit_nofollow->Divide(3,1);
    cnhit_nofollow->cd(1);
    histsNominal->nhit_nofollow_tot->Draw("e1");
    cnhit_nofollow->cd(2);
    histsNominal->nhit_nofollow_sring->Draw("e1");
    cnhit_nofollow->cd(3);
    histsNominal->nhit_nofollow_mring->Draw("e1");

    TCanvas* cnfollowmean_eeffenergy = new TCanvas("cnfollowmean_eeffenergy","Number of Followers vs Prompt Energy",700,900);
    gStyle->SetPalette(56);
    cnfollowmean_eeffenergy->Divide(2,3);
    cnfollowmean_eeffenergy->cd(1);
    cnfollowmean_eeffenergy->cd(1)->SetLogx();
    this->plotMeanFollowerHistogram(histsNominal->nfollowersmean_eeffenergy,
                                    histsUpper->nfollowersmean_eeffenergy,
                                    histsLower->nfollowersmean_eeffenergy,
                                    histsNominal->nfollowersmean_eeffenergy_norm,
                                    histsUpper->nfollowersmean_eeffenergy_norm,
                                    histsLower->nfollowersmean_eeffenergy_norm,
                                    cnfollowmean_eeffenergy);

    cnfollowmean_eeffenergy->cd(2);
    cnfollowmean_eeffenergy->cd(2)->SetLogx();
    histsNominal->nfollowers_eeffenergy->Draw("colz");
    cnfollowmean_eeffenergy->cd(3);
    cnfollowmean_eeffenergy->cd(3)->SetLogx();
    this->plotMeanFollowerHistogram(histsNominal->nfollowersmean_sr_eeffenergy,
                                    histsUpper->nfollowersmean_sr_eeffenergy,
                                    histsLower->nfollowersmean_sr_eeffenergy,
                                    histsNominal->nfollowersmean_sr_eeffenergy_norm,
                                    histsUpper->nfollowersmean_sr_eeffenergy_norm,
                                    histsLower->nfollowersmean_sr_eeffenergy_norm,
                                    cnfollowmean_eeffenergy);
    cnfollowmean_eeffenergy->cd(4);
    cnfollowmean_eeffenergy->cd(4)->SetLogx();
    histsNominal->nfollowers_sr_eeffenergy->Draw("colz");
    cnfollowmean_eeffenergy->cd(5);
    cnfollowmean_eeffenergy->cd(5)->SetLogx();
    this->plotMeanFollowerHistogram(histsNominal->nfollowersmean_mr_eeffenergy,
                                    histsUpper->nfollowersmean_mr_eeffenergy,
                                    histsLower->nfollowersmean_mr_eeffenergy,
                                    histsNominal->nfollowersmean_mr_eeffenergy_norm,
                                    histsUpper->nfollowersmean_mr_eeffenergy_norm,
                                    histsLower->nfollowersmean_mr_eeffenergy_norm,
                                    cnfollowmean_eeffenergy);
    cnfollowmean_eeffenergy->cd(6);
    cnfollowmean_eeffenergy->cd(6)->SetLogx();
    histsNominal->nfollowers_mr_eeffenergy->Draw("colz");

    TCanvas* neutron_efficiency = new TCanvas("neutronEfficiency", "Neutron Efficiency", 1000, 600);
    histsNominal->nfollowers_radial->Sumw2();
    histsNominal->nfollowers_detected_radial->Sumw2();
    histsNominal->nfollowers_detected_radial->Divide(histsNominal->nfollowers_radial);
    histsNominal->nfollowers_detected_radial->Draw("e1");

    _outFile->WriteTObject(cfitpos);
    _outFile->WriteTObject(cseedpos);
    _outFile->WriteTObject(cprompt);
    _outFile->WriteTObject(cmichele);
    _outFile->WriteTObject(cfollowers);
    _outFile->WriteTObject(cnfollow);
    _outFile->WriteTObject(cnhit_nofollow);
    _outFile->WriteTObject(cnfollowmean_eeffenergy);

    // Write histograms to output file for later use - i.e. superimposing data and MC
    histsNominal->writeAllToFileVec(_outFile);
    histsUpper->removeAll();
    histsLower->removeAll();

    _outFile->Close();

}
