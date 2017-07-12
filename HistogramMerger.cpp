//
// Created by Jyotirmai Singh on 7/11/17.
//

#include <iostream>
#include <cmath>
#include "HistogramMerger.h"
#include "HistogramMaker.h"

HistogramMerger::HistogramMerger(TFile* nomFile, TFile* upFile, TFile* lowFile, TFile* outFile) {
    _nomFile = nomFile;
    _upFile = upFile;
    _lowFile = lowFile;
    _outFile = outFile;
}


void HistogramMerger::applyBinErrorCorrection(TH1F* nominal, TH1F* upper, TH1F* lower)
{
    int nBins = nominal->GetSize();
    for (int i = 0; i < nBins; i++) {
        double upperErr = upper->GetBinContent(i);
        double lowerErr = lower->GetBinContent(i);
        double error = fabs(upperErr - lowerErr);
        nominal->SetBinError(i, error);
    }
}

void HistogramMerger::mergeHistograms(TH1F* nominal, TH1F* upper, TH1F* lower,
                                      TCanvas* canv)
{
    TCanvas* canvas = canv;
    TH1F* clonedHisto = nominal->Clone();
    clonedHisto->Draw("e1");
    this->applyBinErrorCorrection(nominal, upper, lower);
    nominal->SetFillColor(kRed);
    nominal->Draw("e2 same");
}

void HistogramMerger::mergeNormDivisionHistograms(TH1F* nominal, TH1F* upper, TH1F* lower,
                                                  TH1F* nominalNorm, TH1F* upperNorm, TH1F* lowerNorm,
                                                  TCanvas* canv)
{
    nominal->Sumw2();
    std::cout << "Done nominal sumw2\n";
    nominalNorm->Sumw2();
    std::cout <<"Done nominal sumw2\n";
    upper->Sumw2();
    std::cout << "Done upper sumw2\n";
    upperNorm->Sumw2();
    std::cout << "Done upper sumw2\n";
    lower->Sumw2();
    std::cout << "Done lower sumw2\n";
    lowerNorm->Sumw2();
    std::cout << "Done lower sumw2\n";

    nominal->Divide(nominalNorm);
    upper->Divide(upperNorm);
    lower->Divide(lowerNorm);
    std::cout << "Done division\n";

    this->mergeHistograms(nominal, upper, lower, canv);
}


void HistogramMerger::makeHistograms()
{
    HistogramMaker* histsNominal = new HistogramMaker(_nomFile);
    HistogramMaker* histsUpper = new HistogramMaker(_upFile);
    HistogramMaker* histsLower = new HistogramMaker(_lowFile);

    std::cout << "Got all histograms ########################" << std::endl;

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
    mergeHistograms(histsNominal->hprompt_eeffenergy,
                          histsUpper->hprompt_eeffenergy,
                          histsLower->hprompt_eeffenergy, cprompt);
    cprompt->cd(5);
    mergeHistograms(histsNominal->hprompt_ueffenergy,
                          histsUpper->hprompt_ueffenergy,
                          histsLower->hprompt_ueffenergy, cprompt);
    cprompt->cd(6);
    mergeHistograms(histsNominal->hprompt_meffenergy,
                          histsUpper->hprompt_meffenergy,
                          histsLower->hprompt_meffenergy, cprompt);
    cprompt->Update();
    std::cout << "Made first set of merged histos ###########################" << std::endl;

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


    std::cout << "Made everything except follower histos ###########################################" << std::endl;

    TCanvas* cnfollowmean_eeffenergy = new TCanvas("cnfollowmean_eeffenergy","Number of Followers vs Prompt Energy",700,900);
    gStyle->SetPalette(56);
    cnfollowmean_eeffenergy->Divide(2,3);
    cnfollowmean_eeffenergy->cd(1);
    cnfollowmean_eeffenergy->cd(1)->SetLogx();
    this->mergeNormDivisionHistograms(histsNominal->nfollowersmean_eeffenergy,
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
    std::cout << "Done cd\n";
    cnfollowmean_eeffenergy->cd(3)->SetLogx();
    std::cout << "Done SetLogx\n";
    this->mergeNormDivisionHistograms(histsNominal->nfollowersmean_sr_eeffenergy,
                                histsUpper->nfollowersmean_sr_eeffenergy,
                                histsLower->nfollowersmean_sr_eeffenergy,
                                histsNominal->nfollowersmean_sr_eeffenergy_norm,
                                histsUpper->nfollowersmean_sr_eeffenergy_norm,
                                histsLower->nfollowersmean_sr_eeffenergy_norm,
                                cnfollowmean_eeffenergy);
    std::cout << "Done MERGE\n";
    cnfollowmean_eeffenergy->cd(4);
    cnfollowmean_eeffenergy->cd(4)->SetLogx();
    histsNominal->nfollowers_sr_eeffenergy->Draw("colz");
    cnfollowmean_eeffenergy->cd(5);
    cnfollowmean_eeffenergy->cd(5)->SetLogx();
    this->mergeNormDivisionHistograms(histsNominal->nfollowersmean_mr_eeffenergy,
                                histsUpper->nfollowersmean_mr_eeffenergy,
                                histsLower->nfollowersmean_mr_eeffenergy,
                                histsNominal->nfollowersmean_mr_eeffenergy_norm,
                                histsUpper->nfollowersmean_mr_eeffenergy_norm,
                                histsLower->nfollowersmean_mr_eeffenergy_norm,
                                cnfollowmean_eeffenergy);
    cnfollowmean_eeffenergy->cd(6);
    cnfollowmean_eeffenergy->cd(6)->SetLogx();
    histsNominal->nfollowers_mr_eeffenergy->Draw("colz");

    _outFile->WriteTObject(cfitpos);
    _outFile->WriteTObject(cseedpos);
    _outFile->WriteTObject(cprompt);
    _outFile->WriteTObject(cmichele);
    _outFile->WriteTObject(cfollowers);
    _outFile->WriteTObject(cnfollow);
    _outFile->WriteTObject(cnhit_nofollow);
    _outFile->WriteTObject(cnfollowmean_eeffenergy);

    _outFile->Close();

}
