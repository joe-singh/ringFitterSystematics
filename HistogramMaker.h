//
// Created by Jyotirmai Singh on 7/7/17.
//

#ifndef HSTOGRAMMAKER_H
#define HISTOGRAMMAKER_H


class HistogramMaker {

public:
    TH1F* hseedpos[3];
    TH1F* hseedposdiff[3];
    TH1F* hfitpos[3];
    TH1F* hfitposdiff[3];

    // Prompts
    TH1F* hprompt_nhits;
    TH1F* hprompt_nrings;
    TH1F* hprompt_pid;
    TH1F* hprompt_eeffenergy;
    TH1F* hprompt_ueffenergy;
    TH1F* hprompt_meffenergy;

    // Michel Electrons
    TH1F* hmichele_nhits;
    TH1F* hmichele_deltat;
    TH1F* hmichele_energy0;
    TH1F* hmichele_energy1;
    TH1F* hmichele_energy2;

    // Neutrons
    TH1F* hfollowers_nhits;
    TH1F* hfollowers_deltat;
    TH1F* hfollowers_energy0;
    TH1F* hfollowers_energy1;
    TH1F* hfollowers_energy2;
    TH1F* hfollowers_dist;
    TH1F* nfollowers_tot;
    TH1F* nfollowers_sring;
    TH1F* nfollowers_mring;
    TH1F* nhit_nofollow_tot;
    TH1F* nhit_nofollow_sring;
    TH1F* nhit_nofollow_mring;
    TH2F* nfollowers_eeffenergy;
    TH1F* nfollowersmean_eeffenergy;
    TH1F* nfollowersmean_eeffenergy_norm;
    TH2F* nfollowers_sr_eeffenergy;
    TH1F* nfollowersmean_sr_eeffenergy;
    TH1F* nfollowersmean_sr_eeffenergy_norm;
    TH2F* nfollowers_mr_eeffenergy;
    TH1F* nfollowersmean_mr_eeffenergy;
    TH1F* nfollowersmean_mr_eeffenergy_norm;
    TH1F* nfollowers_radial;
    TH1F* nfollowers_detected_radial;

    std::vector<TH1*> hists;

    HistogramInitializer();
    HistogramMaker(TFile* file);
    void writeAllToFile(TFile* file);
    void writeAllToFileVec(TFile* file);
    void removeAll();
private:
    void binLogX(TH1* h, bool logLimits);
    void binLogX(TH2* h, bool logLimits);
    TH1F* getHisto(TH1F* hist, TFile* file);
    TH2F* get2DHisto(TH2F* hist, TFile* file);
};


#endif
