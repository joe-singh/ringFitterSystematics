/** @class RingFitterEventCombined.cc
 *  Analyses data in ROOT file and produces lower, nominal, or upper
 *  histograms.
 *  @authors Jyotrimai Singh, Max Smiley, Javier Caravaca */

#define RingFitterEvent_cxx
#include "RingFitterEvent.h"
#include <iostream>
#include <TH1F.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>
#include "HistogramMaker.h"
#include "HighScaleSystematics.h"
#include "LowScaleSystematics.h"

/* Systematic Parameters */

// Low Energy
const double lowEnergyShift = 0.0041;
const double lowEnergyDelta = 0.0;
const double sysUpXYZ = 0.10;
const double sysDownXYZ = -0.57;
const double sysUpScZ = 0.40;
const double sysDownScZ = 0.0;
const double sysUpOffX = 1.15 * 10;
const double sysDownOffX = -0.13 * 10;
const double sysUpOffY = 2.87 * 10;
const double sysDownOffY = -0.17 * 10;
const double sysUpOffZ = 2.58 * 10;
const double sysDownOffZ = -0.15 * 10;
const double sysEdepFidVolUp = 0.85;
const double sysEdepFidVolDown = -0.49;

/* Systematic Flags */
// Edep. Fid. Vol. flag
const bool edepFidVol = true;
// High Energy Flag
const bool highEnergy = true;
// High Energy Shift
const double highEnergyShift = TMath::Abs(37.23 - 38.75) / 37.23;

// TODO Other systematic flags are in LowScaleSystematics.h/.cpp. These two should go there also.

/* Helper Functions */

/** Calculates the norm of the point (x,y,z)
 *  @param x x coordinate.
 *  @param y y coordinate.
 *  @param z z coordinate.
 *  @return x^2 + y^2 + z^2 */
double norm(double x, double y, double z)
{
    return (x * x) + (y * y) + (z * z);
}

/** Applies the Michel Electron Cut.
 *  @param ifol The follower iteration which is being run.
 *  @param rFitEv The place to extract all relevant data, just the
 *                instance of this class running at the moment.
 *  @return true if the conditions of the Michel E cut are met. */
bool micheleCut(int ifol, RingFitterEvent* rFitEv)
{
    return rFitEv->followers_nhit[ifol]>100 &&
           rFitEv->followers_deltat[ifol] > 800e-9 &&
           rFitEv->followers_deltat[ifol] < 20e-6 &&
           (rFitEv->followers_datacleaning2[ifol] & 0x100000) == 0x0;
}

/** Applies the Neutron Cuts
 *  @param ifol The follower iteration which is being run.
 *  @param rFitEv The place to extract all relevant data, just the
 *                instance of this class running at the moment.
 *  @return true if the conditions of the Neutron cut are met.*/
bool neutronCut(int ifol, RingFitterEvent* rFitEv)
{
    double ftkDist = sqrt(norm(rFitEv->followers_ftkpos_fX[ifol], rFitEv->followers_ftkpos_fY[ifol], rFitEv->followers_ftkpos_fZ[ifol]));
    return (rFitEv->followers_deltat[ifol] <= 20e-6  ||
            rFitEv->followers_deltat[ifol] >= 250e-3 ||
            (rFitEv->followers_datacleaning2[ifol] & 0xB56DE1) != 0x0 ||
            ftkDist > 6000.0 ||
            rFitEv->followers_energy[ifol][2] <= 4.0 ||
            rFitEv->followers_ftkitr[ifol] <= 0.5 ||
            !(rFitEv->followers_fit[0]));
}


/** Main Function. Loops over all data and populates all histograms
 *  made in HistogramMaker. Applies low and high scale systematics as
 *  required by flag settings. Saves plots to MC_CombinedSystematicXXX.root,
 *  where XXX is either Nominal, Lower, or Upper depending on dir.
 *  @param USEWATER Whether to use D2O settings.
 *  @param dir The direction of the systematics applied.
 *             1 for upper, -1 for lower, 0 for nominal.
 *  @param data True if using data, false if using MC. */
void RingFitterEvent::Loop(int USEWATER, int dir, bool data)
{
    std::string fileName = "MC_CombinedSystematic";
    std::string dirError = "Incorrect Direction For Loop, must be +1, -1, or 0";
    switch (dir) {
        case 0:
            fileName += "Nominal"; break;
        case 1:
            fileName += "Upper"; break;
        case -1:
            fileName += "Lower"; break;
        default:
            throw new std::invalid_argument(dirError);
    }
    fileName += ".root";

    TFile* outf = new TFile(fileName.c_str(), "recreate");
    HistogramMaker* histograms = new HistogramMaker(0);

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    int nfol_tot = 0;
    Long64_t nbytes = 0;
    Long64_t nb = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        // Simple tracker to monitor progress.
        if (jentry % 100 == 0 && jentry > 0) {
            std::cout << "Entry " << jentry << " complete" << std::endl;
        }

        int seed;
        double spos_fX;
        double spos_fY;
        double spos_fZ;
        double fpos_fX;
        double fpos_fY;
        double fpos_fZ;
        double effen_e;
        double effen_m;


        if(USEWATER==1){
            seed = 0;
            spos_fX = wpos_fX;
            spos_fY = wpos_fY;
            spos_fZ = wpos_fZ;
            fpos_fX = cpos4_fX;
            fpos_fY = cpos4_fY;
            fpos_fZ = cpos4_fZ;
            effen_e = effenergy[0];
            effen_m = effenergy[1];
        }
        else{
            seed = 1;
            spos_fX = qpos_fX;
            spos_fY = qpos_fY;
            spos_fZ = qpos_fZ;
            if(USEWATER==2){
                spos_fX = hpos_fX;
                spos_fY = hpos_fY;
                spos_fZ = hpos_fZ;
            }
            fpos_fX = cpos5_fX;
            fpos_fY = cpos5_fY;
            fpos_fZ = cpos5_fZ;
            effen_e = effenergy[2];
            effen_m = effenergy[3];
        }

        //Seed tests
        histograms->hseedpos[0]->Fill(spos_fX);
        histograms->hseedpos[1]->Fill(spos_fY);
        histograms->hseedpos[2]->Fill(spos_fZ);
        histograms->hseedposdiff[0]->Fill(spos_fX - wpos_fX);
        histograms->hseedposdiff[1]->Fill(spos_fY - wpos_fY);
        histograms->hseedposdiff[2]->Fill(spos_fZ - wpos_fZ);
        histograms->hfitpos[0]->Fill(fpos_fX);
        histograms->hfitpos[1]->Fill(fpos_fY);
        histograms->hfitpos[2]->Fill(fpos_fZ);
        histograms->hfitposdiff[0]->Fill(fpos_fX - cpos4_fX);
        histograms->hfitposdiff[1]->Fill(fpos_fY - cpos4_fY);
        histograms->hfitposdiff[2]->Fill(fpos_fZ - cpos4_fZ);

        // Prompt Cuts
        if(fit[seed] == false) continue;
        if(double(intime[seed]) / double(nhit) <= 0.3) continue;
        if(neck != 0) continue;
        if(owl >= 3) continue;
        if(sqrt(norm(spos_fX, spos_fY, spos_fZ)) >= 8500.) continue;
        if(sqrt(norm(fpos_fX, fpos_fY, fpos_fZ)) <= 0 ||
           sqrt(norm(fpos_fX, fpos_fY, fpos_fZ)) >= 5500) continue;
        if((datacleaning & 0xB62CC103000000) != 0) continue;
        if((trigtype & 0x7) == 0x0) continue;

        // Applying high scale systematics to effen_e
        HighScaleSystematics* HSS = new HighScaleSystematics(highEnergyShift);
        double shiftUpPromptEnergy = HSS->applyEnergyShift(effen_e, 1);
        double shiftDownPromptEnergy = HSS->applyEnergyShift(effen_e, -1);
        double shiftPromptEnergy;

        switch (dir) {
            case 0:
                shiftPromptEnergy = effen_e; break;
            case 1:
                shiftPromptEnergy = shiftUpPromptEnergy; break;
            case -1:
                shiftPromptEnergy = shiftDownPromptEnergy; break;
        }
	
	if (!highEnergy) {
	    shiftPromptEnergy = effen_e;
	}

        histograms->hprompt_nhits->Fill(nhit);
        histograms->hprompt_nrings->Fill(nrings[seed]);
        if (nrings[seed] == 1){
          //Only fill PID where meaningful (i.e. one ring)
          histograms->hprompt_pid->Fill(pid[seed]);

          if (pid[seed] == 0)
             histograms->hprompt_eeffenergy->Fill(shiftPromptEnergy);
          else if (pid[seed] == 1)
             histograms->hprompt_ueffenergy->Fill(effen_m); //why is this not shifted?
        } else if (nrings[seed] == 2){
            histograms->hprompt_meffenergy->Fill(shiftPromptEnergy);
        }

        int nfol = 0;
        double nFolNom;
        double nFolSystematic = 0;
	
        for (int ifol = 0; ifol < followers_; ifol++) {

            // Michel E Cuts
            if (micheleCut(ifol, this)) {
                histograms->hmichele_nhits->Fill(followers_nhit[ifol]);
                histograms->hmichele_deltat->Fill(followers_deltat[ifol]);
                histograms->hmichele_energy0->Fill(followers_energy[ifol][0]);
                histograms->hmichele_energy1->Fill(followers_energy[ifol][1]);
                histograms->hmichele_energy2->Fill(followers_energy[ifol][2]);
            }

            double xPos = followers_ftkpos_fX[ifol];
            double yPos = followers_ftkpos_fY[ifol];
            double zPos = followers_ftkpos_fZ[ifol];
            double followerFTK = followers_energy[ifol][2];

            // Fill radial histogram with neutron radius
            histograms->nfollowers_radial->Fill(sqrt(norm(xPos, yPos, zPos)));

            // Neutron Cuts
            if (neutronCut(ifol, this)) continue;

            histograms->hfollowers_nhits->Fill(followers_nhit[ifol]);
            histograms->hfollowers_deltat->Fill(followers_deltat[ifol]);
            histograms->hfollowers_energy0->Fill(followers_energy[ifol][0]);
            histograms->hfollowers_energy1->Fill(followers_energy[ifol][1]);
            histograms->hfollowers_energy2->Fill(followerFTK);
            histograms->hfollowers_dist->Fill(sqrt(norm((xPos - wpos_fX),
                                                   (yPos - wpos_fY),
                                                   (zPos - wpos_fZ))));


            LowScaleSystematics* LSS = new LowScaleSystematics(lowEnergyShift, lowEnergyDelta, sysUpXYZ, sysDownXYZ,
                                                               sysUpScZ, sysDownScZ, sysUpOffX, sysDownOffX, sysUpOffY,
                                                               sysDownOffY, sysUpOffZ, sysDownOffZ,
                                                               sysEdepFidVolUp, sysEdepFidVolDown);

            double xSystematic;
            double ySystematic;
            double zSystematic;

            TF1* nominalEnergyFunc = new TF1("X", "TMath::Gaus(x, 0.0, 0.155)",
                                             lowEnergyDelta - 10, lowEnergyDelta + 10);
            TF1* energySmearingFunc;

            switch (dir) {
                case 0:
                    xSystematic = xPos;
                    ySystematic = yPos;
                    zSystematic = zPos;
                    energySmearingFunc = nominalEnergyFunc;
                    break;
                case 1:
                    xSystematic = LSS->applyPositionSystematics(xPos, 1, 1);
                    ySystematic = LSS->applyPositionSystematics(yPos, 1, 2);
                    zSystematic = LSS->applyPositionSystematics(zPos, 1, 3);
                    energySmearingFunc = new TF1("X", "TMath::Gaus(x, 0.0, 0.196)",
                                                 lowEnergyDelta - 10, lowEnergyDelta + 10);
                    break;
                case -1:
                    xSystematic = LSS->applyPositionSystematics(xPos, -1, 1);
                    ySystematic = LSS->applyPositionSystematics(yPos, -1, 2);
                    zSystematic = LSS->applyPositionSystematics(zPos, -1, 3);
                    energySmearingFunc = new TF1("X", "TMath::Gaus(x, 0.0, 0.085)",
                                                 lowEnergyDelta - 10, lowEnergyDelta + 10);
                    break;
            }

            double normalSquareDistance = norm(xPos, yPos, zPos);
            double systematicDistance = norm(xSystematic, ySystematic, zSystematic);

            double systematicEnergy = LSS->applySmearing(followerFTK, energySmearingFunc);
            double nominalEnergy = LSS->applySmearing(followerFTK, nominalEnergyFunc);

            // The MC energy correction should not be applied to data.
            if (!data) {

                double r3Systematic = sqrt((systematicDistance))/(600 * 10);
                double r3Nom = sqrt((normalSquareDistance))/(600 * 10);

                double cr3Systematic = 1.01159 - 0.0389943 * r3Systematic + 0.0250065 * pow(r3Systematic, 2);
                double cr3Nom = 1.01159 - 0.0389943 * r3Nom + 0.0250065 * pow(r3Nom, 2);

                systematicEnergy = LSS->correctFollowerEnergy(systematicEnergy, cr3Systematic);
                nominalEnergy = LSS->correctFollowerEnergy(nominalEnergy, cr3Nom);
            }

            double wPlus = LSS->getwPlus(nominalEnergy);
            double wMinus = LSS->getwMinus(nominalEnergy);

           if (sqrt(systematicDistance) <= 6000.0 && systematicEnergy > 4.0) {
               histograms->nfollowers_detected_radial->Fill(sqrt(normalSquareDistance));
               switch (dir) {
                   case 0:
                       nFolSystematic++;
                       break;
                   case 1:
                   {
                       if (edepFidVol) {
                           nFolSystematic += wPlus;
                       } else {
                           nFolSystematic++;
                       }
                       break;
                   }
                   case -1:
                   {
                       if (edepFidVol) {
                           nFolSystematic += wMinus;
                       } else {
                           nFolSystematic++;
                       }
                       break;
                   }
               }
           }

            nfol++;
            nfol_tot++;
        }

        histograms->nfollowers_tot->Fill(nfol);
        histograms->nfollowers_eeffenergy->Fill(shiftPromptEnergy, nfol);
        histograms->nfollowersmean_eeffenergy->Fill(shiftPromptEnergy, nFolSystematic);
        histograms->nfollowersmean_eeffenergy_norm->Fill(shiftPromptEnergy);

        if (nrings[seed] == 1){
            histograms->nfollowers_sring->Fill(nfol);
            histograms->nfollowers_sr_eeffenergy->Fill(shiftPromptEnergy, nfol);
            histograms->nfollowersmean_sr_eeffenergy->Fill(shiftPromptEnergy, nFolSystematic);
            histograms->nfollowersmean_sr_eeffenergy_norm->Fill(shiftPromptEnergy);
        } else if (nrings[seed] == 2){
            histograms->nfollowers_mring->Fill(nfol);
            histograms->nfollowers_mr_eeffenergy->Fill(shiftPromptEnergy,nfol);
            histograms->nfollowersmean_mr_eeffenergy->Fill(shiftPromptEnergy,nFolSystematic);
            histograms->nfollowersmean_mr_eeffenergy_norm->Fill(shiftPromptEnergy);
        }
        if(nfol == 0) {
            histograms->nhit_nofollow_tot->Fill(nhit);
            if(nrings[seed] == 1){
                histograms->nhit_nofollow_sring->Fill(nhit);
            } else if(nrings[seed] == 2){
                histograms->nhit_nofollow_mring->Fill(nhit);
            }
        }

    }

    histograms->writeAllToFile(outf);
    outf->Close();
}
