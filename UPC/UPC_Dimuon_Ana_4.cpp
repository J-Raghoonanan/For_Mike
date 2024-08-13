#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include <TFile.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphPainter.h>

#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include "TF1.h"
#include <bitset>


#include "TTree.h"
#include "TFile.h"
#include "math.h"
#include "TMath.h"

#include <cmath>

#include <Math/Vector4D.h>
#include "VAProgress.h"

using namespace std;

typedef struct { Float_t a;} BRANCH;


TChain *chain[1][5];  //To Load a File
int colors[12] = {4, 6, 8, 28, 1, 2, 38, 42, 45, 46, 5, 30};
int colors2[4][2] = {{4, 6}, {8, 28}, {1, 2}, {38, 42}};
TLegend     *leg;
TCanvas     *c0;
TCanvas     *c1;
TCanvas     *c2;
TCanvas     *c3;
TCanvas     *c4;
TCanvas     *c5;
TCanvas     *c6;
TCanvas     *c7;
TCanvas     *c8;
TCanvas     *c9;
TCanvas     *c10;

char name[100];
int ModuleMask_req[2] = {15, 240};
TString leg_key[4] = {"w. 0n0n", "w. XnXn", "w. 0nXn/Xn0n"};//Which zdc side is triggered
TString leg_key2[4][4] = {
    {"0n0n", "0n1n", "0n2n", "0n3/4n"},
    {"1n0n", "1n1n", "1n2n", "1n3/4n"},
    {"2n0n", "2n1n", "2n2n", "2n3/4n"},
    {"3/4n0n", "3/4n1n", "3/4n2n", "3/4n3/4n"}
};

TString comp_key[2] = {"New", "Old"};
TString TG[2] = {"EMG Sigma Parameter", "EMG Tau Parameter"};
TString zdc[2] = {"ZdcEtA", "ZdcEtC"};


float sigma_parameter[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};
float tau_parameter[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};
float exp_slope_parameter[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};

float dsigma[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};
float dtau[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};
float sigma_y[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};
float sigma_exp[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};
float mean[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};

float exp_ratio[4][4] = {
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1},
    { -1, -1, -1, -1}
};




int unique_counter[16] = {0, 1, 2, 3, -1, 4, 5, 6, -1, -1, 7, 8, -1, -1, -1, 9};

float log_binning[28] = {1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 1e-1};

double xBins[28];


auto pi = TMath::Pi();






#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);

    if (percentage == 1.0)
    {
        cout << endl;
    }
}








void create_chain()
{
    for (int i = 0; i < 2; i++)
    {
        // chain[0][i] = new TChain("zdcNtuple"); //Might have to change this
        chain[0][i] = new TChain("HeavyIonD3PD");

        if (i == 0)
        {
            // sprintf(name, "/Users/jonathanraghoonanan/Desktop/Data_Files/UPC_Dimuons/UPC_dimuon_ntuple.root"); //UPC dimuon small sample
            sprintf(name, "/Users/jonathanraghoonanan/Desktop/CU_docs/Research/Data_Files/UPC_Dimuons/New_All_Data2015_24Dec2020.root"); //
        }
        else if (i == 1)
        {
            sprintf(name, "/Users/jonathanraghoonanan/Desktop/CU_docs/Research/Data_Files/UPC_Dimuons/New_All_Data2018_24Dec2020.root"); //
        }
        else if (i == 2)
        {
            sprintf(name, "/Users/jonathanraghoonanan/Desktop/CU_docs/Research/Data_Files/UPC_Dimuons/New_All_Data2015.root"); //
        }
        else if (i == 3)
        {
            sprintf(name, "/Users/jonathanraghoonanan/Desktop/CU_docs/Research/Data_Files/UPC_Dimuons/New_All_Data2018.root"); //
        }
        else if (i == 4)
        {
            sprintf(name, "/Users/jonathanraghoonanan/Desktop/CU_docs/Research/Data_Files/"); //
        }





        ifstream inFile;
        inFile.open(name);

        if (!inFile)
        {
            cout << "--File " << name << " does not exist!" << endl;
            continue;

        }
        else
        {
            chain[0][i]->Add(TString(name));
            cout << "Adding " << name << " to chain." << endl;
        }
    }
}


unsigned int CountBits(int data)
{
    std::bitset<4> bs(data);
    return bs.count();
}






void ZdcEt_Plots(TString title, int ind) //compare old and new energy plot
{


    char hist[100];
    char con[100];
    c1->cd(ind + 1);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit();
    leg = new TLegend(.55, .55, .70, .70);
    leg->SetTextSize(.05);
    leg->SetBorderSize(0);
    leg->SetTextFont(12);
    leg->SetFillColor(0);
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);
    vector<float> *muon_pair_acop = 0;
    vector<int> *vtx_ntrk = 0;
    vector<float> *muon_pair_muon1_pt = 0;
    vector<float> *muon_pair_muon2_pt = 0;
    float ZdcEtA = 0.0;//contains ZDC energies
    float ZdcEtC = 0.0;
    TString xtitle = "";

    for (int comp = 0; comp < 2; comp++)
    {
        sprintf(hist, "hist%d_%d", ind, comp);
        //interested in <1 GeV?
        TH1* histo = new TH1F(hist, hist, 40, 0, 4);
        // TH1* histo = new TH1F(hist, hist, 39, 0.1, 4);
        // TH1* histo = new TH1F(hist, hist, 20, 0.15, .3);

        for (int compare = 0; compare < 2; compare++)  //
        {
            chain[0][compare]->SetBranchAddress("vtx_ntrk", &vtx_ntrk);
            chain[0][compare]->SetBranchAddress("ZdcEtA", &ZdcEtA);
            chain[0][compare]->SetBranchAddress("ZdcEtC", &ZdcEtC);
            chain[0][compare]->SetBranchAddress("muon_pair_muon1_pt", &muon_pair_muon1_pt);
            chain[0][compare]->SetBranchAddress("muon_pair_muon2_pt", &muon_pair_muon2_pt);
            chain[0][compare]->SetBranchAddress("muon_pair_acop", &muon_pair_acop);

            cout << "Number of events: " << chain[0][compare]->GetEntries() << endl; //why is this 0

            for (int i = 0; i < chain[0][compare]->GetEntries(); i++) //events
            {
                chain[0][compare]->GetEntry(i);
                if ((*vtx_ntrk)[0] == 2) //selects UPC dimuon events
                {
                    int  NPairs = muon_pair_acop->size();
                    for (int j = 0; j < NPairs; j++)
                    {
                        float p1 = 0.0;
                        float p2 = 0.0;

                        p1 = muon_pair_muon1_pt->at(j) / 1000; //pt in MeV, scale to GeV
                        p2 = muon_pair_muon2_pt->at(j) / 1000;

                        if (muon_pair_muon1_pt->at(j) * muon_pair_muon2_pt->at(j) <= 0) //select opposite sign muons
                        {
                            if (comp == 0)
                            {
                                histo->Fill(ZdcEtA / 2.51); //energy in TeV
                                // xtitle = "ZdcEtA / 2.51 TeV";
                            }
                            else if (comp == 1)
                            {
                                histo->Fill(ZdcEtC / 2.51);
                                // xtitle = "ZdcEtC / 2.51 TeV";
                            }
                        }
                        //
                    }



                    // cout << "HEREEEEEEE" << endl;

                }
            } //end of event loop
        } //end of compare

        //moved this from inside the compare loop
        // histo->Scale(1. / histo->Integral());
        histo->Draw( (comp == 0) ? "" : "same" );
        histo->SetLineColor(colors[comp]);
        histo->SetTitle(title);
        histo->GetXaxis()->SetTitle("ZdcEt / 2.51 TeV");
        histo->GetYaxis()->SetTitle("Number of Events");
        // histo->GetYaxis()->SetRangeUser(0.01, 70000.);
        gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn
        leg->AddEntry(histo, zdc[comp]);
        leg->Draw();
        if (ind == 1)
        {
            gPad->SetLogy(1);
        }
        gPad->SetTicks(1, 1);
    } // end of comp
}




void k_perp_Tgraph(TString title, int ind = 0, int c_ind = 0) //ind 0 for Kp, 1 for acop
{
    if (c_ind == 0 || c_ind == 5)
    {
        c6->cd();
    }
    else if (c_ind == 1)
    {
        c7->cd();
    }
    else if (c_ind == 2)
    {
        c8->cd();
    }
    else if (c_ind == 3)
    {
        c9->cd();
    }
    else if (c_ind == 4)
    {
        c10->cd();
    }
    // c6->cd();
    gStyle->SetOptStat(1);
    gStyle->SetOptFit();
    leg = new TLegend(.40, .40, .55, .60);
    leg->SetTextSize(.03);
    leg->SetBorderSize(0);
    leg->SetTextFont(12);
    leg->SetFillColor(0);
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.09);
    gPad->SetLeftMargin(0.15);

    TString ylab;
    // int n = 9;
    int n = 10;

    double x[n];
    double y[n];
    double ey[n];
    double ex[n];

    double y3[n];
    double x3[n];
    double ey3[n];
    double ex3[n];
    auto mg = new TMultiGraph("mg", title);

    cout << "Starting 1D TGraph with errors " << endl;

    for (int Aen = 0; Aen < 4; Aen++)
    {
        for (int Cen = 0; Cen < 4; Cen++)
        {
            if ((4 * Aen + Cen == 0) || (4 * Aen + Cen == 1) || (4 * Aen + Cen == 2) || (4 * Aen + Cen == 3) || (4 * Aen + Cen == 5) || (4 * Aen + Cen == 6) || (4 * Aen + Cen == 7) ||  (4 * Aen + Cen == 10) || (4 * Aen + Cen == 11) || (4 * Aen + Cen == 15))
            {
                float s_sq = pow(sigma_parameter[Aen][Cen], 2);
                float t_sq = pow(tau_parameter[Aen][Cen], 2);
                y[unique_counter[4 * Aen + Cen]] = sqrt((s_sq) + 2 * (t_sq)); //units of GeV ; mult by 1000 to get MeV
                y3[unique_counter[4 * Aen + Cen]] = mean[Aen][Cen];

                if (ind == 0)
                {
                    ey[unique_counter[4 * Aen + Cen]] = 1000 * (1 / (y[unique_counter[4 * Aen + Cen]])) * sigma_y[Aen][Cen];
                    y[unique_counter[4 * Aen + Cen]] = 1000 * y[unique_counter[4 * Aen + Cen]];
                    // units of MeV for K_perp
                }
                else if (ind == 1)
                {
                    ey[unique_counter[4 * Aen + Cen]] = (1 / (y[unique_counter[4 * Aen + Cen]])) * sigma_y[Aen][Cen]; //includes covariance
                    // don't need scaling factor for acoplanarity
                }



                // ey[4 * Aen + Cen] = 1000 * ey[4 * Aen * Cen]; // * 1000 [MeV];


                x[unique_counter[4 * Aen + Cen]] = unique_counter[4 * Aen + Cen];
                ex[unique_counter[4 * Aen + Cen]] = 0.0;

                x3[unique_counter[4 * Aen + Cen]] = unique_counter[4 * Aen + Cen];
                ex3[unique_counter[4 * Aen + Cen]] = 0.0;
                ey3[unique_counter[4 * Aen + Cen]] = 0.0;

                cout << "Regular: (" << x[unique_counter[4 * Aen + Cen]] << ", " << y[unique_counter[4 * Aen + Cen]] << ") with error: " << ey[unique_counter[4 * Aen + Cen]] << endl;

                cout << "EMG Mean: (" << x[unique_counter[4 * Aen + Cen]] << ", " << y3[unique_counter[4 * Aen + Cen]] << ")" << endl;
            }
        }

    }



    // double x2[6] = {0,1,0nXn, 4, 1nXn, XnXn };
    double x2[6] = {0, 1, 2, 4, 5, 8};
    double y2[6] = {1.16e-3, 1.28e-3, 1.30e-3, 1.37e-3, 1.42e-3, 1.48e-3}; //Guesstimated CMS Values
    double ey2[6] = {0, 0, 0, 0, 0, 0};
    double ex2[6] = {0, 0, 0, 0, 0, 0};

    auto gr1 = new TGraphErrors(n, x, y, ex, ey);
    gr1->SetName("gr1");
    gr1->SetTitle("EMG RMS");
    gr1->SetMarkerStyle(21);
    gr1->SetDrawOption("ALP");
    gr1->SetLineColor(2);
    gr1->SetLineWidth(4);
    gr1->SetFillStyle(0);

    auto gr2 = new TGraphErrors(6, x2, y2, ex2, ey2);
    gr2->SetName("gr2");
    gr2->SetTitle("CMS <#alpha>");
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(2);
    gr2->SetDrawOption("ALP");
    gr2->SetLineColor(3);
    gr2->SetLineWidth(4);
    gr2->SetFillStyle(0);

    auto gr3 = new TGraphErrors(n, x3, y3, ex3, ey3);
    gr2->SetName("gr3");
    gr2->SetTitle("EMG <#alpha>");
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(2);
    gr2->SetDrawOption("ALP");
    gr2->SetLineColor(3);
    gr2->SetLineWidth(4);
    gr2->SetFillStyle(0);

    mg->Add( gr1 );
    if (ind == 1)
    {
        mg->Add( gr2 );
        mg->Add( gr3 );
    }

    TAxis* xax = mg->GetXaxis();
    // ChangeLabel(labNum, labAngle, labSize, labAlign, labColor, labFont, labText)
    // xax->ChangeLabel(1, 30., .5, -1, -1, -1, "0n0n");
    // xax->ChangeLabel(2, 30., .5, -1, -1, -1, "0n1n");
    // xax->ChangeLabel(3, 30., .5, -1, -1, -1, "0n2n");
    // xax->ChangeLabel(4, 30., .5, -1, -1, -1, "0n3/4n");
    // xax->ChangeLabel(5, 30., .5, -1, -1, -1, "1n1n");
    // xax->ChangeLabel(6, 30., .5, -1, -1, -1, "1n2n");
    // xax->ChangeLabel(7, 30., .5, -1, -1, -1, "1n3/4n");
    // xax->ChangeLabel(8, 30., .5, -1, -1, -1, "2n2n");
    // xax->ChangeLabel(9, 30., .5, -1, -1, -1, "2n3/4n");
    // xax->ChangeLabel(10, 30., .5, -1, -1, -1, "3/4n3/4n");

    // TText t;
    // t.SetTextSize(0.02);
    // t.SetTextAlign(22);
    // Double_t xt = 1;
    // for (Int_t i = 1; i <= 6; i++)
    // {
    //     t.DrawText(mg->GetBinCenter(i), xt, "0n0n");
    // }

    gPad->Modified(); gPad->Update();
    leg->AddEntry(gr1, "EMG RMS with Errors", "lep");
    leg->AddEntry(gr2, "CMS <#alpha>", "lep");
    if (ind == 1)
    {
        leg->AddEntry(gr3, "EMG <#alpha>", "lep");
    }

    if (c_ind == 0)
    {
        leg->AddEntry(gr1, "3.75 < Muon P_t < 5");
    }
    else if (c_ind == 1)
    {
        leg->AddEntry(gr1, "5 < Muon P_t < 10");
    }
    else if (c_ind == 2)
    {
        leg->AddEntry(gr1, "10 < Muon P_t < 15");
    }
    else if (c_ind == 3)
    {
        leg->AddEntry(gr1, "15 < Muon P_t < 20");
    }
    else if (c_ind == 4)
    {
        leg->AddEntry(gr1, "20 < Muon P_t < 25");
    }


    mg->SetTitle(title);
    mg->GetYaxis()->SetRange(0, 50);
    gPad->Modified(); gPad->Update();
    mg->GetXaxis()->SetTitle("Neutron Index");
    if (ind == 0)
    {
        ylab = "RMS = #sqrt{2#tau^{2} + #sigma^{2}} [MeV]";
    }
    else if (ind == 1)
    {
        ylab = "RMS = #sqrt{2#tau^{2} + #sigma^{2}}";
    }

    mg->GetYaxis()->SetTitle(ylab);
    gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn
    gPad->SetTicks(1, 1);
    mg->Draw("ALP");
    leg->Draw();
    gPad->Modified(); gPad->Update();

    cout << "Finished with 1D TGraph " << endl << endl;

}










void k_perp_plot_comp3_dup(TString title, int c_ind) //compare old and new energy plot // Primary fitting function
{
    double P1 = .0008164;
    double P2 = .02022;
    double P3 = 105.019;
    double P4 = 0.29888;

    // Plotting tail region, getting exp fit for initial values

    for (int comp = 2; comp < 3; comp++) // 0-2 normally; 0-8 for specific energies
    {
        c0->cd();
        leg = new TLegend(.40, .40, .55, .60);
        leg->SetTextSize(.03);
        leg->SetBorderSize(0);
        leg->SetTextFont(12);
        leg->SetFillColor(0);
        char hist[100];
        char con[100];
        vector<float> *muon_pair_acop = 0;
        bool b_L1_ZDC_A_C_L1TBP;
        bool b_L1_ZDC_A_L1TBP;
        bool b_L1_ZDC_C_L1TBP;
        vector<int> *vtx_ntrk = 0;
        vector<float> *muon_pair_muon1_pt = 0;
        vector<float> *muon_pair_muon2_pt = 0;
        float ZdcEtA = 0.0;//contains ZDC energies
        float ZdcEtC = 0.0;


        sprintf(hist, "hist%d", comp);
        TH1* histo0 = new TH1F(hist, hist, 25, 0.15, .4);
        TF1* f2 = new TF1("f2", " [0]*exp(-x/[1])", .15, .4);
        f2->SetParLimits(0, .1, 100000);
        f2->SetParLimits(1, 1e-6, 1);
        cout << "Computing exponential fit to tail region [0.15,0.4]" << endl;
        for (int compare = 0; compare < 2; compare++)  //
        {
            chain[0][compare]->SetBranchAddress("muon_pair_acop", &muon_pair_acop);
            chain[0][compare]->SetBranchAddress("muon_pair_muon1_pt", &muon_pair_muon1_pt);
            chain[0][compare]->SetBranchAddress("muon_pair_muon2_pt", &muon_pair_muon2_pt);
            chain[0][compare]->SetBranchAddress("vtx_ntrk", &vtx_ntrk);
            chain[0][compare]->SetBranchAddress("b_L1_ZDC_A_C_L1TBP", &b_L1_ZDC_A_C_L1TBP);
            chain[0][compare]->SetBranchAddress("b_L1_ZDC_A_L1TBP", &b_L1_ZDC_A_L1TBP);
            chain[0][compare]->SetBranchAddress("b_L1_ZDC_C_L1TBP", &b_L1_ZDC_C_L1TBP);
            chain[0][compare]->SetBranchAddress("b_L1_ZDC_C_L1TBP", &b_L1_ZDC_C_L1TBP);
            chain[0][compare]->SetBranchAddress("ZdcEtA", &ZdcEtA);
            chain[0][compare]->SetBranchAddress("ZdcEtC", &ZdcEtC);

            cout << "Number of events: " << chain[0][compare]->GetEntries() << endl; //why is this 0
            // for (int i = 0; i < 100000; i++)
            for (int i = 0; i < chain[0][compare]->GetEntries(); i++) //events
            {
                chain[0][compare]->GetEntry(i);

                bool a = ((*vtx_ntrk)[0] == 2) && (b_L1_ZDC_A_C_L1TBP == 0);
                bool b = ((*vtx_ntrk)[0] == 2) && (b_L1_ZDC_A_C_L1TBP == 1);
                bool c = ((*vtx_ntrk)[0] == 2) && (b_L1_ZDC_A_L1TBP == 0) && (b_L1_ZDC_C_L1TBP == 1);
                bool d = ((*vtx_ntrk)[0] == 2) && (b_L1_ZDC_A_L1TBP == 1) && (b_L1_ZDC_C_L1TBP == 0);
                bool e = c || d;

                bool conds[3] = {a, b, e};
                if (conds[comp % 3])
                {
                    // cout << "HEREEEEEEE" << endl;
                    int  NPairs = muon_pair_acop->size();
                    for (int j = 0; j < NPairs; j++)
                    {
                        float aco = -1.0;
                        float p1 = 0.0;
                        float p2 = 0.0;
                        float avg_p = 0.0;
                        float kperp = 0.0;

                        p1 = abs(muon_pair_muon1_pt->at(j) / 1000); //pt in MeV, scale to GeV
                        p2 = abs(muon_pair_muon2_pt->at(j) / 1000);
                        aco = abs(muon_pair_acop->at(j));

                        avg_p = (p1 + p2) / 2;
                        kperp = avg_p * aco * pi;

                        // cout << p << endl;
                        if (muon_pair_muon1_pt->at(j) * muon_pair_muon2_pt->at(j) <= 0)
                        {
                            if ( (p1 > 3.75) && (p2 > 3.75) )
                            {
                                histo0->Fill(kperp);
                            }

                        }
                        //
                    }
                }

            } //end of event loop

        } //end of compare

        //moved this from inside the compare loop
        float po = histo0->GetBinContent(1);
        f2->SetParameters(po, .02);
        histo0->Fit("f2");

        P3 = f2->GetParameter(0);
        P4 = f2->GetParameter(1);
        cout << "P4 is given as: " << P4 << endl;

        histo0->Draw();
        gPad->Modified(); gPad->Update();
        // histo->SetLineColor(colors[comp]);
        histo0->SetTitle("Exponential Tail Region");
        histo0->GetXaxis()->SetTitle("K_perp");
        histo0->GetYaxis()->SetTitle("Number of Events");
        // histo0->GetYaxis()->SetRangeUser(0.01, 70000.);
        leg->AddEntry(histo0, "UPC Dimuon K_perp " + leg_key[comp % 3]);
        gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn
        leg->Draw();
        // // gPad->SetLogy(1);
        gPad->SetTicks(1, 1);



        cout << "Finished computing exponential tail values " << endl;
        cout << "(P3,P4) = (" << P3 << ", " << P4 << ")" << endl << endl;
    } //end of comp







    // Plotting full 4x4 plots
    cout << "Starting 4x4 matrix of plots " << endl;

    TObjArray Hlist(0);


    int nBins = 50; // 51 edges
    // int nBins = 10; // 31 edges

    for (int i = 0; i <= nBins; i++)
    {
        // This calculation gets your log bins, with maximum value of A and minimum bin of A*10^-(nBins*cx). You can modify cx to change the log scaling of the binning and modify A to make linear changes.
        float A = 1.0; // Coefficient of the log-x bins
        // float cx = 3.0 / 27.; // Coefficient of the exponential scaling
        float cx = 4.0 / 50.0;
        xBins[i] = A * pow(10.0, ((float)(i - nBins)) * cx );
        cout << xBins[i] << endl;
    }

    for (int Aen = 0; Aen < 4; Aen++) // A side energy
    {
        for (int Cen = 0; Cen < 4; Cen++) // C side energy (4x4 matrix in total)
        {

            if ((4 * Aen + Cen == 0) || (4 * Aen + Cen == 1) || (4 * Aen + Cen == 2) || (4 * Aen + Cen == 3) || (4 * Aen + Cen == 5) || (4 * Aen + Cen == 6) || (4 * Aen + Cen == 7) ||  (4 * Aen + Cen == 10) || (4 * Aen + Cen == 11) || (4 * Aen + Cen == 15))
            {
                char hist[100];
                char con[100];
                if (c_ind == 0 || c_ind == 5)
                {
                    c1->cd(unique_counter[4 * Aen + Cen] + 1);
                }
                else if (c_ind == 1)
                {
                    c2->cd(unique_counter[4 * Aen + Cen] + 1);
                }
                else if (c_ind == 2)
                {
                    c3->cd(unique_counter[4 * Aen + Cen] + 1);
                }
                else if (c_ind == 3)
                {
                    c4->cd(unique_counter[4 * Aen + Cen] + 1);
                }
                else if (c_ind == 4)
                {
                    c5->cd(unique_counter[4 * Aen + Cen] + 1);
                }
                gStyle->SetOptStat(1);
                gStyle->SetOptFit();
                leg = new TLegend(.25, .40, .40, .60);
                leg->SetTextSize(.03);
                leg->SetBorderSize(0);
                leg->SetTextFont(12);
                leg->SetFillColor(0);
                gStyle->SetOptStat(0);
                gPad->SetRightMargin(0.09);
                gPad->SetLeftMargin(0.15);
                vector<float> *muon_pair_acop = 0;
                bool b_L1_ZDC_A_C_L1TBP;
                bool b_L1_ZDC_A_L1TBP;
                bool b_L1_ZDC_C_L1TBP;
                vector<int> *vtx_ntrk = 0;
                vector<float> *muon_pair_muon1_pt = 0;
                vector<float> *muon_pair_muon2_pt = 0;
                float ZdcEtA = 0.0;//contains ZDC energies
                float ZdcEtC = 0.0;
                vector<float> *muon_pair_muon1_eta = 0;
                vector<float> *muon_pair_muon2_eta = 0;
                vector<float> *muon_pair_muon1_phi = 0;
                vector<float> *muon_pair_muon2_phi = 0;

                sprintf(hist, "hist%d_%d_%d", Aen, Cen, c_ind);
                // TH1* histo = new TH1F(hist, hist, 50, 0, 1);
                TH1* histo = new TH1F(hist, hist, nBins, xBins);
                cout << "Drawing " << hist << endl;
                // TH1* histo = new TH1F(hist, hist, 20, 0.15, .3);

                Hlist.Add(histo); //begin writing the histogram to root file

                TF1* f4 = new TF1("f4", "[0]*(exp(([1]**2)/(2*([2]**2))) * ((.5 * (exp(x/[2]) * erfc((1/sqrt(2)) * ((x/[1]) + ([1]/[2]))))) +  (exp(-x/[2]) * (1- .5* erfc((1/sqrt(2)) * ((x/[1]) - ([1]/[2])))))) ) + [3]*exp(-x/[4])", 0, 1); //exp modified gaussian

                f4->SetParLimits(0, 0, 1000000);
                f4->SetParLimits(1, 1e-6, .1);
                f4->SetParLimits(2, 1e-6, .6);
                f4->SetParLimits(3, 1, 1000);
                f4->SetParLimits(4, .01, 1);

                TF1* f5;

                // if ((Aen == 2 && Cen == 3) || (Aen == 3 && Cen == 3))
                if ((Aen == 3 && Cen == 3))
                {
                    f5 = new TF1("f5", "[0]*(1/(2.0*[2]))*(exp(([1]**2.0)/(2.0*([2]**2.0))) * ((.5 * (exp(x/[2]) * erfc((1.0/sqrt(2.0)) * ((x/[1]) + ([1]/[2]))))) +  (exp(-x/[2]) * (1.0- .5* erfc((1.0/sqrt(2.0)) * ((x/[1]) - ([1]/[2])))))) )", 0, 1); //exp modified gaussian
                }
                else
                {
                    // Using fixed exp slopes
                    f5 = new TF1("f5", "[0]*(exp(([1]**2.0)/(2.0*([2]**2.0))) * ((.5 * (exp(x/[2]) * erfc((1.0/sqrt(2.0)) * ((x/[1]) + ([1]/[2]))))) +  (exp(-x/[2]) * (1.0- .5* erfc((1.0/sqrt(2.0)) * ((x/[1]) - ([1]/[2])))))) ) + [3]*exp(-x/[4])", 0, 1); //exp modified gaussian

                    f5->SetParLimits(3, 1, 1000);
                }

                f5->SetParLimits(0, 0, 1000000);
                f5->SetParLimits(1, 1e-6, .1);
                f5->SetParLimits(2, 1e-4, .6);


                for (int compare = 0; compare < 2; compare++)  //
                {
                    chain[0][compare]->SetBranchAddress("muon_pair_acop", &muon_pair_acop);
                    chain[0][compare]->SetBranchAddress("muon_pair_muon1_pt", &muon_pair_muon1_pt);
                    chain[0][compare]->SetBranchAddress("muon_pair_muon2_pt", &muon_pair_muon2_pt);
                    chain[0][compare]->SetBranchAddress("vtx_ntrk", &vtx_ntrk);
                    chain[0][compare]->SetBranchAddress("b_L1_ZDC_A_C_L1TBP", &b_L1_ZDC_A_C_L1TBP);
                    chain[0][compare]->SetBranchAddress("b_L1_ZDC_A_L1TBP", &b_L1_ZDC_A_L1TBP);
                    chain[0][compare]->SetBranchAddress("b_L1_ZDC_C_L1TBP", &b_L1_ZDC_C_L1TBP);
                    chain[0][compare]->SetBranchAddress("b_L1_ZDC_C_L1TBP", &b_L1_ZDC_C_L1TBP);
                    chain[0][compare]->SetBranchAddress("ZdcEtA", &ZdcEtA);
                    chain[0][compare]->SetBranchAddress("ZdcEtC", &ZdcEtC);
                    chain[0][compare]->SetBranchAddress("muon_pair_muon1_eta", &muon_pair_muon1_eta);
                    chain[0][compare]->SetBranchAddress("muon_pair_muon2_eta", &muon_pair_muon2_eta);
                    chain[0][compare]->SetBranchAddress("muon_pair_muon1_phi", &muon_pair_muon1_phi);
                    chain[0][compare]->SetBranchAddress("muon_pair_muon2_phi", &muon_pair_muon2_phi);

                    cout << "Number of events: " << chain[0][compare]->GetEntries() << endl; //why is this 0

                    // for (int i = 0; i < 100000; i++)
                    for (int i = 0; i < chain[0][compare]->GetEntries(); i++) //events
                    {
                        chain[0][compare]->GetEntry(i);
                        // z zero, o one, t two, tf three four
                        bool ZnZn = (ZdcEtA <= 1) && (ZdcEtC <= 1);
                        bool ZnOn = (ZdcEtA <= 1) && (ZdcEtC > 1) && (ZdcEtC <= 1.5 * 2.51);
                        bool ZnTn = (ZdcEtA <= 1) && (ZdcEtC > 1.5 * 2.51) && (ZdcEtC <= 2.5 * 2.51);
                        bool ZnTFn = (ZdcEtA <= 1) && (ZdcEtC > 2.5 * 2.51) && (ZdcEtC <= 4.5 * 2.51);
                        bool OnZn = (ZdcEtA > 1) && (ZdcEtA <= 1.5 * 2.51) && (ZdcEtC <= 1);
                        bool OnOn = (ZdcEtA > 1) && (ZdcEtA <= 1.5 * 2.51) && (ZdcEtC > 1) && (ZdcEtC <= 1.5 * 2.51);
                        bool OnTn = (ZdcEtA > 1) && (ZdcEtA <= 1.5 * 2.51) && (ZdcEtC > 1.5 * 2.51) && (ZdcEtC <= 2.5 * 2.51);
                        bool OnTFn = (ZdcEtA > 1) && (ZdcEtA <= 1.5 * 2.51) && (ZdcEtC > 2.5 * 2.51) && (ZdcEtC <= 4.5 * 2.51);
                        bool TnZn = (ZdcEtA > 1.5 * 2.51) && (ZdcEtA <= 2.5 * 2.51) && (ZdcEtC <= 1);
                        bool TnOn = (ZdcEtA > 1.5 * 2.51) && (ZdcEtA <= 2.5 * 2.51) && (ZdcEtC > 1) && (ZdcEtC <= 1.5 * 2.51);
                        bool TnTn = (ZdcEtA > 1.5 * 2.51) && (ZdcEtA <= 2.5 * 2.51) && (ZdcEtC > 1.5 * 2.51) && (ZdcEtC <= 2.5 * 2.51);
                        bool TnTFn = (ZdcEtA > 1.5 * 2.51) && (ZdcEtA <= 2.5 * 2.51) && (ZdcEtC > 2.5 * 2.51) && (ZdcEtC <= 4.5 * 2.51);
                        bool TFnZn = (ZdcEtA > 2.5 * 2.51) && (ZdcEtA <= 4.5 * 2.51) && (ZdcEtC <= 1);
                        bool TFnOn = (ZdcEtA > 2.5 * 2.51) && (ZdcEtA <= 4.5 * 2.51) && (ZdcEtC > 1) && (ZdcEtC <= 1.5 * 2.51);
                        bool TFnTn = (ZdcEtA > 2.5 * 2.51) && (ZdcEtA <= 4.5 * 2.51) && (ZdcEtC > 1.5 * 2.51) && (ZdcEtC <= 2.5 * 2.51);
                        bool TFnTFn = (ZdcEtA > 2.5 * 2.51) && (ZdcEtA <= 4.5 * 2.51) && (ZdcEtC > 2.5 * 2.51) && (ZdcEtC <= 4.5 * 2.51);


                        bool ZDCetEn[4][4] = {
                            {ZnZn, ZnOn, ZnTn, ZnTFn},
                            {OnZn, OnOn, OnTn, OnTFn},
                            {TnZn, TnOn, TnTn, TnTFn},
                            {TFnZn, TFnOn, TFnTn, TFnTFn}
                        };

                        if ((ZDCetEn[Aen][Cen] || ZDCetEn[Cen][Aen]) && ((*vtx_ntrk)[0] == 2)) // to manage duplicates
                        {
                            // cout << "HEREEEEEEE" << endl;
                            int  NPairs = muon_pair_acop->size();
                            for (int j = 0; j < NPairs; j++)
                            {
                                float aco = -1.0;
                                float p1 = 0.0;
                                float p2 = 0.0;
                                float avg_p = 0.0;
                                float kperp = 0.0;
                                float e1 = 0.0;
                                float e2 = 0.0;
                                double y = 0.0;
                                float phi1 = 0.0;
                                float phi2 = 0.0;
                                float m = .1056583755; // muon mass in GeV

                                p1 = abs(muon_pair_muon1_pt->at(j) / 1000); //pt in MeV, scale to GeV
                                p2 = abs(muon_pair_muon2_pt->at(j) / 1000);
                                aco = abs(muon_pair_acop->at(j));
                                avg_p = (p1 + p2) / 2;
                                kperp = avg_p * aco * pi;
                                e1 = abs(muon_pair_muon1_eta->at(j));
                                e2 = abs(muon_pair_muon1_eta->at(j));
                                phi1 = muon_pair_muon1_phi->at(j);
                                phi2 = muon_pair_muon2_phi->at(j);
                                // pt eta phi m

                                ROOT::Math::PtEtaPhiMVector y1(p1, muon_pair_muon1_eta->at(j), phi1, m);
                                ROOT::Math::PtEtaPhiMVector y2(p2, muon_pair_muon2_eta->at(j), phi2, m);

                                y = abs((y1 + y2).Rapidity());

                                bool bin1 = (p1 > 3.75) && (p2 > 3.75) && (p1 <= 5) && (p2 <= 5);
                                bool bin2 = (p1 > 5) && (p2 > 5) && (p1 <= 10) && (p2 <= 10);
                                bool bin3 = (p1 > 10) && (p2 > 10) && (p1 <= 15) && (p2 <= 15);
                                bool bin4 = (p1 > 15) && (p2 > 15) && (p1 <= 20) && (p2 <= 20);
                                bool bin5 = (p1 > 20) && (p2 > 20) && (p1 <= 25) && (p2 <= 25);

                                bool PT_bins[6] = {
                                    bin1, bin2, bin3, bin4, bin5, 1
                                };

                                if (muon_pair_muon1_pt->at(j) * muon_pair_muon2_pt->at(j) <= 0)
                                {
                                    if ( PT_bins[c_ind] && ((e1 < 2.4) && (e2 < 2.4)) && (y < 2.4) )
                                    {
                                        histo->Fill(kperp);
                                    }
                                }
                                //

                            }
                        }

                    } //end of event loop



                } //end of compare

                histo->Sumw2();
                histo->Scale(1. / histo->Integral(), "width");
                histo->Draw();
                double po = histo->GetBinContent(1);
                // f4->SetParameters(po, P1, P2, P3, P4);
                // f4->SetNpx(500);
                // histo->Fit("f4");
                // if ((Aen == 2 && Cen == 3) || (Aen == 3 && Cen == 3))
                if ((Aen == 3 && Cen == 3))
                {
                    f5->SetParameters(po, P1, P2);
                }
                else
                {
                    f5->SetParameters(po, P1, P2, P3); //fixed P4
                    f5->FixParameter(4, P4);
                    // f5->ReleaseParameter(P4);
                    // f5->SetParLimits(2, 1e-6, .6);


                    // cout << "P4: " << P4 << endl;
                }

                f5->SetNpx(1500);
                // histo->Fit("f5");
                // histo->Fit("f5", "WL");

                // cout << "P4 from fit: " << f5->GetParameter(4) << endl;


                gPad->Modified(); gPad->Update();
                histo->SetLineColor(colors[0]);
                histo->SetTitle(leg_key2[Aen][Cen]);
                histo->GetXaxis()->SetTitle("K_perp");
                histo->GetYaxis()->SetTitle("Number of Events");
                float height = histo->GetBinContent(1);
                histo->GetYaxis()->SetRangeUser(0.01, height);
                leg->AddEntry(histo, "UPC Dimuon K_perp " + leg_key2[Aen][Cen]);
                if (c_ind == 0)
                {
                    leg->AddEntry(histo, "3.75 < Muon P_t < 5");
                }
                else if (c_ind == 1)
                {
                    leg->AddEntry(histo, "5 < Muon P_t < 10");
                }
                else if (c_ind == 2)
                {
                    leg->AddEntry(histo, "10 < Muon P_t < 15");
                }
                else if (c_ind == 3)
                {
                    leg->AddEntry(histo, "15 < Muon P_t < 20");
                }
                else if (c_ind == 4)
                {
                    leg->AddEntry(histo, "20 < Muon P_t < 25");
                }

                gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn
                leg->Draw();
                gPad->SetTicks(1, 1);

                TFitResultPtr result_ptr = histo->Fit(f5, "SIWL"); // WL = weighted log likelihood
                // TFitResultPtr result_ptr = histo->Fit(f5, "SI");
                float sig = f5->GetParameter(1);
                // if (sig < 5e-5)
                // {
                //     // histo_new = histo->Clone(hist);
                //     f5->FixParameter(1, 0);
                //     result_ptr = histo->Fit(f5, "SI");
                // }

                if (!result_ptr->IsValid())
                {
                    sigma_y[Aen][Cen] = 0.0;
                }

                TMatrixD cov = result_ptr->GetCovarianceMatrix();
                cov.Print();
                float cov_tau_sig = cov[1][2]; //covariance between sigma and tau? 5x5 matrix

                sigma_parameter[Aen][Cen] = f5->GetParameter(1); //
                tau_parameter[Aen][Cen] = f5->GetParameter(2);
                exp_slope_parameter[Aen][Cen] = 0; //f5->GetParameter(4);
                dsigma[Aen][Cen] = f5->GetParError(1); //sigma
                dtau[Aen][Cen] = f5->GetParError(2); //tau

                cout << "(dsigma, dtau) = (" << dsigma[Aen][Cen] << ", " << dtau[Aen][Cen] << ")" << endl;

                float error1_sq = 4 * pow((tau_parameter[Aen][Cen] * dtau[Aen][Cen]), 2);
                float error2_sq = pow((sigma_parameter[Aen][Cen] * dsigma[Aen][Cen]), 2);
                float cov_err = 4 * tau_parameter[Aen][Cen] * sigma_parameter[Aen][Cen] * cov_tau_sig;
                sigma_y[Aen][Cen] = sqrt(error1_sq + error2_sq + cov_err);
                // sigma_y[Aen][Cen] = sqrt(error1_sq + error2_sq);
                cout << "sigma_y: " << sigma_y[Aen][Cen] << endl << endl << endl;

                float prod1 = ((pow(sigma_parameter[Aen][Cen], 3) * tau_parameter[Aen][Cen]) - (sigma_parameter[Aen][Cen] * pow(tau_parameter[Aen][Cen], 3)) - pow(sigma_parameter[Aen][Cen], 2)) * erf(sigma_parameter[Aen][Cen] / (tau_parameter[Aen][Cen] * sqrt(2)));
                float prod2 = (sigma_parameter[Aen][Cen] * tau_parameter[Aen][Cen] / sqrt(2 * pi)) * exp(-1 * pow(sigma_parameter[Aen][Cen], 2) / (2 * pow(tau_parameter[Aen][Cen], 2)));
                float prod3 = pow(tau_parameter[Aen][Cen], 2) * exp(pow(sigma_parameter[Aen][Cen], 2) / (2 * pow(tau_parameter[Aen][Cen], 2)));

                mean[Aen][Cen] = prod1 + prod2 + prod3;

            }
        } // end of C side energy
    } // end of A side energy

    cout << "Finished with 4x4 matrix of plots " << endl << endl;

    TFile f("Kt_EMG.root", "RECREATE");
    Hlist.Write();
    f.Close();

}




























void UPC_Dimuon_Ana_4()
{
    create_chain();



    int ana_type = 0;

    cout << "\n\t                            ****  Analysis Types **** "                                << endl << endl;
    cout << "   1.\t    -                          ZdcEt Plots                                      - " << endl << endl;
    cout << "   2.\t    -                          K_perp w/ Muon_Pt Bin                            - " << endl << endl;
    cout << "   3.\t    -                          K_perp w/o Muon_Pt Bin                           - " << endl << endl;
    cout << "   4.\t    -                                                                           - " << endl << endl;

    cout << "Please choose an analysis option: ";
    cin >> ana_type;
    cout << "-----------------------------------" << endl << endl << endl;


    if (ana_type == 1)
    {

        c1 = new TCanvas("canvas1", "canvas1", 1200, 1000);
        c1->Divide(2, 1);

        ZdcEt_Plots("ZdcEt", 0);
        ZdcEt_Plots("ZdcEt", 1);
    }


    else if (ana_type == 2)
    {

        c0 = new TCanvas("canvas0", "canvas0", 800, 600);

        // c1 = new TCanvas("canvas1", "canvas1", 1600, 1400);
        c2 = new TCanvas("canvas2", "canvas2", 1600, 1400);
        // c3 = new TCanvas("canvas3", "canvas3", 1600, 1400);
        // c4 = new TCanvas("canvas4", "canvas4", 1600, 1400);
        // c5 = new TCanvas("canvas5", "canvas5", 1600, 1400);
        // c1->Divide(4, 4);
        // c1->Divide(4, 3);
        c2->Divide(4, 3);
        // c3->Divide(4, 3);
        // c4->Divide(4, 3);
        // c5->Divide(4, 3);

        // c6 = new TCanvas("canvas6", "canvas6", 800, 600);
        c7 = new TCanvas("canvas7", "canvas7", 800, 600);
        // c8 = new TCanvas("canvas8", "canvas8", 800, 600);
        // c9 = new TCanvas("canvas9", "canvas9", 800, 600);
        // c10 = new TCanvas("canvas10", "canvas10", 800, 600);


        // k_perp_plot_comp3_dup("K_perp", 0); //0-4
        // k_perp_Tgraph("K_perp EMG RMS Values", 0, 0); // 10 pt unique RMS_reg
        k_perp_plot_comp3_dup("K_perp", 1);
        k_perp_Tgraph("K_perp EMG RMS Values", 0, 1);
        // k_perp_plot_comp3_dup("K_perp", 2);
        // k_perp_Tgraph("K_perp EMG RMS Values", 0, 2);
        // k_perp_plot_comp3_dup("K_perp", 3);
        // k_perp_Tgraph("K_perp EMG RMS Values", 0, 3);
        // k_perp_plot_comp3_dup("K_perp", 4);
        // k_perp_Tgraph("K_perp EMG RMS Values", 0, 4);

    }

    else if (ana_type == 3)
    {
        c0 = new TCanvas("canvas0", "canvas0", 800, 600);
        c1 = new TCanvas("canvas1", "canvas1", 1600, 1400);
        c1->Divide(4, 3);
        c6 = new TCanvas("canvas6", "canvas6", 800, 600);

        k_perp_plot_comp3_dup("K_perp", 5); //0-5
        k_perp_Tgraph("K_perp EMG RMS Values", 0, 5); // 10 pt unique RMS_reg


    }

    else if (ana_type == 4)
    {
        return;
    }


    else
    {
        cout << "Invalid input; Please try again." << endl;
    }


}
