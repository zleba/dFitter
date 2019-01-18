#include "dplotter.h"

#include "plottingHelper.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include "TStyle.h"

#include <vector>
#include <map>




using namespace PlottingHelper;//pollute the namespace!
using namespace std;

TString rn() {return Form("%d",rand());}


void dPlotter::plotBeta(double xpom)
{
    map<double, vector<point>> dataMap;

    for(point &p : data)
        if(p.xp == xpom)
            dataMap[p.q2].push_back(p);

    int nRows = ceil(dataMap.size() / 4);
    //cout << dataMap.size() << endl;

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 4*200/0.8, nRows*200/0.8);
    SetLeftRight(0.1, 0.1);
    SetTopBottom(0.1, 0.1);

    DivideTransparent(group(1, 0, 4), group(1, 0, nRows));

    int i = 0;
    for(auto item : dataMap) {
        double q2    = item.first;
        auto &points = item.second;
        can->cd(i+1);

        cout << "Radek " << i << endl;
        TGraphErrors *gData = new TGraphErrors(points.size());

        for(int j = 0; j < points.size(); ++j) {
            gData->SetPoint(j, points[j].beta, points[j].xpSig);

            double er = points[j].xpSig * points[j].errTot;
            gData->SetPointError(j, 0, er);
        }

        map<double,double> zMin = { {0.0003, 0.02},
                                    {0.001, 0.015},
                                    {0.003, 0.0015},
                                    {0.01, 0.0015},
                                    {0.03, 0.0015}
                                   };

        TH1D *hFr = new TH1D(rn(), "", 1, zMin[xpom], 1);
        hFr->Draw("axis");
        
        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");

        TGraph *gTh = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gTh->SetPoint(j, points[j].beta, points[j].th);
        }
        gTh->SetLineColor(kBlue);
        gTh->Draw("l");

        DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2));

        gPad->SetLogx();

        GetYaxis()->SetNdivisions(503);

        GetYaxis()->SetRangeUser(0, 0.085);

        SetFTO({24}, {14}, {1.4, 2.2, 0.4, 3.9});

        if(i < dataMap.size() - 4) {
            GetXaxis()->SetTickSize(0);
            GetXaxis()->SetLabelOffset(5000);
        }

        if(i % 4 != 0) {
            GetYaxis()->SetLabelOffset(5000);
            GetYaxis()->SetTickSize(0);
        }

        if(i == 0)
            GetYaxis()->SetTitle("x_{IP} #sigma_{r}^{D(3)}");
        if(i == dataMap.size() -1)
            GetXaxis()->SetTitle("#beta");

        ++i;
    }
    can->SaveAs(Form("plots/xpom%g.pdf", xpom));

}
