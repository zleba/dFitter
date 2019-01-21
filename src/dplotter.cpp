#include "dplotter.h"

#include "plottingHelper.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include "TStyle.h"

#include <vector>
#include <map>

#include "pdf.h"


using namespace PlottingHelper;//pollute the namespace!
using namespace std;

TString rn() {return Form("%d",rand());}


void dPlotter::plotBeta(double xpom)
{
    map<double, vector<point>> dataMap;

    for(point &p : data)
        if(p.xp == xpom)
            dataMap[p.q2].push_back(p);

    int nRows = ceil(dataMap.size() / 4.);
    //cout << dataMap.size() << endl;

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", (1+4)*200, (1+nRows)*200);
    SetLeftRight(0.5/(1+4), 0.5/(1+4));
    SetTopBottom(0.5/(1+nRows), 0.5/(1+nRows));

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
        

        TGraph *gThOrg = new TGraph(points.size());
        TGraph *gTh    = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gThOrg->SetPoint(j, points[j].beta, points[j].thOrg);
            gTh->SetPoint(j, points[j].beta, points[j].th);
        }
        gThOrg->SetLineColor(kBlue);
        gThOrg->SetMarkerColor(kBlue);
        gThOrg->Draw("l* same");

        //MyTheory
        gTh->SetLineColor(kGreen);
        gTh->SetMarkerColor(kGreen);
        gTh->Draw("l* same");


        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");



        gPad->SetLogx();

        GetYaxis()->SetNdivisions(503);

        GetYaxis()->SetRangeUser(0, 0.085);

        SetFTO({24}, {14}, {1.4, 2.2, 0.4, 3.9});

        DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2));


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

    DrawLatexUp(can->GetPad(1), 1.1, Form("x_{IP} = %g", xpom), -1, "l");

    can->SaveAs(Form("dPlots/xpom%g.pdf", xpom));

}

//To print the curve tag
pair<double,double> getRight(TGraphErrors *gr)
{
    pair<double,double> pRight(-100,0);
    for(int i = 0; i < gr->GetN(); ++i) {
        double x, y;
        gr->GetPoint(i, x, y);
        if(pRight.first < x) {
            pRight = {x,y};
        }
    }
    return pRight;
}


void dPlotter::plotQ2(double xpom)
{
    map<double, vector<point>> dataMap;

    for(point &p : data)
        if(p.xp == xpom)
            dataMap[p.beta].push_back(p);

    int nPoints = dataMap.size();

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 500, 800);
    SetLeftRight(0.15, 0.07);
    SetTopBottom(0.1, 0.1);

    gPad->SetLogx();
    gPad->SetLogy();


    /*
       map<double,double> zMin = { {0.0003, 0.02},
       {0.001, 0.015},
       {0.003, 0.0015},
       {0.01, 0.0015},
       {0.03, 0.0015}
       };
       */


    //Plotting style
    TH1D *hFr = new TH1D(rn(), "", 1, 2, 1e4);
    hFr->Draw("axis");

    SetFTO({16}, {14}, {1.4, 2.2, 0.4, 3.9});
    GetYaxis()->SetRangeUser(1e-2, 1e5);
    GetYaxis()->SetTitle("3^{i} * x_{IP} #sigma_{r}^{D(3)}");
    GetXaxis()->SetTitle("Q^{2} [GeV]");



    int i = 0;
    for(auto item : dataMap) {
        double beta  = item.first;
        auto &points = item.second;

        double fac = pow(3, dataMap.size() - 1 - i);
        cout << "Radek " << i << endl;
        TGraphErrors *gData = new TGraphErrors(points.size());

        for(int j = 0; j < points.size(); ++j) {
            gData->SetPoint(j, points[j].q2, fac*points[j].xpSig);

            double er = fac*points[j].xpSig * points[j].errTot;
            gData->SetPointError(j, 0, er);
        }
        

        TGraph *gTh = new TGraph(points.size());
        TGraph *gThOrg = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gTh->SetPoint(j, points[j].q2, fac*points[j].th);
            gThOrg->SetPoint(j, points[j].q2, fac*points[j].thOrg);
        }

        gThOrg->SetLineColor(kBlue);
        gThOrg->SetMarkerColor(kBlue);
        gThOrg->Draw("l* same");

        gTh->SetLineColor(kGreen);
        gTh->SetMarkerColor(kGreen);
        gTh->Draw("l* same");

        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");



        double xT, yT;
        tie(xT,yT) =  getRight(gData);
        TLatex *lat = new TLatex;
        lat->SetTextSize(PxFontToRel(16));
        lat->SetTextAlign(12);
        lat->DrawLatex(xT*1.3, yT, Form("#beta=%g (i=%d)", beta, dataMap.size() - 1 - i));


        ++i;
    }

    DrawLatexUp(1, Form("x_{IP} = %g",xpom), -1, "l");


    can->SaveAs(Form("dPlots/q2%g.pdf", xpom));

}


void dPlotter::plotXpom()
{
    map<pair<double,double>, vector<point>> dataMap;

    //vector<double> betas = {0.01, 0.04, 0.1, 0.2, 0.4, 0.65, 0.9};
    vector<double> betas = { 0.011, 0.043, 0.11, 0.2, 0.43, 0.67, 0.8};

    vector<double> q2s   = {3.5, 5, 6.5, 8.5, 12, 15, 20, 25, 35, 45, 60, 90, 200, 400, 800, 1600};


    for(point &p : data)
        dataMap[{p.q2, p.beta}].push_back(p);


    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 800);
    SetLeftRight(0.1, 0.1);
    SetTopBottom(0.1, 0.1);

    DivideTransparent(group(1, 0, betas.size()), group(1, 0, q2s.size()));

    for(int iq = 0; iq < q2s.size();   ++iq)
    for(int ib = 0; ib < betas.size(); ++ib) {
        double q2   = q2s[iq];
        double beta = betas[ib];

        can->cd(iq*betas.size() + ib + 1);

        TH1D *hFr = new TH1D(rn(), "", 1, 1e-4, 4e-2);
        hFr->Draw("axis");
        GetYaxis()->SetRangeUser(0, 0.085);

        gPad->SetLogx();
        GetYaxis()->SetNdivisions(303);
        SetFTO({14}, {6}, {1.4, 2.2, 0.4, 3.9});


        if(iq != q2s.size() -1)
            GetXaxis()->SetTickSize(0);
        if(ib != 0)
            GetYaxis()->SetTickSize(0);



        if(!dataMap.count({q2,beta})) {
            cout << "Not existing " << q2 <<" "<< beta << endl;
            continue;
        }
        else {
            cout << "Yes existing " << q2 <<" "<< beta << endl;
        }
        auto &points = dataMap.at({q2,beta});


        TGraphErrors *gData = new TGraphErrors(points.size());

        for(int j = 0; j < points.size(); ++j) {
            gData->SetPoint(j, points[j].xp, points[j].xpSig);

            double er = points[j].xpSig * points[j].errTot;
            gData->SetPointError(j, 0, er);
        }


        
        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerSize(0.3);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");

        TGraph *gTh = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gTh->SetPoint(j, points[j].xp, points[j].th);
        }
        gTh->SetLineColor(kBlue);
        gTh->SetMarkerColor(kBlue);
        gTh->Draw("*l same");




        DrawLatex(0.5, 0.5, Form("%d", points.size()));

    }

    //DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2));


    can->SaveAs("dPlots/xpomGrid.pdf");

}





void dPlotter::plotPDFs(bool inLog)
{
    vector<double> q2s = {1.75, 8.5, 20, 90, 800};


    //q0^2 = 1.75
    static PDF myFitA;
    if(myFitA.ag == -99) { 
        myFitA.setPars(0.14591, 0, -0.94705,    1.0587, 2.2964, 0.56894);
        myFitA.evolve();
    }

    //myFitA.checkConvolution();
    //exit(0);

    /*
    myFitA.checkSumRulesInput();
    myFitA.checkSumRules();
    return;
    */

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.1, 0.14);
    SetTopBottom(0.1, 0.1);

    double zMin = 4e-3;

    DivideTransparent(group(1, 0.5, 2), group(1, 0, q2s.size()));

    for(int i = 0; i < q2s.size(); ++i) {

        //Fill Graph
        TGraph *grS = new TGraph();
        TGraph *grG = new TGraph();

        TGraph *grSmy = new TGraph();
        TGraph *grGmy = new TGraph();


        for(double z = zMin; z < 1; z *= 1.001) {
        //for(double z = zMin; z < 1; z += 0.001) {
            double g, qS;
            tie(g, qS) =  PDF::evalFitA(z, q2s[i]);
            //tie(g, qS) =  myFitA.eval(z, q2s[i]);
            qS *= 6;

            grS->SetPoint(grS->GetN(), z, qS);
            grG->SetPoint(grG->GetN(), z, g);

            double gMy, qSMy;
            tie(gMy, qSMy) = myFitA.eval(z, q2s[i]);
            grSmy->SetPoint(grSmy->GetN(), z, (qSMy*6));
            grGmy->SetPoint(grGmy->GetN(), z, gMy );

        }
        grS->SetLineColor(kBlue);
        grG->SetLineColor(kBlue);


        can->cd(2*i + 1);
        TH1D *hFrS = new TH1D(rn(), "", 1, zMin, 1);
        hFrS->Draw("axis");

        grS->Draw("l same");
        grSmy->Draw("l same");
        GetYaxis()->SetRangeUser(0, 0.27);
        //GetYaxis()->SetRangeUser(0.9, 1.1);
        GetYaxis()->SetNdivisions(503);
        GetXaxis()->SetNdivisions(404);
        SetFTO({15}, {6}, {1.4, 2.2, 0.4, 3.9});

        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Singlet");
            GetYaxis()->SetTitle("z #Sigma(z,Q^{2})");
        }
        if(i == q2s.size()-1) {
            GetXaxis()->SetTitle("z");
        }
        else {
            GetXaxis()->SetLabelOffset(50000);
        }

        can->cd(2*i + 2);
        TH1D *hFrG = new TH1D(rn(), "", 1, zMin, 1);
        hFrG->Draw("axis");
        grG->Draw("l same");
        grGmy->Draw("l same");
        GetYaxis()->SetRangeUser(0, 2.25);
        //GetYaxis()->SetRangeUser(0.9, 1.1);
        GetYaxis()->SetNdivisions(503);
        GetXaxis()->SetNdivisions(404);
        SetFTO({15}, {6}, {1.4, 2.2, 0.4, 3.9});
        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Gluon");
            GetYaxis()->SetTitle("z g(z,Q^{2})");
        }
        if(i == q2s.size()-1) {
            GetXaxis()->SetTitle("z");
        }
        else {
            GetXaxis()->SetLabelOffset(50000);
        }

        DrawLatexRight(1, Form("Q^{2}=%g",q2s[i]), -1, "l");

        if(i == 3) {
            auto *leg = newLegend(kPos9);
            leg->AddEntry(grG,   "Published");
            leg->AddEntry(grGmy, "Recalculated");
            DrawLegends({leg});
        }
    }

    if(inLog)
        can->SaveAs("dPlots/pdfsLog.pdf");
    else
        can->SaveAs("dPlots/pdfsLin.pdf");

}



void dPlotter::plotF2FLcharm(bool inLog)
{
    vector<double> q2s = {1.75, 8.5, 20, 90, 800};
    vector<double> rangeF2 = {0.004, 0.02,  0.03, 0.06, 0.13};
    vector<double> rangeFL = {0.0002, 0.002,  0.006, 0.013, 0.025};


    //q0^2 = 1.75
    static PDF myFitA;
    if(myFitA.ag == -99) { 
        myFitA.setPars(0.14591, 0, -0.94705,    1.0587, 2.2964, 0.56894);
        myFitA.evolve();
        myFitA.initConv();
    }

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.1, 0.14);
    SetTopBottom(0.1, 0.1);

    double zMin = 4e-3;

    DivideTransparent(group(1, 0.5, 2), group(1, 0, q2s.size()));

    for(int i = 0; i < q2s.size(); ++i) {

        //Fill Graph
        TGraph *grF2c = new TGraph(); //charm components
        TGraph *grFLc = new TGraph(); 

        TGraph *grF2cmy = new TGraph();//charm components
        TGraph *grFLcmy = new TGraph();

        TGraph *grF2bmy = new TGraph();//bottom components
        TGraph *grFLbmy = new TGraph();


        //Loop over z points
        for(double z = zMin; z < 1; z *= 1.1) {

            //Original -- Charm
            double f2c, flc;
            tie(f2c, flc) = PDF::getF2FL(z, q2s[i], 'C');
            grF2c->SetPoint(grF2c->GetN(), z, f2c);
            grFLc->SetPoint(grFLc->GetN(), z, flc);

            //My -- Charm
            double f2cm, flcm;
            tie(f2cm, flcm) = myFitA.getF2FLmy(z, q2s[i], 'C');
            grF2cmy->SetPoint(grF2cmy->GetN(), z, f2cm);
            grFLcmy->SetPoint(grFLcmy->GetN(), z, flcm);

            //My -- Bottom
            double f2bm, flbm;
            tie(f2bm, flbm) = myFitA.getF2FLmy(z, q2s[i], 'B');
            grF2bmy->SetPoint(grF2bmy->GetN(), z, f2bm);
            grFLbmy->SetPoint(grFLbmy->GetN(), z, flbm);
        }



        grF2c->SetLineColor(kRed);
        grFLc->SetLineColor(kRed);
        grF2cmy->SetLineColor(kBlue);
        grFLcmy->SetLineColor(kBlue);
        grF2bmy->SetLineColor(kBlue);
        grFLbmy->SetLineColor(kBlue);

        grF2cmy->SetLineStyle(2);
        grFLcmy->SetLineStyle(2);
        grF2bmy->SetLineStyle(2);
        grFLbmy->SetLineStyle(2);

        can->cd(2*i + 1);
        TH1D *hFrS = new TH1D(rn(), "", 1, zMin, 1);
        hFrS->Draw("axis");

        grF2c->Draw("l same");
        grF2cmy->Draw("l same");
        grF2bmy->Draw("l same");

        GetYaxis()->SetRangeUser(0, rangeF2[i]);
        GetYaxis()->SetNdivisions(503);
        SetFTO({15}, {6}, {1.2, 2.2, 0.4, 4.3});

        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Charm F_{2}");
            GetYaxis()->SetTitle("F_{2}");
        }
        if(i == q2s.size()-1) {
            GetXaxis()->SetTitle("z");
        }
        else {
            GetXaxis()->SetLabelOffset(50000);
        }


        can->cd(2*i + 2);
        TH1D *hFrG = new TH1D(rn(), "", 1, zMin, 1);
        hFrG->Draw("axis");

        grFLc->Draw("l same");
        grFLcmy->Draw("l same");
        grFLbmy->Draw("l same");


        GetYaxis()->SetRangeUser(0, rangeFL[i]);
        GetYaxis()->SetNdivisions(503);
        SetFTO({15}, {6}, {1.2, 2.2, 0.4, 4.3});
        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Charm F_{L}");
            GetYaxis()->SetTitle("F_{L}");
        }
        if(i == q2s.size()-1) {
            GetXaxis()->SetTitle("z");
        }
        else {
            GetXaxis()->SetLabelOffset(50000);
        }

        DrawLatexRight(1, Form("Q^{2}=%g",q2s[i]), -1, "l");

        if(i == 3) {
            auto *leg = newLegend(kPos9);
            leg->AddEntry(grFLc,   "Published");
            leg->AddEntry(grFLcmy, "Recalculated");
            DrawLegends({leg});
        }



    }

    if(inLog)
        can->SaveAs("dPlots/f2flCharmLog.pdf");
    else
        can->SaveAs("dPlots/f2flCharmLin.pdf");

}


void dPlotter::plotF2FL(bool inLog)
{
    vector<double> q2s = {1.75, 8.5, 20, 90, 800};
    vector<double> rangeF2 = {0.06, 0.06,  0.1, 0.16, 0.3};
    vector<double> rangeFL = {0.012, 0.02,  0.03, 0.05, 0.05};



    //q0^2 = 1.75
    static PDF myFitA;
    if(myFitA.ag == -99) { 
        myFitA.setPars(0.14591, 0, -0.94705,    1.0587, 2.2964, 0.56894);
        myFitA.evolve();
        myFitA.initConv();
    }

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.1, 0.14);
    SetTopBottom(0.1, 0.1);

    double zMin = 4e-3;

    DivideTransparent(group(1, 0.5, 2), group(1, 0, q2s.size()));

    for(int i = 0; i < q2s.size(); ++i) {

        //Fill Graph
        TGraph *grF2  = new TGraph();
        TGraph *grFL  = new TGraph();
        TGraph *grF2c = new TGraph(); //charm components
        TGraph *grFLc = new TGraph(); 

        TGraph *grF2my  = new TGraph();//all components
        TGraph *grFLmy  = new TGraph();
        TGraph *grF2cmy = new TGraph();//charm components
        TGraph *grFLcmy = new TGraph();
        TGraph *grF2clmy = new TGraph();//charm+light components
        TGraph *grFLclmy = new TGraph();


        //Loop over z points
        for(double z = zMin; z < 1; z *= 1.1) {

            //cout << "Helenka " << z << " "<< grF2->GetN() << endl;
            //Original -- All
            double f2, fl;
            tie(f2, fl) = PDF::getF2FL(z, q2s[i], 'A');
            grF2->SetPoint(grF2->GetN(), z, f2);
            grFL->SetPoint(grFL->GetN(), z, fl);

            //Original -- Charm
            double f2c, flc;
            tie(f2c, flc) = PDF::getF2FL(z, q2s[i], 'C');
            grF2c->SetPoint(grF2c->GetN(), z, f2c);
            grFLc->SetPoint(grFLc->GetN(), z, flc);



            //My -- Charm
            double f2cm, flcm;
            tie(f2cm, flcm) = myFitA.getF2FLmy(z, q2s[i], 'C');
            grF2cmy->SetPoint(grF2cmy->GetN(), z, f2cm);
            grFLcmy->SetPoint(grFLcmy->GetN(), z, flcm);

            //cout << "Helenka " << q2s[i] <<" "<< z <<" : " <<" "<<f2<<" "<<   f2c <<" | "<< f2cm << endl;

            //My -- Charm + Light
            double f2lm, fllm;
            tie(f2lm, fllm) = myFitA.getF2FLmy(z, q2s[i], 'L');
            grF2clmy->SetPoint(grF2clmy->GetN(), z, f2cm  + f2lm);
            grFLclmy->SetPoint(grFLclmy->GetN(), z, flcm + fllm);

            //My -- Charm + Light + bottom
            double f2bm, flbm;
            tie(f2bm, flbm) = myFitA.getF2FLmy(z, q2s[i], 'B');
            grF2my->SetPoint(grF2my->GetN(), z, f2cm  + f2lm + f2bm);
            grFLmy->SetPoint(grFLmy->GetN(), z, flcm + fllm  + flbm);
        }


        grF2->SetLineColor(kRed);
        grFL->SetLineColor(kRed);
        grF2my->SetLineColor(kBlue) ;
        grFLmy->SetLineColor(kBlue);

        grF2c->SetLineColor(kRed);
        grFLc->SetLineColor(kRed);
        grF2cmy->SetLineColor(kBlue) ;
        grFLcmy->SetLineColor(kBlue);




        grF2my->SetLineStyle(2) ;
        grFLmy->SetLineStyle(2);
        grF2cmy->SetLineStyle(2) ;
        grFLcmy->SetLineStyle(2);


        can->cd(2*i + 1);
        TH1D *hFrS = new TH1D(rn(), "", 1, zMin, 1);
        hFrS->Draw("axis");

        grF2->Draw("l same");
        grF2my->Draw("l same");
        grF2c->Draw("l same");
        grF2cmy->Draw("l same");

        GetYaxis()->SetRangeUser(0, rangeF2[i]);
        GetYaxis()->SetNdivisions(503);
        SetFTO({15}, {6}, {1.2, 2.2, 0.4, 4.3});

        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Inclusive F_{2}");
            GetYaxis()->SetTitle("F_{2}");
        }
        if(i == q2s.size()-1) {
            GetXaxis()->SetTitle("z");
        }
        else {
            GetXaxis()->SetLabelOffset(50000);
        }


        can->cd(2*i + 2);
        TH1D *hFrG = new TH1D(rn(), "", 1, zMin, 1);
        hFrG->Draw("axis");

        grFL->Draw("l same");
        grFLmy->Draw("l same");
        grFLc->Draw("l same");
        grFLcmy->Draw("l same");


        GetYaxis()->SetRangeUser(0, rangeFL[i]);
        //GetYaxis()->SetRangeUser(0.9, 1.1);
        GetYaxis()->SetNdivisions(503);
        SetFTO({15}, {6}, {1.2, 2.2, 0.4, 4.3});
        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Inclusive F_{L}");
            GetYaxis()->SetTitle("F_{L}");
        }
        if(i == q2s.size()-1) {
            GetXaxis()->SetTitle("z");
        }
        else {
            GetXaxis()->SetLabelOffset(50000);
        }

        DrawLatexRight(1, Form("Q^{2}=%g",q2s[i]), -1, "l");

        if(i == 3) {
            auto *leg = newLegend(kPos9);
            leg->AddEntry(grFL,   "Published");
            leg->AddEntry(grFLmy, "Recalculated");
            DrawLegends({leg});
        }


    }

    if(inLog)
        can->SaveAs("dPlots/f2flLog.pdf");
    else
        can->SaveAs("dPlots/f2flLin.pdf");

}

