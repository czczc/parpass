// convert momentum to KE
double KE(double mom, double mass)
{
    double bg = mom / mass;           // beta*gamma.
    double gamma = sqrt(1. + bg*bg);  // gamma.
    return (gamma-1) * mass;  // KE in GeV
}

// calculate range given initial E0 and dx/dE curve
double range(double E0, TGraph *f)
{
    // range = \int dx/dE * dE
    double dE = 0.1; // MeV
    double Emin = 0.1;  // MeV
    double E = E0;
    double range = 0;
    while (E>Emin) {
        double dx = f->Eval(E) * dE;
        range += dx;
        E -= dE;
        // cout << "E: " << E << " dx: " << dx << endl;
    }
    return range;
}

void SetStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleOffset(1.0, "x");
    gStyle->SetTitleOffset(1.0, "y");
    gStyle->SetTitleFont(42, "P");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);
    gStyle->SetLabelSize(0.05, "x");
    gStyle->SetLabelSize(0.05, "y");
    gStyle->SetHistLineWidth(2);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetMarkerSize(0.3);
}

void plot_range()
{
    SetStyle();
    LArProperties prop = LArProperties();
    prop.PrintInfo();

    TDatabasePDG *dbPDG = new TDatabasePDG();
    TParticlePDG *proton = dbPDG->GetParticle(2212);
    TParticlePDG *kaon = dbPDG->GetParticle(321);
    TParticlePDG *pion = dbPDG->GetParticle(211);
    TParticlePDG *muon = dbPDG->GetParticle(13);

    cout << proton->GetName() << " mass: " << proton->Mass() << " GeV" << endl;
    cout << kaon->GetName() << " mass: " << kaon->Mass() << " GeV" << endl;
    cout << pion->GetName() << " mass: " << pion->Mass() << " GeV" << endl;
    cout << muon->GetName() << " mass: " << muon->Mass() << " GeV" << endl;

    const int N=10000;
    double p[N];
    double x_proton[N], x_kaon[N], x_pion[N], x_muon[N];
    double y_proton[N], y_kaon[N], y_pion[N], y_muon[N];
    double start = 0.; //GeV
    double step = 1e-4; //GeV
    double tcut = 0;
    for (int i=0; i!=N; i++) {
        p[i] = start + step*i;
        x_proton[i] = KE(p[i], proton->Mass()) * 1000; // KE in MeV
          x_kaon[i] = KE(p[i], kaon->Mass()) * 1000; // KE in MeV
          x_pion[i] = KE(p[i], pion->Mass()) * 1000; // KE in MeV
          x_muon[i] = KE(p[i], muon->Mass()) * 1000; // KE in MeV
        //   cout << "p: " << p[i] << " , " << "KE: " << x_muon[i] << endl; 
        y_proton[i] = 1./prop.Eloss(p[i], proton->Mass(), tcut);
          y_kaon[i] = 1./prop.Eloss(p[i], kaon->Mass(), tcut);
          y_pion[i] = 1./prop.Eloss(p[i], pion->Mass(), tcut);
          y_muon[i] = 1./prop.Eloss(p[i], muon->Mass(), tcut);
    }
    // dx/dE curve as a function of KE.
    TGraph *f_proton = new TGraph(N, x_proton, y_proton);
    TGraph *f_kaon = new TGraph(N, x_kaon, y_kaon);
    TGraph *f_pion = new TGraph(N, x_pion, y_pion);
    TGraph *f_muon = new TGraph(N, x_muon, y_muon);

    // const int M = 10000+1;
    const int M = 2000+1;

    double e[M];
    double r_proton[M], r_kaon[M], r_pion[M], r_muon[M];
    e[0] = r_proton[0] = r_kaon[0] = r_pion[0] = r_muon[0] = 0;
    double dE = 0.1; // MeV
    for (int i=1; i<M; i++) {
        e[i] = i * dE;
        r_proton[i] = r_proton[i-1] + f_proton->Eval(e[i]) * dE; // cm
	r_kaon[i] =   r_kaon[i-1] +   f_kaon->Eval(e[i]) * dE;
	r_pion[i] =   r_pion[i-1] +   f_pion->Eval(e[i]) * dE;
	r_muon[i] =   r_muon[i-1] +   f_muon->Eval(e[i]) * dE;
        // r_proton[i] = range(e[i], f_proton);
        // cout << "proton E: " << e[i] << " , range: " << r_proton[i] << endl;

    }

    TH2F *h2 = new TH2F("h2","", 100, 0, r_muon[M-1]+1, 100, 0, e[M-1]+dE);
    h2->GetYaxis()->SetTitle("Kinetic Energy [MeV]");
    h2->GetXaxis()->SetTitle("Range [cm]");
    h2->Draw();

    TGraph *g_proton = new TGraph(M, r_proton, e);
    g_proton->SetLineWidth(2);
    g_proton->SetLineColor(kRed);
    g_proton->Draw("L,same");

    TGraph *g_kaon = new TGraph(M, r_kaon, e);
    g_kaon->SetLineWidth(2);
    g_kaon->SetLineColor(kCyan);
    g_kaon->Draw("L,same");

    TGraph *g_pion = new TGraph(M, r_pion, e);
    g_pion->SetLineWidth(2);
    g_pion->SetLineColor(kViolet);
    g_pion->Draw("L,same");

    TGraph *g_muon = new TGraph(M, r_muon, e);
    g_muon->SetLineWidth(2);
    g_muon->SetLineColor(kGreen);
    g_muon->Draw("L,same");

    TLegend *leg = new TLegend(0.60,0.20, 0.87,0.47);
    leg->AddEntry(g_proton, " proton", "l");
    leg->AddEntry(g_kaon,   " kaon", "l");
    leg->AddEntry(g_pion,   " pion", "l");
    leg->AddEntry(g_muon,   " muon", "l");
    leg->SetFillColor(kWhite);
    leg->Draw();

    gPad->SetGridx();
    gPad->SetGridy();

    // dedx as function of residual range 
    TCanvas *c2 = new TCanvas;
    const int M2 = 1000;
    double rx[M2];
    double dedx_proton[M2], dedx_kaon[M2], dedx_pion[M2], dedx_muon[M2], dedx_mip[M2];
    // rx[0] = dedx_proton[0] = dedx_kaon[0] = dedx_pion[0] = dedx_muon[0] = 0;
    double dx = 0.05; // cm
    for (int i=0; i<M2; i++) {
        rx[i] = (i+1) * dx;
        dedx_proton[i] =  1. / f_proton->Eval( g_proton->Eval(rx[i]) ); // MeV/cm
          dedx_kaon[i] =  1. /   f_kaon->Eval(   g_kaon->Eval(rx[i]) );
          dedx_pion[i] =  1. /   f_pion->Eval(   g_pion->Eval(rx[i]) );
          dedx_muon[i] =  1. /   f_muon->Eval(   g_muon->Eval(rx[i]) );
        //   cout << "rx: " << rx[i] << " , proton dedx: " << dedx_proton[i] << endl;

    }

    TH2F *h22 = new TH2F("h22","", 100, 0, 50, 100, 0, 20);
    h22->GetXaxis()->SetTitle("Residual Range [cm]");
    h22->GetYaxis()->SetTitle("dE/dx (MeV/cm) ");
    h22->Draw();

    TGraph *g2_proton = new TGraph(M2, rx, dedx_proton);
    g2_proton->SetLineWidth(2);
    g2_proton->SetLineColor(kRed);
    g2_proton->Draw("L,same");

    TGraph *g2_kaon = new TGraph(M2, rx, dedx_kaon);
    g2_kaon->SetLineWidth(2);
    g2_kaon->SetLineColor(kCyan);
    g2_kaon->Draw("L,same");


    TGraph *g2_pion = new TGraph(M2, rx, dedx_pion);
    g2_pion->SetLineWidth(2);
    g2_pion->SetLineColor(kViolet);
    g2_pion->Draw("L,same");

    TGraph *g2_muon = new TGraph(M2, rx, dedx_muon);
    g2_muon->SetLineWidth(2);
    g2_muon->SetLineColor(kGreen);
    g2_muon->Draw("L,same");


    TLegend *leg2 = new TLegend(0.65,0.60, 0.87,0.87);
    leg2->AddEntry(g2_proton, " proton", "l");
    leg2->AddEntry(g2_kaon,   " kaon", "l");
    leg2->AddEntry(g2_pion,   " pion", "l");
    leg2->AddEntry(g2_muon,   " muon", "l");
    leg2->SetFillColor(kWhite);
    leg2->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
}


