// #include "SiProperties.h"

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

void plot_time()
{
    SetStyle();
    SiProperties prop = SiProperties();
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
    double step = 1e-5; //GeV
    double tcut = 0;
    for (int i=0; i!=N; i++) {
        p[i] = start + step*i;
        x_proton[i] = KE(p[i], proton->Mass()) * 1000; // KE in MeV
	x_kaon[i] = KE(p[i], kaon->Mass()) * 1000; // KE in MeV
	x_pion[i] = KE(p[i], pion->Mass()) * 1000; // KE in MeV
	x_muon[i] = KE(p[i], muon->Mass()) * 1000; // KE in MeV

	if (i==0){
	  y_proton[i] = 0;
	  y_kaon[i] = 0.;
	  y_pion[i] = 0.;
	  y_muon[i] = 0.;
	}else{
	  y_proton[i] = 1./prop.Eloss(p[i], proton->Mass(), tcut);
	  y_kaon[i] = 1./prop.Eloss(p[i], kaon->Mass(), tcut);
	  y_pion[i] = 1./prop.Eloss(p[i], pion->Mass(), tcut);
	  y_muon[i] = 1./prop.Eloss(p[i], muon->Mass(), tcut);
	}
	//  cout << "p: " << p[i] << " , " << "KE: " << x_muon[i] << " " << y_muon[i] << endl;
	//  break;
    }
    // dx/dE curve as a function of KE.
    TGraph *f_proton = new TGraph(N, x_proton, y_proton);
    TGraph *f_kaon = new TGraph(N, x_kaon, y_kaon);
    TGraph *f_pion = new TGraph(N, x_pion, y_pion);
    TGraph *f_muon = new TGraph(N, x_muon, y_muon);


    // f_muon->Draw("AL");

    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    c1->Divide(2,2);
    c1->cd(1);
    
    {
      // const int M = 10000+1;
      const int M = 1000+1;
      
      double e[M];
      double r_proton[M], r_kaon[M], r_pion[M], r_muon[M];
      e[0] = r_proton[0] = r_kaon[0] = r_pion[0] = r_muon[0] = 0;
      double dE = 4.12/M; // MeV
      double t_muon[M], t_muon_ar[M];
      t_muon[0] = 0.; t_muon_ar[0] = 0;
      for (int i=1; i<M; i++) {
        e[i] = i * dE;
        r_proton[i] = r_proton[i-1] + f_proton->Eval(e[i]) * dE * 1e4; // um
    	r_kaon[i] =   r_kaon[i-1] +   f_kaon->Eval(e[i]) * dE * 1e4;
    	r_pion[i] =   r_pion[i-1] +   f_pion->Eval(e[i]) * dE * 1e4;
    	r_muon[i] =   r_muon[i-1] +   f_muon->Eval(e[i]) * dE * 1e4;
	
	// muon velocity
	double muon_vel = sqrt(1-pow(muon->Mass()*1000.,2)/pow(e[i] +muon->Mass()*1000. ,2))*3e5; // um/ns ...
	t_muon[i] = t_muon[i-1] + f_muon->Eval(e[i]) * dE * 1e4/muon_vel;
	t_muon_ar[i] = t_muon_ar[i-1] + f_muon->Eval(e[i]) * dE * 1e4/muon_vel * pow(muon->Mass()*1000.,2)/pow(e[i] +muon->Mass()*1000. ,2);
	
        // r_proton[i] = range(e[i], f_proton);
        // cout << "proton E: " << e[i] << " , range: " << r_proton[i] << endl;
      }
      double rr = r_muon[M-1];
      double tt = t_muon[M-1];
      double tt_ar = t_muon_ar[M-1];
      double ee = e[M-1];
      double contam[M];
      for (int i=0;i!=M;i++){
	r_muon[i] = rr - r_muon[i];
	t_muon[i] = tt - t_muon[i];
	t_muon_ar[i] = tt_ar - t_muon_ar[i];
	e[i] = (ee - e[i])/(3.875*120e-4);
	contam[i] = (1-exp(-t_muon_ar[i]/2.2e3))/1.235e-8;
      }
      
      TGraph *g1 = new TGraph(M, r_muon, t_muon_ar);
      g1->Draw("AL");
      g1->GetXaxis()->SetTitle("Range (#mum)");
      g1->GetYaxis()->SetTitle("Time (ns)");
      g1->SetTitle("MuonDIF");
      g1->GetYaxis()->SetTitleOffset(1.3);
      std::cout << g1->Eval(200) << " " << g1->Eval(100) << std::endl;

      c1->cd(2);
      TGraph *g5 = new TGraph(M, r_muon, contam);
      g5->Draw("AL");
      g5->GetXaxis()->SetTitle("Range (#mum)");
      g5->GetYaxis()->SetTitle("(MuonDIF contam. w.r.t. 1% Pienu tail)");
      g5->SetTitle("MuonDIF");
      g5->GetYaxis()->SetRangeUser(0,50);
      g5->GetXaxis()->SetRangeUser(0,100);
      
      c1->cd(3);
      TGraph *g3 = new TGraph(M, r_muon, e);
      g3->GetXaxis()->SetTitle("Range (#mum)");
      g3->GetYaxis()->SetTitle("E^{#mu}_{depo}/E_{MIP-120-#mum}");
      g3->SetTitle("MuonDIF");
      g3->Draw("AL");

      c1->cd(4);
      TGraph *g4 = new TGraph(M, e, contam);
      g4->Draw("AL");
      g4->SetTitle("(MuonDIF contam. w.r.t. 1% Pienu tail) vs. E^{#mu}_{depo}/E_{MIP-120-#mum}");
      g4->GetXaxis()->SetRangeUser(0,3);
      g4->GetYaxis()->SetRangeUser(0,20);
    }

    /*
    c1->cd(2);
    
    {
       // const int M = 10000+1;
      const int M = 1000+1;
      
      double e[M];
      double r_proton[M], r_kaon[M], r_pion[M], r_muon[M];
      e[0] = r_proton[0] = r_kaon[0] = r_pion[0] = r_muon[0] = 0;
      double dE = 10.44/M; // MeV
      double t_pion[M], t_pion_ar[M];
      t_pion[0] = 0.; t_pion_ar[0] = 0;
      for (int i=1; i<M; i++) {
        e[i] = i * dE;
        r_proton[i] = r_proton[i-1] + f_proton->Eval(e[i]) * dE * 1e4; // um
    	r_kaon[i] =   r_kaon[i-1] +   f_kaon->Eval(e[i]) * dE * 1e4;
    	r_pion[i] =   r_pion[i-1] +   f_pion->Eval(e[i]) * dE * 1e4;
    	r_muon[i] =   r_muon[i-1] +   f_muon->Eval(e[i]) * dE * 1e4;
	
	// muon velocity
	double pion_vel = sqrt(1-pow(pion->Mass()*1000.,2)/pow(e[i] +pion->Mass()*1000. ,2))*3e5; // um/ns ...
	t_pion[i] = t_pion[i-1] + f_pion->Eval(e[i]) * dE * 1e4/pion_vel;
	t_pion_ar[i] = t_pion_ar[i-1] + f_pion->Eval(e[i]) * dE * 1e4/pion_vel * pow(pion->Mass()*1000.,2)/pow(e[i] +pion->Mass()*1000. ,2);
	
        // r_proton[i] = range(e[i], f_proton);
        // cout << "proton E: " << e[i] << " , range: " << r_proton[i] << endl;
      }
      double rr = r_pion[M-1];
      double tt = t_pion[M-1];
      double tt_ar = t_pion_ar[M-1];
      for (int i=0;i!=M;i++){
	r_pion[i] = rr - r_pion[i];
	t_pion[i] = tt - t_pion[i];
	t_pion_ar[i] = tt_ar - t_pion_ar[i];
      }
      
      TGraph *g2 = new TGraph(M, r_pion, t_pion_ar);
      g2->Draw("AL");
      g2->GetXaxis()->SetTitle("Range (#mum)");
      g2->GetYaxis()->SetTitle("Time (ns)");
      g2->SetTitle("PionDIF");
      g2->GetYaxis()->SetTitleOffset(1.3);
      std::cout << g2->Eval(200) << " " << g2->Eval(100) << std::endl;

      
      
    }

    */    

    // TH2F *h2 = new TH2F("h2","", 100, 0, r_muon[M-1]+100, 100, 0, e[M-1]+dE);
    // h2->GetYaxis()->SetTitle("Kinetic Energy [MeV]");
    // h2->GetXaxis()->SetTitle("Range [#mum]");
    // h2->Draw();

    // TGraph *g_proton = new TGraph(M, r_proton, e);
    // g_proton->SetLineWidth(2);
    // g_proton->SetLineColor(kRed);
    // g_proton->Draw("L,same");

    // TGraph *g_kaon = new TGraph(M, r_kaon, e);
    // g_kaon->SetLineWidth(2);
    // g_kaon->SetLineColor(kCyan);
    // g_kaon->Draw("L,same");

    // TGraph *g_pion = new TGraph(M, r_pion, e);
    // g_pion->SetLineWidth(2);
    // g_pion->SetLineColor(kViolet);
    // g_pion->Draw("L,same");

    // TGraph *g_muon = new TGraph(M, r_muon, e);
    // g_muon->SetLineWidth(2);
    // g_muon->SetLineColor(kGreen);
    // g_muon->Draw("L,same");

    // TLegend *leg = new TLegend(0.60,0.20, 0.87,0.47);
    // leg->AddEntry(g_proton, " proton", "l");
    // leg->AddEntry(g_kaon,   " kaon", "l");
    // leg->AddEntry(g_pion,   " pion", "l");
    // leg->AddEntry(g_muon,   " muon", "l");
    // leg->SetFillColor(kWhite);
    // leg->Draw();

    // gPad->SetGridx();
    // gPad->SetGridy();

    // // dedx as function of residual range 
    // TCanvas *c2 = new TCanvas;
    // const int M2 = 1000;
    // double rx[M2];
    // double dedx_proton[M2], dedx_kaon[M2], dedx_pion[M2], dedx_muon[M2], dedx_mip[M2];
    // // rx[0] = dedx_proton[0] = dedx_kaon[0] = dedx_pion[0] = dedx_muon[0] = 0;
    // double dx = 1.; // um
    // double mip = 3.875; // MIP energy
    // for (int i=0; i<M2; i++) {
    //     rx[i] = (i+1) * dx;
    //     dedx_proton[i] =  1. / f_proton->Eval( g_proton->Eval(rx[i]) ) / mip; // MeV/cm
    //       dedx_kaon[i] =  1. /   f_kaon->Eval(   g_kaon->Eval(rx[i]) ) / mip;
    //       dedx_pion[i] =  1. /   f_pion->Eval(   g_pion->Eval(rx[i]) ) / mip;
    //       dedx_muon[i] =  1. /   f_muon->Eval(   g_muon->Eval(rx[i]) ) / mip;
    //       dedx_mip[i] = mip / mip;
    //     //   cout << "rx: " << rx[i] << " , proton dedx: " << dedx_proton[i] << endl;

    // }

    // TH2F *h22 = new TH2F("h22","", 100, 0, 800, 100, 0, 50);
    // h22->GetXaxis()->SetTitle("Residual Range [#mum]");
    // h22->GetYaxis()->SetTitle("dE/dx / (3.875 MeV/cm) ");
    // h22->Draw();

    // // TGraph *g2_proton = new TGraph(M2, rx, dedx_proton);
    // // g2_proton->SetLineWidth(2);
    // // g2_proton->SetLineColor(kRed);
    // // g2_proton->Draw("L,same");

    // // TGraph *g2_kaon = new TGraph(M2, rx, dedx_kaon);
    // // g2_kaon->SetLineWidth(2);
    // // g2_kaon->SetLineColor(kCyan);
    // // g2_kaon->Draw("L,same");

    // for (Int_t i=0;i!=M2;i++){
    //   dedx_pion[i] /= dedx_mip[i];
    //   dedx_muon[i] /= dedx_mip[i];
    //   dedx_mip[i] = 1.;
    // }
    
    // TGraph *g2_pion = new TGraph(M2, rx, dedx_pion);
    // g2_pion->SetLineWidth(2);
    // g2_pion->SetLineColor(kRed);
    // g2_pion->Draw("L,same");

    // TGraph *g2_muon = new TGraph(M2, rx, dedx_muon);
    // g2_muon->SetLineWidth(2);
    // g2_muon->SetLineColor(kBlack);
    // g2_muon->Draw("L,same");

    // TGraph *g2_mip = new TGraph(M2, rx, dedx_mip);
    // g2_mip->SetLineWidth(2);
    // g2_mip->SetLineColor(kMagenta);
    // g2_mip->Draw("L,same");

    // TLegend *leg2 = new TLegend(0.68,0.70, 0.87,0.87);
    // // leg2->AddEntry(g2_proton, " proton", "l");
    // // leg2->AddEntry(g2_kaon,   " kaon", "l");
    // leg2->AddEntry(g2_pion,   " #pi", "l");
    // leg2->AddEntry(g2_muon,   " #mu ~ 0.9#pi", "l");
    // leg2->AddEntry(g2_mip, " e ~ MIP", "l");
    // leg2->SetFillColor(kWhite);
    

    // TLine *l1 = new TLine(120,0,120,50);
    // l1->Draw();
    // TLine *l2 = new TLine(240,0,240,50);
    // l2->Draw();
    // TLine *l3 = new TLine(360,0,360,50);
    // l3->Draw();
    // TLine *l4 = new TLine(480,0,480,50);
    // l4->Draw();
    // TLine *l5 = new TLine(600,0,600,50);
    // l5->Draw();
    // TLine *l6 = new TLine(720,0,720,50);
    // l6->Draw();

    // leg2->Draw();
    // //gPad->SetGridx();
    // //gPad->SetGridy();
}


