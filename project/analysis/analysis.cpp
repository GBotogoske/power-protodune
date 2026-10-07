
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <TLegend.h>

void analysis(std::string filename="../build/Data.root")
{
    // Open file and get tree

    TFile *file = TFile::Open(filename.c_str());
    TTree *tree = (TTree*)file->Get("Primary_Hit");

    // Create 2D histogram for internal region
    TH2D *hVisInt = new TH2D("hvis","VIS(z,y) - Internal", 36, -4.64103, 4.64103, 32, -4, 4);
    TH2D *hUvInt = new TH2D("huv","VUV(z,y) - Internal",  36, -4.64103, 4.64103, 32, -4, 4);

     // Create 2D histogram for external region
    TH2D *hVisExt = new TH2D("hvis","VIS(z,y)  - External", 36, -4.272, 4.272, 32, -4, 4);
    TH2D *hUvExt = new TH2D("huv","VUV(z,y) - External", 36, -4.272, 4.272, 32, -4, 4);

    int internal_bins_z=12;
    int internal_bins_y=14;

    // Fill the histogram
    Double_t z, y, f, f2;
    tree->SetBranchAddress("Z", &z);
    tree->SetBranchAddress("Y", &y);
    tree->SetBranchAddress("PhotonDetectedVIS", &f);
    tree->SetBranchAddress("PhotonDetectedUV", &f2);

    int number_files=1;
    // float gain=40000.0/25000.0;
    // float ratio=1.0;
    float scale=0.1;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) 
    {
        tree->GetEntry(i);

        if((z>(-1.547) && z<1.547) && (y>(-3.5) && y<3.5))
        {            
            if(f>=1)
            {
                hVisInt->Fill(z, y,f);
            }
                
            if(f2>=1) 
            {
                hUvInt->Fill(z,y,f2);
            }
        }
        else
        {
            // ratio = gain/number_files;
            if(f>=1)
            {
                hVisExt->Fill(z, y,f);
            }
                
            if(f2>=1) 
            {
                hUvExt->Fill(z,y,f2);
            }
        }
    }

    // Draw heatmap
    TCanvas *c1 = new TCanvas();
    hVisInt->Draw("COLZ");  // COLZ = color map
    c1->SaveAs("heatmapVISInt.png");

    TCanvas *c2 = new TCanvas();
    hUvInt->Draw("COLZ");  // COLZ = color map
    c2->SaveAs("heatmapUVInt.png");

    TCanvas *c3 = new TCanvas();
    hVisExt->Draw("COLZ");  // COLZ = color map
    c3->SaveAs("heatmapVISExt.png");

    TCanvas *c4 = new TCanvas();
    hUvExt->Draw("COLZ");  // COLZ = color map
    c4->SaveAs("heatmapUVExt.png");

   // --- Converter h2 para matriz ---
    int nx = hVisExt->GetNbinsX();
    int ny = hVisExt->GetNbinsY();

        // Criar histograma para o buffer
    TH2D *h4 = new TH2D("hall","VUV+VIS(z,y)", nx, -4.64103, 4.64103, ny, -4, 4);
    TH2D *hbuffer = new TH2D("hbuffer","UV_{frac}", nx, -4.64103, 4.64103, ny, -4, 4);
    TH2D *hratio = new TH2D("hratio","UV/VIS", nx, -4.64103, 4.64103, ny, -4, 4);


    std::vector<std::vector<double>> mat_uv(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> mat_vis(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> mat_all(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> buffer(nx, std::vector<double>(ny,0));
    std::vector<std::vector<double>> myratio(nx, std::vector<double>(ny,0));

    double lightyield_avg=0;
    double lightyield_min=10000;

    double cut_pe = 0.0;

    for(int i=1; i<=nx; i++)        
    {
        for(int j=1; j<=ny; j++)
        {
            if ((i > internal_bins_z && i <= nx-internal_bins_z) && (j > 2 && j <= ny-2)){
                mat_vis[i-1][j-1] = hVisInt->GetBinContent(i,j);
                mat_uv[i-1][j-1] = hUvInt->GetBinContent(i,j);
            } else {
                mat_vis[i-1][j-1] = hVisExt->GetBinContent(i,j);
                mat_uv[i-1][j-1] = hUvExt->GetBinContent(i,j);
            }
            
            if(mat_vis[i-1][j-1]<cut_pe)
            {
                mat_vis[i-1][j-1]=0;
            }
            
            if(mat_uv[i-1][j-1]<cut_pe)
            {
                mat_uv[i-1][j-1]=0;
            }

            mat_all[i-1][j-1] = mat_uv[i-1][j-1] +  mat_vis[i-1][j-1];
            buffer[i-1][j-1] = mat_uv[i-1][j-1]/(mat_uv[i-1][j-1] + mat_vis[i-1][j-1] );
            myratio[i-1][j-1] = mat_uv[i-1][j-1]/(mat_vis[i-1][j-1] );
        }
    }

    // Preencher com os valores da matriz buffer
    for (int i = 1; i <= nx; i++)
    {
        for (int j = 1; j <= ny; j++) 
        {
            hbuffer->SetBinContent(i, j, buffer[i-1][j-1]);
            hratio->SetBinContent(i, j, myratio[i-1][j-1]);
            h4->SetBinContent(i, j, mat_all[i-1][j-1]);
        }
    }

    int cont=0;
    for(int i=internal_bins_y-1; i<nx-internal_bins_y; i++)        
    {
        for(int j=2; j<(ny/2-1); j++)
        {
            //std::cout << mat_all[i][j] << std::endl;
            lightyield_avg+=mat_all[i][j];
            cont++;
            if(mat_all[i][j]<lightyield_min)
            {
                lightyield_min = mat_all[i][j];
            }
        }
    }
    lightyield_avg=lightyield_avg/cont;
    std::cout << "AVG Light Yield: " << lightyield_avg << std::endl;
    std::cout << "MIN Light Yield: " << lightyield_min << std::endl;

    // Plotar
    TCanvas *c5 = new TCanvas();
    hbuffer->GetZaxis()->SetRangeUser(0.0,1.0); // já que é uma fração
    hbuffer->Draw("COLZ");
    hbuffer->SetStats(0);
    c5->SaveAs("heatmapBuffer.png");

    TCanvas *c6 = new TCanvas();
    hratio->GetZaxis()->SetRangeUser(0.0,1.0); // já que é uma fração
    hratio->Draw("COLZ");
    c6->SaveAs("heatmapRatio.png");

    TCanvas *c7 = new TCanvas();
    h4->Draw("COLZ");  // COLZ = color map
    c7->SaveAs("heatmapAll.png");

    TH1D *h1d = new TH1D("hist_franc","Histogram Francesco", 100, 0 , 210 );
    for(int i=0; i<nx; i++)        
    {
        for(int j=0; j<ny; j++)
        {
            h1d->Fill(mat_all[i][j]);
        }
    }
    TCanvas *c8 = new TCanvas();
    h1d->Draw();  // COLZ = color map
    c8->SaveAs("histFrancesco.png");

   TH1D *h1dR = new TH1D("hist_R","UV_{frac} distribution (normalized);UV_{frac};Entries", 100, -0.1 , 1.1 );
    for(int i=0; i<nx; i++)        
    {
        for(int j=0; j<ny; j++)
        {
            h1dR->Fill(buffer[i][j]);
        }
    }

    TH1D *h1dRin = new TH1D("hist_R_in","R (inner);R value;Entries", 100, -0.1 , 1.1 );
    TH1D *h1dRout = new TH1D("hist_R_out","R (outer);R value;Entries", 100, -0.1 , 1.1 );

    double out_cont =0.0;
    double out_cont_cut=0.0;
    double out_cut=0.2;

    for(int i=0; i<nx; i++)        
    {
        for(int j=0; j<(ny-1); j++)
        {
            if(i>(internal_bins_y-1) && i<=(nx-internal_bins_y) && (j > 1) && (j < (ny/2-1))){
                h1dRin->Fill(buffer[i][j]);
            } else {
                h1dRout->Fill(buffer[i][j]);
                if(buffer[i][j]>=out_cut)
                {
                    out_cont_cut=out_cont_cut+1;
                }
                out_cont=out_cont+1;
            }
        }
    }


    TCanvas *c9 = new TCanvas("c7","UV_{frac} distributions",800,600);
    c9->SetGrid();

    h1dR->SetLineColor(kBlack);
    h1dR->SetLineWidth(2);
    h1dR->SetStats(0);
    h1dRin->SetLineColor(kBlue+1);
    h1dRin->SetLineWidth(2);
    h1dRin->SetStats(0);
    h1dRout->SetLineColor(kRed+1);
    h1dRout->SetLineWidth(2);

    h1dR->Draw("HIST");
    h1dR->Scale(1/h1dR->GetMaximum());
    h1dRin->Scale(1/h1dRin->GetMaximum());
    h1dRin->Draw("HIST SAME");
    h1dRout->Scale(1/h1dRout->GetMaximum());
    h1dRout->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.55,0.7,0.88,0.88);
    leg->AddEntry(h1dR,"All entries","l");
    leg->AddEntry(h1dRin,"Inner region","l");
    leg->AddEntry(h1dRout,"Outer region","l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    c9->SaveAs("histRATIO.png");

    std::cout << "efficiency: " << out_cont_cut/out_cont <<std::endl;

}

