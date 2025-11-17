
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

    // Create 2D histogram
    TH2D *h2 = new TH2D("hvis","VIS(z,y)", 26, -6.5, 6.5, 26, -6.5, 6.5);
    TH2D *h3 = new TH2D("huv","VUV(z,y)", 26, -6.5, 6.5, 26, -6.5, 6.5);

    // Fill the histogram
    Double_t z, y, f, f2;
    tree->SetBranchAddress("Z", &z);
    tree->SetBranchAddress("Y", &y);
    tree->SetBranchAddress("PhotonDetectedVIS", &f);
    tree->SetBranchAddress("PhotonDetectedUV", &f2);

    int number_files=1;
    float gain=40000.0/25000.0;
    float ratio=5.0;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) 
    {
        tree->GetEntry(i);

        if(z<(-6) || z>6)
        {
            ratio=gain/number_files;
        }
        else
        {
            ratio=1.0/number_files;
        }

        if(f>=1)
        {
            h2->Fill(z, y,f*ratio);
        }
            
        if(f2>=1) 
        {
            h3->Fill(z,y,f2*ratio);
        } // fill with weight = f
    }

    // Draw heatmap
    TCanvas *c1 = new TCanvas();
    h2->Draw("COLZ");  // COLZ = color map
    c1->SaveAs("heatmapVIS.png");

    TCanvas *c2 = new TCanvas();
    h3->Draw("COLZ");  // COLZ = color map
    c2->SaveAs("heatmapUV.png");

   // --- Converter h2 para matriz ---
    int nx = h2->GetNbinsX();
    int ny = h2->GetNbinsY();

        // Criar histograma para o buffer
    TH2D *h4 = new TH2D("hall","VUV+VIS(z,y)", nx, -6.5, 6.5, ny, -6.5, 6.5);
    TH2D *hbuffer = new TH2D("hbuffer","UV/(UV+VIS)", nx, -6.5, 6.5, ny, -6.5, 6.5);
    TH2D *hratio = new TH2D("hratio","UV/VIS", nx, -6.5, 6.5, ny, -6.5, 6.5);


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
            mat_vis[i-1][j-1] = h2->GetBinContent(i,j);
            mat_uv[i-1][j-1] = h3->GetBinContent(i,j);
            if(mat_vis[i-1][j-1]<cut_pe)
            {
                mat_vis[i-1][j-1]=0;
            }
             if(mat_uv[i-1][j-1]<cut_pe)
            {
                mat_uv[i-1][j-1]=0;
            }
            mat_all[i-1][j-1] = mat_uv[i-1][j-1] +  mat_vis[i-1][j-1];
            buffer[i-1][j-1] = mat_uv[i-1][j-1]/(  mat_uv[i-1][j-1] + mat_vis[i-1][j-1] );
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
    for(int i=1; i<nx-1; i++)        
    {
        for(int j=0; j<ny; j++)
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
    TCanvas *c3 = new TCanvas();
    hbuffer->GetZaxis()->SetRangeUser(0.0,1.0); // já que é uma fração
    hbuffer->Draw("COLZ");
    c3->SaveAs("heatmapBuffer.png");

    TCanvas *c4 = new TCanvas();
    hratio->GetZaxis()->SetRangeUser(0.0,1.0); // já que é uma fração
    hratio->Draw("COLZ");
    c4->SaveAs("heatmapRatio.png");

    TCanvas *c5 = new TCanvas();
    h4->Draw("COLZ");  // COLZ = color map
    c5->SaveAs("heatmapAll.png");

    TH1D *h1d = new TH1D("hist_franc","Histogram Francesco", 100, 0 , 210 );
    for(int i=0; i<nx; i++)        
    {
        for(int j=0; j<ny; j++)
        {
            h1d->Fill(mat_all[i][j]);
        }
    }
    TCanvas *c6 = new TCanvas();
    h1d->Draw();  // COLZ = color map
    c6->SaveAs("histFrancesco.png");

   TH1D *h1dR = new TH1D("hist_R","R;R value;Entries", 200, -0.1 , 1.1 );
    for(int i=0; i<nx; i++)        
    {
        for(int j=0; j<ny; j++)
        {
            h1dR->Fill(buffer[i][j]);
        }
    }

    TH1D *h1dRin = new TH1D("hist_R_in","R (inner);R value;Entries", 200, -0.1 , 1.1 );
    TH1D *h1dRout = new TH1D("hist_R_out","R (outer);R value;Entries", 200, -0.1 , 1.1 );

    double out_cont =0.0;
    double out_cont_cut=0.0;
    double out_cut=0.2;

    for(int i=0; i<nx; i++)        
    {
        for(int j=0; j<ny; j++)
        {
            if(i!=0 && i!=(nx-1))
                h1dRin->Fill(buffer[i][j]);
            else
            {
                h1dRout->Fill(buffer[i][j]);
                if(buffer[i][j]>=out_cut)
                {
                    out_cont_cut=out_cont_cut+1;
                }
                out_cont=out_cont+1;
            }       
        }
    }


    TCanvas *c7 = new TCanvas("c7","R distributions",800,600);
    c7->SetGrid();

    h1dR->SetLineColor(kBlack);
    h1dR->SetLineWidth(2);
    h1dRin->SetLineColor(kBlue+1);
    h1dRin->SetLineWidth(2);
    h1dRout->SetLineColor(kRed+1);
    h1dRout->SetLineWidth(2);

    h1dR->Draw("HIST");
    h1dRin->Draw("HIST SAME");
    h1dRout->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.65,0.7,0.88,0.88);
    leg->AddEntry(h1dR,"All entries","l");
    leg->AddEntry(h1dRin,"Inner region","l");
    leg->AddEntry(h1dRout,"Outer region","l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    c7->SaveAs("histRATIO.png");

    std::cout << "efficiency: " << out_cont_cut/out_cont <<std::endl;

}

