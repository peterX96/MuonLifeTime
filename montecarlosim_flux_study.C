#include "TMath.h"
#include <iostream>
#include <cmath>
#include "TTree.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include <TVector3.h>
#include "generator.h"
#include "constants.h"
#include <map>
#include <string>
#include <vector>  // Include vector header
#include "TH1D.h"

// Position along z of surface up and down of slabs, z->-z 
Double_t s1p = 0, s1m = s1p + slab_height;
Double_t s2p = s1m + 2, s2m = s2p + slab_height;
Double_t s3p = s2m + 2, s3m = s3p + slab_height;
Double_t s4p = s3m + 2, s4m = s4p + slab_height;
Double_t s5p = s4m + 2, s5m = s5p + slab_height;
Double_t sp[] = {s1p,s2p,s3p,s4p,s5p};
Double_t sm[] = {s1m,s2m,s3m,s4m,s5m};

// Information variables
Double_t a[25];     // Array to store the information
int type[5];          // Array to store the type of intersection for the particle
int type1[5]; 

int b, c;
Double_t x, y, phi, theta;  // Variables for the line that a particle travels
Double_t x_int1, y_int1;
Double_t x_int2, y_int2;

// Compute distance that particle travels inside the slab
Double_t distance3D(Double_t x1, Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2) {
    TVector3 point1(x1, y1, z1);
    TVector3 point2(x2, y2, z2);
    return (point1 - point2).Mag();
}

// Line starting from the point (x, y) with angle (phi, theta)
// Given a certain height, this function gives the intersection with a plane parallel to the slabs
void line(double x, double y, double phi, double theta, double height, double& x_intersection, double& y_intersection) {
    double t = height / TMath::Cos(theta);
    x_intersection = x + t * TMath::Sin(theta) * TMath::Cos(phi);
    y_intersection = y + t * TMath::Sin(theta) * TMath::Sin(phi);
}

// Check if particle crosses a surface of the slab
int cross(double x, double y) {
    return (std::abs(x) < slab_length / 2 && std::abs(y) < slab_width / 2) ? 1 : 0;
}



void study(TString namefile,int p) {
    // Read the flux parameter generator
    TFile *inputFile = new TFile(namefile);

    // Open the tree
    TTree *Mtree = (TTree*)inputFile->Get("Tree_muons");
    Mtree->SetBranchAddress("x", &x);
    Mtree->SetBranchAddress("y", &y);
    Mtree->SetBranchAddress("phi", &phi);
    Mtree->SetBranchAddress("theta", &theta);

    // Create a file for muon intersection with detector
    TFile *file = new TFile(("Muon_intersection_%"+std::to_string(p)+".root").c_str(), "UPDATE");

    // Fill branch with intersection and path travelled by particle for each slab
    TTree *tree = new TTree("Intersection", "Intersection");

    // s[number = 1,..,5][p = up, m = down]_[x,y coord]
    // s[number = 1,..,5]_[d = distance , type and type1 is to follow]

    tree->Branch("s1p_x", &a[0]);  
    tree->Branch("s1p_y", &a[1]);
    tree->Branch("s1m_x", &a[2]);
    tree->Branch("s1m_y", &a[3]);
    tree->Branch("s1_d", &a[4]);
    tree->Branch("s1_type", &type[0]);
    tree->Branch("s1_type1", &type1[0]);

    tree->Branch("s2p_x", &a[5]);
    tree->Branch("s2p_y", &a[6]);
    tree->Branch("s2m_x", &a[7]);
    tree->Branch("s2m_y", &a[8]);
    tree->Branch("s2_d", &a[9]);
    tree->Branch("s2_type", &type[1]);
    tree->Branch("s2_type1", &type1[1]);

    tree->Branch("s3p_x", &a[10]);
    tree->Branch("s3p_y", &a[11]);
    tree->Branch("s3m_x", &a[12]);
    tree->Branch("s3m_y", &a[13]);
    tree->Branch("s3_d", &a[14]);
    tree->Branch("s3_type", &type[2]);
    tree->Branch("s3_type1", &type1[2]);

    tree->Branch("s4p_x", &a[15]);
    tree->Branch("s4p_y", &a[16]);
    tree->Branch("s4m_x", &a[17]);
    tree->Branch("s4m_y", &a[18]);
    tree->Branch("s4_d", &a[19]);
    tree->Branch("s4_type", &type[3]);
    tree->Branch("s4_type1", &type1[3]);

    tree->Branch("s5p_x", &a[20]);
    tree->Branch("s5p_y", &a[21]);
    tree->Branch("s5m_x", &a[22]);
    tree->Branch("s5m_y", &a[23]);
    tree->Branch("s5_d", &a[24]);
    tree->Branch("s5_type", &type[4]);
    tree->Branch("s5_type1", &type1[4]);

    for (Int_t iev = 0; iev < Mtree->GetEntries(); iev++) {
        Mtree->GetEntry(iev);   // loop over all particle

        for (Int_t j = 0; j < 5; j++) {
            
            // loop over all slabs

            line(x, y, phi, theta, sp[j], x_int1, y_int1);  // point in which iev-th particle cross the up surface of j-slab

            a[0 + j * 5] = x_int1;
            a[1 + j * 5] = y_int1;

            line(x, y, phi, theta, sm[j], x_int2, y_int2);  // point in which iev-th particle cross the down surface of j-slab

            a[2 + j * 5] = x_int2;
            a[3 + j * 5] = y_int2;

            a[4 + j * 5] = distance3D(x, y,sm[j], x_int1, y_int1, sp[j]);  // distance that particle travel inside j-slab

            b = cross(x_int1, y_int1);  // This particle cross upper surface of j-slab ? 1-True, 0-False
            c = cross(x_int2, y_int2);  // This particle cross down surface of j-slab ? 1-True, 0-False

            
            if ( b == 1 && c == 1){
                type[j]=0;  // This particle cross both up and down surface
            }
            else if( b == 1){
                type[j]=1;  // This particle cross only up surface

            }
            else if( c == 1){
                type[j]=2; // This particle cross only down surface
                
            }
            else{
                 type[j]=3; // This particle don't cross down or up surface
            }
            
           if ( b == 1 || c == 1){
                type1[j]=0;  // This particle cross up or down surface
                }
                else{
                type1[j]=1;
                }
        }
        tree->Fill();
        
    }
    
    // Write and close files
    tree->Write();
    file->Close();
    inputFile->Close();
}


void montecarlosim_flux_study(){

int p1[]={0,2,5,10,20,30,50,80,100,200,300,400,500,600,700,800,900,1000};

for(int i=0;i<20;i++){

int n = p1[i];

TString namefile = "Flux_100s_%"+std::to_string(n)+".root";

study(namefile,n);

}

}