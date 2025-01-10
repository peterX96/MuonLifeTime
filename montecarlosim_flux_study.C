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
#include <vector>  
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
Double_t a[25];     
int type[5];          
int type1[5]; 
int sum,cross;

int b, c;
Double_t x, y, phi, theta;  // Variables for the line that a particle travels
Double_t x_int1, y_int1;
Double_t x_int2, y_int2;

// boh

Double_t y_xl,y_xr,x_yf,x_yd;
double x3,y3,z3,x4,y4,z4;

// Compute distance that particle travels inside the slab
Double_t distance3D(Double_t x_a, Double_t y_a, Double_t z_a, Double_t x_b, Double_t y_b, Double_t z_b) {
    TVector3 point1(x_a, y_a, z_a);
    TVector3 point2(x_b, y_b, z_b);
    return (point1 - point2).Mag();
}

// Line starting from the point (x, y) with angle (phi, theta)
// Given a certain height, this function gives the intersection with a plane parallel to the slabs
void line(double x, double y, double phi, double theta, double height, double& x_intersection, double& y_intersection) {
    double t = height / TMath::Cos(theta);
    x_intersection = x + t * TMath::Sin(theta) * TMath::Cos(phi);
    y_intersection = y + t * TMath::Sin(theta) * TMath::Sin(phi);
}

// Given a certain line, this function gives the intersection with the planes that define the box of slabs
// we check the result to understand if there is an intersection with the lateral surface of each slab
void boundary(double x0, double y0, double phi, double theta, int n_slab,int& sum, double& x1,double& y1,double& z1,double& x2,double& y2,double& z2) {
    double t_xl = ((-slab_length / 2)-x0)/(TMath::Sin(theta) * TMath::Cos(phi));
    double t_xr = ((slab_length / 2)-x0)/(TMath::Sin(theta) * TMath::Cos(phi));
    double t_yf = ((slab_width / 2)-y0)/(TMath::Sin(theta) * TMath::Sin(phi));
    double t_yd = ((-slab_width / 2)-y0)/(TMath::Sin(theta) * TMath::Sin(phi));
    y_xl = y0 + t_xl * TMath::Sin(theta) * TMath::Sin(phi);
    y_xr = y0 + t_xr * TMath::Sin(theta) * TMath::Sin(phi);
    x_yf = x0 + t_yf * TMath::Sin(theta) * TMath::Cos(phi);
    x_yd = x0 + t_yd * TMath::Sin(theta) * TMath::Cos(phi);

    if (t_xl > 0){   
        if (std::abs(y_xl) > slab_width / 2){
            t_xl = 0;}
    }
    else {
        t_xl = 0;
    }
    if (t_xr > 0)
    {   
        if (std::abs(y_xr) > slab_width / 2){
            t_xr = 0;}
    }
    else {
        t_xr = 0;
    }

    if (t_yf > 0)
    {  
        if (std::abs(x_yf) > slab_length / 2){
           t_yf  = 0;}
    }
    else {
        t_yf = 0;
    }

    if (t_yd > 0)
    {
        if (std::abs(x_yd) > slab_length / 2){
          t_yd = 0;}
    }
    else {
       t_yd = 0;
    }
    // z is 0 if there is no intersection
    double z_xl = t_xl*TMath::Cos(theta);
    double z_xr = t_xr*TMath::Cos(theta);
    double z_yf = t_yf*TMath::Cos(theta);
    double z_yd = t_yd*TMath::Cos(theta);

    int left = ( z_xl >sp[n_slab]  &&  z_xl <sm[n_slab]) ? 1 : 0;
    int right = ( z_xr >sp[n_slab]  &&  z_xr <sm[n_slab]) ? 1 : 0;
    int forward = ( z_yf >sp[n_slab]  &&  z_yf <sm[n_slab]) ? 1 : 0;
    int downward = ( z_yd >sp[n_slab]  &&  z_yd <sm[n_slab]) ? 1 : 0;
    sum = forward + downward + left + right;

    Double_t v_l[] = {-slab_length / 2,y_xl,z_xl};
    Double_t v_r[] = {slab_length / 2,y_xr,z_xr};
    Double_t v_f[] = {x_yf,slab_width / 2,z_yf};
    Double_t v_d[] = {x_yd,-slab_width / 2,z_yd};

    if (sum == 0)
    {
        x1 = 0;
        y1 = 0;
        z1 = 0;
        x2 = 0;
        y2 = 0;
        z2 = 0;

    }
    else if(sum == 1){

        if (v_l[2] == 0){}

        else  {
        x1 = v_l[0];
        y1 = v_l[1];
        z1 = v_l[2];
        }

        if (v_r[2] == 0){}

        else  {
        x1 = v_r[0];
        y1 = v_r[1];
        z1 = v_r[2];
        }

        if (v_d[2] == 0){}

        else  {
        x1 = v_d[0];
        y1 = v_d[1];
        z1 = v_d[2];
        }

        if (v_f[2] == 0){}

        else  {
        x1 = v_f[0];
        y1 = v_f[1];
        z1 = v_f[2];
        }
        x2=0;
        y2=0;
        z2=0;
        }

    else{ 
                    // find first point
                    if (v_l[2] == 0 && v_r[2] == 0) {      
                // Both v_l[2] and v_r[2] are zero, use v_f and v_d for x1, y1, z1, x2, y2, z2
                x1 = v_f[0];
                y1 = v_f[1];
                z1 = v_f[2];
                x2 = v_d[0];
                y2 = v_d[1];
                z2 = v_d[2];
            }
            else if (v_l[2] == 0 && v_f[2] == 0) {
                // Both v_l[2] and v_f[2] are zero, use v_r and v_d for x1, y1, z1, x2, y2, z2
                x1 = v_r[0];
                y1 = v_r[1];
                z1 = v_r[2];
                x2 = v_d[0];
                y2 = v_d[1];
                z2 = v_d[2];
            }
            else if (v_l[2] == 0 && v_d[2] == 0) {
                // Both v_l[2] and v_d[2] are zero, use v_r and v_f for x1, y1, z1, x2, y2, z2
                x1 = v_r[0];
                y1 = v_r[1];
                z1 = v_r[2];
                x2 = v_f[0];
                y2 = v_f[1];
                z2 = v_f[2];
            }
            else if (v_r[2] == 0 && v_f[2] == 0) {
                // Both v_r[2] and v_f[2] are zero, use v_l and v_d for x1, y1, z1, x2, y2, z2
                x1 = v_l[0];
                y1 = v_l[1];
                z1 = v_l[2];
                x2 = v_d[0];
                y2 = v_d[1];
                z2 = v_d[2];
            }
            else if (v_r[2] == 0 && v_d[2] == 0) {
                // Both v_r[2] and v_d[2] are zero, use v_l and v_f for x1, y1, z1, x2, y2, z2
                x1 = v_l[0];
                y1 = v_l[1];
                z1 = v_l[2];
                x2 = v_f[0];
                y2 = v_f[1];
                z2 = v_f[2];
            }
            else if (v_f[2] == 0 && v_d[2] == 0) {
                // Both v_f[2] and v_d[2] are zero, use v_l and v_r for x1, y1, z1, x2, y2, z2
                x1 = v_l[0];
                y1 = v_l[1];
                z1 = v_l[2];
                x2 = v_r[0];
                y2 = v_r[1];
                z2 = v_r[2];
            }
            else {
            }
        }
    //std::cout <<"SLAB NR:"<< n_slab <<"\t"<< up <<"\t"<< down <<"\t"<<  left <<"\t"<< right <<"\t"<< std::endl;
    //std::cout << t_xl <<"\t"<< t_xr <<"\t"<<  t_yu <<"\t"<< t_yd <<"\t"<< std::endl;
    //std::cout << y_xl <<"\t"<< y_xr <<"\t"<<  x_yu <<"\t"<< x_yd <<"\t"<< std::endl;   
    //std::cout << z_xl <<"\t"<< z_xr <<"\t"<<  z_yu <<"\t"<< z_yd <<"\t"<< std::endl;
}

// Check if particle crosses a surface of the slab
int cross0(double x, double y) {
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

    tree->Branch("cross", &cross);

    for (Int_t iev = 0; iev < Mtree->GetEntries(); iev++) {
        Mtree->GetEntry(iev);   // loop over all particle

        n_cross_slab = 0;

        for (Int_t j = 0; j < 5; j++) {

            //boundary(double x0, double y0, double phi, double theta, int n_slab, int& x1,int& y1,int& z1,int& x2,int& y2,int& z2)

            boundary(x, y, phi, theta,j,sum,x3,y3,z3,x4,y4,z4); // (sum = total number of lateral side) of the j-slab that the particle line of the iev-th particle cross

            /*
            if (sum == 2 && j == 0)
            {
                cout <<"SLAB NR:"<< j << "\t" << x3 << "\t" << y3 << "\t" << z3 << "\t" << x4 << "\t" << y4 << "\t" << z4 << endl;
            }
            */

            

            // loop over all slabs

            line(x, y, phi, theta, sp[j], x_int1, y_int1);  // point in which iev-th particle cross the up surface of j-slab

            a[0 + j * 5] = x_int1;
            a[1 + j * 5] = y_int1;

            line(x, y, phi, theta, sm[j], x_int2, y_int2);  // point in which iev-th particle cross the down surface of j-slab

            a[2 + j * 5] = x_int2;
            a[3 + j * 5] = y_int2;

            b = cross0(x_int1, y_int1);  // This particle cross upper surface of j-slab ? 1-True, 0-False
            c = cross0(x_int2, y_int2);  // This particle cross down surface of j-slab ? 1-True, 0-False

            if (sum == 0){
            a[4 + j * 5] = distance3D(x_int2 ,y_int2,sm[j], x_int1, y_int1, sp[j]);  // distance that particle travel inside j-slab
            }
            else if(sum == 1 && b == 1){

            a[4 + j * 5] = distance3D(x_int1,y_int1,sp[j],x3,y3,z3);  // distance that particle travel inside j-slab
            }
            else if(sum == 1 && c == 1){
            a[4 + j * 5] = distance3D(x_int2,y_int2,sm[j],x3,y3,z3);  // distance that particle travel inside j-slab
            }
            else{
                 a[4 + j * 5] = distance3D(x4,y4,z4,x3,y3,z3);  // distance that particle travel inside j-slab
            }

            // different category of particle
            /*
            if ( b == 1 && c == 1){
                type[j]=0;  // This particle cross both up and down surface
            }
            else if( b == 1){
                type[j]=1;  // This particle  only up surface

            }
            else if( c == 1){
                type[j]=2; // This particle cross only down surface
                
            }
            else{
                 type[j]=3; // This particle don't cross down or up surface
            }
            */

            // b & c are the case in which the particle travel crossing up or down surface
            // sum = 2 is the case in which the particle enter and exit lateraly - pretty rare
           if ( b == 1 || c == 1 || sum == 2){
                type1[j]=0;  // This particle cross the slab
                }
                else{
                type1[j]=1;  // This particle don't cross the slab
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
