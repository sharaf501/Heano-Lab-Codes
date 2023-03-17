//
//  Spatial Structure Fitting of Clinical Data (For Figure 6 in Abubakar et al.)
//  

//  Executes 2 file names. fp is list of parameter combination each with its mean of squared logarithmic residual while fp2 is list of recurrence times for clinical data and simulated data.

#include<iostream>
#include<fstream>
#include<random>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <algorithm>

// Constant Variables
#define N 2500 // Cell number or tissue size
#define SP 50 //  Spatial structure grid size
#define detect_size 5000
#define actual_size 1000000000 // tissue size for cancer diagnosis
#define r0 1.0 // fitness of type 0 cell
#define SRUN 25 // number of entries for clinical data relative to 100
//#define d 1.0
//#define d3 1.0
// death rate is undefined. To be calculated for each parameter combination

using namespace std;

int main(void){
    // new introduction for spatial structure
    int s[SP][SP],ns[SP][SP]; // a grid of SP x SP cells
    int i,j; // row number and column number (position) of a cell in a I x J grid
    
    //Modifiable conditions. To be deduced from clinical dataset
    double r1,r2,r3,u1,u2,u3,po,po1,po2,po3,d,d3;
   
    // Computer randomized parameters
    double op,fit,mu1,mu2,mu3,dd,rnb;
    int rni,rnj,rn4,rn3,rn2;
   
    //  Variables and Types
    int x0,x1,x2,x3;  // Number of Type 0, Type 1, Types S-1 and Type S cells
    int cycle,t1,t0,t2,a;
    double time;
    double g1,g2;
    double Et,Er,sVt,sVr,Vt,Vr,St,Sr,it[10000],ir[10000],tt[10000];
    int sumtt;
    double sumtr;
    int pro[3],bpro[3];
    double prop[3],bprop[3];
    double databox[100];
    int alpha,beta,gamma,delta;
    double totalrate; //for cell division
    double b1rate,b2rate,b3rate,b4rate;
    int trial;
    double Q,SQ,epsilon;
    int counter;
    
    // Input file name (ends with .dat). Specific for each clinical dataset
    FILE *fp;
    FILE *fp2;
    fp=fopen("4pmodel_sp_fit_OSCC.dat","w");
    fp2=fopen("pars_4pmodel_sp_fit_OSCC.dat","w");
    
    std::random_device rnd; //Random number function (rnd is variable name）
    std::mt19937_64 mt(rnd());// Random number generator
    std::uniform_real_distribution<> rand_dist(0,1); // possibility of a choice
    std::uniform_real_distribution<> rand_mu_dist(-4.5,-2); // range of values for mutation rates (in Log10 values)
    std::uniform_real_distribution<> rand_fit_dist(0.9,1.1); // range of values for cell fitness
    std::uniform_real_distribution<> rand_grow_dist(1.1,2.0); // range of values for growth rate
    std::uniform_real_distribution<> rand_tissue_dist(1.0,4.0); // range of values for tissue distribution
    
    std::uniform_int_distribution<> rand50(0, 49);  // Choice for cell death
    std::uniform_int_distribution<> rand4(0, 3); // 1st choice for cell division
    std::uniform_int_distribution<> rand3(0, 2); // 2nd choice for cell division
    std::uniform_int_distribution<> rand2(0, 1); // 2nd choice for cell division
    
    // number of trials to be run for parameter estimation
    for(trial=0;trial<200;trial++){
        
        po1=rand_mu_dist(mt); // random log10 value for u1
        po2=rand_mu_dist(mt); // random log10 value for u2
        po3=rand_mu_dist(mt); // random log10 value for u3
        r1=rand_fit_dist(mt); // random value for fitness of Type 1 cells
        r2=rand_fit_dist(mt); // random value for fitness of Type S-1 cells
        d=rand_tissue_dist(mt); // random value for death rate
        d3=d; // assume death rate of Type 0 and Type S are the same
        r3=d3*rand_grow_dist(mt);
        u1=pow(10.0,po1); // mutation rate from Type 0 to Type 1
        u2=pow(10.0,po2); // mutation rate from Type 1 to Type S-1
        u3=pow(10.0,po3); // mutation rate from Type S-1 to Type S
        
   //     po1=rand_mu_dist(mt);
   //     po2=rand_mu_dist(mt);
   //     po3=rand_mu_dist(mt);
        /*   r1=0.912;
         r2=1.035;
         d=2.989;
         d3=d;
         r3=5.946;
         
         u1=0.00021;
         u2=0.000476;
         u3=0.000047;
         */
        for(a=0;a<100;a++){
            databox[a]=0.0;
        }
        /*    databox[0]=1.74;
         databox[2]=15.0;
         databox[4]=36.7;
         databox[6]=70.0;
         */
        // input proportion of patients with cancer recurrence relative to 100 for each clinical dataset
        // values MUST be same with SRUN. For example, if SRUN is 25, put proportion for every 4th percentile - 4%, 8%, 16%, ........ 100%
        
        databox[0]=0.2;
        databox[1]=3.0;
        databox[2]=3.5;
        databox[3]=4.4;
        databox[4]=4.9;
        databox[5]=6.1;
        databox[6]=6.9;
        databox[7]=7.3;
        databox[8]=8.6;
        databox[9]=9.7;
        databox[10]=26.8;
        databox[16]=81.7;
        databox[18]=82.9;
        databox[20]=98.1;
        databox[21]=130.9;
        
        fprintf(fp2,"N=%d detect_size=%d actual_size=%d r0=%f r1=%f r2=%f r3=%f d=%f d3=%f u1=%f u2=%f u3=%f\n",N,detect_size,actual_size,r0,r1,r2,r3,d,d3,u1,u2,u3);
        
        //    for(po=-5;po<-4.9;po+=0.3){
        // for(r2=0.90;r2<1.11;r2+=0.02){
        // for(r3=1.05;r3<1.51;r3+=0.05){
        //   u1=(double)pow(10.0,po);
        
        sumtt=0;
        sumtr=0;
        sVt=0;
        sVr=0;
        
        for(i=0;i<3;i++){
            pro[i]=0;
            bpro[i]=0;
            prop[i]=0.0;
            bprop[i]=0.0;
            
        }
        
        for(i=0;i<SP;i++){
            for(j=0;j<SP;j++){
                if(s[i][j]==0){
                    x0++;
                }
                else if(s[i][j]==1){
                    x1++;
                }
                else if(s[i][j]==2)
                {
                    x2++;
                }
            }
        }
        for(a=0;a<10000;a++){
            tt[a]=0.0;
            it[a]=0.0;
            ir[a]=0.0;
        }
        
        double tr[SRUN];
        for(a=0;a<SRUN;a++){
            tr[a]=0.0;
        }
        for(cycle=0;cycle<SRUN;cycle++){
            
            for(i=0;i<SP;i++){
                for(j=0;j<SP;j++){
                    s[i][j]=0;
                    ns[i][j]=0;
                }
            }
            
            x3=0;
            time=0;
            
            for(t1=0;t1<1000000000;t1++){
                
                
                op=rand_dist(mt); //decides whilh trial to run
                fit=rand_dist(mt); //Moran proces, decides which cell divides
                mu1=rand_dist(mt);  //decides whether to mutate from Type 0 to Type 1
                mu2=rand_dist(mt); //decides whether to mutate from Type 0 to Type S-1
                mu3=rand_dist(mt); //decides whether to mutate from Type S-1 to Type S
                dd=rand_dist(mt);   //Moran process, decided which cell dies
                rni=rand50(mt); //randomly selects position i for cell death
                rnj=rand50(mt); //radomly selects position j for cell death
                rnb=rand_dist(mt); //random number that determines replacement cell
                rn4=rand4(mt); //randomly selects replacement cell option 1
                rn3=rand3(mt); //randomly selects replacement cell option 2
                rn2=rand2(mt); //randomly selects replacement cell option 3
                
                g2=d*(double)N+r3*(double)x3+d3*(double)x3; //total rate
                time=time+1/g2;
                
                
                if(op>=0 && op<d*N/g2){
                    if(rni==0){
                        if(rnj==0){
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else if(rnj==49){
                            totalrate=0.0;
                            if(s[rni][rnj-1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b1rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else{
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni][rnj-1]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b3rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                    }
                    else if(rni==49){
                        if(rnj==0){
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else if(rnj==49){
                            totalrate=0.0;
                            if(s[rni][rnj-1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else{
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni][rnj-1]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b3rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                    }
                    else{
                        if(rnj==0){
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b3rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else if(rnj==49){
                            totalrate=0.0;
                            if(s[rni][rnj-1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b3rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else{
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni][rnj-1]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b3rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b3rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b4rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b4rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b4rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate+b4rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=(b1rate+b2rate)/totalrate && rnb<(b1rate+b2rate+b3rate)/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                    }
                }
                
                else if(op>=d*N/g2 && op<(d*N/g2+r3*x3/g2)){
                    x3++;
                }
                else if(op>=(d*N/g2+r3*x3/g2) && op<=1){
                    x3--;
                }
                for(i=0;i<SP;i++){
                    for(j=0;j<SP;j++){
                        s[i][j]=ns[i][j];
                    }
                }

                if(x3>=detect_size){
                    x3=0;
                    break;
                }
                if(t1==999999999){
                    break;
                }
            }//End.
            
            x0=0;
            x1=0;
            x2=0;
            for(i=0;i<SP;i++){
                for(j=0;j<SP;j++){
                    if(s[i][j]==0){
                        x0++;
                    }
                    else if(s[i][j]==1){
                        x1++;
                    }
                    else if(s[i][j]==2)
                    {
                        x2++;
                    }
                }
            }

            //Moran process continues for remaining premalignant cells
            
            time=d*(double)N*(1/(r3-d3))*log((double)actual_size/(double)detect_size);
            
            for(t0=0; t0<time ; t0++){
                
                op=rand_dist(mt); //decides whilh trial to run
                fit=rand_dist(mt); //Moran proces, decides which cell divides
                mu1=rand_dist(mt);  //decides whether to mutate from Type 0 to Type 1
                mu2=rand_dist(mt); //decides whether to mutate from Type 0 to Type S-1
                mu3=rand_dist(mt); //decides whether to mutate from Type S-1 to Type S
                dd=rand_dist(mt);   //Moran process, decided which cell dies
                rni=rand50(mt); //randomly selects position i for cell death
                rnj=rand50(mt); //radomly selects position j for cell death
                rnb=rand_dist(mt); //random number that determines replacement cell
                rn4=rand4(mt); //randomly selects replacement cell option 1
                rn3=rand3(mt); //randomly selects replacement cell option 2
                rn2=rand2(mt); //randomly selects replacement cell option 3
                
                            
                if(rni==0){
                    if(rnj==0){
                        totalrate=0.0;
                        if(s[rni][rnj+1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj+1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj+1]==2){
                            b1rate=r2;
                        }
                        if(s[rni+1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni+1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni+1][rnj]==2){
                            b2rate=r2;
                        }
                        totalrate=b1rate+b2rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj+1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj+1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj+1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni+1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni+1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni+1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                    else if(rnj==49){
                        totalrate=0.0;
                        if(s[rni][rnj-1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj-1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj-1]==2){
                            b1rate=r2;
                        }
                        if(s[rni+1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni+1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni+1][rnj]==2){
                            b2rate=r2;
                        }
                        totalrate=b1rate+b2rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj-1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj-1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj-1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni+1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni+1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni+1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                    else{
                        totalrate=0.0;
                        if(s[rni][rnj+1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj+1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj+1]==2){
                            b1rate=r2;
                        }
                        if(s[rni+1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni+1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni+1][rnj]==2){
                            b2rate=r2;
                        }
                        if(s[rni][rnj-1]==0){
                            b3rate=1.0;
                        }
                        else if(s[rni][rnj-1]==1){
                            b3rate=r1;
                        }
                        else if(s[rni][rnj-1]==2){
                            b3rate=r2;
                        }
                        totalrate=b1rate+b2rate+b3rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj+1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj+1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj+1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                            if(s[rni+1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni+1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni+1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni][rnj-1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj-1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj-1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                }
                else if(rni==49){
                    if(rnj==0){
                        totalrate=0.0;
                        if(s[rni][rnj+1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj+1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj+1]==2){
                            b1rate=r2;
                        }
                        if(s[rni-1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni-1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni-1][rnj]==2){
                            b2rate=r2;
                        }
                        totalrate=b1rate+b2rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj+1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj+1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj+1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni-1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni-1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni-1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                    else if(rnj==49){
                        totalrate=0.0;
                        if(s[rni][rnj-1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj-1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj-1]==2){
                            b1rate=r2;
                        }
                        if(s[rni-1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni-1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni-1][rnj]==2){
                            b2rate=r2;
                        }
                        totalrate=b1rate+b2rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj-1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj-1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj-1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni-1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni-1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni-1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                    else{
                        totalrate=0.0;
                        if(s[rni][rnj+1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj+1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj+1]==2){
                            b1rate=r2;
                        }
                        if(s[rni-1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni-1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni-1][rnj]==2){
                            b2rate=r2;
                        }
                        if(s[rni][rnj-1]==0){
                            b3rate=1.0;
                        }
                        else if(s[rni][rnj-1]==1){
                            b3rate=r1;
                        }
                        else if(s[rni][rnj-1]==2){
                            b3rate=r2;
                        }
                        totalrate=b1rate+b2rate+b3rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj+1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj+1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj+1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                            if(s[rni-1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni-1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni-1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni][rnj-1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj-1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj-1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                }
                else{
                    if(rnj==0){
                        totalrate=0.0;
                        if(s[rni][rnj+1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj+1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj+1]==2){
                            b1rate=r2;
                        }
                        if(s[rni-1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni-1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni-1][rnj]==2){
                            b2rate=r2;
                        }
                        if(s[rni+1][rnj]==0){
                            b3rate=1.0;
                        }
                        else if(s[rni+1][rnj]==1){
                            b3rate=r1;
                        }
                        else if(s[rni+1][rnj]==2){
                            b3rate=r2;
                        }
                        totalrate=b1rate+b2rate+b3rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj+1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj+1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj+1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                            if(s[rni-1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni-1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni-1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni+1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni+1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni+1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                    else if(rnj==49){
                        totalrate=0.0;
                        if(s[rni][rnj-1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj-1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj-1]==2){
                            b1rate=r2;
                        }
                        if(s[rni-1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni-1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni-1][rnj]==2){
                            b2rate=r2;
                        }
                        if(s[rni+1][rnj]==0){
                            b3rate=1.0;
                        }
                        else if(s[rni+1][rnj]==1){
                            b3rate=r1;
                        }
                        else if(s[rni+1][rnj]==2){
                            b3rate=r2;
                        }
                        totalrate=b1rate+b2rate+b3rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj-1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj-1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj-1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                            if(s[rni-1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni-1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni-1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni+1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni+1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni+1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                    else{
                        totalrate=0.0;
                        if(s[rni][rnj+1]==0){
                            b1rate=1.0;
                        }
                        else if(s[rni][rnj+1]==1){
                            b1rate=r1;
                        }
                        else if(s[rni][rnj+1]==2){
                            b1rate=r2;
                        }
                        if(s[rni-1][rnj]==0){
                            b2rate=1.0;
                        }
                        else if(s[rni-1][rnj]==1){
                            b2rate=r1;
                        }
                        else if(s[rni-1][rnj]==2){
                            b2rate=r2;
                        }
                        if(s[rni][rnj-1]==0){
                            b3rate=1.0;
                        }
                        else if(s[rni][rnj-1]==1){
                            b3rate=r1;
                        }
                        else if(s[rni][rnj-1]==2){
                            b3rate=r2;
                        }
                        if(s[rni+1][rnj]==0){
                            b4rate=1.0;
                        }
                        else if(s[rni+1][rnj]==1){
                            b4rate=r1;
                        }
                        else if(s[rni+1][rnj]==2){
                            b4rate=r2;
                        }
                        totalrate=b1rate+b2rate+b3rate+b4rate;
                        if(rnb<b1rate/totalrate){
                            if(s[rni][rnj+1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj+1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj+1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                            if(s[rni-1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni-1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni-1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else if(rnb>=(b1rate+b2rate)/totalrate && rnb<(b1rate+b2rate+b3rate)/totalrate){
                            if(s[rni][rnj-1]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni][rnj-1]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if(s[rni][rnj-1]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                        else{
                            if(s[rni+1][rnj]==0){
                                if(mu1<u1){
                                    ns[rni][rnj]=1;
                                }
                                else{
                                    ns[rni][rnj]=0;
                                }
                            }
                            else if(s[rni+1][rnj]==1){
                                if(mu2<u2){
                                    ns[rni][rnj]=2;
                                }
                                else{
                                    ns[rni][rnj]=1;
                                }
                            }
                            else if (s[rni+1][rnj]==2){
                                if(mu3<u3){
                                    x3++;
                                }
                                else{
                                    ns[rni][rnj]=2;
                                }
                            }
                        }
                    }
                }
                
                
                for(i=0;i<SP;i++){
                    for(j=0;j<SP;j++){
                        s[i][j]=ns[i][j];
                    }
                }

            }//End of Moran Process for remaining premalignant cells
            
            x0=0;
            x1=0;
            x2=0;
            for(i=0;i<SP;i++){
                for(j=0;j<SP;j++){
                    if(s[i][j]==0){
                        x0++;
                    }
                    else if(s[i][j]==1){
                        x1++;
                    }
                    else if(s[i][j]==2)
                    {
                        x2++;
                    }
                }
            }
            tt[cycle]=t1+time;

            time=0;
            x3=0;
            for(t2=0;t2<100000000;t2++){
                
                // premalignant cells eventually becomes Type S cells
                
                op=rand_dist(mt); //decides whilh trial to run
                fit=rand_dist(mt); //Moran proces, decides which cell divides
                mu1=rand_dist(mt);  //decides whether to mutate from Type 0 to Type 1
                mu2=rand_dist(mt); //decides whether to mutate from Type 0 to Type S-1
                mu3=rand_dist(mt); //decides whether to mutate from Type S-1 to Type S
                dd=rand_dist(mt);   //Moran process, decided which cell dies
                rni=rand50(mt); //randomly selects position i for cell death
                rnj=rand50(mt); //radomly selects position j for cell death
                rnb=rand_dist(mt); //random number that determines replacement cell
                rn4=rand4(mt); //randomly selects replacement cell option 1
                rn3=rand3(mt); //randomly selects replacement cell option 2
                rn2=rand2(mt); //randomly selects replacement cell option 3
                
                g2=d*(double)N+r3*(double)x3+d3*(double)x3; //total rate
                
                time=time+1/g2;
                if(op>=0 && op<d*N/g2){
                    if(rni==0){
                        if(rnj==0){
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else if(rnj==49){
                            totalrate=0.0;
                            if(s[rni][rnj-1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b1rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else{
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni][rnj-1]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b3rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                    }
                    else if(rni==49){
                        if(rnj==0){
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else if(rnj==49){
                            totalrate=0.0;
                            if(s[rni][rnj-1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            totalrate=b1rate+b2rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else{
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni][rnj-1]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b3rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                    }
                    else{
                        if(rnj==0){
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b3rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else if(rnj==49){
                            totalrate=0.0;
                            if(s[rni][rnj-1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b3rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b3rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                        else{
                            totalrate=0.0;
                            if(s[rni][rnj+1]==0){
                                b1rate=1.0;
                            }
                            else if(s[rni][rnj+1]==1){
                                b1rate=r1;
                            }
                            else if(s[rni][rnj+1]==2){
                                b1rate=r2;
                            }
                            if(s[rni-1][rnj]==0){
                                b2rate=1.0;
                            }
                            else if(s[rni-1][rnj]==1){
                                b2rate=r1;
                            }
                            else if(s[rni-1][rnj]==2){
                                b2rate=r2;
                            }
                            if(s[rni][rnj-1]==0){
                                b3rate=1.0;
                            }
                            else if(s[rni][rnj-1]==1){
                                b3rate=r1;
                            }
                            else if(s[rni][rnj-1]==2){
                                b3rate=r2;
                            }
                            if(s[rni+1][rnj]==0){
                                b4rate=1.0;
                            }
                            else if(s[rni+1][rnj]==1){
                                b4rate=r1;
                            }
                            else if(s[rni+1][rnj]==2){
                                b4rate=r2;
                            }
                            totalrate=b1rate+b2rate+b3rate+b4rate;
                            if(rnb<b1rate/totalrate){
                                if(s[rni][rnj+1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj+1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj+1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=b1rate/totalrate && rnb<(b1rate+b2rate)/totalrate){
                                if(s[rni-1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni-1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni-1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else if(rnb>=(b1rate+b2rate)/totalrate && rnb<(b1rate+b2rate+b3rate)/totalrate){
                                if(s[rni][rnj-1]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni][rnj-1]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if(s[rni][rnj-1]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                            else{
                                if(s[rni+1][rnj]==0){
                                    if(mu1<u1){
                                        ns[rni][rnj]=1;
                                    }
                                    else{
                                        ns[rni][rnj]=0;
                                    }
                                }
                                else if(s[rni+1][rnj]==1){
                                    if(mu2<u2){
                                        ns[rni][rnj]=2;
                                    }
                                    else{
                                        ns[rni][rnj]=1;
                                    }
                                }
                                else if (s[rni+1][rnj]==2){
                                    if(mu3<u3){
                                        x3++;
                                    }
                                    else{
                                        ns[rni][rnj]=2;
                                    }
                                }
                            }
                        }
                    }
                }
                else if(op>=d*N/g2 && op<(d*N/g2+r3*x3/g2)){
                    x3++;
                }
                else if(op>=(d*N/g2+r3*x3/g2) && op<=1){
                    x3--;
                }
                for(i=0;i<SP;i++){
                    for(j=0;j<SP;j++){
                        s[i][j]=ns[i][j];
                    }
                }

                if(x3>=detect_size){
                    x3=0;
                    break;
                }
                
            }//End.
            
            time=time+(1/(r3-d3))*log((double)actual_size/(double)detect_size);
            tr[cycle]=time; //time to recurence
            //            fprintf(fp,"%f %d\n",time,abort);
            
        }
        std::sort(tr,tr+SRUN);
        Q=0.0;
        SQ=0.0; // mean squared logarithmic residuals
        epsilon=0.0; // minimum value for mean squared logarithmic residuals
        counter=0;
        for(a=0;a<SRUN;a++){
            if(databox[a]>0.01){
                Q+=std::fabs(log2(databox[a]/tr[a]));
                counter++;
            }
            
        }
        epsilon=1.5; // maximum threshold for mean squared logarithmic residuals
        SQ=Q/(double)counter;
        if(SQ<epsilon){
            fprintf(fp,"N=%d detect_size=%d actual_size=%d r0=%f r1=%f r2=%f r3=%f d=%f d3=%f u1=%f u2=%f u3=%f SQ=%f\n",N,detect_size,actual_size,r0,r1,r2,r3,d,d3,u1,u2,u3,SQ);
        }
        for(a=0;a<SRUN;a++){
            fprintf(fp2,"%f %f\n",databox[a],tr[a]);
        }
        
        
    }
    fclose(fp);
    fclose(fp2);
    return 0;
}
