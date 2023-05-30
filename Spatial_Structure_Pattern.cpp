//
// Spatial Structure Initiation Pattern (For Figure 4 in Abubakar et al. . Frontiers in Oncology, 2023))
//

// Executes 2 file names. fp is cell numbers for each cell type at specific time point while fp2 is parameter combination

#include<iostream>
#include<fstream>
#include<random>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>

// Constant Variables
#define N 2500 // Cell number or tissue size
#define SP 50 //  Spatial structure grid size
#define detect_size 5000
#define actual_size 1000000000 // tissue size for cancer diagnosis
#define r0 1.0 // fitness of type 0 cell
#define d 1.0 // Death rate of Type 0 cell
#define d3 1.0 // Death rate of Type S cell

using namespace std;

int main(void){
    // new introduction for spatial structure
    int s[SP][SP],ns[SP][SP]; // a grid of SP x SP cells
    int i,j; // row number and column number (position) of a cell in a I x J grid
    for(i=0;i<SP;i++){
        for(j=0;j<SP;j++){
            s[i][j]=0;
            ns[i][j]=0;
        }
    }
    // Modifiable conditions
    double r1,r2,r3,u1,u2,u3,po;
    
    // Computer randomized parameters
    double op,fit,mu1,mu2,mu3,dd,rnb;
    int rni,rnj,rn4,rn3,rn2;
    
    // Variables and Types
    int x0,x1,x2,x3; // Number of Type 0, Type 1, Types S-1 and Type S cells
    int cycle,t1,t0,t2,a;
    double time;
    double timecount=1;
    double g1,g2;
    double Et,Er,sVt,sVr,Vt,Vr,St,Sr,it[10000],ir[10000],tt[10000],tr[10000];
    int sumtt,sumtr;
    int pro[3];
    double prop[4];
    int alpha,beta,gamma,delta;
    double totalrate; //for cell division
    double b1rate,b2rate,b3rate,b4rate;
    
    // Name output file (ends with .dat)
    FILE *fp;
    FILE *fp2;
    fp=fopen("4pmodel1sp_zzbbb1.dat","w");
    fp2=fopen("pars_4pmodelsp_zzbbb.dat","w");

    std::random_device rnd; // Random number function (rnd is variable name）
    std::mt19937_64 mt(rnd());// Random number generator
    std::uniform_real_distribution<> rand_dist(0,1); // possibility of a choice
    std::uniform_int_distribution<> rand50(0, 49); // Choice for cell death
    std::uniform_int_distribution<> rand4(0, 3); // 1st choice for cell division
    std::uniform_int_distribution<> rand3(0, 2); // 2nd choice for cell division
    std::uniform_int_distribution<> rand2(0, 1); // 3rd choice for cell division
    
    
    r1=1.25; // Fitness of Type 1 cell
    r2=1.25; // Fitness of Type S-1 cell
    r3=1.50; // Fitness of Type S cell
    
    u1=pow(10,-2); // mutation rate from Type 0 to Type 1 (in log 10 values)
    u2=pow(10,-2); // mutation rate from Type 1 to Type S-1 (in log 10 values)
    u3=pow(10,-2); // mutation rate from Type S-1 to Type S (in log 10 values)
    
    fprintf(fp2,"N=%d detect_size=%d actual_size=%d r0=%f r1=%f r2=%f r3=%f d=%f d3=%f u1=%f u2=%f u3=%f\n",N,detect_size,actual_size,r0,r1,r2,r3,d,d3,u1,u2,u3);
    
    for(i=0;i<3;i++){
        pro[i]=0;
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
    
    x3=0;
    time=0;

    for(t1=0;t1<100000000;t1++){
        op=rand_dist(mt); //decides which trial to run
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
        else if(op>=d*N/g2 && op<d*N/g2+r3*x3/g2){
            x3++;
        }
        else if(op>=d*N/g2+r3*x3/g2 && op<=1){
            x3--;
        }
        for(i=0;i<SP;i++){
            for(j=0;j<SP;j++){
                s[i][j]=ns[i][j];
            }
        }
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

        if(timecount>=50000){
            break;
        }
        if(x3>=20000000){
            break;
        }
        if(time>=timecount){
            fprintf(fp,"%f %d %d %d %d\n",timecount,x0,x1,x2,x3);
            timecount+=1.0;
        }
    }
            
    fclose(fp);
    fclose(fp2);
    return 0;
}
