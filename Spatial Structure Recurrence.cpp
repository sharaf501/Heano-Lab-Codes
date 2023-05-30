//
//  Spatial Structure Recurrence (For Figure 5 in Abubakar et al. . Frontiers in Oncology, 2023))
//  

//  Executes 2 file names. fp is mean +/- SD of recurrence and proportion of Type S-1 cells for each parameter while fp2 is parameter combination

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
#define tCYCLE 1000 // number of iterations

using namespace std;

int main(void){
    // new introduction for spatial structure
    int s[SP][SP],ns[SP][SP]; // a grid of SP x SP cells
    int i,j; // row number and column number (position) of a cell in a I x J grid
   
    // Modifiable conditions
    double r1,r2,r3,u1,u2,u3,po;
   
    // Computer randomized parameters
    double op,fit,mu1,mu2,mu3,dd,rnb;
    int rni,rnj,rn4,rn3,rn2;
    
    // Variables and Types
    int x0,x1,x2,x3;  // Number of Type 0, Type 1, Types S-1 and Type S cells
    int cycle,t1,t0,t2,a;
    double time;
    double g1,g2;
    double Et,Er,sVt,sVr,Vt,Vr,St,Sr,it[10000],ir[10000],tt[10000],tr[10000];
    int sumtt;
    double sumtr;
    int pro[3],bpro[3];
    double prop[3],bprop[3];
    int alpha,beta,gamma,delta;
    double totalrate; //for cell division
    double b1rate,b2rate,b3rate,b4rate;

    // Name output file (ends with .dat)
   FILE *fp;
   FILE *fp2;
         fp=fopen("r1_4pmodel_sp_2.5k.dat","w");
   fp2=fopen("pars_r1_4pmodel_sp_2.5k.dat","w");
    
    std::random_device rnd; // Random number function (rnd is variable name）
    std::mt19937_64 mt(rnd()); // Random number generator
    std::uniform_real_distribution<> rand_dist(0,1); // possibility of a choice
    std::uniform_int_distribution<> rand50(0, 49); // Choice for cell death
    std::uniform_int_distribution<> rand4(0, 3); // 1st choice for cell division
    std::uniform_int_distribution<> rand3(0, 2); // 2nd choice for cell division
    std::uniform_int_distribution<> rand2(0, 1); // 2nd choice for cell division
    
    
    r1=1.0; // Fitness of Type 1 cell
    r2=1.0; // Fitness of Type S-1 cell
    r3=1.25; // Fitness of Type S cell
    
    u1=pow(10,-3); // mutation rate from Type 0 to Type 1 (in log 10 values)
    u2=pow(10,-3); // mutation rate from Type 1 to Type S-1 (in log 10 values)
    u3=pow(10,-3); // mutation rate from Type S-1 to Type S (in log 10 values)
    
    fprintf(fp2,"N=%d detect_size=%d actual_size=%d r0=%f r1=%f r2=%f r3=%f d=%f d3=%f u1=%f u2=%f u3=%f\n",N,detect_size,actual_size,r0,r1,r2,r3,d,d3,u1,u2,u3);
    

    // adjust range of parameter for dependence on recurrence time
    // Run ONLY one parameter at a time
    // Ensure output variable is same as currently analyzed variable
    
    for(r1=0.90;r1<1.11;r1+=0.02){
    // for(r2=0.90;r2<1.11;r2+=0.02){
    // for(r3=1.05;r3<1.51;r3+=0.05){
    // for(po=-5;po<-1.9;po+=0.3){
    // u1=(double)pow(10.0,po);
    // u2=(double)pow(10.0,po);
    // u3=(double)pow(10.0,po);
        
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
            tr[a]=0.0;
            tt[a]=0.0;
            it[a]=0.0;
            ir[a]=0.0;
        }
        for(cycle=0;cycle<tCYCLE;cycle++){
            
            for(i=0;i<SP;i++){
                for(j=0;j<SP;j++){
                    s[i][j]=0;
                    ns[i][j]=0;
                }
            }

            x3=0;
            time=0;
            
            for(t1=0;t1<1000000000;t1++){
                
                
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
                
 //               g1=r0*(double)x0+r1*(double)x1+r2*(double)x2; //total fitness
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

                if(x3>=detect_size){
                    x3=0;
                    break;
                }
                if(t1==999999999){
                    break;
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

            //Moran process continues for remaining premalignant cells
            
            time=d*(double)N*(1/(r3-d))*log((double)actual_size/(double)detect_size);
            
            for(t0=0; t0<(int)time ; t0++){
                
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
            it[cycle]=x2;
            
            if(it[cycle]>=0 && it[cycle]<=(0.1*N)){
                pro[0]++;
            }
            else if(it[cycle]>(0.1*N) && it[cycle]<=(0.9*N)){
                pro[1]++;
            }
            else if(it[cycle]>(0.9*N) && it[cycle]<=N){
                pro[2]++;
            }
            
            time=0;
            x3=0;
            for(t2=0;t2<1000000000;t2++){
                
                // t2 is event count from detect size to actual size
                // premalignant cells eventually becomes Type S cells
                
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

                if(x3>=detect_size){
                    x3=0;
                    break;
                }
                
            }//End. 
            
            time=time+(1/(r3-d))*log((double)actual_size/(double)detect_size);
            tr[cycle]=time; //time to recurence

            
        }
        
        sumtr=0.0;
        for(a=0;a<tCYCLE;a++){
            sumtr=sumtr+tr[a];
         }
         Er=sumtr/(double)tCYCLE; // Mean of recurrence time for a given number of cycles
        sVr=0.0;
        for(a=0;a<tCYCLE;a++){
         sVr=sVr+(tr[a]-Er)*(tr[a]-Er);
         }
         
         Vr=sVr/(double)tCYCLE;
         
         Sr=sqrt(Vr); // SD of recurrence time for a given number of cycles
        
       
        prop[0]=(double)pro[0]/(double)cycle; // proportion of Type 0 cells
        prop[1]=(double)pro[1]/(double)cycle; // proportion of Type 1 cells
        prop[2]=(double)pro[2]/(double)cycle; // proportion of Type S-1 cells

       
        
        fprintf(fp,"%f %f %f %f %f %f\n",r1,Er,Sr,prop[0],prop[1],prop[2]);
        
    
    
    }
    fclose(fp);
    fclose(fp2);
    return 0;
    
    
    
}
