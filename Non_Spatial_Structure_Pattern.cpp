//
// Non Spatial Structure Initiation Pattern (For Figure 2 in Abubakar et al.)
//

//Executes 2 file names. fp is cell numbers for each cell type at specific time point while fp2 is parameter combination

#include<iostream>
#include<fstream>
#include<random>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>

// Constant Variables
#define N 1000 // Cell number or tissue size
#define detect_size 2000 
#define actual_size 1000000000 // tissue size for cancer diagnosis
#define r0 1.0 // fitness of type 0 cell
#define d 1.0 // Death rate of Type 0 cell
#define d3 1.0 // Death rate of Type S cell


using namespace std;

int main(void){
        
   // Modifiable conditions
    double r1,r2,r3,u1,u2,u3,po;
    
    // Computer randomized parameters
    double op,fit,mu1,mu2,mu3,dd;
    
    // Variables and Types
    int x0,x1,x2,x3; // Number of Type 0, Type 1, Types S-1 and Type S cells
    int i,cycle,t1,t0,t2,a;
    double time;
    double timecount=1;
    double g1,g2;
    double Et,Er,sVt,sVr,Vt,Vr,St,Sr,it[10000],ir[10000],tt[10000],tr[10000];
    int sumtt,sumtr;
    int pro[3];
    double prop[4];
    
    // Input file name below (ends with .dat)
   FILE *fp;
   FILE *fp2;
         fp=fopen("4pmodelaaxux_1.dat","w");
   fp2=fopen("pars_4pmodelaaxux.dat","w");
    
    std::random_device rnd; //Random number function (rnd is variable name）
    std::mt19937_64 mt(rnd());//Random number generator
    std::uniform_real_distribution<> rand_dist(0,1);
    
    r1=0.75; // Fitness of Type 1 cell
    r2=0.75; // Fitness of Type S-1 cell
    r3=1.5; // Fitness of Type S cell
    
    u1=pow(10,-3); // mutation rate from Type 0 to Type 1 (in log 10 values)
    u2=pow(10,-4); // mutation rate from Type 1 to Type S-1 (in log 10 values)
    u3=pow(10,-3); // mutation rate from Type S-1 to Type S (in log 10 values)

    
    fprintf(fp2,"N=%d detect_size=%d actual_size=%d r0=%f r1=%f r2=%f r3=%f d=%f d3=%f u1=%f u2=%f u3=%f\n",N,detect_size,actual_size,r0,r1,r2,r3,d,d3,u1,u2,u3);

        
        for(i=0;i<3;i++){
            pro[i]=0;
        }
        
        

            
            x0=N;
            x1=0;
            x2=0;
            x3=0;
            time=0;
            
            for(t1=0;t1<1000000000;t1++){
                
                
                op=rand_dist(mt); //decides which trial to run
                fit=rand_dist(mt); //Moran proces, decides which cell divides
                mu1=rand_dist(mt); //decides whether to mutate from Type 0 to Type 1
                mu2=rand_dist(mt); //decides whether to mutate from Type 0 to Type S-1
                mu3=rand_dist(mt); //decides whether to mutate from Type S-1 to Type S
                dd=rand_dist(mt);  //Moran process, decided which cell dies
                
                g1=r0*(double)x0+r1*(double)x1+r2*(double)x2; //total fitness
                g2=d*(double)N+r3*(double)x3+d3*(double)x3; //total rate
                
                time=time+1/g2;
                
                if(op>=0 && op<d*N/g2){
                    
                    
                    if(fit>=0 && fit<=r0*(double)x0/g1){
                        
                        if(mu1>=0 && mu1<u1 && dd>=0 && dd<=(double)x0/N){
                            
                            x1++;
                            x0--;
                            
                        }
                        
                        else if(mu1>=0 && mu1<u1 && dd>=(double)x0/N && dd<(double)x0/N+(double)x2/N){
                            
                            x1++;
                            x2--;
                            
                        }
                        
                        else if(mu1<=1 && mu1>=u1 && dd<=1 && dd>(double)x0/N && dd<(double)x0/N+(double)x2/N){
                            
                            x0++;
                            x2--;
                            
                        }
                        
                        else if(mu1<1 && mu1>u1 && dd<=1 && dd>(double)x0/N+(double)x1/N){
                            
                            x0++;
                            x1--;
                            
                        }
                        
                    }
                    
                    else if(fit>r0*(double)x0/g1 && fit<=r0*(double)x0/g1+r1*(double)x1/g1 && fit<=1){
                        
                        if(mu2>=0 && mu2<u2 && dd>=0 && dd<=(double)x0/N){
                            
                            x2++;
                            x0--;
                            
                        }
                        
                        else if(mu2>=0 && mu2<u2 && dd<=1 && dd>=(double)x0/N+(double)x2/N){
                            
                            x2++;
                            x1--;
                            
                        }
                        
                        else if(mu2<1 && mu2>u2 && dd>0 && dd<=(double)x0/N){
                            
                            x1++;
                            x0--;
                            
                        }
                        
                        else if(mu2<1 && mu2>u2 && dd>=(double)x0/N && dd<(double)x0/N+(double)x2/N){
                            
                            x1++;
                            x2--;
                            
                        }
                        
                    }
                    
                    else if(fit<=1 && fit>r0*(double)x0/g1+r1*(double)x1/g1){
                        
                        if(mu3<1 && mu3>u3 && dd>0 && dd<=(double)x0/N){
                            
                            x2++;
                            x0--;
                            
                        }
                        
                        else if(mu3<1 && mu3>u3 && dd<=1 && dd>=(double)x0/N+(double)x2/N){
                            
                            x2++;
                            x1--;
                            
                        }
                       
                        else if(mu3>0 && mu3<u3){
                            
                            x3++;
                            
                        }
                    }
                    
                }
                
                else if(op>=d*N/g2 && op<d*N/g2+r3*x3/g2){
                    
                    x3++;
                    
                }
                else if(op>=d*N/g2+r3*x3/g2 && op<=1){
                    
                    x3--;
                    
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
