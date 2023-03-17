//
//  Non Spatial Structure Recurrrence (For Figure 3 in Abubakar et al.)
//

//Executes 2 file names. fp is mean +/- SD of recurrence and proportion of Type S-1 cells for each parameter while fp2 is parameter combination

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
#define tCYCLE 1000 //  number of iterations


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
    double g1,g2;
    double Et,Er,sVt,sVr,Vt,Vr,St,Sr,it[10000],ir[10000],tt[10000],tr[10000];
    int sumtt,sumtr;
    int pro[3];
    double prop[4];
    
    // Input file name below (ends with .dat)
    FILE *fp;
   FILE *fp2;
         fp=fopen("u1_4pmodel1k.dat","w");
   fp2=fopen("pars_u1_4pmodel1k.dat","w");
    
    std::random_device rnd; // Random number function (rnd is variable name）
    std::mt19937_64 mt(rnd());// Random number generator
    std::uniform_real_distribution<> rand_dist(0,1);
    
    r1=1.0; // Fitness of Type 1 cell
    r2=1.0; // Fitness of Type S-1 cell
    r3=1.25; // Fitness of Type S cell
    
    u1=pow(10,-3); // mutation rate from Type 0 to Type 1 (in log 10 values)
    u2=pow(10,-3); // mutation rate from Type 1 to Type S-1 (in log 10 values)
    u3=pow(10,-3); // mutation rate from Type S-1 to Type S (in log 10 values)
    
    
    fprintf(fp2,"N=%d detect_size=%d actual_size=%d r0=%f r1=%f r2=%f r3=%f d=%f d3=%f u1=%f u2=%f u3=%f\n",N,detect_size,actual_size,r0,r1,r2,r3,d,d3,u1,u2,u3);
    
    
    // adjust range of parameter for dependence on recurrence time
    
    for(po=-5;po<-1.9;po+=0.3){
    // for(r2=0.90;r2<1.11;r2+=0.02){
    // for(r3=1.05;r3<1.51;r3+=0.05){
    u1=(double)pow(10.0,po);
        
        for(i=0;i<tCYCLE;i++){
            it[i]=0.0;
            tr[i]=0.0;
        }

        sumtt=0;
        sumtr=0;
        sVt=0;
        sVr=0;
        
        for(i=0;i<3;i++){
            pro[i]=0;
        }
        
        
        for(cycle=0;cycle<tCYCLE;cycle++){
            
            x0=N;
            x1=0;
            x2=0;
            x3=0;
            time=0;
            
            for(t1=0;t1<1000000000;t1++){
                
                
                op=rand_dist(mt); //decides which trial to run
                fit=rand_dist(mt); //Moran proces, decides which cell divides
                mu1=rand_dist(mt);  //decides whether to mutate from Type 0 to Type 1
                mu2=rand_dist(mt); //decides whether to mutate from Type 0 to Type S-1
                mu3=rand_dist(mt); //decides whether to mutate from Type S-1 to Type S
                dd=rand_dist(mt);   //Moran process, decided which cell dies
                
                g1=r0*(double)x0+r1*(double)x1+r2*(double)x2; //total fitness
                g2=d*(double)N+r3*(double)x3+d3*(double)x3; //total rate
                
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
                
                if(x3>=detect_size){
                    
                    x3=0;
                   break;
                    
                }
                if(t1==999999999){
                    x3=0;
                }
            }
            
            
            //Moran process continues for remaining premalignant cells
            
            time=(1/(r3-d))*log((double)actual_size/(double)detect_size);
            
            for(t0=0; t0<(int)time*N*(int)d ; t0++){
                
                op=rand_dist(mt); //decides which trial to run
                fit=rand_dist(mt); //Moran proces, decides which cell divides
                mu1=rand_dist(mt);  //decides whether to mutate from Type 0 to Type 1
                mu2=rand_dist(mt); //decides whether to mutate from Type 0 to Type S-1
                mu3=rand_dist(mt); //decides whether to mutate from Type S-1 to Type S
                dd=rand_dist(mt);   //Moran process, decided which cell dies
                
                g1=r0*(double)x0+r1*(double)x1+r2*(double)x2; //total fitness
                    
                    //Type0 divides
                
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
                    
                    //Type1 divides
                
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
                    
                //Type 2 divides
                
                else if(fit<=1 && fit>r0*(double)x0/g1+r1*(double)x1/g1){
                    
                    if(mu3<1 && mu3>u3 && dd>0 && dd<=(double)x0/N){
                        
                        x2++;
                        x0--;
                        
                    }
                    
                    else if(mu3<1 && mu3>u3 && dd<=1 && dd>=(double)x0/N+(double)x2/N){
                        
                        x2++;
                        x1--;
                        
                    }
                   
                    
                }
                
            }//End of Moran Process for remaining premalignant cells
            
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
            for(t2=0;t2<1000000000;t2++){
              
                // t2 is event count after surgery
                
                op=rand_dist(mt); //decides which trial to run
                fit=rand_dist(mt); //Moran proces, decides which cell divides
                mu1=rand_dist(mt);  //decides whether to mutate from Type 0 to Type 1
                mu2=rand_dist(mt); //decides whether to mutate from Type 0 to Type S-1
                mu3=rand_dist(mt); //decides whether to mutate from Type S-1 to Type S
                dd=rand_dist(mt);   //Moran process, decided which cell dies
                
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
                
                if(x3>=detect_size){
                    
                    x3=0;
                    break;
                }
                
            }//End. 
            
            time=time+(1/(r3-d))*log((double)actual_size/(double)detect_size);
            tr[cycle]=time; //time to recurence

            
        }
        
        for(a=0;a<tCYCLE;a++){
         sumtr=sumtr+tr[a];
         }
         Er=(double)sumtr/(double)tCYCLE; // Mean of recurrence time for a given number of cycles
        
        for(a=0;a<tCYCLE;a++){
         sVr=sVr+(tr[a]-Er)*(tr[a]-Er);
         }
         
         Vr=sVr/(double)tCYCLE;
         
         Sr=sqrt(Vr); // SD of recurrence time for a given number of cycles
        
        prop[0]=(double)pro[0]/(double)tCYCLE; // proportion of Type 0 cells
        prop[1]=(double)pro[1]/(double)tCYCLE; // proportion of Type 1 cells
        prop[2]=(double)pro[2]/(double)tCYCLE; // proportion of Type S-1 cells
        
        
        fprintf(fp,"%f %f %f %f %f %f\n",u1,Er,Sr,prop[0],prop[1],prop[2]);
        
    
    
    }
    fclose(fp);
    fclose(fp2);
    
    return 0;
    
    
}
