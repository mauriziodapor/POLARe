/* SPAS.c
 Spin-polarization after scattering
 Maurizio Dapor, May 2022 */
    
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_sf.h>

#define pi 3.14159265359

double
theta[200000],S[200000],T[200000],U[200000],Pfx[200000],Pfy[200000],Pfz[200000],Pf[200000];
double c;
double phi,Px,Py,Pz,P;
double cphi,sphi;
long enne,nmax;
char inS[50],inT[50],inU[50],inphi[50],inPx[50],inPy[50],inPz[50];
FILE *cid,*cidS,*cidT,*cidU,*codPfx,*codPfy,*codPfz,*codPf;

int main ()
{
    //Input: File names, azimuth, and initial spin-polarization
    
    cid=fopen("INPUT/input.txt","r");
    fscanf(cid,"%s\n %s\n %s\n %s %lg\n %s %lg\n %s %lg\n %s %lg\n",&inS,&inT,&inU,&inphi,&phi,&inPx,&Px,inPy,&Py,inPz,&Pz);
    fclose(cid);
    
    //Initial spin polarization
        
        P=sqrt(Px*Px+Py*Py+Pz*Pz);
        if(P>1.)
        {
            printf("\n|P| must be < 1\n\n");
            exit(0);
        }
    
    //Input: S function)
    
    enne=0;
    cidS=fopen(inS,"r");
    while( (c=fgetc(cidS)) !=EOF)
    {
        enne=enne+1;
        fscanf(cidS,"%lg%lg",&theta[enne],&S[enne]);
    }
    fclose(cidS);
    
    //Input: T function)
    
    enne=0;
    cidT=fopen(inT,"r");
    while( (c=fgetc(cidT)) !=EOF)
    {
        enne=enne+1;
        fscanf(cidT,"%lg%lg",&theta[enne],&T[enne]);
    }
    fclose(cidT);
    
    //Input: U function)
    
    enne=0;
    cidU=fopen(inU,"r");
    while( (c=fgetc(cidU)) !=EOF)
    {
        enne=enne+1;
        fscanf(cidU,"%lg%lg",&theta[enne],&U[enne]);
    }
    fclose(cidU);
    nmax=enne;
    
    phi=pi*phi/180.;
    cphi=cos(phi);
    sphi=sin(phi);
    
    for(enne=1;enne<nmax;enne++)
    {
        Pfx[enne]=((-S[enne]+Px*sphi-Py*cphi)*sphi+Px*T[enne]+(Py*cphi-Px*sphi)*sphi*T[enne]+Pz*cphi*U[enne])/(1+S[enne]*(Py*cphi-Px*sphi));
        Pfy[enne]=((S[enne]-Px*sphi+Py*cphi)*cphi+Py*T[enne]+(Px*sphi-Py*cphi)*cphi*T[enne]+Pz*U[enne]*sphi)/(1+S[enne]*(Py*cphi-Px*sphi));
        Pfz[enne]=(Pz*T[enne]-(Px*cphi+Py*sphi)*U[enne])/(1+S[enne]*(Py*cphi-Px*sphi));
        Pf[enne]=sqrt(Pfx[enne]*Pfx[enne]+Pfy[enne]*Pfy[enne]+Pfz[enne]*Pfz[enne]);
    }
    
    codPfx=fopen("RESULTS/Pfx.txt","w");
    for(enne=1;enne<nmax;enne++)
    {
        fprintf(codPfx,"%lg %lg\n",theta[enne],Pfx[enne]);
    }
    fclose(codPfx);
    
    codPfy=fopen("RESULTS/Pfy.txt","w");
    for(enne=1;enne<nmax;enne++)
    {
        fprintf(codPfy,"%lg %lg\n",theta[enne],Pfy[enne]);
    }
    fclose(codPfy);
    
    codPfz=fopen("RESULTS/Pfz.txt","w");
    for(enne=1;enne<nmax;enne++)
    {
        fprintf(codPfz,"%lg %lg\n",theta[enne],Pfz[enne]);
    }
    fclose(codPfz);
    
    codPf=fopen("RESULTS/Pf.txt","w");
    for(enne=1;enne<nmax;enne++)
    {
        fprintf(codPf,"%lg %lg\n",theta[enne],Pf[enne]);
    }
    fclose(codPf);
    
    printf("%lg %lg %lg %lg %lg\n",cphi,sphi,Px,Py,Pz);
    
    //Gnuplot figures
        
    FILE *gnuplotPfx = fopen("FIGURES/commandsPfx.txt","w");
    fprintf(gnuplotPfx, "set terminal pdf\n");
    fprintf(gnuplotPfx, "set output 'FIGURES/Pfx.pdf' \n");
    fprintf(gnuplotPfx, "set title 'Spin-polarization after scattering' \n");
    fprintf(gnuplotPfx, "unset key \n");
    fprintf(gnuplotPfx, "set termoption enhanced\n");
    fprintf(gnuplotPfx, "set xlabel '{/Symbol:Italic q} (deg) \n");
    fprintf(gnuplotPfx, "set ylabel 'P^f_x({/Symbol:Italic q})' \n");
    fprintf(gnuplotPfx, "set xrange [0:180] \n");
    fprintf(gnuplotPfx, "plot 'RESULTS/Pfx.txt' w l \n");
    fflush(gnuplotPfx);
    fclose(gnuplotPfx);
    system("gnuplot 'FIGURES/commandsPfx.txt'");
    
    FILE *gnuplotPfy = fopen("FIGURES/commandsPfy.txt","w");
    fprintf(gnuplotPfy, "set terminal pdf\n");
    fprintf(gnuplotPfy, "set output 'FIGURES/Pfy.pdf' \n");
    fprintf(gnuplotPfy, "set title 'Spin-polarization after scattering' \n");
    fprintf(gnuplotPfy, "unset key \n");
    fprintf(gnuplotPfy, "set termoption enhanced\n");
    fprintf(gnuplotPfy, "set xlabel '{/Symbol:Italic q} (deg) \n");
    fprintf(gnuplotPfy, "set ylabel 'P^f_y({/Symbol:Italic q})' \n");
    fprintf(gnuplotPfy, "set xrange [0:180] \n");
    fprintf(gnuplotPfy, "plot 'RESULTS/Pfy.txt' w l \n");
    fflush(gnuplotPfy);
    fclose(gnuplotPfy);
    system("gnuplot 'FIGURES/commandsPfy.txt'");
    
    FILE *gnuplotPfz = fopen("FIGURES/commandsPfz.txt","w");
    fprintf(gnuplotPfz, "set terminal pdf\n");
    fprintf(gnuplotPfz, "set output 'FIGURES/Pfz.pdf' \n");
    fprintf(gnuplotPfz, "set title 'Spin-polarization after scattering' \n");
    fprintf(gnuplotPfz, "unset key \n");
    fprintf(gnuplotPfz, "set termoption enhanced\n");
    fprintf(gnuplotPfz, "set xlabel '{/Symbol:Italic q} (deg) \n");
    fprintf(gnuplotPfz, "set ylabel 'P^f_z({/Symbol:Italic q})' \n");
    fprintf(gnuplotPfz, "set xrange [0:180] \n");
    fprintf(gnuplotPfz, "plot 'RESULTS/Pfz.txt' w l \n");
    fflush(gnuplotPfz);
    fclose(gnuplotPfz);
    system("gnuplot 'FIGURES/commandsPfz.txt'");
    
    FILE *gnuplotPf = fopen("FIGURES/commandsPf.txt","w");
    fprintf(gnuplotPf, "set terminal pdf\n");
    fprintf(gnuplotPf, "set output 'FIGURES/Pf.pdf' \n");
    fprintf(gnuplotPf, "set title 'Spin-polarization after scattering' \n");
    fprintf(gnuplotPf, "unset key \n");
    fprintf(gnuplotPf, "set termoption enhanced\n");
    fprintf(gnuplotPf, "set xlabel '{/Symbol:Italic q} (deg) \n");
    fprintf(gnuplotPf, "set ylabel 'P^f({/Symbol:Italic q})' \n");
    fprintf(gnuplotPf, "set xrange [0:180] \n");
    fprintf(gnuplotPf, "plot 'RESULTS/Pf.txt' w l \n");
    fflush(gnuplotPf);
    fclose(gnuplotPf);
    system("gnuplot 'FIGURES/commandsPf.txt'");
}
    
    
    
