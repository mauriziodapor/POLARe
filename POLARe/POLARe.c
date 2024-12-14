/* POLARe.c

 Relativistic partial wave solution of quantum scattering. Dirac equation. Electrons and positrons projectiles.
 Screening function: Cox and Bonham (1967) (Z<=54), Salvat et al. (1987) (Z>54)
 Exchange potential energy (only for electron projectiles): Furness and McCarty (1973)
 Statistical mixture of spin states: Kessler (1985), Burke and Joachain (1995)
 
 Maurizio Dapor, March 2022 */

#include<stdio.h>
#include<math.h>
#include <gsl/gsl_sf.h>

#define azero 0.529177 //A - azero=hbar^2/(m e^2)
#define ebye 14.3998 // eV A - ebye=e^2
#define mcbyc 511004.06 //eV - mcbyc\,=\,m c^2
#define hbarc 1973.2698 //eV A - hbarc=\,hbar c
#define hbardbymc 0.0038616 //A - hbardbymc=hbar/(m c)
#define pi 3.14159265359
#define h 0.01

char in1[30],in2[30],in3[30],in4[30],in5[30],in6[30],in7[30];
double A1,A2,A3,A4,A5,A6,A7,A8,A9,A10;
double alp1,alp2,alp3,alp4,alp5,alp6,alp7,alp8,alp9,alp10;
double Z,T,rho,Z0,Z1,Z2,Z3,rd,rmax,zz;
double x,j,jj,n,nn,E,Er,W,Z0,Z1,Z2,Z3,p0,p1,p2,p3,K,KF,KL,pl,rmax,r[1000000],V[1000000],T1[1000000],T2[1000000],th[100000],sd[100000];
double c,td1,td2;
double azimuth,Px,Py,Pz,P,Pn;
long enne,nmax,emme,mmax,mtheta,jei,l,spin;
double ct,st,pol,pol1,thetamin,thetamax,dtheta,theta,sigmad,sigmadp,Sf,Tf,Uf;
double cd1,cd2,sd1,sd2,sum1,sum2,sum3,sum4;
double dth,st,stra,stra2;
double pel;
double rcp,vrcp,DELTA,epsilon,ptb,ptco;
FILE *cid1,*cid2,*cid3,*codSD,*codSf,*codTf,*codUf,*codST,*codSDp,*codP,*codAD;
long electron;
char inputfile[50],inputfilep[50],descs[50],Sfunction[50],Tfunction[50],Ufunction[50],descsp[50],total[50],probability[50],density[50];

//Potential energy

double vi(double s)
{
	double phi,pt,phi2,ptexch,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;
	e1=exp(-alp1*s);
	e2=exp(-alp2*s);
	e3=exp(-alp3*s);	
	e4=exp(-alp4*s);
	e5=exp(-alp5*s);
	e6=exp(-alp6*s);
    e7=exp(-alp7*s);
    e8=exp(-alp8*s);
    e9=exp(-alp9*s);
    e10=exp(-alp10*s);
    
//Screening function (Z<=54: Cox and Bonham, Z>54 Salvat et al.)
    
    phi=A1*e1+A2*e2+A3*e3+A4*e4+A5*e5+A6*e6+A7*e7+A8*e8+A9*e9+A10*e10;
phi2=(alp1*alp1*A1*e1+alp2*alp2*A2*e2+alp3*alp3*A3*e3+alp4*alp4*A4*e4+alp5*alp5*A5*e5+alp6*alp6*A6*e6+alp7*alp7*A7*e7+alp8*alp8*A8*e8+alp9*alp9*A9*e9+alp10*alp10*A10*e10)/(hbardbymc*hbardbymc);
    
//Atomic density in A^-3
    
    rho=Z*phi2/(4.*pi*s*hbardbymc);
    
//Electrostatic potential energy in eV (positrons)
    
	pt=Z*ebye*phi/(s*hbardbymc);
    
//Exchange potential energy in eV (Furness and McCarty, 1973)
    
	ptexch=0.5*((E+pt)-sqrt((E+pt)*(E+pt)+4.*pi*rho*azero*ebye*ebye));
    
//Potential energy in eV
    
    if (electron==1) //electron beam
    {
        return (-pt+ptexch)/mcbyc;
        
    }
    else //positron beam
    {
        return pt/mcbyc;
    }
}

//Atomic density

double ad(double s)
{
    double phi2,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;
    e1=exp(-alp1*s);
    e2=exp(-alp2*s);
    e3=exp(-alp3*s);
    e4=exp(-alp4*s);
    e5=exp(-alp5*s);
    e6=exp(-alp6*s);
    e7=exp(-alp7*s);
    e8=exp(-alp8*s);
    e9=exp(-alp9*s);
    e10=exp(-alp10*s);
phi2=(alp1*alp1*A1*e1+alp2*alp2*A2*e2+alp3*alp3*A3*e3+alp4*alp4*A4*e4+alp5*alp5*A5*e5+alp6*alp6*A6*e6+alp7*alp7*A7*e7+alp8*alp8*A8*e8+alp9*alp9*A9*e9+alp10*alp10*A10*e10)/(hbardbymc*hbardbymc);
    rho=Z*phi2/(4.*pi*s*hbardbymc); //atomic density
    return rho*4*pi*s*s*(hbardbymc*hbardbymc);
}

//Definition of the differential equation

double fn(long mm,double t)
{
    double s,vi,f;
    s=r[mm];
    vi=V[mm];
    f=(spin/s)*sin(2.*t)+W-vi-cos(2.*t);
    return f;
}

//Solution of the differential equation (by using the fourth-order Runge-Kutta method)

double phaseshift()
{
    double k1,k2,k3,k4,num,den,tandelta;
    enne=-1;
    while (enne<nmax-4)
    {
        enne=enne+2;
        k1=h*fn(enne,pl);
        k2=h*fn(enne+1,pl+k1/2.);
        k3=h*fn(enne+1,pl+k2/2.);
        k4=h*fn(enne+2,pl+k3);
        pl=pl+((k1+2*k2+2*k3+k4)/6.);
    }
    x=K*rmax;
    j=gsl_sf_bessel_jl(l,x);
    jj=gsl_sf_bessel_jl(l+1,x);
    n=gsl_sf_bessel_yl(l,x);
    nn=gsl_sf_bessel_yl(l+1,x);
    num=K*jj-j*((W+1)*tan(pl)+((1+l+spin)/rmax));
    den=K*nn-n*((W+1)*tan(pl)+((1+l+spin)/rmax));
    tandelta=num/den;
    return tandelta;
}

int main ()
{
    
//Input: Atomic number, electron energy (in eV), electron/positron, atomic polarization, initial spin polarization vector, azimuth
    
    cid1=fopen("INPUT/input.txt","r");
    fscanf(cid1,"%s %lg\n %s %lg\n %s %ld\n %s %lg\n %s %lg\n %s %lg\n %s %lg\n",in1,&Z,in2,&E,in3,&electron,in4,&Px,in5,&Py,in6,&Pz,in7,&azimuth);
    if(Z<1)
    {
        Z=1;
    };
    if(Z>92)
    {
        Z=92;
    }
    fclose(cid1);
    azimuth=azimuth*pi/180.;
    
//Input: Screening function (Cox and Bonham, 1967 - Salvat et al., 1987)
    
	cid2=fopen("INPUT/screening.txt","r");
	while(zz!=Z)
	{
	fscanf(cid2,"%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg",
	       &zz,&A1,&A2,&A3,&A4,&A5,&A6,&A7,&A8,&A9,&A10,&alp1,&alp2,&alp3,&alp4,&alp5,&alp6,&alp7,&alp8,&alp9,&alp10);
	}
	fclose(cid2);
    
    printf("\nZ = %lg\ngamma_i = %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg\nlambda_i = %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg\n",
           zz,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,alp1,alp2,alp3,alp4,alp5,alp6,alp7,alp8,alp9,alp10);
    
//Units
    
    alp1=alp1*hbardbymc/azero;
	alp2=alp2*hbardbymc/azero;
	alp3=alp3*hbardbymc/azero;
	alp4=alp4*hbardbymc/azero;
	alp5=alp5*hbardbymc/azero;
	alp6=alp6*hbardbymc/azero;
    alp7=alp7*hbardbymc/azero;
    alp8=alp8*hbardbymc/azero;
    alp9=alp9*hbardbymc/azero;
    alp10=alp10*hbardbymc/azero;
    
//Initial spin polarization
    
    P=sqrt(Px*Px+Py*Py+Pz*Pz);
    if(P>1.)
    {
        printf("\n|P| must be < 1\n\n");
        exit(0);
    }
    
//Scalar product of the initial polarization by the versor normal to the plane of scattering
    
    Pn=-Px*sin(azimuth)+Py*cos(azimuth);

//Calculation of the potential energy close to the centre of the atom
    
	Z0=Z*ebye/(hbardbymc*mcbyc);
	Z1=-Z0*(A1*alp1+A2*alp2+A3*alp3+A4*alp4+A5*alp5+A6*alp6+A7*alp7+A8*alp8+A9*alp9+A10*alp10);
	Z2=Z0*(A1*alp1*alp1+A2*alp2*alp2+A3*alp3*alp3+A4*alp4*alp4+A5*alp5*alp5
	       +A6*alp6*alp6+A7*alp7*alp7+A8*alp8*alp8+A9*alp9*alp9+A10*alp10*alp10)/2.;
	Z3=-Z0*(A1*alp1*alp1*alp1+A2*alp2*alp2*alp2+A3*alp3*alp3*alp3+A4*alp4*alp4*alp4+A5*alp5*alp5*alp5+A6*alp6*alp6*alp6+A7*alp7*alp7*alp7+A8*alp8*alp8*alp8+A9*alp9*alp9*alp9+A10*alp10*alp10*alp10)/6.;

//Radial atomic density file name
    
    sprintf(density,"RESULTS/AD-Z=%ld.txt",(long int)Z);
            
//Differential elastic scattering cross-section file name
    
    sprintf(descs,"RESULTS/D-Z=%ld-E=%ldeV.txt",(long int)Z,(long int)E);
    
//Sherman function file name
    
    sprintf(Sfunction,"RESULTS/S-Z=%ld-E=%ldeV.txt",(long int)Z,(long int)E);
    
//T function file name
    
    sprintf(Tfunction,"RESULTS/T-Z=%ld-E=%ldeV.txt",(long int)Z,(long int)E);
    
//U function file name
    
    sprintf(Ufunction,"RESULTS/U-Z=%ld-E=%ldeV.txt",(long int)Z,(long int)E);
    
//Differential elastic scattering cross-section (polarized beams) file name
    
    sprintf(descsp,"RESULTS/Dp-Z=%ld-E=%ldeV.txt",(long int)Z,(long int)E);
    
//Total and transport cross-sections file name
    
    sprintf(total,"RESULTS/I-Z=%ld-E=%ldeV.txt",(long int)Z,(long int)E);
    
//Cumulative probability file name
    
    sprintf(probability,"RESULTS/P-Z=%ld-E=%ldeV.txt",(long int)Z,(long int)E);
    
//Energy, momentum, rmax
    
    Er=E/mcbyc;
    W=Er+1;
    K=sqrt(W*W-1);
	rmax=5.;
    rmax=rmax/hbardbymc;
    
//Calculation of the potential V[enne] as a function of r[enne]
    
	rd=h;
    
	while (rd<rmax)
	{
        V[enne]=vi(rd);
        r[enne]=rd;
		rd=rd+h/2.;
        enne=enne+1;
	}
    nmax=enne;
    
//Calculation of the phase shifts T1[l],T2[l]
    
    while(l<2.*K*rmax)
    {
        spin=-(l+1);
        p0=0.5*asin(-Z0/spin);
        p1=(W+Z1-cos(2.*p0))/(1.-2.*spin*cos(2.*p0));
        p2=((1.-spin*p1)*2.*p1*sin(2.*p0)+Z2)/(2.-2.*spin*cos(2.*p0));
        p3=(1.-2.*spin*p1)*2.*p2*sin(2.*p0);
        p3=p3+(1.-(2.*spin*p1/3.))*2.*p1*p1*cos(2.*p0)+Z3;
        p3=p3/(3.-2.*spin*cos(2.*p0));
        pl=p0+p1*r[1]+p2*r[1]*r[1]+p3*r[1]*r[1]*r[1];
        td1=phaseshift();
        spin=l;
        if (spin!=0)
        {
            p0=0.5*(pi-asin(-Z0/spin));
            p1=(W+Z1-cos(2.*p0))/(1.-2.*spin*cos(2.*p0));
            p2=((1.-spin*p1)*2.*p1*sin(2.*p0)+Z2);
            p2=p2/(2.-2.*spin*cos(2.*p0));
            p3=(1.-2.*spin*p1)*2.*p2*sin(2.*p0);
            p3=p3+(1.-(2.*spin*p1/3.))*2.*p1*p1*cos(2.*p0)+Z3;
            p3=p3/(3.-2.*spin*cos(2.*p0));
            pl=p0+p1*r[1]+p2*r[1]*r[1]+p3*r[1]*r[1]*r[1];
            td2=phaseshift();
        }
        if (spin==0) td2=1.;
        T1[l]=td1;
        T2[l]=td2;
        l=l+1;
    }
    
    codSD=fopen(descs,"w");
    codSf=fopen(Sfunction,"w");
    codTf=fopen(Tfunction,"w");
    codUf=fopen(Ufunction,"w");
    codST=fopen(total,"w");
    codSDp=fopen(descsp,"w");
    codP=fopen(probability,"w");
    codAD=fopen(density,"w");

//Calculation of the scattering amplitudes f=sum1+isum2 and g=sum3+isum4 as a function of the scattering angle
    
    emme=1;
    thetamin=0.;
    thetamax=180.;
    dtheta=0.1;
    thetamin=thetamin*pi/180.;
    thetamax=thetamax*pi/180.;
    dtheta=dtheta*pi/180.;
    theta=thetamin;
    while (theta<thetamax-dtheta)
    {
        l=0;
        sum1=0.;
        sum2=0.;
        sum3=0.;
        sum4=0.;
        theta=theta+dtheta;
        ct=cos(theta);
        while(l<2.*K*rmax)
        {
            cd1=(1.-T1[l]*T1[l])/(1.+T1[l]*T1[l]);
            sd1=2.*T1[l]/(1.+T1[l]*T1[l]);
            cd2=(1.-T2[l]*T2[l])/(1.+T2[l]*T2[l]);
            sd2=2.*T2[l]/(1.+T2[l]*T2[l]);
            pol=gsl_sf_legendre_Pl(l,ct);
            if (l==0)
            {
                pol1=0;
            }
            else
            {
                pol1=gsl_sf_legendre_Plm(l,1,ct);
            }
            sum1=sum1+((l+1)*sd1+l*sd2)*pol;
            sum2=sum2+((2*l+1)-(l+1)*cd1-l*cd2)*pol;
            sum3=sum3+(sd1-sd2)*pol1;
            sum4=sum4+(cd2-cd1)*pol1;
            l=l+1;
        }

//Calculation of the differential elastic scattering cross-section as a function of the scattering angle (unpolarized beams)
        
        sigmad=(sum1*sum1+sum2*sum2+sum3*sum3+sum4*sum4)/(4.*K*K/(hbardbymc*hbardbymc));
        fprintf(codSD,"%lg %lg \n",theta*180./pi,sigmad);
        
//Calculation of the Sherman function as a function of the scattering angle
        
        Sf=2.*(sum1*sum4-sum3*sum2)/(sum1*sum1+sum2*sum2+sum3*sum3+sum4*sum4);
        fprintf(codSf,"%lg %lg \n",theta*180./pi,Sf);

//Calculation of the differential elastic scattering cross-section as a function of the scattering angle (polarized beams)
        
        sigmadp=(1.+Pn*Sf)*sigmad;
        fprintf(codSDp,"%lg %lg \n",theta*180./pi,sigmadp);
        
//Calculation of the T function as a function of the scattering angle
        
        Tf=(sum1*sum1+sum2*sum2-sum3*sum3-sum4*sum4)/(sum1*sum1+sum2*sum2+sum3*sum3+sum4*sum4);
        fprintf(codTf,"%lg %lg \n",theta*180./pi,Tf);
        
//Calculation of the U function as a function of the scattering angle
        
        Uf=2.*(sum1*sum3+sum2*sum4)/(sum1*sum1+sum2*sum2+sum3*sum3+sum4*sum4);
        fprintf(codUf,"%lg %lg \n",theta*180./pi,Uf);
        
// TEST: printf("%lg %lg \n",theta*180./pi,Sf*Sf+Tf*Tf+Uf*Uf);
        
        th[emme]=theta;
        sd[emme]=sigmad;
        emme=emme+1;
    }

// Calculation of the total, first transport and second transport elastic scattering cross-sections (by using Bode quadrature rule)
    
    mmax=emme;
    dth=th[2]-th[1];
    st=0.;
    stra=0.;
    stra2=0.;
    emme=1;
    while (emme<=mmax-4)
    {
        st=st+7.*sin(th[emme])*sd[emme];
        st=st+32.*sin(th[emme+1])*sd[emme+1];
        st=st+12.*sin(th[emme+2])*sd[emme+2];
        st=st+32.*sin(th[emme+3])*sd[emme+3];
        st=st+7.*sin(th[emme+4])*sd[emme+4];
        stra=stra+7.*(1-cos(th[emme]))*sin(th[emme])*sd[emme];
        stra=stra+32.*(1-cos(th[emme+1]))*sin(th[emme+1])*sd[emme+1];
        stra=stra+12.*(1-cos(th[emme+2]))*sin(th[emme+2])*sd[emme+2];
        stra=stra+32.*(1-cos(th[emme+3]))*sin(th[emme+3])*sd[emme+3];
        stra=stra+7.*(1-cos(th[emme+4]))*sin(th[emme+4])*sd[emme+4];
        stra2=stra2+7.*(1-(cos(th[emme]))*(cos(th[emme])))*sin(th[emme])*sd[emme];
        stra2=stra2+32.*(1-(cos(th[emme+1]))*(cos(th[emme+1])))*sin(th[emme+1])*sd[emme+1];
        stra2=stra2+12.*(1-(cos(th[emme+2]))*(cos(th[emme+2])))*sin(th[emme+2])*sd[emme+2];
        stra2=stra2+32.*(1-(cos(th[emme+3]))*(cos(th[emme+3])))*sin(th[emme+3])*sd[emme+3];
        stra2=stra2+7.*(1-(cos(th[emme+4]))*(cos(th[emme+4])))*sin(th[emme+4])*sd[emme+4];
        emme=emme+4;
    }
    
// Calculation of the cumulative probability (by using Bode quadrature rule)
    
    pel=0;
    emme=1;
    jei=0;
    for (mtheta=1;mtheta<=mmax;++mtheta)
    {
        while (emme<=mtheta-4)
        {
            pel=pel+7*sin(th[emme])*sd[emme];
            pel=pel+32*sin(th[emme+1])*sd[emme+1];
            pel=pel+12*sin(th[emme+2])*sd[emme+2];
            pel=pel+32*sin(th[emme+3])*sd[emme+3];
            pel=pel+7*sin(th[emme+4])*sd[emme+4];
            emme=emme+4;
        }
        jei=jei+1;
        if (jei==4)
        {
            fprintf(codP,"%lg %lg \n",mtheta*dth*180/pi,pel/st);
            jei=0;
        }
    }
    
//Calculation of the radial atomic density
    
    rd=h;
    while (rd<rmax)
    {
        fprintf(codAD,"%lg %lg \n",rd*hbardbymc,ad(rd)); // radius in A; density in A^-1
        rd=rd+h/2.;
    }
    
    fprintf(codST,"%s = %lg\n%s = %lg\n%s = %ld\n%s = %lg\n\n",in1,Z,in2,E,in3,electron,in7,azimuth*180./pi);
    fprintf(codST,"%s = %lg\n%s = %lg\n%s = %lg\n\n",in4,Px,in5,Py,in6,Pz);
    fprintf(codST,"P = %lg\n\n",P);
    fprintf(codST,"Pn = %lg \n\n",Pn);
    fprintf(codST,"sT(A*A) = %.4f \n",4.*pi*st*dth/45.);
    fprintf(codST,"sTRA(A*A) = %.4f \n",4.*pi*stra*dth/45.);
    fprintf(codST,"sTRA2(A*A) = %.4f \n\n",4.*pi*stra2*dth/30.);
    
    fclose(codSD);
    fclose(codSf);
    fclose(codTf);
    fclose(codUf);
    fclose(codST);
    fclose(codSDp);
    fclose(codP);
    fclose(codAD);
    
    printf("\n%s = %lg\n%s = %lg\n%s = %ld\n%s = %lg\n\n",in1,Z,in2,E,in3,electron,in7,azimuth*180./pi);
    printf("%s = %lg\n%s = %lg\n%s = %lg\n\n",in4,Px,in5,Py,in6,Pz);
    printf("P = %lg\n\n",P);
    printf("Pn = %lg \n\n",Pn);
    printf("sT(A*A) = %.4f \n",4.*pi*st*dth/45.);
    printf("sTRA(A*A) = %.4f \n",4.*pi*stra*dth/45.);
    printf("sTRA2(A*A) = %.4f \n\n",4.*pi*stra2*dth/30.);
    
//Gnuplot figures
    
    FILE *gnuplotAD = fopen("FIGURES/commandsAD.txt","w");
    fprintf(gnuplotAD, "set terminal pdf\n");
    fprintf(gnuplotAD, "set output 'FIGURES/AD-Z=%ld.pdf' \n",(long int)Z);
    fprintf(gnuplotAD, "set title 'Z=%ld' \n",(long int)Z);
    fprintf(gnuplotAD, "unset key \n");
    fprintf(gnuplotAD, "set termoption enhanced\n");
    fprintf(gnuplotAD, "set xlabel 'r (10^{-8}cm)' \n");
    fprintf(gnuplotAD, "set ylabel '4 {/Symbol p} r^2 {/Symbol r} (10^{8}cm^{-1})' \n");
    fprintf(gnuplotAD, "set xrange [0:2] \n");
    fprintf(gnuplotAD, "plot 'RESULTS/AD-Z=%ld.txt' w l \n",(long int)Z);
    fflush(gnuplotAD);
    fclose(gnuplotAD);
    system("gnuplot 'FIGURES/commandsAD.txt'");
    
    FILE *gnuplot = fopen("FIGURES/commandsD.txt","w");
    fprintf(gnuplot, "set terminal pdf\n");
    fprintf(gnuplot, "set output 'FIGURES/D-Z=%ld-E=%ldeV.pdf' \n",(long int)Z,(long int)E);
    fprintf(gnuplot, "set logscale y \n");
    fprintf(gnuplot, "set title 'Z=%ld, E=%ldeV' \n",(long int)Z,(long int)E);
    fprintf(gnuplot, "unset key \n");
    fprintf(gnuplot, "set xlabel 'Scattering Angle (deg)' \n");
    fprintf(gnuplot, "set ylabel 'DESCS (10^{-16}cm^2/sr)' \n");
    fprintf(gnuplot, "plot 'RESULTS/D-Z=%ld-E=%ldeV.txt' w l \n",(long int)Z,(long int)E);
    fflush(gnuplot);
    fclose(gnuplot);
    system("gnuplot 'FIGURES/commandsD.txt'");
    
    FILE *gnuplotS = fopen("FIGURES/commandsS.txt","w");
    fprintf(gnuplotS, "set terminal pdf\n");
    fprintf(gnuplotS, "set output 'FIGURES/S-Z=%ld-E=%ldeV.pdf' \n",(long int)Z,(long int)E);
    fprintf(gnuplotS, "set title 'Z=%ld, E=%ldeV' \n",(long int)Z,(long int)E);
    fprintf(gnuplotS, "unset key \n");
    fprintf(gnuplotS, "set xlabel 'Scattering Angle (deg)' \n");
    fprintf(gnuplotS, "set ylabel 'S' \n");
    fprintf(gnuplotS, "plot 'RESULTS/S-Z=%ld-E=%ldeV.txt' w l \n",(long int)Z,(long int)E);
    fflush(gnuplotS);
    fclose(gnuplotS);
    system("gnuplot 'FIGURES/commandsS.txt'");
    
    FILE *gnuplotT = fopen("FIGURES/commandsT.txt","w");
    fprintf(gnuplotT, "set terminal pdf\n");
    fprintf(gnuplotT, "set output 'FIGURES/T-Z=%ld-E=%ldeV.pdf' \n",(long int)Z,(long int)E);
    fprintf(gnuplotT, "set title 'Z=%ld, E=%ldeV' \n",(long int)Z,(long int)E);
    fprintf(gnuplotT, "unset key \n");
    fprintf(gnuplotT, "set xlabel 'Scattering Angle (deg)' \n");
    fprintf(gnuplotT, "set ylabel 'T' \n");
    fprintf(gnuplotT, "plot 'RESULTS/T-Z=%ld-E=%ldeV.txt' w l \n",(long int)Z,(long int)E);
    fflush(gnuplotT);
    fclose(gnuplotT);
    system("gnuplot 'FIGURES/commandsT.txt'");
    
    FILE *gnuplotU = fopen("FIGURES/commandsU.txt","w");
    fprintf(gnuplotU, "set terminal pdf\n");
    fprintf(gnuplotU, "set output 'FIGURES/U-Z=%ld-E=%ldeV.pdf' \n",(long int)Z,(long int)E);
    fprintf(gnuplotU, "set title 'Z=%ld, E=%ldeV' \n",(long int)Z,(long int)E);
    fprintf(gnuplotU, "unset key \n");
    fprintf(gnuplotU, "set xlabel 'Scattering Angle (deg)' \n");
    fprintf(gnuplotU, "set ylabel 'U' \n");
    fprintf(gnuplotU, "plot 'RESULTS/U-Z=%ld-E=%ldeV.txt' w l \n",(long int)Z,(long int)E);
    fflush(gnuplotU);
    fclose(gnuplotU);
    system("gnuplot 'FIGURES/commandsU.txt'");
    
    FILE *gnuplotDp = fopen("FIGURES/commandsDp.txt","w");
    fprintf(gnuplotDp, "set terminal pdf\n");
    fprintf(gnuplotDp, "set output 'FIGURES/Dp-Z=%ld-E=%ldeV.pdf' \n",(long int)Z,(long int)E);
    fprintf(gnuplotDp, "set logscale y \n");
    fprintf(gnuplotDp, "set title 'Z=%ld, E=%ldeV' \n",(long int)Z,(long int)E);
    fprintf(gnuplotDp, "unset key \n");
    fprintf(gnuplotDp, "set xlabel 'Scattering Angle (deg)' \n");
    fprintf(gnuplotDp, "set ylabel 'DESCS (10^{-16}cm^2/sr)' \n");
    fprintf(gnuplotDp, "plot 'RESULTS/Dp-Z=%ld-E=%ldeV.txt' w l \n",(long int)Z,(long int)E);
    fflush(gnuplotDp);
    fclose(gnuplotDp);
    system("gnuplot 'FIGURES/commandsDp.txt'");
    
    FILE *gnuplotCOMP = fopen("FIGURES/commandsCOMP.txt","w");
    fprintf(gnuplotCOMP, "set terminal pdf\n");
    fprintf(gnuplotCOMP, "set output 'FIGURES/COMP-Z=%ld-E=%ldeV.pdf' \n",(long int)Z,(long int)E);
    fprintf(gnuplotCOMP, "set logscale y \n");
    fprintf(gnuplotCOMP, "set title 'Z=%ld, E=%ldeV' \n",(long int)Z,(long int)E);
    fprintf(gnuplotCOMP, "set xlabel 'Scattering Angle (deg)' \n");
    fprintf(gnuplotCOMP, "set ylabel 'DESCS (10^{-16}cm^2/sr)' \n");
    fprintf(gnuplotCOMP, "plot 'RESULTS/D-Z=%ld-E=%ldeV.txt' w l title 'unpolarized beam','RESULTS/Dp-Z=%ld-E=%ldeV.txt' w l title 'polarized beam' \n",(long int)Z,(long int)E,(long int)Z,(long int)E);
    fflush(gnuplotCOMP);
    fclose(gnuplotCOMP);
    system("gnuplot 'FIGURES/commandsCOMP.txt'");
    
    FILE *gnuplotP = fopen("FIGURES/commandsP.txt","w");
    fprintf(gnuplotP, "set terminal pdf\n");
    fprintf(gnuplotP, "set output 'FIGURES/P-Z=%ld-E=%ldeV.pdf' \n",(long int)Z,(long int)E);
    fprintf(gnuplotP, "set title 'Z=%ld, E=%ldeV' \n",(long int)Z,(long int)E);
    fprintf(gnuplotP, "unset key \n");
    fprintf(gnuplotP, "set xlabel 'Scattering Angle (deg)' \n");
    fprintf(gnuplotP, "set ylabel 'Cumulative probability' \n");
    fprintf(gnuplotP, "plot 'RESULTS/P-Z=%ld-E=%ldeV.txt' w l \n",(long int)Z,(long int)E);
    fflush(gnuplotP);
    fclose(gnuplotP);
    system("gnuplot 'FIGURES/commandsP.txt'");
}
