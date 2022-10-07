#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;
const double k=1.9872036*1e-3; //in kcal/(K¡Pmol)

class MODEL
{
   public:
   MODEL();              //class constructor
   ~MODEL();            //class destructor

   //model description
   int npart;               //number of particles
   double (*x)[3];      //position of particles (in Angstroms)
   double (*m);          //mass (in g/mol)
   double cell[3];      //unit cell length (in A)
   bool   period;       //periodicity 0:nonperiodic 1:periodic
   bool   widom; 		//flag for Widom insertion test


   //MC parameters
   double Tset;         //temperature used in velocity initialization
   int ncycl;              //total steps to run  
   int icycl;               //current step 
   double Rcut;         //cutoff distance for vdw interaction
   double rc2;            //Rcut*Rcut;
   double delx;           //maxmum moving distance
   double beta;           //1/kT
   int accpt;             //number of accepted move
   int Mtest; 			  //number of insertion tests
   double wtest;          //variable for Widom insertion test
   double muex;           //"instantanous" excess chemical potential in kcal/mol
   
   //force field parameters
   double Do;  //in kcal/mol
   double Ro;  //in Angstrom
   double p0;  //p0= 12*Do*pow(Ro,12) in unit kcal/mol*A12
   double p1;  //p1= 12*Do*pow(Ro,6) in unit kcal/mol*A6
   double Ro_Rc3;
   double Ro_Rc6;
   double Ro_Rc9;
   double Ro3;
   double unitc;
    
   //model properties
   double T;                 //temperature (in K)
   int df;                       //degree of freedom
   double Ek;               //kinetic energy in kcal/mol
   double Ep;               //potential energy in kcal/mol
   double Etot;             //Ek+Ep
   double P;                 //pressure    (in GPa)
   double V;                 //volume      (in A3)
   double Ptail;             //tail correction for pressure (GPa)
   double Eptail;           //tail correction for P.E. (kcal/mol)
   double strs[6]; //stress tensor in kcal/mol-A3

   //class functions
   int init();                    //initialization of variables
   int mcmov();              //attempt for Monte Carlo move
   int ener(int ,double &);  //evaluation of interaction energy with one particle
   int tote(double &);               //evaluate Ek, Ep, and Etot, stress
   int sample();             //calculation of properties
   int ener(double*,double &); //calc of interaction energy of test particle
   int Widom(); 					//perform Widom insertion
   double myrand01(); //generate a random number between 0 and 1 [0,1]
   int dump_trj();           //output trajectory data to a file   
   ofstream outf;           //file stream of the trajectory
   ofstream logf;           //file stream of the log
   
   //FEP
   int i;
   int step;
   double lambda_ini;
   double lambda_final;
   double lambda_1;
   double lambda_2;
   double delA;
   int FEP();
   int FEP_ini();
   ofstream FEP_result;           //file stream of the log
   
   //average
   int avgcycl;		//number of cycles to average
   double Tavg;
   double Pavg;
   double Vavg;
   double Ekavg;
   double Epavg;
   double Etotavg;
   double delA_tot;
   int count;
   int avgstep;
   int avg();
};

int main()
{
    MODEL model;                                 //declare a system
    char null[1024];
    int i;
    double temp;
	
	model.outf.open("MC_NVT.gro",ios::out);   //output filename for trajectory
    model.logf.open("MC_NVT.log",ios::out);   //output filename for trajectory
  	model.FEP_result.open("MC_NVT_FEP.log",ios::out);   //output filename for trajectory
  	
  	model.step=300;				//number of steps from A to B
	model.lambda_ini=0.0006;					//state A
    model.lambda_final=1;					//state B
	
	for(i=0;i<2;i++){
		
		model.delA_tot=0;
		
		sprintf(null,"Start to conduct from lambda %f ~ %f divided by %d steps¡G",model.lambda_ini,model.lambda_final,model.step);
    	cout<<null<<endl<<endl;
    	model.logf<<null<<endl<<endl;
		model.FEP_result<<null<<endl;
		
		for(model.i=0;model.i<model.step;model.i++){
			model.init();  								 //initialization
			model.FEP_ini();
			
			if(model.lambda_1==0) {
    			sprintf(null,"From lambda %f ~ %f¡Gskip",model.lambda_1,model.lambda_2);
    			cout<<null<<endl<<endl;
    			model.logf<<null<<endl<<endl;
				model.FEP_result<<null<<endl;
				continue;
			}

    		model.sample();                               //initial property 
  		
  		  	//MC loop
  		  	while(model.icycl<model.ncycl) {
    	    	model.mcmov();                           //attempts to displace a particle
			
				//if( model.icycl%(model.npart*10)==0 ) model.sample();                            //determine system property
				if( model.icycl%100000==0 ) model.sample();                            //determine system property
    	    	//model.sample(); 																   //determine system property 
    	    
    	    	if((model.icycl>model.ncycl-model.avgcycl)&&(model.icycl%model.avgstep==0)) model.FEP();
			}
			model.avg();
		}
		temp=model.lambda_ini;
		model.lambda_ini=model.lambda_final;
		model.lambda_final=temp;
	}		
    	return 0;
}
    	
int MODEL::FEP_ini(){
	char null[1024];
	
	lambda_1=lambda_ini+i*((lambda_final-lambda_ini)/step);			//state 1
    lambda_2=lambda_ini+(i+1)*((lambda_final-lambda_ini)/step);			//state 2
    delA=0;
    
    avgcycl=ncycl/10;		//number of cycles to average
	Tavg=0;
	Pavg=0;
	Vavg=0;
	Ekavg=0;
	Epavg=0;
	Etotavg=0;
	count=0;
	avgstep=1000;
	
	if(lambda_1!=0) {
		sprintf(null,"From lambda %f ~ %f¡G",lambda_1,lambda_2);
		cout<<null<<endl;
    	logf<<null<<endl;
    	FEP_result<<null;
	}
	
	return 0;	
}
			
int MODEL::init()
{   
    //simulation parameters     

    icycl=0;                                  //current step
    npart=512;                              //number of particles
    Tset=100;                                //target temperature in K
    //ncycl=npart*100000;                        //steps to run
    ncycl=500000;                        //steps to run
    widom=0; 			//flag for Widom insertion test
    Mtest=100000; 		//number of insertion tests
    
    T=Tset;

    //allocation of memerory    
    x=new double [npart][3];        //position in Angstrom
    m=new double [npart];             //mass in g/mol

    //assign mass of particles
    int i,j,k;
    for(i=0;i<npart;i++) m[i]=39.948; //molecular weight of Argon        
    //One simple way to place particles in space
    double sep=4;  //separation distance between two particles
    int nt=0;

    //give force field parameters
    Do= 0.185;    // in kcal/mol
    Ro= 3.868;    // in Angstrom
    p0= 12*Do*pow(Ro,12); //in unit kcal/mol*A12
    p1= 12*Do*pow(Ro,6);   //in unit kcal/mol*A6
    cell[0]=cell[1]=cell[2]=pow(30720,1.0/3.0);      //length of unit cell in Angstroms
    period=1;                                //flag for periodicity
    Rcut=4.5*Ro;                          //cutoff distance for vdw interaction
    rc2=Rcut*Rcut;                       //Rcut2 
    
    delx=0.89; // angstrom
    accpt=0;
    beta= 4184/(8.314*T); // in unit of 1/(kcal/mol)
    
    if(period) { //periodic system
      //A rough method to place particles in the box
      double delta[3];
      int pps;
      pps= int(pow(npart,1.0/3.0))+1;            //particles per side
      for(k=0;k<3;k++) delta[k]=cell[k]/pps; //spacing of particles
      for(i=0;i<pps;i++) {
          for(j=0;j<pps;j++) {
              for(k=0;k<pps;k++) {
                  x[nt][0]=i*delta[0];
                  x[nt][1]=j*delta[1];
                  x[nt][2]=k*delta[2];
                  nt++;                                       //number of particles placed
                  if(nt==npart) i=j=k=pps;        //nt has reached npart, exit loops
              }
          }
      }    
    } else { // previous code for non-periodic system
      for(i=0;i<npart;i++) {
          for(j=0;j<=i;j++) {
              for(k=0;k<=j;k++) {
                  x[nt][0]=i*sep;
                  x[nt][1]=j*sep;
                  x[nt][2]=k*sep;
                  nt++;
                  if(nt==npart) i=j=k=npart;  //stop when nt==npart
              }
          }
      }
    }

    df= 3*npart;                                    //degree of freedom

   //calculate volume
    V=cell[0]*cell[1]*cell[2];   //in A^3
    
    //Tail correction for Ep
    Ro_Rc3=pow(Ro/Rcut,3.0);
    Ro_Rc6=Ro_Rc3*Ro_Rc3;
    Ro_Rc9=Ro_Rc6*Ro_Rc3;
    Ro3=Ro*Ro*Ro;
    Eptail=lambda_1*((4.0/3.0)*3.1415926*npart*npart/V*Do*Ro3*(Ro_Rc9/6.0-Ro_Rc3));

    //Tail correction for pressure
    unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
    Ptail=lambda_1*((8.0/3.0)*3.1415926*npart*npart/V/V*Do*Ro3*(Ro_Rc9/3.0-Ro_Rc3)*unitc); 
    
    return 0;
}

int MODEL::mcmov()
{  //attempts to displace a particle
    int k,o,n;
    double eno,enn,xo[3];

    o= int((rand()/(RAND_MAX+1.0))*npart);  //a random # btwn 0 and npart-1
    ener(o,eno);  //calculates interaction energy of o with other particles

    for(k=0;k<3;k++) {
        xo[k]=x[o][k]; //store the original coordinates
        x[o][k] = x[o][k] + (myrand01()-0.5)*delx; //give particle random movement
    }
    ener(o,enn);  //interaction energy of particle o with other particles
    
    if( myrand01() > exp(-beta*(enn-eno)) ) {   //reject the move
        for(k=0;k<3;k++) x[o][k] = xo[k];
    } else accpt++;                 //accpt records the number of successful attempts
    icycl++;                             //icycl records the total number of attempts
    
    return 0;
} 

int MODEL::ener(int i,double &energy)
{   //this is almost the same as the energy calculation in force()
    int j,k,nbox;
    double r2,r2i,r6i,xr[3],redu;    
    energy=0;
    
        for(j=0;j<npart;j++) {
            if(j==i) continue;
            for(k=0;k<3;k++) xr[k]= x[i][k] - x[j][k];  //distance vector
            
            if(period) { //periodic system, find dist within one cell  
                 for(k=0;k<3;k++) {
                     redu= (xr[k]/cell[k]);              //reduced coordinates
                     redu= redu - round (redu);   //between -0.5 and 0.5
                     xr[k] = redu*cell[k];               //real coordinates 
                 } 
            }                         
            r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];  //distance square
            if(r2<rc2) {  //within cutoff distance
               r2i= 1/r2;
               r6i= r2i*r2i*r2i;     
               energy += lambda_1*((p0*r6i - p1*2)*r6i/12.0);
            }
        }
       
    return 0;
}

int MODEL::tote(double &lambda)
{   //this is almost the same as the energy calculation in force()
    int i,j,k,nbox;
    double r2,r2i,r6i,xr[3],redu,ff;    
    Ep=0;
    
    for(i=0;i<6;i++) strs[i]=0;       //set stress to zero 
        
    for(i=0;i<npart;i++) {
        for(j=i+1;j<npart;j++) {

            for(k=0;k<3;k++) xr[k]= x[i][k] - x[j][k];  //distance vector
            
            if(period) { //periodic system, find dist within one cell  
                 for(k=0;k<3;k++) {
                     redu= (xr[k]/cell[k]);              //reduced coordinates
                     redu= redu - round (redu);   //between -0.5 and 0.5
                     xr[k] = redu*cell[k];               //real coordinates 
                 } 
            }                         
            r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];  //distance square
            if(r2<rc2 || period==0) {  //within cutoff distance
               r2i= 1/r2;
               r6i= r2i*r2i*r2i;
               ff = lambda*(r6i*(p0*r6i-p1)*r2i);          //in unit of kcal/mol/A2
               //the stress tensor 
               strs[0]+= ff*xr[0]*xr[0];  //xx in unit of kcal/mol
               strs[1]+= ff*xr[1]*xr[1];  //yy
               strs[2]+= ff*xr[2]*xr[2];  //zz
               strs[3]+= ff*xr[0]*xr[1];  //xy
               strs[4]+= ff*xr[0]*xr[2];  //xz
               strs[5]+= ff*xr[1]*xr[2];  //yz          
               Ep += lambda*((p0*r6i - p1*2)*r6i/12.0);  //in unit of kcal/mol
            }
        }
    }
    
    //Ek=df/2*RT
    Ek= 0.5*df*8.314*T/4184; //in kcal/mol    
    Eptail=lambda*((4.0/3.0)*3.1415926*npart*npart/V*Do*Ro3*(Ro_Rc9/6.0-Ro_Rc3));
    Ep += Eptail;    
    Etot=Ek+Ep;   
	
    //calcualate pressure
    Ptail=lambda*((8.0/3.0)*3.1415926*npart*npart/V/V*Do*Ro3*(Ro_Rc9/3.0-Ro_Rc3)*unitc);
    P=((Ek*2+strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
    P+= Ptail;	
	 
    return 0;
}

MODEL::MODEL()
{
    
};

MODEL::~MODEL()
{
    delete [] x;
    delete [] m;
    outf.close();
    logf.close();
    FEP_result.close();
};

double MODEL::myrand01()
{
    return rand()/double(RAND_MAX);  
     //returns a number between 0 (inclusive) and 1 (inclusive)
}


int MODEL::Widom()
{
	int i,j,k;
	double test[3],energy;
	//determine the position of test particle randomly inside the box
	for(k=0;k<3;k++) test[k]=myrand01()*cell[k];
	ener(test,energy); //interaction energy of test particle
	wtest += exp( -beta*energy);
	return 0;
} 

int MODEL::ener(double *test,double &energy)
{ //this function is similar to the ener() for interaction energy
	int j,k,nbox;
	double r2,r2i,r6i,xr[3];
	energy=0;
	for(j=0;j<npart;j++) {
		for(k=0;k<3;k++) xr[k]= test[k] - x[j][k]; //distance vector
		if(period) { //periodic system, find dist within one cell
			for(k=0;k<3;k++) {
				xr[k] += 0.5*cell[k]; //shift j to box center
				if(xr[k]<0) nbox = int( (xr[k]/cell[k]) -1); //determine image
				else nbox=int(xr[k]/cell[k]);
				xr[k] -= ( nbox+0.5)*cell[k]; //map to original box and shift
			}
		}
		r2= xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //distance square
		if(r2<rc2) { //within cutoff distance
			r2i= 1/r2;
			r6i= r2i*r2i*r2i;
			energy += (p0*r6i - p1*2)*r6i/12.0;
		}
	}
	return 0;
}

int MODEL::sample()
{
    char null[1024];	
    //calculation of system temperature and kinetic energy
    int i,j,k;
    
    tote(lambda_1);
    
    //perform Widom insertion test
	if(widom) {
		wtest=0;
		muex=0;
		for(i=0;i<Mtest;i++) {
			Widom();
		}
		muex = - log(wtest/Mtest)/beta;
		muex += 2*Eptail/npart; //tail correction
	}



    //calcualate pressure
//    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
//    P=((Ek*2+strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
//    P+= Ptail;

    //Display current information
    sprintf(null,"Lambda %f¡Gcurrent Step %d, accepted fraction %f, ",lambda_1,icycl,double(accpt)/icycl);
    cout<<null;
    logf<<null;
    sprintf(null,"T %f K, P %f MPa, V %f A3, ",T,P*1000,V/npart);
    cout<<null;
    logf<<null;    
    sprintf(null,"E(kcal/mol) Ek %f Ep %f Etot %f muex %f ",Ek/npart,Ep/npart,Etot/npart,muex/npart);
	cout<<null<<endl;
	logf<<null<<endl; 
    //printf("Tail contribution %.0f%% in Ep %.0f%% in P \n",Eptail/Ep*100,Ptail/P*100);


    dump_trj();                     //output trajectory file
    return 0;
}

int MODEL::dump_trj()
{
    char null[1024];
    int i;    
    sprintf(null,"Current Step %d, accepted fraction %f, ",icycl,double(accpt)/icycl);
    outf<<null<<endl; //comment in trj file
    sprintf(null,"%5d",npart);
    outf<<null<<endl;    //number of particles
    for(i=0;i<npart;i++) { //position (nm) and velocity (nm/ps)
        sprintf(null,"%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f",
        1,"Ar","Ar",i+1,x[i][0]/10,x[i][1]/10,x[i][2]/10,0,0,0);
        outf<<null<<endl;
    }
    sprintf(null,"%10.5f%10.5f%10.5f",cell[0]/10,cell[1]/10,cell[2]/10);
    outf<<null<<endl;  //box information
    return 0;
}

int MODEL::FEP(){
	double E1,E2;
	
	tote(lambda_2);
    E2=Etot;
	tote(lambda_1);
    E1=Etot;
    
    delA+=exp(-((E2-E1)/npart)/(k*T));
    
    Tavg+=T;
	Pavg+=P;
	Vavg+=V;
	Ekavg+=Ek;
	Epavg+=Ep;
	Etotavg+=Etot;
	
	count++;
	
	return 0;	
}

int MODEL::avg(){
	char null[1024];
	
	Tavg=Tavg/count;
	Pavg=Pavg/count;
	Vavg=Vavg/count;
	Ekavg=Ekavg/count;
	Epavg=Epavg/count;
	Etotavg=Etotavg/count;
			
	delA=-k*Tavg*log(delA/count);
	delA_tot+=delA;
	
	sprintf(null,"From lambda %f ~ %f(average from step %d ~ %d by every %d step)¡G",lambda_1,lambda_2,ncycl-avgcycl+1,ncycl,avgstep);
	cout<<endl<<null<<endl;
	logf<<endl<<null<<endl;
	sprintf(null,"T %f K, P %f MPa, V %f A3, ",Tavg,Pavg*1000,Vavg/npart);
	cout<<null;
	logf<<null;
    FEP_result<<null;
	sprintf(null,"E(kcal/mol) Ek %f Ep %f Etot %f ",Ekavg/npart,Epavg/npart,Etotavg/npart);
    cout<<null;
	logf<<null;
    FEP_result<<null;    
    sprintf(null,"delA = %f, delA_tot = %f",delA,delA_tot); 
    cout<<null<<endl<<endl;
	logf<<null<<endl<<endl;
    FEP_result<<null<<endl; 
    
    return 0;
}




