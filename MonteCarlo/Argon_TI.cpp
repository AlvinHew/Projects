#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

//double lamb[3] = { 0.1127, 0.5, 0.88729 };

double lamb[64] = {
0.00034748
,0.00182994
,0.00449331
,0.00833187
,0.01333659
,0.0194956
,0.02679431
,0.03521541
,0.04473893
,0.05534228
,0.0670003
,0.07968535
,0.09336734
,0.10801382
,0.12359005
,0.14005908
,0.15738184
,0.17551726
,0.19442232
,0.21405218
,0.23436027
,0.25529843
,0.27681699
,0.29886492
,0.32138992
,0.34433856
,0.36765642
,0.39128818
,0.41517779
,0.43926859
,0.46350344
,0.48782485
,0.51217515
,0.53649656
,0.56073141
,0.58482221
,0.60871182
,0.63234358
,0.65566144
,0.67861008
,0.70113508
,0.72318301
,0.74470157
,0.76563973
,0.78594782
,0.80557768
,0.82448274
,0.84261816
,0.85994092
,0.87640995
,0.89198618
,0.90663266
,0.92031465
,0.9329997
,0.94465772
,0.95526107
,0.96478459
,0.97320569
,0.9805044
,0.98666341
,0.99166813
,0.99550669
,0.99817006
,0.99965252
};


class MODEL
{
public:
    MODEL();              //class constructor
    ~MODEL();            //class destructor

    //model description
    int npart;               //number of particles
    double(*x)[3];      //position of particles (in Angstroms)
    double(*m);          //mass (in g/mol)
    double cell[3];      //unit cell length (in A)
    bool   period;       //periodicity 0:nonperiodic 1:periodic

    //MC parameters
    double Tset;         //temperature used in velocity initialization
    int ncycl;              //total steps to run  
    int icycl;               //current step 
    double Rcut;         //cutoff distance for vdw interaction
    double rc2;            //Rcut*Rcut;
    double delx;           //maxmum moving distance
    double beta;           //1/kT
    int accpt;             //number of accepted move

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
    double lambda;

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
    int ener(int, double&);  //evaluation of interaction energy with one particle
    int tote();               //evaluate Ek, Ep, and Etot, stress
    int sample();             //calculation of properties
    int sample2();             
    double myrand01(); //generate a random number between 0 and 1 [0,1]
    int dump_trj();           //output trajectory data to a file   
    ofstream outf;           //file stream of the trajectory
    ofstream logf;           //file stream of the log

    double num_sample;
    double cum_Ep;
};

int main()
{
    for (int i = 0; i < 64; i++) {
        
        //cout << sizeof(lamb) << endl;

        MODEL model;                                 //declare a system

        model.lambda = lamb[i];

        model.init();                                      //initialization
        //model.sample();                               //initial property 
        model.tote();

        model.num_sample = 0.0;
        model.cum_Ep = 0.0;

        //MC loop
        while (model.icycl < model.ncycl) {
            model.mcmov();                           //attempts to displace a particle
            if ((model.icycl > 500000) && (model.icycl % (model.npart * 10) == 0)) model.sample();                            //determine system property
            //model.num_sample++;
            //model.sample();
        }
        model.sample2();
    }

    return 0;
}

int MODEL::init()
{
    //simulation parameters     

    icycl = 0;                                  //current step
    npart = 216;                              //number of particles
    Tset = 100;                                //target temperature in K
    ncycl = npart * 2500;                        //steps to run
    //lambda = 0.99078;

    T = Tset;

    //allocation of memerory    
    x = new double[npart][3];        //position in Angstrom
    m = new double[npart];             //mass in g/mol

    //assign mass of particles
    int i, j, k;
    for (i = 0; i < npart; i++) m[i] = 39.948; //molecular weight of Argon        
    //One simple way to place particles in space
    double sep = 4;  //separation distance between two particles
    int nt = 0;

    //give force field parameters
    Do = 0.185;    // in kcal/mol
    Ro = 3.868;    // in Angstrom
    p0 = 12 * Do * pow(Ro, 12); //in unit kcal/mol*A12
    p1 = 12 * Do * pow(Ro, 6);   //in unit kcal/mol*A6
    cell[0] = cell[1] = cell[2] = pow(10612, 1.0 / 3.0);      //length of unit cell in Angstroms
    period = 1;                                //flag for periodicity
    Rcut = 4.5 * Ro;                          //cutoff distance for vdw interaction
    rc2 = Rcut * Rcut;                       //Rcut2 

    delx = 0.89; // angstrom
    accpt = 0;
    beta = 4184 / (8.314 * T); // in unit of 1/(kcal/mol)

    if (period) { //periodic system
      //A rough method to place particles in the box
        double delta[3];
        int pps;
        pps = int(pow(npart, 1.0 / 3.0)) + 1;            //particles per side
        for (k = 0; k < 3; k++) delta[k] = cell[k] / pps; //spacing of particles
        for (i = 0; i < pps; i++) {
            for (j = 0; j < pps; j++) {
                for (k = 0; k < pps; k++) {
                    x[nt][0] = i * delta[0];
                    x[nt][1] = j * delta[1];
                    x[nt][2] = k * delta[2];
                    nt++;                                       //number of particles placed
                    if (nt == npart) i = j = k = pps;        //nt has reached npart, exit loops
                }
            }
        }
    }
    else { // previous code for non-periodic system
        for (i = 0; i < npart; i++) {
            for (j = 0; j <= i; j++) {
                for (k = 0; k <= j; k++) {
                    x[nt][0] = i * sep;
                    x[nt][1] = j * sep;
                    x[nt][2] = k * sep;
                    nt++;
                    if (nt == npart) i = j = k = npart;  //stop when nt==npart
                }
            }
        }
    }

    df = 3 * npart;                                    //degree of freedom

   //calculate volume
    V = cell[0] * cell[1] * cell[2];   //in A^3

    //Tail correction for Ep
    Ro_Rc3 = pow(Ro / Rcut, 3.0);
    Ro_Rc6 = Ro_Rc3 * Ro_Rc3;
    Ro_Rc9 = Ro_Rc6 * Ro_Rc3;
    Ro3 = Ro * Ro * Ro;
    Eptail = lambda * (4.0 / 3.0) * 3.1415926 * npart * npart / V * Do * Ro3 * (Ro_Rc9 / 6.0 - Ro_Rc3);

    //Tail correction for pressure
    unitc = 4184 * 1E21 / 6.0221367E23; //conver from kcal/mol/A3 to GPa
    Ptail = lambda * (8.0 / 3.0) * 3.1415926 * npart * npart / V / V * Do * Ro3 * (Ro_Rc9 / 3.0 - Ro_Rc3) * unitc;

    //outf.open("ArMC_0.99078.gro", ios::out);   //output filename for trajectory
    //logf.open("ArMC_0.99078.log", ios::out);   //output filename for trajectory

    return 0;
}

int MODEL::mcmov()
{  //attempts to displace a particle
    int k, o, n;
    double eno, enn, xo[3];

    o = int((rand() / (RAND_MAX + 1.0)) * npart);  //a random # btwn 0 and npart-1
    ener(o, eno);  //calculates interaction energy of o with other particles

    for (k = 0; k < 3; k++) {
        xo[k] = x[o][k]; //store the original coordinates
        x[o][k] = x[o][k] + (myrand01() - 0.5) * delx; //give particle random movement
    }
    ener(o, enn);  //interaction energy of particle o with other particles

    if (myrand01() > exp(-beta * (enn - eno))) {   //reject the move
        for (k = 0; k < 3; k++) x[o][k] = xo[k];
    }
    else accpt++;                 //accpt records the number of successful attempts
    icycl++;                             //icycl records the total number of attempts

    return 0;
}

int MODEL::ener(int i, double& energy)
{   //this is almost the same as the energy calculation in force()
    int j, k, nbox;
    double r2, r2i, r6i, xr[3], redu;
    energy = 0;

    for (j = 0; j < npart; j++) {
        if (j == i) continue;
        for (k = 0; k < 3; k++) xr[k] = x[i][k] - x[j][k];  //distance vector

        if (period) { //periodic system, find dist within one cell  
            for (k = 0; k < 3; k++) {
                redu = (xr[k] / cell[k]);              //reduced coordinates
                redu = redu - round(redu);   //between -0.5 and 0.5
                xr[k] = redu * cell[k];               //real coordinates 
            }
        }
        r2 = xr[0] * xr[0] + xr[1] * xr[1] + xr[2] * xr[2];  //distance square
        if (r2 < rc2) {  //within cutoff distance
            r2i = 1 / r2;
            r6i = r2i * r2i * r2i;
            energy += lambda * (p0 * r6i - p1 * 2) * r6i / 12.0;
        }
    }

    return 0;
}

//for sampling and output info, doesnt participate in simulation
int MODEL::tote()
{   //this is almost the same as the energy calculation in force()
    int i, j, k, nbox;
    double r2, r2i, r6i, xr[3], redu, ff;
    Ep = 0;

    for (i = 0; i < 6; i++) strs[i] = 0;       //set stress to zero 

    for (i = 0; i < npart; i++) {
        for (j = i + 1; j < npart; j++) {

            for (k = 0; k < 3; k++) xr[k] = x[i][k] - x[j][k];  //distance vector

            if (period) { //periodic system, find dist within one cell  
                for (k = 0; k < 3; k++) {
                    redu = (xr[k] / cell[k]);              //reduced coordinates
                    redu = redu - round(redu);   //between -0.5 and 0.5
                    xr[k] = redu * cell[k];               //real coordinates 
                }
            }
            r2 = xr[0] * xr[0] + xr[1] * xr[1] + xr[2] * xr[2];  //distance square
            if (r2 < rc2 || period == 0) {  //within cutoff distance
                r2i = 1 / r2;
                r6i = r2i * r2i * r2i;
                /*
                ff = lambda * r6i * (p0 * r6i - p1) * r2i;          //in unit of kcal/mol/A2
                //the stress tensor 
                strs[0] += ff * xr[0] * xr[0];  //xx in unit of kcal/mol
                strs[1] += ff * xr[1] * xr[1];  //yy
                strs[2] += ff * xr[2] * xr[2];  //zz
                strs[3] += ff * xr[0] * xr[1];  //xy
                strs[4] += ff * xr[0] * xr[2];  //xz
                strs[5] += ff * xr[1] * xr[2];  //yz    
                */
                Ep += lambda * (p0 * r6i - p1 * 2) * r6i / 12.0;  //in unit of kcal/mol
            }
        }
    }
    
    //Ek=df/2*RT
    
    //Ek = 0.5 * df * 8.314 * T / 4184; //in kcal/mol    
    Eptail = lambda * (4.0 / 3.0) * 3.1415926 * npart * npart / V * Do * Ro3 * (Ro_Rc9 / 6.0 - Ro_Rc3);
    Ep += Eptail;
    /*
    Etot = Ek + Ep;

    //calcualate pressure
    Ptail = lambda * (8.0 / 3.0) * 3.1415926 * npart * npart / V / V * Do * Ro3 * (Ro_Rc9 / 3.0 - Ro_Rc3) * unitc;
    P = ((Ek * 2 + strs[0] + strs[1] + strs[2]) * unitc) / (3.0 * V); //in GPa
    P += Ptail;
    */
    return 0;
}

MODEL::MODEL()
{

};

MODEL::~MODEL()
{
    delete[] x;
    delete[] m;
    outf.close();
    logf.close();
};

double MODEL::myrand01()
{
    return rand() / double(RAND_MAX);
    //returns a number between 0 (inclusive) and 1 (inclusive)
}


int MODEL::sample()
{
    char null[1024];
    //calculation of system temperature and kinetic energy
    //int i, j, k;

    tote();

    //calcualate pressure
//    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
//    P=((Ek*2+strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
//    P+= Ptail;

    //Display current information
    /*
    sprintf_s(null, "Current Step %d, accepted fraction %f, ", icycl, double(accpt) / icycl);
    cout << null;
    logf << null;
    sprintf_s(null, "T %.2f K, P %.2f MPa, V %.0f A3, ", T, P * 1000, V);
    cout << null;
    logf << null;
    sprintf_s(null, "Energy (kcal/mol) Ek %.2f lambda_Ep %.2f Etot %.8f Ep %.8f", Ek, Ep, Etot, Ep / lambda);
    cout << null << endl;
    logf << null << endl;
    */
    //sprintf_s(null, "lambda %.8f Ep %.8f", lambda, Ep / lambda);
    //cout << null << endl;
    //logf << null << endl;
    cum_Ep += Ep / lambda;
    num_sample++;
    //printf("Tail contribution %.0f%% in Ep %.0f%% in P \n",Eptail/Ep*100,Ptail/P*100);
    
    //dump_trj();                     //output trajectory file
    return 0;
}

int MODEL::sample2()
{
    char null[1024];
    //calculation of system temperature and kinetic energy
    int i, j, k;

    tote();

    //calcualate pressure
//    double unitc=4184*1E21/6.0221367E23; //conver from kcal/mol/A3 to GPa
//    P=((Ek*2+strs[0]+strs[1]+strs[2])*unitc)/(3.0*V); //in GPa
//    P+= Ptail;

    //Display current information
    /*
    sprintf_s(null, "Current Step %d, accepted fraction %f, ", icycl, double(accpt) / icycl);
    cout << null;
    logf << null;
    sprintf_s(null, "T %.2f K, P %.2f MPa, V %.0f A3, ", T, P * 1000, V);
    cout << null;
    logf << null;
    sprintf_s(null, "Energy (kcal/mol) Ek %.2f lambda_Ep %.2f Etot %.8f Ep %.8f", Ek, Ep, Etot, Ep / lambda);
    cout << null << endl;
    logf << null << endl;
    */
    sprintf_s(null, "lambda %.8f Ep %.8f", lambda, cum_Ep / num_sample);
    cout << null << endl;
    logf << null << endl;
    //printf("Tail contribution %.0f%% in Ep %.0f%% in P \n",Eptail/Ep*100,Ptail/P*100);


    //dump_trj();                     //output trajectory file
    return 0;
}

/*
int MODEL::dump_trj()
{
    char null[1024];
    int i;
    sprintf_s(null, "Current Step %d, accepted fraction %f, ", icycl, double(accpt) / icycl);
    outf << null << endl; //comment in trj file
    sprintf_s(null, "%5d", npart);
    outf << null << endl;    //number of particles
    for (i = 0; i < npart; i++) { //position (nm) and velocity (nm/ps)
        sprintf_s(null, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f",
            1, "Ar", "Ar", i + 1, x[i][0] / 10, x[i][1] / 10, x[i][2] / 10, 0, 0, 0);
        outf << null << endl;
    }
    sprintf_s(null, "%10.5f%10.5f%10.5f", cell[0] / 10, cell[1] / 10, cell[2] / 10);
    outf << null << endl;  //box information
    return 0;
}
*/
