#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#define pi acos(-1)

using namespace std;

double lamb[3] = { 0.1127, 0.5, 0.88729 };

class MODEL
{
public:
    MODEL();              //class constructor
    ~MODEL();            //class destructor

    //model description
    int npart;               //number of particles
    double(*x)[3];      //position of particles (in Angstroms)
    double(*m);          //mass (in g/mol)
    double(*chg);		//charge
    double cell[3];      //unit cell length (in A)
    bool   period;       //periodicity 0:nonperiodic 1:periodic

    //MC parameters
    double Tset;         //set temperature 
    double Pset;         //set Pressure
    int ncycl;              //total steps to run  
    int icycl;               //current step 
    int avol;                //accepted volume change moves
    int ivol;                //total volume change attempts
    double Rcut;         //cutoff distance for vdw interaction
    double rc2;            //Rcut*Rcut;
    double delx;           //maxmum moving distance
    double delV;           //maximum volume change
    double beta;           //1/kT
    int accpt;             //number of accepted move

  //force field parameters
    double(*Do_atm);  //in kcal/mol
    double(*Ro_atm);  //in Angstrom
    double Do;
    double Ro;
    double p0;  //p0= 12*Do*pow(Ro,12) in unit kcal/mol*A12
    double p1;  //p1= 12*Do*pow(Ro,6) in unit kcal/mol*A6
    double Ro_Rc3;
    double Ro_Rc6;
    double Ro_Rc9;
    double Ro3;
    double unitc;
    double lambda;

    //bonding parameters
    double len0[3];
    double Kb;
    double Eb;
    double fb;
    double th0;
    double K;

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
    int mcvol();              //attempt for Monte Carlo Volume change
    int ener(int, double&);  //evaluation of interaction energy with one particle
    int tote();               //evaluate Ek, Ep, and Etot, stress
    int sample();             //calculation of properties
    int bond(int, double&);
    double myrand01(); //generate a random number between 0 and 1 [0,1]
    int dump_trj();           //output trajectory data to a file   
    ofstream outf;           //file stream of the trajectory
    ofstream logf;           //file stream of the log
    double cnt;
    int acptO;
    int acptH;
    int sample2();

    double num_sample;
    double cum_Ep;
};

int MODEL::bond(int i, double& eng)
{
    // calculate bond energy and force, ref. Lec.4 p.39
    // Ebond=0.5*Kb*(r-r0)^2
    // fbond=Kb*(r-r0)*x/r
    int j, k;
    double xr[3], redu, r2, len, d[3][3], dd[3], costh, h1, fbond;
    Kb = 547.5;
    //Kb=500;
    len0[0] = len0[1] = 1.0;
    len0[2] = 1.795;

    K = 120;

    th0 = 109.47;
    i = int(i / 3) * 3;

    //distance between O and H1
    d[0][0] = x[i + 1][0] - x[i][0];
    d[0][1] = x[i + 1][1] - x[i][1];
    d[0][2] = x[i + 1][2] - x[i][2];
    //distance between O and H2
    d[1][0] = x[i + 2][0] - x[i][0];
    d[1][1] = x[i + 2][1] - x[i][1];
    d[1][2] = x[i + 2][2] - x[i][2];
    //distance between H1 and H2
    d[2][0] = x[i + 1][0] - x[i + 2][0];
    d[2][1] = x[i + 1][1] - x[i + 2][1];
    d[2][2] = x[i + 1][2] - x[i + 2][2];

    for (j = 0; j < 2; j++) {
        dd[j] = d[j][0] * d[j][0] + d[j][1] * d[j][1] + d[j][2] * d[j][2];
        dd[j] = sqrt(dd[j]);
        //calculate bond stretching energy
        eng += 0.5 * Kb * (dd[j] - len0[j]) * (dd[j] - len0[j]);
        //cout<<0.5*Kb*(dd[j]-len0[j])*(dd[j]-len0[j])<<endl;
    }

    //calculate angle bending energy
    double dot = 0;
    for (int i = 0; i < 3; i++) dot = dot + d[0][i] * d[1][i];
    costh = dot / dd[0] / dd[1];
    double costh0, sinth0;
    costh0 = cos(th0 / 180 * pi);
    sinth0 = sin(th0 / 180 * pi);
    eng += 0.5 * K * (costh - costh0) * (costh - costh0) / sinth0 / sinth0;
    //cout<<"Ea = "<<0.5*K*(costh-costh0)*(costh-costh0)/sinth0/sinth0<<endl<<endl;
    return 0;
}

int main()
{
    for (int i = 0; i < 3; i++) {

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
    npart = 3 * 64;                              //number of particles
    Tset = 300;                                //target temperature in K
    ncycl = npart * 10000;//1e9;                        //steps to run
    Pset = 0.000101325;                                //GPa
    T = Tset;

    //allocation of memerory    
    x = new double[npart][3];        //position in Angstrom
    m = new double[npart];             //mass in g/mol
    chg = new double[npart];
    Do_atm = new double[npart];
    Ro_atm = new double[npart];

    //assign mass, charge, Do and Ro of particles
    int i, j, k, mm;
    for (i = 0; i < npart; i++) {
        if (i % 3 == 0) {		//oxygen properties
            m[i] = 15.9994;
            chg[i] = -0.82;
            Do_atm[i] = 0.1848;    // in kcal/mol
            Ro_atm[i] = 3.5532;    // in Angstrom
        }
        else {			//hydrogen properties
            m[i] = 1.00784;
            chg[i] = 0.41;
            Do_atm[i] = 0.01;    // in kcal/mol
            Ro_atm[i] = 0.9;    // in Angstrom
        }
    }
    //One simple way to place particles in space
    double sep = 4;  //separation distance between two particles
    int nt = 0;

    //give force field parameters

    cell[0] = cell[1] = cell[2] = pow(npart / 3 * 29.40, 1.0 / 3.0);      //length of unit cell in Angstroms
    period = 1;                                //flag for periodicity
    Rcut = 10;                          //cutoff distance for vdw interaction
    rc2 = Rcut * Rcut;                       //Rcut2 

    delx = 0.07; // angstrom
    accpt = 0;
    beta = 4184 / (8.314 * T); // in unit of 1/(kcal/mol)

    if (period) { //periodic system
      //A rough method to place particles in the box
        double delta[3];
        int pps;
        pps = int(pow(npart / 3, 1.0 / 3.0)) + 1;            //particles per side
        for (k = 0; k < 3; k++) delta[k] = cell[k] / pps; //spacing of particles
        cout << delta[0] << endl;
        for (i = 0; i < pps; i++) {
            for (j = 0; j < pps; j++) {
                for (k = 0; k < pps; k++) {
                    for (mm = 0; mm < 3; mm++) {//in every block, 3 particles
                        if (nt % 3 == 0) {
                            x[nt][0] = i * delta[0];
                            x[nt][1] = j * delta[1];
                            x[nt][2] = k * delta[2];
                        }
                        else if (nt % 3 == 1) {
                            x[nt][0] = i * delta[0];
                            x[nt][1] = cos(-19.47 / 180 * pi) + j * delta[1];
                            x[nt][2] = sin(-19.47 / 180 * pi) + k * delta[2];
                        }
                        else {
                            x[nt][0] = i * delta[0];
                            x[nt][1] = j * delta[1];
                            x[nt][2] = 1.0 + k * delta[2];
                        }
                        nt++;
                    }
                    if (nt == npart) i = j = k = pps;        //nt has reached npart, exit loops
                }
            }
        }
    }
    df = 3 * npart;                                    //degree of freedom

   //calculate volume
    V = cell[0] * cell[1] * cell[2];   //in A^3

    //outf.open("MC_NVTwater.gro", ios::out);   //output filename for trajectory
    //logf.open("MC_NVTwater.log", ios::out);   //output filename for trajectory

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
    double r2, r2i, r6i, xr[3], redu;//, lambda = 1;
    energy = 0;

    for (j = 0; j < npart; j++) {
        if (j == i) continue;			//skip atom itself
        //combination rule
        Ro = 0.5 * (Ro_atm[i] + Ro_atm[j]);
        Do = sqrt(Do_atm[i] * Do_atm[j]);
        //exclusion
        if (i % 3 == 0) {
            if (j == i + 1 || j == i + 2) {
                continue;
            }
        }
        else if (i % 3 == 1) {
            if (j == i - 1 || j == i + 1) {
                continue;
            }
        }
        else {
            if (j == i - 2 || j == i - 1) {
                continue;
            }
        }

        p0 = 12 * Do * pow(Ro, 12); //in unit kcal/mol*A12
        p1 = 12 * Do * pow(Ro, 6);   //in unit kcal/mol*A6

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
            energy += lambda * ((p0 * r6i - p1 * 2) * r6i / 12.0 + chg[i] * chg[j] * 332.0637 * sqrt(r2i));

            /*
            cout<<"for "<<i<<" & "<<j<<endl
             <<"vdw = "<<(p0*r6i - p1*2)*r6i/12.0<<endl
             <<"cou = "<<chg[i]*chg[j]*332.0637*sqrt(r2i)<<endl<<endl;*/
        }
    }
    //cout<<"ener result"<<endl;
    bond(i, energy);
    return 0;
}

int MODEL::tote()
{   //this is almost the same as the energy calculation in force()
    int i, j, k, nbox;
    double r2, r2i, r6i, xr[3], redu, ff, cou;//, lambda = 1;
    Ep = 0;

    for (i = 0; i < npart; i++) {
        if (i % 3 == 0) {
            bond(i, Ep);
        }
        for (j = i + 1; j < npart; j++) {
            //combination rule
            Ro = 0.5 * (Ro_atm[i] + Ro_atm[j]);
            Do = sqrt(Do_atm[i] * Do_atm[j]);

            //exclusion
            if (i % 3 == 0) {
                if (j == i + 1 || j == i + 2) {
                    continue;
                }
            }
            else if (i % 3 == 1) {
                if (j == i - 1 || j == i + 1) {
                    continue;
                }
            }
            else {
                if (j == i - 2 || j == i - 1) {
                    continue;
                }
            }
            p0 = 12 * Do * pow(Ro, 12); //in unit kcal/mol*A12
            p1 = 12 * Do * pow(Ro, 6);   //in unit kcal/mol*A6

            for (k = 0; k < 3; k++) {
                xr[k] = x[i][k] - x[j][k];  //distance vector
                //cout<<x[i][k]<<"	"<<x[j][k]<<endl;
            }
            //cout<<"original xr = "<<xr[0]<<", "<<xr[1]<<", "<<xr[2]<<endl;

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
                // consider atom j at origin, the force on atom i at some position r
                // fx = - dU/dx = -(dU/dr)(dr/dx)= - (x/r)(dU/dr)
                // U = Do ( (Ro/r)^12 - 2 (Ro/r)^6) )
                // dU/dr = -12 Do/r ( (Ro/r)^12 - (Ro/r)^6) )
                // fx = 12 x Do/r^2 ( (Ro/r)^12 - (Ro/r)^6) )
                //    =  x ( 12DoRo^12/r^6 - 12DoRo^6 )/r^6 /r^2
                ff = r6i * (p0 * r6i - p1) * r2i;          //in unit of kcal/mol/A2
                //Ecou = kQQ/eR, f = -dU/dx  = -(x/r)(dU/dr) = -(x/r)*(-kQQ/er^2) = kQQ/e/r^3*x
                cou = chg[i] * chg[j] * 332.0637 * sqrt(r6i);

                Ep += lambda * ((p0 * r6i - p1 * 2) * r6i / 12.0 + chg[i] * chg[j] * 332.0637 * sqrt(r2i));

                /*
                cout<<"for "<<i<<" & "<<j<<endl
                <<"vdw = "<<(p0*r6i - p1*2)*r6i/12.0<<endl
                <<"cou = "<<chg[i]*chg[j]*332.0637*sqrt(r2i)<<endl<<endl; */
            }
        }
    }
    //Ek=df/2*RT
    Ek = 0.5 * df * 8.314 * T / 4184; //in kcal/mol   

    //Tail correction for Ep
    Ro_Rc3 = pow(Ro / Rcut, 3.0);
    Ro_Rc6 = Ro_Rc3 * Ro_Rc3;
    Ro_Rc9 = Ro_Rc6 * Ro_Rc3;
    Ro3 = Ro * Ro * Ro;
    Eptail = (4.0 / 3.0) * 3.1415926 * npart * npart / V * Do * Ro3 * (Ro_Rc9 / 6.0 - Ro_Rc3);
    //Ep += Eptail;    
    Etot = Ek + Ep;

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
    sprintf_s(null, "Energy (kcal/mol) kinetic %.2f potential %.2f total %.8f", Ek, Ep / npart * 3, Etot);
    cout << null << endl;
    logf << null << endl;
    */
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
        if (i % 3 == 0) {
            sprintf_s(null, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f",
                1, "O", "O", i + 1, x[i][0] / 10, x[i][1] / 10, x[i][2] / 10, 0, 0, 0);
        }
        else {
            sprintf_s(null, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f",
                1, "H", "H", i + 1, x[i][0] / 10, x[i][1] / 10, x[i][2] / 10, 0, 0, 0);
        }
        outf << null << endl;
    }
    sprintf_s(null, "%10.5f%10.5f%10.5f", cell[0] / 10, cell[1] / 10, cell[2] / 10);
    outf << null << endl;  //box information
    return 0;
}
*/