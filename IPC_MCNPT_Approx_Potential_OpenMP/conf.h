#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<fstream>
#include<iomanip>


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-14
#define RNMX (1.0-EPS)


//#define vectadd(a,b,c) c.x = b.x+a.x; c.y = b.y+a.y; //c.z = b.z+a.z;
//#define vectsub(a,b,c) c.x = b.x-a.x; c.y = b.y-a.y; //c.z = b.z-a.z;
#define vectinnerpdct(a,b) (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
//#define vectimes(a,b,c) c.x = a*b.x; c.y = a*b.y; //c.z=a*b.z;
using namespace std;

const int maxcolloids = 550;
const int maxbin = 500;

//code parameters
int  ind, start1, start2, anncyc;
long *idum;

double Q[4];
double quaternion(double[]);
double rotation(double[][3], double);
double ran2(long*);

class vector{
	public: double x,y,z;
		 vector();
		 vector(double, double, double);
		 void   setrandom(const double& box);
		 void   printvec();
		 void   pbc(double& boxi);
		 double norm();
		 double dot();
		 double dot(const vector&);
		 vector operator-() const;
		 vector operator-(const vector&);
		 vector operator+(const vector&);
		 vector vectimes(const double&, const double&, const double&);
		 vector vectimes(const double&);
		 //vector& operator=(const vector& v);
		 static vector zero();
};

//vector class function definitions
vector::vector() : x(0.0), y(0.0), z(0.0) { }

vector::vector(double x, double y, double z) : x(x), y(y), z(z) { }

void vector::setrandom(const double& box){
	this->x = (ran2(idum)-0.5)*box;	//Coordinates of first particle from -hbox to hbox
	this->y = (ran2(idum)-0.5)*box;
	//this->z = (ran2(idum)-0.5)*box;
	this->z = 0.0;
}
void vector::printvec(){
	cout << this->x << " " << this->y << " " << this->z << endl;
}
double vector::norm(){
	return sqrt(this->x*this->x + this->y*this->y + this->z*this->z);
}
double vector::dot(){
	return (this->x*this->x + this->y*this->y + this->z*this->z);
}
double vector::dot(const vector& v){
	return (this->x*v.x + this->y*v.y + this->z*v.z);
}
vector vector::operator-() const {
	return vector(-this->x, -this->y, -this->z);
}
vector vector::operator-(const vector& v) {
	return vector(this->x - v.x, this->y - v.y, this->z - v.z);
}
vector vector::operator+(const vector& v) {
	return vector(this->x + v.x, this->y + v.y, this->z + v.z);
}
vector vector::vectimes(const double& a, const double& b, const double& c){
	return vector(a*this->x, b*this->y, c*this->z);
}
vector vector::vectimes(const double& a){
	return vector(a*this->x, a*this->y, a*this->z);
}
//vector& vector::operator=(const vector& v) {
//	this->x = v.x; this->y = v.y; this->z = v.z;
//	return *this;
//}

//Class vector static definitions
vector vector::zero(){
	return vector(0.,0.,0.);
}
void vector::pbc(double& boxi)
{
	double hbox = boxi/2.0;
	if (this->x >  hbox) this->x -= boxi;
	if (this->x < -hbox) this->x += boxi;
	if (this->y >  hbox) this->y -= boxi;
	if (this->y < -hbox) this->y += boxi;
	//if (this->z >  hbox) this->z -= boxi;
	//if (this->z < -hbox) this->z += boxi;
}


vector colloid[maxcolloids];			//colloid centers
vector patch[maxcolloids];				//patch unit vector
vector rp[maxcolloids], rc[maxcolloids];	//patch n colloid off centers
int	ibox[maxcolloids];


void	readparameters();
void	readparameterscheck();
void	initialcondition(int, int);
void	equilibration(int, int);
void	production(int, int);
void	mcmove();
void	mcvol_npt();
void	mcvol_gemc();
void	mcswap();
void	writeconfig(long int, int, int);
void	writeconfig_inst(long int, int, int);
void	writexyz();
void	restart();
void	finalparameters();
double	totalenergy(int);
//double	totalenergy_vol(int, vector r[]);
//double	potenergy(int, vector, int);
double potenergychk(int, vector, vector, int);
double pair_energy(int,vector, vector, vector, vector, vector);
//double	potenergychk(int, vector, vector, int);
//double	ucor(int,int);
//double	pressurecalc(int);
vector mic(vector, int);
vector newpatch(vector, double);
//vector newpatch2d(vector);


bool initial, runrestart, runcontinue;

//simulation parameters

long int eqcycle, procycle;
long int f_cyc, f_mov, i_run, p_run; //f_cyc - first cycle, f_mov - first move and so on
long int nind;
int	 count1, count2, scalefact;
int	 nmov, volmov, swapmove, runstage;
int	 npart[2], npart_tot;
int	 eqframes, proframes;
int	 movetry[2], voltry, swaptry, swtry[2];
int	 moveaccept[2], volaccept, swapaccept;
double	 drmax[2], dvmax;

double	 box[2], rho[2];
double	 cutoff, temp, pressure;
double	 U[2], V[2], W[2], Wi;
double	 aven[2], avdens[2], simuP[2];
double	 mu[2], mucalc[2];

int maxit;
double chi, delta;
double tot_charge, zp, zc, ap, ac, term1, term2, constant, const_term2;
double kappa;
//double bessel_kappa[100], l_pow[100];
double ac_pow[100], ap_pow[100];

//for analytical
double legendre(int,double);
double bessel(int, double);
double factorial(int);
double n_plus_half_k(int,int);

//rho dist initialise
int npr, binindex;
double delrho;
double p[maxbin];

vector newp;

//runrestart starts a new set of simulation from a previously completed set
//runcontinue continues an interrupted simulation run

