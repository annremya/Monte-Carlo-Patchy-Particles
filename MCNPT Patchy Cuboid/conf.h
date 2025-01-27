#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<fstream>
#include<iomanip>

#define vectadd(a,b,c) c.x = b.x+a.x; c.y = b.y+a.y; c.z = b.z+a.z;
#define vectsub(a,b,c) c.x = b.x-a.x; c.y = b.y-a.y; c.z = b.z-a.z;
#define vectinnerpdct(a,b) (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
#define vectimes(a,b,c) c.x = a*b.x; c.y = a*b.y; c.z=a*b.z;
using namespace std;

const int nhis        = 200;

cube	*colloid;

void	readparameters();
void	initialcondition();
void	equilibration(int, int);
void	production(int, int);
void	mcmove(int);
void	mcvol();
void	writeconfig(long int, int, int);
void	writeconfig_inst(long int, int, int);
void	writexyz();
void	restart();
void	finalparameters();
void	widom(int,int);
void	widom_run();
void   distribution(int, int, int, int, double);
double	totalenergy();
double	potenergy(int, const cube&);
double pair_energy(const cube&, const cube&);

//code parameters
int  ind, start;

bool initial, runrestart, runcontinue;

//simulation parameters

long int eqcycle, procycle;
long int f_cyc, f_mov, i_run, p_run; //f_cyc - first cycle, f_mov - first move and so on
long int nind, count1, count2;
int	  nmov, volmov, runstage;
int	  npart, nA, nB;
int	  eqframes, proframes;
int	  movetry, voltry;
int	  moveaccept, volaccept;
int	  ighost;
double	  box, rho, volfrac;
double	  drmax, dvmax;
double	  rcut, rcutlj, r3, sigB3;	//cutoff-yukawa(patchy_SALR), rcutlj-lj(iso)
double   temp, epsilon_pi, epsilon_xb, epsilon_hb, pressure;
double	  U, V;
double   scalefact;
double	  aven, avdens, wdens, wdensA, wdensB;
double   wtest[2], muex[2];
int      itest[2];
//widom and overlapping dist
int      ntest[2], nreal[2];
double   umax, umin, factu;
double   p0[nhis][2], p1[nhis][2];

void vector::pbc(){
	double hbox = box/2.0;
	if      (this->x >  hbox) this->x -= box;
	else if (this->x < -hbox) this->x += box;
	if      (this->y >  hbox) this->y -= box;
	else if (this->y < -hbox) this->y += box;
	if      (this->z >  hbox) this->z -= box;
	else if (this->z < -hbox) this->z += box;
}

//runrestart starts a new set of simulation from a previously completed set
//runcontinue continues an interrupted simulation run
//start - type of initialcondition

