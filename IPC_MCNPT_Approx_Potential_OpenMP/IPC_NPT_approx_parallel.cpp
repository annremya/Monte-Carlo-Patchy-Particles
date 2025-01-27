//--------------------------------------------------------------------------------------------------------------------------
//	GEMC + NPT code written for the system of IPC particles based modeling done for Single patch IPC
//	Single Patch IPC approximated model - parallel code
//	Based on algorithms from UMS (Frenkel and Smit)
//	Written by Remya Ann on 1/11/2020
//	Includes a parameter input file and configuartion header file
//--------------------------------------------------------------------------------------------------------------------------

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string>
#include<iomanip>
#include"conf.h"

#include<omp.h>
//set number of threads in job file or before compiling as "export OMP_NUM_THREADS=8"

using namespace std;

int	 nadj      	= 5000;
int	 nsample_e 	= 500; 
int	 nsample_p 	= 500; 
int 	 eqconfig  	= 10000;
int	 proconfig	= 10000;

//int	 shift		= 0;
//int 	 tailcorr	= 1;

//for constant
double qe    = 1.60217662e-19;
double kb    = 1.3807e-23;
double epso  =  8.85e-12;
double epsr  = 80;
double T     = 298.;
double sigma = 1.5e-6;

//Main program
//--------------------------------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])

{
	srand(time(0));
	if(argc > 1){
		ind    = atoi(argv[1]);
		anncyc = atoi(argv[2]);
	}
	readparameters();
	cout << "chi= " << chi << "  " <<	M_PI << endl;

	//check for psi calc
	colloid[0] = vector(0.,0.,0.);
	colloid[1] = vector(2.,0.,0.);
	patch[0] = vector(1.,0.,0.);
	patch[1] = vector(-1.,0.,0.);
	for(int i=0;i<2;i++){
		rp[i]   = colloid[i] +  patch[i].vectimes(ap);
		rc[i]   = colloid[i] + -patch[i].vectimes(ac);
	}
	cout << "pair energy pp = " << pair_energy(0,colloid[0],patch[0],colloid[1],rp[1],rc[1]) << endl;


//	cut2	= cutoff*cutoff;
//	cut3	= cut2*cutoff;
//	invcut3 = 1./cut3;
//	invcut6 = invcut3*invcut3;
//	Ucut	= 4.*invcut6*(invcut6-1.);

	for(int iann=0; iann<anncyc; iann++)
	{

		volmov   = 1;
		swapmove = 1000;//10;
//		nmov     = npart_tot + swapmove + volmov;		//GEMC//Number of moves in each 
		nmov     = npart_tot + volmov;			//NPT //Number of moves in each 
		readparameterscheck();
		cout << "Annealing cycle : " << iann << endl;
		//cout << "Scaling factor is : " << scalefact << " and epsilon is : " << epsilon << " " << epsilon1 << endl;
		count1 = 0;
		count2 = 0;

		char folder[100];
		sprintf(folder,"mkdir run%d_%d_%.2f", ind, scalefact, temp);
		system(folder);
		sprintf(folder,"mkdir run%d_%d_%.2f/eq", ind, scalefact, temp);
		system(folder);
		sprintf(folder,"mkdir run%d_%d_%.2f/pro", ind, scalefact, temp);
		system(folder);
		sprintf(folder,"mv parameter_check_%d_%.2f.dat run%d_%d_%.2f/parameter_check_%d_%.2f.dat", ind, temp, ind, scalefact, temp, ind, temp);
		system(folder);

		long y    = -rand();			//seed for random number generator
		     idum = &y;


		if (iann == 0)
		{
			//------------initial condition---------------------------------
			cout << initial << runrestart << runcontinue << endl;

			if      (initial) 
			{
				cout << "initial" << endl;
				initialcondition(start1, 0);
				initialcondition(start2, 1);
				U[0] = totalenergy(0);
				U[1] = totalenergy(1);
				cout << "Total energy of box 1 = " << U[0]/(double(npart[0])) << endl;
				cout << "Total energy of box 2 = " << U[1]/(double(npart[1])) << endl;
				equilibration(0,0);
				production(0,0);
			}
			else if (runrestart) 
			{
				cout << "runrestart" << endl;
				restart();
				U[0] = totalenergy(0);
				U[1] = totalenergy(1);
				cout << "Total energy of box 1 = " << U[0]/(double(npart[0])) << endl;
				cout << "Total energy of box 2 = " << U[1]/(double(npart[1])) << endl;
				equilibration(0,0);
				production(0,0);

			}
			else if (runcontinue)
			{
				cout << "runcontinue" << endl;
				restart();
				U[0] = totalenergy(0);
				U[1] = totalenergy(1);
				cout << "Total energy of box 1 = " << U[0]/(double(npart[0])) << endl;
				cout << "Total energy of box 2 = " << U[1]/(double(npart[1])) << endl;

				if(runstage == 1)
				{
					equilibration(f_cyc,f_mov);
					production(0,0);
				}
				else
				{
					production(f_cyc,f_mov);
				}
			}
		}

		else{
			U[0] = totalenergy(0);
			U[1] = totalenergy(1);
			cout << "Total energy of box 1 = " << U[0]/(double(npart[0])) << endl;
			cout << "Total energy of box 2 = " << U[1]/(double(npart[1])) << endl;
			writexyz();
			writeconfig(0,0,0);
			writeconfig_inst(0,0,0);
			equilibration(0,0);
			production(0,0);
		}

		//writeconfig((eqcycle+procycle)*nmov);
		writexyz();
		finalparameters();

		//mucalc[0] = (-log(mu[0]/double(swtry[0]))*temp);
		//mucalc[1] = (-log(mu[1]/double(swtry[1]))*temp);
		if((swtry[0]  != 0)&&(mu[0]  != 0)) {mucalc[0] = (-log(mu[0]/double(swtry[0])))*temp;} else {mucalc[0] = 0;}
		if((swtry[1]  != 0)&&(mu[1]  != 0)) {mucalc[1] = (-log(mu[1]/double(swtry[1])))*temp;} else {mucalc[1] = 0;}
		ofstream log;

		char filename[100];
		sprintf(filename, "run%d_%d_%.2f/log.dat", ind, scalefact, temp);
		log.open(filename, ofstream::out | ofstream::app);
		log << std::fixed;
		log << std::setprecision(4);
		if (log.is_open())
		{
			log << "\n\t\t---------Box 1------------\n";
			if(npart[0] != 0)
			{
			log << "Average Energy per particle      = " << aven[0]/(double(npart[0])*double(nind))	<< "\n"
			    << "Running Energy per particle      = " << U[0]/double(npart[0])				<< "\n";
			}
			else
			{
			log << "Average Energy per particle      = " << "0."						<< "\n"
			    << "Running Energy per particle      = " << "0."						<< "\n";
			}

			log << "Average density of the system    = " << avdens[0]/double(nind)			<< "\n"
			    << "Chemical Potential of the system = " << mucalc[0]					<< "\n";

			log << "\n\t\t---------Box 2------------\n";
			if(npart[1] != 0)
			{
			log << "Average Energy per particle      = " << aven[1]/(double(npart[1])*double(nind))	<< "\n"
			    << "Running Energy per particle      = " << U[1]/double(npart[1])				<< "\n";
			}
			else
			{
			log << "Average Energy per particle      = " << "0."						<< "\n"
			    << "Running Energy per particle      = " << "0."						<< "\n";
			}
			log << "Average density of the system    = " << avdens[1]/double(nind)			<< "\n"
			    << "Chemical Potential of the system = " << mucalc[1]					<< "\n";
		}
		else{cerr   << "unable to open log file. \n";}

		//scalefact++;
		//epsilon = epsilon1*pow(1.05,scalefact);
		drmax[0]= 0.1; drmax[1] = 0.1; dvmax = 0.1;
	}
	return 0;

} //End of main program
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void readparameters()
{
	char parafile[100], dum[100];
	sprintf(parafile, "para.%d.inp",ind);
	ifstream infile;
	infile.open(parafile, ifstream::in);
	if (infile.is_open())
	{
	  infile >> npart[0] 	>> dum
		  >> npart[1] 	>> dum
		  >> chi		>> dum
		  >> rho[0]  		>> dum
		  >> rho[1]  		>> dum
		  >> kappa 		>> dum
		  >> temp   		>> dum
		  >> pressure		>> dum
		  >> drmax[0]  	>> dum
		  >> drmax[1]  	>> dum
		  >> dvmax		>> dum
		  >> eqcycle		>> dum
		  >> procycle		>> dum
		  >> runcontinue	>> dum
		  >> initial		>> dum
		  >> runrestart	>> dum
		  >> eqframes   	>> dum
		  >> proframes  	>> dum
		  >> start1		>> dum
		  >> start2		>> dum
		  >> maxit		>> dum
		  >> scalefact	>> dum;
	}
	else { cerr << "unable to open parameter file.\n";}
	infile.close();

	delta     = acos(1.0-(2.*chi));		//M_PI/2.;	//chi*M_PI;	//acos(1.0-(2.*chi)); //chi is patch coverage

	V[0]      = double(8.*npart[0])/rho[0];	//pow(box[0],2.);	//multiply with 8 since 2sigma is the dia
	V[1]      = double(8.*npart[1])/rho[1];	//pow(box[1],2.);
	//box[0]    = pow(V[0],1./3.);
	//box[1]    = pow(V[1],1./3.);
	box[0]    = pow(V[0],1./2.);
	box[1]    = pow(V[1],1./2.);
	npart_tot = npart[0] + npart[1];
//	cutoff    = 1.+ 5.*lambda;
	cutoff    = 5.;
//	epsilon   = epsilon1*pow(1.05,scalefact);
	cout   << "Densities of simulation boxes are   : " << rho[0] << " and " << rho[1] << endl;
	cout   << "Side lengths of simulation boxes are: " << box[0] << " and " << box[1] << endl;

	rho[0]    = double(npart[0])/V[0];
	rho[1]    = double(npart[1])/V[1];
	if(V[0]==0) rho[0] = 0;
	if(V[1]==0) rho[1] = 0;
	cout   << "Densities of simulation boxes are   : " << rho[0] << " and " << rho[1] << endl;
	cout   << "cutoff = " << cutoff << endl;

	tot_charge  = 0.02*2.*sigma/(0.715*1e-9);
	ac          = chi; //R=1
	ap          = 1.-chi;
	zc          = -(1.-chi)*tot_charge;
	zp          = chi*tot_charge;
	term2       = exp(kappa)/(1.+kappa);
	constant    = (qe*qe)/(4.*M_PI*epso*epsr*kb*T*sigma);
	const_term2 = constant*term2*term2;

	for(int l=0; l<=maxit; l++){
//		bessel_kappa[l] = bessel(l,kappa);
//		ac_pow[l]       = pow((-1.*ac),double(l));
//		ap_pow[l]       = pow((ap),double(l));
//		l_pow[l]        = double(2*l+1);
		ac_pow[l] = pow((-1.*ac),double(l))*double(2*l+1);
		ap_pow[l] = pow((ap),double(l))*double(2*l+1);
	}

}

void readparameterscheck()
{
	char parafile[100];
	ofstream outfile;
	sprintf(parafile, "parameter_check_%d_%.2f.dat", ind,temp);
	outfile.open(parafile, ofstream::out | ofstream::trunc);
	if (outfile.is_open())
	{
	  outfile<< "index       = "	<< ind		<< "\n"
		  << "npart1      = " 	<< npart[0]	<< "\n"
		  << "npart2      = " 	<< npart[1]	<< "\n"
		  << "npart_tot   = " 	<< npart_tot	<< "\n"
		  << "chi         = "  	<< chi		<< "\n"
		  << "rho1        = "  	<< rho[0]	<< "\n"
		  << "rho2        = "  	<< rho[1]	<< "\n"
		  << "box1        = "	<< box[0]	<< "\n"
		  << "box2        = "	<< box[1]	<< "\n"
		  << "Volume1     = "	<< V[0]	<< "\n"
		  << "Volume2     = "	<< V[1]	<< "\n"
		  << "kappa       = "	<< kappa 	<< "\n"
		  << "cutoff      = "	<< cutoff 	<< "\n"
		  << "temp        = "	<< temp	<< "\n"
		  << "pressure    = "	<< pressure	<< "\n"
		  << "drmax1      = " 	<< drmax[0]	<< "\n"
		  << "drmax2      = " 	<< drmax[1]	<< "\n"
		  << "dvmax       = "	<< dvmax	<< "\n"
		  << "eqcycle     = "	<< eqcycle	<< "\n"
		  << "procycle    = "	<< procycle	<< "\n"
		  << "runcontinue = "	<< runcontinue<< "\n"
		  << "initial     = "	<< initial	<< "\n"
		  << "runrestart  = "	<< runrestart	<< "\n"
		  << "eqframes    = "	<< eqframes 	<< "\n"
		  << "proframes   = "	<< proframes	<< "\n"
		  << "start1      = "	<< start1	<< "\n"
		  << "start2      = "	<< start2	<< "\n"
		  << "maxit       = "	<< maxit	<< "\n"
		  << "scalefact   = "	<< scalefact	<< "\n";
	}
	else{cerr << "Unable to open file for parameter check output. \n";}
	outfile.close();
	
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void initialcondition(int start, int a) // a = box index
{

	vector t, rij, patchi;
	double	    r2, r1;
	int	    first, last;
	
	if      (a == 0) {first = 0;	 last = npart[0];}
	else if (a == 1) {first = npart[0]; last = npart_tot;}

	if (start == 1)					//Random initial configuration
	{
		cout << "Random initial Configuration for box " << a+1 << endl;
		colloid[first].setrandom(box[a]);		//Coordinates of first particle from -hbox to hbox
		for (int i=first+1; i<last; i++)		//Coordinates of rest of the particles with overlap check
		{
		    do	{
				t.setrandom(box[a]);
				for (int j=first; j<i; j++)
				{
					rij = colloid[j]-t;
					//rij.printvec();
					//vectsub(colloid[j],t,rij);
					//rij.printvec();
					rij.pbc(box[a]);	//Minimum image conversion
					//r2  = vectinnerpdct(rij,rij);
					//r1  = sqrt(r2);
					r1 = rij.norm();
					if (r1<=2.0) break;
				}   
			}   while (r1<2.0);		//Hard Sphere diameter is 2.0

			colloid[i] = t;
			//colloid[i].printvec(); t.printvec();
			cout << i << endl;
		}
	}

	else if (start == 2)				//Simple cubic 3D
	{
		cout << "Simple Cubic Configuration for box " << a+1 << endl;
		int	index;
//		int 	n_unitcell	= ceil(pow(npart[a],(1./3.)));
		int 	n_unitcell	= ceil(pow(npart[a],(1./2.)));
		double	lcell		= box[a]/n_unitcell;
		double	hbox		= box[a]/2.;

		cout << "Unit cell length = \t" << lcell << endl;

		if      (a == 0) {index = 0;}
		else if (a == 1) {index = npart[0];}

		//for(int iz=0; iz<n_unitcell; iz++)
		//{   
		    for(int iy=0; iy<n_unitcell; iy++)
		    {   for(int ix=0; ix<n_unitcell; ix++)
			{
				if (index >= last) {break;}		// make sure only N particles are assigned positions
				colloid[index].x = double(ix)*lcell;
				colloid[index].y = double(iy)*lcell;
				//colloid[index].z = double(iz)*lcell;
				colloid[index].z = 0.0;
				index	         += 1;
			}
		    }
		//}
		for(int i=first; i<last; i++)		//Shifting centre of box to origin
		{
			colloid[i].x	-= hbox;
			colloid[i].y	-= hbox;
			colloid[i].z	-= hbox;
		}
	}

	else if (start == 3)				//Triangular 2d
	{
		cout << "Triangular lattice Configuration \n";
		int 	n_unitcell	= sqrt(2.1)*sqrt(double(npart[a])/2.);
		double	lcell		= box[a]/n_unitcell;
		double  lcell2	= lcell/2.;
		double	hbox		= box[a]/2.;

		cout << "Unit cell length = \t" << lcell << endl;

		//sublattice A
		colloid[first].x = 0.0;
		colloid[first].y = 0.0;

		//sublattice B
		colloid[first+1].x = lcell2;
		colloid[first+1].y = sqrt(3.)*lcell2;

		int a=0, index=0;
		for(int j=0; j<n_unitcell; j++)
		{
		    for(int i=0; i<n_unitcell; i++)
		    {
			for(int ind=first; ind<2; ind++)
			{
				if (index >= last) {break;}		// make sure only N particles are assigned positions
				colloid[ind+a].x = colloid[ind].x + (lcell*double(i));
				colloid[ind+a].y = colloid[ind].y + (sqrt(3.)*lcell*double(j));
				index 		+= 1;
			}
			a += 2;
		    }
		}
		for(int i=first; i<last; i++)		//Shifting centre of box to origin
		{
				colloid[i].x	-= hbox;
				colloid[i].y	-= hbox;
		}
	}

	for(int i=first; i<last; i++)			//Orientations of the particles
	{
		patch[i] = vector(1.,0.,0.);
		//patch[i] = newpatch2d(patch[i]);
		patch[i] = newpatch(patch[i], 0.7854);		//45 degress
		rp[i]   = colloid[i] +  patch[i].vectimes(ap);
		rc[i]   = colloid[i] + -patch[i].vectimes(ac);
		//cout << patch[i].x << "\t" << patch[i].y << "\t" << patch[i].z << endl;
	}
	for(int i=first; i<last; i++)			//Box index of particles
	{
		ibox[i] = a;
	}

	if (a == 1)
	{
		writexyz();
		writeconfig(0,0,0);
	}
	return;

}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void restart()
{
	long int  nrun;
	char filename[100];
	ifstream restartconfig;

	sprintf(filename, "config_inst%d.dat", ind);
	restartconfig.open(filename, ifstream::in);
	if(restartconfig.is_open())
	{
	   restartconfig >> nrun >> f_cyc >> f_mov >> runstage >> npart[0] >> npart[1] >> box[0] >> box[1] >> drmax[0] >> drmax[1] >> dvmax >> U[0];
	   for (int j=0; j<npart_tot; j++)
	   {
		restartconfig >> colloid[j].x >> colloid[j].y >> colloid[j].z
			       >> patch[j].x   >> patch[j].y   >> patch[j].z
			       >> ibox[j];
	   }
	}
	else{cerr << "unable to open config_inst file for input. \n";}
	restartconfig.close();
//	V[0]      = pow(box[0],3.);
//	V[1]      = pow(box[1],3.);
	V[0]      = pow(box[0],2.);
	V[1]      = pow(box[1],2.);
	rho[0]	   = double(npart[0])/V[0];
	rho[1]	   = double(npart[1])/V[1];
	if(V[0]==0) rho[0] = 0;
	if(V[1]==0) rho[1] = 0;
	npart_tot = npart[0] + npart[1];
	cout   << "Densities of simulation boxes are: " << rho[0] << " and " << rho[1] << endl;

	for(int i=0; i<npart_tot; i++)
	{
		rp[i]   = colloid[i] +  patch[i].vectimes(ap);
		rc[i]   = colloid[i] + -patch[i].vectimes(ac);
	}

	double r2, r1;
	vector rij;

//	for (int i=0; i<npart_tot; i++)
//	{
//		for (int j=i+1; j<npart_tot; j++)
//		{
//		     if(ibox[i] == ibox[j])
//		     {
//			//vectsub(colloid[i],colloid[j],rij);
//			rij = colloid[i] - colloid[j];
//			rij.pbc(box[ibox[i]]);		//Minimum image conversion
//			//r2  = vectinnerpdct(rij,rij);
//			r1  = rij.norm();
//			if (r1 <= 1.0)
//			{
//				cout << "Overlap  " << i << "\t" << j << "\t" << r1 << " in box " << ibox[i]+1 << endl;
//				cout << colloid[i].x << "  " << colloid[j].x << endl;
//				cout << "Press Enter to continue" << endl;
//				cin.ignore(); 
//			}
//		     }					
//		}
//	}

	writeconfig(nrun,f_cyc,f_mov);
	writexyz();
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void equilibration(int f_cyc, int f_mov)
{

	ofstream sample, accept, probplot;
	char filename[100];					//output file for sample 
	sprintf(filename, "run%d_%d_%.2f/sample.dat", ind, scalefact, temp);
	sample.open(filename, ofstream::out | ofstream::app);
	sample << std::fixed;
	sample << std::setprecision(4);
	if (sample.is_open())
	{
		sample  << left
			<< setw(10) << "ncycle"		<< "\t"
			<< setw(10) << "toten1/N*eps"	<< "\t"
			<< setw(10) << "toten1/N*eps"	<< "\t"
			<< setw(10) << "avg_dens1"		<< "\t"
			<< setw(10) << "avg_dens2"		<< "\t"
			<< setw(10) << "mu_calc1"		<< "\t"
			<< setw(10) << "mu_calc2"		<< "\n";
	}
	else{cerr	<< "unable to open sample file. \n";}


	sprintf(filename, "run%d_%d_%.2f/accept.dat", ind, scalefact, temp);	//output file for acceptance
	accept.open(filename, ofstream::out | ofstream::app);
	accept << std::fixed;
	accept << std::setprecision(4);
	if (accept.is_open())
	{
		accept  << left
			<< setw(10)<< "ncycle"		<< "\t"
			<< setw(8) << "moveaccept1"  	<< "\t"
			<< setw(7) << "drmax1"      	<< "\t"
			<< setw(8) << "moveaccept2" 	<< "\t"
			<< setw(7) << "drmax2"       	<< "\t"
			<< setw(8) << "volaccept"  		<< "\t"
			<< setw(7) << "dvmax" 		<< "\t"
			<< setw(7) << "swaptry"		<< "\t"
			<< setw(7) << "swapaccept" 		<< "\n";
	}
	else{cerr	<< "unable to open acceptance file. \n";}

	sprintf(filename, "run%d_%d_%.2f/probplot_eq.dat", ind, scalefact, temp);
	probplot.open(filename, ofstream::out | ofstream::app);
	probplot << std::fixed;
	probplot << std::setprecision(4);

	nind = 0;
	int k, ran;
	double insE1[5], insE2[5], sumE[2], insD1[5], insD2[5], sumD[2], avenergy[2];
	double acceptance[2] = {0.,0.}, acceptancev = 0., acceptancesw = 0.;
	avenergy[0] = 0., avenergy[1] = 0;
	mu[0] 	  = 0.;
	mu[1] 	  = 0.;
	swtry[0] = 0.;
	swtry[1] = 0.;

	for(int icyc=f_cyc; icyc<eqcycle; icyc++)
	{
	    if (icyc%100 == 0) {cout << "progress(eq) " << icyc << "\t" << double(icyc)/double(eqcycle) << endl; }
	    //nmov     = npart_tot + swapmove + volmov;
	    for(int imov=f_mov; imov<nmov; imov++)
	    {	
		i_run = (icyc*nmov)+(imov+1);
		ran   = round(ran2(idum)*(nmov));
		if      (ran < npart_tot)
		{
			//cout << i_run << "\t" << U1 << endl;	
			mcmove();
		}
		else if (ran < (npart_tot+volmov))
		{
			mcvol_npt();
		}
//		else
//		{
//			mcswap();
//		}


	    	if( i_run % ((eqcycle*nmov)/eqframes) == 0 ){writexyz();}	//writing a frame for visualization
	    	if( i_run % ((eqcycle*nmov)/eqconfig) == 0 ){writeconfig(i_run,icyc,imov);}
	
	    	if( i_run % nadj == 0 )						//updating drmax to adjust acceptance ratio	
	    	{		
			if (movetry[0]!=0) {acceptance[0] = double(moveaccept[0])/double(movetry[0]);}
			if (acceptance[0] < 0.4)      {drmax[0] /= 1.05;}
			if (acceptance[0] > 0.6)      {drmax[0] *= 1.05;}
			if (drmax[0]      > 0.5)      {drmax[0]  = 0.5;}
			if (drmax[0]     < 1.0e-03)   {drmax[0]  = 1.0e-02;}

			if (movetry[1]!=0) {acceptance[1] = double(moveaccept[1])/double(movetry[1]);}
			if (acceptance[1] < 0.4)      {drmax[1] /= 1.05;}
			if (acceptance[1] > 0.6)      {drmax[1] *= 1.05;}
			if (drmax[1]      > 0.5)      {drmax[1]  = 0.5;}
			if (drmax[1]      < 1.0e-03)  {drmax[1]  = 1.0e-02;}

			if (voltry!=0) {acceptancev = double(volaccept)/double(voltry);}
			if (acceptancev   < 0.4)      {dvmax    /= 1.05;}
			if (acceptancev   > 0.6)      {dvmax    *= 1.05;}
			if (dvmax         > 0.5)      {dvmax     = 0.5;}
			if (dvmax         < 1.0e-03)  {dvmax     = 1.0e-02;}

			if (swaptry!=0) {acceptancesw = double(swapaccept)/double(swaptry);}
//			if (acceptancesw < 0.02)      {swapmove+=100;}
//			if (acceptancesw > 0.02)      {swapmove-=100;}
//			if (swapmove < 0 )            {swapmove =100;}

			if (accept.is_open())
			{
				accept  << left
					<< setw(10) << i_run		<< "\t"
					<< setw(8) << acceptance[0] << "\t"
					<< setw(7) << drmax[0]      << "\t"
					<< setw(8) << acceptance[1] << "\t"
					<< setw(7) << drmax[1]      << "\t"
					<< setw(8) << acceptancev  	<< "\t"
					<< setw(7) << dvmax     	<< "\t"
					<< setw(7) << swaptry     	<< "\t"
					<< setw(7) << acceptancesw	<< "\n";
			}
			else{cerr	<< "unable to open acceptance file for output. \n";}
			movetry[0]	= 0;
			movetry[1]	= 0;
			voltry		= 0;
			swaptry	= 0;
			moveaccept[0] = 0.0;
			moveaccept[1] = 0.0;	
			volaccept	= 0.0;
			swapaccept	= 0.0;
	    	}

	    	if( i_run > (nind*nsample_e - 5) )
	    	{
			if (npart[0] != 0)
			{
				k        = i_run%5;
				insE1[k] = U[0]; 
				insD1[k] = double(npart[0])/V[0];
			}
			else
			{
				k        = i_run%5;
				insE1[k] = 0.0; 
				insD1[k] = 0.0;
			}
			if (npart[1] != 0)
			{
				k        = i_run%5;
				insE2[k] = U[1]; 
				insD2[k] = double(npart[1])/V[1];
			}
			else
			{
				k        = i_run%5;
				insE2[k] = 0.0; 
				insD2[k] = 0.0;
			}

			
	    	}
	    	if( i_run % nsample_e == 0 )
	    	{
			//cout << "sampling" << endl;
			sumE[0]   = 0.; sumE[1]  = 0.;
			for (int j=0; j<5; j++) {sumE[0] += insE1[j]; sumE[1] += insE2[j];}
			sumE[0]  /= 5.; sumE[1] /= 5.;

			sumD[0]   = 0.; sumD[1]  = 0.;
			for (int j=0; j<5; j++) {sumD[0] += insD1[j]; sumD[1] += insD2[j];}
			sumD[0]  /= 5.; sumD[1] /= 5.;

			writeconfig_inst(i_run,icyc,imov);
			nind += 1;

			if((swtry[0]  != 0)&&(mu[0]  != 0)) {mucalc[0] = (-log(mu[0]/double(swtry[0])))*temp;} else {mucalc[0] = 0;}
			if((swtry[1]  != 0)&&(mu[1]  != 0)) {mucalc[1] = (-log(mu[1]/double(swtry[1])))*temp;} else {mucalc[1] = 0;}

			if(npart[0]  != 0) {avenergy[0] = (sumE[0])/double(npart[0]);}
			if(npart[1]  != 0) {avenergy[1] = (sumE[1])/double(npart[1]);}

			if (sample.is_open())
			{
				sample  << left
					<< setw(10) << i_run			<< "\t"
					<< setw(10) << avenergy[0]		<< "\t"
					<< setw(10) << avenergy[1]		<< "\t"
					<< setw(10) << sumD[0]		<< "\t"
					<< setw(10) << sumD[1]		<< "\t"
					<< setw(10) << mucalc[0]		<< "\t"
					<< setw(10) << mucalc[1]		<< "\n";
//					<< swapmove << "   " << nmov 	<< "\n";
			}
			else{cerr 	<< "unable to open sample file for eqm. \n";}

			//prob_plot for analysis of results---------------------------
			if (probplot.is_open())
			{
			       probplot << left
					<< setw(10) << double(npart[0])/double(npart_tot)	<< "\t"
			       	<< setw(10) << V[0]/(V[0]+V[1])				<< "\t"
			       	<< setw(10) << double(npart[1])/double(npart_tot)	<< "\t"
					<< setw(10) << V[1]/(V[0]+V[1])				<< "\n";

			}
			else{cerr 	<< "unable to open probplot file for prod. \n";}
			//------------------------------------------------------------
	    	}
	    }
	}
	accept.close(); sample.close(); probplot.close();
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void production(int f_cyc, int f_mov)
{

	ofstream sample, probplot, rhodist;
	char filename[100];					//output file for sample 

	sprintf(filename, "run%d_%d_%.2f/sample_prod.dat", ind, scalefact, temp);
	sample.open(filename, ofstream::out | ofstream::app);
	sample << std::fixed;
	sample << std::setprecision(4);
	if (sample.is_open())
	{
		sample  << left
			<< setw(10) << "ncycle"		<< "\t"
			<< setw(10) << "toten1/N*eps"	<< "\t"
			<< setw(10) << "toten1/N*eps"	<< "\t"
			<< setw(10) << "avg_dens1"		<< "\t"
			<< setw(10) << "avg_dens2"		<< "\t"
			<< setw(10) << "mu_calc1"		<< "\t"
			<< setw(10) << "mu_calc2"		<< "\n";
	}
	else{cerr	<< "unable to open sample file. \n";}

	sprintf(filename, "run%d_%d_%.2f/probplot.dat", ind, scalefact, temp);
	probplot.open(filename, ofstream::out | ofstream::app);
	probplot << std::fixed;
	probplot << std::setprecision(4);

	sprintf(filename, "run%d_%d_%.2f/rhodist.dat", ind, scalefact, temp);
	rhodist.open(filename, ofstream::out | ofstream::app);
	rhodist << std::fixed;
	rhodist << std::setprecision(4);

	//rho_dist initialise--------------
	npr = 0;
	delrho = 1.2/double(maxbin);
	
	for(int i=0; i<maxbin; i++)
	{
		p[i] = 0.;
	}
	//---------------------------------

	nind = 0, count1 = 0, count2 = 0;
	int k, ran;
	double insE1[5], insE2[5], sumE[2], insD1[5], insD2[5], sumD[2], avenergy[2];
	aven[0]  = 0.; aven[1] = 0.; avdens[0] = 0.; avdens[1] = 0., avenergy[0] =0., avenergy[1]=0.;
	//mu[0] 	 = 0.;
	//mu[1] 	 = 0.;
	//swtry[0] = 0.;
	//swtry[1] = 0.;

	for(int icyc=f_cyc; icyc<procycle; icyc++)
	{
	    if (icyc%100 == 0) {cout << "progress(pro) " << icyc << "\t" << double(icyc)/double(procycle) << endl; }
	    for(int imov=f_mov; imov<nmov; imov++)
	    {	
		//cout << (eqcycle+icyc)*nmov+(imov+1) << endl;
		i_run = (icyc+eqcycle)*nmov+(imov+1);
		p_run = (icyc*nmov)+(imov+1); 
		ran   = round(ran2(idum)*(nmov));
		if      (ran < npart_tot)
		{
			//cout << i_run << "\t" << U1 << endl;	
			mcmove();
		}
		else if (ran < (npart_tot+volmov))
		{
			mcvol_npt();
		}
//		else
//		{
//			mcswap();
//		}

	    	if( p_run % ((procycle*nmov)/proframes) == 0){writexyz();} 		//writing a frame for visualization
	    	if( p_run % ((procycle*nmov)/proconfig) == 0){writeconfig(i_run,icyc,imov);}

	    	if( p_run > (nind*nsample_p - 5) )
	    	{
			if (npart[0] != 0)
			{
				k        = i_run%5;
				insE1[k] = U[0]; 
				insD1[k] = double(npart[0])/V[0];
			}
			else
			{
				k        = i_run%5;
				insE1[k] = 0.0; 
				insD1[k] = 0.0;
			}
			if (npart[1] != 0)
			{
				k        = i_run%5;
				insE2[k] = U[1]; 
				insD2[k] = double(npart[1])/V[1];
			}
			else
			{
				k        = i_run%5;
				insE2[k] = 0.0; 
				insD2[k] = 0.0;
			}
	    	}
	    	if( p_run % nsample_p == 0)
	    	{
			sumE[0]   = 0.; sumE[1]  = 0.;
			for (int j=0; j<5; j++) {sumE[0] += insE1[j]; sumE[1] += insE2[j];}
			sumE[0]  /= 5.; sumE[1] /= 5.;

			sumD[0]   = 0.; sumD[1]  = 0.;
			for (int j=0; j<5; j++) {sumD[0] += insD1[j]; sumD[1] += insD2[j];}
			sumD[0]  /= 5.; sumD[1] /= 5.;

			writeconfig_inst(i_run,icyc,imov);
			aven[0]   += sumE[0];
			aven[1]   += sumE[1];
			avdens[0] += sumD[0];
			avdens[1] += sumD[1];
			nind      += 1;

			if((swtry[0]  != 0)&&(mu[0]  != 0)) {mucalc[0] = (-log(mu[0]/double(swtry[0])))*temp;} else {mucalc[0] = 0;}
			if((swtry[1]  != 0)&&(mu[1]  != 0)) {mucalc[1] = (-log(mu[1]/double(swtry[1])))*temp;} else {mucalc[1] = 0;}

			if(npart[0]  != 0) {avenergy[0] = (sumE[0])/double(npart[0]);}
			if(npart[1]  != 0) {avenergy[1] = (sumE[1])/double(npart[1]);}
	
			if (sample.is_open())
			{
				sample  << left
					<< setw(10) << i_run			<< "\t"
					<< setw(10) << avenergy[0]		<< "\t"
					<< setw(10) << avenergy[1]		<< "\t"
					<< setw(10) << sumD[0]		<< "\t"
					<< setw(10) << sumD[1]		<< "\t"
					<< setw(10) << mucalc[0]		<< "\t"
					<< setw(10) << mucalc[1]		<< "\n";

			}
			else{cerr 	<< "unable to open sample file for prod. \n";}

			//prob_plot for analysis of results---------------------------
			if (probplot.is_open())
			{
				probplot << left
					  << setw(10) << double(npart[0])/double(npart_tot)	<< "\t"
					  << setw(10) << V[0]/(V[0]+V[1])				<< "\t"
					  << setw(10) << double(npart[1])/double(npart_tot)	<< "\t"
					  << setw(10) << V[1]/(V[0]+V[1])				<< "\n";

			}
			else{cerr 	<< "unable to open probplot file for prod. \n";}
			//------------------------------------------------------------

			//rho_dist sampling-------------------------------------------
			npr += 1;
			rho[0]	     = npart[0]/V[0];
			rho[1]	     = npart[1]/V[1];
			if(V[0]==0) rho[0] = 0;
			if(V[1]==0) rho[1] = 0;
			binindex    = int(rho[0]/delrho);
			p[binindex]+= 1;
			binindex    = int(rho[1]/delrho);
			p[binindex]+= 1;
			//rho_dist sampling-------------------------------------------
	    	}
	    }
	}
	sample.close();
	probplot.close();
	
	//rho_dist print-------------------------------------------
	double dist;
	if (rhodist.is_open())
	{
		for(int i=0; i<maxbin; i++)
		{
			dist = delrho*(double(i)+0.5);
			rhodist << left
				 << setw(10) << dist			<< "\t"
				 << setw(10) << p[i]/double(npr)	<< "\n";
		}
	}
	else{cerr 	<< "unable to open rhodist file for prod. \n";}
	rhodist.close();
	
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double totalenergy(int boxind)
{
	double r1;
	double toten = 0., psi_ij, psi_ji;;
	vector rij;

//	W[boxind] = 0.;
	for (int i=0; i<npart_tot-1; i++)
	{
	    if (ibox[i] == boxind)
	    {
		for (int j=i+1; j<npart_tot; j++)
		{
		    if (ibox[j] == boxind)
		    {
			rij = colloid[i] - colloid[j];
			rij.pbc(box[boxind]);		//Minimum image conversion
//			r2  = vectinnerpdct(rij,rij);
			r1  = rij.norm();
			if (r1 < cutoff)
			{
				if (r1 < 2.0){
					toten = 1.0e+050;
					return toten;
				}
				else{
					psi_ij = pair_energy(boxind,colloid[i], patch[i], colloid[j],rp[j], rc[j]);
					psi_ji = pair_energy(boxind,colloid[j], patch[j], colloid[i],rp[i], rc[i]);
					toten += (psi_ij+psi_ji)/2.;
				}
			}
		    }
		}
	    }
	}
//	if (tailcorr)
//	{
//		rhoins	 = double(npart[boxind])/V[boxind];
//		toten	+= double(npart[boxind]*ucor(npart[boxind],boxind));
//	}
	//cout << toten << endl;
	return toten/temp;

}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double potenergy(int j, vector rn, vector pn, int boxind)
{
	double r1;
	double en = 0.0, psi_ij, psi_ji;
	vector rij;

	vector rpn, rcn;

	for (int i=0; i<npart_tot; i++)
	{
	    if (ibox[i] ==  boxind)
	    {
		if (i != j)
		{
			rij = colloid[i] - rn;
			rij.pbc(box[boxind]);		//Minimum image conversion
			r1  = rij.norm();
			if (r1 < cutoff)
			{
				if (r1 < 2.0){
					en = 1.0e+050;
					return en;
				}
				else{
					rpn   = rn +  pn.vectimes(ap);
					rcn   = rn + -pn.vectimes(ac);
					psi_ij = pair_energy(boxind,colloid[i], patch[i], rn, rpn, rcn);
					psi_ji = pair_energy(boxind,rn, pn, colloid[i],rp[i], rc[i]);
					en    += (psi_ij+psi_ji)/2.;
				}
			}
		}
	    }
	}
	return en;

}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double pair_energy(int boxind, vector ri, vector pi, vector rj, vector rjp, vector rjc)
{
	vector rij, rijp, rijc;
	double etap, etac, x;
	double thetai, thetaic, thetaip;
	double sump, sumc, r, kr_term;
	double phi, psi=0.;
	double termc, termp;

	//distance calculation
	rij  = rj-ri;		rij.pbc(box[boxind]);
	rijp = rjp - ri;	rijp.pbc(box[boxind]);
	rijc = rjc - ri;	rijc.pbc(box[boxind]);

	etap = rijp.norm()/rij.norm();
	etac = rijc.norm()/rij.norm();

	//angle calculation
	thetai  = acos(rij.dot(pi)/rij.norm());
	thetaip = acos(rijp.dot(pi)/rijp.norm());
	thetaic = acos(rijc.dot(pi)/rijc.norm());

	//potential of i at point p of j
	sumc    = 0., sump = 0.;
	r       = rij.norm();
	//kr_term = 1./(kappa*sqrt(r));
	term1   = exp(-kappa*r)/r;
	x       = cos(thetaip);
	#pragma omp parallel reduction (+:sumc,sump)
	{
		#pragma omp for
		for(int l=0; l<=maxit; l++){
			sumc += ac_pow[l]*legendre(l,x);
			sump += ap_pow[l]*legendre(l,x);
		}
	}
	termp = (zc*zp*sumc + zp*zp*sump)*exp(-kappa*r*(etap-1.))/etap;
	//phi = (constant*kr_term*(zc*sumc + zp*sump))+((zc+zp)*constant*term1*term2);
	//cout << phi << endl;

	//energy
	//psi += zp*term2*phi;

	//potential of i at point c of j
	sumc    = 0., sump = 0.;
	//r       = rijc.norm();
	//kr_term = 1./(kappa*sqrt(r));
	//term1   = exp(-kappa*r)/r;
	x       = cos(thetaic);
	#pragma omp parallel reduction (+:sumc,sump)
	{
		#pragma omp for
		for(int l=0; l<=maxit; l++){
			sumc += ac_pow[l]*legendre(l,x);
			sump += ap_pow[l]*legendre(l,x);
		}
	}
	termc = (zc*zc*sumc + zp*zc*sump)*exp(-kappa*r*(etac-1.))/etac;
	//phi = (constant*kr_term*(zc*sumc + zp*sump))+((zc+zp)*constant*term1*term2);

	//energy
	//psi += zc*term2*phi;
	//psi = constant*term1*term2*term2*(termc+termp);
	psi = const_term2*term1*(termc+termp);

	return psi;
}

//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
//double ucor(int N, int boxind)
//{
//	double rhoins, invcut9, UTC;
//	invcut9 = invcut3*invcut3*invcut3;
//	rhoins  = double(N)/V[boxind];
//	UTC	= (8./3.)*M_PI*rhoins*((invcut9/3.) - invcut3);
//	return UTC;
//}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
//double pressurecalc(int boxind)
//{
//	double rhoins, PTC, invcut9;
//	rhoins	= double(npart[boxind])/V[boxind];
//	invcut9 = invcut3*invcut3*invcut3;
//	simuP[boxind]	= (rhoins*temp) + (W[boxind]/(3.*V[boxind]));
//	if (tailcorr) 
//	{
//		simuP[boxind] += (16./3.)*M_PI*rhoins*rhoins*((2.*invcut9/3.) - invcut3);
//	}
//	return simuP[boxind];
//}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
//double potenergychk(int j, vector rn, vector pn, int boxind)
//{
//	double en = 0.0;
//	double r2, r1, angi, angj;
//	vector rij, rdum, rdum1;
//	

//	for (int i=0; i<npart_tot; i++)
//	{
//	    if (ibox[i] ==  boxind)
//	    {
//		if (i != j)
//		{
//			vectsub(colloid[i],rn,rij);
//			rdum = rij;

//			rij.pbc(box[boxind]);		//Minimum image conversion
//			rdum1 = rij;

//			r2  = vectinnerpdct(rij,rij);
//			r1  = sqrt(r2);

//			if (r1 < cutoff)
//			{
//				if (r1 <= 1.0)
//				{
//					cout << std::fixed;
//					cout << std::setprecision(15);
//			cout << rdum.x << endl;
//			cout << rdum.y << endl;
//			//cout << rdum.z << endl;

//			cout << rdum1.x << endl;
//			cout << rdum1.y << endl;
//			//cout << rdum1.z << endl;
//			cout << i << " " << j << " " << r2 << " " << r1  << endl;
//			cout << ibox[i] << " " << ibox[j] << " " << boxind  << endl;
//			cout << colloid[i].x << " "<< colloid[i].y << " \n";//<< colloid[i].z << endl;
//			cout << colloid[j].x << " "<< colloid[j].y << " \n";//<< colloid[j].z << endl;
//					cout << "Press enter to continue... \n";
//					cin.ignore();
//					return en;
//				}
//				else if (r1 <= cutoff)
//				{
//					angi = vectinnerpdct(rij,patch[i]);
//					angi = acos(angi/r1);
//					angj = vectinnerpdct(-rij,pn);
//					angj = acos(angj/r1);
//					if      ((angi <= delta) && (angj <= delta)) { en += epsilon;}	//pp rep
//					else if ((angi >= delta) && (angj >= delta)) { en += epsilon;}	//npnp rep
//					else if ((angi <= delta) && (angj >= delta)) { en += epsilon2;}	//p-np att
//					else if ((angi >= delta) && (angj <= delta)) { en += epsilon2;}	//np-p att		
//				}
//			}
//		}
//	    }
//	}
//	return en;

//}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void mcmove()
{
	int	ipart, boxind;
	double	Uo, Un, arg, rand;
	vector dr, ro, rn, po, pn;

	Uo	= 0.;
	Un	= 0.;

	ipart    = int(ran2(idum)*npart_tot); 	//Choose a particle randomly
	ro       = colloid[ipart];			// Assign old configuration to a variable
	po       = patch[ipart];			// Assign old orientation to a variable
	boxind   = ibox[ipart];
	Uo       = potenergy(ipart, ro, po, boxind);
//	Uo       = potenergy(ipart, ro, boxind);
//	if (Uo > 1000) {cout << "move \n"; potenergychk(ipart, ro, po, boxind);}

	movetry[boxind] += 1;

//	dr.x     = (ran2(idum)-0.5)*drmax[boxind]; 		//Choose a random displacement in each direction
//	dr.y     = (ran2(idum)-0.5)*drmax[boxind];	
//	dr.z     = (ran2(idum)-0.5)*drmax[boxind];
	dr.setrandom(drmax[boxind]);

//	vectadd(ro, dr ,rn) 				// rn = New configuration
	rn       = ro + dr;
//	pn       = newpatch(po, 0.1745329);	// pn = New orientation// 10deg
	pn       = newpatch(po, 0.05);	// pn = New orientation// 0.1rad
	Un       = potenergy(ipart, rn, pn, boxind);
//	Un       = potenergy(ipart, rn, boxind);

	arg      = exp(-(Un-Uo)/temp);
	rand     = ran2(idum);
	if(rand < arg)			
	{
		//accepted
		moveaccept[boxind]+= 1;
		rn.pbc(box[boxind]);		//Periodic boundary condition
		colloid[ipart]     = rn;
		patch[ipart]       = pn;
		rp[ipart]          = rn +  pn.vectimes(ap);
		rc[ipart]          = rn + -pn.vectimes(ac);
		U[boxind]	    += (Un-Uo)/temp;
		if(Uo> 1000) cout << "move accepted\t" << Uo << " " <<Un << " "<< U[boxind] << endl; 
	}
	return;

}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void mcvol_gemc()
{
	double U1o, U1n,U2o, U2n, V1o, V1n, V2o, V2n, lnvn, V_tot;
	double box1o, box1n, box2o, box2n;
	double fact1, fact2, arg1, arg2, rand;

	voltry	+= 1;
	U1o	 = U[0];
	U2o	 = U[1];
	V1o	 = V[0];
	V2o	 = V[1];
	V_tot	 = V[0] + V[1];
	lnvn	 = log(V1o/V2o) + (ran2(idum)-0.5)*dvmax;		//random walk in ln(v1/v2)
	V1n	 = V_tot*exp(lnvn)/(1.+exp(lnvn));
	V2n	 = V_tot - V1n;
	box1o	 = box[0];
	box2o	 = box[1];
	//box[0]  = pow(V1n,(1./3.));
	//box[1]  = pow(V2n,(1./3.));
	box[0]  = pow(V1n,(1./2.));
	box[1]  = pow(V2n,(1./2.));

	fact1 = box[0]/box1o;
	fact2 = box[1]/box2o;
	#pragma omp parallel
	{
		#pragma omp for
		for(int  i=0; i<npart_tot; i++)
		{
			 if(ibox[i] == 0) {colloid[i].vectimes(fact1);}		//vectimes(fact1,colloid[i],colloid[i]);
			 if(ibox[i] == 1) {colloid[i].vectimes(fact2);}		//vectimes(fact2,colloid[i],colloid[i]);
		}
		#pragma omp for
		for(int i=0;i<npart_tot;i++){
			rp[i]   = colloid[i] +  patch[i].vectimes(ap);
			rc[i]   = colloid[i] + -patch[i].vectimes(ac);
		}
	}
	U1n	 = totalenergy(0);
	U2n	 = totalenergy(1);

//	arg1	 = exp(-(1./temp)*((U1n-U1o)-(double(npart[0]+1)*temp*log(V1n/V1o))));
//	arg2	 = exp(-(1./temp)*((U2n-U2o)-(double(npart[1]+1)*temp*log(V2n/V2o))));
	arg1	 = pow((V1n/V1o),double(npart[0]+1)) * exp(-(U1n-U1o)/temp);
	arg2	 = pow((V2n/V2o),double(npart[1]+1)) * exp(-(U2n-U2o)/temp);

	rand	 = ran2(idum);
	if (rand < (arg1*arg2))
	{
		 volaccept +=1;
		 U[0]       = U1n;
		 U[1]       = U2n;
		 V[0]       = V1n;
		 V[1]       = V2n;
		 rho[0]     = npart[0]/V[0];
		 rho[1]     = npart[1]/V[1];
		if(V[0]==0) rho[0] = 0;
		if(V[1]==0) rho[1] = 0;
	} else {
		fact1 = box1o/box[0];
		fact2 = box2o/box[1];
		#pragma omp parallel
		{
			#pragma omp for
			for(int  i=0; i<npart_tot; i++)
			{
				 if(ibox[i] == 0) {colloid[i].vectimes(fact1);}
				 if(ibox[i] == 1) {colloid[i].vectimes(fact2);}
			}
			#pragma omp for
			for(int i=0;i<npart_tot;i++){
				rp[i]   = colloid[i] +  patch[i].vectimes(ap);
				rc[i]   = colloid[i] + -patch[i].vectimes(ac);
			}
		}
		box[0]	     = box1o;
		box[1]	     = box2o;
	}
	return;
	
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void mcvol_npt()
{
	double Uo, Un, Vo, Vn, lnvn, boxo, boxn;
	double fact, arg, rand;

	voltry += 1;
	Uo      = U[0];
	Vo      = V[0];
	lnvn    = log(Vo) + (ran2(idum)-0.5)*dvmax;		//perform a random walk in ln V
	Vn      = exp(lnvn);
//	boxn    = pow(Vn, 1./3.);
	boxn    = pow(Vn, 1./2.);
	boxo    = box[0];
	fact    = boxn/boxo;
	box[0]  = boxn;

	#pragma omp parallel
	{
		#pragma omp for
		for (int i=0; i<npart[0]; i++)
		{
			colloid[i].vectimes(fact);
			//colloid[i].x *= fact;
			//colloid[i].y *= fact;
			//colloid[i].z *= fact;
		}
		#pragma omp for
		for(int i=0;i<npart[0];i++){
			rp[i]   = colloid[i] +  patch[i].vectimes(ap);
			rc[i]   = colloid[i] + -patch[i].vectimes(ac);
		}
	}

	Un      = totalenergy(0);
	arg     = -((Un-Uo) + (pressure*(Vn-Vo)) - (double(npart[0]+1)*log(Vn/Vo)*temp)) / temp;
	rand    = ran2(idum);
	if ( rand < exp(arg) )
	{
		//Accepted - update totalenergy, volume and density
		volaccept += 1;
		U[0] 	   = Un;
		//V[0]	   = pow(box[0],3.);
		V[0]	   = pow(box[0],2.);
		rho[0]	   = npart[0]/V[0];
	}
	else
	{
		//Rejected - restore old positions and box size
		fact = 	boxo/boxn;
	#pragma omp parallel
	{
		#pragma omp for
		for (int i=0; i<npart[0]; i++)
		{
			colloid[i].vectimes(fact);
			//colloid[i].x *= fact;
			//colloid[i].y *= fact;
			//colloid[i].z *= fact;
		}
		#pragma omp for
		for(int i=0;i<npart[0];i++){
			rp[i]   = colloid[i] +  patch[i].vectimes(ap);
			rc[i]   = colloid[i] + -patch[i].vectimes(ac);
		}
	}
		box[0] = boxo;
	}
	return;
	
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void mcswap()
{
	int	ipart, in, out;
	double	Uo, Un, arg, rand;
	vector ro, rn, po, pn;
	double del_e, dtail_in, dtail_out;

	swaptry += 1;
	rand   = ran2(idum);
	if (rand < 0.5)
	{
	       in  = 0;
	       out = 1;
	}
	else
	{
	       in  = 1;
	       out = 0;
	}

 	if (npart[out] == 0) 	//condition to add particle from empty box to full box
	{
		ipart = npart_tot+1; ibox[ipart] = out;
	}
	else
	{
		do
		{
		  ipart = int(ran2(idum)*npart_tot);
		} while (ibox[ipart] != out);
	}

//	rn.x 	   = (ran2(idum)-0.5)*box[in];			//New position in inbox
//	rn.y 	   = (ran2(idum)-0.5)*box[in];
//	rn.z 	   = (ran2(idum)-0.5)*box[in];
	rn.setrandom(box[in]);
	pn	   = vector(1.,0.,0.);				// pn = New orientation
	pn	   = newpatch(pn,0.7854);
	Un	   = potenergy(ipart, rn, pn, in);
//	Un	   = potenergy(ipart, rn, in);

	mu[in]    += V[in]*exp(-Un/temp)/double(npart[in]+1);	//update chemical potential
	swtry[in] += 1;

	if (npart[out] == 0) {return;}

	ro	   = colloid[ipart];
	po	   = patch[ipart];
	Uo	   = potenergy(ipart, ro, po, ibox[ipart]);
//	Uo	   = potenergy(ipart, ro, ibox[ipart]);
//	if (Uo > 1000) {cout << "swap \n"; potenergychk(ipart, ro, po, ibox[ipart]);}
 
	arg	   = exp( -( (Un-Uo) + (temp*log( (V[out]*double(npart[in]+1))/(V[in]*double(npart[out])) )) )/temp );

//	del_e     = (Un-Uo) + temp*log( (V[out]*double(npart[in]+1)) / (V[in]*double(npart[out])) );

//	if(tailcorr){
//		dtail_in  = (npart[in]+1)*ucor(npart[in]+1,in) - npart[in]*ucor(npart[in], in);
//		dtail_out = (npart[out]-1)*ucor(npart[out]-1, out) - npart[out]*ucor(npart[out], out);
//		del_e += dtail_in + dtail_out;
//	}

//	arg = exp(-del_e/temp);
	rand	   = ran2(idum);
	if (rand < arg)
	{
		//accepted - add new particle to in box
		swapaccept	+= 1;
		colloid[ipart] = rn;
		patch[ipart]	 = pn;
		rp[ipart]      = rn +  pn.vectimes(ap);
		rc[ipart]      = rn + -pn.vectimes(ac);
		ibox[ipart]	 = in;
		U[in]		+= Un;
		U[out]		-= Uo;
		npart[in]	+= 1;
		npart[out]	-= 1;
		rho[0]	        = npart[0]/V[0];
		rho[1]	        = npart[1]/V[1];
		if(V[0]==0) rho[0] = 0;
		if(V[1]==0) rho[1] = 0;

//		if(tailcorr){
//			totalenergy(in);
//			totalenergy(out);
////			U[in] += dtail_in;
////			U[out] = U[out] + dtail_out;
//		}
	}
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//			Minimum image conversion is same as pbc	
//--------------------------------------------------------------------------------------------------------------------------
vector mic(vector rn, int boxind)
{
	double boxi;

	boxi = box[boxind];
	rn.x = rn.x - boxi*round(rn.x/boxi);
	rn.y = rn.y - boxi*round(rn.y/boxi);
	//rn.z = rn.z - boxi*round(rn.z/boxi);
    	return rn;	
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
vector newpatch(vector en, double maxrad)
{
	vector e;
	e = en;
	double R[3][3];
	rotation(R, maxrad);
	en.x	= R[0][0]*e.x + R[0][1]*e.y + R[0][2]*e.z;
	en.y	= R[1][0]*e.x + R[1][1]*e.y + R[1][2]*e.z;
	en.z	= R[2][0]*e.x + R[2][1]*e.y + R[2][2]*e.z;
	return en;	
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
//vector newpatch2d(vector en)
//{
//	vector e, unitv, t;
//	e = en;
//	double Q[2], unit_vec, t_vec,t_vec1, dtheta = 0.5;	//unit vector generated and resultant vector
//	//[rand(0)*(y-x)]+x  --> TO GET A RANDOM NUMBER IN [x,y]
//	double check = 1.;
//	do
//	{
//	  Q[0] 	  = ran2(idum)*(2.*M_PI);			//random angle from 0-2*pi
//	  Q[1] 	  = (180./M_PI)*Q[0];
//	  unitv.x = cos(Q[0]);	
//	  unitv.y = sin(Q[0]);
//	  check   = acos((unitv.x*e.x) + (unitv.y*e.y));
//	} while (check > 0.1745329);			//10 deg tolerance
//							//(check .gt. 1.570796d0)90 deg tolerance
//	unit_vec  = vectinnerpdct(unitv,unitv);
//	
//	t.x	  = (dtheta*unitv.x) + e.x;
//	t.y	  = (dtheta*unitv.y) + e.y;
//	t_vec	  = vectinnerpdct(t,t);
//	t_vec1	  = sqrt(t_vec);
//	t.x 	  = t.x/t_vec1;
//	t.y 	  = t.y/t_vec1;
//	en 	  = t;
//	return en;
//}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double ran2(long *idum) // ran2 algorithm form Numecial Recipes 2nd Edition
{
	int 	    j;
	long 	    k;
	static long idum2 = 123456789;
	static long iy    = 0.0;
	static long iv[NTAB];
	double 	    temp;

	if (*idum <= 0)				//Initialize
	{
	   if   (-(*idum) < 1) {*idum = 1;}	//Prevent idum = 0
	   else                {*idum = -(*idum);}
	   idum2  = (*idum);
	   
	   for (j = NTAB+7 ; j>=0 ; j--)
	   {
		k     = (*idum)/IQ1;
		*idum = IA1*(*idum-k*IQ1)-k*IR1;
		if (*idum < 0) {*idum += IM1;}
		if (j < NTAB)  {iv[j]  = *idum;}
	   }
	   iy = iv[0];
	}

	k     = (*idum)/IQ1;			//Start here when not initializing
	*idum = IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) {*idum +=IM1;}
	k     = idum2/IQ2;
	idum2 = IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) {idum2 += IM2;}
	j     = iy/NDIV;
	iy    = iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp = AM*iy) > RNMX) {return RNMX;}
	else {return temp;}

}
//--------------------------------------------------------------------------------------------------------------------------
//	Subroutine to generate a random vector uniformly on a 4D unit sphere based on Marsaglia algorithm modified
//	by Vesley (Journal of computation Physics 47,291-296 (1982)). (Ref: UMS, Frenkel and Smit, pg-578.)
//--------------------------------------------------------------------------------------------------------------------------
double quaternion(double Q[4])  
{
	double S1 = 2.0, S2 =2.0;
	double rand1, rand2, rand3, rand4, ranh;

	do
	{
	  rand1 = 1.0 - (2.0*ran2(idum));
	  rand2 = 1.0 - (2.0*ran2(idum));
	  S1    = (rand1*rand1) + (rand2*rand2);
	} while (S1 >= 1.0 );
	do
	{
	  rand3 = 1.0 - (2.0*ran2(idum));
	  rand4 = 1.0 - (2.0*ran2(idum));
	  S2    = (rand3*rand3) + (rand4*rand4);
	} while (S2 >= 1.0 );

	ranh    = sqrt((1.0-S1)/S2);
	Q[0]	= rand1;
	Q[1]	= rand2;
	Q[2]	= rand3*ranh;
	Q[3]	= rand4*ranh;	

	return Q[4];
//	note that :(Q[0]*Q[0])+(Q[1]*Q[1])+(Q[2]*Q[2])+(Q[3]*Q[3]) = 1  
//	such a quaternion can be thought of as a unit-vector in 4D space
}
//--------------------------------------------------------------------------------------------------------------------------
//	Quaternion method to find new orientation in 3d (Ref: UMS, Frenkel and Smit, pg:49) 
//	Read with newpatch subroutine in con.h file for full algorithm
//--------------------------------------------------------------------------------------------------------------------------
double rotation(double R[3][3], double maxrad)
{
	//double R[3][3]; //Rotation matrix
	double check = 1.;
	do
	{
	  quaternion(Q);
	  check = acos(Q[0]);
	} while ( check > maxrad);
	  //while (check > 0.05); // 0.1 rad tolerance, alpha = 2*acos(q0), check if acos(q0) > alpha/2.
	  //while (check > 0.7854); // (pi/2)rad tolerance, alpha = 2*acos(q0), check if acos(q0) > alpha/2.
	
	double a = (Q[0]*Q[0])+(Q[1]*Q[1])+(Q[2]*Q[2])+(Q[3]*Q[3]);
	
	R[0][0] = Q[0]*Q[0] + Q[1]*Q[1] - Q[2]*Q[2] - Q[3]*Q[3];
	R[0][1] = 2.*(Q[1]*Q[2] - Q[0]*Q[3]);
	R[0][2] = 2.*(Q[1]*Q[3] + Q[0]*Q[2]);

	R[1][0] = 2.*(Q[1]*Q[2] + Q[0]*Q[3]);
	R[1][1] = Q[0]*Q[0] - Q[1]*Q[1] + Q[2]*Q[2] - Q[3]*Q[3];
	R[1][2] = 2.*(Q[2]*Q[3] - Q[0]*Q[1]);

	R[2][0] = 2.*(Q[1]*Q[3] - Q[0]*Q[2]);
	R[2][1] = 2.*(Q[2]*Q[3] + Q[0]*Q[1]);
	R[2][2] = Q[0]*Q[0] - Q[1]*Q[1] - Q[2]*Q[2] + Q[3]*Q[3];	 

	return R[3][3];
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void writeconfig(long int i_run, int icyc, int imov)
{
	bool eqm;	
	char filename[100], run[5];
	if (i_run <= (nmov*eqcycle)) {sprintf(run, "eq");	eqm = true;}
	else		               {sprintf(run, "pro");	eqm = false;}
	ofstream config;

	sprintf(filename, "run%d_%d_%.2f/%s/config_%d.dat", ind, scalefact, temp, run, count2 );
	config.open(filename, ofstream::out | ofstream::app);
  	config << std::fixed;
  	config << std::setprecision(4);
	if (config.is_open())
	{
	   config << i_run  << "\t" << icyc   << "\t" << imov << "\t" << npart[0] << "\t" << npart[1] << "\t" 
		   << box[0] << "\t" << box[1] << endl;
	   for (int j=0; j<npart_tot; j++)
	   {
		config << left
			<< setw(10) << colloid[j].x	<< "   "
			<< setw(10) << colloid[j].y	<< "   "
			<< setw(10) << colloid[j].z	<< "   " 
			<< setw(10) << patch[j].x	<< "   "
			<< setw(10) << patch[j].y	<< "   "
			<< setw(10) << patch[j].z	<< "   "
			<< setw(10) << ibox[j]	<< "\n ";
	   }
	   config   << endl;
	   count1 += 1;
	   if (count1%100 == 0) {count2++;}
	}
	else {cerr << "unable to open config_" << run << " file for output. \n";}
	config.close();
	
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void writeconfig_inst(long int i_run, int icyc, int imov)
{
	bool eqm;
	if (i_run <= (nmov*eqcycle)) {eqm = true;}
	else		               {eqm = false;}
	char filename[100];
	ofstream instconfig;

	sprintf(filename, "run%d_%d_%.2f/config_inst%d.dat", ind, scalefact, temp, ind );
	instconfig.open(filename, ofstream::out | ofstream::trunc);
  	instconfig << std::fixed;
  	instconfig << std::setprecision(10);
	if (instconfig.is_open())
	{
	   instconfig << i_run  << "\t" << icyc   << "\t" << imov     << "\t" << eqm      << "\t" << npart[0] << "\t" << npart[1] << "\t" 
		       << box[0] << "\t" << box[1] << "\t" << drmax[0] << "\t" << drmax[1] << "\t" << dvmax    << "\t" 
			<< U[0]   << endl;
	   for (int j=0; j<npart_tot; j++)
	   {
	    instconfig  << left
			  << setw(15) << colloid[j].x	<< "   "
			  << setw(15) << colloid[j].y	<< "   "
			  << setw(15) << colloid[j].z	<< "   " 
			  << setw(15) << patch[j].x		<< "   "
			  << setw(15) << patch[j].y		<< "   "
			  << setw(15) << patch[j].z		<< "   "
			  << setw(15) << ibox[j]		<< "\n ";
	   }

	}
	else {cerr << "unable to open config_inst file for output. \n";}
	instconfig.close();	
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void writexyz()
{
	char filename[100];
	ofstream config;
	sprintf(filename, "run%d_%d_%.2f/config.xyz", ind, scalefact, temp );
	config.open(filename, ofstream::out | ofstream::app);
  	config << std::fixed;
  	config << std::setprecision(4);
	if (config.is_open())
	{
	   config << npart_tot*2 + 2 << endl;
	   config 	       	     << endl;
	   config <<  "C" << "  "    << box[0]/2. 		  << "  " << box[0]/2. << "  " << box[0]/2. << endl;
	   config <<  "C" << "  "    << box[0] + 10.0 + box[1]/2. << "  " << box[1]/2. << "  " << box[1]/2. << endl;
	   for (int j=0; j<npart_tot; j++)
	   {
		if      (ibox[j] == 0)
		{
			config << "He1"						<< "  "
				<< colloid[j].x					<< "  "
				<< colloid[j].y					<< "  "
				<< colloid[j].z					<< "\n";
			config << "He2"						<< "  "
				<< colloid[j].x + patch[j].x/20.			<< "  "
				<< colloid[j].y + patch[j].y/20.			<< "  "
				<< colloid[j].z + patch[j].z/20.			<< "\n";
		}
		else if (ibox[j] == 1)
		{
			config << "He3"						<< "  "
				<< colloid[j].x + box[0] + 10.0			<< "  "
				<< colloid[j].y					<< "  "
				<< colloid[j].z					<< "\n";
			config << "He4"						<< "  "
				<< colloid[j].x + box[0] + 10.0 + patch[j].x/20.	<< "  "
				<< colloid[j].y + patch[j].y/20.			<< "  "
				<< colloid[j].z + patch[j].z/20.			<< "\n";
		}
	   }
	}
	else {cerr << "unable to open xyz file for output. \n";}
	config.close();
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void finalparameters()
{
	char parafile[100];
	ofstream outfile;
	sprintf(parafile, "run%d_%d_%.2f/para.%d.%d.out", ind, scalefact, temp, ind, scalefact);
	outfile.open(parafile, ofstream::out | ofstream::trunc);
	if (outfile.is_open())
	{
	  outfile << ind		<< "\t\t index"	<< "\n"
		  << npart[0]		<< "\t\t npart1" 	<< "\n"
		  << npart[1]		<< "\t\t npart2" 	<< "\n"
		  << chi		<< "\t\t chi"	<< "\n"
		  << rho[0]		<< "\t\t rho1"	<< "\n"
		  << rho[1]		<< "\t\t rho2"	<< "\n"
		  << box[0]		<< "\t\t box1"	<< "\n"
		  << box[1]		<< "\t\t box2"	<< "\n"
		  << V[0]		<< "\t\t volume1"	<< "\n"
		  << V[1]		<< "\t\t volume2"	<< "\n"
		  << kappa		<< "\t\t kappa"	<< "\n"
		  << cutoff 		<< "\t\t cutoff"	<< "\n"
		  << temp 		<< "\t\t temp"	<< "\n"
		  << pressure		<< "\t\t pressure"	<< "\n"
		  << drmax[0]		<< "\t\t drmax1" 	<< "\n"
		  << drmax[1]		<< "\t\t drmax2" 	<< "\n"
		  << dvmax		<< "\t\t dvmax" 	<< "\n"
		  << swapmove		<< "\t\t swapmove"	<< "\n"
		  << eqcycle		<< "\t\t eqcycle"	<< "\n"
		  << procycle		<< "\t\t procycle"	<< "\n"
		  << runcontinue	<< "\t\t runcontinue"<< "\n"
		  << initial		<< "\t\t initial"	<< "\n"
		  << runrestart	<< "\t\t runrestart"	<< "\n"
		  << eqframes 	<< "\t\t eqframes"   << "\n"
		  << proframes	<< "\t\t proframes"  << "\n"
		  << maxit		<< "\t\t maxit"  << "\n"
		  << scalefact	<< "\t\t scalefact"  << "\n"
		  << "\n" "Change runrestart to true for restart run. \n";
	}
	else{cerr << "Unable to open file for final parameter output. \n";}
	outfile.close();

	ofstream outfile1;
	sprintf(parafile, "run%d_%d_%.2f/para.%d.%d.inp", ind, scalefact, temp, ind, scalefact+1);
	outfile1.open(parafile, ofstream::out | ofstream::trunc);
	if (outfile1.is_open())
	{
	 outfile1 << npart[0]	<< "\t\t npart1" 	<< "\n"
		  << npart[1]		<< "\t\t npart2" 	<< "\n"
		  << chi		<< "\t\t chi"		<< "\n"
		  << rho[0]		<< "\t\t rho1"	<< "\n"
		  << rho[1]		<< "\t\t rho2"	<< "\n"
		  << kappa 		<< "\t\t kappa"	<< "\n"
		  << temp 		<< "\t\t temp"	<< "\n"
		  << pressure		<< "\t\t pressure"	<< "\n"
		  << drmax[0]		<< "\t\t drmax1" 	<< "\n"
		  << drmax[1]		<< "\t\t drmax2" 	<< "\n"
		  << dvmax		<< "\t\t dvmax" 	<< "\n"
		  << eqcycle		<< "\t\t eqcycle"	<< "\n"
		  << procycle		<< "\t\t procycle"	<< "\n"
		  << runcontinue	<< "\t\t runcontinue"<< "\n"
		  << "0"		<< "\t\t initial"	<< "\n"
		  << "1"		<< "\t\t runrestart"	<< "\n"
		  << eqframes 	<< "\t\t eqframes"   << "\n"
		  << proframes	<< "\t\t proframes"  << "\n"
		  << "0"		<< "\t\t start1"	<< "\n"
		  << "0"		<< "\t\t start2"	<< "\n"
		  << maxit		<< "\t\t maxit"	<< "\n"
		  << scalefact+1	<< "\t\t scalefact"	<< "\n";
	}
	else{cerr << "Unable to open file for final parameter output. \n";}
	outfile1.close();

	return;
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
double bessel(int n, double x)
{
	double bessel = 0., sum = 0.;
	for(int k=0; k<=n; k++){
		sum += n_plus_half_k(n,k)*pow((2*x),-k);
	}
	bessel = sqrt(M_PI/(2.*x))*exp(-x)*sum;
	return bessel;
}
//----------------------------------------------------------------------
double factorial(int x)
{
	double fact=1.;
	if(x==0) {
		 return fact;
	}
	else if (x>0){
		for(int i=1; i<=x; i++){
			fact = fact*double(i);
		}
	}
	return fact;
}
//----------------------------------------------------------------------
double n_plus_half_k(int n, int k)
{
	double value=0;
	value = factorial(n+k)/(factorial(k)*factorial(n-k));
	return value;
}
//----------------------------------------------------------------------
double legendre(int l, double x)
{
	double Pl[l+1];
	Pl[0] = 1.;
	Pl[1] = x;
	for(int i=1; i<=l; i++){
		Pl[i+1] = ((double((2*i)+1)*x*Pl[i]) - (double(i)*Pl[i-1]))/double(i+1);
//		cout << Pl[i] << " " << Pl[i+1] << " " << i << endl;
	}
	return Pl[l];
}
//----------------------------------------------------------------------
