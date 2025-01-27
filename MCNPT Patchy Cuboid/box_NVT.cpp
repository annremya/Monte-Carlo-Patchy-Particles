//--------------------------------------------------------------------------------------------------------------------------
//	MC NPT code written for a mixture of Inverse Patchy Colloids and isotropic colloids in 2d
//	Orientation dependent square-well potential (Kern-Frenkel model)
//	Modified to MCNVT and widom particle insertion method  and overlapping distribution method added
//	Based on algorithms from UMS (Frenkel and Smit) and under guidance of Dr.Ethayaraja Mani
//	Written by Remya Ann on 27/04/2018 for the third work 
//	Includes a parameter input file and configuartion header file				
//--------------------------------------------------------------------------------------------------------------------------

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string>
#include<iomanip>
using namespace std;
#include"GJK.cpp"
#include"conf.h"

double  Lx = 2., Ly = 1., Lz = 1.;
int	 nadj      	= 2500;
int	 nsample_e 	= 1000; 
int	 nsample_p 	= 1000; 
int 	 eqconfig	= 1000;
int	 proconfig	= 2000;
int	 widom_try	= 0;

double  Picut = 0.5;
double  XBcut = 0.5;
double  delta = 0.2, delta2; //width of patch

//Main program
//--------------------------------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])

{
	srand(time(0));
	if(argc > 1) { ind = atoi(argv[1]);}
	//ind = 1;
	
	readparameters();
	volmov = 0;
	nmov   = npart + volmov;	//Number of moves in each cycle
	count1 = 0;	
	count2 = 0;	
//	r3     = 1.0/(rcutlj*rcutlj*rcutlj);
//	sigB3  = sigma[1]*sigma[1]*sigma[1];
	delta2 = delta*delta;

	char folder[100];
	sprintf(folder,"mkdir run%d_%.2f", ind, epsilon_xb);
	system(folder);
	sprintf(folder,"mkdir run%d_%.2f/eq", ind, epsilon_xb);
	system(folder);
	sprintf(folder,"mkdir run%d_%.2f/pro", ind, epsilon_xb);
	system(folder);
	sprintf(folder,"mv parameter_check_%d_%.2f.dat run%d_%.2f/parameter_check_%d_%.2f.dat", ind,epsilon_xb, ind, epsilon_xb, ind, epsilon_xb);
	system(folder);

	long y    = -rand();					//seed for random number generator
	     idum = &y;


	//------------initial condition---------------------------------
	cout << initial << runrestart << runcontinue << endl;

	if      (initial) 
	{
		cout << "initial" << endl;
		initialcondition();
		U = totalenergy();
		cout << "Total energy = " << U/double(npart) << endl;
		equilibration(0,0);
		production(0,0);
		//widom_run();
	}
	else if (runrestart) 
	{
		cout << "runrestart" << endl;
		restart();
		U = totalenergy();
		cout << "Total energy = " << U/double(npart) << endl;
		equilibration(0,0);
		production(0,0);
		//widom_run();

	}
	else if (runcontinue)
	{
		cout << "runcontinue" << endl;
		restart();
		U = totalenergy();
		cout << "Total energy = " << U/double(npart) << endl;

		if(runstage == 1)
		{
			equilibration(f_cyc,f_mov);
			production(0,0);
			//widom_run();
		}
		else
		{
			production(f_cyc,f_mov);
		}
	}

	//writeconfig((eqcycle+procycle)*nmov);
	writexyz();
	finalparameters();

//	if (itest[0] != 0) {muex[0] = (-log(wtest[0]/double(itest[0])))*temp;} else {muex[0] = 0.;}
//	if (itest[1] != 0) {muex[1] = (-log(wtest[1]/double(itest[1])))*temp;} else {muex[1] = 0.;}
//	cout << "wtest = " << wtest << endl;
//	cout << "muex  = " << muex  << endl;

	ofstream log;
	char filename[100];
	sprintf(filename, "run%d_%.2f/log.dat", ind, epsilon_xb);
	log.open(filename, ofstream::out | ofstream::app);
	log << std::fixed;
	log << std::setprecision(4);
	if (log.is_open())
	{
		log << "Average Energy per particle   = " << aven/(double(npart)*double(nind))	<< "\n"
		    << "Running Energy per particle   = " << U/double(npart)				<< "\n"
		    << "Average density of the system = " << avdens/double(nind)			<< "\n"
		    << "Rho_A                         = " << double(nA)/V				<< "\n"
		    << "Rho_B                         = " << double(nB)/V				<< "\n\n";

//		    << "Density for Mu_ex calculation                    = " << wdens			<< "\n"
//		    << "Rho_A for Mu_ex calculation                      = " << wdensA		<< "\n"
//		    << "Rho_B for Mu_ex calculation                      = " << wdensB		<< "\n"
//		    << "Excess chemical potential of patchy particles    = " << muex[0]		<< "\n"
//		    << "Excess chemical potential of isotropic particles = " << muex[1]		<< "\n";
	}
	else{cerr   << "unable to open log file. \n";}
//	distribution(3,0,0,0,0.);
	free(colloid);
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
	  infile >> nA		>> dum
		  >> nB 		>> dum
		  >> volfrac  	>> dum
		  >> rcut 		>> dum
		  >> temp   		>> dum
		  >> pressure		>> dum
		  >> epsilon_xb	>> dum
		  >> epsilon_pi	>> dum
  		  >> epsilon_hb	>> dum
		  >> drmax  		>> dum
		  >> dvmax		>> dum
		  >> eqcycle		>> dum
		  >> procycle		>> dum
		  >> runcontinue	>> dum
		  >> initial		>> dum
		  >> runrestart	>> dum
		  >> eqframes   	>> dum
		  >> proframes  	>> dum
		  >> start		>> dum
		  >> ighost		>> dum
		  >> scalefact	>> dum;
	}
	else { cerr << "unable to open parameter file.\n";}
	infile.close();

	npart    = nA + nB;
	box      = cbrt(double(npart)*Lx*Ly*Lz/volfrac);
	V        = pow(box,3.);
	rho      = double(npart)/V;
	colloid  = (cube*) malloc(npart*sizeof(cube));
//	epsilon  = epsilon_xb*pow(1.05,scalefact);
//	cout << "Edge length of simulation box is: " << box1    << endl;
//	cout << "Value of epsilon is             : " << epsilon << endl;
	cout << "Edge length of simulation box is: " << box << endl;

	ofstream outfile;
	sprintf(parafile, "parameter_check_%d_%.2f.dat", ind,epsilon_xb);
	outfile.open(parafile, ofstream::out | ofstream::trunc);
	if (outfile.is_open())
	{
	 outfile << "index       = "	<< ind		<< "\n"
		  << "npart       = "	<< npart	<< "\n"
		  << "nA          = "	<< nA		<< "\n"
		  << "nB          = "	<< nB		<< "\n"
		  << "volfrac     = "	<< volfrac	<< "\n"
		  << "rho         = "	<< rho		<< "\n"
		  << "box         = "	<< box		<< "\n"
		  << "Volume      = "	<< V		<< "\n"
		  << "rcut        = "	<< rcut 	<< "\n"
		  << "temp        = "	<< temp	<< "\n"
		  << "pressure    = "	<< pressure	<< "\n"
		  << "epsilon_xb  = "	<< epsilon_xb	<< "\n"
		  << "epsilon_pi  = "	<< epsilon_pi	<< "\n"
		  << "epsilon_hb   = "	<< epsilon_hb	<< "\n"
		  << "drmax       = "	<< drmax	<< "\n"
		  << "dvmax       = "	<< dvmax	<< "\n"
		  << "eqcycle     = "	<< eqcycle	<< "\n"
		  << "procycle    = "	<< procycle	<< "\n"
		  << "runcontinue = "	<< runcontinue<< "\n"
		  << "initial     = "	<< initial	<< "\n"
		  << "runrestart  = "	<< runrestart	<< "\n"
		  << "eqframes    = "	<< eqframes 	<< "\n"
		  << "proframes   = "	<< proframes	<< "\n"
		  << "start       = "	<< start	<< "\n"
		  << "ighost      = "	<< ighost	<< "\n"
		  << "scalefact   = "   	<< scalefact  << "\n";
	}
	else{cerr << "Unable to open file for parameter check output. \n";}
	outfile.close();
	
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void initialcondition()
{
//	for(int i=0; i<npart; i++)
//	{
//		if(i<nA) {moltp[i] = 0;}	//patchy
//		else	  {moltp[i] = 1;}	//isotropic
//	}
	if (start == 1)					//Random initial configuration
	{
		cout << "Random initial Configuration \n";
		cube tmp;
		bool k1,k2;
		GJKalgorithm gjkvar; 
		colloid[0].setrandomcenter(box);
		colloid[0].setaxes();
		colloid[0].setvertices(Lx/2.,Ly/2.,Lz/2.);
		colloid[0].setpatches(Lx/2.,Ly/2.,Lz/2.);
		//colloid[0].center.printvec();
		for(int i=1; i<npart; i++){
			do{
				tmp.setrandomcenter(box);
				tmp.setaxes();
				tmp.rotatecube(abs(ran2(idum)*3.),(ran2(idum)-0.5)*M_PI, Lx/2.,Ly/2.,Lz/2.);
				//tmp.setvertices(Lx/2.,Ly/2.,Lz/2.);
				for(int j=0; j<i; j++){
					k1 = gjkvar.CollisionDetection(colloid[j],tmp);
					k2 = gjkvar.CollisionDetection(tmp,colloid[j]);
					if ((k1==1) || (k2==1)) break;
				}
			} while ((k1==1) || (k2==1));
			colloid[i] = tmp; 
		}
//		for(int i=0; i<npart; i++){
//			colloid[i].setpatches(Lx/2.,Ly/2.,Lz/2.);
//		}
	}

	else if (start == 2)				//Cubic 3D
	{
		bool test = 0;
		cout << "Simple Cubic Configuration \n";
		double AR  = Lx/Ly;
		double ncellx, ncelly, ncellz;
		double lcellx, lcelly, lcellz;
		ncellx = ceil(cbrt(double(npart)/(AR*AR)));
		ncelly = ceil(cbrt(double(npart)*AR));
		ncellz = ncelly;
		lcellx = box/ncellx; lcelly = box/ncelly; lcellz = box/ncellz;
		cout << "Number of cells = " << ncellx << " " << ncelly << " " << ncellz << endl;
		cout << "Unit cell lengths " << lcellx << " " << lcelly << " " << lcellz << endl;
		int index = 0;
		for(int iz=0; iz<ncellz; iz++){
			for(int iy=0; iy<ncelly; iy++){
				for(int ix=0; ix<ncellx; ix++){
					if(index>=npart) break;
					colloid[index].center.x = double(ix)*lcellx + lcellx/2.;
					colloid[index].center.y = double(iy)*lcelly + lcelly/2.;
					colloid[index].center.z = double(iz)*lcellz + lcellz/2.;
					index++;
				}
			}
		}
		double hbox = box/2.;
		for(int i=0; i<npart; i++){
			colloid[i].center.x -= hbox; colloid[i].center.y -= hbox; colloid[i].center.y -= hbox; 
		}
		for(int i=0; i<npart; i++){
			colloid[i].setaxes();
			colloid[i].setvertices(Lx/2.,Ly/2.,Lz/2.);
			colloid[i].setpatches(Lx/2.,Ly/2.,Lz/2.);
		}
		if(test){	//test for overlap in cubic initial configuration
			bool k;
			GJKalgorithm gjkvar; 
			for(int i=0; i<npart-1;i++){
				for(int j=i+1; j<npart; j++){
					k = gjkvar.CollisionDetection(colloid[i],colloid[j]);
					if (k == 1) {
						cout << i <<  " " << j << " overlap\n";
					}
					k = gjkvar.CollisionDetection(colloid[j],colloid[i]);
					if (k == 1) {
						cout << i <<  " " << j << " overlap\n";
					}
				}
			}
		}
	}

//	for(int i=0; i<npart; i++)			//Box index of particles
//	{
//		ibox[i] = 1;
//	}
	writexyz();
	writeconfig(0,0,0);
	return;

}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void restart()
{
	long int  nrun;
	char filename[100];
	int dum;
	ifstream restartconfig;

	sprintf(filename, "config_inst%d.dat", ind);
	restartconfig.open(filename, ifstream::in);
	if(restartconfig.is_open())
	{
	   restartconfig >> nrun >> f_cyc >> f_mov >> runstage >> box >> drmax >> dvmax;
	   for (int j=0; j<npart; j++)
	   {
		restartconfig >> dum
				>> colloid[j].center.x >> colloid[j].center.y >> colloid[j].center.z
				>> colloid[j].ei.x     >> colloid[j].ei.y     >> colloid[j].ei.z
				>> colloid[j].ej.x     >> colloid[j].ej.y     >> colloid[j].ej.z
				>> colloid[j].ek.x     >> colloid[j].ek.y     >> colloid[j].ek.z;
	   }
	}
	else{cerr << "unable to open config_inst file for input. \n";}
	restartconfig.close();

	for(int i=0; i<npart; i++){
		colloid[i].setvertices(Lx/2.,Ly/2.,Lz/2.);
		colloid[i].setpatches(Lx/2.,Ly/2.,Lz/2.);
	}

	box = cbrt(npart/rho);
	V  = pow(box,3.);
	cout << "Edge length of simulation box is: " << box << endl;

	bool test = 1;
	if(test){
		bool k;
		GJKalgorithm gjkvar; 
		for(int i=0; i<npart-1;i++){
			for(int j=i+1; j<npart; j++){
				k = gjkvar.CollisionDetection(colloid[i],colloid[j]);
				if (k == 1) {
					cout << i <<  " " << j << " overlap\n";
				}
				k = gjkvar.CollisionDetection(colloid[j],colloid[i]);
				if (k == 1) {
					cout << i <<  " " << j << " overlap\n";
				}
			}
		}
	}

	writeconfig(nrun,f_cyc,f_mov);
	writexyz();
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void equilibration(int f_cyc, int f_mov)
{

	ofstream sample, accept;
	char filename[100];					//output file for sample 
	sprintf(filename, "run%d_%.2f/sample.dat", ind, epsilon_xb);
	sample.open(filename, ofstream::out | ofstream::app);
	sample << std::fixed;
	sample << std::setprecision(4);
	if (sample.is_open())
	{
		sample	<< left
			<< setw(10) << "ncycle" << "\t"
		      	<< setw(10) << "toten/N"<< "\n";
			//<< setw(10) << "avg_rho"<< "\n";
	}
	else{cerr	<< "unable to open sample file. \n";}


	sprintf(filename, "run%d_%.2f/accept.dat", ind, epsilon_xb);	//output file for acceptance
	accept.open(filename, ofstream::out | ofstream::app);
	accept << std::fixed;
	accept << std::setprecision(4);
	if (accept.is_open())
	{
		accept	<< left
			<< setw(10) << "ncycle"    << "\t"
		      	<< setw(10) << "moveaccept"<< "\t"
		      	<< setw(10) << "drmax"     << "\t"
			<< setw(10) << "volaccept" << "\t"
			<< setw(10) << "dvmax"     << "\n";
	}
	else{cerr	<< "unable to open acceptance file. \n";}

	nind = 0;
	int k, ran;
	double insE[5], sumE, insD[5], sumD;
	double acceptance = 0., acceptancev = 0.;

	for(int icyc=f_cyc; icyc<eqcycle; icyc++)
	{
	    if (icyc%1000 == 0) {cout << "progress(eq) " << icyc << "\t" << double(icyc)/double(eqcycle) << endl; }
	    for(int imov=f_mov; imov<nmov; imov++)
	    {	
		i_run = (icyc*nmov)+(imov+1);   
		ran   = round(ran2(idum)*(nmov));
//		if (ran < npart)
//		{
//			//cout << i_run << "\t" << U1 << endl;	
			mcmove(0);
//		}
//		else
//		{
//			mcvol();
//		}

	    	if( i_run % ((eqcycle*nmov)/eqframes) == 0 ){writexyz();}	//writing a frame for visualization
	    	if( i_run % ((eqcycle*nmov)/eqconfig) == 0 ){writeconfig(i_run,icyc,imov);}
	
	    	if( i_run % nadj == 0 )						//updating drmax to adjust acceptance ratio	
	    	{		
			if (movetry!=0) {acceptance = double(moveaccept)/double(movetry);}
			if (acceptance < 0.4) {drmax /= 1.05;}
			if (acceptance > 0.6) {drmax *= 1.05;}
			if (drmax > 0.5)      {drmax  = 0.5;}
			if (drmax < 1.0e-03)  {drmax  = 1.0e-03;}

			if (voltry!=0) {acceptancev = double(volaccept)/double(voltry);}
			if (acceptancev < 0.4) {dvmax /= 1.05;}
			if (acceptancev > 0.6) {dvmax *= 1.05;}
			if (dvmax > 0.5)       {dvmax  = 0.5;}
			if (dvmax < 1.0e-03)   {dvmax  = 1.0e-03;}

			if (accept.is_open())
			{
				accept  << left
					<< setw(10) << i_run			<< "\t"
					<< setw(10) << acceptance  		<< "\t"
					<< setw(10) << drmax       		<< "\t"
					<< setw(10) << acceptancev 		<< "\t"
					<< setw(10) << dvmax       		<< "\n";
			}
			else{cerr	<< "unable to open acceptance file for output. \n";}
			movetry     = 0;
			voltry      = 0;
			moveaccept  = 0.0;
			volaccept   = 0.0;
	    	}

	    	if( i_run > (nind*nsample_e - 5) )
	    	{
			k       = i_run%5;
			insE[k] = U; 
			insD[k] = npart/V;
	    	}
	    	if( i_run % nsample_e == 0 )
	    	{
			//cout << "sampling" << endl;
			sumE  = 0.;
			for (int j=0; j<5; j++) {sumE += insE[j];}
			sumE /= 5.;

			sumD  = 0.;
			for (int j=0; j<5; j++) {sumD += insD[j];}
			sumD /= 5.;

			writeconfig_inst(i_run,icyc,imov);
			nind += 1;
	
			if (sample.is_open())
			{
				sample << left
					<< setw(10) << i_run				<< "\t"
			       	<< setw(10) << sumE/double(npart)		<< "\n";
					//<< setw(10) << sumD				<< "\n";
			}
			else{cerr 	<< "unable to open sample file for eqm. \n";}

	    	}
	    }
	}
	accept.close();
	sample.close();
	return;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void production(int f_cyc, int f_mov)
{

	ofstream ofmuex, sample;
	char filename[100];					//output file for sample 

	sprintf(filename, "run%d_%.2f/sample_prod.dat", ind, epsilon_xb);
	sample.open(filename, ofstream::out | ofstream::app);
	sample << std::fixed;
	sample << std::setprecision(4);

	if (sample.is_open())
	{
		sample	<< left
			<< setw(10) << "ncycle" << "\t"
		      	<< setw(10) << "toten/N"<< "\n";
			//<< setw(10) << "avg_rho"<< "\n";
	}
	else{cerr	<< "unable to open sample file. \n";}

//	sprintf( filename, "run%d_%.2f/muex.txt", ind, epsilon_xb);
//	ofmuex.open (filename, ofstream::out | ofstream::trunc);
//	ofmuex << std::fixed;
//	ofmuex << std::setprecision(4);
//	if (ofmuex.is_open())
//	{
//		ofmuex	<< left
//			<< setw(10) << "itestA" << "\t"
//		      	<< setw(10) << "wtestA" << "\t"
//			<< setw(10) << "muexA"  << "\t"
//			<< setw(10) << "itestB" << "\t"
//		      	<< setw(10) << "wtestB" << "\t"
//			<< setw(10) << "muexB"  << "\n";
//	}
//	else{cerr	<< "unable to open ofmuex file. \n";}

//	distribution(1,0,0,0,0.);
//	wtest[0] = 0.; wtest[1] = 0.;
//	itest[0] = 0;  itest[1] = 0;
//	muex[0]  = 0.; muex[1]  = 0.;

	nind = 0, count1 = 0, count2 = 0;
	int k, ran;
	double insE[5], sumE, insD[5], sumD;
	aven = 0., avdens = 0.;

	for(int icyc=f_cyc; icyc<procycle; icyc++)
	{
	    if (icyc%1000 == 0) {cout << "progress(pro) " << icyc << "\t" << double(icyc)/double(procycle) << endl; }
	    for(int imov=f_mov; imov<nmov; imov++)
	    {	
		//cout << (eqcycle+icyc)*nmov+(imov+1) << endl;
		i_run = (icyc+eqcycle)*nmov+(imov+1);
		p_run = (icyc*nmov)+(imov+1); 
		ran   = round(ran2(idum)*(nmov));
//		if (ran < npart)
//		{
//			//cout << p_run << "\t" << U1 << endl;
			mcmove(0);
			//mcmove(1);
//		}
//		else
//		{
//			mcvol();
//		}

	    	if( p_run % ((procycle*nmov)/proframes) == 0){writexyz();} 		//writing a frame for visualization
	    	if( p_run % ((procycle*nmov)/proconfig) == 0){writeconfig(i_run,icyc,imov);}

	    	if( p_run > (nind*nsample_p - 5) )
	    	{
			k       = p_run%5;
			insE[k] = U; 
			insD[k] = npart/V;
	    	}
	    	if( p_run % nsample_p == 0)
	    	{
			//widom(ighost,0);
			//widom(ighost,1);
			sumE  = 0.;
			for (int j=0; j<5; j++) {sumE += insE[j];}
			sumE /= 5.;

			sumD  = 0.;
			for (int j=0; j<5; j++) {sumD += insD[j];}
			sumD /= 5.;

			writeconfig_inst(i_run,icyc,imov);
			aven   += sumE;
			avdens += sumD;
			nind   += 1;
	
			if (sample.is_open())
			{
				sample << left
					<< setw(10) << i_run				<< "\t"
			       	<< setw(10) << sumE/double(npart)		<< "\n";
					//<< setw(10) << sumD				<< "\n";
			}
			else{cerr 	<< "unable to open sample file for prod. \n";}
//			if (ofmuex.is_open())
//			{
//			   ofmuex     << left
//					<< setw(10) << itest[0]					<< "\t"
//					<< setw(10) << wtest[0]					<< "\t"
//					<< setw(10) << (-log(wtest[0]/double(itest[0])))*temp 	<< "\t"
//					<< setw(10) << itest[1]					<< "\t"
//					<< setw(10) << wtest[1]					<< "\t"
//					<< setw(10) << (-log(wtest[1]/double(itest[1])))*temp	<< "\n";
//			}
//			else{cerr 	<< "unable to open muex file for muex_sample. \n";}
	    	}
	    }
	}
	return;
	sample.close();
	//ofmuex.close();
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
//void widom_run()
//{
//	ofstream ofmuex, sample;
//	char filename[100];	
//	sprintf( filename, "run%d_%.2f/muex.txt", ind, epsilon_xb);
//	ofmuex.open (filename, ofstream::out | ofstream::trunc);
//	ofmuex << std::fixed;
//	ofmuex << std::setprecision(4);
//	if (ofmuex.is_open())
//	{
//		ofmuex	<< left
//			<< setw(10) << "itestA" << "\t"
//		      	<< setw(10) << "wtestA" << "\t"
//			<< setw(10) << "muexA"  << "\t"
//			<< setw(10) << "itestB" << "\t"
//		      	<< setw(10) << "wtestB" << "\t"
//			<< setw(10) << "muexB"  << "\n";
//	}
//	else{cerr	<< "unable to open ofmuex file. \n";}

//	sprintf(filename, "run%d_%.2f/sample_w.dat", ind, epsilon_xb);
//	sample.open(filename, ofstream::out | ofstream::app);
//	sample << std::fixed;
//	sample << std::setprecision(4);

//	if (sample.is_open())
//	{
//		sample	<< left
//			<< setw(10) << "ncycle" << "\t"
//		      	<< setw(10) << "toten/N"<< "\n";
//	}
//	else{cerr	<< "unable to open sample file. \n";}

//       distribution(1,0,0,0,0.);
//	wtest[0] = 0.; wtest[1] = 0.;
//	itest[0] = 0;  itest[1] = 0;
//	muex[0]  = 0.; muex[1]  = 0.;

//	wdens  = double(npart)/V;
//	wdensA = double(nA)/V;
//	wdensB = double(nB)/V;

//	for(int icyc=0; icyc<widom_try; icyc++)
//	{
//	    if (icyc%1000 == 0) {cout << "progress(widom) " << icyc << "\t" << double(icyc)/double(widom_try) << endl; }
//	    for(int imov=0; imov<nmov; imov++)
//	    {		
//	       i_run = (icyc+eqcycle+procycle)*nmov+(imov+1);
//		mcmove(1);
//		if(imov%50 == 0)
//		{
//			widom(ighost,0);
//			widom(ighost,1);
//		}
//		if(icyc%500 == 0)
//		{
//			if (sample.is_open())
//			{
//				sample << left
//					<< setw(10) << i_run				<< "\t"
//					<< setw(10) << U/double(npart)		<< "\n";
//			}
//			else{cerr 	<< "unable to open sample file for widom_sample. \n";}
//			if (ofmuex.is_open())
//			{
//			   ofmuex     << left
//					<< setw(10) << itest[0]					<< "\t"
//					<< setw(10) << wtest[0]					<< "\t"
//					<< setw(10) << (-log(wtest[0]/double(itest[0])))*temp 	<< "\t"
//					<< setw(10) << itest[1]					<< "\t"
//					<< setw(10) << wtest[1]					<< "\t"
//					<< setw(10) << (-log(wtest[1]/double(itest[1])))*temp	<< "\n";
//			}
//			else{cerr 	<< "unable to open muex file for muex_sample. \n";}
//		}
//	    }
//	}
//	sample.close();
//	ofmuex.close();

//	return;
//}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double pair_energy(const cube& A, const cube& B)
{
	double en=0.;
	vector r;
	double angi, angj, angmini, angminj, ang, r1; 
	//pi-pi interaction
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			vectsub(A.patPi[i], B.patPi[j],r);
			r.pbc();
			r1 = r.norm();
			if(r1 <= Picut){
				if(i%2 == 0) angi = r.dot(A.ej);
				else         angi = r.dot(-A.ej);
				if(j%2 == 0) angj = -r.dot(B.ej);
				else         angj = -r.dot(-B.ej);
				angi = acos(angi/r1);
				angj = acos(angj/r1);
				ang = (angi*angi + angj*angj)/delta2;
				en -= epsilon_pi*exp(-ang);
				//cout << "PI int " << i << " " << j << " " << angi << " " << angj << endl;
				//cout << r1 << " "; r.printvec();
			}
		}
	}
	//XB interaction - A-B
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			vectsub(A.patA[i],B.patB[j],r);
			r.pbc();
			r1 = r.norm();
			if(r1 <= XBcut){
				//ei unit vector
				if(i%2 == 0) angi = r.dot(A.ei);
				else         angi = r.dot(-A.ei);
				if(j%2 == 0) angj = -r.dot(B.ei);
				else         angj = -r.dot(-B.ei);
				angi    = acos(angi/r1);
				angj    = acos(angj/r1);
				angmini = angi; angminj = angj;
				//ej unit vector
				if(i<=1) angi = r.dot(A.ej);
				else     angi = r.dot(-A.ej);
				if(j<=1) angj = -r.dot(B.ej);
				else     angj = -r.dot(-B.ej);
				angi    = acos(angi/r1);
				angj    = acos(angj/r1);
				if(angi < angmini) angmini = angi;
				if(angj < angminj) angminj = angj;

				ang = (angi*angi + angj*angj)/delta2;
				en -= epsilon_xb*exp(-ang);
				//cout << "AB int " << i << " " << j << " " << angi << " " << angj << endl;
				//cout << r1 << " "; r.printvec();
			}
		}
	}
	//XB interaction - B-A
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			vectsub(A.patB[i],B.patA[j],r);
			r.pbc();
			r1 = r.norm();
			if(r1 <= XBcut){
				//ei unit vector
				if(i%2 == 0) angi = r.dot(A.ei);
				else         angi = r.dot(-A.ei);
				if(j%2 == 0) angj = -r.dot(B.ei);
				else         angj = -r.dot(-B.ei);
				angi = acos(angi/r1);
				angj = acos(angj/r1);
				angmini = angi; angminj = angj;
				//ej unit vector
				if(i<=1) angi = r.dot(A.ej);
				else     angi = r.dot(-A.ej);
				if(j<=1) angj = -r.dot(B.ej);
				else     angj = -r.dot(-B.ej);
				angi    = acos(angi/r1);
				angj    = acos(angj/r1);
				if(angi < angmini) angmini = angi;
				if(angj < angminj) angminj = angj;

				ang = (angi*angi + angj*angj)/delta2;
				en -= epsilon_xb*exp(-ang);
				//cout << "BA int " << i << " " << j << " " << angi << " " << angj << endl;
				//cout << r1 << " "; r.printvec();
			}
		}
	}
	return en;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double totalenergy()
{
	double toten = 0., rr1;
	bool   k1, k2;
	vector rr;
	GJKalgorithm gjkvar;
	for(int i=0; i<npart-1;i++){
		for(int j=i+1; j<npart; j++){
			k1 = gjkvar.CollisionDetection(colloid[i],colloid[j]);
			k2 = gjkvar.CollisionDetection(colloid[j],colloid[i]);
			if ((k1 == 1) || (k2 == 1)) {
				//cout << i <<  " " << j << " overlap\n";
				toten = 1.0e+050;
				return toten;
			} else{
				vectsub(colloid[i].center, colloid[j].center,rr);
				rr.pbc();
				rr1 = rr.norm();
				if(rr1 <= rcut){
					toten += pair_energy(colloid[i], colloid[j]);
				}
			}
		}
	}
	return toten/temp;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double potenergy(int j, const cube& A)
{
	double en = 0., rr1;
	bool   k1, k2;
	vector rr;
	GJKalgorithm gjkvar;
	for(int i=0; i<npart;i++){
		if(i!=j){
			k1 = gjkvar.CollisionDetection(colloid[i],A);
			k2 = gjkvar.CollisionDetection(A,colloid[i]);
			if ((k1 == 1) || (k2 == 1)) {
				//cout << i <<  " " << j << " overlap\n";
				en = 1.0e+050;
				return en;
			} else{
				vectsub(colloid[i].center, A.center,rr);
				rr.pbc();
				rr1 = rr.norm();
				if(rr1 <= rcut){
					en += pair_energy(colloid[i], A);
				}
			}
		}
	}
	return en;
}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
void mcmove(int prod)
{
	int 	     ipart;
	double      Uo, Un, arg, rand;
	vector      dr;
	cube        Ao, An;

	movetry += 1;
	Uo	 = 0.;
	Un	 = 0.;

	ipart    = int(ran2(idum)*npart); 		//Choose a particle randomly
	Ao       = colloid[ipart];
	An       = Ao;
	Uo       = potenergy(ipart, Ao);

	dr.x     = (ran2(idum)-0.5)*drmax;
	dr.y     = (ran2(idum)-0.5)*drmax;
	dr.z     = (ran2(idum)-0.5)*drmax;

	vectadd(Ao.center, dr, An.center);
	An.center.pbc();
	An.rotatecube(abs(ran2(idum)*3.),(ran2(idum)-0.5)*0.1745329, Lx/2.,Ly/2.,Lz/2.);
	//0.1745329 rad = 10 deg (-5deg to +5 deg rotation)
	Un       = potenergy(ipart, An);

//	if (prod) {distribution(2,moltp[ipart],0,1,Uo);}
	arg      = exp(-(Un-Uo)/temp);
	rand     = ran2(idum);
	if(rand < arg)			
	{
		//accepted
		moveaccept    += 1;
		colloid[ipart] = An;
		colloid[ipart].setvertices(Lx/2.,Ly/2.,Lz/2.);
		colloid[ipart].setpatches(Lx/2.,Ly/2.,Lz/2.);
		U		+= (Un-Uo)/temp;
	}
	return;

}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
//void mcvol()
//{
//	double Uo, Un, Vo, Vn, lnvn, boxo, boxn;
//	double fact, arg, rand;

//	voltry += 1;
//	Uo      = U;
//	Vo      = V;
//	lnvn    = log(Vo) + (ran2(idum)-0.5)*dvmax;		//perform a random walk in ln V
//	Vn      = exp(lnvn);
//	boxn    = pow(Vn, 1./3.);
//	boxo    = box;
//	fact    = boxn/boxo;
//	box     = boxn;
//	for (int i=0; i<npart; i++)
//	{
//		colloid[i].x *= fact;
//		colloid[i].y *= fact;
//		colloid[i].z *= fact;
//	}

//	Un      = totalenergy();
//	arg     = -((Un-Uo) + (pressure*(Vn-Vo)) - (double(npart+1)*log(Vn/Vo)*temp)) / temp;
//	rand    = ran2(idum);
//	if ( rand < exp(arg) )
//	{
//		//Accepted - update totalenergy, volume and density
//		volaccept += 1;
//		U 	   = Un;
//		V	   = pow(box,3.);
//		rho	   = npart/V;
//	}
//	else
//	{
//		//Rejected - restore old positions and box size
//		fact	   = boxo/boxn;
//		for (int i=0; i<npart; i++)
//		{
//			colloid[i].x *= fact;
//			colloid[i].y *= fact;
//			colloid[i].z *= fact;
//		}
//		box 	   = boxo;
//	}
//	return;
//}
//--------------------------------------------------------------------------------------------------------------
//	Widom
//--------------------------------------------------------------------------------------------------------------
//void widom(int ighost, int type)
//{
////     int type;
//	vector test, testp;
//	double Utest, rhoins, ecor;
//	
//	//yukawa correction
//	if (type == 0)
//	{
//		rhoins = double(nA)/(box*box*box);
//		ecor   = 2.*3.1415*rhoins*A_yuk*Xi*Xi*(rcut+Xi)*exp(-rcut/Xi);
//	}
//	//LJ correction
//	if (type == 1)
//	{
//		rhoins = double(nB)/(box*box*box);
//		ecor	= 8.*3.1415*rhoins*epsilon_hb*sigB3*r3*(1./9*r3-1./3.);
//	}

//    	for(int i=0; i<ighost; i++)
//	{
////	    	if (ran2(idum) < 0.5) {type = 0;}
////	    	else		     {type = 1;}	    
//	    	itest[type] += 1;
//	    	test.x       = (ran2(idum) - 0.5)*box;
//	  	test.y       = (ran2(idum) - 0.5)*box;            
//		test.z       = (ran2(idum) - 0.5)*box;
//	    	testp.x      = 1.0;
//	    	testp.y      = 0.0;
//		testp.z      = 0.0;
//	    	testp        = newpatch(testp);
//	    	moltp[npart] = type;
//	    
//	    	Utest        = potenergy(npart, test, testp);
//		Utest       += 2*ecor;
//	    	wtest[type] += exp(-Utest/temp);
////	    	cout << Utest << " " <<  wtest[type] << " " << type << endl;
//	    	distribution(2,type,1,0,Utest);
//    }
//    return;
//	
//}
//--------------------------------------------------------------------------------------------------------------
//	Distributiom_Real and test particle energies
//--------------------------------------------------------------------------------------------------------------
//void distribution(int i, int k, int test, int tail, double en)
//{
//    int    iu;
//    double enhelp, rhoins, ecor;
//    if (i == 1)		//initialise
//    {
//	umax   =  20.;
//	umin   = -20.;
//	factu  = double(nhis)/(umax-umin);
//	for (int j=0; j<nhis; j++)
//	{
//		p0[j][0] = 0.;	//fuA
//		p0[j][1] = 0.;	//fuB
//		p1[j][0] = 0.;	//guA
//		p1[j][1] = 0.;	//guB
//	}
//	ntest[0]  = 0; ntest[1]  = 0;
//	nreal[0]  = 0; nreal[1]  = 0;

//    }
//    else if (i == 2)	//sample histogram only during production cycle
//    {
//	if (tail)
//	{
//		//yukawa correction
//		if (k == 0)
//		{
//			rhoins = double(nA)/(box*box*box);
//			ecor   = 2.*3.1415*rhoins*A_yuk*Xi*Xi*(rcut+Xi)*exp(-rcut/Xi);
//		}
//		//LJ correction
//		if (k == 1)
//		{
//			rhoins = double(nB)/(box*box*box);
//			ecor	= 8.*3.1415*rhoins*epsilon_hb*sigB3*r3*(1./9*r3-1./3.);
//		}
//		enhelp = en + (2.*ecor);
//	}
//	else
//	{
//		enhelp = en;
//	}
//	iu     = int((enhelp-umin)*factu);
//	if (test)
//	{
//		ntest[k] += 1;
//		if((iu>0) && (iu<=nhis)) p0[iu][k] += 1.;
//	}
//	else
//	{
//		nreal[k] += 1;
//		if((iu>0) && (iu<=nhis)) p1[iu][k] += 1.;
//	}

//    }
//    else
//    {
//	char filename[100];
//	double f0[2], f1[2];
//	ofstream dist;
//	sprintf(filename, "run%d_%.2f/dist.dat", ind, epsilon_xb);
//	dist << std::fixed;
//	dist << std::setprecision(4);
//	dist.open(filename, ofstream::out | ofstream::trunc);
//	if(dist.is_open())
//	{
//		dist 	<< "#Number of test patchy particle samples    = " << ntest[0] << "\n"
//			<< "#Number of real patchy particle samples    = " << nreal[0] << "\n"
//			<< "#Number of test isotropic particle samples = " << ntest[1] << "\n"
//			<< "#Number of real isotropic particle samples = " << nreal[1] << "\n\n"
//			<< left
//			<< setw(10) << "#enhelp" 		    		 			<< "\t"
//			<< setw(10) << "p0A"  		<< "\t" << setw(7)  << "f0A"	<< "\t"
//			<< setw(10) << "p1A" 		<< "\t" << setw(7)  << "f1A" 	<< "\t"
//			<< setw(10) << "f0A-f1A/beta"	<< "\t" << setw(10) << "muexA"	<< "\t"
//			<< setw(10) << "p0B"  		<< "\t" << setw(7)  << "f0B"	<< "\t"
//			<< setw(10) << "p1B"  		<< "\t" << setw(7)  << "f1B" 	<< "\t"
//			<< setw(10) << "f0B-f1B/beta"	<< "\t" << setw(10) << "muexB"	<< "\n";
//		for (int j=0; j<nhis; j++)
//		{
//			enhelp    = umin + ((double(j)+0.5)/factu);
//			f0[0]     = log(p0[j][0]/double(ntest[0])) - (enhelp/(2.*temp));
//			f1[0]     = log(p1[j][0]/double(nreal[0])) + (enhelp/(2.*temp));
//			f0[1]     = log(p0[j][1]/double(ntest[1])) - (enhelp/(2.*temp));
//			f1[1]     = log(p1[j][1]/double(nreal[1])) + (enhelp/(2.*temp));
//			dist << left
//			     << setw(10) << enhelp 		    		 				<< "\t"
//			     << setw(10) << p0[j][0]/double(ntest[0])  << "\t" << setw(7)  << f0[0]	<< "\t"
//			     << setw(10) << p1[j][0]/double(nreal[0])  << "\t" << setw(7)  << f1[0] 	<< "\t"
//			     << setw(10) << (f1[0]-f0[0])*temp	     << "\t" << setw(10) << muex[0]<< "\t"
//			     << setw(10) << p0[j][1]/double(ntest[1])  << "\t" << setw(7)  << f0[1]	<< "\t"
//			     << setw(10) << p1[j][1]/double(nreal[1])  << "\t" << setw(7)  << f1[1] 	<< "\t"
//			     << setw(10) << (f1[1]-f0[1])*temp	     << "\t" << setw(10) << muex[1]<< "\n";
//		}
//	}
//	else { cerr << "Unable to open dist file\n";}
//	dist.close();
//	//plot f0 vs enhelp, f1 vs enhelp and (f1-f0)*temp vs enhelp
//	//(f1-f0) = beta*muex
//    }

//}
//--------------------------------------------------------------------------------------------------------------------------
//					
//--------------------------------------------------------------------------------------------------------------------------
double ran2(long *idum) // ran2 algorithm form Numecial Recipes 2nd Edition
{
	int 	     j;
	long 	     k;
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
//					
//--------------------------------------------------------------------------------------------------------------------------
void writeconfig(long int i_run, int icyc, int imov)
{	
	char filename[100], run[5];
	if (i_run <= (nmov*eqcycle)) {sprintf(run, "eq");}
	else		               {sprintf(run, "pro");}
	ofstream config;

	sprintf(filename, "run%d_%.2f/%s/config_%ld.dat", ind, epsilon_xb, run, count2 );
	config.open(filename, ofstream::out | ofstream::app);
	config << std::fixed;
	config << std::setprecision(4);
	if (config.is_open())
	{
		config << i_run << "\t" << icyc << "\t" << imov << "\t" << box << endl;
		for (int j=0; j<npart; j++)
		{
			config << j << endl;
			colloid[j].printcube(config);
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
	char filename[100], run[5];
	if (i_run <= (nmov*eqcycle)) {eqm = true;}
	else		               {eqm = false;}
	ofstream instconfig;

	sprintf(filename, "run%d_%.2f/config_inst%d.dat", ind, epsilon_xb, ind );
	instconfig.open(filename, ofstream::out | ofstream::trunc);
  	instconfig << std::fixed;
  	instconfig << std::setprecision(10);
	if (instconfig.is_open())
	{
		instconfig << i_run << "\t" << icyc << "\t" << imov << "\t" << eqm  << "\t" << box << "\t" << drmax << "\t" << dvmax << endl;
		for (int j=0; j<npart; j++)
		{
			instconfig << j << endl;
			colloid[j].printcube(instconfig);
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
	sprintf(filename, "run%d_%.2f/config.xyz", ind, epsilon_xb );
	config.open(filename, ofstream::out | ofstream::app);
  	config << std::fixed;
  	config << std::setprecision(4);
	if (config.is_open())
	{
		config << (npart*19)+ 1  << endl;
		config 	       	    << endl;
		config <<  "C1"  << "  "  << box/2. << "  " << box/2. << "  " << box/2. << endl;
		for(int i=0; i<npart; i++){
			colloid[i].printcubexyz(config);
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
	sprintf(parafile, "run%d_%.2f/para.%d.out", ind, epsilon_xb, ind);
	outfile.open(parafile, ofstream::out | ofstream::trunc);
	if (outfile.is_open())
	{
	  outfile<< ind		<< "\t\t index"		<< "\n"
		  << npart		<< "\t\t npart" 		<< "\n"
		  << nA		<< "\t\t npartA"		<< "\n"
		  << nB		<< "\t\t npartB"		<< "\n"
		  << volfrac		<< "\t\t volfrac"		<< "\n"
		  << rho		<< "\t\t rho"			<< "\n"
		  << box		<< "\t\t box"			<< "\n"
		  << V			<< "\t\t volume"		<< "\n"
		  << rcut 		<< "\t\t rcut"		<< "\n"
		  << temp 		<< "\t\t temp"		<< "\n"
		  << pressure		<< "\t\t pressure"		<< "\n"
		  << epsilon_xb	<< "\t\t epsilon_xb"		<< "\n"
		  << epsilon_pi	<< "\t\t epsilon_pi"		<< "\n"
		  << epsilon_hb	<< "\t\t epsilon_hb"		<< "\n"
		  << drmax		<< "\t\t drmax" 		<< "\n"
		  << dvmax		<< "\t\t dvmax" 		<< "\n"
		  << eqcycle		<< "\t\t eqcycle"		<< "\n"
		  << procycle		<< "\t\t procycle"		<< "\n"
		  << runcontinue	<< "\t\t runcontinue"	<< "\n"
		  << initial		<< "\t\t initial"		<< "\n"
		  << runrestart	<< "\t\t runrestart"		<< "\n"
		  << eqframes 	<< "\t\t eqframes"   	<< "\n"
		  << proframes	<< "\t\t proframes"  	<< "\n"
		  << ighost		<< "\t\t ighost"		<< "\n"
		  << scalefact	<< "\t\t scalefact"		<< "\n"
		  << "\n" "Change runrestart to true for restart run. \n";
	}
	else{cerr << "Unable to open file for final parameter output. \n";}
	outfile.close();

	ofstream outfile1;
	sprintf(parafile, "run%d_%.2f/para.%d.inp", ind, epsilon_xb, ind+1);
	outfile1.open(parafile, ofstream::out | ofstream::trunc);
	if (outfile1.is_open())
	{
	 outfile1<< nA		<< "\t\t nA"			<< "\n"
		  << nB		<< "\t\t nB"			<< "\n"
		  << volfrac		<< "\t\t volfrac"		<< "\n"
		  << rcut 		<< "\t\t rcut"		<< "\n"
		  << temp 		<< "\t\t temp"		<< "\n"
		  << pressure		<< "\t\t pressure"		<< "\n"
		  << epsilon_xb	<< "\t\t epsilon_xb"		<< "\n"
		  << epsilon_pi	<< "\t\t epsilon_pi"		<< "\n"
		  << epsilon_hb	<< "\t\t epsilon_hb"		<< "\n"
		  << drmax		<< "\t\t drmax" 		<< "\n"
		  << dvmax		<< "\t\t dvmax" 		<< "\n"
		  << eqcycle		<< "\t\t eqcycle"		<< "\n"
		  << procycle		<< "\t\t procycle"		<< "\n"
		  << runcontinue	<< "\t\t runcontinue"	<< "\n"
		  << "0"		<< "\t\t initial"		<< "\n"
		  << "1"		<< "\t\t runrestart"		<< "\n"
		  << eqframes 	<< "\t\t eqframes"   	<< "\n"
		  << proframes	<< "\t\t proframes"  	<< "\n"
		  << "0"		<< "\t\t start"		<< "\n"
		  << ighost		<< "\t\t ighost"		<< "\n"
		  << scalefact	<< "\t\t scalefact"		<< "\n";

	}
	else{cerr << "Unable to open file for final parameter output. \n";}
	outfile1.close();
	return;
}
