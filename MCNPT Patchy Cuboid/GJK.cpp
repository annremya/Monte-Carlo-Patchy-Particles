#include"GJK.h"

//Class Cube function definitions
void cube::setrandomcenter(const double& box)
{
	this->center = vector((ran2(idum)-0.5)*box,(ran2(idum)-0.5)*box,(ran2(idum)-0.5)*box);
}
//void cube::setrandomcenter(const double& box)
//{
//	this->center = vector(ran2(idum)*box,ran2(idum)*box,ran2(idum)*box);
//}
void cube::setaxes()
{
	this->ei = vector(1,0,0); this->ej = vector(0,1,0); this->ek = vector(0,0,1);
}

void cube::setvertices(const double& Lx, const double& Ly, const double& Lz)
{
	vector cent;
	cent = this->center+this->ek.vectimes(Lz);
	this->vert[0] = (cent+ this->ei.vectimes(Lx)) + (cent+ this->ej.vectimes(Ly));
	this->vert[1] = (cent+ this->ei.vectimes(Lx)) + (cent+-this->ej.vectimes(Ly));
	this->vert[2] = (cent+-this->ei.vectimes(Lx)) + (cent+ this->ej.vectimes(Ly));
	this->vert[3] = (cent+-this->ei.vectimes(Lx)) + (cent+-this->ej.vectimes(Ly));
	for(int j=0; j<4; j++){
		this->vert[j] = this->vert[j] - cent;
	}
	cent = this->center+-this->ek.vectimes(Lz);
	this->vert[4] = (cent+ this->ei.vectimes(Lx)) + (cent+ this->ej.vectimes(Ly));
	this->vert[5] = (cent+ this->ei.vectimes(Lx)) + (cent+-this->ej.vectimes(Ly));
	this->vert[6] = (cent+-this->ei.vectimes(Lx)) + (cent+ this->ej.vectimes(Ly));
	this->vert[7] = (cent+-this->ei.vectimes(Lx)) + (cent+-this->ej.vectimes(Ly));
	for(int j=4; j<8; j++){
		this->vert[j] = this->vert[j] - cent;
	}
}

void cube::setpatches(const double& Lx, const double& Ly, const double& Lz)
{
	this->patPi[0] = this->center +  this->ej.vectimes(Ly);	//+y
	this->patPi[1] = this->center + -this->ej.vectimes(Ly);	//-y
	//+x
	this->patA[0]  = this->patPi[0] + this->ei.vectimes(Lx) + -this->ek.vectimes(Lz/2.);
	this->patB[0]  = this->patPi[0] + this->ei.vectimes(Lx) +  this->ek.vectimes(Lz/2.);
	//-x
	this->patA[1]  = this->patPi[0] + -this->ei.vectimes(Lx) + this->ek.vectimes(Lz/2.);
	this->patB[1]  = this->patPi[0] + -this->ei.vectimes(Lx) +  -this->ek.vectimes(Lz/2.);
	//+x
	this->patA[2]  = this->patPi[1] + this->ei.vectimes(Lx) + this->ek.vectimes(Lz/2.);
	this->patB[2]  = this->patPi[1] + this->ei.vectimes(Lx) + -this->ek.vectimes(Lz/2.);
	//-x
	this->patA[3]  = this->patPi[1] + -this->ei.vectimes(Lx) + -this->ek.vectimes(Lz/2.);
	this->patB[3]  = this->patPi[1] + -this->ei.vectimes(Lx) + this->ek.vectimes(Lz/2.);
}

vector cube::fartherstPointInDirection(const vector& dir)
{
	vector r = this->vert[0];
	double maxdot = r.dot(dir);
	for(int i=1; i<8; i++){
		double dot = this->vert[i].dot(dir);
		//find biggest dot product
		if(dot > maxdot){
			maxdot = dot; r = this->vert[i];
		}
	}
	return r;
}

void cube::printcubexyz(ofstream& file)
{
	file << "C  " << this->center.x  << "  " << this->center.y  << "  " << this->center.z  << endl;
	for(int i=0; i<8; i++){
		file << "He " << this->vert[i].x  << "  " << this->vert[i].y  << "  " << this->vert[i].z  << endl;
	}
	for(int i=0; i<4; i++){
		file << "A  " << this->patA[i].x  << "  " << this->patA[i].y  << "  " << this->patA[i].z  << endl;
	}
	for(int i=0; i<4; i++){
		file << "B  " << this->patB[i].x  << "  " << this->patB[i].y  << "  " << this->patB[i].z  << endl;
	}
	for(int i=0; i<2; i++){
		file << "P  " << this->patPi[i].x << "  " << this->patPi[i].y << "  " << this->patPi[i].z << endl;
	}
}

void cube::printcubeverticesxyz(ofstream& file)
{
	file << "C  " << this->center.x  << "  " << this->center.y  << "  " << this->center.z  << endl;
	for(int i=0; i<8; i++){
		file << "He " << this->vert[i].x  << "  " << this->vert[i].y  << "  " << this->vert[i].z  << endl;
	}
}

void cube::printcube(ofstream& file){
	file << this->center.x << " " << this->center.y << " " << this->center.z << endl;
	file << this->ei.x     << " " << this->ei.y     << " " << this->ei.z     << endl;
	file << this->ej.x     << " " << this->ej.y     << " " << this->ej.z     << endl;
	file << this->ek.x     << " " << this->ek.y     << " " << this->ek.z     << endl;
}

void cube::rotatecube(int j, double theta, const double& Lx, const double& Ly, const double& Lz)
{
	this->ei = this->ei.rotation(j, theta);
	this->ej = this->ej.rotation(j, theta);
	this->ek = this->ek.rotation(j, theta);
	//cout << this->ei.dot() << " " << this->ej.dot() << " " << this->ek.dot() << endl;
	setvertices(Lx,Ly,Lz);
	setpatches(Lx,Ly,Lz);
}

//Class GJK function definitions
GJKalgorithm::GJKalgorithm()
{
	a = b = c = d = vector::zero();
}

GJKalgorithm::~GJKalgorithm()
{

}

vector GJKalgorithm::support(cube A, cube B, const vector& dir) const
{
	vector p1 = A.fartherstPointInDirection(dir);
	vector p2 = B.fartherstPointInDirection(-dir);
	vector p3 = p1-p2;
	return p3;
}

bool GJKalgorithm::line(vector& dir)
{
	vector ab = b - a;
	vector ao = -a;
	//can't be behind B
	//new direction
	dir = vector::doubleCross(ab, ao);
	c = b;
	b = a;
	nPtSimplex = 2;
	return false;
}

bool GJKalgorithm::triangle(vector& dir)
{
	vector ao  = -a;	//(origin-a)
	vector ab  = b - a;
	vector ac  = c - a;
	vector abc = vector::cross(ab,ac);
	//point can't be behind/in the direction of B,C or BC
	
	vector ab_abc = vector::cross(ab,abc);
	//is the origin away from ab edge? in the same plane
	//if ao is in that direction then
	if(ab_abc.dot(ao)>0){
		//change points
		c = b;
		b = a;
		//direction is not ab_abc bcz it dsnt point towards origin
		dir = vector::doubleCross(ab,ao);
		//direction changed, cant build tetrahedron
		return false;
	}
	
	vector abc_ac = vector::cross(abc,ac);
	//is origin away from ac edge? or is it in abc?
	//if ao is in that direction then
	if(abc_ac.dot(ao)>0){
		//keep c the same
		b = a;
		//direction is not abc_ac as it dsnt point towards origin
		dir = vector::doubleCross(ac,ao);
		//direction changed, cant build tetrahedron
		return false;
	}
	
	//now can build tetrahedron; check above or below
	if(abc.dot(ao)>0){
		//base of tetrahedron
		d   = c;
		c   = b;
		b   = a;
		dir = abc; // new direction
	} 
	else{
		//upside down tetrahedron
		d   = b;
		b   = a;
		dir = -abc;
	}
	nPtSimplex = 3;
	return false;
}

bool GJKalgorithm::tetrahedron(vector& dir)
{
	vector ao  = -a; //0-a
	vector ab  = b - a;
	vector ac  = c - a;

	//Case I - infront of triangle abc
	//dont have to change ao, ab, ac, abc
	vector abc =  vector::cross(ab,ac); // build triangle abc
	if(abc.dot(ao)>0){
		checkTetrahedron(ao,ab,ac,abc,dir);
	}

	//Case II - in front of triangle acd
	vector ad  = d - a;
	vector acd = vector::cross(ac,ad); // build triangle acd
	if(acd.dot(ao)>0){
		b   = c;
		c   = d;
		ab  = ac;
		ac  = ad;
		abc = acd;
		checkTetrahedron(ao,ab,ac,abc,dir);
	}

	//Case III - Infront of triangle adb
	vector adb = vector::cross(ad,ab); //build triangle adb
	if(adb.dot(ao)>0){
		c   = b;
		b   = d;
		ac  = ab;
		ab  = ad;
		abc = adb;
		checkTetrahedron(ao,ab,ac,abc,dir);
	}
	//origin in Tetrahedron
	return true;
}

bool GJKalgorithm::checkTetrahedron(const vector& ao, const vector& ab, const vector& ac, const vector& abc, vector& dir)
{
	vector ab_abc = vector::cross(ab,abc);
	if(ab_abc.dot(ao)>0){
		c   = b;
		b   = a;
		dir = vector::doubleCross(ab,ao);
		//dir is not ab_abc bcz it dsnt point towards origin
		//abxaoxab is the direction we r looking for
		//build new triangle, d will be lost
		nPtSimplex = 2;
		return false;
	}

	vector acp = vector::cross(abc,ac);
	if(acp.dot(ao)>0){
		b   = a;
		dir = vector::doubleCross(ac,ao);
		//acxaoxac is the direction we r looking for
		nPtSimplex = 2;
		return false;
	}
	//build new tetrahedron with new base
	d = c;
	c = b;
	b = a;
	dir =  abc;
	nPtSimplex = 3;
	return false;
}

bool GJKalgorithm::ContainsOrigin(vector& dir)
{
	if      (nPtSimplex == 2) {return triangle(dir);}
	else if (nPtSimplex == 3) {return tetrahedron(dir);}
	return false;
}

bool GJKalgorithm::CollisionDetection(const cube& A, const cube& B)
{
	//cout << "new------\n";
	vector dir = vector(1.,1.,1.);
	c   = support(A,B,dir);
	dir = -c; //negative direction
	
	b   = support(A,B, dir);
	if(b.dot(dir)<0) {return false;}

	dir = vector::doubleCross(c-b,-b);
	nPtSimplex = 2;	//begin with 2 points in simplex

	int steps = 0;	//avoid infinite loop
	while (steps<50){
		a = support(A,B,dir);
		//cout << "a.dot(dir) " << a.dot(dir) << " " << nPtSimplex << endl;
		//dir.printvec();
		if(a.dot(dir)<0) {
			return false;
		}
		else {
			if(ContainsOrigin(dir)) {
				//cout << a.dot(dir) << " " << steps << " true" << endl; 
				//cout << "true\n";
				return true;
			}
		}
		steps++;
	}
	return false;
}


