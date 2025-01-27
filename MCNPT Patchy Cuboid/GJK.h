#ifndef GJK_algorithm
#define GJK_algorithm
//the above statement is to avoid redeclaration of variables

//Random number generator
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
double ran2(long*);
long *idum;

//CLass Vector
class vector{
	public: double x,y,z;
		 vector();
		 vector(double, double, double);
		 void   printvec();
		 void   pbc();
		 double norm();
		 double dot();
		 double dot(const vector&);
		 vector operator-() const;
		 vector operator-(const vector&);
		 vector operator+(const vector&);
		 vector vectimes(const double&, const double&, const double&);
		 vector vectimes(const double&);
		 vector rotation(int, double);
		 
		 static vector zero();
		 static vector cross(const vector&, const vector&);
		 static vector doubleCross(const vector&, const vector&);
};

//vector class function definitions
vector::vector() : x(0.0), y(0.0), z(0.0) { }

vector::vector(double x, double y, double z) : x(x), y(y), z(z) { }

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
vector vector::rotation(int i, double theta){
	vector a;
	//along x
	if(i == 0){
		return vector(this->x, cos(theta)*this->y + sin(theta)*this->z, -sin(theta)*this->y + cos(theta)*this->z);
	}
	//along y
	if(i == 1){
		return vector(cos(theta)*this->x + sin(theta)*this->z, this->y, -sin(theta)*this->x + cos(theta)*this->z);
	}
	//along z
	if(i == 2){
		return vector(cos(theta)*this->x + sin(theta)*this->y, -sin(theta)*this->x + cos(theta)*this->y, this->z);

	}

}


//Class vector static definitions
vector vector::zero(){
	return vector(0.,0.,0.);
}
vector vector::cross(const vector& v1, const vector& v2){
	return vector(v1.y*v2.z - v2.y*v1.z, v1.z*v2.x-v2.z*v1.x, v1.x*v2.y-v2.x*v1.y);
}
vector vector::doubleCross(const vector& v1, const vector& v2){
	return cross(cross(v1,v2),v1);
}


//Class Cube
class cube{
	public:
	vector center;
	vector ei, ej, ek;
	vector vert[8], patPi[2], patA[4], patB[4];
	vector fartherstPointInDirection(const vector&);
	void   setaxes();
	void   setrandomcenter(const double& box);
	void   setvertices(const double& Lx, const double& Ly, const double& Lz);
	void   setpatches(const double& Lx, const double& Ly, const double& Lz);
	void   printcube(ofstream&);
	void   printcubexyz(ofstream&);
	void   printcubeverticesxyz(ofstream&);
	void   rotatecube(int, double, const double& Lx, const double& Ly, const double& Lz);
};

//Class GJK
class GJKalgorithm{
	public:
		 GJKalgorithm();
		~GJKalgorithm();
		bool CollisionDetection(const cube&, const cube&);
	private:
		int    nPtSimplex = 0;
		vector a, b, c, d;
		vector support(cube, cube, const vector&) const;
		bool   ContainsOrigin(vector&);
		bool   line(vector&);
		bool   triangle(vector&);
		bool   tetrahedron(vector&);
		bool   checkTetrahedron(const vector&, const vector&, const vector&, const vector&, vector&);
};

#endif


