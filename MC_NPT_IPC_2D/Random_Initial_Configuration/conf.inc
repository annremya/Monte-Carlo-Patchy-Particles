!	conf.inc
	
        Integer N, Np
	  Integer Ng, Nt
	  Parameter (N = 500, Np = 1)
	  Parameter (Ng = 0, Nt = N + Ng )
        Double precision rx(N), ry(N), rz(N)
        Double precision rgx(Ng), rgy(Ng), rgz(Ng)
        Double precision ex(N, Np), ey(N, Np), ez(N, Np)
       
	common /conf1/ rx, ry, rz, ex, ey, ez, rgx, rgy, rgz
!       rx(i), ry(i), rz(i)          : position of particle i (2.2 mu)
!       rgx(i), rgy(i), rgz(i)          : position of particle i (3 mu)
!       ex(i,j), ey(i,j), ez(i,j)   : orientation of patch j on particle i
!       N = NUmber of 2.2mu m particles, Np = no ofpatches
!       Ng = number of negatively charged 3 mu m particles
