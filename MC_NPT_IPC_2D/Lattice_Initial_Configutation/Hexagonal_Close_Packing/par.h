!     Control parameters      
      Integer  nsample, nadj, nxyze, nxyzp, nsamplep
      Integer*8  eqcycle, procycle
      Integer nann
      Parameter (eqcycle=1d+9, procycle=6d+6)
      Parameter (nsample = 250000, nadj = 25000, nsamplep = 100 )
      Parameter (nann = 1, nxyze = 2000000, nxyzp = 2000)
!     Simulation parameters	
      Double precision  pres, rcut, T0, pi, alpha
      Parameter(pres = 1.0d0, rcut = 2.0d0, alpha = 1.0)
      Parameter (T0 = 1.0d0, pi = 3.14159d0)

!     Potential parameters	
      Double precision sig, zeta, kappa, bjerrum
      Double precision epsiln1, deg, epsiln2
      Double precision elp, elpnp, elnp
      Parameter (sig = 1.0d0)
      Parameter (zeta = 0.0d0, kappa = 0.0d0, bjerrum = 0.28d0)
      Parameter (elp = 1.5d0, elpnp = 1.28d0,elnp = 1.46d0)
	Parameter (epsiln1 = -6.0d0,epsiln2 = 6.0d0)
	Parameter (deg = 90.0d0)  !patch coverage = 50
      Logical restart
      Parameter (restart = .true.)     

! nann             : no. of annealing (Temperature) steps
! bjerrum          : Bjerrum length
! alpha            : cooling-schedule factor
!patch coverage % = 1 - cos(deg)/2
!deg = patch angle in deg

