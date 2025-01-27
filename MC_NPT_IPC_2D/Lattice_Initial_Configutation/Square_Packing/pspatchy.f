!     Single 50% patch 2-d model for oppositely charged patchy particles (Inverse Patchy Colloids)
!     Interactions: Kern-Frenkel Potential - H-boding (square-well with angular term)(attraction and repulsion) 
!     Simulation Technique: MC NPT
!     Last modified 25-02-2016
!     Code written by M. Ethayaraja modified to IPC by Remya Ann Mathews KalapurakalAnn   
!-------------------------------------------------------------------------
       Program patchypom
!	 use IFPORT
       Implicit None
       Include 'par.h'
       Include 'conf.inc'


       Integer*8 ii, j, a, jj, ncycle,icycle, imove
       Integer k, nind, attemptv,Naccv, attempt, Nacc,fileind
       Double Precision  En, Ecut, T, AA, prefact, sumE, sumrho
       Double Precision  ratiom, ratiov, dr,dv, insE(5),insrho(5)
       Double Precision rseed, randf, box, ran, rho, rhoavg
       Character atom*4, lable*4,dumbox*3

       Common /pot/ AA, Ecut
 !      box = 50.0d0
       rseed = -20021105.0d0  !(always initialize the random no. with
!                             negative signed seed)  
       
	open (12, file='dr.dat')
	read(12,*) dr, dv
!	read(12,*) box
	close(12)
	open (10, file='initial.txt')
	read(10,*) dumbox, box
	close(10)
	print*,box
!	dr = 0.1d0
!	dv = 0.1d0
       prefact = zeta**2 * exp(kappa) / (4.0d0 * bjerrum)

       Call Initial(rseed,box)
       !write the initial configuration file for visualizationon
       Call configwrite(0,box)

       fileind = 50
       open (13, file='acceptance.dat')

       Do jj = 1, nann
          T  = (alpha**jj)*T0 
          AA = prefact*T
          Ecut = AA * exp(-kappa*rcut) / rcut
          Call totalenergy(box,En)

          write(fileind, *) 'Temperature:', T, 'Equilibration'
          write(fileind+1, *) 'Temperature:', T, 'Production'
          Write(13, *) 'Temperature:', T
          Write(13,*)  'MCcyle:     transaccep ratio:      dr:
     &                     volaccep ratio:	dv:'

          Do ii=1, 2
             if (ii .eq. 1) then   ! equilibration steps
                ncycle = eqcycle
             else                  ! production steps
                ncycle = procycle
             endif
             attempt = 0
	      attemptv = 0
             Nacc    = 0
             Naccv    = 0
             nind    = 0
	      rhoavg = 0.0d0
             Do icycle = 1, ncycle

	         ran = randf(rseed)*(N+1)+1
	         If (ran .le. N) then
		     Call mcmove(box,T, En, attempt, Nacc, dr, rseed)
	         Else
		     Call mcvol(box, T, En, attemptv, Naccv, dv, rseed)
	         Endif

                ! write configuration file for visualization
                If (ii .eq. 1) then
                   If(mod (icycle, nxyze) .eq. 0) then   
                      Call configwrite(3,box)     
                   Endif
                ! write instantaneous energy and density (averaged over 5
                ! configuration)   
                   if ( icycle .gt. (nind*nsample - 5) ) then
                        k = mod(icycle, 5)
                        insE(k+1) = En
			   insrho(k+1) = N/(box*box)
                   Endif
                   If ( mod(icycle, nsample) .eq. 0 ) then
                      sumE = 0.0d0
			 sumrho = 0.0
                      Do k=1, 5
                         sumE = sumE + insE(k)
			    sumrho = insrho(k) + sumrho
                      Enddo
                         sumE = sumE / 5.0
			    sumrho = sumrho / 5.0
                      Write(fileind, 6666) dble( (jj-1)*(eqcycle+
     $                procycle)+icycle), sumE/dble(Nt), sumrho 
                      nind = nind + 1
		      open (14, file='intmconfig.txt')
       			Do j=1, Ng
	   		 Write(14,*) rgx(j), rgy(j)!, rgz(j)
      			Enddo 
       			Do j=1, N                                   
          		 Write(14,*) rx(j), ry(j)!, rz(j)
          		 Do a = 1, Np
            		    Write(14, *) ex(j,a),ey(j,a)!,ez(j,a)
          		 Enddo
	 		Enddo
		  close(14)
                   Endif

                ! reset the displacement and volume   
                   If ( mod(icycle, nadj) .eq. 0) then
                      If(attempt .ne. 0 .and. attemptv .ne. 0)then
                         ratiom = real(Nacc)/real(attempt)
                         ratiov = real(Naccv)/real(attemptv)
                         write(13, 1313)icycle, ratiom, dr, ratiov, dv
                         if (ratiom .lt. 0.4d0) then
                            dr = dr / 1.05d0
                         elseif(ratiom .gt. 0.6d0) then
                            dr = dr * 1.05d0
                         endif
                         if (ratiov .lt. 0.4d0) then
                            dv = dv * 0.95     !/ 1.05d0
                         elseif(ratiov .gt. 0.5d0) then
                           dv = dv * 1.05d0
                         endif
                         if (dr .gt. 0.5d0) dr = 0.5d0
                         if (dr .lt. 1.0d-5) dr = 1.0d-5
                         if (dv .gt. 0.5d0) dv = 0.5d0
                         if (dv .lt. 1.0d-4) dv = 1.0d-4
                         attempt = 0
                         Nacc = 0
                         attemptv = 0
                         Naccv = 0
                      Endif
                   Endif
                Endif   
                
                If (ii. eq. 2) then
                   if ( icycle .gt. (nind*nsamplep - 5) ) then
                        k = mod(icycle, 5)
                        insE(k+1) = En
			   insrho(k+1) = N/(box*box)
                   Endif

                   If ( mod (icycle, nsamplep ) .eq. 0) then
                      sumE = 0.0
			 sumrho = 0.0
                      Do k = 1, 5
                         sumE = insE(k) + sumE
			    sumrho = insrho(k) + sumrho
                      Enddo
                         sumE = sumE / 5.0
			    sumrho = sumrho / 5.0
                      Write( fileind+1, 6666) dble((jj-1)*(eqcycle+
     $                procycle)+eqcycle+icycle), sumE/dble(Nt), sumrho
			 rho = N/(box*box) 
			 rhoavg = rhoavg + sumrho
                      open(66, file='rawdata.dat', ACCESS='append')
                      WRITE (66, *) icycle, N, rho, sumrho
                      close(66)
                      nind = nind + 1
		      open (14, file='intmconfig.txt')
       			Do j=1, Ng
	   		 Write(14,*) rgx(j), rgy(j)!, rgz(j)
      			Enddo 
       			Do j=1, N                                   
          		 Write(14,*) rx(j), ry(j)!, rz(j)
          		 Do a = 1, Np
            		    Write(14, *) ex(j,a),ey(j,a)!,ez(j,a)
          		 Enddo
	 		Enddo
		  close(14)
                   Endif
                   If(mod (icycle, nxyzp) .eq. 0)then   
                      Call configwrite(1,box)     
                   Endif
                Endif
             Enddo  ! end of (equilibration or production)
          Enddo  ! end of ii-loop (equilibrium and production)

          Open(11,file='config_final.txt')
          write(11,*) 'Temperature = ', 1.0/epsiln2
       Do j=1, Ng
	    Write(11,*) rgx(j), rgy(j)!, rgz(j)
       Enddo 
       Do j=1, N                                   
          Write(11,*) rx(j), ry(j)!, rz(j)
          Do a = 1, Np
             Write(11, *) ex(j,a),ey(j,a)!,ez(j,a)
          Enddo
	 Enddo


      write(26,*)
      write(26,*) 'Temperature=', 1.0/epsiln2
      write(26,*) 'running energy=', En/dble(Nt)
      Call totalenergy(box,En)
      write(26,*) 'final total energy=', En/dble(Nt)
      Write(26,*) 'Average density = ', rhoavg/dble(nind)
      write(26,*)
      open (12, file='dr.dat')
      write(12,*) dr, dv
      write(12,*) box
      close(12)
      fileind = fileind+2
      Enddo !end of temperature scaling

    
      Call configwrite(2,box)
6666  Format(2X, E14.7, 2(2X, F20.7))
1313  Format(2X, i15, 4(F20.7, 2X))
      Stop 
      End
!-----------------subroutine to write configuration-------------
       Subroutine configwrite(x,box)
       Implicit None
       Include 'par.h'
       Include 'conf.inc'
       Integer j, x, a, fileind,y
       Double precision box2, box3,box

       box2 = 0.5d0*box
	 box3 = 0.0000000
       If (x .eq. 0) then
           fileind = 30
           open (fileind, file ='initial.xyz', ACCESS = 'append')
       elseif(x .eq. 1) then
           fileind = 21
           open (fileind, file = 'runningprod.xyz')!, ACCESS= 'append')
       elseif(x .eq. 2) then
           fileind = 22
           open (fileind, file = 'final.xyz', ACCESS = 'append')
       elseif(x .eq. 3) then
           fileind = 23
           open (fileind, file = 'runningeq.xyz')!, ACCESS = 'append')
       Endif    
       write(fileind,*) ((N*(Np+1))+1+Ng)
       Write(fileind,*) 
!       Write(fileind, 2020) 'C', -box2, -box2, -box2
!       Write(fileind, 2020) 'C',  box2, -box2, -box2
!       Write(fileind, 2020) 'C',  box2,  box2, -box2
!       Write(fileind, 2020) 'C',  box2, -box2,  box2
!       Write(fileind, 2020) 'C', -box2,  box2, -box2
!       Write(fileind, 2020) 'C', -box2,  box2,  box2
!       Write(fileind, 2020) 'C', -box2, -box2,  box2
!       Write(fileind, 2020) 'C',  box2,  box2,  box2

 
!       Write(fileind, 2020) 'C', -box2, -box2, -box3
       Write(fileind, 2020) 'C',  box2, -box2, -box3
!       Write(fileind, 2020) 'C',  box2,  box2, -box3
!       Write(fileind, 2020) 'C',  box2, -box2,  box3
!       Write(fileind, 2020) 'C', -box2,  box2, -box3
!       Write(fileind, 2020) 'C', -box2,  box2,  box3
!       Write(fileind, 2020) 'C', -box2, -box2,  box3
!       Write(fileind, 2020) 'C',  box2,  box2,  box3



	 Do j=1,Ng
          Write(fileind,2020) 'Hf', rgx(j), rgy(j), box3
	 Enddo
       Do j=1, N                                    
          Write(fileind,2020) 'He1', rx(j), ry(j), box3
          Do a = 1, Np
             Write(fileind,2020) 'He2', rx(j)+ex(j,a)/20.0,ry(j)+ 
     &                                ey(j,a)/20.0, box3
          Enddo
       Enddo 

	 
2020   Format(2X, 1A, 3(F14.7, 2X))  
       Return
       End
!-----------------------------------------------------------------------
       ! This subroutine generates the initial configuration
!-----------------------------------------------------------------------
       Subroutine Initial(rseed,box)
       Implicit None
       Include 'par.h'
       Include 'conf.inc'
       Integer i, j, k, x
       Double precision a1, a2, a3, a4, a5, xr(3)
       Double precision rxij, ryij, rzij, rij2, box2
       Double precision rgxij, rgyij, rgzij, rgij2
       Double Precision rseed, randf, box
       Double precision temp1,dum1
	Character eq, temp*12, rho*4,dum*3
       
       box2 = 0.5d0 * box  
       If (restart) Then
           Open(11, file='initial.txt') 
	     read(11,*) dum, dum1, rho, dum1
 	     Do j = 1, Ng
              Read(11,*) rgx(j), rgy(j)!, rgz(j)
           Enddo       
           Do j = 1, N
              Read(11,*) rx(j), ry(j)!, rz(j)
              Do i = 1, Np
                 Read(11,*) ex(j,i), ey(j,i)!, ez(j,i)		 
              Enddo       
           Enddo 
    
      
       Else
       Do j = 1, Ng
!		print *, j
98           xr(1) = (1.d0 - 2.d0*randf(rseed)) * box2
             xr(2) = (1.d0 - 2.d0*randf(rseed)) * box2
!             xr(3) = (1.d0 - 2.d0*randf(rseed)) * box2  
!		print *, xr(1), xr(2), xr(3)  
             Do k=1, j - 1
                rgxij   = rgx(k) - xr(1)
                rgyij   = rgy(k) - xr(2)
!                rgzij   = rgz(k) - xr(3)
              rgxij  = rgxij - box * Anint (rgxij/box)
	        rgyij  = rgyij - box * Anint (rgyij/box)
!	        rgzij  = rgzij - box * Anint (rgzij/box)
                rgij2   = rgxij*rgxij + rgyij*rgyij! + rgzij*rgzij                          
                If ( sqrt(rgij2) .lt. 1.41d0 ) Then     !1.36+0.05
                     Goto 98
                Endif
             Enddo ! loop over k
                   rgx(j) = xr(1)
                   rgy(j) = xr(2)
!                   rgz(j) = xr(3)
!			 print *, rx(j)
                   write(6,99002) j               
          Enddo   !1st loop over j

          Do j = 1, N
!		print *, j
99           xr(1) = (1.d0 - 2.d0*randf(rseed)) * box2
             xr(2) = (1.d0 - 2.d0*randf(rseed)) * box2
!            xr(3) = (1.d0 - 2.d0*randf(rseed)) * box2  
!		print *, xr(1), xr(2), xr(3)  
            Do k=1, Ng
                rgxij   = rgx(k) - xr(1)
                rgyij   = rgy(k) - xr(2)
 !               rgzij   = rgz(k) - xr(3)
              rgxij  = rgxij - box * Anint (rgxij/box)
	        rgyij  = rgyij - box * Anint (rgyij/box)
!	        rgzij  = rgzij - box * Anint (rgzij/box)
                rgij2   = rgxij*rgxij + rgyij*rgyij! + rgzij*rgzij                          
                If ( sqrt(rgij2) .lt. 1.232d0 ) Then         !1.18+0.05
                     Goto 99
                Endif
             Enddo ! loop over k
             Do k=1, j - 1
                rxij   = rx(k) - xr(1)
                ryij   = ry(k) - xr(2)
 !               rzij   = rz(k) - xr(3)
              rxij  = rxij - box * Anint (rxij/box)
	        ryij  = ryij - box * Anint (ryij/box)
!	        rzij  = rzij - box * Anint (rzij/box)
                rij2   = rxij*rxij + ryij*ryij !+ rzij*rzij                          
                If ( sqrt(rij2) .lt. 1.05d0 ) Then          !1.0+0.05
                     Goto 99
                Endif
             Enddo ! loop over k
                   rx(j) = xr(1)
                   ry(j) = xr(2)
!                   rz(j) = xr(3)
!			 print *, rx(j)
                   write(6,99003) j              
          Enddo   !2nd loop over j

!---------fixing patches on the particles (components of unit vectors)
          a1 = 0.0d0
          a2 = 1.0d0
          a3 = 0.809015d0
          a4 = 0.309018d0
          a5 = 0.5d0
          ex(1,1) =  a2
          ey(1,1) =  a1
 !         ez(1,1) =  a1
          

          Do j = 2, N
            Do i = 1, Np
	         ex(j,i) = ex(1,i)
               ey(j,i) = ey(1,i)
!	         ez(j,i) = ez(1,i) 
            Enddo
	    Enddo	 

	 Endif
       
99002  Format('Random initialization:', /, i5,
     &         ' non-patchy particle placed randomly')
99003  Format('Random initialization:', /, i5,
     &         ' patchy particle placed randomly')
       Return
       End
!-----------------------------------------------------------------------
!     Constrained Monte Carlo simulation
!     Headgroups are held fixed; only counterions are moved
!-----------------------------------------------------------------------
       Subroutine mcmove(box,T, En, attempt, Nacc,  dr, rseed)
       Implicit None
       Include 'par.h'
       Include 'conf.inc'      
       Integer i, newN, Nacc, Nrot, attempt,k
       Double precision exi(Np), eyi(Np), Q(4)!,ezi(Np)
       Double precision vxi(Np), vyi(Np), vzi(Np), rxij, ryij, rzij, rij
       Double precision AA, T, Enew, Eold, En, Ecut, delE, delEB, dr
       Double precision rxn, ryn, rzn, check, checkmin
       Double precision q0, q1, q2, q3, R11, R12, R13, R21, R22, R23
       Double precision R31, R32, R33, arge, arg, box2
       Double precision rxio, ryio, rzio
       Double Precision rseed, randf, dtheta, PI2
       Double precision txi(Np), tyi(Np) , ti(Np), vi(Np)
       Double precision fact, fact1, box
       
       Common /pot/ AA, Ecut

       box2 = 0.5d0*box
       attempt = attempt + 1
	 dtheta = 0.5d0
	 PI2=4.D0*DATAN(1.D0)    !DATAN = double precision arc tangent

!------select a particle at random       
       newN = int(randf(rseed)*Nt) + 1
       If(newN .gt. Nt) newN = Nt
!------calculate energy old configuration

!------start move of patchy particle-----------------------------------
       If(newN .gt. Ng) then
	    k = newN - Ng
          rxio = rx(k)
          ryio = ry(k)
!          rzio = rz(k)
	    Do i = 1, Np
	       exi(i) = ex(k, i)
	       eyi(i) = ey(k, i)
!	       ezi(i) = ez(k, i)
          Enddo

!          Call energy(newN, rxio, ryio, rzio, exi, eyi, ezi, Eold)
          Call energy(box,newN, rxio, ryio, exi, eyi, Eold)
          rxn = rxio + ( 1.0d0 - 2.0d0*randf(rseed) ) * dr
          ryn = ryio + ( 1.0d0 - 2.0d0*randf(rseed) ) * dr
 !         rzn = rzio + ( 1.0d0 - 2.0d0*randf(rseed) ) * dr
!-----------------------Rotation-----------------------------------    
!          check = 0.0d0
!          Do while (check .lt. 0.995d0)
!             Call sphere (Q, rseed)
!             check = acos( Q(1) )
!          Enddo   
!          q0 = Q(1)
!          q1 = Q(2)
!          q2 = Q(3)
!          q3 = Q(4)
!          R11 = q0*q0 + q1*q1 - q2*q2 - q3*q3
!          R12 = 2.0*(q1*q2 + q0*q3)
!          R13 = 2.0*(q1*q3 - q0*q2)
!          R21 = 2.0*(q1*q2 - q0*q3)
!          R22 = q0*q0 - q1*q1 + q2*q2 - q3*q3
!          R23 = 2.0*(q2*q3 + q0*q1)
!          R31 = 2.0*(q1*q3 + q0*q2)
!          R32 = 2.0*(q2*q3 - q0*q1)
!          R33 = q0*q0 - q1*q1 - q2*q2 + q3*q3
!          Do i=1, Np
!          vxi(i) = R11*exi(i) + R12*eyi(i)! + R13*ezi(i)
!          vyi(i) = R21*exi(i) + R22*eyi(i)! + R23*ezi(i)
!!          vzi(i) = R31*exi(i) + R32*eyi(i) + R33*ezi(i)
!          exi(i) = vxi(i)
!          eyi(i) = vyi(i)
!!          ezi(i) = vzi(i)
!          Enddo

	    check = 1.0d0
          Do while (check .gt. 0.1745329d0)		!10 deg tolerance
!          Do while (check .gt. 1.570796d0)		!90 deg tolerance
!          Do while (check .gt. PI)
            Call Rotation (Q, rseed)
            Do i=1, Np
	      vxi(i) = cos(Q(1)) 
	      vyi(i) = sin(Q(1)) 
		fact = (vxi(i)*exi(i))+(vyi(i)*eyi(i))
            check = acos(fact)
		fact1 = (PI/180.0d0)*fact
		Enddo
          Enddo 
          Do i=1, Np
	       vxi(i) = cos(Q(1))
	       vyi(i) = sin(Q(1))
	       vi(i) = vxi(i)*vxi(i)+vyi(i)*vyi(i)

	       txi(i) = (dtheta*vxi(i))+exi(i)
	       tyi(i) = (dtheta*vyi(i))+eyi(i)

	       ti(i) = txi(i)*txi(i)+tyi(i)*tyi(i)
	       ti(i) = sqrt(ti(i))

	       txi(i) = txi(i)/ti(i)
	       tyi(i) = tyi(i)/ti(i)
 
!	       ti(i) = txi(i)*txi(i)+tyi(i)*tyi(i)
!	       ti(i) = sqrt(ti(i))
	
		 exi(i) = txi(i)
		 eyi(i) = tyi(i)
	    Enddo
!--------------------------------------------------------------------
!          Call energy(newN, rxn, ryn, rzn, exi, eyi, ezi, Enew)
          Call energy(box,newN, rxn, ryn, exi, eyi, Enew)
          delE  = Enew - Eold
          delEB = delE / T
          If ( delEB .lt. 100.0d0) Then
               arg = exp(-delEB)
               If (randf(rseed) .lt. arg) then
!---------------accepted
                Nacc = Nacc + 1
                En   = En + delE
                if (rxn .lt. -box2) rxn = rxn + box
                if (rxn .gt.  box2) rxn = rxn - box
                if (ryn .lt. -box2) ryn = ryn + box
                if (ryn .gt.  box2) ryn = ryn - box
!                if (rzn .lt. -box2) rzn = rzn + box
!                if (rzn .gt.  box2) rzn = rzn - box
                rx(k) = rxn
                ry(k) = ryn
!                rz(k) = rzn   
                Do i=1, Np
                   ex(k, i) = exi(i)
                   ey(k, i) = eyi(i)
!                   ez(k, i) = ezi(i)
                Enddo     
               Endif
          Endif
!--------start move of non patchy particle---------------------
	 Else
	    k = newN
	    rxio = rgx(k)
          ryio = rgy(k)
!          rzio = rgz(k)

!          Call energy(newN, rxio, ryio, rzio, exi, eyi, ezi, Eold)
          Call energy(box,newN, rxio, ryio, exi, eyi, Eold)
          rxn = rxio + ( 1.0d0 - 2.0d0*randf(rseed) ) * dr
          ryn = ryio + ( 1.0d0 - 2.0d0*randf(rseed) ) * dr
!          rzn = rzio + ( 1.0d0 - 2.0d0*randf(rseed) ) * dr

!          Call energy(newN, rxn, ryn, rzn, exi, eyi, ezi, Enew)
          Call energy(box,newN, rxn, ryn, exi, eyi, Enew)
          delE  = Enew - Eold
          delEB = delE / T
          If ( delEB .lt. 100.0d0) Then
               arg = exp(-delEB)
               If (randf(rseed) .lt. arg) then
!---------------accepted
                Nacc = Nacc + 1
                En   = En + delE
                if (rxn .lt. -box2) rxn = rxn + box
                if (rxn .gt.  box2) rxn = rxn - box
                if (ryn .lt. -box2) ryn = ryn + box
                if (ryn .gt.  box2) ryn = ryn - box
 !               if (rzn .lt. -box2) rzn = rzn + box
!                if (rzn .gt.  box2) rzn = rzn - box
                rgx(k) = rxn
                rgy(k) = ryn
!                rgz(k) = rzn      
               Endif
          Endif
	 Endif
!-------end move of non-patchy particle---------------------------------
       Return
       End

!-----------------------------------------------------------------------
!     Calculation of Energy
!-----------------------------------------------------------------------
      Subroutine totalenergy(box,En)
      
      Implicit None
      Include 'par.h'
      Include 'conf.inc'
      integer i, j, a, k, nlist(N), list(N,N)
      Double precision rxi, ryi, rzi, rxij, ryij, rzij, rij2, rij 
      Double precision rijm, rij2m, Uij, En, Ecut, hpi, sig2, rcut2, AA
      Double precision angmin, angmini, angminj, angi, angj,del,fij
	Double precision fact1, fact2
!      Double precision exi, eyi, ezi
      Double precision exi(Np), eyi(Np), ezi(Np),box
      Double precision rgxij, rgyij, rgzij, rgij2, rgij

      Common /pot/AA, Ecut
!    ------ Outer loops begins ---------------
!	print*, box
      hpi = pi
      sig2 = 2.0d0*sig*sig
      rcut2 = rcut*rcut
      del = (pi/180)*deg
      En = 0.0d0

	!non-patch--non-patch
	Do i = 1, Ng-1
!	   print *, i, '--------'
	rxi = rgx(i)
	ryi = rgy(i)
!	rzi = rgz(i)

        Do j = i+1, Ng      
!	     print *, j, '--------'
           rgxij = rgx(j) - rxi             	!distance between 1 and 2 from 1 to 2  is (x2-x1)
           rgyij = rgy(j) - ryi 
 !          rgzij = rgz(j) - rzi
          rgxij = rgxij - box * anint(rgxij/box)
          rgyij = rgyij - box * anint(rgyij/box)
 !         rgzij = rgzij - box * anint(rgzij/box)

          rgij2 = rgxij*rgxij + rgyij*rgyij! + rgzij*rgzij
          rgij  = sqrt(rgij2)
!	    fact = exp(1.0d0 - rgij)            !exp(-k*(rij-1))
	  
	    If(rgij .lt. rcut) then
	      If (rgij .lt. 1.36d0) then !if2
              En = 1.0d+50
		  Return
	      Elseif(rgij .lt. elnp) then    !repulsive
!		  En = En + ((A*fact)/rgij)
		  En = En + epsiln2
	      Endif
	    Endif
	  Enddo
	Enddo

      !non-patch--patch
	Do i = 1, Ng
	rxi = rgx(i)
	ryi = rgy(i)
!	rzi = rgz(i)

        Do j = 1, N		
           rgxij = rx(j) - rxi             	
           rgyij = ry(j) - ryi 
!           rgzij = rz(j) - rzi
          rgxij = rgxij - box * anint(rgxij/box)
          rgyij = rgyij - box * anint(rgyij/box)
!          rgzij = rgzij - box * anint(rgzij/box)

          rgij2 = (rgxij*rgxij) + (rgyij*rgyij) !+ (rgzij*rgzij)
          rgij  = sqrt(rgij2)
!	    fact = exp(1.0d0 - rgij)            !exp(-k*(rij-1))
	  
	    If(rgij .lt. rcut) then
	      If (rgij .lt. 1.18d0) then !if2
              En = 1.0d+50
		  Return
	      Elseif(rgij .lt. elpnp) then
		  Do a = 1, Np
 	         angj = (-rgxij*ex(j,a))+(-rgyij*ey(j,a))!+
    ! &                 (-rgzij*ez(j,a))
               angj = acos(angj/rgij)
		  Enddo
		  If(angj .lt. del) then        !repulsive
!		     En = En + ((B*fact)/rgij)
		     En = En + epsiln2
		  Elseif (angj .gt. del) then   !attractive
		     En = En + epsiln1
		  Endif
	      Endif
	    Endif
	  Enddo
	Enddo
  
      !patch--patch
	Do i = 1, N-1
	rxi = rx(i)
	ryi = ry(i)
	rzi = rz(i)

 	Do a = 1, Np
	 exi(a) = ex(i,a)
	 eyi(a) = ey(i,a)
!	 ezi(a) = ez(i,a)
	Enddo

        Do j = i+1, N	
           rxij = rx(j) - rxi             	
           ryij = ry(j) - ryi 
!           rzij = rz(j) - rzi
            rxij = rxij - box * anint(rxij/box)
            ryij = ryij - box * anint(ryij/box)
!            rzij = rzij - box * anint(rzij/box)

           rij2 = (rxij*rxij) + (ryij*ryij)! + (rzij*rzij)
           rij  = sqrt(rij2)
!	     fact = exp(1.0d0 - rij)            !exp(-k*(rij-1))
         
	   If (rij .lt. rcut) then          
	      If (rij .lt. 1.0d0) then 
              En = 1.0d+50
		  Return
             Elseif (rij .lt. elp) then
		  Do a = 1, Np
              angi = (rxij*exi(a))+(ryij*eyi(a))!+(rzij*ezi(a))  !rxij.eai
              angi = acos(angi/rij)
		  Enddo
	 	  Do a =1,Np
               angj = (-rxij*ex(j,a)) + (-ryij*ey(j,a))! + 
    ! &                 (-rzij*ez(j,a))
               angj = acos(angj/rij)
		  Enddo
              If ((abs(angi) .lt. del) .and. 
     &            (abs(angj) .lt. del))then     !repulsive-case1
	               En = En + epsiln2
              Elseif ((abs(angi) .lt. del) .and. 
     &                (abs(angj) .gt. del))then !attractive-case3
	               En = En + epsiln1
              Elseif ((abs(angi) .gt. del) .and. 
     &                (abs(angj) .lt. del))then !attractive-case4
	               En = En + epsiln1
              Elseif((abs(angi) .gt. del) .and. 
     &               (abs(angj) .gt. del))then  !repulsive-case2
	               En = En + epsiln2
	   	Endif
            Endif 
          Endif
	  Enddo
	Enddo
!	print*, 'En', En
!----------------------------------------------------------------------------------------
	
      Return
      End

!----------------------------------------------------------------
!       Energy calculation 
!----------------------------------------------------------------
!      Subroutine energy(newN, rxi, ryi, rzi, exi, eyi, ezi, Eni)
      Subroutine energy(box,newN, rxi, ryi, exi, eyi, Eni)
      Implicit none
      Include 'par.h'
      Include 'conf.inc'
      Integer newN, i, j, a, k, neighbor, nlist(N), list(N,N)
      Double precision exi(Np), eyi(Np), ezi(Np)
!      Double precision exi(Np), eyi(Np), ezi(Np)
      Double precision rxi, ryi, rzi, fij
      Double precision rxij, ryij, rzij, rij2, rij, rijm, rij2m
      Double precision Eni, Ecut, Uij, hpi, sig2, rcut2, AA,del
      Double precision angmin, angmini, angminj, angi, angj
	Double precision fact1, fact2,box
!      Double precision exi, eyi, ezi
      Double precision rgxij, rgyij, rgzij, rgij2, rgij
    
      Common /pot/AA, Ecut

      Eni = 0.0d0
	del = (pi/180)*deg 
      hpi = pi
      sig2 = 2.0d0*sig*sig
      rcut2 = rcut*rcut

	If(newN .gt. Ng) then
!	print*, 'patchy paticle chosen'
	goto 101
	Else 
!	print*, 'non patchy paticle chosen'
	goto 100
	Endif
	

 !-------------------------non patchy particle--------------------------------------
	!non-patch--non-patch
100	k = newN

        Do j = 1, Ng      
	   If(k .ne. j) then
           rgxij = rgx(j) - rxi             	!distance between 1 and 2 from 1 to 2  is (x2-x1)
           rgyij = rgy(j) - ryi 
!           rgzij = rgz(j) - rzi
          rgxij = rgxij - box * anint(rgxij/box)
          rgyij = rgyij - box * anint(rgyij/box)
!          rgzij = rgzij - box * anint(rgzij/box)

          rgij2 = rgxij*rgxij + rgyij*rgyij! + rgzij*rgzij
          rgij  = sqrt(rgij2)
!	    fact = exp(1.0d0 - rgij)            !exp(-k*(rij-1))
	  
	    If(rgij .lt. rcut) then
	      If (rgij .lt. 1.36d0) then
              Eni = 1.0d+50
		  Return
	      Elseif(rgij .lt. elnp) then   !repulsive
!		  Eni = Eni + ((A*fact)/rgij)
		  Eni = Eni + epsiln2
	      Endif
	    Endif
	   Endif
	Enddo
 
      !non-patch--patch

        Do j = 1, N		
           rgxij = rx(j) - rxi             	
           rgyij = ry(j) - ryi 
!           rgzij = rz(j) - rzi
          rgxij = rgxij - box * anint(rgxij/box)
          rgyij = rgyij - box * anint(rgyij/box)
!          rgzij = rgzij - box * anint(rgzij/box)

          rgij2 = (rgxij*rgxij) + (rgyij*rgyij)! + (rgzij*rgzij)
          rgij  = sqrt(rgij2)
!	    fact = exp(1.0d0 - rgij)            !exp(-k*(rij-1))
	  
	    If(rgij .lt. rcut) then
	      If (rgij .lt. 1.18d0) then
              Eni = 1.0d+50
		  Return
	      Elseif(rgij .lt. elpnp) then
		  Do a = 1, Np
 	         angj = (-rgxij*ex(j,a))+(-rgyij*ey(j,a))!+
!     &                 (-rgzij*ez(j,a))
               angj = acos(angj/rgij)
		  Enddo
		  If(angj .lt. del) then       !repulsive-case1
!		     Eni = Eni + ((B*fact)/rgij)
		     Eni = Eni + epsiln2
		  Elseif (angj .gt. del) then  !attractive-case2
		     Eni = Eni + epsiln1
		  Endif
	      Endif
	    Endif
	  Enddo
	Return
!----------------------------------------------------------------------------------------

!--------------------------------patchy-particle-----------------------------------------

      !patch--non-patch
101	k =  newN -Ng

        Do j = 1, Ng		
!	     print *, j, '--------'
           rgxij = rgx(j) - rxi             	
           rgyij = rgy(j) - ryi 
!           rgzij = rgz(j) - rzi
          rgxij = rgxij - box * anint(rgxij/box)
          rgyij = rgyij - box * anint(rgyij/box)
!          rgzij = rgzij - box * anint(rgzij/box)

          rgij2 = (rgxij*rgxij) + (rgyij*rgyij) !+ (rgzij*rgzij)
          rgij  = sqrt(rgij2)
!	    fact = exp(1.0d0 - rgij)            !exp(-k*(rij-1))
	  
	    If(rgij .lt. rcut) then
	      If (rgij .lt. 1.18d0) then !if2
              Eni = 1.0d+50
		  Return
	      Elseif(rgij .lt. elpnp) then
		  Do a = 1, Np
 	         angi = (rgxij*exi(a))+(rgyij*eyi(a))!+(rgzij*ezi(a))
               angi = acos(angi/rgij)
		  Enddo
		  If(angi .lt. del) then       !repulsive-case1
!		     Eni = Eni + ((B*fact)/rgij)
		     Eni = Eni + epsiln2
		  Elseif (angi .gt. del) then  !attractive-case2
		     Eni = Eni + epsiln1
		  Endif
	      Endif
	    Endif
	  Enddo

      !patch--patch

        Do j = 1, N	
	   If (k .ne. j) then
           rxij = rx(j) - rxi             	
           ryij = ry(j) - ryi 
!           rzij = rz(j) - rzi
          rxij = rxij - box * anint(rxij/box)
          ryij = ryij - box * anint(ryij/box)
!          rzij = rzij - box * anint(rzij/box)

          rij2 = (rxij*rxij) + (ryij*ryij) !+ (rzij*rzij)
          rij  = sqrt(rij2)
!	    fact = exp(1.0d0 - rij)            !exp(-k*(rij-1))
       
	    If (rij .lt. rcut) then          
	      If (rij .lt. 1.0d0) then 
              Eni = 1.0d+50
		  Return
            Elseif (rij .lt. elp) then  
		  Do a = 1, Np
               angi = (rxij*exi(a))+(ryij*eyi(a))!+(rzij*ezi(a))  !rxij.eai
               angi = acos(angi/rij)
		  Enddo
		  Do a = 1, Np
               angj = (-rxij*ex(j,a)) + (-ryij*ey(j,a))! + 
!     &                 (-rzij*ez(j,a))
               angj = acos(angj/rij)
		  Enddo
              If ((abs(angi) .lt. del) .and. 
     &            (abs(angj) .lt. del))then     !repulsive-case1
	               Eni = Eni + epsiln2
              Elseif ((abs(angi) .lt. del) .and. 
     &                (abs(angj) .gt. del))then !attractive-case3
	               Eni = Eni + epsiln1
              Elseif ((abs(angi) .gt. del) .and. 
     &                (abs(angj) .lt. del))then !attractive-case4
	               Eni = Eni + epsiln1
              Elseif((abs(angi) .gt. del) .and. 
     &               (abs(angj) .gt. del))then  !repulsive-case2
	               Eni = Eni + epsiln2
	   	  Endif
            Endif
          Endif
	   Endif
	Enddo
	Return
      End
!-------------------------------------------------------------
!-----------random number generator---------------------------
!-------------------------------------------------------------
       Double Precision Function randf(rseed)
       Implicit None

       Integer IM1, IM2, IA1, IA2, IQ1, IQ2, IR1, IR2, NTAB, IV
       Integer IMM1, Ndiv, IY, IDUM, IDUM2, j, k
       Double Precision EPS, RNMX, AM, rseed
!       Common /randseed/ rseed

       Parameter (IM1=2147483563,IM2=2147483399,
     &            IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     &            IR2=3791,NTAB=32)
        Parameter (EPS=1.2D-14,RNMX=1.0D00-EPS)
!    C RAN2 OF NUMERICAL RECIPES 2ND ED.
        Dimension IV(NTAB)
        Save IV,IY,IDUM2
        Data IDUM2/123456789/,IV/NTAB*0/,IY/0/

        IDUM=int(rseed)
        AM=1.d0/im1
        IMM1=IM1-1
        Ndiv=1+IMM1/NTAB

        IF (IDUM.LE.0) THEN
            IDUM=MAX(-IDUM,1)
            IDUM2=IDUM
            DO j = NTAB+8,1,-1
               K=IDUM/IQ1
               IDUM=IA1*(IDUM-K*IQ1)-K*IR1
               IF(IDUM.LT.0)IDUM=IDUM+IM1
               IF(J.LE.NTAB)IV(J)=IDUM
            ENDDO
            IY=IV(1)
        ENDIF

        K=IDUM/IQ1
        IDUM=IA1*(IDUM-K*IQ1)-K*IR1
        IF(IDUM.LT.0)IDUM=IDUM+IM1
        K=IDUM2/IQ2
        IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
        IF(IDUM2.LT.0)IDUM2=IDUM2+IM2
        J=1+IY/NDIV
        IY=IV(J)-IDUM2
        IV(J)=IDUM
        IF(IY.LT.1)IY=IY+IMM1

        randf=MIN(AM*IY,RNMX)

        rseed=float(IDUM)

        RETURN
        END
!--------------------------------------------------

!-----------------------------------------------------------------------
!      This subroutine generates a point on a sphere of radius R-Re
!      using Marsaglia Algorithm (Ref: UMS, Frenkel & Smit, pp-578).
!-----------------------------------------------------------------------
       Subroutine Sphere (Q, rseed)     
       Implicit None
       
       Double precision  Q(4), S12, S34, ran1, ran2, ran3, ran4, ranh  
       Double Precision rseed, randf

       S12  = 2.0d0
       Do while (S12 .ge. 1.0d0) 
          ran1 = 1.d0 - 2.d0 * randf(rseed)
          ran2 = 1.d0 - 2.d0 * randf(rseed)
          S12  = ran1*ran1 + ran2*ran2
       Enddo
       S34  = 2.0d0
       Do while (S34 .ge. 1.0d0) 
          ran3 = 1.d0 - 2.d0 * randf(rseed)
          ran4 = 1.d0 - 2.d0 * randf(rseed)
          S34  = ran3*ran3 + ran4*ran4
       Enddo
          ranh  = sqrt((1.0d0 - S12)/S34) 
          Q(1) = ran1
          Q(2) = ran2
          Q(3) = ran3*ranh 
          Q(4) = ran4*ranh
       Return
       End
!---------------------------------------------------------------------
!-----------------------------------------------------------------------
!      This subroutine generates a random orientation from 0-2pi
!-----------------------------------------------------------------------
       Subroutine Rotation (Q,rseed)  
!	 use IFPORT   
       Implicit None

       
       Double precision  Q(2),PI	
       Double precision  randf, rseed	  !comment out use in gfortran and use rand(0), use random(0) or rand(0) with use ifport in intel fortran

	 PI=4.D0*DATAN(1.D0)    		  !DATAN = double precision arc tangent

!       Q(1) = (rand(0)*(2.0*PI))		  ![rand(0)*(y-x)]+x  --> TO GET A RANDOM NUMBER IN [x,y]
!	 Q(1) = (random(0)*(2.0*PI))
	 Q(1) = (randf(rseed)*(2.0*PI))
	 Q(2) = (180.d0/PI)*Q(1)
	
!	 print *, 'Q', Q(1),Q(2)

       Return
       End
!-----------------------------------------------------------------------
!     Attempt to change the volume
!-----------------------------------------------------------------------
       Subroutine mcvol(box, T,En, attemptv, Naccv,  dv, rseed)
       Implicit None
       Include 'par.h'
       Include 'conf.inc'      
       Integer i, Naccv, attemptv
       Double precision rxij, ryij, rzij, rij
       Double precision arg, T, AA, Ecut, En
       Double Precision rseed, randf, dv
       Double precision box, boxn, fact,ran
       Double precision Uo, Un, vo, vn, lnvn

       Common /pot/ AA, Ecut

       attemptv = attemptv + 1

!------calculate total old energy
       Call totalenergy(box,Uo)
	vo = box*box		!box**3 in case of 3d
	lnvn = log(vo)+(randf(rseed)-0.5)*dv
	vn = exp(lnvn)
	
	boxn = sqrt(vn)	!vn**(1/3) in case of 3d
	fact = boxn/box

	Do i = 1, N
	   rx(i) = rx(i)*fact
	   ry(i) = ry(i)*fact
!	   rz(i) = rz(i)*fact
	Enddo

       Call totalenergy(boxn,Un)

	arg = -((Un-Uo)+pres*(vn-vo)-(N+1)*log(vn/vo)*T)/T
	ran = randf(rseed)
	If (ran .lt. exp(arg)) Then
!----------accepted
           Naccv = Naccv + 1
           En   = Un
	    box = boxn
	Elseif(ran .gt. exp(arg)) Then
	    fact = box/boxn
	    Do i = 1, N
	       rx(i) = rx(i)*fact
	       ry(i) = ry(i)*fact
	    Enddo
	Endif
       Return
       End

!-----------------------------------------------------------------------
