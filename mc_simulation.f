*****************************************************************************
*          Molecular Dynamics code to simulate at NVE collectivity 
*          a system of N atoms of Ar inside a cubic box, interacting
*          through a truncated Lennard-Jones force field at 2.5 sigma. 
*          Leap-frog Verlet integration algorithm.
*
*****************************************************************************
      PROGRAM mc_simulation
      implicit double precision(a-h,o-z)

c     1. Defining dimensions
 
      dimension r(3,1000), U_vals(1000), P_vals(1000), g(1000)
      dimension Ur(1000), dUdr(1000), rhis(1000)

c     2. Reading data and computing related quantities

      open(1,file='leap-lj.data',status='old')
      read(1,*) ncycles
      read(1,*) npart
      read(1,*) tref
      read(1,*) pref 
      read(1,*) sigma,epsil
      read(1,*) deltax
      read(1,*) nhis 


      close(1)

      nf = 3*npart-3 !number of degrees of freedom 
      rc = 2.5d0 !sigma = range of the potential in reduced units
      beta = 1.d0 / tref

c     3. Reading initial configuration (positions, velocities) in A and A/ps

      open(2,file='leap-conf.data',status='old')
      do is = 1,npart
         read(2,*) (r(l,is),l=1,3)
         read(2,*) 
      end do         
      read(2,*) box
      close(2)
c     


      do is = 1, nhis
            g(is) = 0
      enddo
c     Opening files to write results
c
      open(3,file='energy-leap.dat',status='unknown')
      open(4,file='temp-leap.dat',status='unknown')
c
c     4. Change to reduced units

      call reduced(npart,r,box,sigma)

c     5. Start the loop to generate new configurations 


      isuccess = 0
      isamp = 0
      nsamp = int(ncycles/999) ! Interval entre samples
      pi = 4.d0*DATAN(1.d0)
      rho = npart / (box**3)



      do icycle = 1,ncycles
         call moveparticle(npart, r, rc, box, deltax, beta, isuccess)

         if (mod(icycle, nsamp) == 0) then
            call sampleVals(isample, npart, nhis, r, rc, box, beta,
     &                  U_vals, P_vals, g)
         endif
         if (isample>1000) then
            print *, "Too many samples!"
            exit
         endif
      enddo
      close(3)
      close(4)

      isample = isample - 1
      print *, 'Success rate: ', int(abs(100*isuccess/ncycles)), "%"
      print *, 'Number of samples:', isample
      print *, 'Number of cycles:', icycle
      print *, 'Reduced density:', rho

      uAvg = 0.d0
      pAvg = 0.d0
      uStd = 0.d0
      pStd = 0.d0

      print *, "====== Using samplig of MC simulation ======"

c     Pressure tail correction assuming constant g(r) = 1 for LJ potential
      dP = 16.d0*pi*(rho ** 2)*(1**6)*1/(3.d0 * (rc**3))
     & * ((1.d0/3.d0)*(1/rc)**6 - 1.d0)

      dU = 8.d0*pi/3.d0 * rho  * npart *(1.d0/3.d0 * (1.d0/rc)**6 - 1) 
     & * (1.d0/rc)**3

      print *, 'Energy correction: ', dU
      print *, 'Pressure correction: ', dP


      call stats(U_vals, isample, uAvg, uStd)
      call stats(P_vals, isample, pAvg, pStd)
      uAvg = uAvg + dU
      pAvg = pAvg + dP

      print *, 'U mean: ', uAvg, "std: ", uStd, " (", 
     &      int(100*uStd/uAvg), "% )"

      print *, 'P mean: ', pAvg, "std: ", pStd, " (", 
     &      int(100*pStd/pAvg), "% )"


c     6. Saving last configuration in A and A/ps 

      open(11,file='newconf.data',status='unknown')
      do is = 1,npart
         write(11,*) (r(l,is)*sigma,l=1,3)
         write(11,*) ! Write empty line
      end do         
      write(11,*) box*sigma
      close(11)

      delg = box/(2.d0*nhis) ! bin size
      do is = 1, nhis
            rhis(is) = delg*(dfloat(is-1)+0.5d0)
            vb = ((is+1)**3 - is**3)*delg**3
            rnid = (4.d0/3.d0)*pi*vb*rho
            g(is) = g(is) /(isample * npart*rnid)

            Ur(is) = g(is) * 4.d0/rhis(is)**4 * (1/rhis(is)**6 - 1) ! El integrando de la U, es decir, U(r)*r^2 * g
            dUdr(is) = g(is) * 24.d0/rhis(is)**4 * (-2/rhis(is)**6 + 1) ! El integrando de la P, es decir, dU/dr * r^3 * g
      end do

      ! Hacemos la integrales

      U_int = 2.d0*pi*rho*npart*tegrate_simpson(Ur, delg, nhis)
      P_int = rho/beta
     & - 2.d0*pi/3.d0 * rho**2 * tegrate_simpson(dUdr, delg, nhis)
      
      ! Tail corrections para la U, P integradas (muy parecidas a las de antes pero con un rc distinto)
      rc = rhis(nhis)
      dP_int = 16.d0*pi*(rho ** 2)*(1**6)*1/(3.d0 * (rc**3)) 
     & * ((1.d0/3.d0)*(1/rc)**6 - 1.d0)

      dU_int = 8.d0*pi/3.d0 * rho* npart *(1.d0/3.d0 * (1.d0/rc)**6 - 1) 
     & * (1.d0/rc)**3

      print *, "====== Using integration of g(r) ======"

      print *, 'Energy correction for integrated U: ', dU_int
      print *, 'Pressure correction for integrated P: ', dP_int


      U_int = U_int + dU_int
      P_int = P_int + dP_int
      print*, "Total energy using integration, U = ", U_int 
      print*, "Total pressure using integration, P = ", P_int 

      open(5,file='g.data',status='unknown')
      do is = 1,nhis
         write(5,*) delg*dfloat(is-1), g(is)
      end do         
      close(5)

      open(7,file='U.data',status='unknown')
      do is = 1,1000
         write(7,*) is, U_vals(is)
      end do
      close(7)

      open(8,file='P.data',status='unknown')
      do is = 1,1000
         write(8,*) is, P_vals(is)
      end do         
      close(8)

      stop
      end

*********************************************************
*********************************************************
c              subroutine reduced
*********************************************************
*********************************************************
      subroutine reduced(npart,r,box,sigma)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)

      rgas = 8.314472673d0 !J/(mol*K) 
      box = box/sigma
      do is = 1,npart
         do l = 1,3
            r(l,is) = r(l,is)/sigma
         end do
      end do

      return
      end

*********************************************************
*********************************************************
c              subroutine moveparticle
*********************************************************
*********************************************************
      subroutine moveparticle(npart, r0, rc, box, deltax,
     & beta, isuccess)
      implicit double precision(a-h,o-z)
      dimension r0(3,1000), rn(3,1000), g(1000)
c     This subroutine choses a particle at random and attempts a move

c     Notation: variable0 -> Old values
c               variablen -> new attempted values

        
      ! Chose particle at random
      iSel = 1 + int(rand()*npart)
      if (iSel > npart) iSel = npart

      do i = 1, npart
        do l = 1,3
            rn(l,i) = r0(l,i)
         enddo
      enddo

      ! Move a particle at random
      do l = 1,3
         rn(l,iSel) = r0(l,iSel) + deltax*(rand()-0.5d0)
      enddo

c      call boundaryConds(npart, rn, box)


      ! Energy difference
      call deltaEnergy(npart, r0, rn, iSel, box, rc, deltaU)
      acc = exp(-beta*deltaU) ! Acceptance probability.

      roll = rand()
      if(roll<= acc) then
         do i = 1, npart ! We accept
            do l = 1,3
               r0(l,i) = rn(l,i)
            enddo
         enddo
         call boundaryConds(npart, r0, box)
         isuccess = isuccess + 1
      endif ! We do nothing

c      write(6,*) isel,  U0, Un, deltaU, deltau2, acc, roll, isuccess

      end

*********************************************************
*********************************************************
c              subroutine deltaEnergy
*********************************************************
*********************************************************
      subroutine deltaEnergy(npart, r0, rn, iSel,box, rc, deltaU)
      implicit double precision(a-h,o-z)
      dimension rn(3,1000), r0(3,1000), g(1000)
c     This subroutine calculates the energy difference between two
c     configurations where the only difference is particle i
      deltaU = 0.d0 
      U0 = 0.d0
      Un = 0.d0

c      atom-atom interactions

      do js = 1,npart
         if (js == iSel)  then
            cycle
         else
            call lj(iSel, js, r0, rr, box, rc, U0, Pkin)
            call lj(iSel, js, rn, rr, box, rc, Un, Pkin)
            deltaU = deltaU + Un - U0

         end if
      end do
      return
      end

*********************************************************
*********************************************************
c              subroutine getValues
*********************************************************
*********************************************************
      subroutine getValues(npart, nhis, r, box, rc, Utot, Pkin, g)
      implicit double precision(a-h,o-z)
      dimension r(3,1000), g(1000)
c     This subroutine calculates the energy of a configuration


      Utot = 0.d0
      Pkin = 0.d0
      rr = 0.d0
      delg = box/(2.d0*nhis) ! bin size


c      atom-atom interactions
      do is = 1,npart-1
c         print *, 'Atom ', is, ' coords = ', r(1,is), r(2,is), r(3,is)
         do js = is+1,npart
            call lj(is,js,r, rr, box,rc,Uij,rFij)
            Utot = Utot + Uij
            Pkin = Pkin + rFij/(3*box**3)

            ig = int(rr/delg)

            if(rr < box/2.d0) then
               if (ig>0) then
                     g(ig) = g(ig) + 2
               endif
            endif
         end do
      end do
      end
*********************************************************
*********************************************************
      subroutine lj(is,js,r, rr, box,rc,Uij, rFij)
      implicit double precision(a-h,o-z)
      dimension r(3,1000), rij(3), g(1000)

      rr2 = 0.d0
      Uij = 0.d0
      rijl = 0.d0
      rFij = 0.d0
      do l = 1,3
         rijl = r(l,js) - r(l,is)
         rij(l) = rijl - box*dnint(rijl/box)
         rr2 = rr2 + rij(l)**2
      end do

      rr = dsqrt(rr2)

      if (rr < rc) then
         ynvrr2 = 1.d0/rr2
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
         Uij = 4.d0*(ynvrr12-ynvrr6)  
         rFij = forcedist*rr2
      end if



      return
      end

*********************************************************
*********************************************************
c              subroutine boundaryConds
*********************************************************
*********************************************************
      subroutine boundaryConds(npart, r, box)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
c     This subroutine applies periodic boundary conditions
c     (box with a corner at the origin)

      do i = 1, npart
         do l = 1,3
            if(r(l,i) < 0 ) r(l,i) = r(l,i) + box*ceiling(-r(l,i)/box)
            if(r(l,i) > 0 ) r(l,i) = r(l,i) - box*floor(r(l,i)/box) 
         enddo
      end do

      return
      end

*********************************************************
*********************************************************
c              subroutine sample
*********************************************************
*********************************************************
      subroutine sampleVals(isample, npart, nhis, r, rc, box, beta,
     &                  U_vals, P_vals, g)

      implicit double precision(a-h,o-z)
      dimension r(3,1000), U_vals(1000), P_vals(1000), g(1000)
      ! write(6, *) "Sample: ", isample
      Utot = 0.d0
      Pkin = 0.d0
      rho = npart/box**3

      call getValues(npart, nhis, r, box, rc, Utot, Pkin, g)
      


      U_vals(isample) = Utot
      P_vals(isample) = Pkin + rho/beta

      isample = isample + 1

      end



*********************************************************
*********************************************************
c              functions
*********************************************************
*********************************************************
      subroutine stats(A, n, avg, std)
      implicit double precision(a-h,o-z)
      dimension A(1000)
      avg = 0.d0
      std = 0.d0
      x2 = 0.d0

      do i = 1, n
            avg = avg + A(i) / n
            x2 = x2 + A(i)**2 / n
      enddo
      std = sqrt(x2 - avg**2)
      end 
*********************************************************
*********************************************************
c              functions
*********************************************************
*********************************************************
      double precision function tegrate_simpson(x, dt, N) ! Formula de Simpson
      implicit none

      integer, intent(in) :: N
      integer :: i
      double precision, dimension(1000), intent(in) :: x(1:1000)
      double precision, intent(in) :: dt

    

      tegrate_simpson = 0.0

      do i = 2, N-1, 2
            tegrate_simpson = tegrate_simpson + 
     &      dt * (1.0/3.0* x(i-1) + 4.0/3.0 * x(i)+ 1.0/3.0 * x(i+1))
      end do

      return
      end