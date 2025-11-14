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
      rho = npart / box**3

c     Pressure tail correction assuming constant g(r) = 1 for LJ potential
      dP = 16.d0*pi*(rho ** 2)*(1**6)*1/(3.d0 * (rc**3))
     & * ((1.d0/3.d0)*(1/rc)**6 - 1.d0)
      print *, 'Pressure correction: ', dP

      dU = 8.d0*pi/3.d0 * rho  * npart * (1.d0/3.d0 * (1.d0/rc)**6 - 1)
     & * (1.d0/rc)**3
      print *, 'Energy correction: ', dU


      do icycle = 1,ncycles
         call moveparticle(npart, r, rc, box, deltax, beta, isuccess)

         if (mod(icycle, nsamp) == 0) then
            call sampleVals(isample, npart, nhis, r, rc, box, beta,
     &                  U_vals, P_vals, g, dU, dP)
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

      uAvg = 0.d0
      pAvg = 0.d0
      uStd = 0.d0
      pStd = 0.d0
      call stats(U_vals, isample, uAvg, uStd)
      call stats(P_vals, isample, pAvg, pStd)

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
            vb = ((is+1)**3 - is**3)*delg**3
            rnid = (4.d0/3.d0)*pi*vb*rho
            g(is) = g(is) /(isample * npart*rnid)
      end do

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
c      call getValues(npart, nhis, r0, box, rc, U0, Pkin, g)
c      call getValues(npart, nhis, rn, box, rc, Un, Pkin, g)
      call Energy(npart, r0, iSel, box, rc, U0)
      call Energy(npart, rn, iSel, box, rc, Un)
      call deltaEnergy(npart, r0, rn, iSel, box, rc, deltaU)
      deltaU2=Un-U0


c      deltaU = U0 - Un
      acc = exp(-beta*deltaU2) ! Acceptance probability.
c     If deltaU<0 -> -beta*deltaU >0 -> acc > 0

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

      write(6,*) isel,  U0, Un, deltaU, deltau2, acc, roll, isuccess

      end

*********************************************************
*********************************************************
c              subroutine deltaEnergy
*********************************************************
*********************************************************
      subroutine deltaEnergy(npart, r0, rn, iSel,box, rc, deltaU)
      implicit double precision(a-h,o-z)
      dimension rn(3,1000), r0(3,1000)
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
c              subroutine deltaEnergy
*********************************************************
*********************************************************
      subroutine Energy(npart, r0, iSel,box, rc, u)
      implicit double precision(a-h,o-z)
      dimension rn(3,1000), r0(3,1000)
c     This subroutine calculates the energy difference between two
c     configurations where the only difference is particle i


      U = 0.d0

c      atom-atom interactions

      do js = 1,npart
         if (js == iSel)  then
            cycle
         else
            call lj(iSel, js, r0, rr, box, rc, dU, Pkin)
            U = U + dU

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
      delg = box/(2.d0*nhis) ! bin size


c      atom-atom interactions
      do is = 1,npart-1
c         print *, 'Atom ', is, ' coords = ', r(1,is), r(2,is), r(3,is)
         do js = is+1,npart
            call lj(is,js,r, rr, box,rc,Uij,rFij)
            Utot = Utot + Uij
            Pkin = Pkin + rFij/(3*box**3)

            if(rr < box/2.d0) then
                  ig = int(rr/delg)
                  if ((ig>0)) then
                     g(ig) = g(ig) + 2
                  endif
               endif
         end do
      end do


      return
      end

*********************************************************
*********************************************************
c              subroutine Lennard-Jones
*********************************************************
*********************************************************
      subroutine lj(is,js,r, rr, box,rc,Uij, rFij)
      implicit double precision(a-h,o-z)
      dimension r(3,1000), rij(3)

      rr2 = 0.d0
      Uij = 0.d0
      rijl = 0.d0
      rFij = 0.d0
      do l = 1,3
         rijl = r(l,js) - r(l,is)
         rij(l) = rijl - box*dnint(rijl/box)
         rr2 = rr2 + rijl**2
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
     &                  U_vals, P_vals, g, dU, dP)
      implicit double precision(a-h,o-z)
      dimension r(3,1000), U_vals(1000), P_vals(1000), g(1000)
      
      Utot = 0.d0
      Pkin = 0.d0
      rho = npart/box**3

      call getValues(npart, nhis, r, box, rc, Utot, Pkin, g)
      


      U_vals(isample) = Utot + dU
      P_vals(isample) = Pkin + rho/beta + dP

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