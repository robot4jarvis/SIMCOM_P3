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
 
      dimension r(3,1000)

c     2. Reading data and computing related quantities

      open(1,file='leap-lj.data',status='old')
      read(1,*) ncycles
      read(1,*) npart
      read(1,*) tref
      read(1,*) pref 
      read(1,*) sigma,epsil
      read(1,*) deltax 

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
c     Opening files to write results
c
      open(3,file='energy-leap.dat',status='unknown')
      open(4,file='temp-leap.dat',status='unknown')
c
c     4. Change to reduced units

      call reduced(npart,r,box,sigma)

c     5. Start the loop to generate new configurations 


      isuccess = 0
      do icycle = 1,ncycles
         call moveparticle(npart, r, rc, box, deltax, beta, isuccess)
c         call boundaryConds(npart, r, box)

      enddo
      print *, 'Success rate: ', int(100*isuccess/ncycles), "%"


      
      close(3)
      close(4)

c     6. Saving last configuration in A and A/ps 

      open(11,file='newconf.data',status='unknown')
      do is = 1,npart
         write(11,*) (r(l,is)*sigma,l=1,3)
         write(11,*) ! Write empty line
      end do         
      write(11,*) box*sigma
      close(11)

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
      dimension r0(3,1000), rn(3,1000) 
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
         rn(l,iSel) = r0(l,iSel) + deltax*(rand()-0.d5)
      enddo

      ! Energy difference
c      call totEnergy(npart, r0, box, rc, Utot0)
c      call totEnergy(npart, rn, box, rc, Utotn)
      call deltaEnergy(npart, r0, rn, iSel, box, rc, deltaU)

c      Udelta = Utotn - Utot0
      acc = exp(-beta*deltaU)

      if(rand() >= acc) then
         do i = 1, npart
            do l = 1,3
               r0(l,i) = rn(l,i)
            enddo
         enddo
         isuccess = isuccess + 1
      endif ! We do nothing

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
         if (js == iSel)  cycle
         call lj(iSel, js, r0, box, rc, U0)
         call lj(iSel, js, rn, box, rc, Un)
         deltaU = deltaU + Un - U0
      end do

      return
      end

*********************************************************
*********************************************************
c              subroutine totEnergy
*********************************************************
*********************************************************
      subroutine totEnergy(npart, r, box, rc, Utot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
c     This subroutine calculates the energy of a configuration


      Utot = 0.d0 

c      atom-atom interactions
      do is = 1,npart-1
c         print *, 'Atom ', is, ' coords = ', r(1,is), r(2,is), r(3,is)
         do js = is+1,npart
            call lj(is,js,r,box,rc,Uij)
            Utot = Utot + Uij
         end do
      end do


      return
      end

*********************************************************
*********************************************************
c              subroutine Lennard-Jones
*********************************************************
*********************************************************
      subroutine lj(is,js,r,box,rc,Uij)
      implicit double precision(a-h,o-z)
      dimension r(3,1000), rij(3)

      rr2 = 0.d0
      Uij = 0.d0
      rijl = 0
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
            if(r(l,i) < 0 ) r(l,i) = r(l,i) + box
            if(r(l,i) > 0 ) r(l,i) = r(l,i) - box

         enddo
      end do

      return
      end

