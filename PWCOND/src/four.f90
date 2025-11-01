subroutine four(w0, z0, dz, tblm, taunew, r, rab, betar)
!
! This routine computes the bidimensional fourier transform of the
! beta function. It has been implemented for s, p, d-orbitals.
!
!   w0(z,g,m)=1/S * \int w(r) \exp{-ig r_\perp} dr_\perp
!   where w(r) - beta function of the alpha's orbital.
!
!   (see Gradshtein "Tables of integrals")
! For a fixed l it computes w0 for all m (for l<=2). For l=3, only the requested m.
!
! The order of spherical harmonics used:
!             s ;
!             p_z, p_{-x}, p_{-y} ;
!             d_{z^2-1}, d_{-xz}, d_{-yz}, d_{x^2-y^2}, d_{xy}
!
  USE kinds, ONLY: DP
  USE constants, ONLY : tpi, fpi
  USE radial_grids, only : ndmx
  USE cell_base, ONLY : alat, tpiba
  USE cond, ONLY : sarea, nz1, ngper, gper, ninsh, gnsh, ngpsh
  USE io_global, ONLY: stdout, ionode

implicit none

  integer :: kz, ig, ign, igphi, &
             indexr, iz, lb, ir, nmesh, nmeshs, tblm(4), m_in, msel
  real(DP), parameter :: eps=1.d-8
  complex(DP), parameter :: cim=(0.d0, 1.d0)
  real(DP) :: gn, s1, s2, s3, s4, cs, sn, cs2, sn2, cs3, sn3, rz, dz1, zr, &
                   dr, z0, dz,  bessj, taunew(4), r(ndmx),         &
                   rab(ndmx), betar(ndmx), beta_boundary
  real(DP), allocatable :: x1(:), x2(:), x3(:), x4(:)
  real(DP), allocatable :: fx1(:), fx2(:), fx3(:), fx4(:), zsl(:)
  complex(DP) :: w0(nz1, ngper, 5)
  complex(DP), allocatable :: wadd(:,:)

  ! f-only storage
  real(DP), allocatable :: x5(:), x6(:)
  real(DP), allocatable :: fx5(:), fx6(:)
  complex(DP), allocatable :: wadd2(:,:), wadd3(:,:)

  ! Optional incoming f-order mapping (UPF -> analytical):
  ! Analytical internal order for f (1..7): m=0, cosφ, sinφ, cos2φ, sin2φ, cos3φ, sin3φ
  ! Set map_f_in(k_in) = k_target in [1..7]. Default identity.
  integer, parameter :: map_f_in(7) = (/ 1,2,3,4,5,6,7 /)

  ! Debug controls (env-driven, ionode-only)
  logical :: ldbg_f
  integer :: debug_max, debug_hits, istat, lstr, slot_sel
  character(len=32) :: env

  ldbg_f = .false.
  debug_max = 5
  debug_hits = 0
  call get_environment_variable('PWCOND_DEBUG_F', env, length=lstr, status=istat)
  if (istat==0 .and. lstr>0) then
     if (env(1:1)=='1' .or. env(1:1)=='T' .or. env(1:1)=='t' .or. env(1:1)=='Y' .or. env(1:1)=='y') ldbg_f=.true.
  end if
  call get_environment_variable('PWCOND_DEBUG_F_MAX', env, length=lstr, status=istat)
  if (istat==0 .and. lstr>0) then
     read(env,*,err=10,end=10) debug_max
  end if
10 continue

  allocate( x1(0:ndmx) )
  allocate( x2(0:ndmx) )
  allocate( x3(0:ndmx) )
  allocate( x4(0:ndmx) )
  allocate( fx1( nz1 ) )
  allocate( fx2( nz1 ) )
  allocate( fx3( nz1 ) )
  allocate( fx4( nz1 ) )
  allocate( zsl( nz1) )
  allocate( wadd( nz1, ngper ) )

  lb = tblm(3)
  m_in = tblm(4)
  msel = m_in
  if (lb.eq.3) then
     if (m_in<1 .or. m_in>7) call errore('four','invalid m for l=3',1)
     msel = map_f_in(m_in)
     allocate( x5(0:ndmx), x6(0:ndmx) )
     allocate( fx5(nz1), fx6(nz1) )
     allocate( wadd2(nz1,ngper), wadd3(nz1,ngper) )
     ! Compute the selected output slot for logging
     select case (msel)
     case (1); slot_sel=1
     case (2); slot_sel=2
     case (3); slot_sel=3
     case (4); slot_sel=4
     case (5); slot_sel=5
     case (6); slot_sel=2
     case (7); slot_sel=3
     end select
     if (ionode .and. ldbg_f .and. debug_hits<debug_max) then
        write(stdout,'(A,3(I0,1X),A,I0)') 'four.f90: l=3 m_in=',m_in,' msel=',msel,' -> slot=',slot_sel
        debug_hits = debug_hits + 1
     end if
  end if

  nmesh=indexr(taunew(4)*alat,ndmx,r)
  dz1=dz/nz1
  zsl(1)=(z0+dz1*0.5d0-taunew(3))*alat
  do kz = 2, nz1
    zsl(kz) = zsl(kz-1)+dz1*alat
  enddo

  ! Zero outputs for safety
  w0(:,:,:) = (0.d0,0.d0)
  fx1(:) = 0.d0 ; fx2(:) = 0.d0 ; fx3(:) = 0.d0 ; fx4(:) = 0.d0
  if (lb.eq.3) then
     fx5(:) = 0.d0 ; fx6(:) = 0.d0
  end if
  wadd(:,:) = (0.d0,0.d0)
  if (lb.eq.3) then
     wadd2(:,:) = (0.d0,0.d0)
     wadd3(:,:) = (0.d0,0.d0)
  end if

  ig=0
  do ign=1, ngpsh

     gn=gnsh(ign)
     do kz=1, nz1
       if (abs(zsl(kz))+eps.le.taunew(4)*alat) then
         iz=indexr(zsl(kz),nmesh,r)
         if ((nmesh-iz)/2*2.eq.nmesh-iz) then
            nmeshs=nmesh
         else
            nmeshs=nmesh+1
         endif
         do ir=iz, nmeshs
            rz=sqrt(r(ir)**2-zsl(kz)**2)
            if (lb.eq.0) then
               x1(ir)=betar(ir)*bessj(0,gn*rz)
            elseif (lb.eq.1) then
               x1(ir)=betar(ir)*bessj(1,gn*rz)/r(ir)*rz
               x2(ir)=betar(ir)*bessj(0,gn*rz)/r(ir)
            elseif (lb.eq.2) then
               x1(ir)=betar(ir)*bessj(2,gn*rz)*rz**2/r(ir)**2
               x2(ir)=betar(ir)*bessj(1,gn*rz)*rz/r(ir)**2
               x3(ir)=betar(ir)*bessj(0,gn*rz)/r(ir)**2
               x4(ir)=betar(ir)*bessj(0,gn*rz)
            elseif (lb.eq.3) then
               ! f kernels
               x1(ir)=betar(ir)*bessj(3,gn*rz)*rz**3/r(ir)**3  ! ρ^3/r^3
               x2(ir)=betar(ir)*bessj(2,gn*rz)*rz**2/r(ir)**3  ! ρ^2/r^3
               x3(ir)=betar(ir)*bessj(1,gn*rz)*rz/r(ir)**3     ! ρ/r^3
               x4(ir)=betar(ir)*bessj(1,gn*rz)*rz**3/r(ir)**3  ! ρ^3/r^3
               x5(ir)=betar(ir)*bessj(0,gn*rz)/r(ir)**3        ! 1/r^3
               x6(ir)=betar(ir)*bessj(0,gn*rz)*rz**2/r(ir)**3  ! ρ^2/r^3
            else
               call errore ('four','ls not programmed ',1)
            endif
         enddo
         call simpson(nmeshs-iz+1,x1(iz),rab(iz),fx1(kz))
         if (iz.eq.1) then
            dr=r(iz)
         else
            dr=r(iz)-r(iz-1)
         endif
         zr=r(iz)-abs(zsl(kz))
         if (lb.eq.0) then
            if (iz.eq.1) then
               x1(iz-1)=betar(iz)-betar(iz)/dr*zr
            else
               x1(iz-1)=betar(iz)-(betar(iz)-betar(iz-1))/dr*zr
            endif
            fx1(kz)=fx1(kz)+(x1(iz-1)+x1(iz))*0.5d0*zr
         else
            fx1(kz)=fx1(kz)+x1(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x2(iz),rab(iz),fx2(kz))
         endif
         if (lb.eq.1) then
            if(iz.eq.1) then
              x2(iz-1)=0.d0
            else
              x2(iz-1)=(betar(iz)-(betar(iz)-betar(iz-1))/dr*zr)/abs(zsl(kz))
            endif
            fx2(kz)=fx2(kz)+(x2(iz-1)+x2(iz))*0.5d0*zr
         endif
         if (lb.eq.2) then
            fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x3(iz),rab(iz),fx3(kz))
            call simpson(nmeshs-iz+1,x4(iz),rab(iz),fx4(kz))
            if(iz.eq.1) then
               x3(iz-1)=0.d0
               x4(iz-1)=0.d0
            else
               x3(iz-1)=(betar(iz)-(betar(iz)-betar(iz-1))/dr*zr)/abs(zsl(kz))**2
               x4(iz-1)=betar(iz)-(betar(iz)-betar(iz-1))/dr*zr
            endif
            fx3(kz)=fx3(kz)+(x3(iz-1)+x3(iz))*0.5d0*zr
            fx4(kz)=fx4(kz)+(x4(iz-1)+x4(iz))*0.5d0*zr
         elseif (lb.eq.3) then
            ! f: complete kernels with boundary treatment
            fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x3(iz),rab(iz),fx3(kz))
            fx3(kz)=fx3(kz)+x3(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x4(iz),rab(iz),fx4(kz))
            fx4(kz)=fx4(kz)+x4(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x5(iz),rab(iz),fx5(kz))
            call simpson(nmeshs-iz+1,x6(iz),rab(iz),fx6(kz))
            ! ρ=0 boundary: only x5 survives; extrapolate it
            if(iz.eq.1) then
               beta_boundary = 0.d0
            else
               beta_boundary = betar(iz)-(betar(iz)-betar(iz-1))/dr*zr
            endif
            x5(iz-1) = 0.d0
            x6(iz-1) = 0.d0
            if (iz.gt.1) x5(iz-1) = beta_boundary / (abs(zsl(kz))**3)
            fx5(kz)=fx5(kz)+(x5(iz-1)+x5(iz))*0.5d0*zr
            fx6(kz)=fx6(kz)+(x6(iz-1)+x6(iz))*0.5d0*zr
         endif
       else
          fx1(kz)=0.d0
          fx2(kz)=0.d0
          fx3(kz)=0.d0
          fx4(kz)=0.d0
          if (lb.eq.3) then
             fx5(kz)=0.d0
             fx6(kz)=0.d0
          endif
       endif
     enddo

     do igphi=1, ninsh(ign)
        ig=ig+1
        if (gn.gt.eps) then
          cs=gper(1,ig)*tpiba/gn
          sn=gper(2,ig)*tpiba/gn
        else
          cs=0.d0
          sn=0.d0
        endif
        cs2=cs**2-sn**2
        sn2=2*cs*sn
        if (lb.eq.3) then
           cs3 = cs * (4.d0*cs**2 - 3.d0)
           sn3 = sn * (3.d0 - 4.d0*sn**2)
        endif

        do kz=1, nz1
            if (lb.eq.0) then
               w0(kz,ig,1)=fx1(kz)
            elseif (lb.eq.1) then
               w0(kz,ig,2)=cs*fx1(kz)
               w0(kz,ig,1)=fx2(kz)
               w0(kz,ig,3)=sn*fx1(kz)
            elseif (lb.eq.2) then
               w0(kz,ig,5)=sn2*fx1(kz)
               w0(kz,ig,2)=cs*fx2(kz)
               w0(kz,ig,1)=fx3(kz)
               w0(kz,ig,3)=sn*fx2(kz)
               w0(kz,ig,4)=cs2*fx1(kz)
               wadd(kz,ig)=fx4(kz)
            elseif (lb.eq.3) then
               ! For f: fill only the slot of the selected m (others remain zero)
               select case (msel)
               case (1)  ! m=0 -> slot 1
                  w0(kz,ig,1)=fx5(kz)       ! t1
                  wadd(kz,ig)=fx6(kz)       ! wa1
               case (2)  ! cos φ (|m|=1) -> slot 2
                  w0(kz,ig,2)=cs*fx3(kz)    ! t2
                  wadd2(kz,ig)=cs*fx4(kz)   ! wa2
               case (3)  ! sin φ (|m|=1) -> slot 3
                  w0(kz,ig,3)=sn*fx3(kz)    ! t3
                  wadd3(kz,ig)=sn*fx4(kz)   ! wa3
               case (4)  ! cos 2φ (|m|=2) -> slot 4
                  w0(kz,ig,4)=cs2*fx2(kz)   ! t4
               case (5)  ! sin 2φ (|m|=2) -> slot 5
                  w0(kz,ig,5)=sn2*fx2(kz)   ! t5
               case (6)  ! cos 3φ (|m|=3) -> slot 2
                  w0(kz,ig,2)=cs3*fx1(kz)   ! t6
               case (7)  ! sin 3φ (|m|=3) -> slot 3
                  w0(kz,ig,3)=sn3*fx1(kz)   ! t7
               end select
               if (ionode .and. ldbg_f .and. debug_hits<debug_max) then
                  write(stdout,'(A,2I0,A,2I0,A,I0)') 'four.f90: f fill kz=',kz,' ig=',ig,' l=',lb,' msel=',msel,' slot=',slot_sel
                  debug_hits = debug_hits + 1
               end if
            endif
        enddo
     enddo

  enddo

  if (lb.eq.0) then
     s1=tpi/sarea/sqrt(fpi)
  elseif (lb.eq.1) then
     s1=tpi/sarea*sqrt(3.d0/fpi)
  elseif (lb.eq.2) then
     s1=-tpi/2.d0/sarea*sqrt(15.d0/fpi)
     s2=tpi/sarea*sqrt(5.d0/tpi/8.d0)
  elseif (lb.eq.3) then
     s1 = tpi/sarea * sqrt(35.d0/(fpi*4.d0))  ! |m|=3
     s2 = tpi/sarea * sqrt(105.d0/(fpi*4.d0)) ! |m|=2
     s3 = tpi/sarea * sqrt(21.d0/(fpi*8.d0))  ! |m|=1
     s4 = tpi/sarea * sqrt(7.d0/(fpi*4.d0))   ! m=0
     if (ionode .and. ldbg_f .and. debug_hits<debug_max) then
        write(stdout,'(A,1X,4(ES12.4,1X))') 'four.f90: f norms s1(|m|=3), s2(|m|=2), s3(|m|=1), s4(m=0)=', s1, s2, s3, s4
        debug_hits = debug_hits + 1
     end if
  endif
  do ig=1, ngper
    do kz=1, nz1
      if (lb.eq.0) then
        w0(kz,ig,1)=s1*w0(kz,ig,1)
      elseif (lb.eq.1) then
        w0(kz,ig,2)=cim*s1*w0(kz,ig,2)
        w0(kz,ig,1)=s1*zsl(kz)*w0(kz,ig,1)
        w0(kz,ig,3)=cim*s1*w0(kz,ig,3)
      elseif (lb.eq.2) then
        w0(kz,ig,5)=s1*w0(kz,ig,5)
        w0(kz,ig,2)=-2.d0*cim*s1*zsl(kz)*w0(kz,ig,2)
        w0(kz,ig,1)=3.d0*zsl(kz)**2*s2*w0(kz,ig,1)-s2*wadd(kz,ig)
        w0(kz,ig,3)=-2.d0*cim*s1*zsl(kz)*w0(kz,ig,3)
        w0(kz,ig,4)=s1*w0(kz,ig,4)
      elseif (lb.eq.3) then
        select case (msel)
        case (1)  ! m=0 -> slot 1
          w0(kz,ig,1)=s4*(2.d0*zsl(kz)**3*w0(kz,ig,1)-3.d0*zsl(kz)*wadd(kz,ig))
        case (2)  ! cos φ (|m|=1) -> slot 2
          w0(kz,ig,2)=cim*s3*(4.d0*zsl(kz)**2*w0(kz,ig,2)-wadd2(kz,ig))
        case (3)  ! sin φ (|m|=1) -> slot 3
          w0(kz,ig,3)=cim*s3*(4.d0*zsl(kz)**2*w0(kz,ig,3)-wadd3(kz,ig))
        case (4)  ! cos 2φ (|m|=2) -> slot 4
          w0(kz,ig,4)=-s2*zsl(kz)*w0(kz,ig,4)
        case (5)  ! sin 2φ (|m|=2) -> slot 5
          w0(kz,ig,5)=-s2*zsl(kz)*w0(kz,ig,5)
        case (6)  ! cos 3φ (|m|=3) -> slot 2
          w0(kz,ig,2)=cim*s1*w0(kz,ig,2)
        case (7)  ! sin 3φ (|m|=3) -> slot 3
          w0(kz,ig,3)=cim*s1*w0(kz,ig,3)
        end select
        if (ionode .and. ldbg_f .and. debug_hits<debug_max) then
           write(stdout,'(A,2I0,A,1P,ES12.4)') 'four.f90: f final kz=',kz,' ig=',ig,' |w0|=',abs(w0(kz,ig,slot_sel))
           debug_hits = debug_hits + 1
        end if
      endif
    enddo
  enddo

  deallocate(x1)
  deallocate(x2)
  deallocate(x3)
  deallocate(x4)
  deallocate(fx1)
  deallocate(fx2)
  deallocate(fx3)
  deallocate(fx4)
  deallocate(zsl)
  deallocate(wadd)
  if (lb.eq.3) then
     deallocate(x5, x6)
     deallocate(fx5, fx6)
     deallocate(wadd2, wadd3)
  endif

  return
end subroutine four

function indexr(zz, ndim, r)
  USE kinds, only : DP
  implicit none

  integer :: iz, ndim, indexr
  real(DP) :: zz, r(ndim)
!
!     abs(zz)<r(indexr)
!
  iz = 1
  do while(r(iz).le.abs(zz)+1.d-10)
    iz=iz+1
  enddo
  indexr=iz
  return
end function indexr
