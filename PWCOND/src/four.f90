!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Modified by M. Pourfath (2025)
! Extended to include f-orbital (l = 3) calculations.
!
! IMPORTANT: Physics invariants and phase conventions are documented in
!            .github/copilot-instructions.md
!            Any modifications MUST preserve these invariants.
!
subroutine four(w0, z0, dz, tblm, taunew, r, rab, betar)
!
! This routine computes the bidimensional fourier transform of the
! beta function. It has been implemented for s, p, d, f-orbitals.
!
!   w0(z,g,m)=1/S * \int w(r) \exp{-ig r_\perp} dr_\perp
!   where w(r) - beta function of the alpha's orbital.
!
!   (see Gradshtein "Tables of integrals")
! For a fixed l it computes w0 for all m.
!
! Physics background:
! - UPF (Unified Pseudopotential Format) pseudopotentials define projectors using
!   complex spherical harmonics Y_l^m(theta,phi) in the Condon-Shortley convention
! - Internally, this routine uses REAL spherical harmonics (cos/sin combinations)
! - Phase conventions: The 2D Fourier transform uses the plane-wave expansion:
!   exp(-i g·r_perp) = sum_m (-i)^m J_m(g r_perp) exp(i m(phi - phi_g))
! - This leads to phase pattern: m=0→+1, m=1→-i, m=2→-1, m=3→+i
!
! The order of spherical harmonics used (real harmonics, PWCOND ordering):
!             s ;
!             p_z, p_{-x}, p_{-y} ;
!             d_{z^2-1}, d_{-xz}, d_{-yz}, d_{x^2-y^2}, d_{xy}
!             f_{z(5z^2-3r^2)}, f_{-x(5z^2-r^2)}, f_{-y(5z^2-r^2)},
!             f_{z(x^2-y^2)}, f_{xyz}, f_{-x(x^2-3y^2)}, f_{-y(3x^2-y^2)}
!
! input:  tblm   -  array characterizing the orbital.
!         taunew -  coordinates and radius of the orbital.
!         z0     -  the initial z
!         dz     -  the slab width
!
! output: w0(z, g, m), where
!                      z0< z <z0+dz
!                      g - 2D g-vector
!
  USE kinds, ONLY: DP
  USE constants, ONLY : tpi, fpi
  USE radial_grids, only : ndmx
  USE cell_base, ONLY : alat, tpiba
  USE cond, ONLY : sarea, nz1, ngper, gper, ninsh, gnsh, ngpsh

implicit none

  integer :: kz, ig, ign, igphi, &
             indexr, iz, lb, ir, nmesh, nmeshs, tblm(4)
  real(DP), parameter :: eps=1.d-8
  complex(DP), parameter :: cim=(0.d0, 1.d0)
  real(DP) :: gn, s1, s2, s3, s4, cs, sn, cs2, sn2, cs3, sn3, rz, dz1, zr, &
                   dr, z0, dz,  bessj, taunew(4), r(ndmx),         &
                   rab(ndmx), betar(ndmx)
  real(DP), allocatable :: x1(:), x2(:), x3(:), x4(:), x5(:), x6(:)
  real(DP), allocatable :: fx1(:), fx2(:), fx3(:), fx4(:), fx5(:), fx6(:), zsl(:)
  complex(DP) :: w0(nz1, ngper, 7)
  complex(DP), allocatable :: wadd(:,:), wadd2(:,:), wadd3(:,:)
  complex(DP) :: t1, t2, t3, t4, t5, t6, t7, wa1, wa2, wa3


  allocate( x1(0:ndmx) )
  allocate( x2(0:ndmx) )
  allocate( x3(0:ndmx) )
  allocate( x4(0:ndmx) )
  allocate( x5(0:ndmx) )
  allocate( x6(0:ndmx) )
  allocate( fx1( nz1 ) )
  allocate( fx2( nz1 ) )
  allocate( fx3( nz1 ) )
  allocate( fx4( nz1 ) )
  allocate( fx5( nz1 ) )
  allocate( fx6( nz1 ) )
  allocate( zsl( nz1) )
  allocate( wadd( nz1, ngper ) )
  allocate( wadd2( nz1, ngper ) )
  allocate( wadd3( nz1, ngper ) )

  lb = tblm(3)
  nmesh=indexr(taunew(4)*alat,ndmx,r)
  dz1=dz/nz1
  zsl(1)=(z0+dz1*0.5d0-taunew(3))*alat
  do kz = 2, nz1
    zsl(kz) = zsl(kz-1)+dz1*alat
  enddo


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
         ! RADIAL INTEGRATION: Compute radial integrals x1..x6 using
         ! Bessel functions J_m(g*r_perp) with appropriate powers of r and r_perp.
         ! These are pure RADIAL quantities - NO phase factors here.
         ! Phase corrections applied later in angular assembly (see lines 250-290).
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
               x1(ir)=betar(ir)*bessj(3,gn*rz)*rz**3/r(ir)**3
               x2(ir)=betar(ir)*bessj(2,gn*rz)*rz**2/r(ir)**3
               x3(ir)=betar(ir)*bessj(1,gn*rz)*rz/r(ir)**3
               x4(ir)=betar(ir)*bessj(1,gn*rz)*rz**3/r(ir)**3
               x5(ir)=betar(ir)*bessj(0,gn*rz)/r(ir)**3
               x6(ir)=betar(ir)*bessj(0,gn*rz)*rz**2/r(ir)**3
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
              x2(iz-1)=(betar(iz)-(betar(iz)-   &
                        betar(iz-1))/dr*zr)/abs(zsl(kz))
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
               x3(iz-1)=(betar(iz)-(betar(iz)-   &
                         betar(iz-1))/dr*zr)/abs(zsl(kz))**2
               x4(iz-1)=betar(iz)-(betar(iz)-       &
                         betar(iz-1))/dr*zr
            endif
            fx3(kz)=fx3(kz)+(x3(iz-1)+x3(iz))*0.5d0*zr
            fx4(kz)=fx4(kz)+(x4(iz-1)+x4(iz))*0.5d0*zr
         elseif (lb.eq.3) then
            fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x3(iz),rab(iz),fx3(kz))
            fx3(kz)=fx3(kz)+x3(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x4(iz),rab(iz),fx4(kz))
            fx4(kz)=fx4(kz)+x4(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x5(iz),rab(iz),fx5(kz))
            call simpson(nmeshs-iz+1,x6(iz),rab(iz),fx6(kz))
            if(iz.eq.1) then
               x5(iz-1)=0.d0
            else
               x5(iz-1)=(betar(iz)-(betar(iz)-betar(iz-1))/dr*zr)/(abs(zsl(kz))**3)
            endif
            x6(iz-1)=0.d0
            fx5(kz)=fx5(kz)+(x5(iz-1)+x5(iz))*0.5d0*zr
            fx6(kz)=fx6(kz)+(x6(iz-1)+x6(iz))*0.5d0*zr
         endif
       else
          fx1(kz)=0.d0
          fx2(kz)=0.d0
          fx3(kz)=0.d0
          fx4(kz)=0.d0
          fx5(kz)=0.d0
          fx6(kz)=0.d0
       endif
     enddo
     ! ANGULAR ASSEMBLY: Combine radial integrals (fx1..fx6) with angular factors
     ! (cs, sn, cs2, sn2, cs3, sn3) to form real spherical harmonics.
     ! This implements the cos/sin combinations that convert complex Y_l^m
     ! to real harmonics in the PWCOND ordering.
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
        cs3 = cs * (4.d0*cs**2 - 3.d0)
        sn3 = sn * (3.d0 - 4.d0*sn**2)

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
               ! f-orbitals (l=3): 7 channels for m = 0, ±1, ±2, ±3
               ! m=0: w0(:,:,1) from fx5 (J_0 integral)
               ! m=1: w0(:,:,2),w0(:,:,3) from fx3,fx4 with cs,sn (J_1 integral, cos/sin)
               ! m=2: w0(:,:,4),w0(:,:,5) from fx2 with cs2,sn2 (J_2 integral, cos/sin)
               ! m=3: w0(:,:,6),w0(:,:,7) from fx1 with cs3,sn3 (J_3 integral, cos/sin)
               ! M-BLOCK CONSISTENCY: Each |m| uses same Bessel J_m, same radial integral,
               ! differs only by cos/sin angular factor. Do NOT change this ordering.
               w0(kz,ig,1)=fx5(kz)
               wadd(kz,ig)=fx6(kz)
               w0(kz,ig,2)=cs*fx3(kz)
               wadd2(kz,ig)=cs*fx4(kz)
               w0(kz,ig,3)=sn*fx3(kz)
               wadd3(kz,ig)=sn*fx4(kz)
               w0(kz,ig,4)=cs2*fx2(kz)
               w0(kz,ig,5)=sn2*fx2(kz)
               w0(kz,ig,6)=cs3*fx1(kz)
               w0(kz,ig,7)=sn3*fx1(kz)
            endif
        enddo
     enddo

  enddo

  ! NORMALIZATION AND PHASE FACTORS: Apply normalization constants (s1..s4) and
  ! complex phase factors (cim = i) to complete the transformation.
  ! Phase pattern: m=0→real, m=1→i, m=2→-1, m=3→-i (PWCOND gauge)
  ! This differs from canonical (-i)^m pattern by global sign for odd m.
  ! These phase factors come from the plane-wave expansion and must be
  ! consistent with the UPF Ylm definitions (see copilot-instructions.md).
  if (lb.eq.0) then
     s1=tpi/sarea/sqrt(fpi)
  elseif (lb.eq.1) then
     s1=tpi/sarea*sqrt(3.d0/fpi)
  elseif (lb.eq.2) then
     s1=-tpi/2.d0/sarea*sqrt(15.d0/fpi)
     s2=tpi/sarea*sqrt(5.d0/tpi/8.d0)
  elseif (lb.eq.3) then
     s1 = tpi/sarea * sqrt(35.d0/(fpi*8.d0))
     s2 = tpi/sarea * sqrt(105.d0/(fpi*4.d0))
     s3 = tpi/sarea * sqrt(21.d0/(fpi*8.d0))
     s4 = tpi/sarea * sqrt(7.d0/(fpi*4.d0))
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
        ! f-orbitals final phase and normalization
        ! Phase factors: m=0→real, m=1→+i, m=2→-1, m=3→-i (PWCOND gauge)
        ! Implementation: m=0 uses s4, m=1 uses cim*s3, m=2 uses -s2, m=3 uses -cim*s1
        ! Note: m=2,3 have negative signs, consistent with PWCOND gauge
        t1=w0(kz,ig,1);wa1=wadd(kz,ig)
        t2=w0(kz,ig,2);wa2=wadd2(kz,ig)
        t3=w0(kz,ig,3);wa3=wadd3(kz,ig)
        t4=w0(kz,ig,4);t5=w0(kz,ig,5)
        t6=w0(kz,ig,6);t7=w0(kz,ig,7)
        w0(kz,ig,1)=s4*(2.d0*zsl(kz)**3*t1-3.d0*zsl(kz)*wa1)
        w0(kz,ig,2)=cim*s3*(4.d0*zsl(kz)**2*t2-wa2)
        w0(kz,ig,3)=cim*s3*(4.d0*zsl(kz)**2*t3-wa3)
        w0(kz,ig,4)=-s2*zsl(kz)*t4
        w0(kz,ig,5)=-s2*zsl(kz)*t5
        w0(kz,ig,6)=-cim*s1*t6
        w0(kz,ig,7)=-cim*s1*t7
      endif
    enddo
  enddo

  ! DEBUG OUTPUT: Write w0 data for phase checking with tools/check_pwcond_phases.py
  ! Set environment variable PWCOND_DEBUG_PHASES=1 to enable this output
  call write_w0_debug(w0, lb, nz1, ngper)

  deallocate(x1)
  deallocate(x2)
  deallocate(x3)
  deallocate(x4)
  deallocate(x5)
  deallocate(x6)
  deallocate(fx1)
  deallocate(fx2)
  deallocate(fx3)
  deallocate(fx4)
  deallocate(fx5)
  deallocate(fx6)
  deallocate(zsl)
  deallocate(wadd)
  deallocate(wadd2)
  deallocate(wadd3)

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

subroutine write_w0_debug(w0, lb, nz1, ngper)
  !
  ! Debug output subroutine for phase checking
  ! Appends w0(kz,ig,m) data to 'w0_debug.dat' for analysis with
  ! tools/check_pwcond_phases.py
  !
  ! Enable by setting environment variable: export PWCOND_DEBUG_PHASES=1
  ! Data from multiple iterations is appended to allow comparison across orbitals.
  !
  USE kinds, ONLY: DP
  implicit none
  
  integer, intent(in) :: lb, nz1, ngper
  complex(DP), intent(in) :: w0(nz1, ngper, 7)
  
  integer :: kz, ig, m, nm, debug_unit, ios
  integer, external :: find_free_unit
  character(len=10) :: env_val
  logical :: debug_enabled, file_exists
  
  ! Check if debug output is enabled via environment variable
  call get_environment_variable('PWCOND_DEBUG_PHASES', env_val)
  debug_enabled = (trim(env_val) == '1')
  
  if (.not. debug_enabled) return
  
  ! Determine number of m channels based on lb (angular momentum quantum number)
  ! Number of m channels: s(l=0)→1, p(l=1)→3, d(l=2)→5, f(l=3)→7
  if (lb.eq.0) then
     nm = 1  ! s-orbital
  elseif (lb.eq.1) then
     nm = 3  ! p-orbitals
  elseif (lb.eq.2) then
     nm = 5  ! d-orbitals
  elseif (lb.eq.3) then
     nm = 7  ! f-orbitals
  else
     return  ! Unknown lb value
  endif
  
  ! Check if file exists to determine whether to write header
  inquire(file='w0_debug.dat', exist=file_exists)
  
  ! Find free unit and open debug file for appending
  call find_free_unit(debug_unit)
  if (file_exists) then
     ! Append to existing file
     open(unit=debug_unit, file='w0_debug.dat', status='old', &
          position='append', action='write', iostat=ios)
  else
     ! Create new file and write header
     open(unit=debug_unit, file='w0_debug.dat', status='new', &
          action='write', iostat=ios)
     if (ios == 0) then
        write(debug_unit, '(A)') '# w0 debug output for phase checking'
        write(debug_unit, '(A)') '# Format: kz ig m Re(w0) Im(w0)'
        write(debug_unit, '(A)') '# All indices are zero-based for Python compatibility'
        write(debug_unit, '(A)') '# Each block below corresponds to one call to four()'
     endif
  endif
  
  if (ios /= 0) return  ! Failed to open file
  
  ! Write separator and metadata for this iteration
  write(debug_unit, '(A)') '#'
  write(debug_unit, '(A,I2)') '# lb (angular momentum) = ', lb
  write(debug_unit, '(A,I4)') '# nz1 = ', nz1
  write(debug_unit, '(A,I4)') '# ngper = ', ngper
  write(debug_unit, '(A,I2)') '# nm (number of m channels) = ', nm
  
  ! Write a subset of data (to keep file size manageable)
  ! Write first few kz and ig values, all m channels
  do kz = 1, min(3, nz1)
     do ig = 1, min(5, ngper)
        do m = 1, nm
           write(debug_unit, '(3I6,2ES24.15)') kz-1, ig-1, m-1, &
                real(w0(kz,ig,m)), aimag(w0(kz,ig,m))
        enddo
     enddo
  enddo
  
  close(debug_unit)
  
  write(*,'(A)') '*** PWCOND DEBUG: w0 data written to w0_debug.dat'
  write(*,'(A,I2,A,I2,A)') '***              lb=', lb, ', nm=', nm, ' channels'
  
end subroutine write_w0_debug
