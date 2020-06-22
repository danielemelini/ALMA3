! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'write_log'
! Writes various information to the log file opened on unit 99
!
! Initial version DM February 24, 2020
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine write_log(imode)
use general_parameters
use rheological_profile
implicit none
!
 integer :: imode
 integer :: i
 integer :: onefound, densiv
!
 character(10) :: lab(0:6) = (/ 'Fluid     ', 'Elastic   ', 'Maxwell   ', 'Newton    ', &
                                'Kelvin    ', 'Burgers   ', 'Andrade   ' /)
!
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! | MODE=0 : General parameters                         |
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
 if( imode==0 ) then
!
 write(99,*) ' '
 write(99,*) ' -------------'
 write(99,*) ' ALMA log file'
 write(99,*) ' -------------'
!
 write(99,*) ' '
 write(99,*) ' Number of significant digits set to: ',nd
 write(99,*) ' '
 write(99,*) ' Order of the Salzer accelerated sequence: ',order
 write(99,*) ' '
 write(99,*) ' Harmonic degrees (min/max/step): ',lmin, lmax, lstep
 write(99,*) ' '
! 
 if (iload==0) write(99,*) ' Requested Love numbers are of the TIDAL type'
 if (iload==1) write(99,*) ' Requested Love numbers are of the LOADING type'
!
!
 end if
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! | MODE=1 : Un-normalized rheological profile          |
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
 if( imode==1 ) then
!
 write(99,*) ' '
 write(99,*) ' Radii of the interfaces from bottom to top, (km)'
!
 i=0
 write(99,*) ' ', i, ':', to_sp(r(i))/1000, '(Planet center)'
!
 i=1
 write(99,*) ' ', i, ':', to_sp(r(i))/1000, '(CMB)'	
! 
 do i=2, nla+1
     write(99,*) ' ', i, ':', to_sp(r(i))/1000
 enddo
!
 i=nla+2
 write(99,*) ' ', i, ':', to_sp(r(i))/1000, '(External surface)' 
!
!
 write(99,*) ' '
 write(99,*) ' Rigidity of the layers from bottom to top, (*10^11 Pa)' 
!
 i=0
 write(99,*) ' ', i, ':', to_sp(mu(i))/1.e11, '(Core)'
!
 do i=1, nla
     write(99,*) ' ', i, ':', to_sp(mu(i))/1.e11   
 enddo
! 
 i=nla+1
 write(99,*) ' ', i, ':', to_sp(mu(i))/1.e11, '(Lithosphere)'	
!
!
 write(99,*) ' '
 write(99,*) ' Density of the layers from bottom to top, (kg/m^3)' 
!
 i=0
 write(99,*) ' ', i, ':', to_sp(rho(i)), '(Core)'
!
 do i=1, nla
     write(99,*) ' ', i, ':', to_sp(rho(i))    
 enddo	
!
 i=nla+1
 write(99,*) ' ', i, ':', to_sp(rho(i)), '(Lithosphere)'
!
!
 write(99,*) ' '
 densiv=0
 onefound=0
 do i=0,nla
      if(rho(i).lt.rho(i+1)) then 
      onefound=1 
      densiv=1
      write(99,*) ' Density inversion for layers ', i, i+1
      densiv=0
      endif		
 enddo	 
!
 if(onefound==0) write(99,*) ' NO density inversions found!'
!
!
 write(99,*) ' '
 write(99,*) ' Viscosity of the layers from bottom to top, (Pa.s)'	
 i=0
 write(99,*) ' ', i, ':', to_dp(eta(i)), '(Core)'
 do i=1, NLA
     write(99,*) ' ', i, ':', to_dp(eta(i))    
 enddo	
 i=NLA+1
 write(99,*) ' ', i, ':', to_dp(eta(i)), ' (Lithosphere)'
!
!
 write(99,*) ' '
 write(99,*) ' Rheology of the layers from bottom to top'	
 i=0
 write(99,*) ' ', i, ':', lab(irheol(i)), ' (Core)'
 do i=1, NLA
     write(99,*) ' ', i, ':', lab(irheol(i))    
 enddo	
 i=NLA+1
 write(99,*) ' ', i, ':', lab(irheol(i)), ' (Lithosphere)'

!
!
 write(99,*) ' '
 write(99,*) ' Mass of the planet (kg) = ', to_dp(mass) 
!
 write(99,*) ' '
 write(99,*) ' Gravity at the interfaces from bottom to top, (m/s^2)'
 i=0
 write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i)), '(Earth center)'
 i=1
 write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i)), '(CMB)'	
 do i=2, NLA+1
     write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i))
 enddo
 i=NLA+2
 write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i)), '(External surface)'	 
!
!
 end if
!
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! | MODE=2 : Normalized rheological profile             |
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
 if( imode==2 ) then
!
!
 write(99,*) ' '
 write(99,*) ' Reference radius (m) =',       to_sp(r0) 
 write(99,*) ' Reference density (kg/m^3) =', to_sp(rho0) 
 write(99,*) ' Reference rigidity (Pa) =',    to_sp(mu0) 
 write(99,*) ' Reference time (s) =',         to_sp(t0) 
 write(99,*) ' Reference mass (kg) =',        to_sp(mass0) 
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Radii of the interfaces from bottom to top'
 i=0
 write(99,*) ' ', i, ':', to_sp(r(i)), '(Planet center)'
 i=1
 write(99,*) ' ', i, ':', to_sp(r(i)), '(CMB)'	
 do i=2, nla+1
     write(99,*) ' ', i, ':', to_sp(r(i))
 enddo
 i=nla+2
 write(99,*) ' ', i, ':', to_sp(r(i)), '(External surface)' 
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Rigidity of the layers from bottom to top' 
 i=0
 write(99,*) ' ', i, ':', to_sp(mu(i)), '(Core)'
 do i=1, nla
     write(99,*) ' ', i, ':', to_sp(mu(i))  
 enddo
 i=nla+1
 write(99,*) ' ', i, ':', to_sp(mu(i)), '(Lithosphere)'	
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Density of the layers from bottom to top' 
 i=0
 write(99,*) ' ', i, ':', to_sp(mu(i)), '(Core)'
 do i=1, NLA
     write(99,*) ' ', i, ':', to_sp(mu(i))    
 enddo	
 i=NLA+1
 write(99,*) ' ', i, ':', to_sp(mu(i)), '(Lithosphere)'
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Viscosity of the layers from bottom to top, (Pa.s)'	
 i=0
 write(99,*) ' ', i, ':', to_dp(eta(i)), '(Core)'
 do i=1, NLA
     write(99,*) ' ', i, ':', to_dp(eta(i))    
 enddo	
 i=NLA+1
 write(99,*) ' ', i, ':', to_dp(eta(i)), '(Lithosphere)'
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Mass of the Earth = ', to_dp(mass) 
 write(99,*) ' '
 write(99,*) ' NORMALIZED Newton constant = ', to_dp(G) 
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Gravity at the interfaces from bottom to top'
 i=0
 write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i)), '(Planet center)'
 i=1
 write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i)), '(CMB)'	
 do i=2, NLA+1
     write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i))
 enddo
 i=NLA+2
 write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i)), '(External surface)'	 
!
!
 end if
!
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! | MODE=3 : Time steps                                 |
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
 if (imode==3) then
!
 write(99,*) ' '
!
    if( itime==0 ) then
        write(99,*) ' Time scale is LINEAR' 
    elseif( itime==1 ) then
        write(99,*) ' Time scale is LOGARITHMIC' 
    end if
    write(99,*) ' '
    write(99,*) ' Time points (kyr):' 
    do i=1,p+1
         write(99,*) ' ', i, ':', to_sp(t(i))
    end do
 
end if
!
!
!
!
!
!
end subroutine write_log