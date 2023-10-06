! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'write_log'
! Writes various information to the log file opened on unit 99
!
! Initial version DM February 24, 2020
! Modified by DM Aug 11, 2021 - added some messages in the MODE=0 section
! Modified Oct  5, 2023       - changed the numbering scheme for the layers
!                               fixed some bugs and wrong labels
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
 if (itype==0) write(99,*) ' REAL (time-domain) LNs will be computed'
 if (itype==1) write(99,*) ' COMPLEX (frequency-domain) LNs will be computed'
 if (itype==2) write(99,*) ' REAL (time-domain) DERIVATIVES of LNs will be computed'
!
 if (itype==0) then
    if (ihist==1) write(99,*) ' A Heaviside load time-history is assumed'
    if (ihist==2) write(99,*) ' A linear ramp load time-history is assumed'
    if (ihist==2) write(99,*) ' (loading phase length is ',tau,' kyrs)'
 end if
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
 do i=0,nla
    if(i==0) then
       write(99,*) ' ', i, ':', to_sp(r(i))/1000, '(Planet center)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp(r(i))/1000, '(External surface)' 
    elseif(i==1) then
       write(99,*) ' ', i, ':', to_sp(r(i))/1000, '(CMB)'
    else
       write(99,*) ' ', i, ':', to_sp(r(i))/1000
    end if
 end do
!
!
!
 write(99,*) ' '
 write(99,*) ' Rigidity of the layers from bottom to top, (*10^11 Pa)' 
!
 do i=1,nla
    if(i==1) then
       write(99,*) ' ', i, ':', to_sp(mu(i))/1.e11, '(Core)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp(mu(i))/1.e11, '(Lithosphere)'
    else
       write(99,*) ' ', i, ':', to_sp(mu(i))/1.e11   
    end if
 end do
! 
!
!
 write(99,*) ' '
 write(99,*) ' Density of the layers from bottom to top, (kg/m^3)' 
!
 do i=1,nla
    if(i==1) then
       write(99,*) ' ', i, ':', to_sp(rho(i)), '(Core)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp(rho(i)), '(Lithosphere)'
    else
       write(99,*) ' ', i, ':', to_sp(rho(i))    
    end if
 end do
!
!
!
 write(99,*) ' '
 densiv=0
 onefound=0
 do i=1,nla-1
      if(rho(i).lt.rho(i+1)) then 
      onefound=1 
      densiv=1
      write(99,*) ' Density inversion for layers ', i, i+1
      densiv=0
      endif
 end do
!
 if(onefound==0) write(99,*) ' NO density inversions found!'
!
!
 write(99,*) ' '
 write(99,*) ' Viscosity of the layers from bottom to top, (Pa.s)'	
!
 do i=1,nla
    if(i==1) then
       write(99,*) ' ', i, ':', to_dp(eta(i)), '(Core)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_dp(eta(i)), ' (Lithosphere)'
    else
       write(99,*) ' ', i, ':', to_dp(eta(i))    
    end if
 end do
!
!
!
 write(99,*) ' '
 write(99,*) ' Rheology of the layers from bottom to top'
!
 do i=1,nla
    if(i==1) then
       write(99,*) ' ', i, ':', lab(irheol(i)), ' (Core)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', lab(irheol(i)), ' (Lithosphere)'
    else
       write(99,*) ' ', i, ':', lab(irheol(i))    
    end if
 end do
!
!
!
 write(99,*) ' '
 write(99,*) ' Mass of the planet (kg) = ', to_dp(mass) 
!
 write(99,*) ' '
 write(99,*) ' Gravity at the interfaces from bottom to top, (m/s^2)'
!
 do i=0,nla
    if(i==0) then
       write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i)), '(Planet center)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i)), '(External surface)'	 
    elseif(i==1) then
       write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i)), '(CMB)'
    else
       write(99,*) ' ', i, ':', to_sp(r(i))/1000, to_sp(gra(i))
    end if
 end do
!
!
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
!
 do i=0,nla
    if(i==0) then
       write(99,*) ' ', i, ':', to_sp(r(i)), '(Planet center)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp(r(i)), '(External surface)' 
    elseif(i==1) then
       write(99,*) ' ', i, ':', to_sp(r(i)), '(CMB)'
    else
       write(99,*) ' ', i, ':', to_sp(r(i))
    end if
 end do
!
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Rigidity of the layers from bottom to top' 
!
!
 do i=1,nla
    if(i==1) then
       write(99,*) ' ', i, ':', to_sp(mu(i)), '(Core)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp(mu(i)), '(Lithosphere)'
    else
       write(99,*) ' ', i, ':', to_sp(mu(i))   
    end if
 end do
!
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Density of the layers from bottom to top' 
!
!
 do i=1,nla
    if(i==1) then
       write(99,*) ' ', i, ':', to_sp(rho(i)), '(Core)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp(rho(i)), '(Lithosphere)'
    else
       write(99,*) ' ', i, ':', to_sp(rho(i))    
    end if
 end do
!
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Viscosity of the layers from bottom to top'
!
 do i=1,nla
    if(i==1) then
       write(99,*) ' ', i, ':', to_dp(eta(i)), '(Core)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_dp(eta(i)), ' (Lithosphere)'
    else
       write(99,*) ' ', i, ':', to_dp(eta(i))    
    end if
 end do
!
!
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Mass of the Earth = ', to_dp(mass) 
 write(99,*) ' '
 write(99,*) ' NORMALIZED Newton constant = ', to_dp(G) 
!
 write(99,*) ' '
 write(99,*) ' NORMALIZED Gravity at the interfaces from bottom to top'
!
 do i=0,nla
    if(i==0) then
       write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i)), '(Planet center)'
    elseif(i==nla) then
       write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i)), '(External surface)'	 
    elseif(i==1) then
       write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i)), '(CMB)'	
    else
       write(99,*) ' ', i, ':', to_sp((r(i)*r0))/1000, to_sp(gra(i))
    end if
 end do
!
!
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
