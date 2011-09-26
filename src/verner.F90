
! photoionization cross-sections from Verner
!===========================================================

  !> HI photo ionization x-section (Verner) [cm^2]
  !-------------------------------------------------------------------------
  function Verner_HI_photo_cs( Ry ) result( sigma )
    implicit none
    real, intent(in) :: Ry  !< energy [Rydbergs]
    real :: sigma           !< cross section [cm^2]
    real, parameter :: Eth = 13.6d0
    real, parameter :: Emax = 5.0d4
    real, parameter :: E0 = 4.298d-1
    real, parameter :: sig0 = 5.475d4
    real, parameter :: ya = 3.288d1
    real, parameter :: P = 2.963d0

    real :: eV
    real :: x
    real :: y

    if (Ry < 1.) then
       sigma = 0.0
       return
    endif
    eV = Ry * 13.6d0
    x = eV / E0
    y = x
  
    sigma = sig0 * (x-1)**2 * y**(0.5d0 * P - 5.5d0) * (1 + sqrt(y/ya))**(-P)
    sigma = sigma * 1.0d-18
    if (sigma < 0) sigma = 0.0
  end function Verner_HI_photo_cs

!------------------------------------------
!> HeI photo ionization (Osterbrok)  [cm^2]
    function Osterbrok_HeI_photo_cs(Ry) result(sigma)    
        implicit none
        real, intent(in) :: Ry !< frequency [Ry]
        real :: sigma !< cross section
        real :: scaled_freq
        real, parameter :: nu_HeI = 24.587
        if (Ry * 13.6 < nu_HeI) then
           sigma = 0.0e0
           return
        end if
        scaled_freq = Ry * 13.6 / nu_HeI
        sigma = 7.2e-18 * & 
                ( 1.66e0 * ( scaled_freq )**(-2.05e0) + &
                  0.66e0 * ( scaled_freq )**(-3.05e0) )
    end function Osterbrok_HeI_photo_cs

!------------------------------------------
!> HeII photo ionization (Osterbrok)  [cm^2]
    function Osterbrok_HeII_photo_cs(Ry) result(sigma)    
        implicit none
        real, intent(in) :: Ry !< frequency (Ry)
        real :: sigma !< cross section
        real, parameter :: nu_HeII = 54.416
        if (Ry * 13.6 < nu_HeII) then
           sigma = 0.0e0
           return
        end if
        sigma = 1.58e-18 * ( Ry * 13.6 / nu_HeII ) ** (-3)
    end function Osterbrok_HeII_photo_cs
