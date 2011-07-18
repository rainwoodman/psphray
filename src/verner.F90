! photoionization cross-sections from Verner
!===========================================================

  !> HI photo ionization x-section (Verner) [cm^2]
  !-------------------------------------------------------------------------
  function Verner_HI_photo_cs( Ry ) result( sigma )
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

    eV = Ry * 13.6d0
    x = eV / E0
    y = x
  
    sigma = sig0 * (x-1)**2 * y**(0.5d0 * P - 5.5d0) * (1 + sqrt(y/ya))**(-P)
    sigma = sigma * 1.0d-18

  end function Verner_HI_photo_cs
