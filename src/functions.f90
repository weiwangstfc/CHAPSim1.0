!**********************************************************************************************************************************
!> @brief
!>        some functions
!> @details
!> FUNCTION: DGAL
!> FUNCTION: GAL
!> FUNCTION: ISNAN1
!> FUNCTION: ISINF1
!> FUNCTION: logbaseb
!> @note
!> @toDO
! REVISION HISTORY:
! 06/ 2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
FUNCTION DGAL(AL,ZITAF)
    USE WPRECISION
    IMPLICIT NONE
    !      IMPLICIT DOUBLE PRECISION(A- H,O-Z)
    REAL(WP), INTENT(IN) :: AL
    REAL(WP), INTENT(IN) :: ZITAF
    REAL(WP) :: DGAL

    REAL(WP) :: DG1
    REAL(WP) :: DG2
    REAL(WP) :: DG3
    REAL(WP) :: ALZITF


    ALZITF = AL * ZITAF
    DG1 = 1.0_WP / (DSINH(ALZITF) * DCOSH(ALZITF))
    DG2 = - ALZITF / DSINH(ALZITF)**2
    DG3 = - ALZITF / DCOSH(ALZITF)**2
    DGAL = DG1 + DG2 + DG3

    RETURN
END

!**********************************************************************************************************************************
FUNCTION GAL(AL,ZITAF)
    USE WPRECISION
    IMPLICIT NONE
    !      IMPLICIT DOUBLE PRECISION (A- H,O-Z)
    REAL(WP), INTENT(IN) :: AL
    REAL(WP), INTENT(IN) :: ZITAF
    REAL(WP) :: GAL
    REAL(WP) :: ALZITF

    ALZITF = AL * ZITAF
    GAL = AL / (DCOSH(ALZITF) * DSINH(ALZITF))

    RETURN
END

!**********************************************************************************************************************************
LOGICAL FUNCTION ISNAN1(a)
    USE WPRECISION
    IMPLICIT NONE

    REAL(WP) :: a
    !IF (a /= a) THEN  !commented by WW
    !IF ((a + 1.0) == a) THEN ! for pgi pgf90
    !ISNAN1 = .TRUE.
    !ELSE IF (a/= A) THEN
    !ISNAN1 = .TRUE.   !for intel IFort
    !ELSE
    !ISNAN1 = .FALSE.
    !endIF

    IF( ((a + 1.0_WP) == a) .OR. (a /= a) ) THEN
        ISNAN1 = .TRUE.
    ELSE
        ISNAN1 = .FALSE.
    END IF


    RETURN
end FUNCTION ISNAN1

!**********************************************************************************************************************************
LOGICAL FUNCTION ISINF1(a)
    USE WPRECISION
    IMPLICIT NONE

    REAL(WP) :: a
    REAL(WP) :: b
    !IF (a /= a) THEN  !commented by WW
    !IF ((a + 1.0) == a) THEN ! for pgi pgf90
    !ISNAN1 = .TRUE.
    !ELSE IF (a/= A) THEN
    !ISNAN1 = .TRUE.   !for intel IFort
    !ELSE
    !ISNAN1 = .FALSE.
    !endIF
    b = 0.0_WP
    b = HUGE(b)

    IF( DABS(a) >   DABS(b) ) THEN
        ISINF1 = .TRUE.
    ELSE
        ISINF1 = .FALSE.
    END IF


    RETURN
end FUNCTION ISINF1
!**********************************************************************************************************************************

REAL(WP) FUNCTION logbaseb(a, b)
    USE WPRECISION
    IMPLICIT NONE

    REAL(WP) :: a
    REAL(WP) :: b

    logbaseb = log10(a) / log10(b)

    RETURN
END FUNCTION logbaseb
