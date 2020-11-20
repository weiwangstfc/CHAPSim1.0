!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> SUBROUTINE: random_initialize
!> SUBROUTINE: rvec_random
!> SUBROUTINE: r_random
!> @note
!> @todo
! REVISION HISTORY:
! 10/2013 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE random_initialize ( seed )
    USE WPRECISION
    !
    !*******************************************************************************
    !
    !! random_INITIALIZE initializes the FORTRAN 90 random number seed.
    !
    !
    !  DIScussion:
    !
    !    IF you DOn't initialize the random number generator, its behavior
    !    IS not specified.  IF you initialize it simply by:
    !
    !      CALL random_seed
    !
    !    its behavior IS not specified.  On the DEC ALPHA, IF that's all you
    !    DO, the same random number sequence IS RETURNed.  In order to actually
    !    try to scramble up the random number generator a bit, thIS routine
    !    goes through the tedious process of getting the size of the random
    !    number seed, making up values based on the current time, and setting
    !    the random number seed.
    !
    !  ModIFied:
    !
    !    19 December 2001
    !
    !  Author:
    !
    !    John BurkARdt
    !
    !  parameters:
    !
    !    Input/output, INTEGER(4) SEED.
    !    IF SEED IS zero on input, THEN you're asking thIS routine to come up
    !    with a seed value, whICh IS RETURNed as output.
    !    IF SEED IS nonzero on input, THEN you're asking thIS routine to
    !    USE the input value of SEED to initialize the random number generator,
    !    and SEED IS not changed on output.
    !
    IMPLICIT none
    !
    INTEGER(4) :: count
    INTEGER(4) :: count_max
    INTEGER(4) :: count_rate
    LOGICAL, PARAMETER :: debug = .FALSE.
    INTEGER(4) :: i
    INTEGER(4) :: seed
    INTEGER(4), ALLOCATABLE :: seed_vector(:)
    INTEGER(4) :: seed_size
    REAL(WP) :: t
    !
    !  Initialize the random number seed.
    !
    CALL random_seed
    !
    !  Determine the size of the random number seed.
    !
    CALL random_seed ( size = seed_size )
    !
    !  Allocate a seed of the right size.
    !
    allocate ( seed_vector(seed_size) ); seed_vector = 0

    IF ( seed /= 0 ) THEN

        IF ( debug ) THEN
            WRITE ( *, '(a)' ) ' '
            WRITE ( *, '(a)' ) 'random_INITIALIZE'
            WRITE ( *, '(a, I20)' ) '  Initialize random_NUMBER, USEr SEED = ', seed
        END IF

    ELSE

        CALL system_Clock ( count, count_rate, count_max )

        seed = count

        IF ( debug ) THEN
            WRITE ( *, '(a)' ) ' '
            WRITE ( *, '(a)' ) 'random_INITIALIZE'
            WRITE ( *, '(a, I20)' ) '  Initialize random_NUMBER, ARbitrARy SEED = ', &
            seed
        END IF

    END IF
    !
    !  Now set the seed.
    !
    seed_vector(1:seed_size) = seed

    CALL random_seed ( put = seed_vector(1:seed_size) )
    !
    !  Free up the seed space.
    !
    DEALLOCATE ( seed_vector )
    !
    !  CALL the random number routine a bunch of times.
    !random_initialize
    DO I = 1, 100
        CALL random_number ( harvest = t )
    END DO

    RETURN
end

!**********************************************************************************************************************************
SUBROUTINE rvec_random ( alo, ahi, n, a )
    USE WPRECISION
    !
    !*******************************************************************************
    !
    !! RVEC_random RETURNs a random REAL(WP) vector in a given range.
    !
    !
    !  ModIFied:
    !
    !    04 FebruARy 2001
    !
    !  Author:
    !
    !    John BurkARdt
    !
    !  parameters:
    !
    !    Input, REAL(WP) ALO, AHI, the range allowed for the entries.
    !
    !    Input, INTEGER(4) N, the number of entries in the vector.
    !
    !    Output, REAL(WP) A(N), the vector of randomly chosen values.
    !
    IMPLICIT none
    !
    INTEGER(4) n
    !
    REAL(WP) a(N)
    REAL(WP) ahi
    REAL(WP) alo
    INTEGER(4) i
    !
    DO I = 1, n
        CALL r_random ( alo, ahi, a(I) )
    END DO

    RETURN
end

!**********************************************************************************************************************************
SUBROUTINE r_random ( rlo, rhi, r )
    USE WPRECISION
    !
    !*******************************************************************************
    !
    !! R_random RETURNs a random REAL(WP) in a given range.
    !
    !
    !  ModIFied:
    !
    !    06 April 2001
    !
    !  Author:
    !
    !    John BurkARdt
    !
    !  parameters:
    !
    !    Input, REAL(WP) RLO, RHI, the minimum and maximum values.
    !
    !    Output, REAL(WP) R, the randomly chosen value.
    !
    IMPLICIT none
    !
    REAL(WP) :: r
    REAL(WP) :: rhi
    REAL(WP) :: rlo
    REAL(WP) :: t
    !
    !  PICk T, a random number in (0, 1).
    !
    CALL random_number ( harvest = t )
    !
    !  Set R in ( RLO, RHI ).
    !
    r = ( 1.0E+00 - t ) * rlo + t * rhi

    RETURN
end
