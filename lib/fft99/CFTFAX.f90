!> @par  Subroutine CFTFAX : Continuous Fourier Transform (CFT)
!>     THIS ROUTINE WAS MODIFIED FROM TEMPERTON"S ORIGINAL
!>     BY DAVE FULKER.  IT NO LONGER PRODUCES FACTORS IN ASCENDING
!>     ORDER, AND THERE ARE NONE OF THE ORIGINAL 'MODE' OPTIONS.
!
!> @param N [in]
!>        THE LENGTH OF EACH COMPLEX TRANSFORM TO BE PERFORMED
!>        N MUST BE GREATER THAN 1 AND CONTAIN NO PRIME
!>        FACTORS GREATER THAN 5.
!
!> @param IFAX[out]
!               IFAX(1)
!                 THE NUMBER OF FACTORS CHOSEN OR -99 IN CASE OF ERROR
!               IFAX(2) THRU IFAX( IFAX(1)+1 )
!                 THE FACTORS OF N IN THE FOLLOWIN ORDER:  APPEARING
!                 FIRST ARE AS MANY FACTORS OF 4 AS CAN BE OBTAINED.
!                 SUBSEQUENT FACTORS ARE PRIMES, AND APPEAR IN
!                 ASCENDING ORDER, EXCEPT FOR MULTIPLE FACTORS.
!
!> @param TRIGS[out]
!               2N SIN AND COS VALUES FOR USE BY THE TRANSFORM ROUTINE
!

      SUBROUTINE CFTFAX(N,IFAX,TRIGS)
      USE WPRECISIONFFT
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      INCLUDE 'mpif.h'   
!      DIMENSION IFAX(13),TRIGS(*)

      REAL(WP),INTENT(INOUT)      :: TRIGS(*)
      INTEGER(4),INTENT(INOUT)  :: IFAX(13) 
      INTEGER(4),INTENT(IN)      :: N
      
      INTEGER(4) :: K

!************************************************************
      CALL FACT(N,IFAX)
      K = IFAX(1)
      IF (K .LT. 1 .OR. IFAX(K+1) .GT. 5) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 ) WRITE(6,*) '# CFTFAX -- INVALID N'
      
      CALL CFTRIG (N, TRIGS)
      
      RETURN
      END
      

      

