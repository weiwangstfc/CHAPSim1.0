!> @par CFTRIG :  Continuous Fourier Transform (CFT)
!> @param N     [in]
!> @param TRIGS [out] TRIGS(i)=cos(pi*(i-1)/N); TRIGS(i+1)=sin(pi*(i-1)/N)
      
    SUBROUTINE CFTRIG(N,TRIGS)
        USE WPRECISIONFFT
        IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      INCLUDE 'mpif.h'   
!      DIMENSION TRIGS(*)

        REAL(WP),INTENT(INOUT)  :: TRIGS(*) 
        INTEGER(4),INTENT(IN)  :: N
      
        REAL(WP)    :: PI
        REAL(WP)    :: DEL
        REAL(WP)    :: ANGLE
        INTEGER(4) :: I, L
      
      
!**********************************************************************      
      PI=2.0_WP*DASIN(1.0_WP)
      DEL=(PI+PI)/DBLE(N)  !DEL=2Pi/N
      L=N+N
      
      DO I=1,L,2
         ANGLE=0.50_WP*DBLE(I-1)*DEL  !ANGLE=Pi(I-1)/N
         TRIGS(I)=DCOS(ANGLE)       !
         TRIGS(I+1)=DSIN(ANGLE)
      END DO
   
      RETURN
      END
