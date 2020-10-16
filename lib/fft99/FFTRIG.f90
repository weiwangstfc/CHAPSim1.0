!> @par Subroutine FFTRIG
!>      trigonometric variables/constants
!> @param N     [in]
!> @param MODE  [in]  =3
!> @param TRIGS [out] AN ARRAY OF TRIGNOMENTRIC FUNCTION VALUES SUBSEQUENTLY
!>                    USED BY THE FFT ROUTINES.   
!>                    A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
!>                    EVEN, OR 3*N/2+1 IF N/2 IS ODD.

      SUBROUTINE FFTRIG(TRIGS,N,MODE)
      USE WPRECISIONFFT
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     INCLUDE 'mpif.h'   
!      DIMENSION TRIGS(*)
        !INTEGER,PARAMETER  :: WP=KIND(0.0D0) !WORKING PRECESION
        REAL(WP),INTENT(INOUT) :: TRIGS(*)
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: MODE 
      
      REAL(WP)    :: PI
      REAL(WP)    :: DEL
      REAL(WP)    :: ANGLE
      INTEGER(4) :: IMODE 
      INTEGER(4) :: NN, NH
      INTEGER(4) :: L, LA
      INTEGER(4) :: I
      
!******************************************************************      
      PI=2.0_WP*DASIN(1.0_WP)
      IMODE=IABS(MODE)
      NN=N
      
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/DBLE(NN)
      L=NN+NN
 
!>    @note I=1:(NN*2) . I=1,2*NN-1,2 for cos. I=2,2*NN,2 for sin.     
      DO 10 I=1,L,2
         ANGLE=0.50_WP*DBLE(I-1)*DEL
         TRIGS(I)=DCOS(ANGLE)
         TRIGS(I+1)=DSIN(ANGLE)
   10 CONTINUE
   
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
    
!>    @note I=(2*NN+1):( 2*NN+2*(NN+1)/2 )   
!>          if NN is even, I=(2*NN+1):(3*NN)
!>          if NN is odd,  I=(2*NN+1):(3*NN+1)   
      DEL=0.50_WP*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
         ANGLE=0.50_WP*DBLE(I-1)*DEL
         TRIGS(LA+I)=DCOS(ANGLE)
         TRIGS(LA+I+1)=DSIN(ANGLE)
   20 CONTINUE
   
!>    @note FOR IMODE==3, no use.   
      IF (IMODE.LE.3) RETURN
      
      DEL=0.50_WP*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      
      DO 30 I=2,NN
         ANGLE=DBLE(I-1)*DEL
         TRIGS(LA+I)=2.0_WP*DSIN(ANGLE)
   30 CONTINUE
   
      RETURN
      
   40 CONTINUE
      DEL=0.50_WP*DEL
      DO 50 I=2,N
         ANGLE=DBLE(I-1)*DEL
         TRIGS(LA+I)=DSIN(ANGLE)
   50 CONTINUE
   
      RETURN
      
      END SUBROUTINE FFTRIG
