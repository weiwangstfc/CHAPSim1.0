!> @par subroutine FACT
!> @param N [in]
!> @param IFAX [out]
!> @note What's the difference between FACT and FAX? (Same N-->Different IFAX)

      SUBROUTINE FACT(N,IFAX)
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      INCLUDE 'mpif.h'   
!     FACTORIZATION ROUTINE THAT FIRST EXTRACTS ALL FACTORS OF 4
!      DIMENSION IFAX(13)
      INTEGER(4),INTENT(INOUT)  :: IFAX(13)
      INTEGER(4),INTENT(IN)      :: N
      
      INTEGER(4)   :: NN
      INTEGER(4)   :: K
      INTEGER(4)   :: L
      INTEGER(4)   :: INC
      INTEGER(4)   :: MAXD
!***************************************************************      
      IF (N.GT.1) GO TO 10
      IFAX(1) = 0
      IF (N.LT.1) IFAX(1) = -99
      RETURN
      
   10 NN=N
      K=1
!     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
!     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
!     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
!     NOW FIND REMAINING FACTORS
   50 L=5
      MAXD = INT( DSQRT(DBLE(NN)) )
      INC=2
!     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 IF (L.GT.MAXD) GO TO 75
      L=L+INC
      INC=6-INC
      GO TO 60
   75 K = K+1
      IFAX(K) = NN
   80 IFAX(1)=K-1
!     IFAX(1) NOW CONTAINS NUMBER OF FACTORS
      RETURN
      END
