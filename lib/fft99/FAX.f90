!> @par Subroutine FAX
!> @param IFAX [out]   \n 
!>               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
!>               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
!>               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.
!>        IFAX(1) = how many prime factors in N.
!>        IFAX(2:IFAX(1)+1) prime factors  based on 3,4,5 
!>        N=2*IFAX(2)*IFAX(3)*IFAX(4)...IFAX(IFAX(1)+1)
!>
!> @param N    [in]   cell no. in one direction
!> @param MODE [in]      

      SUBROUTINE FAX(IFAX,N,MODE)
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     INCLUDE 'mpif.h'   
!      DIMENSION IFAX(13)
      INTEGER(4),INTENT(INOUT)  :: IFAX(13)
      INTEGER(4),INTENT(IN)      :: N
      INTEGER(4),INTENT(IN)      :: MODE
      
      INTEGER(4)   :: NN
      INTEGER(4)   :: I, II, ITEM, ISTOP
      INTEGER(4)   :: K, L
      INTEGER(4)   :: INC
      INTEGER(4)   :: NFAX
!*******************************************************************      
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
 
      
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
!>    @note if cell no. N = odd number, below, IFAX(1)=-99 
      IFAX(1)=-99
      
      RETURN
      
   10 K=1
!>    @note  TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
!>    @note real(N/2/4) == integer, below IFAX(2)=4, K=2    
      K=K+1
      IFAX(K)=4
      NN=NN/4      
      IF (NN.EQ.1) GO TO 80      
      GO TO 20   
      
!>    @note TEST FOR EXTRA FACTOR OF 2
!>    @note N/2/2 has remaider. IFAX(2)=2, K=2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
      
!>    @note TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50   
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
      
!>    @note NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
!     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L      
      IF (NN.EQ.1) GO TO 80      
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
      
!>    @note N/2/4==1       
   80 IFAX(1)=K-1
!     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
!     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      
      DO 100 II=2,NFAX
         ISTOP=NFAX+2-II
         DO 90 I=2,ISTOP
            IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
            ITEM=IFAX(I)
            IFAX(I)=IFAX(I+1)
            IFAX(I+1)=ITEM
   90    CONTINUE
  100 CONTINUE
  
  110 CONTINUE
  
      RETURN
      END
      
      
