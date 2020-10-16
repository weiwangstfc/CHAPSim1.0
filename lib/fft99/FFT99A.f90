!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!> @par SUBROUTINE FFT99A - 
!>      PREPROCESSING STEP FOR FFT99, ISIGN=+1 (SPECTRAL TO GRIDPOINT TRANSFORM)

      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      USE WPRECISIONFFT
      IMPLICIT NONE
!C     INCLUDE 'mpif.h'   
!      DOUBLE PRECISION A(*),WORK(*)
!      DIMENSION TRIGS(*)
      REAL(WP),INTENT(INOUT)  :: A(*)
      REAL(WP),INTENT(INOUT)  :: WORK(*)
      REAL(WP),INTENT(IN)      :: TRIGS(*)
      INTEGER(4),INTENT(IN)  :: INC
      INTEGER(4),INTENT(IN)  :: JUMP
      INTEGER(4),INTENT(IN)  :: N
      INTEGER(4),INTENT(IN)  :: LOT
 
      INTEGER(4) :: NH, NX21
      INTEGER(4) :: INK
      INTEGER(4) :: IA, IB, IABASE,IBBASE
      INTEGER(4) :: JA, JB, JABASE, JBBASE
      INTEGER(4) :: L, K
      REAL(WP)     :: C
      REAL(WP)     :: S 
!*****************************************************************
      NH=N/2
      NX21=N+1
      INK=INC+INC
      
!>    @note transfer data A to work space, the ends of work spaces. 
!C     A(0) AND A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
!DIR$ IVDEP
      DO 10 L=1,LOT
         WORK(JA)=A(IA)+A(IB)
         WORK(JB)=A(IA)-A(IB)
         IA=IA+JUMP
         IB=IB+JUMP
         JA=JA+NX21
         JB=JB+NX21
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
!
      DO 30 K=3,NH,2
         IA=IABASE
         IB=IBBASE
         JA=JABASE
         JB=JBBASE
         C=TRIGS(N+K)    !cos
         S=TRIGS(N+K+1)  !sin
!DIR$ IVDEP
         DO 20 L=1,LOT
            WORK(JA)=(A(IA)+A(IB))-                                &
                     (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
            WORK(JB)=(A(IA)+A(IB))+                                &
                     (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
            WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+  &
                       (A(IA+INC)-A(IB+INC))
            WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-  &
                       (A(IA+INC)-A(IB+INC))
            IA=IA+JUMP
            IB=IB+JUMP
            JA=JA+NX21
            JB=JB+NX21
   20    CONTINUE
         IABASE=IABASE+INK
         IBBASE=IBBASE-INK
         JABASE=JABASE+2
         JBBASE=JBBASE-2
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
!DIR$ IVDEP
      DO 40 L=1,LOT
          WORK(JA)=2.0_wp*A(IA)
          WORK(JA+1)=-2.0_wp*A(IA+INC)
          IA=IA+JUMP
          JA=JA+NX21
   40 CONTINUE
!
   50 CONTINUE
      RETURN
      END
