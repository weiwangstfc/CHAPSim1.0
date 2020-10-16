!> @par  SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
!>    (GRIDPOINT TO SPECTRAL TRANSFORM)

      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
      USE WPRECISIONFFT
      IMPLICIT NONE
!     INCLUDE 'mpif.h'   
!      DIMENSION TRIGS(*)
!      DOUBLE PRECISION WORK(*),A(*)
      !INTEGER,PARAMETER  :: WP=KIND(0.0D0) !WORKING PRECESION
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
      INTEGER(4) :: K, L
      REAL(WP)     :: C
      REAL(WP)     :: S
      REAL(WP)     :: SCALE0
       
!******************************************************
      NH=N/2
      NX21=N+1
      INK=INC+INC

!     A(0) AND A(N/2)
      SCALE0=1.0_wp/DBLE(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
!DIR$ IVDEP
      DO 10 L=1,LOT
         A(JA)=SCALE0*(WORK(IA)+WORK(IB))
         A(JB)=SCALE0*(WORK(IA)-WORK(IB))
         A(JA+INC)=0.00_wp
         A(JB+INC)=0.00_wp
         IA=IA+NX21
         IB=IB+NX21
         JA=JA+JUMP
         JB=JB+JUMP
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      SCALE0=0.50_wp*SCALE0
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
!
      DO 30 K=3,NH,2
         IA=IABASE
         IB=IBBASE
         JA=JABASE
         JB=JBBASE
         C=TRIGS(N+K)
         S=TRIGS(N+K+1)
!DIR$ IVDEP
         DO 20 L=1,LOT
            A(JA)=SCALE0*((WORK(IA)+WORK(IB))                         &
                  +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
            A(JB)=SCALE0*((WORK(IA)+WORK(IB))                         &
                  -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
            A(JA+INC)=SCALE0*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)    &
                      +WORK(IB+1)))+(WORK(IB+1)-WORK(IA+1)))
            A(JB+INC)=SCALE0*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)    &
                      +WORK(IB+1)))-(WORK(IB+1)-WORK(IA+1)))
            IA=IA+NX21
            IB=IB+NX21
            JA=JA+JUMP
            JB=JB+JUMP
   20    CONTINUE
         IABASE=IABASE+2
         IBBASE=IBBASE-2
         JABASE=JABASE+INK
         JBBASE=JBBASE-INK
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE0=2.0_wp*SCALE0
!DIR$ IVDEP
      DO 40 L=1,LOT
         A(JA)=SCALE0*WORK(IA)
         A(JA+INC)=-SCALE0*WORK(IA+1)
         IA=IA+NX21
         JA=JA+JUMP
   40 CONTINUE
!
   50 CONTINUE
      RETURN
      END
