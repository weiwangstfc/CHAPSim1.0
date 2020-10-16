!> @par SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"
!>     PERFORMS ONE PASS THROUGH DATA
!>     AS PART OF MULTIPLE COMPLEX FFT ROUTINE
!>     FFT of two sets of real data (combined into a set of complex) in one pass
!>     Procedure 2 in paper Cooley, Lewist and Welch, 1970. P274.
!>
!>@note input data A+Bi ==> C+Di  DFT(A)=C; DFT(B)=D
 
!>@param A [in] IS FIRST REAL INPUT VECTOR
!>@param B [in] IS FIRST IMAGINARY INPUT VECTOR
!>               A+iB
!>@param C [out] IS FIRST REAL OUTPUT VECTOR
!>@param D [out] IS FIRST IMAGINARY OUTPUT VECTOR
!>                C+iD
!>@param TRIGS IS PRECALCULATED TABLE OF SINES " COSINES
!>@param INC1 IS ADDRESSING INCREMENT FOR A AND B
!>@param INC2 IS ADDRESSING INCREMENT FOR C AND D
!>@param INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
!>@param INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
!>@param LOT IS THE NUMBER OF VECTORS
!>@param N IS LENGTH OF VECTORS
!>@param IFAC IS CURRENT FACTOR OF N
!>@param LA IS PRODUCT OF PREVIOUS FACTORS


      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
      USE WPRECISIONFFT
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!    INCLUDE 'mpif.h'   
!      DIMENSION A(*),B(*),C(*),D(*),TRIGS(*)
!      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,     &
!           SIN72/0.951056516295154/,COS72/0.309016994374947/,      &
!           SIN60/0.866025403784437/
      !INTEGER,PARAMETER  :: WP=KIND(0.0D0) !WORKING PRECESION
      REAL(WP),INTENT(IN)     :: A(*)
      REAL(WP),INTENT(IN)     :: B(*)
      REAL(WP),INTENT(INOUT)  :: C(*)
      REAL(WP),INTENT(INOUT)  :: D(*)
      REAL(WP),INTENT(IN)      :: TRIGS(*)
      INTEGER(4),INTENT(IN)  :: INC1
      INTEGER(4),INTENT(IN)  :: INC2
      INTEGER(4),INTENT(IN)  :: INC3
      INTEGER(4),INTENT(IN)  :: INC4
      INTEGER(4),INTENT(IN)  :: LOT
      INTEGER(4),INTENT(IN)  :: N
      INTEGER(4),INTENT(IN)  :: IFAC
      INTEGER(4),INTENT(IN)  :: LA
      
      INTEGER(4)  :: M
      INTEGER(4)  :: I, IA, IB, IC, ID, IE, IBASE, IINK, IGO, IJK
      INTEGER(4)  :: J, JA, JB, JC, JD, JE, JBASE, JINK, JUMP
      INTEGER(4)  :: K, KB, KC, KD, KE  
      INTEGER(4)  :: L, LA1 
      REAL(WP)     :: C1, C2, C3, C4
      REAL(WP)     :: S1, S2, S3, S4
      REAL(WP)     :: SIN36
      REAL(WP)     :: COS36
      REAL(WP)     :: SIN72
      REAL(WP)     :: COS72
      REAL(WP)     :: SIN60
       
      
!************************************************************************

      SIN36=0.587785252292473_WP
      COS36=0.809016994374947_WP
      SIN72=0.951056516295154_WP
      COS72=0.309016994374947_WP
      SIN60=0.866025403784437_WP
      
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
!
!     CODING FOR FACTOR 2
!
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO L=1,LA
         I=IBASE
         J=JBASE
!DIR$ IVDEP
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            D(JA+J)=B(IA+I)+B(IB+I)
            C(JB+J)=A(IA+I)-A(IB+I)
            D(JB+J)=B(IA+I)-B(IB+I)
            I=I+INC3
            J=J+INC4
         END DO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      END DO
   
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         DO L=1,LA
            I=IBASE
            J=JBASE
!DIR$ IVDEP
            DO IJK=1,LOT
               C(JA+J)=A(IA+I)+A(IB+I)
               D(JA+J)=B(IA+I)+B(IB+I)
               C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
               D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
               I=I+INC3
               J=J+INC4
            END DO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          END DO
         JBASE=JBASE+JUMP
      END DO
   
      RETURN
!
!     CODING FOR FACTOR 3
!
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO IJK=1,LOT
             C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
             D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
             C(JB+J)=(A(IA+I)-0.50_WP*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)  &
                    -B(IC+I)))
             C(JC+J)=(A(IA+I)-0.50_WP*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)  &
                    -B(IC+I)))
             D(JB+J)=(B(IA+I)-0.50_WP*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)  &
                    -A(IC+I)))
             D(JC+J)=(B(IA+I)-0.50_WP*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)  &
                    -A(IC+I)))
             I=I+INC3
             J=J+INC4
          END DO
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
      END DO
   
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         KC=KB+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         DO L=1,LA
            I=IBASE
            J=JBASE
!DIR$ IVDEP
            DO IJK=1,LOT
               C(JA+J)= A(IA+I)+(A(IB+I)+A(IC+I))
               D(JA+J)= B(IA+I)+(B(IB+I)+B(IC+I))
               C(JB+J)= C1*((A(IA+I)-0.50_WP*(A(IB+I)+A(IC+I)))      &
                        -(SIN60*(B(IB+I)-B(IC+I))))                &
                        -S1*((B(IA+I)-0.50_WP*(B(IB+I)+B(IC+I)))     &
                        +(SIN60*(A(IB+I)-A(IC+I))))
               D(JB+J)= S1*((A(IA+I)-0.50_WP*(A(IB+I)+A(IC+I)))      &
                        -(SIN60*(B(IB+I)-B(IC+I))))                &
                        +C1*((B(IA+I)-0.50_WP*(B(IB+I)+B(IC+I)))     &
                        +(SIN60*(A(IB+I)-A(IC+I))))

               C(JC+J)= C2*((A(IA+I)-0.50_WP*(A(IB+I)+A(IC+I)))      &
                        +(SIN60*(B(IB+I)-B(IC+I))))                &
                        -S2*((B(IA+I)-0.50_WP*(B(IB+I)+B(IC+I)))     &
                        -(SIN60*(A(IB+I)-A(IC+I))))
     
               D(JC+J)= S2*((A(IA+I)-0.50_WP*(A(IB+I)+A(IC+I)))      &
                        +(SIN60*(B(IB+I)-B(IC+I))))                &
                        +C2*((B(IA+I)-0.50_WP*(B(IB+I)+B(IC+I)))     &
                        -(SIN60*(A(IB+I)-A(IC+I))))

               I=I+INC3
               J=J+INC4
               END DO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         END DO
         JBASE=JBASE+JUMP
      END DO
   
      RETURN
!
!     CODING FOR FACTOR 4
!
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      
      DO 100 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 95 IJK=1,LOT
             C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
             C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
             D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
             D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
             C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
             C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
             D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
             D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
             I=I+INC3
             J=J+INC4
   95    CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
  100 CONTINUE
  
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         KC=KB+KB
         KD=KC+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         DO L=1,LA
            I=IBASE
            J=JBASE
!DIR$ IVDEP
            DO IJK=1,LOT
               C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
               D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
               C(JC+J)= C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))  &
                      -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))     
               D(JC+J)= S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))  &
                      +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
               C(JB+J)= C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))  &
                      -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
               D(JB+J)= S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))  &
                      +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
               C(JD+J)= C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))  &
                      -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
               D(JD+J)= S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))  &
                      +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))    
               I=I+INC3
               J=J+INC4
            END DO
  
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         END DO
  
         JBASE=JBASE+JUMP
      END DO
  
      RETURN
!
!     CODING FOR FACTOR 5
!
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO L=1,LA
         I=IBASE
         J=JBASE
!DIR$ IVDEP
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
            D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
            C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))                   &
                   -COS36*(A(IC+I)+A(ID+I)))                           &
                   -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
            C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))                   &
                   -COS36*(A(IC+I)+A(ID+I)))                           &
                   +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
            D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))                   &
                   -COS36*(B(IC+I)+B(ID+I)))                           &
                   +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
            D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))                   &
                   -COS36*(B(IC+I)+B(ID+I)))                           &
                   -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
            C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))                   &
                   +COS72*(A(IC+I)+A(ID+I)))                           &
                   -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
            C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))                   &
                   +COS72*(A(IC+I)+A(ID+I)))                           &
                   +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
            D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))                   &
                   +COS72*(B(IC+I)+B(ID+I)))                           &
                   +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
            D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))                   &
                   +COS72*(B(IC+I)+B(ID+I)))                           &
                   -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
            I=I+INC3
            J=J+INC4
         END DO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      END DO
  
      IF (LA.EQ.M) RETURN
      
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO K=LA1,M,LA
         KB=K+K-2
         KC=KB+KB
         KD=KC+KB
         KE=KD+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         C4=TRIGS(KE+1)
         S4=TRIGS(KE+2)
         DO L=1,LA
            I=IBASE
            J=JBASE
!DIR$ IVDEP
            DO IJK=1,LOT
               C(JA+J)= A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
               D(JA+J)= B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
               C(JB+J)= C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))        &
                       -COS36*(A(IC+I)+A(ID+I)))                    &
                       -(SIN72*(B(IB+I)-B(IE+I))                    &
                       +SIN36*(B(IC+I)-B(ID+I))))                   &
                       -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))        &
                       -COS36*(B(IC+I)+B(ID+I)))                    &
                       +(SIN72*(A(IB+I)-A(IE+I))                    &
                       +SIN36*(A(IC+I)-A(ID+I))))
               D(JB+J)= S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))        &
                       -COS36*(A(IC+I)+A(ID+I)))                    &
                       -(SIN72*(B(IB+I)-B(IE+I))                    &
                       +SIN36*(B(IC+I)-B(ID+I))))                   &
                       +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))        &
                       -COS36*(B(IC+I)+B(ID+I)))                    &
                       +(SIN72*(A(IB+I)-A(IE+I))                    &
                       +SIN36*(A(IC+I)-A(ID+I))))
               C(JE+J)= C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))        & 
                       -COS36*(A(IC+I)+A(ID+I)))                    &
                       +(SIN72*(B(IB+I)-B(IE+I))                    &
                       +SIN36*(B(IC+I)-B(ID+I))))                   &
                       -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))        &
                       -COS36*(B(IC+I)+B(ID+I)))                    &
                       -(SIN72*(A(IB+I)-A(IE+I))                    &
                       +SIN36*(A(IC+I)-A(ID+I))))
               D(JE+J)= S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))        &
                       -COS36*(A(IC+I)+A(ID+I)))                    &
                       +(SIN72*(B(IB+I)-B(IE+I))                     &
                       +SIN36*(B(IC+I)-B(ID+I))))                    &
                       +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))         &
                       -COS36*(B(IC+I)+B(ID+I)))                     &
                       -(SIN72*(A(IB+I)-A(IE+I))                     &
                       +SIN36*(A(IC+I)-A(ID+I))))
               C(JC+J)= C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))        &
                       +COS72*(A(IC+I)+A(ID+I)))                    &
                       -(SIN36*(B(IB+I)-B(IE+I))                    &
                       -SIN72*(B(IC+I)-B(ID+I))))                   &
                       -S2*((B(IA+I)-COS36*(B(IB+I)                 &
                       +B(IE+I))+COS72*(B(IC+I)+B(ID+I)))           &
                       +(SIN36*(A(IB+I)-A(IE+I))                    &
                       -SIN72*(A(IC+I)-A(ID+I))))
               D(JC+J)= S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))        &
                       +COS72*(A(IC+I)+A(ID+I)))                    &
                       -(SIN36*(B(IB+I)-B(IE+I))                    &
                       -SIN72*(B(IC+I)-B(ID+I))))                   &
                       +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))        &
                       +COS72*(B(IC+I)+B(ID+I)))                    &
                       +(SIN36*(A(IB+I)-A(IE+I))                    &
                       -SIN72*(A(IC+I)-A(ID+I))))
               C(JD+J)= C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))        &
                       +COS72*(A(IC+I)+A(ID+I)))                    &
                       +(SIN36*(B(IB+I)-B(IE+I))                    &
                       -SIN72*(B(IC+I)-B(ID+I))))                   &
                       -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))        &
                       +COS72*(B(IC+I)+B(ID+I)))                    &
                       -(SIN36*(A(IB+I)-A(IE+I))                    &
                       -SIN72*(A(IC+I)-A(ID+I))))
               D(JD+J)= S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))        &
                       +COS72*(A(IC+I)+A(ID+I)))                    &
                       +(SIN36*(B(IB+I)-B(IE+I))                    &
                       -SIN72*(B(IC+I)-B(ID+I))))                   &
                       +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))        &
                       +COS72*(B(IC+I)+B(ID+I)))                    &
                       -(SIN36*(A(IB+I)-A(IE+I))                    &
                       -SIN72*(A(IC+I)-A(ID+I))))
               I=I+INC3
               J=J+INC4
            END DO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         END DO
         JBASE=JBASE+JUMP
      END DO
  
      RETURN
      END
