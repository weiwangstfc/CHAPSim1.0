!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!> @par    Introduction of FFT
!> PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
!>              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
!>              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
!>              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE
!>              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER
!>              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
!>              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS
!>              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.
!>              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE
!>              THAT IS MOSTLY WRITTEN IN CAL.
!>
!>              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES:
!>
!>            SUBROUTINE FFTFAX
!>                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
!>                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
!>                (PROVIDED THAT N IS NOT CHANGED).
!>
!>            SUBROUTINES FFT99 AND FFT991
!>                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT
!>                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE.
!>
!>@par
!> ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH
!>
!>                   *FORTRAN,P=XLIB,SN=FFT99F
!>              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY
!>              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE
!>              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE
!>              ACCESSED USING
!
!>                   FETCH P=CRAYLIB,SN=FFT99
!>                   FETCH P=CRAYLIB,SN=CAL99
!>@par
!> USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1,
!>              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
!>              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH
!>              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS
!>              OF LENGTH N IS
!
!>                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),
!>                  +          WORK(M*(N+1))
!
!>                   CALL FFTFAX (N, IFAX, TRIGS)
!>                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
!
!>              SEE THE INDIVIDUAL WRITE-UPS FOR FFTFAX, FFT99, AND
!>              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE
!>              ARGUMENTS.
!>@par
!> HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
!>              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
!>              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!> @par
!>     SUBROUTINE "FFT99" - MULTIPLE FAST REAL PERIODIC TRANSFORM
!>     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
!>     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!>     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!>     (1970), 315-337)
!
!>@param   A 
!>         IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!>@param   WORK
!>         IS AN AREA OF SIZE (N+1)*LOT
!>@param   TRIGS
!>         IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!>@param   IFAX
!>         IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!>@param   INC
!>         IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!>         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!>@param   JUMP
!>         IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!>@param   N 
!>         IS THE LENGTH OF THE DATA VECTORS
!>@param   LOT
!>         IS THE NUMBER OF DATA VECTORS
!>@param   ISIGN
!>     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!>           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!>     ORDERING OF COEFFICIENTS:
!>         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!>         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!>     ORDERING OF DATA:
!>         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
!>         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
!
!>     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!>     PARALLEL
!
!>     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!>     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!>     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!>         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!>     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*DCOS(2*J*K*PI/N))
!>               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*DSIN(2*J*K*PI/N))
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!> @par SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
!>                       AND
!> SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
!>@par
!> PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
!>              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
!>              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT
!>              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE
!>              GRIDPOINT VALUES (FFT99).  GIVEN A SET
!>              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF
!>              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE
!>              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN
!>              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS
!>              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE
!>              OPTIMIZED FOR USE ON THE CRAY-1.
!>@par
!> ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)
!> DIMENSIONS
!
!> ARGUMENTS
!>@par
!>
!> ON INPUT    
!>@param   A
!>               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA
!>               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY
!>               THE RESULTS.
!>@param  WORK             
!>               A WORK ARRAY OF DIMENSION M*(N+1)
!>@param  TRIGS             
!>               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
!>@param  IFAX            
!>               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
!>@param  INC              
!>               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF
!>               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR
!>               CONSECUTIVELY STORED DATA).
!>@param JUMP            
!>               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF
!>               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,
!>               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
!>               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
!>               INC AND JUMP, SEE THE EXAMPLES BELOW.
!>@param   N
!>               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
!>               TRANSFORMS, BELOW).
!>@param   M
!>               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
!>@param   ISIGN
!>               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
!>                    GRIDPOINT VALUES.
!>               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
!                    COEFFICIENTS.
!>@par
!> ON OUTPUT   
!> @param   A
!>               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED
!>               EACH CONTAINING THE SEQUENCE:
!
!>               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)
!
!>               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH
!>               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:
!
!>               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0.
!>               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).
!                   (EXPLICIT CYCLIC CONTINUITY)
!
!>               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY:
!>                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!>                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!>                 AND I=SQRT (-1)
!
!>               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH
!>               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS
!>               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS
!>               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS
!>               A(K), B(K), 0 .LE. K .LE N/2.
!
!>               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY:
!>                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N))
!>                 WHERE C(K)=A(K)+I*B(K) AND I=DSQRT(-1)
!
!>               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1
!>               (OR VICE VERSA) RETURNS THE ORIGINAL DATA.
!
!>               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL
!>               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1,
!>               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS.
!>@par
!> EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT
!>               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF
!>               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE
!>               ARRANGED LIKE THIS:
!>@par
!> FIRST DATA   A(1)=    . . .                A(66)=             A(70)
!> VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
!>@par
!> SECOND DATA  A(71)=   . . .                                  A(140)
!> VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
!
!>               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1,
!>               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC
!>               CONTINUITY).
!
!>               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS:
!
!>                FIRST         SECOND                          LAST
!>                DATA          DATA                            DATA
!>                VECTOR        VECTOR                          VECTOR
!
!>                 A(1)=         A(2)=                           A(19)=
!
!>                 X(63)         X(63)       . . .               X(63)
!>        A(20)=   X(0)          X(0)        . . .               X(0)
!>        A(39)=   X(1)          X(1)        . . .               X(1)
!                  .             .                               .
!                  .             .                               .
!                  .             .                               .
!
!>               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING
!>               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH
!>               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT
!>               DATA VECTOR.
!>
!-----------------------------------------------------------------------
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!>@note for PHACALC.f 
!>       FFT99(XR,WORK,TRIGXX1,IFXX1,1,M1+1,N1M,N3M,-1)
!>@param A     [in/out] =XR  [in]-array of QCAP with j fixed one cell forward (divergence of velocity divided by time step.)
!>                      XR(M1M+2,M3M)  or size of M*(N+2)
!>@param WORK  [out] size of M*(N+1) 
!>@param TRIGS [in] =TRIGXX1
!>@param IFAX  [in] =IFXX1
!>@param INC   [in] =1      intervals between data index based on the same spacing. 
!>@param JUMP  [in] =M1+1=N1M+2+1 
!>@param N     [in] =N1M
!>@param LOT   [in] =N3M
!>@param ISIGN [in] =-1

!>@warning     Insert a "cdir$ ivdep" or "!dir$ ivdep" statement right before the do-loop vectorize the loop. 
!>             [VERIFY] Make sure that these arrays in the loop do not have unsafe cross-iteration dependencies: foo. 
!>             A cross-iteration dependency exists if a memory location is modified in an iteration of the loop and 
!>             accessed (by a read or a write) in another iteration of the loop. Make sure that there are no such dependencies, 
!>             or that any cross-iteration dependencies can be safely ignored.   
!>             A test with gfortran gives an error message "Error: Rank mismatch in array reference"

      SUBROUTINE FFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      USE WPRECISIONFFT
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     INCLUDE 'mpif.h'   
!      INTEGER MYID
!     COMMON /NODES/SIZE,NUMPROCS,NPN,MYID,SW11,EW11  !"NPN" HAMAN "SIZE-1" AST!  !COMMON/NODES/MY_NODE,NUMPROCS,NPN
!

!      DIMENSION TRIGS(*),IFAX(13)
!      DOUBLE PRECISION A(*),WORK(*)
      REAL(WP),INTENT(INOUT)  :: A(*)
      REAL(WP),INTENT(INOUT)  :: WORK(*)
      REAL(WP),INTENT(IN)      :: TRIGS(*)
      INTEGER(4),INTENT(IN)  :: IFAX(13)
      INTEGER(4),INTENT(IN)  :: INC
      INTEGER(4),INTENT(IN)  :: JUMP
      INTEGER(4),INTENT(IN)  :: N
      INTEGER(4),INTENT(IN)  :: LOT
      INTEGER(4),INTENT(IN)  :: ISIGN
      
      INTEGER(4)  :: NFAX, NX21, NH
      INTEGER(4)  :: INK, IGO      
      INTEGER(4)  :: IBASE, JBASE      
      INTEGER(4)  :: I, J, K, L,M      
      INTEGER(4)  :: IA, IB, LA
            
!***********************************************************************      
!
      NFAX=IFAX(1)
      NX21=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
!
!>    @note IF NECESSARY, TRANSFER DATA TO WORK AREA
!>          Transfer data A to work area WORK 
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=INC+1
      JBASE=1
      DO 20 L=1,LOT  ! numbers of arrays.
         I=IBASE
         J=JBASE
!DIR$ IVDEP          
         DO 10 M=1,N ! numbers in one array
             WORK(J)=A(I)
             I=I+INC
             J=J+1
   10   CONTINUE
        IBASE=IBASE+JUMP
        JBASE=JBASE+NX21
   20 CONTINUE
!
      IGO=60
      GO TO 40
!
!>    @note PREPROCESSING (ISIGN=+1)
!>          ISIG=+1,inverse FFT, A TRANSFORM FROM FOURIER COEFFICIENTS TO
!>          GRIDPOINT VALUES 
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)  !default ISIGN=+1 (SPECTRAL TO GRIDPOINT TRANSFORM)
      IGO=60
!
!     COMPLEX TRANSFORM
!     -----------------
!
   40 CONTINUE
      IA=INC+1
      LA=1
      DO 80 K=1,NFAX
         IF (IGO.EQ.60) GO TO 60
   50    CONTINUE
         CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,     &
               INK,2,JUMP,NX21,LOT,NH,IFAX(K+1),LA)
         IGO=60
         GO TO 70
   60    CONTINUE
         CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,     &
               2,INK,NX21,JUMP,LOT,NH,IFAX(K+1),LA)
         IGO=50
   70    CONTINUE
          LA=LA*IFAX(K+1)
   80 CONTINUE
!
      IF (ISIGN.EQ.-1) GO TO 130
!
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=IA
      DO 100 L=1,LOT
         I=IBASE
         J=JBASE
!DIR$ IVDEP
         DO 90 M=1,N
            A(J)=WORK(I)
            I=I+1
            J=J+INC
   90   CONTINUE
        IBASE=IBASE+NX21
        JBASE=JBASE+JUMP
  100 CONTINUE
!
!     FILL IN CYCLIC BOUNDARY POINTS
  110 CONTINUE
      IA=1
      IB=N*INC+1
!DIR$ IVDEP
      DO 120 L=1,LOT
         A(IA)=A(IB)
         A(IB+INC)=A(IA+INC)
         IA=IA+JUMP
         IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!
!     POSTPROCESSING (ISIGN=-1):
!     --------------------------
!
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!
  140 CONTINUE
      RETURN
      END

