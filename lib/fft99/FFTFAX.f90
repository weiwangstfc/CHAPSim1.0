!>@par    SUBROUTINE FFTFAX
!>                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
!>                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
!>                (PROVIDED THAT N IS NOT CHANGED).
!-----------------------------------------------------------------------
!>@par    SUBROUTINE FFTFAX (N,IFAX,TRIGS)
!> 
!>@par    PURPOSE
!>              A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE
!>              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT
!>              ROUTINES (PROVIDED THAT N IS NOT CHANGED).
!
!>@par ARGUMENT DIMENSIONS   
!>             IFAX(13),TRIGS(3*N/2+1) 
!
!>@par ARGUMENTS
!>@par ON INPUT
!>        N      \n
!>             AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR
!>             GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
!>             THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE
!>             DEFINITIONS OF THE TRANSFORMS).
!>             theoretically, N=2^N1*3^N2*5^N3, an even number larger than 4. 
!
!>        IFAX    \n
!>               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
!>               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
!>               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.
!>        IFAX(1) = how many prime factors in N.
!>        IFAX(2:IFAX(1)+1) prime factors  based on 3,4,5 
!>        N=2*IFAX(2)*IFAX(3)*IFAX(4)...IFAX(IFAX(1)+1)
!
!>        TRIGS    \n
!>               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
!>               EVEN, OR 3*N/2+1 IF N/2 IS ODD.
!>@par
!> ON OUTPUT    IFAX
!               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE
!               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
!               IN IFAX(2),IFAX(3),...  IF FFTFAX IS CALLED WITH N ODD,
!               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)
!               IS SET TO -99.
!
!              TRIGS
!               AN ARRAY OF TRIGNOMENTRIC FUNCTION VALUES SUBSEQUENTLY
!               USED BY THE FFT ROUTINES.
!
!-----------------------------------------------------------------------
!> @note the original code is given by below link
!>        http://map.nasa.gov/GEOSgcm_f90toHTML/html_code/FVdycore_GridComp/fft99.F90.html
!> @param N     [in]   Cell no. in one direction
!> @param IFAX  [out]  given by subroutine FAX
!> @param TRIGS [out] 

    
      SUBROUTINE FFTFAX(N,IFAX,TRIGS)
      USE WPRECISIONFFT
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     INCLUDE 'mpif.h'   
!      DIMENSION IFAX(13),TRIGS(*)
      !INTEGER,PARAMETER  :: WP=KIND(0.0D0) !WORKING PRECESION
      REAL(WP),INTENT(INOUT)      :: TRIGS(*)
      INTEGER(4),INTENT(INOUT)  :: IFAX(13)
      INTEGER(4),INTENT(IN)      :: N
      
      INTEGER(4)  :: MODE=3
      INTEGER(4)  :: I
! MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
! TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
! DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
! WAS WRITTEN.
!
!     DATA MODE /3/
!**************************************************************      
      CALL FAX (IFAX, N, MODE)
      I = IFAX(1)
      IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99
      
!     IF (IFAX(1) .LE. 0 )CALL ULIBER(33, ' FFTFAX -- INVALID N', 20)
      IF (IFAX(1) .LE. 0 ) THEN
         WRITE(*,*) ' FFTFAX -- INVALID N'
         STOP
      ENDIF
      
!>    @warning if MODE==3, below subroutine is just return without anything.      
      CALL FFTRIG (TRIGS, N, MODE)
      
      RETURN
      
      END
