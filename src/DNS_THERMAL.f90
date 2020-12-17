!**********************************************************************************************************************************
!                         Channel And Pipe flow Simulation (CHAPSim)
!**********************************************************************************************************************************
! History: Code Development History
!          The early version of CHAPSim was developed by Dr Mehdi Seddighi-Moornani in f77 platform for simulations of unsteady 
!          turbulent flows in a channel during his PhD study starting from 2007. The code has later been revamped by Dr Wei Wang 
!          for solving thermo-fluids problems in f90 platform in 2014/2015. All this early development work was carried out while
!          Drs Seddighi-Moornani and Wang were working with Professor Shuisheng He at Universities of Aberdeen and Sheffield. In 
!          2020, CHAPSim was selected by CCP.NTH to be developed as an open-source community code in partnership with STFC 
!          Daresbury Laboratory (Contact: Dr Wei Wang), Liverpool John Moores University (Contact: Dr Mehdi Seddighi) and
!          University of Sheffield (Contact: Professor Shuisheng He) and released under the GNU General Public License, that is a 
!          free, copyleft license.
! License: GNU General Public License
!          This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
!          License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later 
!          version. This program is distributed in the hope that it will be useful, but without any warranty; without even the 
!          implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more 
!          details. You should have received a copy of the GNU General Public License along with this program; If not, write to the 
!          Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add information on 
!          how to contact you by electronic and paper mail.
! Contact: Main contributors
!          Wei Wang (wei.wang@stfc.ac.uk, wei.wang@sheffield.ac.uk), STFC Daresbury Laboratory
!          Mehdi Seddighi (m.seddighi@lJMu.ac.uk), Liverpool John Moors University
!          Shuisheng He (s.he@sheffield.ac.uk), The University Of Sheffield
!**********************************************************************************************************************************
! DESCRIPTION:
!> @brief
!> Code Name: Channel And Pipe flow Simulation (CHAPSim)
!> @version  CHAPSim_noCHT_2020
!> @par        SIMULATION CODE(DOUBLE PRECISION)
!>             PARALLEL
!>                    - MPI                                          \n
!>             MESH                                                  \n
!>                    - STAGGERED GRID                               \n
!>             NONLINEAR                                             \n
!>                    - CONVECTION & DIVERGENCE FORM                 \n
!>                    - 2ND ORDER FINITE DIFFERENCE METHOD           \n
!>                    - RUNGE KUTTA & ADAMS-BASHFORTH METHOD         \n
!>             VISCOUS                                               \n
!>                    - CRANK-NICOLSON METHOD                        \n
!>             PRESSURE                                              \n
!>                    - FFT                                          \n
!>                    - FRACTIONAL STEP METHOD                       \n
!>             THERMODYNAMIC                                         \n
!>                    - QUAS-IMCOMPRESSIBLE FLOW                     \n
!>                    - THE SOLVED THERMAL VARIABLE IS $\RHO * H$    \n
!**********************************************************************************************************************************
!> REVISION HISTORY:
!> 11/2012 - Initial Version (tg domain only), by Mehdi Seddighi
!> 05/2014 - Added cylindrical coordinate for pipe flow, by Kui He.
!> 02/2014 - Added thermodynamics simulation and carried out code optimization, by Wei Wang
!**********************************************************************************************************************************

PROGRAM CHAPSim_DNS_Solver
    USE init_info
    IMPLICIT NONE

    !====PREPARE DIR PATH AND FILE NAME ==========
    CALL MKDIR
    CALL Start_MPI
    IF (MYID == 0) THEN
        CALL date_and_time(DATE = date, TIME = time)
        fllog = date(1:4) // '.' // Date(5:8) // '.' // Time(1:4) // '.log'
        OPEN(logflg_pg, FILE = TRIM(FilePath0) // 'CHAPSim_Solver.' //fllog)
        CALL CHKHDL('**** CHAPSim DNS solver Starts******', MYID)
    END IF

    !=== READ IN USER-GIVEN PARAMETERS ===============
    IF (MYID == 0) CALL CHKHDL('1. Read ini file...', MYID)
    IF (MYID == 0) &
    CALL READINI
    CALL BCAST_READINI

    IF (MYID == 0) CALL CHKHDL('2. Mesh decomposition (both y-dir and z-dir)...', MYID)
    CALL mesh_Ydecomp
    CALL mesh_Zdecomp

    IF (MYID == 0) CALL CHKHDL('3. Allocate variables...', MYID)
    CALL MEM_ALLOCAT

    IF (MYID == 0) CALL CHKHDL('4. Set up global and local index for MPI ...', MYID)
    CALL INDEXL2G

    IF (MYID == 0) CALL CHKHDL('5. Set up + 1 /- 1 index (all three directions) ...', MYID)
    CALL INDEXSET

    IF (MYID == 0) CALL CHKHDL('6. Set up RK coefficients...', MYID)
    IF (MYID == 0) &
    CALL CONS_RKCOEF
    CALL BCAST_RKCOEF

    IF (MYID == 0) CALL CHKHDL('7. Set up mesh information in x and z directions...', MYID)
    IF (MYID == 0) &
    CALL CONSPARA
    CALL BCAST_CONSPARA

    IF (MYID == 0) CALL CHKHDL('8. Set up mesh distribution in y direction...', MYID)
    IF (MYID == 0) &
    CALL COORDJR
    CALL BCAST_COORDJR

    IF (MYID == 0) CALL CHKHDL('9. Set up laminar (initial) velocity profile...', MYID)
    IF (MYID == 0) &
    CALL LAMPOISLPROF
    CALL BCAST_LAMPOISLPROF

    IF (MYID == 0) CALL CHKHDL('10. Set up mesh coefficient for TDMA...', MYID)
    IF (MYID == 0) &
    CALL TDMA_COEF
    CALL BCAST_TDMA_COEF
    !CALL CFL_VISCOUS ! \nu dependent

    !====PREPARE POISSON SOLVER ==========
    IF (MYID == 0) CALL CHKHDL('11. Initialization of the FFT solver...', MYID)
    IF(TgFlowFlg .AND. IoFlowFlg) THEN
        CALL FFT99_POIS3D_INIT(ITG)
        CALL FISHPACK_POIS3D_INIT
    ELSE
        IF(TgFlowFlg) CALL FFT99_POIS3D_INIT(ITG)
        IF(IoFlowFlg) THEN
            CALL FFT99_POIS3D_INIT(IIO) !Method One
            !CALL FISHPACK_POIS3D_INIT  !Metthod Two, good
        END IF
    END IF

    CALL MPI_BARRIER(ICOMM, IERROR)

    !====PREPARE THERMAL INFO============
    IF(IoFlowFlg) THEN
        IF(iThermoDynamics == 1) THEN
            IF (MYID == 0) CALL CHKHDL('12. Initialization of thermal properties...', MYID)
            IF (MYID == 0) CALL thermal_init
            CALL BCAST_THERMAL_INIT

            CALL MPI_BARRIER(ICOMM, IERROR)
        ELSE
            IF (MYID == 0) CALL CHKHDL('12. No thermodynamics is considered.', MYID)
        END IF
    END IF

    IF(iPostProcess == 2) THEN
        CALL POSTPROCESS_INTEGRAL_INSTANS
        CALL MPI_BARRIER(ICOMM, IERROR)
        IF(MYID == 0) CALL CHKHDL('<===Only postprocessed given instantanous results, now the code stops...==>', MYID)
        STOP
    END IF

    !==== The main CFD Solver ==========
    IF (MYID == 0) CALL CHKHDL('13. The main CFD solver Starts...', MYID)
    CALL SOLVE

    CALL MEM_DEALLOCAT
    IF (MYID == 0) CALL CHKHDL('15. The CFD solver finished.', MYID)
    CALL MPI_FINALIZE(IERROR)
    STOP

END PROGRAM CHAPSim_DNS_Solver
