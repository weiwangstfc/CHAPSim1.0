SUBROUTINE CHT_MESH_GEO !master only
  USE CHT_info
  USE mehs_info
  USE init_info

  DY_SOLID_TOP = HY_Solid_Top / DBLE(NCL2_Solid_Top)
  DY_SOLID_BOT = HY_Solid_Bot / DBLE(NCL2_Solid_Bot)

  NCL1s_CHT_TOP = NINT(HXs_Solid_Top * DXI)
  NCL1s_CHT_BOT = NINT(HXs_Solid_Bot * DXI)

  NCL1e_CHT_TOP = NINT(HXe_Solid_Top * DXI)
  NCL1e_CHT_BOT = NINT(HXe_Solid_Bot * DXI)

END SUBROUTINE

!*******************************************************************************
SUBROUTINE CHT_THERMAL_preparation !master only
  USE CHT_info
  USE thermal_info

  !=====nondimensionalization of the thermal properties based on flow info===
  Cp_Solid_Top = Cp_Solid_Top / Cp0
  Cp_Solid_Bot = Cp_Solid_Bot / Cp0

  D_Solid_Top = D_Solid_Top / D0
  D_Solid_Bot = D_Solid_Bot / D0

  K_Solid_Top = K_Solid_Top / K0
  K_Solid_Bot = K_Solid_Bot / K0

  !==== Heat source /sink in the solid ========
  IF (iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux) THEN
      IF(iCase == ICHANL) THEN
          HEAT_SRC_SOLID_TOP = thermalWallBC_Dim(iTopWall) / HY_Solid_Top !unit: W/ M3
      END IF

      IF(iCase == IPIPEC) THEN
          HEAT_SRC_SOLID_TOP = 2.0_WP * thermalWallBC_Dim(iTopWall) * HYT / ((HYT + HY_Solid_Top)**2 - HYT**2)
      END IF

      IF(iCase == IANNUL) THEN
          HEAT_SRC_SOLID_TOP = 2.0_WP * thermalWallBC_Dim(iTopWall) * HYT / ((HYT + HY_Solid_Top)**2 - HYT**2)
      END IF

      HEAT_SRC_SOLID_TOP = HEAT_SRC_SOLID_TOP / (T0 * D0 * Cp0 * U0 / L0) !dimensionless
  END IF

  IF (iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
      IF(iCase == ICHANL) THEN
          HEAT_SRC_SOLID_BOT = thermalWallBC_Dim(iBotWall) / HY_Solid_BOT !unit: W/ M3
      END IF

      IF(iCase == IPIPEC) THEN
          HEAT_SRC_SOLID_BOT = 0.0_WP ! not used.
      END IF

      IF(iCase == IANNUL) THEN
          HEAT_SRC_SOLID_BOT = 2.0_WP * thermalWallBC_Dim(iBotWall) * HYB / (HYB**2 - (HYT - HY_Solid_BOT)**2)
      END IF

      HEAT_SRC_SOLID_BOT = HEAT_SRC_SOLID_BOT / (T0 * D0 * Cp0 * U0 / L0) !dimensionless
  END IF

  IF (iThermalWallType(iTopWall) == IsoPowerDENSITY) THEN
      HEAT_SRC_SOLID_TOP = thermalWallBC_Dim(iTopWall) / (T0 * D0 * Cp0 * U0 / L0) !dimensionless
  END IF

  IF (iThermalWallType(iBotWall) == IsoPowerDENSITY) THEN
      HEAT_SRC_SOLID_BOT = thermalWallBC_Dim(iBotWall) / (T0 * D0 * Cp0 * U0 / L0) !dimensionless
  END IF

END SUBROUTINE

!*******************************************************************************
SUBROUTINE BCAST_CHT

  CALL MPI_BCAST( NCL1s_CHT_TOP,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
  CALL MPI_BCAST( NCL1s_CHT_BOT,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
  CALL MPI_BCAST( NCL1e_CHT_TOP,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
  CALL MPI_BCAST( NCL1e_CHT_BOT,    1, MPI_INTEGER4,         0, ICOMM, IERROR )

  CALL MPI_BCAST( Cp_Solid_Top,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
  CALL MPI_BCAST( Cp_Solid_Bot,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
  CALL MPI_BCAST( D_Solid_Top,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
  CALL MPI_BCAST( D_Solid_Bot,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
  CALL MPI_BCAST( K_Solid_Top, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
  CALL MPI_BCAST( K_Solid_Bot, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

  CALL MPI_BCAST( HEAT_SRC_SOLID_TOP, 1,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
  CALL MPI_BCAST( HEAT_SRC_SOLID_BOT, 1,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

  CALL MPI_BCAST( DY_SOLID_TOP, 1,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
  CALL MPI_BCAST( DY_SOLID_BOT, 1,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

END SUBROUTINE

!*******************************************************************************
SUBROUTINE CHT_pARtitioned_ISS_Coupling
  USE CHT_info
  USE mehs_info
  USE init_info

  !========== Heat Flux from flow to solid at the interface (i, J', K) ============





END SUBROUTINE

SUBROUTINE CHT_Solid_HT_SOLVER



END SUBROUTINE
