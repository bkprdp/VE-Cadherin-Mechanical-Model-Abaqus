C     USER SCRIPT CONTAINING RANDOM AMPLITUDE LOADING APPLIED AT THE CELL-CELL JUNCTIONS
      SUBROUTINE UAMP(
     1     ampName, time, ampValueOld, dt, nProps, props, nSvars,
     2     svars, lFlagsInfo,
     3     nSensor, sensorValues, sensorNames, jSensorLookUpTable,
     4     AmpValueNew,
     5     lFlagsDefine,
     6     AmpDerivative, AmpSecDerivative, AmpIncIntegral,
     7     AmpDoubleIntegral)

      INCLUDE 'ABA_PARAM.INC'

      character*80 ampName
      real*8 AmpValueNew

c      ! ----- variables for seed setting -----
      INTEGER :: i_seed
      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
      INTEGER, DIMENSION(1:8) :: dt_seed
c      ! ----- end of variables for seed setting -----

      if (ampName .eq. "AMP-1") then
        CALL RANDOM_SEED(size=i_seed)
        ALLOCATE(a_seed(1:i_seed))
        CALL RANDOM_SEED(get=a_seed)
        CALL DATE_AND_TIME(values=dt_seed)
        a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
        CALL RANDOM_SEED(put=a_seed)
        DEALLOCATE(a_seed)
        CALL RANDOM_NUMBER(AmpValueNew)
      else
        CALL RANDOM_SEED(size=i_seed)
        CALL RANDOM_NUMBER(AmpValueNew)
      end if

      if (AmpValueNew .ge. 0.5) THEN
        AmpValueNew = 0.5 - AmpValueNew
      ELSE
        AmpValueNew = AmpValueNew
      end if

      RETURN
      END