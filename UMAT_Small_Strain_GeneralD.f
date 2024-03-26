c The ABAQUS UMAT subroutine that describes the cell as a strain-rate dependent material.
c A version of this code is primarliy used to study endothelial permeability as described in the article
c "P.Keshavanarayana, F.Spill, A mechanical modelling framework to study endothelial permeability, Biophysical Journal, 2024." 

c If you spot any bug in the code, please let the authors know. It will help in further development. 
c Authors : Pradeep Keshavanarayana, Fabian Spill
c Contact : p.keshavanarayana@bham.ac.uk, f.spill@bham.ac.uk
c 
c The code is distributed with CC BY license. Please cite the above article if you are using this code. 
c
c      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(NDI),DROT(NDI,NDI),
     4 DFGRD0(NDI,NDI),DFGRD1(NDI,NDI)

C     DEFINE VARIABLES   
      REAL(KIND=8) E_PASSIVE,NU,ID_MAT(NDI,NDI),K_F,K_B,SIG_MAX,K_T
      REAL(KIND=8) SIGMA_MAT(NDI,NDI),DDSDDE_CONS(NTENS,NTENS)
      REAL(KIND=8) C_MAT(NTENS,NTENS),ALAMBDA,AMU
      INTEGER N_FIBRES, I, J, N_FIBRES_EFF 

C     MATERIAL PARAMETERS
      E_PASSIVE=PROPS(1)
      NU=PROPS(2)
      SIG_MAX=PROPS(3)
      K_T = PROPS(4)
      K_F = PROPS(5)
      K_B = PROPS(6)

C     NUMBER OF FIBRES IN CYTOPLASM (N_FIBRES)
      IF ( NDI .EQ. 2) THEN
            N_FIBRES = 20
      ELSE IF (NDI .EQ. 3) THEN
            N_FIBRES = 20*20
            N_FIBRES_EFF = 20
      END IF

C     INITIALISE MATRICES
      DO I=1,NDI
        DO J=1,NDI
          ID_MAT(I,J) = 0.0D0
          SIGMA_MAT(I,J) = 0.0D0
          IF (I .EQ. J) THEN
            ID_MAT(I,J) = 1.0D0
          END IF
        END DO
      END DO

      IF (TIME(1) .EQ. 0.0) THEN
         DO I=1,NSTATV
            STATEV(I) = 0.0D0
         END DO
      END IF

      DO I=1,NTENS
        DO J=1,NTENS
            C_MAT(I,J) = 0.0D0
        END DO
      END DO

C     CONVERTING YOUNGS MODULUS (E) AND POISSON RATIO (NU) TO LAME CONSTANTS (ALAMBDA,AMU)
      ALAMBDA = E_PASSIVE*NU/((1.0D0+NU)*(1.0D0-2.0D0*NU))
      AMU = E_PASSIVE/(2.0D0*(1.0D0+NU))

C     STIFFNESS MATRIX FOR A 2D ELASTIC MATERIAL
      IF ( NDI .EQ. 2) THEN
            C_MAT(1,1) = 2.0D0*AMU + ALAMBDA
            C_MAT(1,2) = ALAMBDA
            C_MAT(1,3) = 0.0D0
            C_MAT(2,1) = ALAMBDA
            C_MAT(2,2) = 2.0D0*AMU + ALAMBDA
            C_MAT(2,3) = 0.0D0
            C_MAT(3,1) = 0.0D0
            C_MAT(3,2) = 0.0D0
            C_MAT(3,3) = AMU
      ELSE
            C_MAT(1,1) = 2.0D0*AMU + ALAMBDA
            C_MAT(1,2) = ALAMBDA
            C_MAT(1,3) = ALAMBDA
            C_MAT(2,1) = ALAMBDA
            C_MAT(2,2) = 2.0D0*AMU + ALAMBDA
            C_MAT(2,3) = ALAMBDA
            C_MAT(3,1) = ALAMBDA
            C_MAT(3,2) = ALAMBDA
            C_MAT(3,3) = 2.0D0*AMU + ALAMBDA
            C_MAT(4,4) = AMU
            C_MAT(5,5) = AMU
            C_MAT(6,6) = AMU
      END IF

C     EVALUATE TOTAL STRESS IN THE CELL
      CALL EVAL_SIGMA(E_PASSIVE,NU,ID_MAT,N_FIBRES,STATEV,NSTATV,
     1     DSTRAN,DTIME,DROT,K_F,K_B,SIG_MAX,K_T,C_MAT,SIGMA_MAT,
     2     NDI,NTENS,N_FIBRES_EFF)

C     RETURN UPDATED TOTAL STRESS TO ABAQUS
      IF ( NDI .EQ. 2) THEN
            STRESS(1) = SIGMA_MAT(1,1)
            STRESS(2) = SIGMA_MAT(2,2)
            STRESS(3) = SIGMA_MAT(1,2)
      ELSE 
            STRESS(1) = SIGMA_MAT(1,1)
            STRESS(2) = SIGMA_MAT(2,2)
            STRESS(3) = SIGMA_MAT(3,3)
            STRESS(4) = SIGMA_MAT(2,3)
            STRESS(5) = SIGMA_MAT(1,3)
            STRESS(6) = SIGMA_MAT(1,2)
      END IF
      
C     RETURN DDSDDE MATRIX TO ABAQUS
      DDSDDE = C_MAT

      RETURN
      END SUBROUTINE

C     DEFINITION OF TOTAL STRESS EVALUATION SUBROUTINE
      SUBROUTINE EVAL_SIGMA(E_PASSIVE,NU,ID_MAT,N_FIBRES,STATEV,
     1           NSTATV,DSTRAN,DTIME,DROT,K_F,K_B,SIG_MAX,K_T,
     2           C_MAT,SIGMA_MAT,NDI,NTENS,N_FIBRES_EFF)

      IMPLICIT NONE
      INTEGER N_FIBRES, I, J, NSTATV, NDI,NTENS, N_FIBRES_EFF
      REAL(KIND=8) SIGMA_MAT(NDI,NDI),SIGMA_ACTIVE_MAT(NDI,NDI)
      REAL(KIND=8) SIGMA_PASSIVE_MAT(NDI,NDI)
      REAL(KIND=8) DSTRAN(NTENS),E_PASSIVE,NU,ID_MAT(NDI,NDI),STATEV(NSTATV)
      REAL(KIND=8) DTIME,K_F,K_B,SIG_MAX,K_T,C_MAT(NTENS,NTENS),DROT(NTENS,NTENS)

      DO I=1,NDI
        DO J=1,NDI
          SIGMA_PASSIVE_MAT(I,J) = 0.0D0
          SIGMA_ACTIVE_MAT(I,J) = 0.00D0
        END DO
      END DO

C     CALL SUBROUTINE TO EVALUATE PASSIVE STRESS
      CALL EVAL_SIGMA_PASSIVE_ELASTIC(E_PASSIVE,NU,DSTRAN,N_FIBRES,
     1                       STATEV,NSTATV,C_MAT,SIGMA_PASSIVE_MAT,
     2                       NDI,NTENS)

C     CALL SUBROUTINE TO EVALUATE ACTIVE STRESS
      CALL EVAL_SIGMA_ACTIVE(N_FIBRES,STATEV,NSTATV,DSTRAN,
     1                       DTIME, DROT, K_F,K_B,SIG_MAX,K_T,
     2                       SIGMA_ACTIVE_MAT,NDI,NTENS,N_FIBRES_EFF)

C     TOTAL STRESS = ACTIVE_STRESS + PASSIVE_STRESS
      SIGMA_MAT = SIGMA_PASSIVE_MAT + SIGMA_ACTIVE_MAT

      RETURN
      END SUBROUTINE

C     DEFINITION OF PASSIVE STRESS EVALUATION SUBROUTINE
      SUBROUTINE EVAL_SIGMA_PASSIVE_ELASTIC(E_PASSIVE,NU,DSTRAN,
     1           N_FIBRES,STATEV,NSTATV,C_MAT,SIGMA_PASSIVE_MAT,
     2           NDI,NTENS)
      
      IMPLICIT NONE  
      INTEGER N_FIBRES,I,J, NSTATV,NDI,NTENS
      REAL(KIND=8) E_PASSIVE, NU, C_MAT(NTENS,NTENS), STATEV(NSTATV)
      REAL(KIND=8) SIGMA_PASSIVE_VEC(NTENS), DSTRAN(NTENS)
      REAL(KIND=8) SIGMA_PASSIVE_MAT(NDI,NDI), AMU, ALAMBDA

C     PASSIVE STRESS FROM PREVIOUS TIME STEP
      IF (NDI .EQ. 2) THEN
            SIGMA_PASSIVE_VEC(1) = STATEV(2*N_FIBRES+1)
            SIGMA_PASSIVE_VEC(2) = STATEV(2*N_FIBRES+2)
            SIGMA_PASSIVE_VEC(3) = STATEV(2*N_FIBRES+3)
      ELSE
            SIGMA_PASSIVE_VEC(1) = STATEV(2*N_FIBRES+1)
            SIGMA_PASSIVE_VEC(2) = STATEV(2*N_FIBRES+2)
            SIGMA_PASSIVE_VEC(3) = STATEV(2*N_FIBRES+3)
            SIGMA_PASSIVE_VEC(4) = STATEV(2*N_FIBRES+4)
            SIGMA_PASSIVE_VEC(5) = STATEV(2*N_FIBRES+5)
            SIGMA_PASSIVE_VEC(6) = STATEV(2*N_FIBRES+6)
      END IF

C     UPDATE PASSIVE STRESS AS PASSIVE_STRESS = PASSIVE_STRESS + C*DELTA_EPSILON
      DO I=1,NTENS
        DO J=1,NTENS
          SIGMA_PASSIVE_VEC(I) = SIGMA_PASSIVE_VEC(I) + 
     1       C_MAT(I,J)*DSTRAN(J)
        END DO
      END DO

C     UPDATE STATEV WITH THE UPDATED PASSIVE STRESS
      IF (NDI .EQ. 2) THEN
            STATEV(2*N_FIBRES+1) = SIGMA_PASSIVE_VEC(1)
            STATEV(2*N_FIBRES+2) = SIGMA_PASSIVE_VEC(2)
            STATEV(2*N_FIBRES+3) = SIGMA_PASSIVE_VEC(3)
      ELSE
            STATEV(2*N_FIBRES+1) = SIGMA_PASSIVE_VEC(1)
            STATEV(2*N_FIBRES+2) = SIGMA_PASSIVE_VEC(2)
            STATEV(2*N_FIBRES+3) = SIGMA_PASSIVE_VEC(3)
            STATEV(2*N_FIBRES+4) = SIGMA_PASSIVE_VEC(4)
            STATEV(2*N_FIBRES+5) = SIGMA_PASSIVE_VEC(5)
            STATEV(2*N_FIBRES+6) = SIGMA_PASSIVE_VEC(6)
      END IF

C     REWRITE PASSIVE STRESS VECTOR IN THE FORM OF A MATRIX
      IF (NDI .EQ. 2) THEN
            SIGMA_PASSIVE_MAT(1,1) = SIGMA_PASSIVE_VEC(1)
            SIGMA_PASSIVE_MAT(2,2) = SIGMA_PASSIVE_VEC(2)
            SIGMA_PASSIVE_MAT(1,2) = SIGMA_PASSIVE_VEC(3)
            SIGMA_PASSIVE_MAT(2,1) = SIGMA_PASSIVE_MAT(1,2)
      ELSE      
            SIGMA_PASSIVE_MAT(1,1) = SIGMA_PASSIVE_VEC(1)
            SIGMA_PASSIVE_MAT(2,2) = SIGMA_PASSIVE_VEC(2)
            SIGMA_PASSIVE_MAT(3,3) = SIGMA_PASSIVE_VEC(3)
            SIGMA_PASSIVE_MAT(2,3) = SIGMA_PASSIVE_VEC(4)
            SIGMA_PASSIVE_MAT(1,3) = SIGMA_PASSIVE_VEC(5)
            SIGMA_PASSIVE_MAT(1,2) = SIGMA_PASSIVE_VEC(6)
            SIGMA_PASSIVE_MAT(2,1) = SIGMA_PASSIVE_MAT(1,2)
            SIGMA_PASSIVE_MAT(3,1) = SIGMA_PASSIVE_MAT(1,3)
            SIGMA_PASSIVE_MAT(3,2) = SIGMA_PASSIVE_MAT(2,3)
      END IF

      RETURN
      END SUBROUTINE

C     DEFINITION OF ACTIVE STRESS EVALUATION SUBROUTINE
      SUBROUTINE EVAL_SIGMA_ACTIVE(N_FIBRES,STATEV,NSTATV,DSTRAN,
     1           DTIME, DROT, K_F,K_B,SIG_MAX,K_T,SIGMA_ACTIVE_MAT,
     2           NDI,NTENS, N_FIBRES_EFF)                           

      IMPLICIT NONE
      INTEGER N_FIBRES, NSTATV, NDI,NTENS
      REAL(KIND=8) SIGMA_ACTIVE_MAT(NDI,NDI),DTIME,K_F,K_B,SIG_MAX
      REAL(KIND=8) STATEV(NSTATV),K_T,DROT(NDI,NDI), TIME(2)
      REAL(KIND=8) DPHI,DOMEGA,PHI_VEC(N_FIBRES),OMEGA_VEC(N_FIBRES)
      REAL(KIND=8) SIGMA_ACTIVE_THETA_OLD(N_FIBRES)
      REAL(KIND=8) EPSDOT_MAT(NDI,NDI)
      REAL(KIND=8) M_VEC(NDI),EPSDOT_THETA(N_FIBRES)
      REAL(KIND=8) ETA_SF_THETA(N_FIBRES),SIGMA_ACTIVE_THETA(N_FIBRES)
      REAL(KIND=8) EPSDOT_VEC(NTENS),DSTRAN(NTENS)
      REAL(KIND=8) PI
      INTEGER I,J,N_FIBRES_EFF,COUNT_FIBRE

      PI=4.D0*ATAN(1.D0)

C     DISTRIBUTION OF FIBRES IN THE SPHERE
      IF (NDI .EQ. 2) THEN
            DO I=1,N_FIBRES
                  PHI_VEC(I) = 0.0D0 + (2.0D0*PI/N_FIBRES)*(I-1)
            END DO
      ELSE 
            DO I=1,N_FIBRES_EFF
                  PHI_VEC(I) = 0.0D0 + (2.0D0*PI/N_FIBRES_EFF)*(I-1)
                  OMEGA_VEC(I) = 0.0D0 + (PI/N_FIBRES_EFF)*(I-1)
            END DO
      END IF

C     INITIALISE ACTIVE STRESS AND STRAIN RATE MATRICES (NEEDED FOR NUMERICAL INTEGRATION)
      DO I=1,NDI
        DO J=1,NDI
          EPSDOT_MAT(I,J) = 0.0D0
          SIGMA_ACTIVE_MAT(I,J) = 0.0D0
        END DO
      END DO

C     ACTIVE STRESS FROM THE PREVIOUS TIME STEP      
      DO I=1,N_FIBRES
        SIGMA_ACTIVE_THETA_OLD(I) = STATEV(N_FIBRES+I)
      END DO

C     EVALUATE STRAIN_RATE = DELTA_EPSILON/DELTA_T
      IF (NDI .EQ. 2) THEN
            EPSDOT_VEC(1) = DSTRAN(1)/DTIME
            EPSDOT_VEC(2) = DSTRAN(2)/DTIME
            EPSDOT_VEC(3) = DSTRAN(3)/DTIME  
      ELSE
            EPSDOT_VEC(1) = DSTRAN(1)/DTIME
            EPSDOT_VEC(2) = DSTRAN(2)/DTIME
            EPSDOT_VEC(3) = DSTRAN(3)/DTIME
            EPSDOT_VEC(4) = DSTRAN(4)/DTIME
            EPSDOT_VEC(5) = DSTRAN(5)/DTIME
            EPSDOT_VEC(6) = DSTRAN(6)/DTIME
      END IF

C     REWRITE STRAIN RATE VECTOR IN THE FORM OF A MATRIX
      IF (NDI .EQ. 2) THEN
            EPSDOT_MAT(1,1) = EPSDOT_VEC(1)
            EPSDOT_MAT(2,2) = EPSDOT_VEC(2)
            EPSDOT_MAT(1,2) = (1.0/2.0)*EPSDOT_VEC(3)
            EPSDOT_MAT(2,1) = EPSDOT_MAT(1,2)
      ELSE
            EPSDOT_MAT(1,1) = EPSDOT_VEC(1)
            EPSDOT_MAT(2,2) = EPSDOT_VEC(2)
            EPSDOT_MAT(3,3) = EPSDOT_VEC(3)
            EPSDOT_MAT(2,3) = (1.0/2.0)*EPSDOT_VEC(4)
            EPSDOT_MAT(1,3) = (1.0/2.0)*EPSDOT_VEC(5)
            EPSDOT_MAT(1,2) = (1.0/2.0)*EPSDOT_VEC(6)
            EPSDOT_MAT(2,1) = EPSDOT_MAT(1,2)
            EPSDOT_MAT(3,1) = EPSDOT_MAT(1,3)
            EPSDOT_MAT(3,2) = EPSDOT_MAT(2,3)
      END IF

C     EVALUATE STRAIN RATE IN EACH FIBRE FROM STRAIN RATE MATRIX
      IF (NDI .EQ. 2) THEN
            DO I=1,N_FIBRES
                  M_VEC(1) = COS(PHI_VEC(I))
                  M_VEC(2) = SIN(PHI_VEC(I))
                  EPSDOT_THETA(I) = EPSDOT_MAT(1,1)*M_VEC(1)*M_VEC(1)
     1                  + EPSDOT_MAT(2,2)*M_VEC(2)*M_VEC(2)
     2                  + EPSDOT_MAT(1,2)*M_VEC(1)*M_VEC(2)
     3                  + EPSDOT_MAT(2,1)*M_VEC(2)*M_VEC(1)
            END DO
      ELSE IF (NDI .EQ. 3) THEN
            COUNT_FIBRE = 1
            DO I=1,N_FIBRES_EFF
                  DO J=1,N_FIBRES_EFF
                        M_VEC(1) = SIN(OMEGA_VEC(J)*COS(PHI_VEC(I)))
                        M_VEC(2) = SIN(OMEGA_VEC(J)*SIN(PHI_VEC(I)))
                        M_VEC(3) = COS(OMEGA_VEC(J))

                        EPSDOT_THETA(COUNT_FIBRE) = 
     1                    EPSDOT_MAT(1,1)*M_VEC(1)*M_VEC(1)
     2                  + EPSDOT_MAT(2,2)*M_VEC(2)*M_VEC(2)
     3                  + EPSDOT_MAT(3,3)*M_VEC(3)*M_VEC(3)
     4                  + EPSDOT_MAT(1,2)*M_VEC(1)*M_VEC(2)
     5                  + EPSDOT_MAT(1,3)*M_VEC(1)*M_VEC(3)
     6                  + EPSDOT_MAT(2,3)*M_VEC(2)*M_VEC(3)
     7                  + EPSDOT_MAT(2,1)*M_VEC(2)*M_VEC(1)
     8                  + EPSDOT_MAT(3,1)*M_VEC(3)*M_VEC(1)
     9                  + EPSDOT_MAT(3,2)*M_VEC(3)*M_VEC(2)

                        COUNT_FIBRE = COUNT_FIBRE + 1
                  END DO
            END DO
      END IF

C     EVALUATE STRESS FIBRE CONCENTRATION IN EACH FIBRE
      CALL EVAL_ETA_SF_THETA(K_F,K_B,SIG_MAX,N_FIBRES,DTIME,NSTATV,
     1            STATEV,SIGMA_ACTIVE_THETA_OLD, ETA_SF_THETA)

C     UPDATE STATEV WITH UPDATED STRESS FIBRE CONCENTRATION IN EACH FIBRE
      DO I=1,N_FIBRES
            STATEV(I) = ETA_SF_THETA(I)
      END DO
      
C     CALL SUBROUTINE TO EVALUATE ACTIVE STRESS IN EACH FIBRE
      CALL EVAL_ACTIVE_STRESS_THETA(K_T,SIG_MAX,N_FIBRES,STATEV,
     1     NSTATV,EPSDOT_THETA,SIGMA_ACTIVE_THETA,PHI_VEC)

C     UPDATE STATEV WITH UPDATED ACTIVE STRESS IN EACH FIBRE
      DO I=1,N_FIBRES
            STATEV(N_FIBRES + I) = SIGMA_ACTIVE_THETA(I)
      END DO

      IF (NDI .EQ. 2) THEN
            DPHI = 2.0D0*PI/N_FIBRES
      ELSE IF (NDI .EQ. 3) THEN
            DPHI = 2.0D0*PI/N_FIBRES_EFF
            DOMEGA = PI/N_FIBRES_EFF
      END IF

C     INTEGRATION OF ACTIVE STRESS IN EACH FIBRE TO EVALUATE HOMOGENISED ACTIVE STRESS MATRIX
      IF ( NDI .EQ. 2) THEN
            DO I=1,N_FIBRES
                  M_VEC(1) = COS(PHI_VEC(I))
                  M_VEC(2) = SIN(PHI_VEC(I))
                  SIGMA_ACTIVE_MAT(1,1) = SIGMA_ACTIVE_MAT(1,1) + 
     1                            SIGMA_ACTIVE_THETA(I)
     2            *M_VEC(1)*M_VEC(1)*DPHI
                  SIGMA_ACTIVE_MAT(2,2) = SIGMA_ACTIVE_MAT(2,2) + 
     1                             SIGMA_ACTIVE_THETA(I)  
     2            *M_VEC(2)*M_VEC(2)*DPHI
                  SIGMA_ACTIVE_MAT(1,2) = SIGMA_ACTIVE_MAT(1,2) + 
     1                            SIGMA_ACTIVE_THETA(I)
     2            *M_VEC(1)*M_VEC(2)*DPHI          
                  
            END DO
            SIGMA_ACTIVE_MAT(1,1) = (1.0D0/PI)*SIGMA_ACTIVE_MAT(1,1)
            SIGMA_ACTIVE_MAT(2,2) = (1.0D0/PI)*SIGMA_ACTIVE_MAT(2,2)
            SIGMA_ACTIVE_MAT(1,2) = (1.0D0/PI)*SIGMA_ACTIVE_MAT(1,2)
            SIGMA_ACTIVE_MAT(2,1) = SIGMA_ACTIVE_MAT(1,2)
      ELSE
            COUNT_FIBRE = 1
            DO I=1,N_FIBRES_EFF  
                  DO J=1,N_FIBRES_EFF
                        M_VEC(1) = SIN(OMEGA_VEC(J)*COS(PHI_VEC(I)))
                        M_VEC(2) = SIN(OMEGA_VEC(J)*SIN(PHI_VEC(I)))
                        M_VEC(3) = COS(OMEGA_VEC(J))

                        SIGMA_ACTIVE_MAT(1,1) = SIGMA_ACTIVE_MAT(1,1) 
     1                           + SIGMA_ACTIVE_THETA(COUNT_FIBRE)
     2            *M_VEC(1)*M_VEC(1)*SIN(OMEGA_VEC(J))*DPHI*DOMEGA
                        SIGMA_ACTIVE_MAT(2,2) = SIGMA_ACTIVE_MAT(2,2)
     1                           + SIGMA_ACTIVE_THETA(COUNT_FIBRE)  
     2            *M_VEC(2)*M_VEC(2)*SIN(OMEGA_VEC(J))*DPHI*DOMEGA
                        SIGMA_ACTIVE_MAT(3,3) = SIGMA_ACTIVE_MAT(3,3)
     1                           + SIGMA_ACTIVE_THETA(COUNT_FIBRE)
     2            *M_VEC(3)*M_VEC(3)*SIN(OMEGA_VEC(J))*DPHI*DOMEGA
                        SIGMA_ACTIVE_MAT(1,2) = SIGMA_ACTIVE_MAT(1,2)
     1                           + SIGMA_ACTIVE_THETA(COUNT_FIBRE)
     2            *M_VEC(1)*M_VEC(2)*SIN(OMEGA_VEC(J))*DPHI*DOMEGA
                        SIGMA_ACTIVE_MAT(1,3) = SIGMA_ACTIVE_MAT(1,3)
     1                           + SIGMA_ACTIVE_THETA(COUNT_FIBRE)
     2            *M_VEC(1)*M_VEC(3)*SIN(OMEGA_VEC(J))*DPHI*DOMEGA
                        SIGMA_ACTIVE_MAT(2,3) = SIGMA_ACTIVE_MAT(2,3)
     1                           + SIGMA_ACTIVE_THETA(COUNT_FIBRE)
     2            *M_VEC(2)*M_VEC(3)*SIN(OMEGA_VEC(J))*DPHI*DOMEGA
                  
                        COUNT_FIBRE = COUNT_FIBRE + 1  
                  END DO
            END DO
            SIGMA_ACTIVE_MAT(1,1) = (3.0D0/(4.0D0*PI))*SIGMA_ACTIVE_MAT(1,1)
            SIGMA_ACTIVE_MAT(2,2) = (3.0D0/(4.0D0*PI))*SIGMA_ACTIVE_MAT(2,2)
            SIGMA_ACTIVE_MAT(3,3) = (3.0D0/(4.0D0*PI))*SIGMA_ACTIVE_MAT(3,3)
            SIGMA_ACTIVE_MAT(1,2) = (3.0D0/(4.0D0*PI))*SIGMA_ACTIVE_MAT(1,2)
            SIGMA_ACTIVE_MAT(1,3) = (3.0D0/(4.0D0*PI))*SIGMA_ACTIVE_MAT(1,3)
            SIGMA_ACTIVE_MAT(2,3) = (3.0D0/(4.0D0*PI))*SIGMA_ACTIVE_MAT(2,3)
            SIGMA_ACTIVE_MAT(2,1) = SIGMA_ACTIVE_MAT(1,2)
            SIGMA_ACTIVE_MAT(3,1) = SIGMA_ACTIVE_MAT(1,3)
            SIGMA_ACTIVE_MAT(3,2) = SIGMA_ACTIVE_MAT(2,3)
      END IF

      RETURN
      END SUBROUTINE

C     DEFITION OF EVALUATUION OF STRESS FIBRE CONCENTRATION IN EACH FIBRE SUBROUTINE
      SUBROUTINE EVAL_ETA_SF_THETA(K_F,K_B,SIG_MAX,N_FIBRES,DTIME,
     1           NSTATV,STATEV,SIGMA_ACTIVE_THETA_OLD, ETA_SF_THETA)

      IMPLICIT NONE
      INTEGER N_FIBRES,I,J, NSTATV
      REAL(KIND=8) ATHETA, CAL_CONC, SIGMA_BAR_THETA, SIG_MAX, DTIME
      REAL(KIND=8) ETA_SF_THETA(N_FIBRES),K_F, K_F_BAR, K_B, K_B_BAR
      REAL(KIND=8) SIGMA_ACTIVE_THETA_OLD(N_FIBRES)
      REAL(KIND=8) STATEV(NSTATV)

C     STRESS FIBRE CONCENTRATION IN EACH FIBRE OBTAINED BY FORWARD EULER INTEGRATION
      ATHETA = 1000.0D0
      K_F_BAR = K_F/ATHETA
      K_B_BAR = K_B/ATHETA

      DO I = 1, N_FIBRES
        ETA_SF_THETA(I) = STATEV(I)
            CAL_CONC = 0.1D0
            IF ( ETA_SF_THETA(I) .EQ. 0.0) THEN
                  SIGMA_BAR_THETA = 0.0D0
            ELSE      
                  SIGMA_BAR_THETA = SIGMA_ACTIVE_THETA_OLD(I)/
     1                        (ETA_SF_THETA(I)*SIG_MAX)
                  IF (SIGMA_BAR_THETA .gt. 1.0) THEN
                        SIGMA_BAR_THETA = 1.0D0
                  END IF
        END IF

        ETA_SF_THETA(I) = ETA_SF_THETA(I) + DTIME*( 
     1              (1.0D0 - ETA_SF_THETA(I))*CAL_CONC*K_F_BAR -
     2              (1.0D0 - SIGMA_BAR_THETA)*ETA_SF_THETA(I)*K_B_BAR)

      END DO

      RETURN
      END SUBROUTINE EVAL_ETA_SF_THETA

C     DEFITION OF EVALUATUION OF ACTIVE STRESS IN EACH FIBRE SUBROUTINE
      SUBROUTINE EVAL_ACTIVE_STRESS_THETA(K_T,SIG_MAX,N_FIBRES,STATEV,
     1           NSTATV,EPSDOT_THETA,SIGMA_ACTIVE_THETA, PHI_VEC)

      IMPLICIT NONE
      INTEGER N_FIBRES,I, NSTATV
      REAL(KIND=8) ETA_SF_THETA_I, SIG_MAX, K_T
      REAL(KIND=8) SIGMA_ACTIVE_THETA(N_FIBRES),EPSDOT_THETA(N_FIBRES)
      REAL(KIND=8) STATEV(NSTATV), PHI_VEC(N_FIBRES)

C     EVALUATE ACTIVE STRESS IN EACH FIBRE
      DO I = 1,N_FIBRES
        ETA_SF_THETA_I = STATEV(I)
        SIGMA_ACTIVE_THETA(I) = ETA_SF_THETA_I*SIG_MAX*( 1.0D0 + 
     1                          (K_T*EPSDOT_THETA(I)/
     2           SQRT(1.0D0 + EPSDOT_THETA(I)*EPSDOT_THETA(I))) )
      END DO

      RETURN
      END SUBROUTINE EVAL_ACTIVE_STRESS_THETA

c     INCLUDE AMPLITUDE SUBROUTINE
c      INCLUDE "UAMP.f"