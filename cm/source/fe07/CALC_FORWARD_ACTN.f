      REAL*8 FUNCTION CALC_FORWARD_ACTN(DERIVATIVE_TYPE,DERIV_PARAM_NUM,
     '  HEART_ny,i,j,NHP,NKH,NPLIST3,NVHP,nx,NYNP,ACTN_TIME_TOL,
     '  T_BH,TIME,YP,INTERPOLATE,ERR,ERROR)

C#### Function: CALC_FORWARD_ACTN
C###  Description:
C###    CALC_FORWARD_ACTN calculates and returns the body surface potential,
C###    PHI, at a particular body surface location, j, at a specific point in
C###    time, TIME, from a supplied activation field using the
C###    activation waveform that has been defined by the user with
C###    define transfer. DERIVATIVE_TYPE specifies the type of
C###    derivative to return in PHI. If DERIVATIVE_TYPE is 0|1|2|3 then
C###    PHI contains PHI_j|dPHI_j/dt|dPhi_j/dtau_i|dPhi_j/du_k where
C###    i is supplied as is HEART_ny, the mesh dof number corresponding
C###    to i. k is supplied in DERIV_PARAM_NUM. If INTERPOLATE is .TRUE.
C###    then body surface potential is calculated from an interpolated
C###    activation field otherwise a discrete activation field is used.
C###    The parameter ACTN_TIME_TOL is used to provide a tolerance for
C###    on/off type waveform types e.g. Heaviside step.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER DERIVATIVE_TYPE,DERIV_PARAM_NUM,ERR,HEART_ny,i,j,
     '  NHP(NPM,0:NRM),NKH(NHM,NPM,NCM,0:NRM),NPLIST3(0:NPM),
     '  NVHP(NHM,NPM,NCM,0:NRM),nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 ACTN_TIME_TOL,T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  TIME,YP(NYM,NIYM)
      LOGICAL INTERPOLATE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER icol,nh,nhx,nk,nolist,np,nv,ny
      REAL*8 HALF_WINDOW,PHI,s

      PHI=0.0d0
      CALC_FORWARD_ACTN=0.0d0 !Initialise so it does return a value on error
      IF(INTERPOLATE) THEN !Use interpolated activation field
        ERR=1
        ERROR='>>Not implemeted yet'
        RETURN

C       LKC 7-JUL-2003 Commenting this out as FTNCHK doesn't like having
C       code that it can't reach and the code is not going to get done
C       in the distant future.
C        
C        IF(DERIVATIVE_TYPE.EQ.0) THEN !Calculate Phi_j
CC         Note the resting potential does not enter the equation
CC         as it is assumed the the row sums of T_BH are zero. At the
CC         moment the transmembrane jump is assumed constant as is thus
CC         included as a multiplicative factor after the sum calculations
C          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !Arctan function
C          ELSE
C            ERR=1
C            ERROR='>>Invalid activation waveform type'
C            RETURN
C          ENDIF
C        ELSE IF(DERIVATIVE_TYPE.EQ.1) THEN !Calculate dPhi
C          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !Arctan function
C          ELSE
C            ERR=1
C            ERROR='>>Invalid activation waveform type'
C            RETURN
C          ENDIF
C        ELSE IF(DERIVATIVE_TYPE.EQ.2) THEN !Calculate
C          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !Arctan function
C          ELSE
C            ERR=1
C            ERROR='>>Invalid activation waveform type'
C            RETURN
C          ENDIF
C        ELSE IF(DERIVATIVE_TYPE.EQ.3) THEN !Calculate
C          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
C          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !Arctan function
C          ELSE
C            ERR=1
C            ERROR='>>Invalid activation waveform type'
C            RETURN
C          ENDIF
C        ELSE
C          ERR=1
C          ERROR='>>Invalid derivative type'
C          RETURN
C        ENDIF
        
      ELSE !Use discrete activation field
C       At the moment assume direct 1-1 mapping between heart nodes
C       and row elements in the transfer matrix and thus loop over
C       to ISIZE_TBH. In the future will have to do this properly.
        IF(DERIVATIVE_TYPE.EQ.0) THEN !Calculate Phi_j
C         Note the resting potential does not enter the equation
C         as it is assumed the the row sums of T_BH are zero. At the
C         moment the transmembrane jump is assumed constant as is thus
C         included as a multiplicative factor after the sum calculations
          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
            icol=0
            DO nolist=1,NPLIST3(0) !list of heart nodes
              np=NPLIST3(nolist)
              DO nhx=1,NHP(np,TRSF_NR_FIRST)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                  DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '              KTYP93(1,TRSF_NR_FIRST),1)
                    ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                    icol=icol+1 !Appropriate column of T_BH
                    s=TIME-YP(ny,1)
                    IF(s.GT.-ACTN_TIME_TOL) THEN
                      PHI=PHI+T_BH(j,icol)
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhx
            ENDDO !nolist
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
            HALF_WINDOW=TRSF_ACTN_WAVE_WIDTH/2.0d0
            icol=0
            DO nolist=1,NPLIST3(0) !list of heart nodes
              np=NPLIST3(nolist)
              DO nhx=1,NHP(np,TRSF_NR_FIRST)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                  DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '              KTYP93(1,TRSF_NR_FIRST),1)
                    ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                    icol=icol+1 !Appropriate column of T_BH
                    s=TIME-YP(ny,1)
                    IF(s.GE.-HALF_WINDOW.AND.s.LE.0.0d0) THEN
                      PHI=PHI+T_BH(j,icol)*
     '                  0.5d0*(s/HALF_WINDOW+1.0d0)**2
                    ELSE IF(s.GT.0.0d0.AND.s.LE.HALF_WINDOW) THEN
                      PHI=PHI+T_BH(j,icol)*
     '                  (1.0d0-0.50d0*(s/HALF_WINDOW-1.0d0)**2)
                    ELSE IF(s.GT.HALF_WINDOW) THEN
                      PHI=PHI+T_BH(j,icol)
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhx
            ENDDO !nolist
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !arctan function
            icol=0
            DO nolist=1,NPLIST3(0) !list of heart nodes
              np=NPLIST3(nolist)
              DO nhx=1,NHP(np,TRSF_NR_FIRST)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                  DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '              KTYP93(1,TRSF_NR_FIRST),1)
                    ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                    icol=icol+1 !Appropriate column of T_BH
                    s=TIME-YP(ny,1)
                    PHI=PHI+T_BH(j,icol)*
     '                0.5d0*(1.0d0+2.0d0/PI*
     '                ATAN(PI*s/TRSF_ACTN_WAVE_WIDTH))
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhx
            ENDDO !nolist
          ELSE
            ERR=1
            ERROR='>>Invalid activation waveform type'
            RETURN
          ENDIF
          PHI=PHI*TRSF_ACTN_WAVE_JUMP
        ELSE IF(DERIVATIVE_TYPE.EQ.1) THEN !Calculate dPhi_j/dt
          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
            icol=0
            DO nolist=1,NPLIST3(0) !list of heart nodes
              np=NPLIST3(nolist)
              DO nhx=1,NHP(np,TRSF_NR_FIRST)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                  DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '              KTYP93(1,TRSF_NR_FIRST),1)
                    ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                    icol=icol+1 !Appropriate column of T_BH
                    s=TIME-YP(ny,1)
                    IF(DABS(s).LE.ACTN_TIME_TOL) THEN
                      PHI=PHI+T_BH(j,icol)
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhx
            ENDDO !nolist
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
            HALF_WINDOW=TRSF_ACTN_WAVE_WIDTH/2.0d0
            icol=0
            DO nolist=1,NPLIST3(0) !list of heart nodes
              np=NPLIST3(nolist)
              DO nhx=1,NHP(np,TRSF_NR_FIRST)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                  DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '              KTYP93(1,TRSF_NR_FIRST),1)
                    ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                    icol=icol+1 !Appropriate column of T_BH
                    s=TIME-YP(ny,1)
                    IF(s.GE.-HALF_WINDOW.AND.s.LE.0.0d0) THEN
                      PHI=PHI-T_BH(j,icol)*
     '                  (s/HALF_WINDOW+1.0d0)/HALF_WINDOW
                    ELSE IF(s.GT.0.0d0.AND.s.LE.HALF_WINDOW) THEN
                      PHI=PHI+T_BH(j,icol)*
     '                  (s/HALF_WINDOW-1.0d0)/HALF_WINDOW
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhx
            ENDDO !nonolist
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !arctan function
            icol=0
            DO nolist=1,NPLIST3(0) !list of heart nodes
              np=NPLIST3(nolist)
              DO nhx=1,NHP(np,TRSF_NR_FIRST)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                  DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '              KTYP93(1,TRSF_NR_FIRST),1)
                    ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                    icol=icol+1 !Appropriate column of T_BH
                    s=TIME-YP(ny,1)
                    PHI=PHI+T_BH(j,icol)*TRSF_ACTN_WAVE_WIDTH/
     '                (TRSF_ACTN_WAVE_WIDTH**2+PI**2*s**2)
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhx
            ENDDO !nolist
          ELSE
            ERR=1
            ERROR='>>Invalid activation waveform type'
            RETURN
          ENDIF
          PHI=PHI*TRSF_ACTN_WAVE_JUMP
        ELSE IF(DERIVATIVE_TYPE.EQ.2) THEN !Calculate dPhi_j/dtau_i
          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
            s=TIME-YP(HEART_ny,1)
            IF(DABS(s).LE.ACTN_TIME_TOL) THEN
              PHI=-TRSF_ACTN_WAVE_JUMP*T_BH(j,i)
            ELSE
              PHI=0.0d0
            ENDIF
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
            HALF_WINDOW=TRSF_ACTN_WAVE_WIDTH/2.0d0
            s=TIME-YP(HEART_ny,1)
            IF(s.GE.-HALF_WINDOW.AND.s.LE.0.0d0) THEN
              PHI=-TRSF_ACTN_WAVE_JUMP*T_BH(j,i)*(s/HALF_WINDOW+
     '          1.0d0)/HALF_WINDOW
            ELSE IF(s.GT.0.0d0.AND.s.LE.HALF_WINDOW) THEN
              PHI=TRSF_ACTN_WAVE_JUMP*T_BH(j,i)*(s/HALF_WINDOW-
     '          1.0d0)/HALF_WINDOW
            ELSE
              PHI=0.0d0
            ENDIF
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !arctan function
            s=TIME-YP(HEART_ny,1)
            PHI=-TRSF_ACTN_WAVE_JUMP*T_BH(j,i)*TRSF_ACTN_WAVE_WIDTH/
     '        (TRSF_ACTN_WAVE_WIDTH**2+PI**2*s**2)
          ELSE
            ERR=1
            ERROR='>>Invalid activation waveform type'
            RETURN
          ENDIF
        ELSE IF(DERIVATIVE_TYPE.EQ.3) THEN !Calculate dPhi_j/du_k
          IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function
            IF(DERIV_PARAM_NUM.EQ.1) THEN !Resting potential
              PHI=0.0d0
            ELSE IF(DERIV_PARAM_NUM.EQ.2) THEN !Transmembrane jump
              PHI=0.0d0
              icol=0
              DO nolist=1,NPLIST3(0) !list of heart nodes
                np=NPLIST3(nolist)
                DO nhx=1,NHP(np,TRSF_NR_FIRST)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                    DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '                KTYP93(1,TRSF_NR_FIRST),1)
                      ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                      icol=icol+1 !Appropriate column of T_BH
                      s=TIME-YP(ny,1)
                      IF(s.GT.-ACTN_TIME_TOL) PHI=PHI+T_BH(j,icol)
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nhx
              ENDDO !nolist
            ELSE
              ERR=1
              ERROR='>>Invalid derivative parameter number'
              RETURN
            ENDIF
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function
            IF(DERIV_PARAM_NUM.EQ.1) THEN !Resting potential
              PHI=0.0d0
            ELSE IF(DERIV_PARAM_NUM.EQ.2) THEN !Transmembrane jump
              PHI=0.0d0
              HALF_WINDOW=TRSF_ACTN_WAVE_WIDTH/2.0d0
              icol=0
              DO nolist=1,NPLIST3(0) !list of heart nodes
                np=NPLIST3(nolist)
                DO nhx=1,NHP(np,TRSF_NR_FIRST)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                    DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '                KTYP93(1,TRSF_NR_FIRST),1)
                      ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                      icol=icol+1 !Appropriate column of T_BH
                      s=TIME-YP(ny,1)
                      IF(s.GE.-HALF_WINDOW.AND.s.LE.0.0d0) THEN
                        PHI=PHI+T_BH(j,icol)*
     '                    0.5d0*(s/HALF_WINDOW+1.0d0)**2
                      ELSE IF(s.GT.0.0d0.AND.s.LE.HALF_WINDOW) THEN
                        PHI=PHI+T_BH(j,icol)*
     '                    (1.0d0-0.50d0*(s/HALF_WINDOW-1.0d0)**2)
                      ELSE IF(s.GT.HALF_WINDOW) THEN
                        PHI=PHI+T_BH(j,icol)
                      ENDIF
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nhx
              ENDDO !nolist
            ELSE IF(DERIV_PARAM_NUM.EQ.3) THEN !Window width
              ERR=1
              ERROR='>>Not implemented yet'
              RETURN
            ELSE
              ERR=1
              ERROR='>>Invalid derivative parameter number'
              RETURN
            ENDIF
          ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !arctan function
            IF(DERIV_PARAM_NUM.EQ.1) THEN !Resting potential
              PHI=0.0d0
            ELSE IF(DERIV_PARAM_NUM.EQ.2) THEN !Transmembrane jump
              icol=0
              DO nolist=1,NPLIST3(0) !list of heart nodes
                np=NPLIST3(nolist)
                DO nhx=1,NHP(np,TRSF_NR_FIRST)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                    DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '                KTYP93(1,TRSF_NR_FIRST),1)
                      ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                      icol=icol+1 !Appropriate column of T_BH
                      s=TIME-YP(ny,1)
                      PHI=PHI+T_BH(j,icol)*0.5d0*(1.0d0+2.0d0/PI*
     '                  ATAN(PI*s/TRSF_ACTN_WAVE_WIDTH))
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nhx
              ENDDO !nolist
            ELSE IF(DERIV_PARAM_NUM.EQ.3) THEN !Window width
              icol=0
              DO nolist=1,NPLIST3(0) !list of heart nodes
                np=NPLIST3(nolist)
                DO nhx=1,NHP(np,TRSF_NR_FIRST)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,TRSF_NR_FIRST)
                    DO nk=1,MAX(NKH(nh,np,1,TRSF_NR_FIRST)-
     '                KTYP93(1,TRSF_NR_FIRST),1)
                      ny=NYNP(nk,nv,nh,np,0,1,TRSF_NR_FIRST)
                      icol=icol+1 !Appropriate column of T_BH
                      s=TIME-YP(ny,1)
                      PHI=PHI-TRSF_ACTN_WAVE_JUMP*T_BH(j,icol)*s/
     '                  (TRSF_ACTN_WAVE_WIDTH**2+PI**2*s**2)
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nhx
              ENDDO !nolist
            ELSE
              ERR=1
              ERROR='>>Invalid derivative parameter number'
              RETURN
            ENDIF
          ELSE
            ERR=1
            ERROR='>>Invalid activation waveform type'
            RETURN
          ENDIF
        ELSE
          ERR=1
          ERROR='>>Invalid derivative type'
          RETURN
        ENDIF
      ENDIF

      CALC_FORWARD_ACTN=PHI

      RETURN
      END

