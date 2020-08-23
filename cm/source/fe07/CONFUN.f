      SUBROUTINE CONFUN(IBT,IDO,NBJ,NEEDCON,NEELEM,NEL,NKH,
     '  NLL,NLNO,NNL,NONL,NONY,NPL,NPNE,NPNODE,NPNY,nr,NVHP,
     '  NVJL,NVJP,NYNO,NYNP,PAOPTY,
     '  CM,CJACM,CONTR,CONJAC,DL,PAOPTI,SE,XP,PARAMTYPE,ERROR,*)

C#### Subroutine: CONFUN
C###  Description:
C###    CONFUN returns the constraints for optimising geometric
C###    and material parameters etc.

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),
     '  NEEDCON(*),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NKH(NHM,NPM,NCM,0:NRM),NLL(12,NEM),NLNO(NOPM),NNL(0:4,12,NBFM),
     '  NONL(NLM),NONY(0:NOYM,NYM,NRCM),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJL(4,NJM,NLM),NVJP(NJM,NPM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),PAOPTY(NOPM)
      REAL*8 CM(*),CJACM(NCOM,*),CONTR(*),CONJAC(NCOM,*),DL(3,NLM),
     '  PAOPTI(*),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER nb,nj,njj2,nk,NKTOT,nl,no,no_aero,nocont,
     '  nonode,noopti,noy,np,nv,NVTOT,ny,nyo
      REAL*8 PERIM,TEMP
      LOGICAL INOPTI

      CALL ENTERS('CONFUN',*9999)

      IF(PARAMTYPE(1:13).EQ.'IPOPTI_params') THEN
        IF(KTYP27.EQ.7) THEN !Aero wake press diff & sail stress
          CONTR(1)=TE_VELOC_DIFF !trailing edge velocity difference

        ELSE IF(KTYP27.EQ.8) THEN !Aerofoil lift & wake press. diff.
c          IF(NEEDCON(1).GT.0) THEN !constraint is evaluated
            PERIM=0.0d0
            DO no_aero=1,NL_AERO(0,1)+NL_AERO(0,2)
              IF(no_aero.LE.NL_AERO(0,1)) THEN !upper surface
                nl=NL_AERO(no_aero,1)
              ELSE                             !lower surface
                nl=NL_AERO(no_aero-NL_AERO(0,1),2)
              ENDIF
!           Calculate the current estimate of DL(3,nl)
              CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,
     '          NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,XP,ERROR,*9999)
              PERIM=PERIM+DL(3,nl)
            ENDDO
            CONTR(1)=PERIM
c         ENDIF
        ENDIF
!       IF(DOP) THEN
          WRITE(OP_STRING,'('' Constraint CONTR(1) = '',D12.4)')
     '      CONTR(1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!       ENDIF

      ELSE IF(PARAMTYPE(1:12).EQ.'DATA_FITTING') THEN
!       Copy in the current estimates of XP and DL from PAOPTI
        DO noopti=1,NTOPTI
          IF(PAOPTY(noopti).EQ.1) THEN !Parameter is geometric dof
            DO nyo=1,NYNO(0,noopti,2)
              ny=NYNO(nyo,noopti,2)
              nk=NPNY(1,ny,0)
              nv=NPNY(2,ny,0)
              nj=NPNY(3,ny,0)
              np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              XP(nk,nv,nj,np)=PAOPTI(noopti)
            ENDDO
          ELSE IF(PAOPTY(noopti).EQ.2) THEN !Parameter is a line
            nl=NLNO(noopti)
            DL(1,nl)=PAOPTI(noopti)
            DL(2,nl)=DL(1,nl)
            IF(JTYP2B.EQ.1.AND.NPL(4,0,nl).GT.0) THEN
              DL(1,NPL(4,0,nl))=-1.0d0*DL(2,nl)
              DL(2,NPL(4,0,nl))=-1.0d0*DL(1,nl)
              DL(3,NPL(4,0,nl))=DL(3,nl)
            ENDIF
          ELSE
            ERROR=' >>Invalid PAOPTY for data fitting'
            GOTO 9999
          ENDIF
        ENDDO
        IF(G1SCALING) THEN
!         Calculate the current estimate of DL(3,nl)
          DO nl=1,NLT
            CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,
     '        NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,XP,ERROR,*9999)
          ENDDO
          DO nb=1,NBFT
            CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,DL,SE,
     '        ERROR,*9999)
          ENDDO
        ENDIF
!       Generate the normality constraints on the arc length derivatives
        nocont=0
        IF(KTYP1B.EQ.1) THEN !Derivative constraints
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
C           Find number of versions and derivatives at the node
            nj=NJ_LOC(NJL_GEOM,1,nr)
            NKTOT=NKH(nj,np,1,nr)
            NVTOT=NVHP(nj,np,1,nr)
            DO njj2=2,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj2,nr)
              IF(NKH(nj,np,1,nr).NE.NKTOT.AND.NKTOT.NE.0) THEN
                WRITE(OP_STRING,
     '            '('' Differing numbers of derivatives at node '''
     '            //',I3)') np
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                NKTOT=0
              ENDIF
              IF(NVHP(nj,np,1,nr).NE.NVTOT.AND.NVTOT.NE.0) THEN
                WRITE(OP_STRING,
     '            '('' Differing numbers of versions at node '''
     '            //',I3)') np
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                NVTOT=0
              ENDIF
            ENDDO
            DO nk=2,NKTOT
              IF(nk.EQ.2.OR.nk.EQ.3.OR.nk.EQ.5) THEN !first deriv
                INOPTI=.FALSE.
                DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj2,nr)
                  DO nv=1,NVHP(nj,np,1,nr)
                    ny=NYNP(nk,nv,nj,np,0,1,nr)
                    IF(ny.GT.0) THEN
                      no=NONY(1,ny,2)
                      IF(no.GT.0) THEN
C                       Only need one constraint if mesh dofs map to
C                       one solution dof.
                        IF(NYNO(1,no,2).EQ.ny) INOPTI=.TRUE.
                      ENDIF !no.GT.0
                    ENDIF !ny.GT.0
                  ENDDO !nv
                ENDDO !nj
                IF(INOPTI) THEN
                  nocont=nocont+1
                  IF(KTYP29.EQ.1) THEN
                    IF(NEEDCON(nocont).GT.0) THEN
                      TEMP=0.0d0
                      DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                        nj=NJ_LOC(NJL_GEOM,njj2,nr)
                        TEMP=TEMP+XP(nk,1,nj,np)**2
                      ENDDO !nj
                      CONTR(nocont)=DSQRT(TEMP)
                      IF(TEMP.EQ.0.0d0) THEN
                        ERROR='>> Node normality constraint is zero'
                        GOTO 9999
                      ENDIF
                      DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                        nj=NJ_LOC(NJL_GEOM,njj2,nr)
                        ny=NYNP(nk,1,nj,np,0,1,nr)
                        DO noy=1,NONY(0,ny,2)
                          noopti=NONY(noy,ny,2)
                          CONJAC(nocont,noopti)=XP(nk,1,nj,np)/
     '                      CONTR(nocont)
                        ENDDO !noy
                      ENDDO !nj
                    ENDIF
                  ELSE IF(KTYP29.EQ.2) THEN
                    CM(nocont)=0.0d0
                    DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj2,nr)
                      CM(nocont)=CM(nocont)+XP(nk,1,nj,np)**2
                    ENDDO
                    CM(nocont)=DSQRT(CM(nocont))
C CPB 2/6/94 Should put in code to test if the jacobian needs to be
C calculated
                    IF(CM(nocont).EQ.0.0d0) THEN
                      ERROR='>> Node normality constraint is zero'
                      GOTO 9999
                    ENDIF
                    DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj2,nr)
                      ny=NYNP(nk,1,nj,np,0,1,nr)
                      DO noy=1,NONY(0,ny,2)
                        noopti=NONY(noy,ny,2)
                        CJACM(nocont,noopti)=XP(nk,1,nj,np)/CM(nocont)
                      ENDDO !noy
                    ENDDO !nj
                  ENDIF !ktyp29
                ENDIF !inopti
              ENDIF !nk
            ENDDO !nk
          ENDDO !npnode (np)
        ENDIF !Derivative constraints
        IF(G1SCALING) THEN
!         Generate the line length constraints
          DO nl=1,NLT
            IF(NONL(nl).GT.0) THEN !The line is an optimising variable
              nocont=nocont+1
              IF(KTYP29.EQ.1) THEN
                IF(NEEDCON(nocont).GT.0) THEN
                  CONTR(nocont)=DL(3,nl)-DL(1,nl)
                   CALL D_CONFUN(nocont,nl,IDO,NBJ,NEL(0,nl),
     '              NONL,NONY,NPL(1,0,nl),NPNE,nr,NVJP,NYNP,
     '              CJACM,CONJAC,DL(1,nl),XP,PARAMTYPE,ERROR,*9999)
                ENDIF
              ELSE IF(KTYP29.EQ.2) THEN
                CM(nocont)=DL(3,nl)-DL(1,nl)
                CALL D_CONFUN(nocont,nl,IDO,NBJ,NEL(0,nl),
     '            NONL,NONY,NPL(1,0,nl),NPNE,nr,NVJP,NYNP,
     '            CJACM,CONJAC,DL(1,nl),XP,PARAMTYPE,ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(CONFUN_1)
          WRITE(OP_STRING,'('' Constraint values:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO nocont=1,NTCNTR
            IF(KTYP29.EQ.1) THEN
              WRITE(OP_STRING,'('' CONTR('',I3,'') = '',D12.4)')
     '          nocont,CONTR(nocont)
            ELSE IF(KTYP29.EQ.2) THEN
              WRITE(OP_STRING,'('' CONTR('',I3,'') = '',D12.4)')
     '          nocont,CM(nocont)
            ENDIF
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'('' Constraint Jacobian values:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO nocont=1,NTCNTR
            IF(KTYP29.EQ.1) THEN
              WRITE(OP_STRING,'('' CONJAC('',I3,'',no=1,ntopti) ='','
     '          //'8(1X,D12.6)/,:(26X,8(1X,D12.6)))')
     '          nocont,(CONJAC(nocont,noopti),noopti=1,NTOPTI)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ELSE IF(KTYP29.EQ.2) THEN
              WRITE(OP_STRING,'('' CJACM('',I3,'',no=1,ntopti) = '','
     '          //'8(1X,D12.6)/,:(25X,8(1X,D12.6)))')
     '          nocont,(CJACM(nocont,noopti),noopti=1,NTOPTI)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
CC$OMP     END CRITICAL(CONFUN_1)
        ENDIF
      ENDIF

      CALL EXITS('CONFUN')
      RETURN
 9999 CALL ERRORS('CONFUN',ERROR)
      CALL EXITS('CONFUN')
      RETURN 1
      END


