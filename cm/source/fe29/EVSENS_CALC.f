      SUBROUTINE EVSENS_CALC(ISTATE,IUSER,IIY,IY,NEELEM,NMNO,NPNODE,nr,
     '  NVJP,NYNP,PAOPTI,PE,PF,PMIN,PMAX,R,RESID,RESJAC,SENSIT,USER,XC,
     '  XC2,XP,YP,TYPE,CONSTR,DEFORMED,FIX,FIXP,OPFILE,OUTPTPRN,ERROR,*)

C#### Subroutine: EVSENS
C###  Description:
C###    EVSENS calculates sensitivity of optimisation parameters
C###    wrt geometry pressure bcs or force bcs.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER ISTATE(*),IUSER(*),IIY,IY,NEELEM(0:NE_R_M,0:NRM),
     '  NMNO(1:2,0:NOPM,NXM),NPNODE(0:NP_R_M,0:NRM),nr,NVJP(NJM,NPM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 PAOPTI(*),PE(2,NEM),PF(2,NEM),PMIN(*),PMAX(*),R(NOPM,*),
     '  RESID(*),RESJAC(NREM,*),SENSIT(NOPM,NJM,*),USER(*),XC(*),XC2(*),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER TYPE*(*),ERROR*(*)
      LOGICAL CONSTR,DEFORMED,FIX(NYM,NIYFIXM,NXM),FIXP(2,NEM),OPFILE,
     '  OUTPTPRN
!     Local Variables
      INTEGER iface,ne,noelem,nj,nk,noopti,nonode,np,nv,nx,ny
      REAL*8 DELTA_POS
      DATA DELTA_POS/1.0D-4/


      CALL ENTERS('EVSENS_CALC',*9999)

      nx=1 !Temporary cpb 22/11/94

      IF(TYPE(1:5).EQ.'NODAL') THEN !sensitivity wrt nodal parameters
        !Loop over all global nodes and nj directions (nc=1,nk=1)
        nk=1  !only perturb nodal values (not derivatives)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            DO nv=1,NVJP(nj,np)
              ny=NYNP(nk,nv,nj,np,0,1,nr)
              IF(.NOT.CONSTR.OR.FIX(ny,IIY,nx)) THEN
                IF(DEFORMED) THEN
                  !Perturb current nodal,nj'th value in YP
                  YP(ny,IY,nx)=YP(ny,IY,nx)+DELTA_POS
                ELSE
                  !Perturb undeformed current nodal,nj'th coord in XP
                  XP(nk,nv,nj,np)=XP(nk,nv,nj,np)+DELTA_POS
                ENDIF

                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$               call mp_setlock()
                  WRITE(OP_STRING,'('' Optimising with perturbed '
     '              //'value at node '',I3,'' nj='',I1)') np,nj
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$               call mp_unsetlock()
                ENDIF

C ***           Optimise mat params for perturbed nodal params
                CALL NAGMINA(ISTATE,PAOPTI,PMIN,PMAX,R,RESID,RESJAC,
     '            XC2,IUSER,USER,OUTPTPRN,ERROR,*9999)

                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$               call mp_setlock()
                  WRITE(OP_STRING,'('' Material param # '','
     '              //'5(I2,10X)/,(20X,5(I2,10X)))')
     '              (NMNO(1,noopti,nx),noopti=1,NTOPTI)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' Solution    '',5D12.4/,'
     '              //'(13X,5D12.4))') (XC2(noopti),noopti=1,NTOPTI)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$               call mp_unsetlock()
                ENDIF

                IF(DEFORMED) THEN
                  !Reset current nodal,nj'th value in YP
                  YP(ny,IY,nx)=YP(ny,IY,nx)-DELTA_POS
                ELSE
                  !Reset undeformed current nodal,nj'th value in XP
                  XP(nk,nv,nj,np)=XP(nk,nv,nj,np)-DELTA_POS
                ENDIF

C ***           Calc parameter sensitivity wrt current nodal param
                DO noopti=1,NTOPTI
                  SENSIT(noopti,nj,np)=(XC2(noopti)-XC(noopti))
     '                                 /DELTA_POS
                ENDDO
              ENDIF
            ENDDO !nv
          ENDDO !nj
        ENDDO !nonode (np)

        !Output parameter sensitivities
        IF(OPFILE) IOFI=IOFILE1 !necessary since EVRESI changes IOFI
        IF(DEFORMED) THEN
          WRITE(OP_STRING,'(/'' Sensitivities wrt nodal values '
     '      //'in YP(IY='',I1,''):'')') IY
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(/'' Sensitivities wrt undeformed '
     '      //'nodal values in XP:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'('' Mat. param. number'',5(I2,10X)/,'
     '    //'(20X,5(I2,10X)))') (NMNO(1,noopti,nx),noopti=1,NTOPTI)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        nk=1 !only perturb force values (not derivatives)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            DO nv=1,NVJP(nj,np)
              ny=NYNP(nk,nv,nj,np,0,1,nr)
              IF(.NOT.CONSTR.OR.FIX(ny,IIY,nx)) THEN
                WRITE(OP_STRING,'('' Node '',I3,'' nj='',I1,'
     '            //'5D12.4/,(14X,5D12.4))')
     '            np,nj,(SENSIT(noopti,nj,np),noopti=1,NTOPTI)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !nv
          ENDDO !nj
        ENDDO !nonode (np)

      ELSE IF(TYPE(1:8).EQ.'PRESSURE') THEN !sensit. wrt pressure bcs
        !Loop over elements
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          !Loop over pressure bcs in the element
          DO iface=1,2
            DO noopti=1,NTOPTI
              SENSIT(noopti,iface,ne)=-9999 !to distinguish from
            ENDDO                           !non-bc in output
            IF(FIXP(iface,ne)) THEN
              !Perturb current pressure bc
              PE(iface,ne)=PE(iface,ne)+DELTA_POS
              PF(iface,ne)=PF(iface,ne)+DELTA_POS

              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$             call mp_setlock()
                WRITE(OP_STRING,'('' Optimising with pressure bc '','
     '            //'I1,'' in element '',I3,'' perturbed'')') iface,ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$             call mp_unsetlock()
              ENDIF

C ***         Optimise material parameters for perturbed pressure bc
              CALL NAGMINA(ISTATE,PAOPTI,PMIN,PMAX,R,RESID,RESJAC,
     '          XC2,IUSER,USER,OUTPTPRN,ERROR,*9999)

              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$             call mp_setlock()
                WRITE(OP_STRING,'('' Material param # '',5(I2,10X)/,'
     '            //'(20X,5(I2,10X)))')
     '            (NMNO(1,noopti,nx),noopti=1,NTOPTI)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Solution    '',5D12.4/,'
     '            //'(13X,5D12.4))') (XC2(noopti),noopti=1,NTOPTI)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$             call mp_unsetlock()
              ENDIF

              !Reset current pressure bc
              PE(iface,ne)=PE(iface,ne)-DELTA_POS
              PF(iface,ne)=PF(iface,ne)-DELTA_POS

C ***         Calculate parameter sensitivity wrt current pressure bc
              DO noopti=1,NTOPTI
                SENSIT(noopti,iface,ne)=(XC2(noopti)-XC(noopti))
     '                                         /DELTA_POS
              ENDDO

            ENDIF
          ENDDO
        ENDDO

        !Output parameter sensitivities
        IF(OPFILE) IOFI=IOFILE1 !necessary since EVRESI changes IOFI
        WRITE(OP_STRING,'(/'' Sensitivities wrt pressure bcs:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Mat param number'',8X,5(I2,10X)/,'
     '    //'(25X,5(I2,10X)))') (NMNO(1,noopti,nx),noopti=1,NTOPTI)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO iface=1,2
            IF(FIXP(iface,ne)) THEN
              WRITE(OP_STRING,'('' Element '',I3,'', Face '',I1,'
     '          //'5D12.4/,(20X,5D12.4))')
     '          ne,iface,(SENSIT(noopti,iface,ne),noopti=1,NTOPTI)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('EVSENS_CALC')
      RETURN
 9999 CALL ERRORS('EVSENS_CALC',ERROR)
      CALL EXITS('EVSENS_CALC')
      RETURN 1
      END


