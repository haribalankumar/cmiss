      SUBROUTINE GET_TNVECTOR(IBT,IDO,INP,NBJ,NENP,NKJE,np,NPF,
     '  NP_INTERFACE,NPNE,nr,NRE,NVJE,nx,NORMAL,SE,TANGENT,XA,XE,
     '  XP,ERROR,*)

C#### Subroutine: GET_TNVECTOR
C###  Description:
C###    GET_TNVECTOR returns the tangent(s) and normal vectors to
C###    a mesh at a node np. It is only implemented for rectangular
C###    cartesian coordinates at the moment.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NENP(NPM,0:NEPM),NKJE(NKM,NNM,NJM,NEM),np,
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),nx
      REAL*8 SE(NSM,NBFM,NEM),NORMAL(3),TANGENT(3,2),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ni,nii,nj,njj,nn,noelem,nonp,nu,NU1(0:3)
      REAL*8 PXI,SUM,TEMP(3),XI(2)
      LOGICAL COLLAPSEDNODE,FEMFEMCOUP,FOUND,HERMSECTOR,ISATCOLLAPSE

      DATA NU1/1,2,4,7/

      CALL ENTERS('GET_TNVECTOR',*9999)

      CALL ASSERT(ITYP10(nr).EQ.1,
     '  '>>Only implemented for rc coordinates',ERROR,*9999)

      DO nj=1,3
        NORMAL(nj)=0.0d0
        TANGENT(nj,1)=0.0d0
        TANGENT(nj,2)=0.0d0
      ENDDO !nj

      ne=0
      DO nonp=1,NENP(np,0)
        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
      ENDDO !nonp
      CALL ASSERT(ne.NE.0,'>>Could not find an element in the region',
     '  ERROR,*9999)
      nb=NBJ(1,ne)
      COLLAPSEDNODE=.FALSE.
      HERMSECTOR=.FALSE.
      FOUND=.FALSE.
      nn=1
      DO WHILE(.NOT.FOUND.OR.nn.GT.NNT(nb))
        IF(NPNE(nn,nb,ne).EQ.NP) THEN
          FOUND=.TRUE.
        ELSE
          nn=nn+1
        ENDIF
      ENDDO
      CALL ASSERT(FOUND,'>>Could not find local node in element',
     '  ERROR,*9999)
      DO ni=1,NIT(nb)
        IF(IBT(1,ni,nb).EQ.3.AND.IBT(2,ni,nb).EQ.4) THEN !Herm Sim
          HERMSECTOR=.TRUE.
        ENDIF
      ENDDO
      COLLAPSEDNODE=ISATCOLLAPSE(IBT(1,1,nb),INP(1,1,nb),nb,nn)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Node : '',I5)') np
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      IF(.NOT.COLLAPSEDNODE) THEN
        DO ni=1,MIN(NIT(nb),2)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(''   Xi direction : '',I1)') ni
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(IBT(1,ni,nb).EQ.1.OR.
     '      ((IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6).AND.
     '      IBT(2,ni,nb).NE.4)) THEN !Lagrange element - get ave tangent
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(''   Lagrangian basis function'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''   Number of elements surrounding '
     '          //'node : '',I5)') NENP(np,0)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            DO nj=1,NJT
              TANGENT(nj,ni)=0.0d0
            ENDDO !nj
            DO noelem=1,NENP(np,0)
              ne=NENP(np,noelem)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'(''     Element '',I2,'', ne = '',I5)')
     '            noelem,ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,
     '          *9999)
              DO nj=1,NJT
                nb=NBJ(nj,ne)
                FOUND=.FALSE.
                nn=1
                DO WHILE(.NOT.FOUND.OR.nn.GT.NNT(nb))
                  IF(NPNE(nn,nb,ne).EQ.NP) THEN
                    FOUND=.TRUE.
                  ELSE
                    nn=nn+1
                  ENDIF
                ENDDO
                CALL ASSERT(FOUND,
     '            '>>Could not find local node in element',ERROR,*9999)
                DO nii=1,NIT(nb)
                  XI(nii)=DBLE(INP(nn,nii,nb)-1)/DBLE(IBT(2,nii,nb))
                ENDDO !nii
                nu=NU1(ni)
                TEMP(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '            nu,XI,XE(1,nj))
                TANGENT(nj,ni)=TANGENT(nj,ni)+TEMP(nj)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'(''       nj = '',I2,'', nb = '',I2,'
     '              //''', nn = '',I2,'', NIT(nb) = '',I1,'', nu = '','
     '              //'I1)') nj,nb,nn,NIT(nb),nu
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(''       XI :'',3(1X,D12.4))')
     '              (XI(nii),nii=1,NIT(nb))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(''       Element tangent :'','
     '              //'3(1X,D12.4))') (TEMP(njj),njj=1,NJT)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ENDDO !nj
            ENDDO !ne
            DO nj=1,NJT
              TANGENT(nj,ni)=TANGENT(nj,ni)/DBLE(NENP(np,0))
            ENDDO !nj
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(''   Tangent :'',3(1X,D12.4))')
     '          (TANGENT(nj,ni),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ELSE IF(IBT(1,ni,nb).EQ.2.OR.
     '        (IBT(1,ni,nb).EQ.3.AND.IBT(2,ni,nb).EQ.4).OR.
     '        ((IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6).AND.
     '        IBT(2,ni,nb).EQ.4)) THEN !Hermite element
            DO nj=1,NJT
              TANGENT(nj,ni)=XP(ni+1,1,nj,np)
            ENDDO !nj
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(''   Hermitian basis function'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''   Tangent :'',3(1X,D12.4))')
     '          (TANGENT(nj,ni),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ELSE
            ERROR='>>Not implemented'
            GOTO 9999
          ENDIF
          SUM=0.0d0
          DO nj=1,NJT
            SUM=SUM+TANGENT(nj,ni)**2
          ENDDO !nj
          SUM=DSQRT(SUM)
          IF(DABS(SUM).GT.ZERO_TOL) THEN
            DO nj=1,NJT
              TANGENT(nj,ni)=TANGENT(nj,ni)/SUM
            ENDDO !nj
          ENDIF
        ENDDO !ni
        IF(NJT.EQ.2) THEN
          NORMAL(1)=-TANGENT(2,1)
          NORMAL(2)=TANGENT(1,1)
        ELSE
          NORMAL(1)=TANGENT(2,1)*TANGENT(3,2)-TANGENT(3,1)*TANGENT(2,2)
          NORMAL(2)=TANGENT(3,1)*TANGENT(1,2)-TANGENT(1,1)*TANGENT(3,2)
          NORMAL(3)=TANGENT(1,1)*TANGENT(2,2)-TANGENT(2,1)*TANGENT(1,2)
        ENDIF
        SUM=0.0d0
        DO nj=1,NJT
          SUM=SUM+NORMAL(nj)**2
        ENDDO !nj
        SUM=DSQRT(SUM)
        IF(DABS(SUM).GT.ZERO_TOL) THEN
          DO nj=1,NJT
            NORMAL(nj)=NORMAL(nj)/SUM
          ENDDO !nj
        ELSE
          ERROR='>>Zero length normal vector'
          GOTO 9999
        ENDIF
      ELSE
        IF(NJT.EQ.3) THEN
          IF(HERMSECTOR) THEN
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(''   Hermite Sector basis function'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(nn.EQ.1) THEN
              XI(2)=0.0d0
            ELSE
              XI(2)=1.0d0
            ENDIF
            DO noelem=1,NENP(np,0)
              ne=NENP(np,noelem)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'(''     Element '',I2,'', ne = '',I5)')
     '            noelem,ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,
     '          *9999)
              XI(1)=0.0d0
              DO nj=1,NJT
                TANGENT(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,4,XI,XE(1,nj))
              ENDDO !nj
              XI(1)=1.0d0
              DO nj=1,NJT
                TANGENT(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,4,XI,XE(1,nj))
              ENDDO !nj
              TEMP(1)=TANGENT(2,1)*TANGENT(3,2)-TANGENT(3,1)*
     '          TANGENT(2,2)
              TEMP(2)=TANGENT(3,1)*TANGENT(1,2)-TANGENT(1,1)*
     '          TANGENT(3,2)
              TEMP(3)=TANGENT(1,1)*TANGENT(2,2)-TANGENT(2,1)*
     '          TANGENT(1,2)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'(''     Tangent 1 :'',3(1X,D12.4))')
     '            (TANGENT(nj,1),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''     Tangent 2 :'',3(1X,D12.4))')
     '            (TANGENT(nj,2),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''     Normal :'',3(1X,D12.4))')
     '            (TEMP(nj),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              NORMAL(1)=NORMAL(1)+TEMP(1)
              NORMAL(2)=NORMAL(2)+TEMP(2)
              NORMAL(3)=NORMAL(3)+TEMP(3)
            ENDDO !ne
            IF(nn.EQ.1) THEN
              DO nj=1,NJT
                NORMAL(nj)=-NORMAL(nj)/DBLE(NENP(np,0))
              ENDDO !nj
            ELSE
              DO nj=1,NJT
                NORMAL(nj)=NORMAL(nj)/DBLE(NENP(np,0))
              ENDDO !nj
            ENDIF
          ELSE
            CALL ASSERT(NIT(nb).EQ.2,'>>Not Implemented',ERROR,*9999)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(''   Sector basis function'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6) THEN !Xi1 collapsed
              IF(IBT(1,1,nb).EQ.5) THEN !Collapsed at xi2=0
                XI(2)=0.0d0
              ELSE !Collapsed at xi2=1
                XI(2)=1.0d0
              ENDIF
              nu=4
            ELSE !Xi2 collapsed
              IF(IBT(1,2,nb).EQ.5) THEN !Collapsed at xi1=0
                XI(1)=0.0d0
              ELSE !Collapsed at xi1=1
                XI(1)=1.0d0
              ENDIF
              nu=2
            ENDIF
            DO noelem=1,NENP(np,0)
              ne=NENP(np,noelem)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'(''     Element '',I2,'', ne = '','
     '            //'I5)') noelem,ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '          XA(1,1,ne),XE,XP,ERROR,*9999)
              IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6) THEN !Xi1 collap
                XI(1)=0.0d0
              ELSE
                XI(2)=0.0d0
              ENDIF
              DO nj=1,NJT
                TANGENT(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,nu,XI,XE(1,nj))
              ENDDO !nj
              IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6) THEN !Xi1 collap
                XI(1)=1.0d0
              ELSE
                XI(2)=1.0d0
              ENDIF
              DO nj=1,NJT
                TANGENT(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,nu,XI,XE(1,nj))
              ENDDO !nj
              TEMP(1)=TANGENT(2,1)*TANGENT(3,2)-TANGENT(3,1)*
     '          TANGENT(2,2)
              TEMP(2)=TANGENT(3,1)*TANGENT(1,2)-TANGENT(1,1)*
     '          TANGENT(3,2)
              TEMP(3)=TANGENT(1,1)*TANGENT(2,2)-TANGENT(2,1)*
     '          TANGENT(1,2)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'(''     Tangent 1 :'',3(1X,D12.4))')
     '            (TANGENT(nj,1),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''     Tangent 2 :'',3(1X,D12.4))')
     '            (TANGENT(nj,2),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''     Normal :'',3(1X,D12.4))')
     '              (TEMP(nj),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              NORMAL(1)=NORMAL(1)+TEMP(1)
              NORMAL(2)=NORMAL(2)+TEMP(2)
              NORMAL(3)=NORMAL(3)+TEMP(3)
            ENDDO !ne
C cpb 2/11/96 Need to reverse the normal if the element is collapsed at
C the xi=0 end.
            IF(IBT(1,1,nb).EQ.5.OR.IBT(1,2,nb).EQ.5) THEN !collapse xi=0
              DO nj=1,NJT
                NORMAL(nj)=-NORMAL(nj)/DBLE(NENP(np,0))
              ENDDO !nj
            ELSE !collapse xi=1
              DO nj=1,NJT
                NORMAL(nj)=NORMAL(nj)/DBLE(NENP(np,0))
              ENDDO !nj
            ENDIF
          ENDIF
          SUM=0.0d0
          DO nj=1,NJT
            SUM=SUM+NORMAL(nj)**2
          ENDDO !nj
          SUM=DSQRT(SUM)
          IF(DABS(SUM).GT.ZERO_TOL) THEN
            DO nj=1,NJT
              NORMAL(nj)=NORMAL(nj)/SUM
            ENDDO !nj
          ELSE
            ERROR='>>Zero length normal vector'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF
      ENDIF

      IF((NP_INTERFACE(np,0).GT.1).AND.(nr.NE.NP_INTERFACE(np,1))) THEN
C       Interface node in slave region
        FEMFEMCOUP=ITYP4(nr,nx).EQ.1.AND.
     '    ITYP4(NP_INTERFACE(np,1),nx).EQ.1
        IF(.NOT.FEMFEMCOUP) THEN
          IF(NJT.EQ.2) THEN
C           Reverse sign of normal at the interface
            DO nj=1,NJT
              NORMAL(nj)=-NORMAL(nj)
            ENDDO !nj
          ELSE IF(NJT.EQ.3) THEN
C           Reverse sign of normal and s2 tangent at the interface
            DO nj=1,NJT
              NORMAL(nj)=-NORMAL(nj)
              TANGENT(nj,2)=-TANGENT(nj,2)
            ENDDO !nj
          ENDIF ! njt
        ENDIF !FEM-FEM coupling
      ENDIF ! np_interface

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        IF(NIT(nb).EQ.1) THEN
          WRITE(OP_STRING,'('' Tangent :'',3(1X,D12.4))')
     '      (TANGENT(nj,1),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(NIT(nb).EQ.2) THEN
          WRITE(OP_STRING,'('' Tangent 1 :'',3(1X,D12.4))')
     '      (TANGENT(nj,1),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Tangent 2 :'',3(1X,D12.4))')
     '      (TANGENT(nj,2),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'('' Normal :'',3(1X,D12.4))')
     '    (NORMAL(nj),nj=1,NJT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('GET_TNVECTOR')
      RETURN
9999  CALL ERRORS('GET_TNVECTOR',ERROR)
      CALL EXITS('GET_TNVECTOR')
      RETURN 1
      END


