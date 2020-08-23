      SUBROUTINE GET_TNVECTOR2(IBT,IDO,INP,NBJ,NBJF,ne,nef,NFF,
     '  NKEF,NKJE,NNF,NPF,NPNE,NPNF,nr,NVJE,NVJF,NORMAL,POS,SE,SF,
     '  TAN_DERI,TANGENT,XA,XE,XI,XP,ERROR,*)

C#### Subroutine: GET_TNVECTOR2
C###  Description:
C###    GET_TNVECTOR2 returns the tangent(s) and normal vectors at
C###    a given xi point. It is only implemented for rectangular
C###    cartesian coordinates at the moment.

C*** JWF 14-4-03 Modified so can calculate normal/tangents to faces of
C*** a volume element


      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NBJF(NJM,NFM),ne,nef,NFF(6,NEM),
     '  NKEF(0:4,16,6,NBFM),NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM)
      REAL*8 POS(3),SE(NSM,NBFM,NEM),SF(NSM,NBFM),NORMAL(3),
     '  TAN_DERI(3,2,2),TANGENT(3,2),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XI(2),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nf,ni,nii,nj,njj,nn,nu,nu_deri1,nu_deri2,NU1(0:3),
     &  NKJF(NKM,NNM,NJM)
      REAL*8 PXI,SUM,TEMP(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('GET_TNVECTOR2',*9999)

      CALL ASSERT(ITYP10(nr).EQ.1,
     '  '>>Only implemented for rc coordinates',ERROR,*9999)

      DO nj=1,3
        NORMAL(nj)=0.0d0
        TANGENT(nj,1)=0.0d0
        TANGENT(nj,2)=0.0d0
C*** 21/02/08 JHC initialise tangent derivative vectors         
        TAN_DERI(nj,1,1)=0.0d0
        TAN_DERI(nj,1,2)=0.0d0
        TAN_DERI(nj,2,1)=0.0d0
        TAN_DERI(nj,2,2)=0.0d0
      ENDDO !nj

C     JWF 3-3-02 Modified for contact_mechanics.

      nf=NFF(nef,ne) ! global face #
      IF (ktyp5G(nr).GE.1) THEN
        nb=NBJF(1,nf)
      ELSE
        nb=NBJ(1,ne)
      ENDIF

      DO ni=1,MIN(NIT(nb),2)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(''   Xi direction : '',I1)') ni
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(IBT(1,ni,nb).EQ.1.OR.IBT(1,ni,nb).EQ.2.OR.
     '    (IBT(1,ni,nb).EQ.3.AND.IBT(2,ni,nb).EQ.4).OR.
     '    ((IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6).AND.
     '    IBT(2,ni,nb).NE.4)) THEN !Lagrange or Hermite element
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(''   Lagrangian or Hermitian basis fn'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DO nj=1,NJT
            TANGENT(nj,ni)=0.0d0
C*** 21/02/08 JHC initialise tangent derivative vectors         
            TAN_DERI(nj,ni,1)=0.0d0
            TAN_DERI(nj,ni,2)=0.0d0
          ENDDO !nj

C         JWF 3-3-02 Modified for contact_mechanics.

          IF (ktyp5G(nr).GE.1) THEN
            CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
     '        NBJF(1,nf),nef,NKJE(1,1,1,ne),
     '        NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
     '        NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
            CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,NVJF,
     '       SF,XA(1,1,1),XE,XP,ERROR,*9999)
          ELSE
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '        ERROR,*9999)
          ENDIF

          DO nj=1,NJT

C           JWF 3-3-02 Modified for contact_mechanics.
            IF (ktyp5G(nr).GE.1) THEN
              nb=NBJF(1,nf)
            ELSE
              nb=NBJ(1,ne)
            ENDIF

            nu=NU1(ni)
            TEMP(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        nu,XI,XE(1,nj))
            TANGENT(nj,ni)=TANGENT(nj,ni)+TEMP(nj)

C*** 21/02/08 JHC computing tangent derivative vectors         
            nu_deri1=nu+1
            TAN_DERI(nj,ni,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,nu_deri1,XI,XE(1,nj))

            nu_deri2=6
            IF (ni.EQ.1) THEN
              TAN_DERI(nj,ni,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '          INP(1,1,nb),nb,nu_deri2,XI,XE(1,nj))
            ELSE
              TAN_DERI(nj,2,1)=TAN_DERI(nj,1,2)
            ENDIF

            IF(ni.EQ.1) THEN !find position of projection
              POS(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '          1,XI,XE(1,nj))
            ENDIF
 
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(''       nj = '',I2,'', nb = '',I2,'
     '          //''', nn = '',I2,'', NIT(nb) = '',I1,'', nu = '','
     '          //'I1)') nj,nb,nn,NIT(nb),nu
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''       XI :'',3(1X,D12.4))')
     '          (XI(nii),nii=1,NIT(nb))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''       Element tangent :'','
     '          //'3(1X,D12.4))') (TEMP(njj),njj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDDO !nj
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(''   Tangent :'',3(1X,D12.4))')
     '        (TANGENT(nj,ni),nj=1,NJT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ELSE
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF

C*** 21/02/08 JHC Do not normalise tangent vectors in case of contact. 
C                 However, for tied contact these tangent vectors are normalised in UPDATA.f in fe23.
        IF (ktyp5G(nr).LT.1) THEN ! not contact
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
        ENDIF
C NEWE JHC
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

      CALL EXITS('GET_TNVECTOR2')
      RETURN
9999  CALL ERRORS('GET_TNVECTOR2',ERROR)
      CALL EXITS('GET_TNVECTOR2')
      RETURN 1
      END


