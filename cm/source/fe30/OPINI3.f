      SUBROUTINE OPINI3(ITHRES,NBH,NBJ,NEELEM,NHE,NHP,
     '  NKH,NPNODE,nr,NVHP,NWQ,nx,NYNP,FIX,FULL,AQ,THRES,YG,
     '  YP,YQ,ERROR,*)

C#### Subroutine: OPINI3
C###  Description:
C###    OPINI3 outputs initial and boundary data for FE30 problems.

      IMPLICIT NONE

      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER ITHRES(3,NGM,NEM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M),NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NPNODE(0:NP_R_M),nr,NVHP(NHM,NPM,NCM),NWQ(8,0:NQM,NAM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 AQ(NMAQM,NQM),THRES(3,NGM,NEM),YG(NIYGM,NGM,NEM),
     '  YP(NYM,NIYM),YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),FULL
!     Local Variables
      INTEGER loop,LOOP_MAX,LOOP_STEP,maqp1t0,maqp1t1,maqp1i,
     '  maqp2t0,maqp2t1,maqp2i,
     '  nc,ne,ng,nh,NHKT,nhx,nk,noelem,nonode,np,nq,nv
      CHARACTER BC_STR*8,BC_TYP*8,CHAR*2,CHAR1*1,
     '  FORMAT*200,NODE_STR*10

      CALL ENTERS('OPINI3',*9999)

      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.8) THEN
! Threshold Modelling
        WRITE(OP_STRING,'(/'' Activation array:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          WRITE(OP_STRING,'('' ITHRES(1,ng,'',I4,''): '',/(50I2))')
     '      ne,(ITHRES(1,ng,ne),ng=1,NGT(NBJ(1,ne)))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(/'' Purkinje tissue array:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          WRITE(OP_STRING,'('' ITHRES(2,ng,'',I4,''): '',/(50I2))')
     '      ne,(ITHRES(2,ng,ne),ng=1,NGT(NBJ(1,ne)))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(/'' Absolute activation times:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          WRITE(OP_STRING,
     '      '('' YG(1,ng,'',I4,''): '',/(10(1X,E12.4)))')
     '      ne,(YG(1,ng,ne),ng=1,NGT(NBJ(1,ne)))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(/'' Time since activation:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          WRITE(OP_STRING,
     '      '('' THRES(1,ng,'',I4,''): '',/(10(1X,E12.4)))')
     '      ne,(THRES(1,ng,ne),ng=1,NGT(NBJ(1,ne)))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

      ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN
! cellular based modelling

        IF (ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7) THEN
!electrical or coupled model
          IF (ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).LT.5) THEN
            WRITE(OP_STRING,'('' Current injection conditions: '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            CALL ASSERT(NMAQM.GE.6,'>>Increase NMAQM, must be >= 6',
     '        ERROR,*9999)

C         Call maq_loc to get the three indicies for each pulse
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t0,
     '        MAQ_START,ERROR,*9999)
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t1,
     '        MAQ_STOP,ERROR,*9999)
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,
     '        MAQ_CURRENT,ERROR,*9999)
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t0,
     '        MAQ_START,ERROR,*9999)
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t1,
     '        MAQ_STOP,ERROR,*9999)
            CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2i,
     '        MAQ_CURRENT,ERROR,*9999)

            DO nq=1,NQT
              IF(AQ(maqp1t1,nq).GT.ZERO_TOL) THEN
                WRITE(OP_STRING,'('' Grid point '',I6,''   activated at'
     '            //' '',F8.2,'' ms (Pulse 1)'')') nq,AQ(maqp1t0,nq)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Grid point '',I6,'' deactivated at'
     '            //' '',F8.2,'' ms (Pulse 1)'')') nq,AQ(maqp1t1,nq)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Current injected: '',F8.4)')
     '            AQ(maqp1i,nq)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(AQ(maqp2t1,nq).GT.ZERO_TOL) THEN
                WRITE(OP_STRING,'('' Grid point '',I6,''   activated at'
     '            //' '',F8.2,'' ms (Pulse 2)'')') nq,AQ(maqp2t0,nq)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Grid point '',I6,'' deactivated at'
     '            //' '',F8.2,'' ms (Pulse 2)'')') nq,AQ(maqp2t1,nq)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Current injected: '',F8.4)')
     '            AQ(maqp2i,nq)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDIF !ITYP3(nr,nx).LT.5
        ENDIF !ITYP19(nr,nx).EQ.1.OR.7

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volumes also
      ELSE IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '        ITYP4(nr,nx).EQ.7) THEN 
! Collocation
        DO nq=1,NQT
          IF(NWQ(1,nq,1).EQ.0) THEN !nq is internal point
            IF(DABS(YQ(nq,1,1)).GT.1.D-6) THEN
              WRITE(OP_STRING,'('' YQ('',I6,'',1,1): '',E12.4,'
     '          //''' Non-zero i.c.'')') nq,YQ(nq,1,1)
            ELSE
              WRITE(OP_STRING,'('' YQ('',I6,'',1,1): '',E12.4)')
     '          nq,YQ(nq,1,1)
            ENDIF
          ELSE IF(NWQ(1,nq,1).GT.0) THEN !nq is boundary point
            IF(NWQ(5,nq,1).EQ.1) THEN      !Dirichlet b.c.
              WRITE(OP_STRING,'('' YQ('',I6,'',1,1): '',E12.4,'
     '          //''' Boundary point: Dirichlet'')') nq,YQ(nq,1,1)
            ELSE IF(NWQ(5,nq,1).EQ.2) THEN !Neumann b.c.
              WRITE(OP_STRING,'('' YQ('',I6,'',1,1): '',E12.4,'
     '          //''' Boundary point: Neumann'')') nq,YQ(nq,1,1)
            ENDIF
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO


C PM 26-JUL-01 :for flow in elastic tube
      ELSE IF (ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.5.AND.
     '  ITYP3(nr,nx).EQ.1) THEN
                                            ! Flow in elastic tube

        WRITE(OP_STRING,'('' initial pressure: '',E12.4)') Po
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF (VENOUS_NETWORK.EQ.'N') THEN
          WRITE(OP_STRING,'('' No venous network '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'('' Venous side is included and '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF (MICRO_NETWORK.EQ.'Y') THEN
            WRITE(OP_STRING,'('' coupled via micro-circulation ' //
     '        'network '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' coupled directly '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE
!All except cardiac activation and collocation
        IF(ITYP6(nr,nx).EQ.1) THEN       !linear
          IF(ITYP5(nr,nx).EQ.1) THEN       !static
            LOOP_MAX=1  !bc's
            LOOP_STEP=1
          ELSE IF(ITYP5(nr,nx).EQ.2) THEN  !time-depentdent
            LOOP_MAX=3  !bc's/initial values
            LOOP_STEP=2
          ENDIF
        ELSE IF(ITYP5(nr,nx).EQ.5) THEN !wavefront paths
          LOOP_MAX=3 !bc's/initial values
          LOOP_STEP=2 !no increments
        ELSE IF(ITYP6(nr,nx).EQ.2) THEN  !nonlinear
          LOOP_STEP=1
          IF(ITYP5(nr,nx).EQ.1) THEN       !static
            LOOP_MAX=2  !bc's/increments
          ELSE IF(ITYP5(nr,nx).EQ.2) THEN  !time-depentdent
            LOOP_MAX=3  !bc's/increments/initial values
          ENDIF
        ENDIF
        DO nonode=1,NPNODE(0)
          np=NPNODE(nonode)
          DO loop=1,LOOP_MAX,LOOP_STEP
            IF(loop.EQ.1) THEN !essential and flux bcs
              WRITE(NODE_STR,'(A6,I4)') ' Node ',np
              BC_STR=' bcs  : '
            ELSE IF(loop.EQ.2) THEN !incremental bcs
              NODE_STR='          '
              BC_STR=' incr.: '
            ELSE IF(loop.EQ.3) THEN !initial conditions
              NODE_STR='          '
              BC_STR=' init.: '
            ENDIF
            DO nc=1,MIN(2,NCM)
              IF(nc.EQ.1) THEN
                BC_TYP=' Essent.'
              ELSE IF(nc.EQ.2) THEN
                NODE_STR='          '
                BC_TYP=' Natural'
              ENDIF
              NHKT=0
              DO nhx=1,NHP(np)
                nh=nh_loc(nhx,nx)
                DO nv=1,NVHP(nh,np,nc)
                  DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                    NHKT=NHKT+1
                  ENDDO
                ENDDO
              ENDDO
              WRITE(CHAR,'(I2)') NHKT
              FORMAT='('''//NODE_STR(1:10)//BC_TYP(1:8)
     '          //BC_STR(1:8)//''','//CHAR(1:2)
     '          //'L1,5E12.4,:,/(26X,'//CHAR(1:2)//'X,5E12.4))'
              WRITE(OP_STRING,FORMAT)
     '          (((FIX(NYNP(nk,nv,nh_loc(nhx,nx),np,0,nc,nr),loop),nk=1,
     '          MAX(NKH(nh_loc(nhx,nx),np,nc)-KTYP93(nc,nr),1)),nv=1,
     '          NVHP(nh_loc(nhx,nx),np,nc)),nhx=1,NHP(np)),
     '          (((YP(NYNP(nk,nv,nh_loc(nhx,nx),np,0,nc,nr),loop),NK=1,
     '          MAX(NKH(nh_loc(nhx,nx),np,nc)-KTYP93(nc,nr),1)),nv=1,
     '          NVHP(nh_loc(nhx,nx),np,nc)),nhx=1,NHP(np))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nc
          ENDDO !loop
        ENDDO !nonode (np)

        WRITE(OP_STRING,'(/'' Nodal parameters''/,1X,16(''=''),'
     '    //'/14X,''#Variables  #Derivs/var  #Equations  #Derivs/eqn'','
     '    //'/14X,''  NHP(np)   NKH(nh,np,1)   NHP(np)   NKH(nh,np,2)'
     '    //''')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nonode=1,NPNODE(0)
          np=NPNODE(nonode)
          WRITE(CHAR1,'(I1)') NHP(np)
          FORMAT='(1X,''Node'',I5,10X,I1,10X,'//CHAR1//'I3,'
     '                         //'10X,I1,10X,'//CHAR1//'I3)'
          WRITE(OP_STRING,FORMAT) np,
     '       NHP(np),(NKH(nh_loc(nhx,nx),np,1),nhx=1,NHP(np)),
     '       NHP(np),(NKH(nh_loc(nhx,nx),np,2),nhx=1,NHP(np))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

        WRITE(OP_STRING,'(/'' Auxillary parameters''/,1X,20(''=''),'
     '    //'/17X,''   Variable       Equation   '','
     '    //'/17X,''(#Parameters)   (#Aux eqtns)'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          WRITE(CHAR1,'(I1)') NHE(ne)
          FORMAT='(1X,''Element'',I5,'
     '      //'8X,'//CHAR1//'(1X,I1,''('',I1,'')''),'
     '      //'8X,'//CHAR1//'(1X,I1,''('',I1,'')''))'
          WRITE(OP_STRING,FORMAT) ne,
     '      (nh,NAT(NBH(nh_loc(nhx,nx),1,ne)),nhx=1,NHE(ne)),
     '      (nh,NAT(NBH(nh_loc(nhx,nx),2,ne)),nhx=1,NHE(ne))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

        IF(FULL) THEN !extended output
        ENDIF

      ENDIF

      CALL EXITS('OPINI3')
      RETURN
 9999 CALL ERRORS('OPINI3',ERROR)
      CALL EXITS('OPINI3')
      RETURN 1
      END

      

