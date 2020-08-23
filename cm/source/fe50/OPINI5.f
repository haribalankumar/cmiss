      SUBROUTINE OPINI5(NBH,NEELEM,NHP,NKH,
     '  NPNODE,nr,NVHP,NW,nx,NYNE,NYNP,FIX,FULL,YP,ERROR,*)

C#### Subroutine: OPINI5
C###  Description:
C###    OPINI5 outputs initial and boundary data.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M),NHP(NPM),
     '  NKH(NHM,NPM,NCM),NPNODE(0:NP_R_M),nr,NVHP(NHM,NPM,NCM),
     '  NW(NEM,3),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),FULL
!     Local Variables
      INTEGER loop,LOOP_MAX,na,NATMAX,nb,nc,ne,nh,NHKMAX,NHKT,nhx,nk,
     '  noelem,nonode,np,nv
      CHARACTER BC_STR(3)*7,BC_TYP(2)*8,CHAR1*2,CHAR2*2,
     '  FIX_STR*200,NUM_STR*10
      LOGICAL AUXIL

      DATA BC_STR /' bcs : ', ' incr: ', ' init: '/
      DATA BC_TYP /' Displac', ' Force  '/

      CALL ENTERS('OPINI5',*9999)

      IF(ITYP6(nr,nx).EQ.1) THEN       !linear
        IF(ITYP5(nr,nx).EQ.1) THEN       !static
          LOOP_MAX=2  !bc's
        ELSE IF(ITYP5(nr,nx).EQ.2) THEN  !time-dependent
          LOOP_MAX=3  !bc's/initial values
        ENDIF
      ELSE IF(ITYP6(nr,nx).EQ.2) THEN  !nonlinear
        LOOP_MAX=3  !bc's/increments/initial values
      ENDIF
C     Find maximum # dof for nodes
      NHKMAX=0
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO loop=1,LOOP_MAX
          DO nc=1,MIN(2,NCM)
            NHKT=0
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,nc)
                DO nk=1,NKH(nh,np,nc)
                  NHKT=NHKT+1
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF(NHKT.GT.NHKMAX) NHKMAX=NHKT
        ENDDO
      ENDDO

      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO loop=1,LOOP_MAX
          DO nc=1,MIN(2,NCM)
            NUM_STR='          '
            IF(nc.EQ.1.AND.loop.EQ.1)
     '        WRITE(NUM_STR,'(A6,I4)') ' Node ',np
            NHKT=0
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,nc)
                DO nk=1,NKH(nh,np,nc)
                  NHKT=NHKT+1
                ENDDO
              ENDDO
            ENDDO
            WRITE(CHAR1,'(I2)') NHKT
            WRITE(CHAR2,'(I2)') NHKMAX-NHKT+1
            IF(loop.EQ.3) THEN !don't print FIX array for init conds
              FORMAT='('//CHAR1(1:2)//'X,' //CHAR2(1:2)//'X)'
              WRITE(FIX_STR,FORMAT)
            ELSE !print FIX array for all others
              FORMAT='('//CHAR1(1:2)//'L1,'//CHAR2(1:2)//'X)'
              WRITE(FIX_STR,FORMAT)
     '          (((FIX(NYNP(nk,nv,NH_LOC(nhx,nx),np,0,nc,nr),loop),nk=1,
     '          NKH(NH_LOC(nhx,nx),np,nc)),
     '          nv=1,NVHP(NH_LOC(nhx,nx),np,nc)),nhx=1,NHP(np))
            ENDIF
            FORMAT='('''//NUM_STR(1:10)//BC_TYP(nc)(1:8)
     '        //BC_STR(loop)(1:7)//FIX_STR(1:NHKMAX)
     '        //''',4D12.4,:,/(24X,'//CHAR1(1:2)//'X,'
     '        //CHAR2(1:2)//'X,4D12.4))'
            WRITE(OP_STRING,FORMAT)
     '          (((YP(NYNP(nk,nv,NH_LOC(nhx,nx),np,0,nc,nr),loop),nk=1,
     '          NKH(NH_LOC(nhx,nx),np,nc)),
     '          nv=1,NVHP(NH_LOC(nhx,nx),np,nc)),nhx=1,NHP(np))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nc
        ENDDO !loop
      ENDDO !nonode (np)

C LC 25/2/97 archived section :
C!!!! MPN 15-Sep-95: not used anymore

      AUXIL=.FALSE.
      DO nb=1,NBFT
        IF(NAT(nb).GT.0) AUXIL=.TRUE.
      ENDDO
      IF(AUXIL) THEN
        NATMAX=0
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          DO nhx=1,NH_LOC(0,nx)
            nh=NH_LOC(nhx,nx)
            IF(NAT(NBH(nh,1,ne)).GT.NATMAX) NATMAX=NAT(NBH(nh,1,ne))
          ENDDO
        ENDDO

        WRITE(OP_STRING,
     '    '(/'' Auxiliary element variables:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          WRITE(OP_STRING,'('' Element '',I5,'', NW='',I1)') ne,NW(ne,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO loop=1,LOOP_MAX
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              WRITE(CHAR1,'(I2)') NAT(NBH(nh,1,ne))
              WRITE(CHAR2,'(I2)') NATMAX-NAT(NBH(nh,1,ne))+1
              IF(loop.EQ.3) THEN !don't print FIX array for init conds
                WRITE(FIX_STR,'('//CHAR1(1:2)//'X,' //CHAR2(1:2)//'X)')
              ELSE !print FIX array for all others
                WRITE(FIX_STR,'('//CHAR1(1:2)//'L1,'//CHAR2(1:2)//'X)')
     '            (FIX(NYNE(na,nh,0,1,ne),loop),na=1,NAT(NBH(nh,1,ne)))
              ENDIF
              FORMAT='(''  Dep var '',I1,'''
     '          //BC_STR(loop)(1:7)//FIX_STR(1:NATMAX)
     '          //''',5D12.4,:,/(17X,'//CHAR1(1:2)//'X,'
     '          //CHAR2(1:2)//'X,5D12.4))'
              WRITE(OP_STRING,FORMAT)
     '          nh,(YP(NYNE(na,nh,0,1,ne),loop),na=1,NAT(NBH(nh,1,ne)))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nh
          ENDDO !loop

C            WRITE(OP_STRING,
C     '        '(/'' Element'',I3,'' Aux param '',I2,'': '','
C     '        //'4D12.4/, (25X,4D12.3))')
C     '        ne,ne,(YP(NYNE(na,1,0,1,ne),1),na=1,NAT(NBH(NH_LOC(1,nx),1,ne)))
C            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            DO nh=1,NH_LOC(0,nx)
C              WRITE(OP_STRING,'(18X,I1,6X,5D11.3/,(25X,5D11.3))')
C     '          nh,(YP(NYNE(na,nh,0,1,ne),1),na=1,NAT(NBH(nh,1,ne)))
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            ENDDO
        ENDDO !noelem (ne)
      ENDIF

      IF(KTYP57(nr).EQ.2) THEN
        WRITE(OP_STRING,
     '    '(/'' Boundary pressure increments have been entered'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP51(nr).EQ.3) THEN !3D
          IF(KTYP5A(nr).EQ.1) THEN
            WRITE(OP_STRING,'('' Free hydrostatic pressure vars are '
     '        //'determined using incompressibility constraints.'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(KTYP5A(nr).EQ.2) THEN
            WRITE(OP_STRING,'('' Free hydrostatic pressure vars are '
     '        //'determined by matching normal Cauchy stress on Xi3 '
     '        //'faces.'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      IF(KTYP58(nr).EQ.2) THEN
        WRITE(OP_STRING,
     '    '(/'' Elements are isochoric: the deformed state metric '
     '    //'tensors ''/'' are evaluated at the Gauss pts with I3=1'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(FULL) THEN !extended output
      ENDIF

      CALL EXITS('OPINI5')
      RETURN
 9999 CALL ERRORS('OPINI5',ERROR)
      CALL EXITS('OPINI5')
      RETURN 1
      END


