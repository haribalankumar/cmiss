      SUBROUTINE OPINI9(NHP,NKH,NPNODE,nr,NVHP,nx,NYNP,FIX,FULL,YP,
     '  ERROR,*)

C#### Subroutine: OPINI9
C###  Description:
C###    OPINI9 outputs initial and boundary data for BE problems.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NHP(NPM),NKH(NHM,NPM,NCM),NPNODE(0:NP_R_M),nr,
     '  NVHP(NHM,NPM,NCM),nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),FULL
!     Local Variables
      INTEGER loop,LOOP_MAX,LOOP_STEP,nc,nh,NHKT,nhx,nk,nonode,np,nv
      CHARACTER BC_STR*7,BC_TYP*8,CHAR*2,CHAR1*1,
     '  FORMAT*200,NODE_STR*10
      CALL ENTERS('OPINI9',*9999)

      IF(ITYP6(nr,nx).EQ.1) THEN       !linear
        IF(ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4) THEN !static or quasi
          LOOP_MAX=1  !bc's
          LOOP_STEP=2
        ELSE IF(ITYP5(nr,nx).EQ.2) THEN  !time-depentdent
          LOOP_MAX=3  !bc's/initial values
          LOOP_STEP=2
        ENDIF
      ENDIF
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO loop=1,LOOP_MAX,LOOP_STEP
          IF(loop.EQ.1) THEN
            WRITE(NODE_STR,'(A6,I4)') ' Node ',np
            BC_STR=' bcs : '
          ELSE IF(loop.EQ.2) THEN
            NODE_STR='          '
            BC_STR=' incr: '
          ELSE IF(loop.EQ.3) THEN
            NODE_STR='          '
            BC_STR=' init: '
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
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,nc)
                DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                  NHKT=NHKT+1
                ENDDO
              ENDDO
            ENDDO
            WRITE(CHAR,'(I2)') NHKT
            FORMAT='('''//NODE_STR(1:10)//BC_TYP(1:8)
     '        //BC_STR(1:7)//''','//CHAR(1:2)
     '        //'L1,5E12.4,:,/(25X,'//CHAR(1:2)//'X,5E12.4))'
            WRITE(OP_STRING,FORMAT)
     '        (((FIX(NYNP(nk,nv,NH_LOC(nhx,nx),np,0,nc,nr),loop),nk=1,
     '        MAX(NKH(NH_LOC(nhx,nx),np,nc)-KTYP93(nc,nr),1)),
     '        nv=1,NVHP(NH_LOC(nhx,nx),np,nc)),
     '        nhx=1,NHP(np)),
     '        (((YP(NYNP(nk,nv,NH_LOC(nhx,nx),np,0,nc,nr),loop),
     '        nk=1,MAX(NKH(NH_LOC(nhx,nx),np,nc)-KTYP93(nc,nr),1)),nv=1,
     '        NVHP(NH_LOC(nhx,nx),np,nc)),nhx=1,NHP(np))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO ! nc
        ENDDO ! loop
      ENDDO ! nonode (np)

      WRITE(OP_STRING,'(/'' Nodal parameters''/,1X,16(''=''),'
     '  //'/14X,''#Variables  #Derivs/var  #Equations  #Derivs/eqn'','
     '  //'/14X,''  NHP(np)   NKH(nh,np,1)   NHP(np)   NKH(nh,np,2)'
     '  //''')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        WRITE(CHAR1,'(I1)') NHP(np)
        FORMAT='(1X,''Node'',I5,10X,I1,10X,'//CHAR1//'I3,'
     '                       //'10X,I1,10X,'//CHAR1//'I3)'
        WRITE(OP_STRING,FORMAT) np,
     '     NHP(np),(NKH(NH_LOC(nhx,nx),np,1),nhx=1,NHP(np)),
     '     NHP(np),(NKH(NH_LOC(nhx,nx),np,2),nhx=1,NHP(np))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      IF(KTYP93(1,nr).GT.0) THEN !no cross derivative information
        WRITE(OP_STRING,'(/'' ***   KTYP93(1,nr) = '',I2,/'
     '    //'''       No cross derivative variables will '
     '    //'be solved for'')') KTYP93(1,nr)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(FULL) THEN !extended output
      ENDIF

      CALL EXITS('OPINI9')
      RETURN
 9999 CALL ERRORS('OPINI9',ERROR)
      CALL EXITS('OPINI9')
      RETURN 1
      END


