      SUBROUTINE OPINI4(NHP,NKH,NPNODE,nr,NVHP,nx,NYNP,YP,FIX,ERROR,*)

C#### Subroutine: OPINI4
C###  Description:
C###    OPINI4 outputs initial and boundary data.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NHP(NPM),NKH(NHM,NPM,NCM),NPNODE(0:NP_R_M),nr,
     '  NVHP(NHM,NPM,NCM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*)
!     Local Variables
      INTEGER loop,LOOP_MAX,nc,nh,nhx,
     '  NHKT,nk,nonode,np,nv
      CHARACTER BC_STR*7,BC_TYP*8,CHAR*2,FORMAT*200,NODE_STR*10

      CALL ENTERS('OPINI4',*9999)

      IF(ITYP6(nr,nx).EQ.1) THEN       !linear
        IF(ITYP5(nr,nx).EQ.1) THEN       !static
          LOOP_MAX=1  !bc's
        ELSE IF(ITYP5(nr,nx).EQ.2) THEN  !time-depentdent
          LOOP_MAX=3  !bc's/initial values
        ENDIF
      ENDIF
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO loop=1,LOOP_MAX
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
              BC_TYP=' Displac'
            ELSE IF(nc.EQ.2) THEN
              NODE_STR='          '
              BC_TYP=' Force  '
            ENDIF
            NHKT=0
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,nc)
                DO nk=1,NKH(nh,np,nc)
                  NHKT=NHKT+1
                ENDDO
              ENDDO
            ENDDO
            WRITE(CHAR,'(I2)') NHKT
            FORMAT='('''//NODE_STR(1:10)//BC_TYP(1:8)
     '        //BC_STR(1:7)//''','//CHAR(1:2)
     '        //'L1,5D12.4,:,/(25X,'//CHAR(1:2)//'X,5D12.4))'
            WRITE(OP_STRING,FORMAT)
     '        (((FIX(NYNP(nk,nv,NH_LOC(nhx,nx),np,0,nc,nr),loop),nk=1,
     '        NKH(NH_LOC(nhx,nx),np,nc)),
     '        nv=1,NVHP(NH_LOC(nhx,nx),np,nc)),nhx=1,NHP(np)),
     '        (((YP(NYNP(nk,nv,NH_LOC(nhx,nx),np,0,nc,nr),loop),
     '        nk=1,NKH(NH_LOC(nhx,nx),np,nc)),
     '        nv=1,NVHP(NH_LOC(nhx,nx),np,nc)),nhx=1,NHP(np))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nc
        ENDDO !loop
      ENDDO !nonode (np)

      CALL EXITS('OPINI4')
      RETURN
 9999 CALL ERRORS('OPINI4',ERROR)
      CALL EXITS('OPINI4')
      RETURN 1
      END


