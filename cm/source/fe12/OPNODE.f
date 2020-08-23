      SUBROUTINE OPNODE(nc,NHP,NKH,NKJ,
     '  NP_INTERFACE,NPLIST,NPNODE,nr,
     '  NVJP,NWP,nx,XP,ZP,EXTEND,INTERFACE,TYPE,ERROR,*)

C#### Subroutine: OPNODE
C###  Description:
C###    OPNODE outputs nodal coordinates for all nodes in region nr.

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER nc,NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM),NKJ(NJM,NPM),
     '  NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),
     '  NPNODE(0:NP_R_M,0:NRM),NWP(NPM,2),nr,
     '  NVJP(NJM,NPM),nx
      REAL*8 XP(NKM,NVM,NJM,NPM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL EXTEND,INTERFACE
!     Local Variables
      INTEGER IB2,IE2,IB6,IE6,IBEG,IEND,nj,NKJMAX,
     '  nonode,np
      REAL*8 DUMMY(1)
      CHARACTER CHAR2*3,CHAR6*6,OPT3(5)*21

      DATA OPT3(1)/'rectangular cartesian'/,
     '     OPT3(2)/'cylindrical polar    '/,
     '     OPT3(3)/'spherical polar      '/,
     '     OPT3(4)/'prolate spheroidal   '/,
     '     OPT3(5)/'oblate  spheroidal   '/

      CALL ENTERS('OPNODE',*9999)


C LKC 2-MAY-1999 Initialise DUMMY(1)
      DUMMY(1)=0.D0
      IF(NPLIST(0).EQ.0) THEN
        WRITE(OP_STRING,'(/'' No nodes defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        GO TO 9998
      ENDIF

c cpb 7/12/94 Needs to be redone for corners

      WRITE(CHAR2,'(I2)') nr
      CALL STRING_TRIM(CHAR2,IB2,IE2)
      WRITE(CHAR6,'(I6)') NPNODE(0,nr)
      CALL STRING_TRIM(CHAR6,IB6,IE6)
      WRITE(OP_STRING,'(/''  Number of nodes in region '//CHAR2(IB2:IE2)
     '  //' is '//CHAR6(IB6:IE6)//'. Listing nodes:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      CALL WRITE_LONG(INTTYPE,1,1,IOFI,NPLIST(0),10,10,NPLIST(1),
     '  DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)

      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(ITYP10(nr).EQ.4) THEN
        WRITE(OP_STRING,'(''  Focus position = '',E11.4)') FOCUS
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(EXTEND) THEN
        NKJMAX=0
        np=NPNODE(1,nr)
        DO nj=1,NJT
          IF(NKJ(nj,np).GT.NKJMAX) NKJMAX=NKJ(nj,np)
        ENDDO
        WRITE(CHAR2,'(I3)') 14+12*MIN(NKJMAX,4)
        CALL STRING_TRIM(TYPE,IBEG,IEND)
        FORMAT='(/''  Nodal coordinates:'','//CHAR2//'X,'''
     '    //TYPE(IBEG:IEND)//':'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'(/''  Nodal coordinates: ('
     '    //OPT3(ITYP10(nr))//')'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      DO nonode=1,NPLIST(0)
        np=NPLIST(nonode)
        CALL OPNODE1(NHP,NKH(1,1,nc),NKJ(1,np),np,
     '    NP_INTERFACE,nr,NVJP,NWP,nx,
     '    XP(1,1,1,np),ZP(1,1,1,np,nc),EXTEND,INTERFACE,
     '    ERROR,*9999)
      ENDDO

C KAT 18Jun98: Element based parameters should be output in OPELEM
C      IF(NQT.GT.0) THEN
C        AUXILLIARY=.FALSE.
C        DO nb=1,NBFT
C          IF(NAT(nb).GT.0) AUXILLIARY=.TRUE.
C        ENDDO
C        IF(AUXILLIARY) THEN
C          WRITE(OP_STRING,'(/''  Element parameters XA(1,nj,nq):'')')
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C          DO nq=1,NQT
C            WRITE(OP_STRING,'(3X,''Element '',I5,'':'',3X,6D13.4)')
C     '        nq,(XA(1,nj,nq),nj=1,NJ_LOC(0,0,0))
C            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C          ENDDO
C        ENDIF
C      ENDIF

      WRITE(OP_STRING,'(/4X,''Highest node # in region #'',I3,'
     '  //''' is '',I6)') nr,NPT(nr)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'( 4X,''Highest node # in all regions is '','
     '  //'I6)') NPT(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

 9998 CALL EXITS('OPNODE')
      RETURN
 9999 CALL ERRORS('OPNODE',ERROR)
      CALL EXITS('OPNODE')
      RETURN 1
      END


