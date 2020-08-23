      SUBROUTINE NX_LOC(ACTION,CLASS,nx,TYPE,ERROR,*)

C#### Subroutine: NX_LOC
C###  Description:
C###    NX_LOC returns nx number.
C###  TYPE is NX_FIT   for fitting problems
C###          NX_OPTI  for optimisation problems
C###          NX_SOLVE for solution problems
C###  CLASS is a user specified number to indicate differentiate
C###    between problems of the same type.
C###  ACTION is NX_ALLOCATE to allocate an nx# for a problem type
C###            NX_ALLOCATE_AND_LOCK to allocate and lock an nx#
C###                        for a problem type
C###            NX_FREE is to free a locked allocated nx#
C###            NX_INQUIRE to inquire about the nx# for a problem type
C###  If ACTION is NX_ALLOCATE_AND_LOCK the nx# is
C###  locked until it is freed with another call to NX_LOC with an
C###  ACTION of NX_FREE

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ACTION,CLASS,nx,TYPE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LOCK,nxx
      LOGICAL FOUND

      CALL ENTERS('NX_LOC',*9999)

      IF(ACTION.EQ.NX_ALLOCATE) THEN
        LOCK=NX_NOLOCK
      ELSE IF(ACTION.EQ.NX_ALLOCATE_AND_LOCK) THEN
        LOCK=NX_LOCK
      ENDIF

      IF(TYPE.EQ.NX_FIT.OR.TYPE.EQ.NX_OPTI.OR.TYPE.EQ.NX_SOLVE) THEN

        IF(ACTION.EQ.NX_ALLOCATE.OR.ACTION.EQ.NX_ALLOCATE_AND_LOCK) THEN
C
C Try and find a free nx and check that the nx has not already been
C allocated for the type
C
          FOUND=.FALSE.
          nx=0
          DO WHILE(nx.LT.NXM.AND.(.NOT.FOUND))
            nx=nx+1
            IF((NX_TYPE(nx).EQ.0).OR.(NX_LOCKS(nx).EQ.NX_NOLOCK))
     '        FOUND=.TRUE.
          ENDDO
          IF(FOUND) THEN
            NX_CLASS(nx)=CLASS
            NX_TYPE(nx)=TYPE
            NX_LOCKS(nx)=LOCK
            NX_LIST(0)=NX_LIST(0)+1
            NX_LIST(NX_LIST(0))=nx
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Allocate nx= '',I2,'' for class= '','
     '          //'I2,'', type= '',I2)') nx,CLASS,TYPE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            nx=0
            ERROR='>>Increase NXM'
            GOTO 9999
          ENDIF

        ELSE IF(ACTION.EQ.NX_FREE) THEN
          FOUND=.FALSE.
          nxx=0
          DO WHILE(nxx.LT.NX_LIST(0).AND.(.NOT.FOUND))
            nxx=nxx+1
            nx=NX_LIST(nxx)
            IF((NX_TYPE(nx).EQ.TYPE).AND.NX_CLASS(nx).EQ.CLASS)
     '        FOUND=.TRUE.
          ENDDO
C CPB 7/6/94 Should a free action about a non-allocated nx for a
C particular type be an error???
          IF(.NOT.FOUND) THEN
            ERROR='>>Cannot free nx, none allocated for the '
     '        //'problem type'
            GOTO 9999
          ELSE
C MPN 20/6/94 Compact NX_LIST's and NX_LOCKS by removing NXth location
            NX_LIST(0)=NX_LIST(0)-1
            DO nxx=nx,NX_LIST(0)
              NX_LIST(nxx)=NX_LIST(nxx+1)
            ENDDO
            NX_LIST(NX_LIST(0)+1)=0
            NX_CLASS(nx)=0
            NX_LOCKS(nx)=NX_NOLOCK
            NX_TYPE(nx)=0
          ENDIF

        ELSE IF(ACTION.EQ.NX_INQUIRE) THEN
          FOUND=.FALSE.
          nxx=0
          DO WHILE(nxx.LT.NX_LIST(0).AND.(.NOT.FOUND))
            nxx=nxx+1
            nx=NX_LIST(nxx)
            IF((NX_TYPE(nx).EQ.TYPE).AND.NX_CLASS(nx).EQ.CLASS)
     '        FOUND=.TRUE.
          ENDDO
          IF(.NOT.FOUND) nx=0
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Inquire: nx= '',I2,'' for class= '','
     '        //'I2,'', type= '',I2)') nx,CLASS,TYPE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE
          nx=0
          ERROR='>>Invalid nx action'
          GOTO 9999
        ENDIF

      ELSE
        nx=0
        ERROR='>>Invalid nx type'
        GOTO 9999
      ENDIF

      CALL EXITS('NX_LOC')
      RETURN
 9999 CALL ERRORS('NX_LOC',ERROR)
      CALL EXITS('NX_LOC')
      RETURN 1
      END


