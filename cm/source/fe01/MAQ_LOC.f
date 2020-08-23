      SUBROUTINE MAQ_LOC(ACTION,TYPE,maq,SUBTYPE,ERROR,*)

C#### Subroutine: MAQ_LOC
C###  Description:
C###    <HTML>
C###    MAQ_LOC returns maq number.
C###    <PRE>
C###    TYPE:    MAQ_I_PULSE1 for first pulse current input conditions
C###             MAQ_I_PULSE2 for second pulse current input conditions
C###             MAQ_TIME for time information
C###             MAQ_COORD for stored grid point coordinates
C###
C###    SUBTYPE: MAQ_START is the start time for current injection
C###             MAQ_STOP is the stop time for current injection
C###             MAQ_CURRENT is the amplitude of the injected current
C###             MAQ_ACTIV_TIME is the activation time of a grid point
C###             MAQ_M_DPOT is the maximum change in potential
C###             MAQ_X_UNDEF is the undeformed x grid point position
C###             MAQ_Y_UNDEF is the undeformed y grid point position
C###             MAQ_Z_UNDEF is the undeformed z grid point position
C###             MAQ_RHO_X is the x component of the dipole strength
C###             MAQ_RHO_Y is the y component of the dipole strength
C###             MAQ_RHO_Z is the z component of the dipole strength
C###             MAQ_NORMAL_X is the x component of the grid point normal
C###             MAQ_NORMAL_Y is the y component of the grid point normal
C###             MAQ_NORMAL_Z is the z component of the grid point normal
C###
C###    ACTION:  MAQ_ALLOCATE to allocate an maq# for a grid property
C###             MAQ_ALLOCATE_AND_LOCK to allocate and lock an maq#
C###               for a grid property
C###             MAQ_FREE is to free a locked allocated maq#
C###             MAQ_INQUIRE to inquire about the maq# for type/subtype
C###  If ACTION is MAQ_ALLOCATE_AND_LOCK the maq# is
C###  locked until it is freed with another call to MAQ_LOC with an
C###  ACTION of MAQ_FREE.
C###  </PRE>
C###  </HTML>
C***  Martin Buist 11 February 1998
C***  Greg Sands   24 July 2003

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'maqloc00.cmn'
      INCLUDE 'maqloc00.inc'
!     Parameter List
      INTEGER ACTION,maq,SUBTYPE,TYPE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LOCK,maqq
      LOGICAL FOUND

      CALL ENTERS('MAQ_LOC',*9999)

      IF(ACTION.EQ.MAQ_ALLOCATE) THEN
        LOCK=MAQ_NOLOCK
      ELSE IF(ACTION.EQ.MAQ_ALLOCATE_AND_LOCK) THEN
        LOCK=MAQ_LOCK
      ENDIF

      IF((TYPE.EQ.MAQ_I_PULSE1).OR.(TYPE.EQ.MAQ_I_PULSE2).OR.
     &  (TYPE.EQ.MAQ_TIME).OR.(TYPE.EQ.MAQ_COORD)) THEN

        IF(ACTION.EQ.MAQ_ALLOCATE.OR.ACTION.EQ.
     &    MAQ_ALLOCATE_AND_LOCK) THEN

          FOUND=.FALSE.
          maq=0
          DO WHILE(maq.LT.NMAQM.AND.(.NOT.FOUND))
            maq=maq+1
            IF((MAQ_TYPE(maq).EQ.0).OR.(MAQ_LOCKS(maq).EQ.MAQ_NOLOCK))
     &        FOUND=.TRUE.
          ENDDO
          IF(FOUND) THEN
            MAQ_SUBTYPE(maq)=SUBTYPE
            MAQ_TYPE(maq)=TYPE
            MAQ_LOCKS(maq)=LOCK
            MAQ_LIST(0)=MAQ_LIST(0)+1
            MAQ_LIST(MAQ_LIST(0))=maq
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Alloc. maq= '',I2,'' for type= '','
     &          //'I2,'', subtype= '',I2)') maq,TYPE,SUBTYPE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            maq=0
            ERROR='>>Increase NMAQM'
            GOTO 9999
          ENDIF

        ELSE IF(ACTION.EQ.MAQ_FREE) THEN

          FOUND=.FALSE.
          maqq=0
          DO WHILE(maqq.LT.MAQ_LIST(0).AND.(.NOT.FOUND))
            maqq=maqq+1
            maq=MAQ_LIST(maqq)
            IF((MAQ_TYPE(maq).EQ.TYPE).AND.MAQ_SUBTYPE(maq).EQ.SUBTYPE)
     &        FOUND=.TRUE.
          ENDDO
          IF(.NOT.FOUND) THEN
            WRITE(OP_STRING,'('' Cannot free maq, none allocated for'
     &        //' this type/subtype'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE
            MAQ_LIST(0)=MAQ_LIST(0)-1
            DO maqq=maq,MAQ_LIST(0)
              MAQ_LIST(maqq)=MAQ_LIST(maqq+1)
            ENDDO
            MAQ_LIST(MAQ_LIST(0)+1)=0
            MAQ_SUBTYPE(maq)=0
            MAQ_LOCKS(maq)=MAQ_NOLOCK
            MAQ_TYPE(maq)=0
          ENDIF

        ELSE IF(ACTION.EQ.MAQ_INQUIRE) THEN

          FOUND=.FALSE.
          maqq=0
          DO WHILE(maqq.LT.MAQ_LIST(0).AND.(.NOT.FOUND))
            maqq=maqq+1
            maq=MAQ_LIST(maqq)
            IF((MAQ_TYPE(maq).EQ.TYPE).AND.MAQ_SUBTYPE(maq).EQ.SUBTYPE)
     &        FOUND=.TRUE.
          ENDDO
          IF(.NOT.FOUND) maq=0
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Inquire: maq= '',I2,'' for type= '','
     &        //'I2,'', subtype= '',I2)') maq,TYPE,SUBTYPE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE
          maq=0
          ERROR='>>Invalid maq action'
          GOTO 9999
        ENDIF

      ELSE
        maq=0
        ERROR='>>Invalid maq type'
        GOTO 9999
      ENDIF

      CALL EXITS('MAQ_LOC')
      RETURN
 9999 CALL ERRORS('MAQ_LOC',ERROR)
      CALL EXITS('MAQ_LOC')
      RETURN 1
      END


