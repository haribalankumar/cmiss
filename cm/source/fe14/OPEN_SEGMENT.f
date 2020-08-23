      SUBROUTINE OPEN_SEGMENT(ISEGNUM,ISEG,iw,CLABEL,index,INDEX_OLD,
     '  NLABEL,IVIS,CSEG,ERROR,*)

C#### Subroutine: OPEN_SEGMENT
C###  Description:
C###    OPEN_SEGMENT opens graphics segment ISEGNUM. If ISEGNUM does
C###    not already exist (has the value 0), then NTSG is incremented
C###    and a new segment is created.  Otherwise the old segment
C###    structure is replaced by the new structure.  All subsequent
C###    call to graphics primitives, etc will be part of this
C###    structure. Note that GKS segments cannot be nested, so a call to
C###    OPEN_SEGMENT must be followed by a call to CLOSE_SEGMENT before
C###    the next OPEN_SEGMENT.

C**** When created CSEG(ISEGNUM) stores a string comprising the first 48
C**** characters of CLABEL together with index(4 chars) and
C**** NLABEL(5 chars) then a '/' and the iw number(2 chars).
C**** index is the bundle table index for the current segment and if
C**** an old segment is being replaced, INDEX_OLD is read from the
C**** CSEG string and returned for possible use in the replacement
C**** segment creation.
C**** IVIS=1 creates a visible segment
C**** IVIS=2 creates an invisible segment

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER index,INDEX_OLD,ISEG(*),ISEGNUM,IVIS,iw,NLABEL
      CHARACTER CLABEL*(*),CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER ERR,IFROMC
      CHARACTER CHAR2*2,CHAR4*4,CHAR5*5,CLABEL2*(52)

      CALL ENTERS('OPEN_SEGMENT',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'' IWKT(iw)='',I2,'
     '    //''' IWKG(iw)='',I2)') iw,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISEGNUM='',I5)') ISEGNUM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(ISEGNUM.EQ.0) THEN !segment does not already exist
        NTSG=NTSG+1
        ISEGNUM=NTSG
      ELSE !recover index, delete old segment and use the old ISEG,CSEG labels
        INDEX_OLD=IFROMC(CSEG(ISEGNUM)(49:52))
        IF(DOP) THEN
          WRITE(OP_STRING,'('' INDEX='',I3)') index
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(IWKT(iw).EQ.1) THEN      !GKS
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' Delete segment number '',I3)') ISEGNUM
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
            CALL DLOBJT(ISEGNUM,ERR)
            IF(err.ne.0) THEN
              ERROR='>>Error from GX'
              GOTO 9999
            ENDIF
          ENDIF
        ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        ENDIF
      ENDIF

      ISEG(ISEGNUM)=2
      IF(IVIS.EQ.2) ISEG(ISEGNUM)=1

C     open structure for subsequent graphic, set insert mode
      IF(IWKT(iw).EQ.1) THEN      !GKS
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Create segment isegnum='',I5,'' on iw='',I2)')
     '      ISEGNUM,iw
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL OPEN_SEGMENT_GX(ISEGNUM,ISEG,ERROR,*9999)
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
      ENDIF

      WRITE(CHAR2,'(I2)') iw
      WRITE(CHAR4,'(I4)') index
      WRITE(CHAR5,'(I5)') NLABEL
      CLABEL2(1:)=CLABEL
      CSEG(ISEGNUM)=CLABEL2(1:48)//CHAR4(1:4)//CHAR5(1:5)//'/'//
     '  CHAR2(1:2)

      CALL EXITS('OPEN_SEGMENT')
      RETURN
 9999 CALL ERRORS('OPEN_SEGMENT', ERROR)
      CALL EXITS('OPEN_SEGMENT')
      RETURN 1
      END


