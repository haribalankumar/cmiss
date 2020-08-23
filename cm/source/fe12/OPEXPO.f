      SUBROUTINE OPEXPO(ERROR,*)

C#### Subroutine: OPEXPO
C###  Description:
C###    OPEXPO outputs export parameters.


C**** LKC 4-SEP-98 Adding export of signal traces

C**** NT_SUB_DIV is the number of local subdivisions given to each
C**** element.

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmgui01.cmn'
      INCLUDE 'emap00.cmn'
      INCLUDE 'expo00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ncol,nonr,nr,nrow
      CHARACTER RIG_TYPES(4)*10

C LKC 3-OCT-98 These are not consistent with those in
C      emap00.cmn which leads to problems
C      PARAMETER(EMAP_MIXED=-1,EMAP_SOCK=0,EMAP_PATCH=1,EMAP_TORSO=2)
C
C      DATA RIG_TYPES / 'SOCK','PATCH','TORSO','MIXED' /

      DATA RIG_TYPES / 'MIXED','SOCK','PATCH','TORSO' /

      CALL ENTERS('OPEXPO',*9999)
      CALL ASSERT(CALL_EXPO,
     '  '>>Export parameters not defined',ERROR,*9999)

      IF(DEFEXPO_TYPE.EQ.1) THEN !Signal
        IF(SIGEXPO_TYPE.EQ.1) THEN !UNEMAP
          WRITE(OP_STRING,'(/'' UNEMAP export signal paramters'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C LKC 3-OCT-98 alter DATA RIG_TYPES
C          WRITE(OP_STRING,'(/'' The main rig type is '',A)')
C     '      RIG_TYPES(EMAP_RIGTYPE(0)+1)
          WRITE(OP_STRING,'(/'' The main rig type is '',A)')
     '      RIG_TYPES(EMAP_RIGTYPE(0)+2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(EMAP_RIGTYPE(0).EQ.1) THEN
            WRITE(OP_STRING,'('' The sock focus is '',D12.4)')
     '        EMAP_SOCKFOCUS(1)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' The rig name is '',A)') EMAP_RIGNAME
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The number of regions is '',I1)')
     '      EMAP_NUMREGIONS
          DO nonr=1,EMAP_NUMREGIONS


C LKC 7-SEP-1999 adding exporting a region list
            nr=EXPORT_NRLIST(nonr)


            WRITE(OP_STRING,'(/'' Region '',I1)') nr
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(EMAP_RIGTYPE(0).EQ.4) THEN
              WRITE(OP_STRING,'('' The region rig type is '',A)')
     '          RIG_TYPES(EMAP_RIGTYPE(nr))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(EMAP_RIGTYPE(0).EQ.1) THEN
                WRITE(OP_STRING,'('' The sock focus is '',D12.4)')
     '            EMAP_SOCKFOCUS(nr)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            WRITE(OP_STRING,'('' The region name is '',A)')
     '        EMAP_REGIONNAME(nr)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nr

C LKC 4-SEPT-98 not used
C
C        ELSE IF(SIGEXPO_TYPE.EQ.2) THEN !Map3d
C          WRITE(OP_STRING,'(/'' MAP3D export geometry paramters'')')
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C          WRITE(OP_STRING,'(/'' Local subdivisions per element = '','
C     '      //'I2)') NT_SUB_DIV
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(SIGEXPO_TYPE.EQ.2) THEN !CMGUI
          WRITE(OP_STRING,'(/'' Exporting to CMGUI (via exgobj)'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(OP_STRING(1),
     '      '('' Signal  Name     : '',A)') NAME_EXGOBJ(1)
          WRITE(OP_STRING(2),
     '      '('' Signal Material  : '',A)') MATERIAL_EXGOBJ(1)
          WRITE(OP_STRING(3),
     '      '('' Pointer Name     : '',A)') NAME_EXGOBJ(2)
          WRITE(OP_STRING(4),
     '      '('' Pointer Material : '',A)') MATERIAL_EXGOBJ(2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)


          WRITE(OP_STRING,'(/'' Transformation Matrix'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          DO nrow=1,4
            WRITE(IOFI,'(4F12.3)')
     '        (TRANSFORM_EXGOBJ(nrow,ncol),ncol=1,4)
          ENDDO
          OP_STRING(1)=' '
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE
          ERROR='Not implemented'
        ENDIF
      ENDIF

      CALL EXITS('OPEXPO')
      RETURN
 9999 CALL ERRORS('OPEXPO',ERROR)
      CALL EXITS('OPEXPO')
      RETURN 1
      END

