      SUBROUTINE OPGAUSG(ERROR,*)

C#### Subroutine: OPGAUSG
C###  Description:
C###    OPGAUSG outputs Gauss point groups.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER index,n,ne,ne_num_points,noelem,nogrga,nog,NOGT
!     Functions
      INTEGER ILISTMBR

      CALL ENTERS('OPGAUSG',*9999)

      WRITE(OP_STRING,'('' Number of Gauss point'
     '  //' groups = '',I2)') NTGRGA
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nogrga=1,NTGRGA
        NOGT=NLIGRGA(nogrga)
        WRITE(OP_STRING,'('' Group '',I2,'' Label= '',A,'' Total '
     '    //'points = '','
     '    //'I6)') nogrga,LAGRGA(nogrga),NOGT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        index=2
        DO noelem=1,ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),1)
          ne=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index)
          WRITE(OP_STRING,'('' Element '',I5)') ne
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          index=index+1
          ne_num_points=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index)
          DO nog=1,ne_num_points/10+1
            WRITE(OP_STRING,'((10(1X,I8)))')
     '        (ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),n),
     '        n=index+(nog-1)*10+1,index+MIN(nog*10,ne_num_points))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
          index=index+ne_num_points+1
        ENDDO
      ENDDO

      CALL EXITS('OPGAUSG')
      RETURN
 9999 CALL ERRORS('OPGAUSG',ERROR)
      CALL EXITS('OPGAUSG')
      RETURN 1
      END


