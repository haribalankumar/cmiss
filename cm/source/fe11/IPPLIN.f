      SUBROUTINE IPPLIN(ERROR,*)

C#### Subroutine: IPPLIN
C###  Description:
C###    IPPLIN defines polyline parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'plin00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IBEG,ICHAR,IEND,INFO,INDEX_PLIN,IPFILE,nj,NL,no_plin,
     '  no_point,NOQUES,NT_POINT
      CHARACTER CHAR1*1,CHAR2*2
      LOGICAL FILEIP

      CALL ENTERS('IPPLIN',*9999)
      IPFILE=1 !is input file version number on 24-Jan-1990
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.3) THEN
        WRITE(UNIT=IFILE,REC=1,FMT='(A,I2)') 'CMISS Version '//CMISS
     '    //' IPPLIN File Version ',IPFILE
        WRITE(UNIT=IFILE,REC=2,FMT='(A)') 'Heading: '//HEADING
        WRITE(UNIT=IFILE,REC=3,FMT='(1X)')
      ELSE IF(IOTYPE.EQ.2.OR.IOTYPE.EQ.4) THEN
        READ(UNIT=IFILE,REC=1,FMT='(39X,I2)') IPFILE
        IF(DOP) THEN
          WRITE(OP_STRING,'('' File version number is '',I2)') IPFILE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ASSERT(IPFILE.EQ.1,
     '    'Old file version number - redo file input',ERROR,*9999)
        READ(UNIT=IFILE,REC=2,FMT='(9X,A)') HEADING
        WRITE(OP_STRING,'('' File heading: '',A)') HEADING
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        READ(UNIT=IFILE,REC=3,FMT='(1X)')
      ENDIF

      FORMAT='(/$,'' Enter number of polylines (<100)[1]: '',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=NT_PLIN
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,100,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NT_PLIN=IDATA(1)

      DO no_plin=1,NT_PLIN
        INDEX_PLIN=no_plin
        WRITE(CHAR2,'(I2)') INDEX_PLIN
        CALL STRING_TRIM(CHAR2,IBEG,IEND)
        FORMAT='('' For polyline '//CHAR2(IBEG:IEND)//
     '    ' enter type [1]:'''//
     '    '/''   (1) Piecewise linear '''//
     '    '/''   (2) Bezier curve     '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=INDEX_PLIN_TYPE(INDEX_PLIN)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) INDEX_PLIN_TYPE(INDEX_PLIN)=IDATA(1)

        FORMAT='($,'' Enter number of sections [1]: '',I2)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=NT_PLIN_SECTIONS(INDEX_PLIN)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NT_PLIN_SECTIONS(INDEX_PLIN)=IDATA(1)

        IF(INDEX_PLIN_TYPE(INDEX_PLIN).EQ.1) THEN !Piecewise linear polyline
          WRITE(CHAR1,'(I1)') NJT
          NT_POINT=NT_PLIN_SECTIONS(INDEX_PLIN)+1
          DO no_point=1,NT_POINT
            WRITE(CHAR2,'(I2)') no_point
            CALL STRING_TRIM(CHAR2,IBEG,IEND)
            FORMAT='($,'' Enter the '//CHAR1
     '        //' r.c. coords for point '
     '        //CHAR2(IBEG:IEND)//' [0]: '',D11.4)'
            IF(IOTYPE.EQ.3) THEN
              DO nj=1,NJT
                RDATA(nj)=PLIN_DATA(nj,no_point,INDEX_PLIN)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nj=1,NJT
                PLIN_DATA(nj,no_point,INDEX_PLIN)=RDATA(nj)
              ENDDO
            ENDIF
          ENDDO

        ELSE IF(INDEX_PLIN_TYPE(INDEX_PLIN).EQ.2) THEN !Bezier polyline
          WRITE(CHAR1,'(I1)') NJT
          DO nl=1,NT_PLIN_SECTIONS(INDEX_PLIN)
            DO I=1,4
              FORMAT='($,'' Enter the '//CHAR1
     '          //' coords for control point '
     '          //CHAR2(1:1)//' [0]: '',3D11.4)'
              IF(IOTYPE.EQ.3) THEN
                DO nj=1,NJT
                  RDATA(nj)=PLIN_DATA(nj,3*(nl-1)+I,INDEX_PLIN)
                ENDDO
              ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          0,IMAX,LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                DO nj=1,NJT
                  PLIN_DATA(nj,3*(nl-1)+I,INDEX_PLIN)=RDATA(nj)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      CALL EXITS('IPPLIN')
      RETURN
 9999 CALL ERRORS('IPPLIN',ERROR)
      CALL EXITS('IPPLIN')
      RETURN 1
      END


