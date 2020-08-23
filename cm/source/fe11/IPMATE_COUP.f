      SUBROUTINE IPMATE_COUP(nr,nx,ERROR,*)

C#### Subroutine: IPMATE_COUP
C###  Description:
C###    IPMAT3 inputs coupling parameters.

C DMAL CREATED: 08-AUG-2002

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'file01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'


!     Parameter List
      INTEGER nr,nx
!      REAL*8
      CHARACTER ERROR*(*)
!      LOGICAL
!     Local Variables
      INTEGER COUP_INSTANCE,COUP_TYPE,COUP_SUBTYPE,
     '  IBEG,IB1,IB2,IB3,IEND,IE1,IE2,IE3,ICHAR,INFO,IPTYPE,
     '  j,no_coup,NOQUES,QTYPE(0:10,2),ques_num,var_num
!      REAL*8
      CHARACTER CHAR,CHAR1*16,CHAR2*16,CHAR3*16,QUESTIONS(10)*64,
     ' QUESTION*64
!      LOGICAL
!     Functions
!      LOGICAL

      CALL ENTERS('IPMATE_COUP',*9999)

      ICHAR=999
      QTYPE(0,1)=0
      NOQUES=0
      var_num=1

      CALL ASSERT(ITYP2(nr,nx).NE.0,
     '  '>>Solution type has not been defined',ERROR,*9999)
      CALL ASSERT(ITYP2(nr,nx).GT.2,'>>Solution type is incorrect',
     '  ERROR,*9999)

C     Loop(1) through COUP_CL
      DO no_coup=1,COUP_CL(0,1,0) !L1

        COUP_TYPE=COUP_CL(no_coup,1,0)
        COUP_SUBTYPE=COUP_CL(no_coup,2,0)
        COUP_INSTANCE=COUP_CL(no_coup,3,0)

        WRITE(CHAR1,'(I2)') COUP_TYPE
        CALL STRING_TRIM(CHAR1,IB1,IE1)
        WRITE(CHAR2,'(I2)') COUP_SUBTYPE
        CALL STRING_TRIM(CHAR2,IB2,IE2)
        WRITE(CHAR3,'(I2)') COUP_INSTANCE
        CALL STRING_TRIM(CHAR3,IB3,IE3)
        FORMAT='(/'' Coupling Model - Type: '//CHAR1(IB1:IE1)//
     '    ', Subtype: '//CHAR2(IB2:IE2)//', Instance: '//
     '    CHAR3(IB3:IE3)//''')'

        IPTYPE=0
        CALL GINOUT(IOTYPE,IPTYPE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.EQ.1)THEN

        ELSEIF(IOTYPE.EQ.2)THEN
C          IF(CBBREV(CO,'DEPVAR',2,noco+1,PART2,N3CO)) THEN
C            nhc1=IFROMC(CO(N3CO+1))

        ENDIF

        FORMAT='(/$,'' Is this a boundary condition coupling model'//
     '    ' [N]: '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(COUP_CL(no_coup,4,0).EQ.1) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3)THEN
          IF(ADATA(1).EQ.'Y')THEN
            COUP_CL(no_coup,4,0)=1
          ELSE
            COUP_CL(no_coup,4,0)=0
          ENDIF
        ENDIF
        COUP_CL(no_coup,5,0)=0 ! number of variables in model

        CALL COUP_QUES(COUP_TYPE,COUP_SUBTYPE,QTYPE,
     '    QUESTIONS,ERROR,*9999) !Get questions
        COUP_CL(no_coup,5,0)=QTYPE(0,1)
        DO ques_num=1,QTYPE(0,1) !L2
          QUESTION=QUESTIONS(ques_num)
          CALL STRING_TRIM(QUESTION,IBEG,IEND)
          IDEFLT(1)=1
          IF(QTYPE(ques_num,1).EQ.0)THEN
            COUP_CL(no_coup,1,ques_num)=ques_num
            COUP_CL(no_coup,2,ques_num)=var_num
            COUP_CL(no_coup,3,ques_num)=QTYPE(ques_num,1)
            COUP_CL(no_coup,4,ques_num)=2 ! field
            FORMAT='(/$,'' The '//QUESTION(IBEG:IEND)//
     '        ' is stored in field [1] : '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=COUP_CL(no_coup,5,ques_num)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '        IDEFLT,1,NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,
     '        RDEFLT,RMIN,RMAX,1,ERROR,*9999)
            IF(IOTYPE.NE.3) COUP_CL(no_coup,5,ques_num)=IDATA(1)
          ELSE
            COUP_CL(no_coup,1,ques_num)=ques_num
            COUP_CL(no_coup,2,ques_num)=var_num
            COUP_CL(no_coup,3,ques_num)=QTYPE(ques_num,1)
            FORMAT='(/'' The '//QUESTION(IBEG:IEND)//' is [1] :'''//
     '        '/''   (1) Constant spatially                   '''//
     '        '/''   (2) Defined by a field (XP)              '''//
     '        '/''   (3) Defined by a dependent variable (YP) '''//
     '       '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=COUP_CL(no_coup,4,ques_num)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '        IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,1,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) COUP_CL(no_coup,4,ques_num)=IDATA(1)

            IF(COUP_CL(no_coup,4,ques_num).EQ.1)THEN ! constant
              RDEFLT(1)=0.0d0
              WRITE(CHAR1,'(E11.4)') RDEFLT(1)
              CALL STRING_TRIM(CHAR1,IB1,IE1)
              FORMAT='($,'' Enter the value of the constant ['//
     '          CHAR1(IB1:IE1)//']: '',E11.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=CCR(COUP_CL(no_coup,5,ques_num))
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '          -RMAX,RMAX,ques_num,ERROR,*9999)
              IF(IOTYPE.NE.3)THEN
                NCCR=NCCR+1
                COUP_CL(no_coup,5,ques_num)=NCCR
                CCR(NCCR)=RDATA(1)
              ENDIF
            ELSEIF(COUP_CL(no_coup,4,ques_num).EQ.2)THEN ! field
              WRITE(CHAR,'(I1)') IDEFLT(1)
              FORMAT='($,'' Enter the field # ['//CHAR//']: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=COUP_CL(no_coup,5,ques_num)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,1,NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,
     '          RDEFLT,RMIN,RMAX,1,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                COUP_CL(no_coup,5,ques_num)=IDATA(1)
              ENDIF
            ELSEIF(COUP_CL(no_coup,4,ques_num).EQ.3)THEN  ! dep. var.
              WRITE(CHAR,'(I1)') IDEFLT(1)
              FORMAT='($,'' Enter the region # ['//CHAR//']: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=COUP_CL(no_coup,5,ques_num)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,1,NH_LOC(0,nx),LDATA,LDEFLT,RDATA,RDEFLT,
     '          RMIN,RMAX,1,ERROR,*9999)
              IF(IOTYPE.NE.3) COUP_CL(no_coup,5,ques_num)=IDATA(1)
              FORMAT='($,'' Enter the class # ['//CHAR//']: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=COUP_CL(no_coup,6,ques_num)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,1,NH_LOC(0,nx),LDATA,LDEFLT,RDATA,RDEFLT,
     '          RMIN,RMAX,1,ERROR,*9999)
              IF(IOTYPE.NE.3) COUP_CL(no_coup,6,ques_num)=IDATA(1)
              FORMAT='($,'' Enter the dependent variable # ['//
     '          CHAR//']: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=COUP_CL(no_coup,7,ques_num)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,1,NH_LOC(0,nx),LDATA,LDEFLT,RDATA,RDEFLT,
     '          RMIN,RMAX,1,ERROR,*9999)
              IF(IOTYPE.NE.3) COUP_CL(no_coup,7,ques_num)=IDATA(1)
            ENDIF
          ENDIF
        ENDDO

C   Ask question
C   Save answers

C        ENDDO !L2
      IF(DOP)THEN
        write(*,*) 'COUPLING',no_coup
        DO j=1,QTYPE(0,1)
          write(*,*) COUP_CL(no_coup,1,j),COUP_CL(no_coup,2,j),
     '      COUP_CL(no_coup,3,j),COUP_CL(no_coup,4,j),
     '      COUP_CL(no_coup,5,j),COUP_CL(no_coup,6,j),
     '      COUP_CL(no_coup,7,j)
        ENDDO
      ENDIF

      ENDDO !L1

      CALL EXITS('IPMATE_COUP')
      RETURN
 9999 CALL ERRORS('IPMATE_COUP',ERROR)
      CALL EXITS('IPMATE_COUP')
      RETURN 1
      END


