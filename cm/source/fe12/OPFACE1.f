      SUBROUTINE OPFACE1(NBJ,NBJF,nf,NLF,NNF,NPF,NPNE,
     '  NPNF,DF,SF,ERROR,*)

C#### Subroutine: OPFACE1
C###  Description:
C###     OPFACE1 outputs one face.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBJF(NJM),nf,NLF(4),NNF(0:17,6,NBFM),NPF(9),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM)
      REAL*8 DF,SF(NSM,NBFM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LOC_NBJF(9),NAF,nb,nbe,ne,nef,nj,nk,nn,nne
      CHARACTER NAME2(7)*10,NAME5(2)*8

      DATA NAME2(1)/'l.Lagrange'/,NAME2(2)/'q.Lagrange'/,
     '  NAME2(3)/'c.Lagrange'/,NAME2(4)/'c.Hermite '/,
     '  NAME2(5)/'          '/,NAME2(6)/'q1.Hermite'/,
     '  NAME2(7)/'q2.Hermite'/,
     '  NAME5 /'External','Internal'/

      CALL ENTERS('OPFACE1',*9999)

      CALL ASSERT(nf.LE.NFT,'>>Face greater than max',ERROR,*9999)

      IF(NIT(NBJ(1,1)).GT.1.AND.NJ_LOC(NJL_GEOM,0,1).EQ.3) THEN
        WRITE(OP_STRING,'(//3X,'' Global face segment data:'',//,6X,'
     '    //'''Face  Xi-dirs'',7X,'' Basis types'',6X,''Basis nos.'','
     '    //'7X,''Area'',6X,'' Face type'',4X,'' Elements'',4X,'
     '    //''' Nodes'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        nb=NBJF(1)
        ne=NPF(6)
        nef=NPF(8)
        nbe=NBJ(1,ne)
        DO nn=1,NNF(0,nef,nbe)
          nne=NNF(1+nn,nef,nbe)
          NPNF(nn,nb)=NPNE(nne,nbe,ne)
        ENDDO
        DO nj=1,9
          LOC_NBJF(nj)=0
        ENDDO
        DO nj=1,NJT
          LOC_NBJF(nj)=NBJF(nj)
        ENDDO
        WRITE(OP_STRING,'(6X,I4,5X,I1,1X,I1,6X,A,1X,A,2X,6I2,1X,'
     '    //'E11.4,3X,A,2X,2I5,2X,8I5)') nf,NPF(1),NPF(3),
     '    NAME2(NPF(2)),NAME2(NPF(4)),(LOC_NBJF(nj),nj=1,6),
     '    DF,NAME5(NPF(5)),NPF(6),NPF(7),(NPNF(nn,nb),nn=1,NNT(nb))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(NKT(0,nb).GT.1) THEN
          WRITE(OP_STRING,'(/6X,''  Face  '',2X,''Scale factors '','
     '      //'''for face no. '',I4,'' :  (Global line nos '',4I5,'
     '      //''')'',/7X,'' Node'',5X,''nk = '',7X,7(I1,10X))')
     '      nf,(NLF(NAF),NAF=1,4),(nk,nk=1,NKT(0,nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nn=1,NNT(nb)
            WRITE(OP_STRING,'(10X,I1,11X,7E11.3)') nn,
     '        (SF(nk+(nn-1)*NKT(0,nb),nb),nk=1,NKT(0,nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('OPFACE1')
      RETURN
 9999 CALL ERRORS('OPFACE1',ERROR)
      CALL EXITS('OPFACE1')
      RETURN 1
      END


