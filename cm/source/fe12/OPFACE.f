      SUBROUTINE OPFACE(NBJ,NBJF,NLF,NNF,NPF,NPNE,NPNF,DF,SF,ERROR,*)

C#### Subroutine: OPFACE
C###  Description:
C###    OPFACE outputs face data.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBJF(NJM,NFM),NLF(4,NFM),
     '  NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM)
      REAL*8 DF(NFM),SF(NSM,NBFM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LOC_NBJF(9),NAF,nb,nbe,ne,nef,nf,nj,nk,n,nn,nne,ns
      CHARACTER NAME2(7)*10,NAME5(2)*8

      DATA NAME2(1)/'l.Lagrange'/,NAME2(2)/'q.Lagrange'/,
     '  NAME2(3)/'c.Lagrange'/,NAME2(4)/'c.Hermite '/,
     '  NAME2(5)/'          '/,NAME2(6)/'q1.Hermite'/,
     '  NAME2(7)/'q2.Hermite'/,
     '  NAME5 /'External','Internal'/

      CALL ENTERS('OPFACE',*9999)
c     IF(NIT(NBJ(1,1)).GT.1.AND.NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
        WRITE(OP_STRING,'(//3X,'' Global face segment data : The '
     '    //'total no. of face segments = '',I5,/,6X,''Face  Xi-dirs'','
     '    //'7X,'' Basis types'',6X,''Basis nos.'',7X,''Area'',6X,'
     '    //''' Face type'',4X,'' Elements'',4X,'' Nodes'')') NFT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nf=1,NFT
          nb=NBJF(1,nf)
          ne=NPF(6,nf)
          nef=NPF(8,nf)
          nbe=NBJ(1,ne)
          DO nj=1,9
            LOC_NBJF(nj)=0
          ENDDO
          DO nj=1,NJT
            LOC_NBJF(nj)=NBJF(nj,nf)
          ENDDO
          DO nn=1,NNF(0,nef,nbe)
            nne=NNF(1+nn,nef,nbe)
            NPNF(nn,nb)=NPNE(nne,nbe,ne)
          ENDDO
          WRITE(OP_STRING,'(6X,I4,5X,I1,1X,I1,6X,A,1X,A,2X,6I2,1X,'
     '      //'D11.4,3X,A,2X,2I5,2X,8I5)') nf,NPF(1,nf),NPF(3,nf),
     '      NAME2(NPF(2,nf)),NAME2(NPF(4,nf)),
     '      (LOC_NBJF(nj),nj=1,6),
     '      DF(nf),NAME5(NPF(5,nf)),NPF(6,nf),
     '      NPF(7,nf),(NPNF(nn,nb),nn=1,NNT(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

        DO nf=1,NFT
          nb=NBJF(1,nf)
          IF(NKT(0,nb).GT.1) THEN
            WRITE(OP_STRING,'(/6X,''  Face  '',2X,''Scale factors '','
     '        //'''for face no. '',I4,'' :  (Global line nos '',4I5,'
     '        //''')'',/7X,'' Node'',5X,''nk = '',7X,7(I1,10X))')
     '        nf,(NLF(NAF,nf),NAF=1,4),(nk,nk=1,NKT(0,nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nn=1,NNT(nb)
              ns=0
              DO n=1,nn-1
                ns=ns+NKT(n,nb)
              ENDDO
              WRITE(OP_STRING,'(10X,I1,11X,7D11.3)') nn,
     '          (SF(nk+ns,nb),nk=1,NKT(nn,nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDDO
c     ENDIF

      CALL EXITS('OPFACE')
      RETURN
 9999 CALL ERRORS('OPFACE',ERROR)
      CALL EXITS('OPFACE')
      RETURN 1
      END


