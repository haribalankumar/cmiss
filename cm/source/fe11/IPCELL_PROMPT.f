      SUBROUTINE IPCELL_PROMPT(IBT,IDO,INP,NEELEM,NELIST,NENQ,NPNE,
     '  NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,CP,CQ,XE,ERROR,*)

C#### Subroutine: IPCELL_PROMPT
C###  Description:
C###    IPCELL_PROMPT prompts for cellular parameters for simple
C###    membrane models.

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b32.cmn'
      INCLUDE 'file01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'titl30.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENQ(0:8,NQM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),nr,nx
      REAL*8 CE(NMM,NEM),CP(NMM,NPM),CQ(NMM,NQM),XE(NSM,NJM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IB,IB5,IBEG,ICHAR,IE,IE5,IEND,II,IJ,IK,il,
     '  ILT_LOCAL(0:1),INFO,nb,NBMIN,ne,neq,nii,nij,nik,nn,noelem,
     '  nonode,NOQUES,np,nq,nqq,nqsc,ns
      REAL*8 PXI,XI(3)
      CHARACTER CHAR*100,CHAR1*1,CHAR11*11,CHAR5*5

      CALL ENTERS('IPCELL_PROMPT',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      ILT_LOCAL(0)=9
      IF(ITYP3(nr,nx).EQ.1) THEN
        ILT_LOCAL(1)=12
      ELSE IF(ITYP3(nr,nx).EQ.2) THEN
        ILT_LOCAL(1)=18
      ELSE IF(ITYP3(nr,nx).EQ.3) THEN
        ILT_LOCAL(1)=14
      ELSE IF(ITYP3(nr,nx).EQ.4) THEN
        ILT_LOCAL(1)=13
      ENDIF
      CALL ASSERT(ILT_LOCAL(1).LE.NMM,'>>Increase NMM',ERROR,*9999)

      DO il=1,2
        CHAR=TITL32(il,ITYP2(nr,nx),ITYP3(nr,nx))
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        IDEFLT(1)=1
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='(/'' '//CHAR(IBEG:IEND)//' is ['//CHAR1//']: '''//
     '    '/''   (1) Constant spatially                      '''//
     '    '/''   (2) Piecewise constant (defined by elements)'''//
     '    '/''   (3) Piecewise linear   (defined by nodes)   '''//
     '    '/''   (4) Defined by grid points (CQ)             '''//
     '    '/''   (5) Defined elsewhere by Gauss pt array (YG)'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=ILP(il,1,nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,il,ERROR,*9999)
        IF(IOTYPE.NE.3) ILP(il,1,nr,nx)=IDATA(1)

        IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
          RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),ITYP3(nr,nx))
          WRITE(CHAR11,'(E11.4)') RDEFLT(1)
          CALL STRING_TRIM(CHAR11,IB,IE)
          FORMAT='($,'' The value is ['//CHAR11(IB:IE)//']: '',E11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=CQ(il,1)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nq=NQR(1,nr),NQR(2,nr)
              CQ(il,nq)=RDATA(1)
            ENDDO
          ENDIF

        ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(ne.EQ.1) RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),
     '        ITYP2(nr,nx),ITYP3(nr,nx))
            IF(ne.GT.1) RDEFLT(1)=CE(il,ne-1)
            WRITE(CHAR5,'(I5)') ne
            CALL STRING_TRIM(CHAR5,IB5,IE5)
            WRITE(CHAR11,'(E11.4)') RDEFLT(1)
            CALL STRING_TRIM(CHAR11,IB,IE)
            FORMAT='($,'' The value in element '//CHAR5(IB5:IE5)//
     '        ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=CE(il,ne)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              CE(il,ne)=RDATA(1)
              DO nqq=1,NQET(NQS(ne))
                nq=NQNE(ne,nqq)
                CQ(il,nq)=RDATA(1)
              ENDDO
            ENDIF
          ENDDO

        ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
          FORMAT='($,'' Enter basis type number [1]: '',I1)'
          NBMIN=1
          IF(IOTYPE.EQ.3) IDATA(1)=NMB(il,1,nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,NBMIN,NBFM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
          IF(IOTYPE.NE.3) NMB(il,1,nx)=IDATA(1)

          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            IF(nonode.EQ.1) RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),
     '        ITYP2(nr,nx),ITYP3(nr,nx))
            IF(nonode.GT.1) RDEFLT(1)=CP(il,
     '        NPNODE(nonode-1,nr))
            WRITE(CHAR5,'(I5)') np
            CALL STRING_TRIM(CHAR5,IB5,IE5)
            WRITE(CHAR11,'(E11.4)') RDEFLT(1)
            CALL STRING_TRIM(CHAR11,IB,IE)
            FORMAT='($,'' The value at node '//CHAR5(IB5:IE5)//
     '        ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=CP(il,np)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) CP(il,np)=RDATA(1)
          ENDDO

          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NMB(il,1,nx)
            ns=0
            DO nn=1,NNT(nb)
              ns=ns+1
C Assuming linear basis!!!! (as above i.e. doesn't ask for derivatives)
              XE(ns,1)=CP(il,NPNE(nn,nb,ne))
            ENDDO !nn
            nqsc=NQS(ne)
            II=MAX(1,NQXI(1,nqsc))
            IJ=1
            IK=1
            IF(NQXI(0,nqsc).GT.1) IJ=MAX(1,NQXI(2,nqsc))
            IF(NQXI(0,nqsc).GT.2) IK=MAX(1,NQXI(3,nqsc))
            DO i=1,3
              XI(i)=0.0d0
            ENDDO !i
            DO nik=1,IK
              DO nij=1,IJ
                DO nii=1,II
                  neq=nii+((nij-1)*NQXI(1,nqsc))
                  IF(NQXI(0,nqsc).GT.1) neq=neq+((nik-1)*
     '              NQXI(1,nqsc)*NQXI(2,nqsc))
                  nq=NQNE(ne,neq)
                  IF(NENQ(1,nq).EQ.ne) THEN
                    IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                    IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                    IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)
                    CQ(il,nq)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,XI,XE(1,1))
                  ENDIF
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
          ENDDO !noelem

        ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by grid points
          CALL ASSERT(NQT.GT.0,'>>No grid points defined',ERROR,*9999)

          IF(IOTYPE.NE.3) THEN
            FORMAT='($,'' Enter grid values by element [Y]? '',A)'
            ADEFLT(1)='Y'
            ADATA(1)='Y'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

            IF(ADATA(1).EQ.'Y') THEN
 802          FORMAT='($,'' Enter element #/s [EXIT]: '',I5)'
              CDATA(1)='ELEMENTS' !for use with group input
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IZERO,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '          INFO,ERROR,*9999)

              IF(IDATA(1).NE.0) THEN !not default exit
                NELIST(0)=IDATA(0)
                DO ne=1,NELIST(0)
                  NELIST(ne)=IDATA(ne)
                ENDDO

                FORMAT='($,'' Is the value constant within '
     '            //'element(s) [Y]? '',A)'
                ADEFLT(1)='Y'
                ADATA(1)='Y'
                CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)

                IF(ADATA(1).EQ.'Y') THEN
                  RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),
     '              ITYP3(nr,nx))
                  WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                  CALL STRING_TRIM(CHAR11,IB,IE)
                  FORMAT='($,'' The value in the selected elements'
     '              //' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     &              RMAX,INFO,ERROR,*9999)

                  DO noelem=1,NELIST(0)
                    ne=NELIST(noelem)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      CQ(il,nq)=RDATA(1)
                    ENDDO
                  ENDDO
                  GO TO 802
                ELSE
                  DO noelem=1,NELIST(0)
                    ne=NELIST(noelem)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      IF(nqq.EQ.1) THEN
                        RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),
     '                    ITYP3(nr,nx))
                      ELSE
                        RDEFLT(1)=CQ(il,NQNE(ne,nqq-1))
                      ENDIF
                      WRITE(CHAR5,'(I5)') nq
                      CALL STRING_TRIM(CHAR5,IB5,IE5)
                      WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                      CALL STRING_TRIM(CHAR11,IB,IE)
                      FORMAT='($,'' The value for grid '
     '                  //CHAR5(IB5:IE5)//
     '                  ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                  IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &                  -RMAX,RMAX,INFO,ERROR,*9999)
                      CQ(il,nq)=RDATA(1)
                    ENDDO
                  ENDDO
                  GOTO 802
                ENDIF
              ENDIF
            ELSE
              DO nq=1,NQT
                IF(nq.EQ.1) THEN
                  RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),
     '              ITYP3(nr,nx))
                ELSE
                  RDEFLT(1)=CQ(il,nq-1)
                ENDIF
                WRITE(CHAR5,'(I5)') nq
                CALL STRING_TRIM(CHAR5,IB5,IE5)
                WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                CALL STRING_TRIM(CHAR11,IB,IE)
                FORMAT='($,'' The value for grid '//CHAR5(IB5:IE5)//
     '            ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                CQ(il,nq)=RDATA(1)
              ENDDO
            ENDIF
          ELSE IF(IOTYPE.EQ.3) THEN
            FORMAT='($,'' Enter grid values by element [Y]? '',A)'
            ADEFLT(1)='N'
            ADATA(1)='N'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            DO nq=1,NQT
              RDEFLT(1)=CQ(il,nq)
              WRITE(CHAR5,'(I5)') nq
              CALL STRING_TRIM(CHAR5,IB5,IE5)
              WRITE(CHAR11,'(E11.4)') RDEFLT(1)
              CALL STRING_TRIM(CHAR11,IB,IE)
              FORMAT='($,'' The value for grid '//CHAR5(IB5:IE5)//
     '          ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
              RDATA(1)=CQ(il,nq)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '          RMAX,INFO,ERROR,*9999)
            ENDDO
          ENDIF

        ELSE IF(ILP(il,1,nr,nx).EQ.5) THEN
          !defined elsewhere by elements

        ENDIF
      ENDDO !il

      DO il=ILT_LOCAL(0),ILT_LOCAL(1)
        CHAR=TITL32(il,ITYP2(nr,nx),ITYP3(nr,nx))
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        IDEFLT(1)=1
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='(/'' '//CHAR(IBEG:IEND)//' is ['//CHAR1//']: '''//
     '    '/''   (1) Constant spatially                      '''//
     '    '/''   (2) Piecewise constant (defined by elements)'''//
     '    '/''   (3) Piecewise linear   (defined by nodes)   '''//
     '    '/''   (4) Defined by grid points (CQ)             '''//
     '    '/''   (5) Defined elsewhere by Gauss pt array (YG)'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=ILP(il,1,nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,il,ERROR,*9999)
        IF(IOTYPE.NE.3) ILP(il,1,nr,nx)=IDATA(1)

        IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
          RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),ITYP3(nr,nx))
          WRITE(CHAR11,'(E11.4)') RDEFLT(1)
          CALL STRING_TRIM(CHAR11,IB,IE)
          FORMAT='($,'' The value is ['//CHAR11(IB:IE)//']: '',E11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=CQ(il,1)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nq=NQR(1,nr),NQR(2,nr)
              CQ(il,nq)=RDATA(1)
            ENDDO
          ENDIF

        ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(ne.EQ.1) RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),
     '        ITYP2(nr,nx),ITYP3(nr,nx))
            IF(ne.GT.1) RDEFLT(1)=CE(il,ne-1)
            WRITE(CHAR5,'(I5)') ne
            CALL STRING_TRIM(CHAR5,IB5,IE5)
            WRITE(CHAR11,'(E11.4)') RDEFLT(1)
            CALL STRING_TRIM(CHAR11,IB,IE)
            FORMAT='($,'' The value in element '//CHAR5(IB5:IE5)//
     '        ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=CE(il,ne)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              CE(il,ne)=RDATA(1)
              DO nqq=1,NQET(NQS(ne))
                nq=NQNE(ne,nqq)
                CQ(il,nq)=RDATA(1)
              ENDDO
            ENDIF
          ENDDO

        ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
          FORMAT='($,'' Enter basis type number [1]: '',I1)'
          NBMIN=1
          IF(IOTYPE.EQ.3) IDATA(1)=NMB(il,1,nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,NBMIN,NBFM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
          IF(IOTYPE.NE.3) NMB(il,1,nx)=IDATA(1)

          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            IF(nonode.EQ.1) RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),
     '        ITYP2(nr,nx),ITYP3(nr,nx))
            IF(nonode.GT.1) RDEFLT(1)=CP(il,
     '        NPNODE(nonode-1,nr))
            WRITE(CHAR5,'(I5)') np
            CALL STRING_TRIM(CHAR5,IB5,IE5)
            WRITE(CHAR11,'(E11.4)') RDEFLT(1)
            CALL STRING_TRIM(CHAR11,IB,IE)
            FORMAT='($,'' The value at node '//CHAR5(IB5:IE5)//
     '        ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=CP(il,np)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) CP(il,np)=RDATA(1)
          ENDDO

          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NMB(il,1,nx)
            ns=0
            DO nn=1,NNT(nb)
              ns=ns+1
C Assuming linear basis!!!! (as above i.e. doesn't ask for derivatives)
              XE(ns,1)=CP(il,NPNE(nn,nb,ne))
            ENDDO !nn
            nqsc=NQS(ne)
            II=MAX(1,NQXI(1,nqsc))
            IJ=1
            IK=1
            IF(NQXI(0,nqsc).GT.1) IJ=MAX(1,NQXI(2,nqsc))
            IF(NQXI(0,nqsc).GT.2) IK=MAX(1,NQXI(3,nqsc))
            DO i=1,3
              XI(i)=0.0d0
            ENDDO !i
            DO nik=1,IK
              DO nij=1,IJ
                DO nii=1,II
                  neq=nii+((nij-1)*NQXI(1,nqsc))
                  IF(NQXI(0,nqsc).GT.1) neq=neq+((nik-1)*
     '              NQXI(1,nqsc)*NQXI(2,nqsc))
                  nq=NQNE(ne,neq)
                  IF(NENQ(1,nq).EQ.ne) THEN
                    IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                    IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                    IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)
                    CQ(il,nq)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,XI,XE(1,1))
                  ENDIF
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
          ENDDO !noelem

        ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by grid points
          CALL ASSERT(NQT.GT.0,'>>No grid points defined',ERROR,*9999)

          IF(IOTYPE.NE.3) THEN
            FORMAT='($,'' Enter grid values by element [Y]? '',A)'
            ADEFLT(1)='Y'
            ADATA(1)='Y'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

            IF(ADATA(1).EQ.'Y') THEN
 804          FORMAT='($,'' Enter element #/s [EXIT]: '',I5)'
              CDATA(1)='ELEMENTS' !for use with group input
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IZERO,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '          INFO,ERROR,*9999)

              IF(IDATA(1).NE.0) THEN !not default exit
                NELIST(0)=IDATA(0)
                DO ne=1,NELIST(0)
                  NELIST(ne)=IDATA(ne)
                ENDDO

                FORMAT='($,'' Is the value constant within '
     '            //'element(s) [Y]? '',A)'
                ADEFLT(1)='Y'
                ADATA(1)='Y'
                CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)

                IF(ADATA(1).EQ.'Y') THEN
                  RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),
     '              ITYP3(nr,nx))
                  WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                  CALL STRING_TRIM(CHAR11,IB,IE)
                  FORMAT='($,'' The value in the selected elements'
     '              //' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     &              RMAX,INFO,ERROR,*9999)

                  DO noelem=1,NELIST(0)
                    ne=NELIST(noelem)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      CQ(il,nq)=RDATA(1)
                    ENDDO
                  ENDDO
                  GO TO 804
                ELSE
                  DO noelem=1,NELIST(0)
                    ne=NELIST(noelem)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      IF(nqq.EQ.1) THEN
                        RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),
     '                    ITYP3(nr,nx))
                      ELSE
                        RDEFLT(1)=CQ(il,NQNE(ne,nqq-1))
                      ENDIF
                      WRITE(CHAR5,'(I5)') nq
                      CALL STRING_TRIM(CHAR5,IB5,IE5)
                      WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                      CALL STRING_TRIM(CHAR11,IB,IE)
                      FORMAT='($,'' The value for grid '
     '                  //CHAR5(IB5:IE5)//
     '                  ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                  IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &                  -RMAX,RMAX,INFO,ERROR,*9999)
                      CQ(il,nq)=RDATA(1)
                    ENDDO
                  ENDDO
                  GOTO 804
                ENDIF
              ENDIF
            ELSE
              DO nq=1,NQT
                IF(nq.EQ.1) THEN
                  RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),
     '              ITYP3(nr,nx))
                ELSE
                  RDEFLT(1)=CQ(il,nq-1)
                ENDIF
                WRITE(CHAR5,'(I5)') nq
                CALL STRING_TRIM(CHAR5,IB5,IE5)
                WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                CALL STRING_TRIM(CHAR11,IB,IE)
                FORMAT='($,'' The value for grid '//CHAR5(IB5:IE5)//
     '            ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                CQ(il,nq)=RDATA(1)
              ENDDO
            ENDIF
          ELSE IF(IOTYPE.EQ.3) THEN
            FORMAT='($,'' Enter grid values by element [Y]? '',A)'
            ADEFLT(1)='N'
            ADATA(1)='N'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            DO nq=1,NQT
              RDEFLT(1)=CQ(il,nq)
              WRITE(CHAR5,'(I5)') nq
              CALL STRING_TRIM(CHAR5,IB5,IE5)
              WRITE(CHAR11,'(E11.4)') RDEFLT(1)
              CALL STRING_TRIM(CHAR11,IB,IE)
              FORMAT='($,'' The value for grid '//CHAR5(IB5:IE5)//
     '          ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
              RDATA(1)=CQ(il,nq)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '          RMAX,INFO,ERROR,*9999)
            ENDDO
          ENDIF

        ELSE IF(ILP(il,1,nr,nx).EQ.5) THEN
          !defined elsewhere by elements

        ENDIF
      ENDDO !il

      CALL EXITS('IPCELL_PROMPT')
      RETURN
 9999 CALL ERRORS('IPCELL_PROMPT',ERROR)
      CALL EXITS('IPCELL_PROMPT')
      RETURN 1
      END


