      SUBROUTINE DERIV_INFO(IDO,NBJ,NEELEM,NENP,NKJE,NKJ,NPNE,
     '  NPNODE,XP,ERROR,*)

C#### Subroutine: DERIV_INFO
C###  Description:
C###    DERIV_INFO adjusts local derivative numbers nkje(nk,nn,nj,ne)
C###    according to the interpolation used in each element sharing
C###    node np.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IRESID,nb,NBNK,ne,NENK,ni,nj,NJ2,nk,nk2,nn,
     '  noelem,nonode,np,nr
      LOGICAL ALL_DERIVATIVES,FOUND

      CALL ENTERS('DERIV_INFO',*9999)

C!!! KAT: This default value does not look good.
C!!! This element may have nothing to do with the node in question.
C!!! Don't we need a default for each node?
      NENK=NEELEM(1,1) !Default value
      DO nr=1,NRT !Added AJP 12-10-92 (Previously nr=1)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nj=1,NJT
            DO noelem=1,NEELEM(0,nr) !Loop over elements
              ne=NEELEM(noelem,nr)
              nb=NBJ(nj,ne)
              IF(nb.NE.0) THEN
                DO nn=1,NNT(nb) !Loop over element nodes
                  IF(np.EQ.NPNE(nn,nb,ne))THEN !Node is in element ne
                    IF(NENP(np,0,nr).EQ.0)THEN
                      NENK=ne
                    ELSE
                      IF(NKT(nn,nb).GT.NKT(0,NBJ(nj,NENK))) NENK=ne
                      !Element sharing np with largest NKT(0,nb)
                    ENDIF
                  ENDIF
                ENDDO !end nn
              ENDIF
            ENDDO !end ne
            DO noelem=1,NENP(np,0,nr)
              ne=NENP(np,noelem,nr)
              nb=NBJ(nj,ne)
              IF(nb.NE.0.AND.NKT(0,nb).NE.NKJ(nj,np))THEN
!               !Adjust derivative numbers
!               !AJP 16-5-93
!               !Check all geometric derivatives at node np.  If they
!               !are all zero then node is part of a distorted element
!               !so set nkje(nk,nn,nj,ne) to be zero for these nk values
                IF(NBC(nb).EQ.5)THEN !BE basis
                  ALL_DERIVATIVES=.FALSE.
                  NJ2=0
                  DO WHILE((NJ2.LT.NJT).AND.(.NOT.ALL_DERIVATIVES))
                    NJ2=NJ2+1
                    nk=1
                    DO WHILE((nk.LT.NKJ(nj2,np)).AND.
     '                (.NOT.ALL_DERIVATIVES))
                      nk=nk+1
                      IF(DABS(XP(nk,1,NJ2,np)).GT.RDELTA)
     '                  ALL_DERIVATIVES=.TRUE.
                    ENDDO
                  ENDDO
                ELSE
                  ALL_DERIVATIVES=.TRUE.
                ENDIF !END AJP 16-5-93

                NBNK=NBJ(nj,NENK) !nb with largest NKT
                nn=1
                DO WHILE(np.NE.NPNE(nn,nb,ne))
                  nn=nn+1  !Find local node nn associated with np and ne
                  IF(nn.GT.NNT(nb))THEN
                    WRITE(OP_STRING,
     '                '('' >>Error in DERIV_INFO - nn too large'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ENDDO
                DO nk=1,NKT(0,nb)
                  nk2=0
                  FOUND=.FALSE.
                  DO WHILE(.NOT.FOUND)
                    nk2=nk2+1
                    IF(nk2.GT.NKT(0,NBNK))THEN
                      WRITE(OP_STRING,*)
     '                  ' Error in DERIV_INFO - NK2 too large'
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      GOTO 9999
                    ENDIF
                    IRESID=0
                    DO ni=1,NIT(nb)
C!!! KAT: nn for NBNK should be the nn for np.
                      IRESID=ABS(IDO(nk2,nn,ni,NBNK)-IDO(nk,nn,ni,nb))+
     '                  IRESID
                    ENDDO
                    IF(IRESID.EQ.0)FOUND=.TRUE.
                  ENDDO
                  IF(ALL_DERIVATIVES.OR.(nk.EQ.1))THEN  !AJP 16-5-93
                    NKJE(nk,nn,nj,ne)=nk2
                  ELSE
                    NKJE(nk,nn,nj,ne)=0
                  ENDIF
                ENDDO !End of loop over nk
              ENDIF
            ENDDO !End of loop over elements sharing np
          ENDDO !End of loop over nj
        ENDDO !End of loop over np
      ENDDO !End of loop over nr

      CALL EXITS('DERIV_INFO')
      RETURN
 9999 CALL ERRORS('DERIV_INFO',ERROR)
      CALL EXITS('DERIV_INFO')
      RETURN 1
      END


