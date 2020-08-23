      SUBROUTINE CALC_NUNK(IDO,NBJ,NENP,NKJE,NKJ,NPNE,NPNODE,
     '  NRLIST,NUNK,ERROR,*)

C#### Subroutine: CALC_NUNK
C###  Description:
C###    CALC_NUNK calculates the mapping array NUNK

C**** Created by Carey Stevens September 1998

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'parameters.inc'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NUNK(NKM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NUIDO(2,2,2),nb,ne,nj,njj,njtype,nk,nkp,nn,noelem,nonode,
     '  no_nrlist,np,nr

      DATA NUIDO/1,2,4,6,7,9,10,11/

      CALL ENTERS('CALC_NUNK',*9999)

      DO no_nrlist=1,NRLIST(0)
        nr=NRLIST(no_nrlist)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO njtype=1,3 !NJL_GEOM,NJL_FIBR & NJL_FIEL
            DO njj=1,NJ_LOC(njtype,0,nr)
              nj=NJ_LOC(njtype,njj,nr)
              DO nk=1,NKJ(nj,np)
                NUNK(nk,nj,np)=0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO no_nrlist=1,NRLIST(0)
        nr=NRLIST(no_nrlist)

        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)

          DO njtype=1,3 !NJL_GEOM,NJL_FIBR & NJL_FIEL
            DO njj=1,NJ_LOC(njtype,0,nr)
              nj=NJ_LOC(njtype,njj,nr)

              noelem=1
              DO WHILE(noelem.LE.NENP(np,0,nr))
                ne=NENP(np,noelem,nr)
                nb=NBJ(nj,ne)
                IF(nb.NE.0) THEN
                  nn=1
                  DO WHILE ((nn.LE.NNT(nb)).AND.
     '              (NPNE(nn,nb,ne).NE.np))
                    nn=nn+1
                  ENDDO

                  DO nk=1,NKT(nn,nb)
                    nkp=NKJE(nk,nn,nj,ne)

                    IF(NUNK(nkp,nj,np).EQ.0) THEN
                      IF(NIT(nb).EQ.1) THEN
                        NUNK(nkp,nj,np)=NUIDO(IDO(nk,nn,1,nb),1,1)
                      ELSE IF(NIT(nb).EQ.2) THEN
                        NUNK(nkp,nj,np)=NUIDO(IDO(nk,nn,1,nb),
     '                    IDO(nk,nn,2,nb),1)
                      ELSE
                        NUNK(nkp,nj,np)=NUIDO(IDO(nk,nn,1,nb),
     '                    IDO(nk,nn,2,nb),IDO(nk,nn,3,nb))
                      ENDIF
                    ENDIF !NUNK(nkp,nj,np).EQ.0
                  ENDDO !nk
                ENDIF
                noelem=noelem+1
              ENDDO !noelem
            ENDDO !njj
          ENDDO !njtype
        ENDDO !no_node
      ENDDO !nr

      CALL EXITS('CALC_NUNK')
      RETURN
 9999 CALL ERRORS('CALC_NUNK',ERROR)
      CALL EXITS('CALC_NUNK')
      RETURN 1
      END


