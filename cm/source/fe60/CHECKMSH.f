      SUBROUTINE CHECKMSH(NBJ,NEELEM,NFVC,NPLIST,NPNE,NPNODE,nr,NVCNODE,
     '  VC,ERROR,*)

C#### Subroutine: CHECKMSH
C###  Description:
C###    Checks the mesh integrity and quality

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NFVC(2,0:NFVCM,NVCM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVCNODE(2,NP_R_M)
      REAL*8 VC(0:NVCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nonode,np,noelem,ne,nb,N_BOUNDARY,N_INTBOUN,N_INTERNAL,
     '  nn,nfvl,nvc,cnonode,cnvc
      REAL*8 VOL_DIFF_AVG,VOL_DIFF_TOT,VOL_DIFF,VOL_DIFF_MAX

      CALL ENTERS('CHECKMSH',*9999)

C     ..See if any elements are made up of boundary and internal nodes
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NPLIST(np)=nonode
      ENDDO
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        nb=NBJ(1,ne)
        N_BOUNDARY=0
        N_INTBOUN=0
        N_INTERNAL=0
        DO nn=1,NNT(nb)
          np=NPNE(nn,nb,ne)
          nonode=NPLIST(np)
          IF(NVCNODE(TYPE,nonode).EQ.BOUNDARY) THEN
            N_BOUNDARY=N_BOUNDARY+1
          ELSEIF(NVCNODE(TYPE,nonode).EQ.INTBOUN) THEN
            N_INTBOUN=N_INTBOUN+1
          ELSEIF(NVCNODE(TYPE,nonode).EQ.INTERNAL) THEN
            N_INTERNAL=N_INTERNAL+1
          ENDIF
        ENDDO
        IF(N_INTERNAL.GT.0.AND.N_BOUNDARY.GT.0) THEN
          WRITE(OP_STRING,'('' Element boundary violation: '//
     '      'ne = '',I6,'' np = '',4I6)') ne,
     '      (NPNE(nn,nb,ne),nn=1,NNT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO

C     ..Obtain volume difference average and maximum
      VOL_DIFF_AVG=0.d0
      VOL_DIFF_MAX=0.d0
      DO nvc=1,NVCT
        VOL_DIFF_TOT=0.d0
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
            cnvc=NVCNODE(MAP,cnonode)
            VOL_DIFF=DABS(VC(nvc)-VC(cnvc))
            VOL_DIFF_TOT=VOL_DIFF_TOT+VOL_DIFF
            VOL_DIFF=VOL_DIFF/VC(nvc)
            IF(VOL_DIFF.GT.VOL_DIFF_MAX) VOL_DIFF_MAX=VOL_DIFF
          ENDIF
        ENDDO
        VOL_DIFF_TOT=VOL_DIFF_TOT/VC(nvc)
        VOL_DIFF_AVG=VOL_DIFF_AVG+VOL_DIFF_TOT/DBLE(NFVC(1,0,nvc))
      ENDDO
      VOL_DIFF_AVG=VOL_DIFF_AVG/(2.d0*DBLE(NVCT))
      WRITE(OP_STRING,'('' Average volume difference = '',D14.4)')
     '  VOL_DIFF_AVG
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Maximum volume difference = '',D14.4)')
     '  VOL_DIFF_MAX
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      CALL EXITS('CHECKMSH')
      RETURN
 9999 CALL ERRORS('CHECKMSH',ERROR)
      CALL EXITS('CHECKMSH')
      RETURN 1
      END


