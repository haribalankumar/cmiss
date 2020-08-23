      SUBROUTINE CHNODE(ISEG,ISELNO,ISLINE,ISLINO,ISNONO,
     '  IIW,MXI,NBJ,NEELEM,NHP,NKH,NLL,NLLIST,np,NPL,NPNE,
     '  NPNODE,nx,NYNP,DL,XP,ZC,ZP,CSEG,BOTH,DEFORM,FIX,FRAME,ERROR,*)

C#### Subroutine: CHNODE
C###  Description:
C###    CHNODE redefines position of individual node.
C###    np is node number of deleted and then recreated node

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER IIW,ISEG(*),ISELNO(NWM,NEM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),
     '  ISNONO(NWM,NPM),MXI(2,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NLL(12,NEM),NLLIST(0:NLM),np,
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM),ZC(NJM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL BOTH,DEFORM,FIX(NYM,NIYFIXM,NXM),FRAME
!     Local Variables
      INTEGER INDEX,INDEX_POLYLINE,INSTAT,iw,nb,ne,nj,nn,nolist,
     '  noelem,nr
      REAL*8 XX(3),ZZ(3)

      CALL ENTERS('CHNODE',*9999)

C *** Redefine node position
      CALL DENODE(1,INSTAT,ISEG,ISNONO,IIW,NHP(1,0,nx),NKH,np,
     '  NPNODE,NYNP,XP,BOTH,DEFORM,FIX,FRAME,CSEG,'OLD',ERROR,*9999)
      NLLIST(0)=NLT
      DO nolist=1,NLLIST(0)
        NLLIST(nolist)=nolist
      ENDDO
      iw=IIW
      CALL ACWK(iw,1,ERROR,*9999)

C *** Update line segments
C??? CS 17/2/97 crashed if no lines were defined
C      IF(ISLINE(iw,NTLINE).GT.0) THEN
      IF(ISLINE(iw,1).GT.0) THEN
        INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
        DO nr=1,NRT
          CALL SGLINE(INDEX,ISEG,ISLINE(iw,NTLINE),ISLINO(iw),iw,NLLIST,
     '      NTLINE,NPL,nr,nx,CSEG,'UNDEFORMED',DL,'SOLID',.TRUE.,XP,ZP,
     '      ERROR,*9999)
        ENDDO !nr
      ENDIF

C *** Update element segments
      IF(ISELNO(iw,1).GT.0) THEN
        DO nr=1,NRT
          DO 650 noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            DO nn=1,NNT(nb)
              IF(NPNE(nn,nb,ne).eq.np) GO TO 606
            ENDDO
            GO TO 650
 606        DO nj=1,NJT
              ZC(nj,ne)=0.d0
            ENDDO
            DO nn=1,NNT(nb)
              DO nj=1,NJT
                nb=NBJ(nj,ne)
                XX(nj)=XP(1,1,nj,NPNE(nn,nb,ne))
              ENDDO
              CALL XZ(ITYP10(nr),XX,ZZ)
              DO nj=1,NJT
                ZC(nj,ne)=ZC(nj,ne)+ZZ(nj)
              ENDDO
            ENDDO
            DO nj=1,NJT
              nb=NBJ(nj,ne)
              ZC(nj,ne)=ZC(nj,ne)/DBLE(NNT(nb))
            ENDDO
            IF(iw.ne.5.and.iw.ne.6) THEN
              CALL SGELEM(0,ISEG,ISELNO(iw,ne),iw,MXI(1,ne),NBJ(1,ne),
     '          ne,NLL(1,ne),NPL,nr,.FALSE.,.FALSE.,CSEG,DL,XP,
     '          ZC(1,ne),ERROR,*9999)
            ENDIF
 650      CONTINUE
        ENDDO !nr
      ENDIF !iselno

      CALL DAWK(iw,1,ERROR,*9999)

      CALL EXITS('CHNODE')
      RETURN
 9999 CALL ERRORS('CHNODE',ERROR)
      CALL EXITS('CHNODE')
      RETURN 1
      END


