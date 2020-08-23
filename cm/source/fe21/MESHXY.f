      SUBROUTINE MESHXY(IBT,IDO,INP,NBJ,nr,NEELEM,NPNE,DEPTH,SE,
     '             XP,WIDTH,Z,ERROR,*)

C#### Subroutine: MESHXY
C###  Description:
C###    MESHXY calculates the width of mesh in x and y directions.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
C      INCLUDE 'cmiss$reference:opti00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),nr,NEELEM(0:NE_R_M,0:NRM),
     '  NPNE(NNM,NBFM,NEM)
      REAL*8 DEPTH,SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),WIDTH,Z
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER max,nb,ne,nj,nj2,nk,nn,nn1,noelem,ns,i,SDE(3,2),zne
      REAL*8 DFN,FN,PTS(4,2),PSI1,TEMP,XI(3)

      CALL ENTERS('MESHXY',*9999)
      zne=0
C     Loop over elem find sign changes
      IF(NJT.EQ.2) THEN
        DO ne=1,NET(nr)
          DO nj=1,NJT
            nb=NBJ(nj,ne)
            IF(XP(1,1,nj,NPNE(1,nb,ne)).LE.zero_tol.AND.
     '         XP(1,1,nj,NPNE(2,nb,ne)).GE.zero_tol) THEN
              SDE(nj,1)=ne
            ENDIF
            IF(XP(1,1,nj,NPNE(2,nb,ne)).LE.zero_tol.AND.
     '            XP(1,1,nj,NPNE(1,nb,ne)).GE.zero_tol) THEN
              SDE(nj,2)=ne
            ENDIF
          ENDDO !nj
        ENDDO !ne
        ELSEIF(NJT.EQ.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(3,ne)
            IF(NBJ(3,ne).NE.1) NPNE(4,nb,ne)=NPNE(3,nb,ne)
            IF((XP(1,1,3,NPNE(3,nb,ne)).LE.Z.AND.
     '       XP(1,1,3,NPNE(1,nb,ne)).GE.Z).OR.(XP(1,1,3,NPNE(2,nb,ne))
     '       .LE.Z.AND.XP(1,1,3,NPNE(4,nb,ne)).GE.Z)) THEN
              DO nj=1,2
                IF(XP(1,1,nj,NPNE(1,nb,ne)).LE.zero_tol.AND.
     '           XP(1,1,nj,NPNE(2,nb,ne)).GT.-zero_tol) THEN
                SDE(nj,1)=ne
                zne=ne
                ENDIF
                IF(XP(1,1,nj,NPNE(2,nb,ne)).LE.zero_tol.AND.
     '           XP(1,1,nj,NPNE(1,nb,ne)).GT.-zero_tol) THEN
                 SDE(nj,2)=ne
                ENDIF
              ENDDO !nj
            ELSEIF((XP(1,1,3,NPNE(3,nb,ne)).GE.Z.AND.
     '       XP(1,1,3,NPNE(1,nb,ne)).LE.Z).OR.(XP(1,1,3,NPNE(2,nb,ne))
     '       .LE.Z.AND.XP(1,1,3,NPNE(4,nb,ne)).GE.Z)) THEN
              DO nj=1,2
                IF(XP(1,1,nj,NPNE(1,nb,ne)).LE.zero_tol.AND.
     '           XP(1,1,nj,NPNE(2,nb,ne)).GT.-zero_tol) THEN
                 SDE(nj,1)=ne
                 zne=ne
                ENDIF
                IF(XP(1,1,nj,NPNE(2,nb,ne)).LE.zero_tol.AND.
     '           XP(1,1,nj,NPNE(1,nb,ne)).GT.-zero_tol) THEN
                 SDE(nj,2)=ne
                ENDIF
              ENDDO !nj
            ENDIF
          ENDDO !noelem
          CALL ASSERT(zne.GT.0,'>>Invalid measurement request',
     '    ERROR,*9999)
        ENDIF
C     Carry out Newton-Raphson within elem to find axes cross pt.
        DO nn=1,4
          DO nj=1,2
            PTS(nn,nj)=0.0d0
          ENDDO
        ENDDO
        DO nj=1,2
          DO nn1=1,2
            ne=SDE(nj,nn1)
            nb=NBJ(3,ne)
            XI(1)=0.5d0
            XI(3)=0.d0
            TEMP=0.d0
            i=1
            DO WHILE(i.LT.6)
            IF(NJT.EQ.3) THEN
            CALL CHMESH_FIND_ZETA(IBT,IDO,INP,max,nb,ne,NPNE,SE,
     '        XP,Z,XI(1),XI(2),ERROR,*9999)
              i=i+1
            ELSE
              XI(2)=0.0d0
              i=6
            ENDIF
            DO WHILE(DABS(XI(1)-TEMP).GT.zero_tol)
              TEMP=XI(1)
              FN=0.d0
              DFN=0.d0
              DO nn=1,NNT(nb)
                DO nk=1,NKT(nn,nb)
                  IF(NJT.EQ.2) ns=(nn-1)*2+nk
                  IF(NJT.EQ.3) ns=(nn-1)*4+nk
                  FN=FN+PSI1(IBT,IDO,INP,nb,1,nk,nn,XI)*
     '              XP(nk,1,nj,NPNE(nn,nb,ne))*SE(ns,nb,ne)
                  DFN=DFN+PSI1(IBT,IDO,INP,nb,2,nk,nn,XI)*
     '              XP(nk,1,nj,NPNE(nn,nb,ne))*SE(ns,nb,ne)
                ENDDO !nk
              ENDDO !nn
              IF(dabs(DFN).GT.zero_tol) THEN
                XI(1)=XI(1)-(FN/DFN)
              ELSE
                XI(1)=0.d0
              ENDIF
            ENDDO
            ENDDO
            nj2=(nj+1)-2*(nj-1)
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                IF(NJT.EQ.2) ns=(nn-1)*2+nk
                 IF(NJT.EQ.3) ns=(nn-1)*4+nk
                 PTS(nn1,nj)=PTS(nn1,nj)+
     '             PSI1(IBT,IDO,INP,nb,1,nk,nn,XI)*
     '             XP(nk,1,nj2,NPNE(nn,nb,ne))*SE(ns,nb,ne)
              ENDDO !nk
            ENDDO !nn
          ENDDO !nn1
        ENDDO !nj


C       Calc mesh height and width
        DEPTH=DABS(PTS(1,1)-PTS(2,1))
        WIDTH=DABS(PTS(1,2)-PTS(2,2))



      CALL EXITS('MESHXY')
      RETURN
 9999 CALL ERRORS('MESHXY',ERROR)
      CALL EXITS('MESHXY')
      RETURN 1
      END

