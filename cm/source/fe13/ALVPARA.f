      SUBROUTINE ALVPARA(nb_host,NELIST,NLL,NPL,NPLIST,NPNE,NXI,centre,
     '  centre2,XP,ERROR,*)

C#### Subroutine: ALVPARA
C###  Description:
C###    ALVPARA calculates alveolar parameters for input into
C###    IPMESH2 to generate the pulmonary alveolar-capillary network.
      
C*** Created by Kelly Burrowes, Jan, 2003.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter variables
      INTEGER nb_host,NELIST(0:NEM),NLL(12,NEM),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 centre(NJT),centre2(NJT),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nae,ne,ne2,ne3,nj,nl,nn,nn2,noelem,noelem2,noelem3,nonode,
     '  np,np1,np2
      REAL*8 TOL
      CHARACTER STRING*(255)
      LOGICAL COLLAPSED,edge,NEW1,NEW2,NODE1,NODE2,
     '  TRIANGULAR

      CALL ENTERS('ALVPARA',*9999)

      IF(NPNE(1,nb_host,NELIST(1)).EQ.NPNE(2,nb_host,NELIST(1))) THEN
        TRIANGULAR=.TRUE.
      ELSE
        TRIANGULAR=.FALSE.
      ENDIF
      DO nj=1,NJT
        centre2(nj)=0.d0 !initialise
      ENDDO 
      NPLIST(0)=0
      DO noelem=1,NELIST(0)
        ne=NELIST(noelem)
        IF(TRIANGULAR) THEN !triangular element
          DO nn=2,NNT(nb_host)
            TOL=1.d-8
            np1=NPNE(nn,nb_host,ne)
            IF(nn.LT.NNT(nb_host)) THEN
              np2=NPNE((nn+1),nb_host,ne) !nodes along each element edge
            ELSE
              np2=NPNE(2,nb_host,ne)
            ENDIF
            COLLAPSED=.FALSE.
            IF(np1.EQ.np2) THEN
              COLLAPSED=.TRUE. !collapsed element
            ELSE IF(XP(1,1,1,np1).LT.XP(1,1,1,np1)+TOL.AND.XP(1,1,1,np1)
     '          .GT.XP(1,1,1,np2)-TOL.AND.XP(1,1,2,np1).LT.XP(1,1,2,np2) 
     '          +TOL.AND.XP(1,1,2,np1).GT.XP(1,1,2,np2)-TOL.AND.
     '          XP(1,1,3,np1).LT.XP(1,1,3,np2)+TOL.AND.XP(1,1,3,np1).GT.
     '          XP(1,1,3,np2)-TOL) THEN
              COLLAPSED=.TRUE.
            ENDIF
            EDGE=.TRUE.
            noelem2=0
            DO WHILE(EDGE.AND.noelem2.LT.NELIST(0))
              noelem2=noelem2+1
              ne2=NELIST(noelem2)
              NODE1=.FALSE.
              NODE2=.FALSE.
              IF(ne.NE.ne2) THEN
                DO nn2=2,NNT(nb_host) 
C... may be mulitple nodes at same position, therefore need to check XP
                  np=NPNE(nn2,nb_host,ne2)
                  IF(np.EQ.np1) THEN
                    NODE1=.TRUE.
                  ELSE IF(XP(1,1,1,np).LT.XP(1,1,1,np1)+TOL.AND.
     '                XP(1,1,1,np).GT.XP(1,1,1,np1)-TOL.AND.XP(1,1,2,np)
     '                .LT.XP(1,1,2,np1)+TOL.AND.XP(1,1,2,np).GT.
     '                XP(1,1,2,np1)-TOL.AND.XP(1,1,3,np).LT.
     '                XP(1,1,3,np1)+TOL.AND.XP(1,1,3,np).GT.
     '                XP(1,1,3,np1)-TOL) THEN
                    NODE1=.TRUE.
                  ELSE IF(np.EQ.np2) THEN
                    NODE2=.TRUE.
                  ELSE IF(XP(1,1,1,np).LT.XP(1,1,1,np2)+TOL.AND.
     '                XP(1,1,1,np).GT.XP(1,1,1,np2)-TOL.AND.XP(1,1,2,np)
     '                .LT.XP(1,1,2,np2)+TOL.AND.XP(1,1,2,np).GT.
     '                XP(1,1,2,np2)-TOL.AND.XP(1,1,3,np).LT.
     '                XP(1,1,3,np2)+TOL.AND.XP(1,1,3,np).GT.
     '                XP(1,1,3,np2)-TOL) THEN
                    NODE2=.TRUE.
                  ENDIF
                ENDDO !nn2
                IF(NODE1.AND.NODE2) EDGE=.FALSE. !not inlet edge
                IF(COLLAPSED) THEN !collapsed ne only 1 np on edge
                  IF(NODE1.OR.NODE2) EDGE=.FALSE.
                ENDIF
              ENDIF !ne.NE.ne2
            ENDDO !WHILE
            IF(EDGE) THEN
              TOL=1.d-2
              NEW1=.TRUE.
              NEW2=.TRUE.
              DO nonode=1,NPLIST(0) !check if np already in list
                np=NPLIST(nonode)
                IF(np.EQ.np1) THEN
                  NEW1=.FALSE.
                ELSE IF(XP(1,1,1,np).LT.XP(1,1,1,np1)+TOL.AND.
     '              XP(1,1,1,np).GT.XP(1,1,1,np1)-TOL.AND.
     '              XP(1,1,2,np).LT.XP(1,1,2,np1)+TOL.AND.
     '              XP(1,1,2,np).GT.XP(1,1,2,np1)-TOL.AND.
     '              XP(1,1,3,np).LT.XP(1,1,3,np1)+TOL.AND.
     '              XP(1,1,3,np).GT.XP(1,1,3,np1)-TOL) THEN
                  NEW1=.FALSE.
                ENDIF
                IF(np.EQ.np2) THEN
                  NEW2=.FALSE.
                ELSE IF(XP(1,1,1,np).LT.XP(1,1,1,np2)+TOL.AND.
     '              XP(1,1,1,np).GT.XP(1,1,1,np2) -TOL.AND.
     '              XP(1,1,2,np).LT.XP(1,1,2,np2)+TOL.AND.
     '              XP(1,1,2,np).GT.XP(1,1,2,np2)-TOL.AND.
     '              XP(1,1,3,np).LT.XP(1,1,3,np2)+TOL.AND.
     '              XP(1,1,3,np).GT.XP(1,1,3,np2)-TOL) THEN
                  NEW2=.FALSE.
                ENDIF
              ENDDO
              IF(NEW1) THEN
                NPLIST(0)=NPLIST(0)+1
                NPLIST(NPLIST(0))=np1
              ENDIF
              IF(NEW2) THEN
                NPLIST(0)=NPLIST(0)+1
                NPLIST(NPLIST(0))=np2
              ENDIF
            ENDIF !EDGE
          ENDDO !nn
        ELSE !NOT.TRIANGULAR
          DO nj=-NIT(nb_host),NIT(nb_host)
            EDGE=.FALSE.
            IF(NXI(nj,0,ne).EQ.0) THEN !element is on edge
              EDGE=.TRUE.
            ELSE
              EDGE=.TRUE.                              
              DO noelem2=1,NXI(nj,0,ne)
                ne2=NXI(nj,noelem2,ne) !adjacent element #
 !if ne2 is not in current group then ne must be on edge
                DO noelem3=1,NELIST(0)
                  ne3=NELIST(noelem3)
                  IF(ne2.EQ.ne3) EDGE=.FALSE.
                ENDDO
              ENDDO
            ENDIF
            IF(EDGE) THEN
              IF(nj.EQ.-2) nae=1 !Xi_-2=1 !not sure if an array has 
              IF(nj.EQ.-1) nae=3 !Xi_-1=3 !this information ??
              IF(nj.EQ.1) nae=4 !Xi_+1=4
              IF(nj.EQ.2) nae=2 !Xi_+2=2
              nl=NLL(nae,ne) !line segment #
              IF(nl.NE.0) THEN
                np1=NPL(2,1,nl)
                np2=NPL(3,1,nl)
                NEW1=.TRUE.
                NEW2=.TRUE.
                IF(np1.eq.np2) THEN
                  NEW1=.FALSE.
                  NEW2=.FALSE.
                ENDIF
                DO nonode=1,NPLIST(0) !check if np already in list
                  np=NPLIST(nonode)
                  IF(np.EQ.np1) THEN
                    NEW1=.FALSE.
                  ENDIF
                  IF(np.EQ.np2) THEN
                    NEW2=.FALSE.
                  ENDIF
                ENDDO
                IF(NEW1) THEN
                  NPLIST(0)=NPLIST(0)+1
                  NPLIST(NPLIST(0))=np1
                ENDIF
                IF(NEW2) THEN
                  NPLIST(0)=NPLIST(0)+1
                  NPLIST(NPLIST(0))=np2
                ENDIF
              ENDIF
            ENDIF !NXI
          ENDDO !nj
        ENDIF
      ENDDO !noelem
      STRING='edge'
      CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)
C... calculates rotation anlge for points of sphere      
      CALL ALVANGLE(NPLIST,centre,centre2,XP,TRIANGULAR,ERROR,*9999)
      
      CALL EXITS('ALVPARA')
      RETURN
 9999 CALL ERRORS('ALVPARA',ERROR)
      CALL EXITS('ALVPARA')
      RETURN 1
      END


