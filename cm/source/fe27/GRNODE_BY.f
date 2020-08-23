      SUBROUTINE GRNODE_BY(NBJ,NEELEM,NPLIST,NPNE,nr,NXI,CE,
     &  radius,XP,xyz,STRING2,CENTRE,ERROR,*)

C#### Subroutine: GRNODE_BY
C###  Description:
C###    GRNODE_BY groups nodes by generation, order, lobe, or terminal.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),nr,NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),radius,XP(NKM,NVM,NJM,NPM),xyz(3)
      CHARACTER ERROR*(*),STRING2*255
      LOGICAL CENTRE
!     Local Variables
      INTEGER nb,ne,ne0,nj,noelem,np
      REAL*8 distance,SUM

      CALL ENTERS('GRNODE_BY',*9999)

      NPLIST(0)=0

      IF(CENTRE)THEN
        DO noelem=1,NEELEM(0,nr)
           ne=NEELEM(noelem,nr)
           IF(NXI(1,0,ne).EQ.0)THEN !is a terminal element
              nb=NBJ(1,ne) !basis function for geometry
              np=NPNE(2,nb,ne) !end (terminal) node number
              distance=0.d0
              DO nj=1,NJT
                 distance=distance+(XP(1,1,nj,np)-xyz(nj))**2
              ENDDO !nj
              distance=DSQRT(distance)
              IF(distance.LE.radius)THEN !is within the specified radius
                 NPLIST(0)=NPLIST(0)+1
                 NPLIST(NPLIST(0))=np
              ENDIF !if within the radius
           ENDIF !if terminal
        ENDDO !noelem
      ELSE
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne)
          IF(NXI(1,0,ne).NE.0)THEN
            ne0=NXI(1,1,ne)
            IF(CE(5,ne).EQ.1.d0.AND.CE(5,ne0).LT.1.d0)THEN
              NPLIST(0)=NPLIST(0)+1
              NPLIST(NPLIST(0))=np
              IF(DOP)THEN
                SUM=CE(1,ne)*PI*CE(4,ne)**2.d0
                ne0=NXI(-1,1,ne)
                DO WHILE(ne0.NE.0)
                  SUM=SUM+CE(1,ne0)*PI*CE(4,ne0)**2.d0
                  ne0=NXI(-1,1,ne0)
                ENDDO
              ENDIF !DOP
            ENDIF
          ELSE
            IF(CE(5,ne).EQ.1.d0)THEN
              NPLIST(0)=NPLIST(0)+1
              NPLIST(NPLIST(0))=np
              SUM=CE(1,ne)*PI*CE(4,ne)**2.d0
              IF(DOP)THEN
                ne0=NXI(-1,1,ne)
                DO WHILE(ne0.NE.0)
                  SUM=SUM+CE(1,ne0)*PI*CE(4,ne0)**2.d0
                  ne0=NXI(-1,1,ne0)
                ENDDO
              ENDIF !DOP
            ENDIF
          ENDIF !NXI
        ENDDO !noelem
      ENDIF

      CALL EXITS('GRNODE_BY')
      RETURN
 9999 CALL ERRORS('GRNODE_BY',ERROR)
      CALL EXITS('GRNODE_BY')
      RETURN 1
      END


