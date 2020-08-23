      SUBROUTINE DENODE_FROM_DATA(LD,NBJ,NDDATA,NDP,NENP,NPLIST,NPNE,
     &  nr,WD,XID,XP,ZD,ERROR,*)
      
C#### Subroutine: DENODE_FROM_DATA
C###  Description:
C###  NODE_DATA calculates data positions to coincide with node list


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),NDP(NDM),
     &  NENP(NPM,0:NEPM,0:NRM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),nr
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,nj,nonode,np

      CALL ENTERS('DENODE_FROM_DATA',*9999)

      DO nonode=1,NPLIST(0) !for each node in list
        np=NPLIST(nonode)
        NDT=NDT+1
        DO nj=1,3
          ZD(nj,NDT)=XP(1,1,nj,np) !single version, coordinates only
          WD(nj,NDT)=1.0d0
          XID(nj,NDT)=0.d0
        ENDDO !nj
        NDP(NDT)=np !stores the node associated with data
        NDDATA(NDT,nr)=NDT !store data point number
        NDDATA(0,nr)=NDDATA(0,nr)+1 !increment # of data points
        ne=NENP(np,1,nr) !first element node is in
        LD(NDT)=ne !store a single element that the node is in
        nb=NBJ(1,ne) !basis function for geometry of element ne

        IF(NIT(nb).EQ.1)THEN !1 Xi direction
          IF(NPNE(2,nb,ne).EQ.np) XID(1,NDT)=1.d0

        ELSE IF(NIT(nb).EQ.2)THEN !2 Xi directions
          IF(NPNE(1,nb,ne).EQ.np)THEN
          ELSE IF(NPNE(2,nb,ne).EQ.np)THEN
            XID(1,NDT)=1.d0
          ELSE IF(NPNE(3,nb,ne).EQ.np)THEN
            XID(2,NDT)=1.d0
          ELSE IF(NPNE(4,nb,ne).EQ.np)THEN
            XID(1,NDT)=1.d0
            XID(2,NDT)=1.d0
          ENDIF

        ELSE IF(NIT(nb).EQ.3)THEN !3 Xi directions
          IF(NPNE(1,nb,ne).EQ.np)THEN
          ELSE IF(NPNE(2,nb,ne).EQ.np)THEN
            XID(1,NDT)=1.d0
          ELSE IF(NPNE(3,nb,ne).EQ.np)THEN
            XID(2,NDT)=1.d0
          ELSE IF(NPNE(4,nb,ne).EQ.np)THEN
            XID(1,NDT)=1.d0
            XID(2,NDT)=1.d0
          ELSE IF(NPNE(5,nb,ne).EQ.np)THEN
            XID(3,NDT)=1.d0
          ELSE IF(NPNE(6,nb,ne).EQ.np)THEN
            XID(1,NDT)=1.d0
            XID(3,NDT)=1.d0
          ELSE IF(NPNE(7,nb,ne).EQ.np)THEN
            XID(2,NDT)=1.d0
            XID(3,NDT)=1.d0
          ELSE IF(NPNE(8,nb,ne).EQ.np)THEN
            XID(1,NDT)=1.d0
            XID(2,NDT)=1.d0
            XID(3,NDT)=1.d0
          ENDIF
        ENDIF
      ENDDO !nonode
        
      CALL EXITS('DENODE_FROM_DATA')
      RETURN
 9999 CALL ERRORS('DENODE_FROM_DATA',ERROR)
      CALL EXITS('DENODE_FROM_DATA')
      RETURN 1
      END

