      SUBROUTINE CAP_IO(nb,NBJ,NEELEM,NKJ,NKJE,NPNE,NPNODE,nr,nr_feed,
     &  NRE,NVJE,NVJP,SE,XP,ERROR,*)

C####  Subroutine: CAP_IO
C###   Description:
C###     CAP_IO creates an inlet arterioles and outlet venules for
C###     the pulmonary capillary network.

C***   Last modified December, 2002.
C***   Created by Kelly Burrowes, 8th March, 2002.

C***   This subroutine still under development, will be tidied up
C***   when larger pulmonary vessels modelled.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter values
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NKJ(NJM,NPM),
     &  NKJE(NKM,NNM,NJM,NEM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,nr_feed,NRE(NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ne,nj,noelem,nonode,nonode2,nonode3,np,np2,np3,np_closest
      REAL*8 ADD,length,min_length
      LOGICAL FIRST

      CALL ENTERS('CAP_IO',*9999)

      nonode=0
      np2=NPT(0) !highest global node number
      nonode2=NPNODE(0,nr) !highest node # in region nr
      ne=NET(0)
      noelem=NEELEM(0,nr) !highest ne # in region nr
      ADD=0.1d0

C... read in inlet & outlet nodes & elements from feed region
C... put these nodes & elements into mesh region (nr)
        FIRST=.TRUE.
        IF(nr_feed.NE.0) THEN
         IF(NPNODE(0,nr_feed).NE.0) THEN !feed nodes read in
          DO nonode=1,NPNODE(0,nr_feed)
            np=NPNODE(nonode,nr_feed)
            nonode2=nonode2+1
            NPNODE(nonode2,nr)=np            
C... Find closest capillary node
            min_length=1.d6
            DO nonode3=1,NPNODE(0,nr)
              np3=NPNODE(nonode3,nr)
              length=0.d0
              DO nj=1,NJT
                length=length+(XP(1,1,nj,np3)-XP(1,1,nj,np))**2.d0
              ENDDO
              length=DSQRT(length)
              IF(length.LT.min_length) THEN
                np_closest=np3
                min_length=length
              ENDIF 
            ENDDO !nonode
C... Now can make feed elements
            noelem=noelem+1
            ne=ne+1
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=np
            NPNE(2,nb,ne)=np_closest
            IF(FIRST) THEN
              CAP_INLET=ne
              FIRST=.FALSE.
            ENDIF
            CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,np,np_closest,nr,NRE,NVJE,
     &        NVJP,SE,ERROR,*9999)
          ENDDO !nonode
        ELSE
          CALL ASSERT(NPNODE(0,nr_feed).NE.0,
     &      '>>No nodes defined in feed region',ERROR,*9999)
        ENDIF
      ELSE
C... if no feed nodes input, just allocate 1st and last nodes in region
C... inlet
          np=NPNE(1,nb,NEELEM(1,nr)) !1st node of 1st elem in region
          nonode2=nonode2+1
          np2=np2+1
          NPNODE(nonode2,nr)=np2 !new node
          DO nj=1,NJT
            XP(1,1,nj,np2)=XP(1,1,nj,np)-ADD
          ENDDO !nj
          noelem=noelem+1
          ne=ne+1
          CAP_INLET=ne !highest global ne # + 1
          WRITE(OP_STRING,'('' Inlet vessel ne #= '',I5)') ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          NEELEM(noelem,nr)=CAP_INLET
          NPNE(1,nb,CAP_INLET)=np2 !Creates arteriole inlet element
          NPNE(2,nb,CAP_INLET)=np
          CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,np2,np,nr,NRE,NVJE,NVJP,SE,
     &      ERROR,*9999)
C... outlet
          np=NPNE(2,nb,NEELEM(NEELEM(0,nr),nr)) !2nd np of last ne in nr
          nonode2=nonode2+1
          np2=np2+1
          NPNODE(nonode2,nr)=np2 !new node
          DO nj=1,NJT
            XP(1,1,nj,np2)=XP(1,1,nj,np)+ADD
          ENDDO !nj
          noelem=noelem+1
          ne=ne+1
          CAP_OUTLET=ne
          WRITE(OP_STRING,'('' Outlet vessel ne #= '',I5)') ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          NEELEM(noelem,nr)=CAP_OUTLET
          NPNE(1,nb,CAP_OUTLET)=np !Creates arteriole inlet element
          NPNE(2,nb,CAP_OUTLET)=np2
          CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,np,np2,nr,NRE,NVJE,NVJP,SE,
     &      ERROR,*9999)
        ENDIF
C      ENDIF !VORONOI_ALV

C     CALL ASSERT(INLET.EQ.1,'>>No inlet vessel defined',ERROR,*9999)
C     CALL ASSERT(OUTLET.EQ.1,'>>No outlet vessel defined',ERROR,*9999)
      NEELEM(0,nr)=noelem
      NEELEM(0,0)=ne
      NET(0)=ne
      NET(nr)=noelem
      NPNODE(0,nr)=nonode2
      NPNODE(0,0)=np2
      NPT(0)=np2
      NPT(nr)=nonode2

      CALL EXITS('CAP_IO')
      RETURN
 9999 CALL ERRORS('CAP_IO',ERROR)
      CALL EXITS('CAP_IO')
      RETURN 1
      END


