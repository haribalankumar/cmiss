      SUBROUTINE DEXI_NONLIN_CONTAIN(IBT,IDO,INP,LD,LDTEMP,NBJ,NBH,
     '  nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,NKHE,NKJE,NPF,
     '  NPNE,nolist,nr,NRE,NVHE,NVJE,NW,nx,NXI,CURVCORRECT,
     '  SE,XA,XE,XI,XID,XIQ,XP,XQ,
     '  ZA,ZD,ZP,DEFORM,FOUND,GRPGRID,PASS2,
     '  SPECIFY,ERROR,*)

C#### Subroutine: DEXI_NONLIN_CONTAIN
C###  Description:
C###    DEXI_NONLIN_CONTAIN

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LDTEMP,NBJ(NJM,NEM),NBH(NHM,NCM,NEM),nd,ND0,ND1,ne,
     '  NELIST(0:NEM),NHE(NEM,NXM),ni,NITB,nj1,
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),nolist,nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XID(NIM,NDM),
     '  XIQ(NIM,NQM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORM,FOUND,GRPGRID,SPECIFY
      CHARACTER ERROR*(*)
      LOGICAL PASS2
!     Local Variables
      INTEGER IT,nb,ne0,ne_close,nj,nk,noelem,nonode,np,nv
      REAL*8 centre(3),closest,distance,LOOSE_TOL,Z_LEN
      LOGICAL INELEM

      CALL ENTERS('DEXI_NONLIN_CONTAIN',*9999)

      LOOSE_TOL = 1d-8

      DO nd=ND0,ND1
! find the closest element to initiate the search
        closest=1.d8

        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(NJ1,ne)
          DO nj=1,NJT
            centre(nj)=0.d0
          ENDDO
!! find the right array to use for number of nodes!!
          DO nonode=1,8
             np=NPNE(nonode,nb,ne)
             nk=1 !coordinates
             nv=1 !base estimate on first version only
             DO nj=1,NJT
                centre(nj)=centre(nj)+XP(nk,nv,nj,np)
             ENDDO
          ENDDO
          distance=0.d0
          DO nj=1,NJT
            distance=distance+(ZD(nj,nd)-centre(nj)/8.d0)**2
          ENDDO
          distance=DSQRT(distance)
          IF(distance.LT.closest)THEN
            ne_close=ne
            closest=distance
          ENDIF
        ENDDO !noelem

        ne=ne_close

        IF(DEFORM)THEN
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     &      NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &      NVHE(1,1,1,ne),NW(ne,1),nx,
     &      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     &      XE,ZP,ERROR,*9999)
        ELSE
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     &      NPF(1,1),NPNE(1,1,ne),
     &      NRE(ne),NVJE(1,1,1,ne),
     &      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        ENDIF
          
        FOUND=.FALSE.
        NITB=NIT(NBJ(NJ1,ne))
        DO ni=1,NITB
           XI(ni)=0.5d0
        ENDDO
        
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,
     &    XID(1,nd),ZD(1,nd),FOUND,INELEM,ERROR,*9999)

        IF(FOUND)THEN
          LD(nd)=ne
        ELSE
! check each of the adjacent elements
          ne0=ne
          closest=1.d8
          DO ni=1,NITB
            DO noelem=1,NXI(-ni,0,ne0)
              ne=NXI(-ni,noelem,ne0)
              DO nonode=1,8
                np=NPNE(nonode,nb,ne)
                nk=1          !coordinates
                nv=1          !base estimate on first version only
                DO nj=1,NJT
                  centre(nj)=centre(nj)+XP(nk,nv,nj,np)
                ENDDO
              ENDDO
              distance=0.d0
              DO nj=1,NJT
                distance=distance+(ZD(nj,nd)-centre(nj)/8.d0)**2
              ENDDO
              distance=DSQRT(distance)
              IF(distance.LT.closest)THEN
                ne_close=ne
                closest=distance
              ENDIF
            ENDDO

            DO noelem=1,NXI(ni,0,ne0)
               ne=NXI(ni,noelem,ne0)
               DO nonode=1,8
                 np=NPNE(nonode,nb,ne)
                 nk=1          !coordinates
                 nv=1          !base estimate on first version only
                 DO nj=1,NJT
                   centre(nj)=centre(nj)+XP(nk,nv,nj,np)
                 ENDDO
               ENDDO
               distance=0.d0
               DO nj=1,NJT
                 distance=distance+(ZD(nj,nd)-centre(nj)/8.d0)**2
               ENDDO
               distance=DSQRT(distance)
               IF(distance.LT.closest)THEN
                 ne_close=ne
                 closest=distance
               ENDIF
            ENDDO 
          ENDDO !ni
           
          ne=ne_close
          IF(DEFORM)THEN
             CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     &            NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &            NVHE(1,1,1,ne),NW(ne,1),nx,
     &            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     &            XE,ZP,ERROR,*9999)
          ELSE
             CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     &            NPF(1,1),NPNE(1,1,ne),
     &            NRE(ne),NVJE(1,1,1,ne),
     &            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          ENDIF
          
          FOUND=.FALSE.
          NITB=NIT(NBJ(NJ1,ne))
          DO ni=1,NITB
             XI(ni)=0.5d0
          ENDDO
          
          INELEM=.FALSE.
          IT=0
          CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,
     &         XID(1,nd),ZD(1,nd),FOUND,INELEM,ERROR,*9999)
          
          IF(FOUND)THEN
             LD(nd)=ne
          ELSE
             write(*,*) 'not found for',nd
          ENDIF
       ENDIF
      ENDDO ! nd=1,NDT loop

      CALL EXITS('DEXI_NONLIN_CONTAIN')
      RETURN
 9999 CALL ERRORS('DEXI_NONLIN_CONTAIN',ERROR)
      CALL EXITS('DEXI_NONLIN_CONTAIN')
      RETURN 1
      END


