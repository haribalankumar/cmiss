
      SUBROUTINE UPXI_POINTS_VERSIONS
     & (FD,from_field,IBT,IDO,INP,LD,NBH,NBHF,NBJ,
     &  NBJF,NCLOSFACE,NDP,NDP_INV,NFF,NFLIST,NHE,NKEF,NKHE,NKJE,
     &  NKJF,NNF,NPF,NPLIST,NPNE,NPNF,NRE,NVHE,NVJE,NVJF,NVJP,nx,NXI,
     &  to_field,FDSQ,SE,SF,SQ,XA,XE,XI_1,
     &  XI_2,XI_3,XID,XP,ZD,ZE,ZP,EDGE,ERROR,*)

C#### Subroutine: UPXI_POINTS
C###  Description:
C###    UPXI_POINTS calculates the location of a data point projection
C###    to the surface of a volume mesh (stored in fields), where the
C###    data points are each mapped to a node in a volume mesh. The 
C###    projection is stored in fields. 

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      
!     Parameter List
      INTEGER FD(NDM),from_field(3),IBT(3,NIM,NBFM),
     &  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),LD(NDM),
     &  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     &  NPLIST(0:NPM),NDP(NDM),NDP_INV(NDM),NFF(6,NEM),NFLIST(0:NFM),
     &  NHE(NEM,NXM),NKEF(0:4,16,6,NBFM),
     &  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),
     &  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NRE(NEM),
     &  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJF(NNM,NBFM,NJM),NVJP(NJM,NPM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  to_field(3)

      INTEGER VD(NDM), FD1(NDM), FD2(NDM), VD1(NDM), VD2(NDM),FDg(NDM)

      REAL*8 SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     &  SQ(NDM),XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),
     &  XI_1,XI_2,XI_3,
     &  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     &  ZD(NJM,NDM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM),Xmod(3,11,3)
      CHARACTER ERROR*(*)
      LOGICAL EDGE

!     Local Variables
      INTEGER nb,nbf,NCLOSFACE(0:NFM),ncf,ncfm,nd,ne,nf,nef,ni,nj,
     &  NKJF(NKM,NNM,NJM),nn,noface,np_nf,nr,N1LIST,IT,ITMAX,
     &  nflast,nflast2,nodata,np,nv,mxif(-3:3),
     &  xidir(1:6),nu,nured

      INTEGER nvd

      REAL*8 PXI,SQND,SQNDlast,X(3),XIlast(3),FDSQ(0:NFM),TEMP,
     &  SUM,DSDXI(3),DLAMCH,SCALEFACTOR(3)

      REAL*8 PXIDERIV

      LOGICAL FOUND,INLIST

      PARAMETER(ITMAX=20)

      CALL ENTERS('UPXI_POINTS',*9999)

      ncfm=2
      !Map to change direction indices from NXI to NFF format
      mxif(-1)=1
      mxif(1)=2
      mxif(-2)=3
      mxif(2)=4
      mxif(-3)=5
      mxif(3)=6
C.....Set up inverse mapping for nd to NDP
      DO nodata=1,NDT
        np=NDP(nodata) !the node number that corresponds to data point
        NDP_INV(np)=nodata !the data point that corresponds to node np
      ENDDO

      write(*,*) 'NPLIST0 =',NPLIST(0)
      write(*,*) 'N1LIST=',N1LIST
      DO nodata=1,NPLIST(0) !NPLIST contains a list of nodes
        np=NPLIST(nodata) !the node number in the list
        nd=NDP_INV(np) !the data point for node np
        SQ(nd)=0.0d0
        LD(nd)=0
        FD(nd)=0
      ENDDO


!=====================================================
      DO nodata=1,NPLIST(0)
!=====================================================
      
        np=NPLIST(nodata) !the node number in the list
        nd=NDP_INV(np) !the data point for node np
        write(*,*) '=================================='
        write(*,*) 'number, np, nd=',nodata,np,nd
        write(*,*) '=================================='
        DO ncf=0,NFM
          NCLOSFACE(ncf)=0
          FDSQ(ncf)=0.0d0
        ENDDO

        ncfm=0
        DO noface=1,NFLIST(0)
          nf=NFLIST(noface) !global face number
          write(*,*) np, 'checking for face=',nf
          ne=NPF(6,nf)      !first element nf is in
C note that this subroutine projects surface nodes, so the faces in
C NFLIST must only appear in a single element. i.e. nf is unique to ne
          nef=NPF(8,nf)     !local face number of nf in ne (1st, 2nd, ...)
          nr=NRE(ne)        !region ne is in
C upon return from CALC_FACE_INFORMATION_IND, NPNF contains the
C global node numbers for face nf in element ne
         CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),
     &    nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
     &    nr,NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
          nbf=NBJF(1,nf)
          DO nn=1,4
            np_nf=NPNF(nn,nbf)
            IF(np_nf.EQ.np)THEN !include this face in the search
              ncfm=ncfm+1
              NCLOSFACE(ncfm)=nf
              FDSQ(ncfm)=1.d7
            ENDIF
          ENDDO
        ENDDO !noface (face counter)
        SQ(nd)=1.d7

        DO ni=1,3
          XID(ni,nd)=0.5d0
        ENDDO

        write(*,*) 'ncfm = ',ncfm
        DO noface=1,ncfm
          nf=NCLOSFACE(noface) ! starting face
          ne=NPF(6,nf)
          IT=0
          FOUND=.FALSE.
          nflast2=-1
          nflast=-1
          
          XIlast(1)=XI(1)
          XIlast(2)=XI(2)
          SQNDlast=SQND

          DO ni=1,3
            XI(ni)=XID(ni,nd)
          ENDDO

          nef=NPF(8,nf)
          nr=NRE(ne)
          nb=NBJ(1,ne)
          xidir(2)=NPF(1,nf)
          xidir(1)=-xidir(2) 
          xidir(4)=NPF(3,nf)
          xidir(3)=-xidir(4) 
          xidir(6)=NNF(1,nef,nb)
          xidir(5)=-xidir(6)
          

          CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),
     &      nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
     &      nr,NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
          CALL XPXE_POINTS(from_field,NBJF(1,nf),NKJF,NPF(1,nf),NPNF,
     &      nr,NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)

          FOUND=.TRUE.

          CALL PROJ_ORTHOG(IBT,IDO,INP,NBJF(1,nf),SQND,
     &      XE,XI,ZD(1,nd),FOUND,ERROR,*9999)

          IF(FOUND.AND.SQND.LT.SQ(nd))THEN

            SQ(nd)=SQND
            FD(nd)=nef
            FDg(nd)=nf
            LD(nd)=ne
             
            write(*,*) 'nef=',nef

             IF(MOD(nef,2).EQ.1) THEN
               XI(3)=0.0d0
             ELSE
               XI(3)=1.0d0
             ENDIF
             IF(EDGE)THEN
                IF(XI_1.NE.2.d0)THEN
                   XID(1,nd)=XI_1
                ELSE
                   XID(1,nd)=XI(1)
                ENDIF
                IF(XI_2.NE.2.d0)THEN
                   XID(2,nd)=XI_2
                ELSE
                   XID(2,nd)=XI(2)
                ENDIF
                IF(XI_3.NE.2.d0)THEN
                   XID(3,nd)=XI_3
                ELSE
                   XID(3,nd)=XI(3)
                ENDIF
             ELSE
                XID(NPF(1,nf),nd)=XI(1)
                XID(NPF(3,nf),nd)=XI(2)
                XID(NNF(1,nef,nb),nd)=XI(3)
             ENDIF
          ENDIF ! found
       ENDDO      !nf LOOP

          if(nd.eq.37.and.FDg(nd).eq.86) then
                       VD(nd)=2
                       FD1(nd)=80
                       VD1(nd)=3
                       FD2(nd)=81
                       VD2(nd)=1
          endif
          if(nd.eq.37.and.FDg(nd).eq.81) then
                       VD(nd)=1
                       FD1(nd)=80
                       VD1(nd)=3
                       FD2(nd)=86
                       VD2(nd)=2
          endif
           if(nd.eq.4.and.FDg(nd).eq.104) then
                        VD(nd)=3
                        FD1(nd)=100
                        VD1(nd)=1
                        FD2(nd)=101
                        VD2(nd)=2
             endif
            if(nd.eq.4.and.FDg(nd).eq.101) then
                        VD(nd)=2
                        FD1(nd)=100
                        VD1(nd)=1
                        FD2(nd)=104
                        VD2(nd)=3
             endif
            if(nd.eq.4.and.FDg(nd).eq.74) then
                         VD(nd)=3
                         FD1(nd)=66
                         VD1(nd)=1
                         FD2(nd)=67
                         VD2(nd)=2
             endif
            if(nd.eq.4.and.FDg(nd).eq.67) then
                         VD(nd)=2
                         FD1(nd)=66
                         VD1(nd)=1
                         FD2(nd)=74
                         VD2(nd)=3
             endif
          if(nd.eq.10.and.FDg(nd).eq.74) then
                         VD(nd)=3
                         FD1(nd)=66
                         VD1(nd)=1
                         FD2(nd)=67
                         VD2(nd)=2
           endif
           if(nd.eq.10.and.FDg(nd).eq.86) then
                         VD(nd)=2
                         FD1(nd)=80
                         VD1(nd)=3
                         FD2(nd)=81
                         VD2(nd)=1
           endif
           if(nd.eq.10.and.FDg(nd).eq.67) then
                         VD(nd)=2
                         FD1(nd)=66
                         VD1(nd)=1
                         FD2(nd)=74
                         VD2(nd)=3
           endif
           if(nd.eq.10.and.FDg(nd).eq.81) then
                         VD(nd)=1
                         FD1(nd)=80
                         VD1(nd)=3
                         FD2(nd)=86
                         VD2(nd)=2
           endif
          if(nd.eq.36.and.FDg(nd).eq.58) then
                         VD(nd)=3
                         FD1(nd)=52
                         VD1(nd)=1
                         FD2(nd)=53
                         VD2(nd)=2
           endif
           if(nd.eq.36.and.FDg(nd).eq.104) then
                         VD(nd)=3
                         FD1(nd)=100
                         VD1(nd)=1
                         FD2(nd)=101
                         VD2(nd)=2
           endif
           if(nd.eq.36.and.FDg(nd).eq.53) then
                          VD(nd)=2
                          FD1(nd)=52
                          VD1(nd)=1
                          FD2(nd)=58
                          VD2(nd)=3
           endif
           if(nd.eq.36.and.FDg(nd).eq.101) then
                         VD(nd)=2
                         FD1(nd)=100
                         VD1(nd)=1
                         FD2(nd)=104
                         VD2(nd)=3
           endif
           write(*,*) '-----',noface,'-------'
           write(*,*) nd,FDg(nd),VD(nd),FD1(nd),VD1(nd)
      ENDDO     ! nodata

!! HARDCODING

      write(*,*) 'next nodata loop-----'


C ##########################################################
C.....Update a field from Xi calculations      
!=====================================================
      DO nodata=1,NPLIST(0)
c     DO nd=1,NDT
!=====================================================

        np=NPLIST(nodata) !the node number in the list
        nd=NDP_INV(np) !the data point for node np

        write(*,*) 'nodata,np,nd,LD(nd)=',nodata,np,nd,LD(nd)

!       nf = FD(nd)
!       CALL XPXE_POINTS(from_field,NBJF(1,nf),NKJF,NPF(1,nf),NPNF,
!    &      nr,NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)

        IF (LD(nd).gt.0) THEN
          ne=LD(nd)
          nr=NRE(ne)

C         CALL XPXE_POINTS(from_field,NBJ(1,ne),NKJE(1,1,1,ne),
C    &      NPF(1,1),NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),SE(1,1,ne),
C    &      XA(1,1,ne),XE,XP,ERROR,*9999)
 
          write(*,*) 'IT(NBJ(1,ne))=',NIT(NBJ(1,ne)) 

          DO ni=1,NIT(NBJ(1,ne))
            XI(ni)=XID(ni,nd)
          ENDDO
          np=NDP(nd)
          write(*,*) 'nd,np=',nd,np

         
         IF(INLIST(np,NPLIST(1),NPLIST(0),N1LIST))THEN

          DO nj=1,NJT
          nvd = NVJP(nj,np)
          write(*,*) 'node versions=',nvd

          IF(nvd.eq.1) THEN !! ------nvd=1

           write(*,*) 'coming inside nvd=1 IF loop'
           nb=NBJ(nj,ne)
           IF(nb.GT.0) THEN
  
           CALL XPXE_POINTS(from_field,NBJF(1,nft),NKJF,NPF(1,nft),NPNF,
     &      nr,NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)

            DO nu=1,11
             Xmod(1,nu,nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     &       nu,XI,XE(1,nj))
            ENDDO !nu
           ENDIF !nb

          ELSE !!------for nvd>1

!          nb=NBJ(nj,ne)
           nbf=NBJF(1,FDg(nd))

           write(*,*) 'nbf for nvd>1',nbf

          IF(nbf.GT.0) THEN
!          DO nv=1,NVJP(nj,np) !for number of versions of geometry
 
           nv = VD(nd)
           nft = FDg(nd)

           write(*,*) 'nv nft = ', nv , nft
           write(*,*) 'NBJF(1,nft),NPF(1,nft)'
           write(*,*) NBJF(1,nft),NPF(1,nft)

          CALL XPXE_POINTS(from_field,NBJF(1,nft),NKJF,NPF(1,nft),NPNF,
     &      nr,NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)

           write(*,*) 'nut nbf=',NUT(nbf)
           write(*,*) 'IBT 1 1 nbf =',IBT(1,1,nbf)
           write(*,*) 'IDO 1 1 0 nbf =',IDO(1,1,0,nbf)
           write(*,*) 'INP 1 1 nbf =',INP(1,1,nbf)
           write(*,*) 'node = ',np
           DO nu=1, 6  !NUT(nbf)
            Xmod(nv,nu,nj)=
     &       PXIDERIV(IBT(1,1,nbf),IDO(1,1,0,nbf),INP(1,1,nbf),nbf,
     &      nu,XI,XE(1,nj))
           write(*,*) nv,nu,nj,'Xmod=',Xmod(nv,nu,nj)
           ENDDO !nu

           write(*,*) 'Next version interpolation'

           nv = VD1(nd)
           nft = FD1(nd)
 
           write(*,*) 'nv nft = ', nv, nft
           write(*,*) 'NBJF(1,nft),NPF(1,nft)'
           write(*,*) NBJF(1,nft),NPF(1,nft)

          CALL XPXE_POINTS(from_field,NBJF(1,nft),NKJF,NPF(1,nft),NPNF,
     &      nr,NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)

           write(*,*) 'nut nbf=',NUT(nbf)
           write(*,*) 'IBT 1 1 nbf =',IBT(1,1,nbf)
           write(*,*) 'IDO 1 1 0 nbf =',IDO(1,1,0,nbf)
           write(*,*) 'INP 1 1 nbf =',INP(1,1,nbf)
           write(*,*) 'node = ',np
           DO nu=1, 6  !NUT(nbf)
            Xmod(nv,nu,nj)=
     &      PXIDERIV(IBT(1,1,nbf),IDO(1,1,0,nbf),INP(1,1,nbf),nbf,
     &      nu,XI,XE(1,nj))
           write(*,*) nft,'----',nv,nu,nj,Xmod(nv,nu,nj)
           ENDDO !nu

           write(*,*) 'next version interpolation'

           nv = VD2(nd)
           nft = FD2(nd)

           write(*,*) 'nv nft = ', nv , nft
           write(*,*) 'NBJF(1,nft),NPF(1,nft)'
           write(*,*) NBJF(1,nft),NPF(1,nft)

          CALL XPXE_POINTS(from_field,NBJF(1,nft),NKJF,NPF(1,nft),NPNF,
     &      nr,NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)

           write(*,*) 'nut nbf=',NUT(nbf)
           write(*,*) 'IBT 1 1 nbf =',IBT(1,1,nbf)
           write(*,*) 'IDO 1 1 0 nbf =',IDO(1,1,0,nbf)
           write(*,*) 'INP 1 1 nbf =',INP(1,1,nbf)
           DO nu=1, 6  !NUT(nbf)
            Xmod(nv,nu,nj)=
     &      PXIDERIV(IBT(1,1,nbf),IDO(1,1,0,nbf),INP(1,1,nbf),nbf,
     &      nu,XI,XE(1,nj))
           ENDDO !nu
!          ENDDO !nv
          ENDIF

          ENDIF !nb
          ENDDO !nj

          write(*,*) 'normalizing derivatives'

         nvd = NVJP(nj,np)
         SUM=0.0d0
         DO ni=1,NIT(NBJ(1,ne))
         if (NVJP(ni,np).gt.1) then

         DO nvd=1,NVJP(ni,np)
           nured=1+ni*(1+ni)/2
           SUM=Xmod(nvd,nured,1)**2+Xmod(nvd,nured,2)**2
           IF(NJT.gt.2) THEN
           SUM=SUM+Xmod(nvd,nured,3)**2
           ENDIF !njt
           DSDXI(ni)=SQRT(SUM)

          IF(DABS(DSDXI(ni)).LT.DSQRT(DLAMCH('EPS'))) THEN
            DSDXI(ni)=1.0d0
            WRITE(OP_STRING,
     '        '('' >>Warning: Arc length derivative is zero'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          IF(NBI(nb).GE.5.AND.NBI(nb).LE.7) THEN
           IF(ni.eq.3) then
!               arc-length or ave. arc-length scale factors
            DSDXI(ni)=DSDXI(ni)
           endif
           ELSE
           ! unit scale factor              
           DSDXI(ni)=1.0d0
          ENDIF !nbi 
         ENDDO !nvd

         endif !------nvjp
         ENDDO !ni

         write(*,*) 'version=',NVJP(nj,np)
         DO nj=1,NJT
           DO nv=1,NVJP(nj,np) !for number of versions of geometry
            XP(1,nv,to_field(nj),np)=Xmod(nv,1,nj)
            XP(2,nv,to_field(nj),np)=Xmod(nv,2,nj)/DSDXI(1) !dXi1
            XP(3,nv,to_field(nj),np)=Xmod(nv,4,nj)/DSDXI(2) !dXi2
            XP(5,nv,to_field(nj),np)=Xmod(nv,7,nj)/DSDXI(3) !dXi3
!     XP(4,nv,to_field(nj),np)=Xmod(nv,6,nj)/(DSDXI(1)*DSDXI(2))!dXi1dXi2
            IF(nv.eq.1) THEN
!     XP(6,nv,to_field(nj),np)=Xmod(nv,9,nj)/(DSDXI(1)*DSDXI(3))!dXi1dXi3
!     XP(7,nv,to_field(nj),np)=Xmod(nv,10,nj)/(DSDXI(1)*DSDXI(3))!dXi2dXi3
            ENDIF
            write(*,*) 'Xmod2:dxi1=',DSDXI(1),XP(2,nv,to_field(nj),np)
            write(*,*) 'Xmod4:dxi2=',DSDXI(2),XP(3,nv,to_field(nj),np)
            write(*,*) 'Xmod7:dxi2=',DSDXI(3),XP(5,nv,to_field(nj),np)
           ENDDO !nv
         ENDDO !nj

        ENDIF
        ENDIF !LD
        stop

!======================================
      ENDDO !nd
!======================================

      CALC_XI = .TRUE.

      CALL EXITS('UPXI_POINTS')
      RETURN
 9999 CALL ERRORS('UPXI_POINTS',ERROR)
      CALL EXITS('UPXI_POINTS')
      RETURN 1

      END
