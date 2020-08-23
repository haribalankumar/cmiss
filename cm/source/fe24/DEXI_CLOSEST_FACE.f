      SUBROUTINE DEXI_CLOSEST_FACE(FD,IBT,IDO,INP,LD,NBH,NBHF,NBJ,
     &  NBJF,ND0,ND1,NFF,NFLIST,NHE,NKEF,
     &  NKHE,NKJE,NNF,NPF,NPNE,NPNF,NRE,NSEED,NSEARCH,
     &  NVHE,NVJE,
     &  NVJF,nx,NXI,SE,SF,SQ,XA,XE,XI,
     &  XI_1,XI_2,XI_3,XID,XP,ZD,ZE,ZP,DEFORM,NEW,ORTHO,
     &  SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,CROSS,ERROR,*)

C#### Subroutine: DEXI_CLOSEST_FACE
C###  Description:
C###    DEXI_CLOSEST_FACE calculates the local coordinates of a
C###    projection of a data point on to the closest face of a
C###    volume mesh.
C**** Added by JWF 11/12/01
C**** Modified by PM 10/09/02
C**** Modified by LKC Mar 05
C**** Modified by GDR Mar 05 - added element crossing search
C**** Modified by GDR Apr 05 - changed initial face selection

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      
!     Parameter List
      INTEGER FD(NDM),IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     &  INP(NNM,NIM,NBFM),LD(NDM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     &  NBJ(NJM,NEM),NBJF(NJM,NFM),ND0,ND1,NFF(6,NEM),NFLIST(0:NFM),
     &  NHE(NEM,NXM),NKEF(0:4,16,6,NBFM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),NPF(9,NFM),
     &  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NRE(NEM),NSEED,NSEARCH,
     &  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJF(NNM,NBFM,NJM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     &  SQ(NDM),XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XI_1,XI_2,XI_3,
     &  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     &  ZD(NJM,NDM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL DEFORM,NEW,ORTHO,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,
     &  CROSS
!     Local Variables
      INTEGER nb,nbface,NCLOSFACE(0:NFM),ncf,ncfm,
     &  nd,ne,nf,nef,nh,ni,ni2,nj,njj,nk,
     &  NKHF(NKM,NNM,NHM),NKJF(NKM,NNM,NJM),
     &  nn,noface,nr,ns,N1LIST,IT,ITMAX,NITB,
     &  nflast,neadj,nelast,nefadj,nflast2,
     &  mxif(-3:3),xidir(1:6),NVHF(NNM,NBFM,NHM)
      REAL*8 SQND,TEMP,SQNDlast,XIlast(3),FDSQ(0:NFM)
      LOGICAL FOUND,INLIST,reject,LDOP
!     Functions
      REAL*8 PXI
      
      PARAMETER(ITMAX=20)

      CALL ENTERS('DEXI_CLOSEST_FACE',*9999)
               
C!!! LKC 1-MAR-2005 I don't think this loop is actually necessary if the data
C!!!   is allowed to cross element boundaries
C GDR 21Mar2005 I have implemented adjacent face searching so this is not
C really needed. However there are at least 2 cases (that I know of) where
C the boundary crossing will fail and the wrong face will be found and
C therefore 2 initial guesses are still needed until these are fixed.
C Both cases are illustrated in example 21h. First is when there is a collapsed
C face, this is not fatal and can be fixed by accounting for this in the
C boundary crossing algorithm. I have attempted this. Second case is when the
C initial guess has an orthogonal projection onto the wrong face because
C of an acute angle between faces. You can see the effect by setting ncfm
C to 1 and looking at the initial projections in 21h. This cannot (AFAIK)
C be detected by the boundary crossing algorithm. However the frustum approach
C to making an initial guess will probably work for this case. I will add
C boundary crossing to that one day so that the only difference between frustum
C and closest_face will be the way that the initial guess is made.
 
      ! GDR 7Apr05
      ! ncfm determines the number of faces which will be used as initial
      ! searching locations based on the distance from the data point to
      ! a seed point on the face. The default is 2 (set in DEXI.f).
      ! All the faces can be searched if required. This takes more time
      ! but for some problems might be needed.
      IF(NSEARCH.LE.NFLIST(0)) THEN
        ncfm=NSEARCH
      ELSE
        ncfm=NFLIST(0)
      ENDIF
                  
      !Map to change direction indices from NXI to NFF format
      mxif(-1)=1
      mxif(1)=2
      mxif(-2)=3
      mxif(2)=4
      mxif(-3)=5
      mxif(3)=6
      
      ! for debugging this subroutine only set LDOP to .TRUE.
      LDOP=DOP
      !LDOP=.TRUE.
      
      DO nd=ND0,ND1
        SQ(nd)=0.0d0
        IF(NEW) THEN
          LD(nd)=0
          FD(nd)=0
        ENDIF
      ENDDO

      DO nd=ND0,ND1

C     Find the distance to each face based on the distance to a seed point

        ! init arrays
        DO ncf=0,NFM
          NCLOSFACE(ncf)=0
          FDSQ(ncf)=0.0d0
        ENDDO

C LKC 1-MAR-2005 Putting section of code to skip the initial face
C search if using OLD xi positions.
        IF(NEW) THEN
          IF(NFLIST(0).EQ.1) THEN !if only one face is specified in the list.
            ncfm=1
          ENDIF

C!!! LKC 1-MAR-2005 This code is crap. It is trying to find the 2 closest faces
C (the reason why is another matter). Why search the full list twice?
C Should really only build up the list of closest faces as we go and only loop
C through the faces once.

C!!! GDR 21Mar05 For the reason 2 faces are searched see my first comment.
          
C GDR 7Apr05 Changed initial face selection so that now it only does
C one loop and stores the distance to the closest seed point on each face.
C Then it sorts the list to determine which faces are closest.
C Any number of these closest faces can be searched for the closest
C projection and the search will start with the closest based on the
C initial seed point. 

          !DO ncf=1,ncfm            
            DO noface=1,NFLIST(0)
              nf=NFLIST(noface)
!               IF((nf.EQ.NCLOSFACE(ncf-1)).AND.(NCLOSFACE(ncf-1).NE.0))
!      &          THEN
!                 GOTO 100
!               ENDIF
              
              ne=NPF(6,nf)
              nef=NPF(8,nf)
              nr=NRE(ne)

C VJ adding comment: If user asks to project data point onto deformed surface interpolate over deformed surfaces
              IF (DEFORM) THEN
                CALL CALC_FACE_INFORMATION_DEP(NBH(1,1,ne),
     &            NBHF(1,1,nf),nef,NHE(ne,nx),NKHE(1,1,1,ne),NKEF,
     &            NKHF,NNF,NPNE(1,1,ne),NPNF,NVHE(1,1,1,ne),NVHF,
     &            nx,SE(1,1,ne),SF,ERROR,*9999)
C news VJ 3Feb2005 setting up arrays to be used in ZPZE_FACE
C to optimise for material parameters using geometric least squares on deformed faces
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  nh=NH_LOC(njj,nx)
                  nbface=NBHF(nh,1,nf)
C                  NBHF(nh,1,nf)=nb
                  DO nn=1,NNT(nbface)
C                    NVHF(nn,nb,nh)=NVJF(nn,nb,nj)
                    DO nk=1,NKT(nn,nbface)
C                      NKHF(nk,nn,nh)=NKJF(nk,nn,nj)
                    ENDDO
                  ENDDO
                ENDDO

                CALL ZPZE_FACE(nf,NBHF,1,NKHF,NPNF,NVHF,nx,SF,ZE,ZP,
     &            ERROR,*9999)
C copying entries in ZE into XE for performing data projections.
                 DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                   nj=NJ_LOC(NJL_GEOM,njj,nr)
                   nh=NH_LOC(njj,nx)
                   nbface=NBHF(nh,1,nf)      
                   DO ns=1,NST(nbface)
                     XE(ns,nj)=ZE(ns,nh)
                   ENDDO
                 ENDDO
C newe VJ
C               CALL ZPZE(NBJF(1,nf),1,NHE(ne,nx),NKJF,NPF(1,nf),
C     '           NPNF,nr,NVJF,NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SF,
C     '           ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
              ELSE
                CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),
     &            nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
     &            nr,NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
                CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,NVJF,
     &            SF,XA(1,1,1),XE,XP,ERROR,*9999)
              ENDIF

C LKC 2-MAR-2005 Search for X times in each xi direction instead of
C just the cell centres [0.5,0.5].
C
C              XI(1)=0.5d0
C              XI(2)=0.5d0
C              SQND=0.0d0
C              DO nj=1,3 !max nj=3, as face fitting is always with volume elements
C                nb=NBJF(nj,nf)
C                TEMP=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     &            XE(1,nj))-ZD(nj,nd)
C                SQND=SQND+TEMP*TEMP
C              ENDDO !nj             
C              IF((NCLOSFACE(ncf).EQ.0).OR.(SQND.LT.SQ(nd))) THEN
C                NCLOSFACE(ncf)=nf
C                SQ(nd)=SQND
C              ENDIF

              DO ni=1,NSEED
                XI(1)=1.d0/DBLE(NSEED+1)*DBLE(ni)
                DO ni2=1,NSEED
                  XI(2)=1.d0/DBLE(NSEED+1)*DBLE(ni2)

                  SQND=0.0d0
                  DO nj=1,3 !max nj=3, as face fitting is always with volume elements
                    nb=NBJF(nj,nf)
                    TEMP=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     &                XI,XE(1,nj))-ZD(nj,nd)
                    SQND=SQND+TEMP*TEMP
                  ENDDO !nj
                  
!                   IF((NCLOSFACE(ncf).EQ.0).OR.(SQND.LT.SQ(nd))) THEN
!                     NCLOSFACE(ncf)=nf
!                     SQ(nd)=SQND
!                   ENDIF
                  IF((NCLOSFACE(noface).EQ.0).OR.(SQND.LT.SQ(nd))) THEN
                    FDSQ(noface)=SQND
                    !NCLOSFACE(ncf)=nf
                    NCLOSFACE(noface)=nf
                    SQ(nd)=SQND
                  ENDIF
                ENDDO !ni
              ENDDO !ni2
                           
! 100        ENDDO !noface (face counter)
            ENDDO !noface (face counter)
          !ENDDO !ncf (list of closest faces)
          
          ! Now sort the list into increasing order
          CALL RSORT(NFLIST(0),FDSQ(1),NCLOSFACE(1))
          
          ! NCLOSFACE() now contains the closest faces to this data point
          ! based on distance to a seed point, sorted in order. 
 
          ! initialise SQ() with closest distance to face based on seed point
          SQ(nd)=FDSQ(1)
          
        ELSE !old

C LKC 1-MAR-2005 Using OLD xi positions as in initial guess.
C the ipxi file should have contained ne, nf_local, xi1, xi2, xi3          

          ncfm=1  
          ne=LD(nd)
          nef=FD(nd) !local face number
          NCLOSFACE(1)=NFF(nef,ne) !global face number & starting search face
          
        ENDIF !not NEW

C*** There is now a list of (currently 2) faces which will be used to
C        search for projections.

        IF(NEW) THEN !otherwise use existing xi positions (OLD)
          DO ni=1,3
            XID(ni,nd)=0.5d0
          ENDDO
        ENDIF !new

        DO noface=1,ncfm
          nf=NCLOSFACE(noface) ! starting face
          
          ne=NPF(6,nf)
          IT=0
          FOUND=.FALSE.
          nflast2=-1
          nflast=-1
          
          IF(LDOP) THEN
            WRITE(OP_STRING,'(''ND '',I4,'' Start face '',
     &        I4)') nd,nf
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
            
          IF(SPECIFY) THEN ! specifies starting xi values
            IF(SET_XI_1) THEN
              XI(1)=XI_1
            ENDIF
            IF(SET_XI_2) THEN
              XI(2)=XI_2
            ENDIF
            IF(SET_XI_3) THEN
              XI(3)=XI_3
            ENDIF
          ELSE
             DO ni=1,3
               XI(ni)=XID(ni,nd)
             ENDDO
          ENDIF

          DO WHILE(.NOT.FOUND.AND.IT.LT.ITMAX.AND.nf.NE.0)
            IT=IT+1            
            nef=NPF(8,nf)
            nr=NRE(ne)
            nb=NBJ(1,ne)
            
            xidir(2)=NPF(1,nf)
            xidir(1)=-xidir(2) 
            xidir(4)=NPF(3,nf)
            xidir(3)=-xidir(4) 
            xidir(6)=NNF(1,nef,nb)
            xidir(5)=-xidir(6)
            
            IF(DEFORM) THEN
              CALL CALC_FACE_INFORMATION_DEP(NBH(1,1,ne),
     &          NBHF(1,1,nf),nef,NHE(ne,nx),NKHE(1,1,1,ne),NKEF,
     &          NKHF,NNF,NPNE(1,1,ne),NPNF,NVHE(1,1,1,ne),NVHF,
     &          nx,SE(1,1,ne),SF,ERROR,*9999)
C news VJ 3Feb2005 setting up arrays to be used in ZPZE_FACE
C to optimise for material parameters using geometric least squares on deformed faces
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                nh=NH_LOC(njj,nx)
                nbface=NBHF(nh,1,nf)
C                NBHF(nh,1,nf)=nb
                DO nn=1,NNT(nbface)
C                 NVHF(nn,nb,nh)=NVJF(nn,nb,nj)
                  DO nk=1,NKT(nn,nbface)
C                   NKHF(nk,nn,nh)=NKJF(nk,nn,nj)
                  ENDDO
                ENDDO
              ENDDO

              CALL ZPZE_FACE(nf,NBHF,1,NKHF,NPNF,NVHF,nx,SF,ZE,ZP,
     &          ERROR,*9999)
C copying entries in ZE into XE for performing data projections.
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                nh=NH_LOC(njj,nx)
                nbface=NBHF(nh,1,nf)      
                DO ns=1,NST(nbface)
                  XE(ns,nj)=ZE(ns,nh)
                ENDDO
              ENDDO
C newe VJ
C             CALL ZPZE(NBJF(1,nf),1,NHE(ne,nx),NKJF,NPF(1,nf),
C     '         NPNF,nr,NVJF,NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SF,
C     '         ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
            ELSE
              CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),
     &          nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
     &          nr,NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
              CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,NVJF,
     &          SF,XA(1,1,1),XE,XP,ERROR,*9999)
            ENDIF
              
            FOUND=.TRUE.
            XIlast(1)=XI(1)
            XIlast(2)=XI(2)
            SQNDlast=SQND
            SQND=0.0d0
                                                      
            CALL PROJ_ORTHOG(IBT,IDO,INP,NBJF(1,nf),SQND,
     &        XE,XI,ZD(1,nd),FOUND,ERROR,*9999)
            
            ! Disable face searching if the CROSS parameter is false (the default)
            IF(CROSS) THEN
              FOUND=.FALSE.
            ELSE
              FOUND=.TRUE.
            ENDIF
              
              IF(LDOP) THEN
                WRITE(OP_STRING,'(''ND '',I4,'' XI '',
     &          F11.3,F11.3,'' SQD '',E11.3,I4,I4,I4)')
     &          nd,XI(1),XI(2),SQND,ne,nef,nf
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              
            
            IF(.NOT.SPECIFY) THEN !try neighbouring face
              nflast2=nflast
              nflast=nf
              nelast=ne
              nefadj=nef
              NITB=NIT(NBJF(1,nf))
                           
              DO ni=1,NITB
                IF(.NOT.FOUND) THEN
                
C GDR 21Mar05 If the projection is on the face edge ie Xi=0.0 or 1.0
C Then find the adjacent face in that direction.
C Note that the face Xi indices may not match the element
C Xi indices so we need to map them to the correct dirs.
C That is what the mxif and xidir arrays are for.
C??? Is there a better way to do this mapping?

C If there is no orthogonal projection and the data point is
C bouncing back and forth between 2 faces then choose the closest
C face. This is what the *last variables are for.

                  IF(LDOP) THEN
                    WRITE(OP_STRING,'(''ND '',I4,'' Element '',
     &                I4,'' Face '',I4,'' IT '',I2,'' XI['',
     &                I1,''<=>''I1,'']='',F11.3,'' SQD '',E11.3)')
     &                nd,ne,nf,IT,ni,xidir(mxif(ni)),XI(ni),SQND
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  
                  IF(XI(ni).LE.0.0d0
     &              .AND.SQND.GT.0.0d0) THEN

                    neadj=NXI(xidir(mxif(-ni)),1,ne)
                    IF(neadj.GT.0) THEN
                      !nefadj=nef
                      ! account for collapsed face
                      IF(NFF(nef,neadj).EQ.0) THEN
                        ! this face is collapsed so choose the next one
                        nefadj = mxif(xidir(mxif(-ni)))
                      ELSE
                        ! same dir as main face
                        nefadj = nef
                      ENDIF
                      ne=neadj
                    ELSE
                      neadj=ne
                      nefadj=mxif(xidir(mxif(-ni)))
                    ENDIF
                    nf=NFF(nefadj,ne)
                    IF(nf.NE.0.AND.neadj.GT.0
     &                .AND.MOD(nefadj,2).EQ.MOD(nef,2)
     &                .AND.MOD(nefadj,2).EQ.0) THEN
                      XI(ni)=1.0d0
                    ENDIF
                    IF(LDOP) THEN
                      WRITE(OP_STRING,'(''Hit edge xi=0 nefadj:'',
     &                  I4,'' nfadj '',I4)') nefadj,nf
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF   
                  ELSE IF(XI(ni).GE.1.0d0
     &                  .AND.SQND.GT.0.0d0) THEN
     
                    neadj=NXI(xidir(mxif(ni)),1,ne)
                    IF(neadj.GT.0) THEN
                      !nefadj=nef
                      ! account for collapsed face
                      IF(NFF(nef,neadj).EQ.0) THEN
                        ! this face is collapsed so choose the next one
                        nefadj = mxif(xidir(mxif(ni)))
                      ELSE
                        ! same dir as main face
                        nefadj = nef
                      ENDIF
                      ne=neadj
                    ELSE
                      neadj=ne
                      nefadj=mxif(xidir(mxif(ni)))
                    ENDIF
                    nf=NFF(nefadj,ne)
                    IF(nf.NE.0.AND.neadj.GT.0
     &                .AND.MOD(nefadj,2).EQ.MOD(nef,2)
     &                .AND.MOD(nefadj,2).EQ.1) THEN
                      XI(ni)=0.0d0
                    ENDIF
                    IF(LDOP) THEN
                      WRITE(OP_STRING,'(''Hit edge xi=1 nefadj:'',
     &                  I4,'' nfadj: '',I4)') nefadj,nf
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(SQND.EQ.0.0d0) THEN
                    FOUND=.TRUE.
                  ENDIF !in bounds
                                            
                ENDIF ! FOUND
              ENDDO
              
              ! check if the adjacent face is in the list
              IF(.NOT.FOUND.AND.
     &            .NOT.INLIST(nf,NFLIST(1),NFLIST(0),N1LIST)) THEN
                IF(LDOP) THEN
                  WRITE(OP_STRING,'(''Face '', I4,
     &            '' not in list, using face '',I4)') nf, nflast
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                nf=nflast
                ne=nelast
                nef=NPF(8,nf)
                !SQND=SQNDlast
                !XI(1)=XIlast(1)
                !XI(2)=XIlast(1)
              ENDIF
              
              FOUND=(nf.EQ.nflast.OR.nf.EQ.nflast2)              
            ENDIF !not specify
            
            ! if neither face has an orthogonal projection
            ! then choose the closest
            IF(nf.EQ.nflast2) THEN
              IF(LDOP) THEN
                WRITE(OP_STRING,'(''No orthogonal proj for face:'',
     &            I4,'' SQD '',F11.3,'' nflast:'',I4,'' SQD '',
     &            F11.3)') nf,SQND,nflast,SQNDlast
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(SQNDlast.LT.SQND) THEN
                nf=nflast
                ne=nelast
                SQND=SQNDlast
              ELSE
                XI(1)=XIlast(1)
                XI(2)=XIlast(2)                
              ENDIF
              IF(XI(1).EQ.0.0d0) XI(1)=1.0d0
              IF(XI(1).EQ.1.0d0) XI(1)=0.0d0
              IF(XI(2).EQ.0.0d0) XI(2)=1.0d0
              IF(XI(2).EQ.1.0d0) XI(2)=0.0d0
              
              nef=NPF(8,nf)
              
              IF(LDOP) THEN
                WRITE(OP_STRING,'(''using face:'',
     &            I4)') nf
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            
            ! Reject points with a non-orthogonal projection
            ! if the ortho flag is set.
            reject=.FALSE.
            IF(ORTHO.AND.SQND.NE.0.0d0
     &         .AND.((XI(1).EQ.0.0d0.OR.XI(1).EQ.1.0d0)
     &         .OR.(XI(2).EQ.0.0d0.OR.XI(2).EQ.1.0d0))) THEN
              reject=.TRUE.
            ENDIF
            
            IF((SQND.LE.SQ(nd)).AND.FOUND.AND..NOT.reject) THEN
            !IF(FOUND.AND..NOT.reject) THEN

              SQ(nd)=SQND
              FD(nd)=nef
              LD(nd)=ne
              !XIone=XI(1)
              !XItwo=XI(2)
              IF(MOD(nef,2).EQ.1) THEN
                XI(3)=0.0d0
              ELSE
                XI(3)=1.0d0
              ENDIF
              
              XID(NPF(1,nf),nd)=XI(1)
              XID(NPF(3,nf),nd)=XI(2)
              XID(NNF(1,nef,nb),nd)=XI(3)
              
              IF(LDOP) THEN
                WRITE(OP_STRING,'(''ND '',I4,'' using XI '',
     &          F11.3,F11.3,F11.3,'' SQD '',E11.3, '' ne:'',I4,
     &          '' nef:'',I4,'' nf:'',I4)')
     &          nd,XI(1),XI(2),XI(3),SQND,ne,nef,nf
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

            ENDIF
            
!            ELSE
!               SQ(nd)=0.0d0
!               FD(nd)=0
!               LD(nd)=0
!            ENDIF !FOUND
     
!           IF((SQND.LT.SQ(nd)).OR.(noface.EQ.1)) THEN
!             SQ(nd)=SQND
!             FD(nd)=nef
!             LD(nd)=ne
!             XIone=XI(1)
!             XItwo=XI(2)
!           ENDIF
        
          ENDDO !not found
        ENDDO !noface
        
        IF(reject) THEN
              WRITE(OP_STRING,'(''>>Warning no orthogonal '',
     &          ''projection found for data point '',I6)') nd
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF    

! C       Determine the value of the 3rd (final) xi direction.
!         ne=LD(nd)
!         nef=FD(nd)
!         nf=NFF(nef,ne)
!         nr=NRE(ne)
!         CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),
!      &    nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
!      &    nr,NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
! 
!         nbface=NBJF(1,nf)
!         DO nn=1,NNT(nbface)
!           A_LIST(nn)=NPNF(nn,nbface)
!         ENDDO
! 
!         nb=NBJ(1,ne)
!         IF(INLIST(NPNE(1,nb,ne),A_LIST,NNT(nbface),N1LIST)) THEN
!           Xithree=0.0d0
!         ELSE
!           Xithree=1.0d0
!         ENDIF

!         XID(NPF(1,nf),nd)=XIone
!         XID(NPF(3,nf),nd)=XItwo
!         XID(NNF(1,nef,nb),nd)=XIthree

C??? GDR wouldn't it be easier to do this to determine the 3rd Xi dir?                
!         ! odd element face numbers (nef) at Xi=0, even at Xi=1
!         ! so divide nef by 2 and use remainder
!         IF(MOD(nef,2).EQ.1) THEN
!           Xithree=0.0d0
!         ELSE
!           Xithree=1.0d0
!         ENDIF

      ENDDO !nd
      
      CALL EXITS('DEXI_CLOSEST_FACE')
      RETURN
 9999 CALL ERRORS('DEXI_CLOSEST_FACE',ERROR)
      CALL EXITS('DEXI_CLOSEST_FACE')
      RETURN 1
      END
