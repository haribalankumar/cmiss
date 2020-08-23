      SUBROUTINE DEXI_IN_FRUSTUM(FD,IBT,IDO,INP,LD,NBH,NBJ,NBJF,
     &  ND0,ND1,NFF,NFLIST,NEELEM,NHE,NKEF,NKHE,NKJE,NNF,NPF,
     &  NPNE,NPNF,NRE,NVHE,NVJE,NVJF,NW,nx,NXI,
     &  CURVCORRECT,SE,SF,SQ,XA,XE,XI,XI_1,XI_2,XI_3,XID,XP,ZA,
     &  ZD,ZP,DEFORM,NEW,SET_XI_1,SET_XI_2,SET_XI_3,
     &  SPECIFY,ERROR,*)

C#### Subroutine: DEXI_IN_FRUSTUM
C###  Description:
C###    DEXI_IN_FRUSTUM calculates the local coordinates of a
C###    projection of a data point on to the specified faces of
C###    a volume mesh. It uses a frustum inclusion technique to
C###    select the data points to project onto each face.
C**** Added by Glenn Ramsey February 2005
C**** Loosely based on DEXI_CLOSEST_FACE

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      

!     Parameter List
      INTEGER FD(NDM),IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     &  INP(NNM,NIM,NBFM),LD(NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     &  NBJF(NJM,NFM),ND0,ND1,NFF(6,NEM),NFLIST(0:NFM),
     &  NEELEM(0:NE_R_M,0:NRM),
     &  NHE(NEM,NXM),NKEF(0:4,16,6,NBFM),NKHE(NKM,NNM,NHM,NEM),
     &  NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),NPF(9,NFM),
     &  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NRE(NEM),
     &  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJF(NNM,NBFM,NJM),NW(NEM,3),nx,NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     &  SQ(NDM),XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XI_1,XI_2,XI_3,
     &  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),
     &  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL DEFORM,NEW,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY
!     Local Variables
      INTEGER nb,nbface,nd,ne,nf,nef,nface,ni,nj,NKJF(NKM,NNM,NJM),
     &  nn,nr,XIDIR(1:6),ADJELEM(0:4),ADJFACE(0:4),
     &  np0(1:4),np(1:4),mxi(-3:3),
     &  ne2,nf2,nef2,nface2,nb2,COPY,t,
     &  tf,NTRI(2),NTRIFACE
      REAL*8 SQND,d1,d2,d3,d4,XIone,XItwo,XIthree,
     &  p1(NJM), p2(NJM), p3(NJM),p4(NJM),ptmp(NJM),
     &  norm1(4),norm2(4),norm3(4),
     &  norm4(4),norm5(4),norm6(4),ntmp(4),
     &  frustum1(4,0:1), frustum2(4,0:1), frustum3(4,0:1),
     &  frustum4(4,0:1),
     &  NORMAL(4),NORMLIST(NJM,NFM,0:1),
     &  flip(0:1),face
      LOGICAL FOUND,DATAUSED(NDM)

      CALL ENTERS('DEXI_IN_FRUSTUM',*9999)
            
      !Map to change direction indices from NXI to NFF format
      mxi(-1)=1
      mxi(1)=2
      mxi(-2)=3
      mxi(2)=4
      mxi(-3)=5
      mxi(3)=6
      
      ! Initialise data related variables
      DO nd=ND0,ND1
        DATAUSED(nd)=.FALSE.
        FD(nd)=0 ! element face no. 
        LD(nd)=0 ! element no.
        SQ(nd)=0.0d0
        IF(NEW) THEN
          DO ni=1,3
            XID(ni,nd)=0.5d0
          ENDDO
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
      ENDDO
      
      DO nj=1,3
        DO ni=1,NFM
          DO t=0,1
            NORMLIST(nj,ni,t)=0.0d0
          ENDDO
        ENDDO
      ENDDO  

      ! for each face
      DO nface=1,NFLIST(0)
        nf=NFLIST(nface) ! the global id of this face
                
        !DOP=.TRUE.
        !find the 4 adjacent faces
                
        nef=NPF(8,nf) !local face #
        ne=NPF(6,nf) !element to which face is associated
        nr=NRE(ne) ! region for this element  
        !nb=NBJ(1,NEELEM(1,nr))
        nb=NBJ(1,ne)
        
        XIDIR(2)=NPF(1,nf)
        XIDIR(1)=-XIDIR(2) 
        XIDIR(4)=NPF(3,nf)
        XIDIR(3)=-XIDIR(4) 
        XIDIR(6)=NNF(1,nef,nb)
        XIDIR(5)=-XIDIR(6)
                
        ! get the global node numbers for this face
        nbface=NBJF(1,nf)
        DO nn=1,NNT(nbface)
          np0(nn)=NPNE(NNF(1+nn,nef,nb),nb,ne)
        ENDDO
        
        ! put the face of interest in the adjacent face list
        ! so it can be processed with the others
        ADJFACE(0)=nef ! ADJFACE has local face #'s of adjacent faces
        ADJELEM(0)=ne  ! ADJELEM global element #'s of adjacent elements
        
        DO ni=1,4
          ADJFACE(ni)=0
          ADJELEM(ni)=NXI(XIDIR(ni),1,ne)
          IF(ADJELEM(ni).NE.0) THEN
            IF(NFF(nef,ADJELEM(ni)).EQ.0) THEN
              ! this face is collapsed so choose the next one
              ADJFACE(ni) = mxi(XIDIR(ni))
            ELSE
              ! same dir as main face
              ADJFACE(ni) = nef
            ENDIF
          ELSE            
            !The adjacent face will be on this element in the
            !same dir for which the element is missing... unless
            !the face is collapsed.
            ADJELEM(ni) = ne
            ADJFACE(ni)=mxi(XIDIR(ni))
          ENDIF
        ENDDO                             
          
        IF(DOP) THEN
          CALL WRITES(IODI,OP_STRING,ERROR,*9999) ! blank line
          WRITE(OP_STRING,'(''Face '',I4,'' Xi dirs '',
     &      I3,I3,I3,I3,I3,I3)')
     &      nf,XIDIR(1),XIDIR(2),XIDIR(3),XIDIR(4),XIDIR(5),XIDIR(6)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF          
        
        ! for each face - ie the face of interest and its adjacent faces
        DO nface2 = 0,4
          nf2=NFF(ADJFACE(nface2), ADJELEM(nface2)) ! global face number
          
           ! reversal of the triangle normals, set to -1 to reverse
          flip(0)=1 ! triangle 0
          flip(1)=1 ! triangle 1
          
          IF(DOP) THEN
            WRITE(OP_STRING,'(''adjacent('',I1,'')='',I4)')nface2,nf2
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
           !if the face exists and is not already calculated
          IF((nf2.NE.0).AND.(NORMLIST(1,nf2,0).EQ.0)) THEN
          
            ! divide the face into two triangle segments, find and cache
            ! the normals for these triangles. The cache is the array NORMLIST.
            
            ! 1st triangle (index 0) has vertices (Xi1,Xi2) = (0,0) (0,1) (1,0) 
            ! 2nd triangle (index 1) has vertices (Xi1,Xi2) = (1,1) (1,0) (0,1)
                    
            ! get the global node numbers for the face
            ne2=NPF(6,nf2) !element to which face is associated
            nb2=NBJ(1,ne)
            nef2=ADJFACE(nface2)
            
            nbface=NBJF(1,nf)
            DO nn=1,NNT(nbface)                    
              np(nn)=NPNE(NNF(1+nn,nef2,nb2),nb2,ne2)
            ENDDO
            
            ! retrieve the geometric position of each node in this face
            DO nj=1,3
              p1(nj)=XP(1,1,nj,np(1))
              p2(nj)=XP(1,1,nj,np(2))
              p3(nj)=XP(1,1,nj,np(3))
              p4(nj)=XP(1,1,nj,np(4))
            ENDDO
            
            IF(DOP) THEN
              WRITE(OP_STRING,'(''p1 node'',I4,f8.3,f8.3,f8.3)')
     &          np(1), p1(1),p1(2),p1(3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''p2 node'',I4,f8.3,f8.3,f8.3)')
     &          np(2),p2(1),p2(2),p2(3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''p3 node'',I4,f8.3,f8.3,f8.3)')
     &          np(3),p3(1),p3(2),p3(3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''p4 node'',I4,f8.3,f8.3,f8.3)')
     &          np(4),p4(1),p4(2),p4(3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            
            ! Get the unit normal to each triangular segment of the face
            
            ! Because of the way that the mesh is structured
            ! the normals calculated will point in the +Xi dir
            ! so reverse them for faces in the -Xi dir.
            
            ! odd element face numbers at Xi=0, even at Xi=1
            ! so divide by 2 and use remainder
                         
            IF((nface2.EQ.0)) THEN ! main face
              ! flip main face (nef) if in -ve dir
              IF((MOD(nef,2).EQ.1)) THEN ! -Xi dir
                flip(0)=-1
                flip(1)=-1
              ELSE ! +Xi dir
                flip(0)=1
                flip(1)=1
              ENDIF
            ELSE ! adjacent faces
              ! adjacent faces: nef2
              IF((MOD(nef2,2).EQ.1)) THEN ! -Xi dir
                  ! don't flip xi on adjacent faces in +Xi dir
                  flip(0)=-1
                  flip(1)=-1                
              ELSE
                  flip(0)=1
                  flip(1)=1              
              ENDIF
            ENDIF
                  
            IF(DOP) THEN
              WRITE(OP_STRING,'(''adjacent('',I4,'')'',
     &          I4,'' dir '',I2,'' flip '',f3.0,'' '',f3.0)')
     &          nface2,nf2,NNF(1,nef2,nb2),flip(0),flip(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            ! If the face has a collapsed node then one of the triangles
            ! will be a line. In this case put a copy of the normal
            ! for the non-collapsed face in the cache for the collapsed
            ! face. This means that testing for a collapsed face can be
            ! avoided for some situations. 
            
            COPY=-1 ! flag used to specify if one normal needs to be a copy of the other
            
            ! lower triangle (index 0)
            IF( (np(1).EQ.np(2)).OR.(np(1).EQ.np(3))) THEN! triangle 0 is collapsed
              COPY=0
            ELSE
            
              IF(DOP) THEN
                WRITE(OP_STRING,'(''adjacent('',I4,'')'',
     &           I4,'' lower'')') nface2,nf2
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

              CALL PLANE_FROM_3_PTS(NORMAL,2,p1,p3,p2,ERROR,*9999)
              IF(DOP) THEN
                WRITE(OP_STRING,'(''triangle 0 normal '')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO nj=1,4
                NORMLIST(nj,nf2,0)=flip(0)*NORMAL(nj)
              
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''NORMAL('',I1,'')='',F8.3)')
     &              nj,NORMLIST(nj,nf2,0)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
     
              ENDDO
            ENDIF
            ! upper triangle (index 1)
            IF((np(4).EQ.np(3)).OR.((np(4).EQ.np(2)))) THEN! triangle 1 is collapsed
              COPY=1
            ELSE
            
              IF(DOP) THEN
                WRITE(OP_STRING,'(''adjacent('',I4,'')='',
     &            I4,'' upper'')') nface2,nf2
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              
              CALL PLANE_FROM_3_PTS(NORMAL,2,p4,p2,p3,ERROR,*9999)
              IF(DOP) THEN
                WRITE(OP_STRING,'(''triangle 1 normal '')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              
              DO nj=1,4
                NORMLIST(nj,nf2,1)=flip(1)*NORMAL(nj)
                
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''NORMAL('',I1,'')='',F8.3)')
     &              nj,NORMLIST(nj,nf2,1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                             
              ENDDO
            ENDIF
            IF(COPY.GT.-1) THEN
              DO nj=1,4
                NORMLIST(nj,nf2,COPY) = flip(COPY)*NORMAL(nj)
              ENDDO
            ENDIF
          ENDIF ! ((nf2.NE.0)
        ENDDO !nface2
        
        ! now that we have the adjacent faces and all their normals
        ! calculate the frustums for the 2 triangles on this face and
        ! find data points that are in the frustum
        
        ! retrieve the geometric position of each node 
        DO nj=1,3
          p1(nj)=XP(1,1,nj,np0(1)) ! reusing p*()
          p2(nj)=XP(1,1,nj,np0(2))
          p3(nj)=XP(1,1,nj,np0(3))
          p4(nj)=XP(1,1,nj,np0(4))
        ENDDO
        
          IF(DOP) THEN
            WRITE(OP_STRING,'(''face='',I4,'' nodes:'')') nf
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''p1 node'',I4,f8.3,f8.3,f8.3)')
     &        np0(1), p1(1),p1(2),p1(3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''p2 node'',I4,f8.3,f8.3,f8.3)')
     &        np0(2),p2(1),p2(2),p2(3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''p3 node'',I4,f8.3,f8.3,f8.3)')
     &        np0(3),p3(1),p3(2),p3(3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''p4 node'',I4,f8.3,f8.3,f8.3)')
     &        np0(4),p4(1),p4(2),p4(3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          
        ! If the face of interest is in the -Xi dir reverse
        ! the direction of all the frustum planes  
          
        flip(0)=1
        IF(MOD(nef,2).EQ.1) THEN !-Xi dir
          flip(0)=-1
          !flip1=-1
        ENDIF
          
        ! calculate average normals for each triangle vertex
        ! and define points on the frustum plane relative to triangle
        ! vertices so we can use these sets of 3 points to generate
        ! the frustum planes
        DO nj=1,3
          norm1(nj)=0
          norm2(nj)=0
          norm3(nj)=0
          norm4(nj)=0
          norm5(nj)=0
          norm6(nj)=0
          
          !Select the triangle adjacent to the face diagonal
          !to calculate the average normal.
          !If there are collapsed nodes then it will be on
          !another face and on the same side of the diagonal -
          ! ie will have the same triangle index
          face=nf
          t=1
          IF((np0(4).EQ.np0(2))
     &      .OR.(np0(4).EQ.np0(3))) THEN ! face is collapsed
            IF(ADJFACE(2).NE.0) THEN
              face=NFF(ADJFACE(2), ADJELEM(2))
            ENDIF
            t=0
          ENDIF
           
          norm1(nj)=(NORMLIST(nj,NINT(face),t)+NORMLIST(nj,nf,0))/2
     &      +p3(nj)
          
          face = NFF(ADJFACE(3), ADJELEM(3))
          IF((ADJFACE(3).NE.0).AND.(face.NE.0)) THEN 
            norm3(nj)=(NORMLIST(nj,NINT(face),1)+NORMLIST(nj,nf,0))/2
     &        +p2(nj)
          ENDIF
          
          face = NFF(ADJFACE(1), ADJELEM(1))  
          IF((ADJFACE(1).NE.0).AND.(face.NE.0)) THEN 
            norm2(nj)=(NORMLIST(nj,NINT(face),1)+NORMLIST(nj,nf,0))/2
     &        +p1(nj)
          ENDIF
          
          face=nf
          t=0
          IF( (np0(1).EQ.np0(2))
     &      .OR.(np0(1).EQ.np0(3))) THEN ! face is collapsed
            IF(ADJFACE(1).NE.0) THEN
              face=NFF(ADJFACE(1), ADJELEM(1))
            ENDIF
            t=1            
          ENDIF
          norm4(nj)=(NORMLIST(nj,NINT(face),t)+NORMLIST(nj,nf,1))/2
     &      +p2(nj)
          
          face = NFF(ADJFACE(4), ADJELEM(4))  
          IF((ADJFACE(4).NE.0).AND.(face.NE.0)) THEN
            norm6(nj)=(NORMLIST(nj,NINT(face),0)+NORMLIST(nj,nf,1))/2
     &        +p3(nj)
          ENDIF
          
          face = NFF(ADJFACE(2), ADJELEM(2))  
          IF((ADJFACE(2).NE.0).AND.(face.NE.0)) THEN
            face = NFF(ADJFACE(2), ADJELEM(2))  
            norm5(nj)=(NORMLIST(nj,NINT(face),0)+NORMLIST(nj,nf,1))/2
     &        +p4(nj)
          ENDIF
        ENDDO ! nj=1,3
          
        ! frustums on triangle 0 - lower 
        IF((np0(1).NE.np0(2)).AND.(np0(1).NE.np0(3))) THEN! triangle 0 not collapsed
          CALL PLANE_FROM_3_PTS(NORMAL,2,p3,norm1,p2,ERROR,*9999)
          DO nj=1,4
            frustum1(nj,0)=flip(0)*NORMAL(nj)
          ENDDO
          CALL PLANE_FROM_3_PTS(NORMAL,2,p1,norm2,p3,ERROR,*9999)
          DO nj=1,4
            frustum2(nj,0)=flip(0)*NORMAL(nj)
          ENDDO
          CALL PLANE_FROM_3_PTS(NORMAL,2,p2,norm3,p1,ERROR,*9999)
          DO nj=1,4
            frustum3(nj,0)=flip(0)*NORMAL(nj)
            ! if tri 1 is collapsed then this will be ignored
            frustum4(nj,1)=flip(0)*NORMAL(nj)
          ENDDO
        ELSE !triangle 0 was collapsed
          DO nj=1,3
            ptmp(nj)=p2(nj)+p3(nj)-p4(nj)
            ntmp(nj)=ptmp(nj)+NORMLIST(nj,nf,1)
          ENDDO
          CALL PLANE_FROM_3_PTS(NORMAL,2,p2,ntmp,ptmp,ERROR,*9999)
          DO nj=1,4
            frustum4(nj,1)=flip(0)*NORMAL(nj)
          ENDDO                  
        ENDIF
                
        ! frustums on triangle 1 - upper 
        IF((np0(4).NE.np0(3)).AND.((np0(4).NE.np0(2)))) THEN! triangle 1 not collapsed          
          CALL PLANE_FROM_3_PTS(NORMAL,2,p2,norm4,p3,ERROR,*9999)
          DO nj=1,4
            frustum1(nj,1)=flip(0)*NORMAL(nj)
          ENDDO
          CALL PLANE_FROM_3_PTS(NORMAL,2,p4,norm5,p2,ERROR,*9999)
          DO nj=1,4
            frustum2(nj,1)=flip(0)*NORMAL(nj)
          ENDDO          
          CALL PLANE_FROM_3_PTS(NORMAL,2,p3,norm6,p4,ERROR,*9999)
          DO nj=1,4
            frustum3(nj,1)=flip(0)*NORMAL(nj)
            frustum4(nj,0)=flip(0)*NORMAL(nj)
          ENDDO
        ELSE
            DO nj=1,3
              ptmp(nj)=p3(nj)+p2(nj)-p1(nj)
              ntmp(nj)=ptmp(nj)+NORMLIST(nj,nf,0)
            ENDDO
            CALL PLANE_FROM_3_PTS(NORMAL,2,p3,ntmp,ptmp,ERROR,*9999)
            DO nj=1,4
              frustum4(nj,0)=flip(0)*NORMAL(nj)
            ENDDO
        ENDIF
        
        IF(DOP) THEN
          DO t=0,1
            WRITE(OP_STRING,'(''frustum1-'',I1)') t
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nj=1,4
              WRITE(OP_STRING,'(''NORMAL('',I1,'')='',F8.3)')
     &        nj,frustum1(nj,t)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO                
            WRITE(OP_STRING,'(''frustum2-'',I1)') t
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nj=1,4
              WRITE(OP_STRING,'(''NORMAL('',I1,'')='',F8.3)')
     &        nj,frustum2(nj,t)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            WRITE(OP_STRING,'(''frustum3-'',I1)') t
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nj=1,4
              WRITE(OP_STRING,'(''NORMAL('',I1,'')='',F8.3)')
     &        nj,frustum3(nj,t)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO                
            WRITE(OP_STRING,'(''frustum4-'',I1)') t
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nj=1,4
              WRITE(OP_STRING,'(''NORMAL('',I1,'')='',F8.3)')
     &        nj,frustum4(nj,t)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO               
        ENDIF        
        
        NTRIFACE=2
        NTRI(1)=0
        NTRI(2)=1
        IF( (np0(1).EQ.np0(2)).OR.(np0(1).EQ.np0(3)) ) THEN! triangle 0 is collapsed
          NTRIFACE=1
          NTRI(1)=1
        ELSEIF ((np0(4).EQ.np0(3)).OR.((np0(4).EQ.np0(2)))) THEN! triangle 1 is collapsed
          NTRIFACE=1
          NTRI(1)=0
        ENDIF
        
        !DOP=.FALSE.
                   
        ! Test to see if the data point lies in the frustum.
        ! Do this by calculating the distance from the data point to
        ! each of the frustum planes (by substituting the point coords
        ! into the plane equation). If the distance is -ve then
        ! the point lies on the 'wrong' side of the plane.
        DO tf=1,NTRIFACE
          t=NTRI(tf)
          DO nd=ND0,ND1
            IF(.NOT.DATAUSED(nd)) THEN
              d1=frustum1(4,t)
              DO nj=1,3
                d1=d1+frustum1(nj,t)*ZD(nj,nd)
              ENDDO
              IF(d1.GE.0.0d0) THEN
                d2=frustum2(4,t)
                DO nj=1,3
                  d2=d2+frustum2(nj,t)*ZD(nj,nd)
                ENDDO
                IF(d2.GE.0.0d0) THEN
                  d3=frustum3(4,t)
                  DO nj=1,3
                    d3=d3+frustum3(nj,t)*ZD(nj,nd)
                  ENDDO
                  IF(d3.GE.0.0d0) THEN
                    d4=frustum4(4,t)
                    DO nj=1,3
                      d4=d4+frustum4(nj,t)*ZD(nj,nd)
                    ENDDO
                    IF(d4.GE.0.0d0) THEN                 
                    ! point is in the +ve halfspace of this plane too
                      DATAUSED(nd)=.TRUE.
                      
                      IF(DOP) THEN                    
                        WRITE(OP_STRING,'(''FOUND gl. face='',I5,
     &                    '' nd='',I5,'' tri= '',I1,'' d='',
     &                    f12.3,f12.3,f12.3,f12.3)')
     &                    nf,nd,t,d1,d2,d3,d4
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                                
                      IF (DEFORM) THEN
C JWF- Call to ZPZE for deformed option has temporarily
C used NBJF for NBHF, NKJF for NKHF, and NVJF for NVHF
                        CALL CALC_FACE_INFORMATION_DEP(NBH(1,1,ne),
     &                    NBJF(1,nf),nef,NHE(ne,nx),NKHE(1,1,1,ne),
     &                    NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,
     &                    NVHE(1,1,1,ne),NVJF,nx,SE(1,1,ne),
     &                    SF,ERROR,*9999)
                        CALL ZPZE(NBJF(1,nf),1,NHE(ne,nx),NKJF,
     &                    NPF(1,nf),NPNF,nr,NVJF,NW(ne,1),nx,
     &                    CURVCORRECT(1,1,1,ne),SF,
     &                    ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
                      ELSE
                        CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
     &                    NBJF(1,nf),nef,NKJE(1,1,1,ne),NKEF,NKJF,
     &                    NNF,NPNE(1,1,ne),NPNF,nr,NVJE(1,1,1,ne),
     &                    NVJF,SE(1,1,ne),SF,ERROR,*9999)
                        CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,
     &                    nr,NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)
                      ENDIF
          
                      FOUND=.TRUE.
                      SQND=0.0d0
                      
                      CALL PROJ_ORTHOG(IBT,IDO,INP,NBJF(1,nf),SQND,
     &                  XE,XI,ZD(1,nd),FOUND,ERROR,*9999)
                      FD(nd)=nef
                      LD(nd)=ne
                      SQ(ND)=SQND
                      XIone=XI(1)
                      XItwo=XI(2)
                      
C       Determine the value of the 3rd xi direction.
  
                      ! odd element face numbers at Xi=0, even at Xi=1
                      ! so divide by 2 and use remainder
                      IF(MOD(nef,2).EQ.1) THEN
                        Xithree=0.0d0
                      ELSE
                        Xithree=1.0d0
                      ENDIF
                      
                      XID(NPF(1,nf),nd)=XIone
                      XID(NPF(3,nf),nd)=XItwo
                      XID(NNF(1,nef,nb),nd)=XIthree                    
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF !.NOT.DATAUSED
          ENDDO !nd
        ENDDO  !t
      ENDDO !nface
            
      CALL EXITS('DEXI_IN_FRUSTUM')
      RETURN
 9999 CALL ERRORS('DEXI_IN_FRUSTUM',ERROR)
      CALL EXITS('DEXI_IN_FRUSTUM')
      RETURN 1
      END
