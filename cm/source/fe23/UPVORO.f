      SUBROUTINE UPVORO(NENFVC,NFVC,NODENVC,NPNODE,NRLIST,NVCNODE,VC,
     '  VC_INIT,XNFV,XP,ZA,STRING,ERROR,*)

C#### Subroutine: UPVORO
C###  Description:
C###    UPVORO updates Voronoi cell volumes.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'voro00.inc'
      INCLUDE 'voro01.inc'
!     Parameter list
      INTEGER NENFVC(0:NFVCM,NFVM),NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NVCNODE(2,NP_R_M)
      REAL*8 VC(0:NVCM),VC_INIT(2,NVCM),XNFV(-(NJM+1):NJM,NFVM),
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local variables
      INTEGER IBEG,IEND,cnonode,cnp,ne,nf,nfvl,nj,njj,nj1,nj_vol,
     '  nonode,np,npl,nr,nvc,N3CO,VERT,VERTTOTAL
      REAL*8 CENTROID(3),FACEAREA,PROJECTION(3),VC_field,
     '  VC_sum,VERTICES(3,60),XPROD(3)
      LOGICAL ALL_REGIONS
!     Functions
      INTEGER IFROMC
      REAL*8 DDOT,DOTPROD
      LOGICAL CBBREV

      CALL ENTERS('UPVORO',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        
C---------------------------------------------------------------------
        
C#### Command: FEM update voronoi
C###  Parameter:      <field FIELD#[1]>
C###  Specify the field number for the cell volumes.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region number the fields and cells are in.
        
        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<field FIELD#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        
C---------------------------------------------------------------------
        
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPVORO',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1) !only for a single region
        
        IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
          nj1=IFROMC(CO(N3CO+1))
        ELSE
          nj1=1
        ENDIF
        nj_vol=NJ_LOC(NJL_FIEL,nj1,nr)
        CALL ASSERT(nj_vol.GT.0,'>>Non-existent field',ERROR,*9999)
      
        VC(0)=0
        VC_field=0
c        nfv=1
        
        DO nvc=1,NVCT !each Voronoi cell
          nonode=NODENVC(nvc) !local node # at cell centre
          np=NPNODE(nonode,nr) !global node # at cell centre
          VC_INIT(1,nvc)=VC(nvc) !store the previous volume
          VC_sum=0.d0
          
C         For each face in cell nvc:
          DO nfvl=1,NFVC(1,0,nvc) !for each face in cell
            nf=NFVC(2,nfvl,nvc) !global face #
            cnonode=NFVC(1,nfvl,nvc) !local node (nr) at centre of adjacent cell
            cnp=NPNODE(cnonode,nr)
            VERT=0
            
C           For each vertex in face nf:
            DO npl=1,NENFVC(0,nf) !for each vertex
              ne=NENFVC(npl,nf) !Delaunay element # for vertex
              VERT=VERT+1
              
c             calculate voronoi cell volumes = acinus volume
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                VERTICES(njj,VERT)=ZA(1,njj,1,ne)-XP(1,1,nj,np)
              ENDDO !njj
              
            ENDDO !npl
            VERTTOTAL=VERT
            
C           ..Make the end VERTEX the start VERTEX to complete the
C           loop
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              VERTICES(njj,VERT+1)=VERTICES(njj,1)
            ENDDO
            
            IF(NVCNODE(TYPE,cnonode).EQ.BOUNDARY)THEN
C             When the opposite node is a boundary node, the volume
C             cannot update correctly using the original approach.
C             i.e.
C             the method below (for non-boundary opposite nodes) uses
C             the vector from centre of cell to opposite node (at
C             centre
C             of adjacent cell) to calculate the unit normal to the
C             face. But for a boundary node, it is not within the host
C             volume, so it does not have valid Xi coordinates, and
C             therefore its global position will not change during
C             deformation. This section updates the boundary node
C             position s.t. it is always opposite the corresponding
C             internal boundary node (must be np).
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                XP(1,1,nj,cnp)=XP(1,1,nj,np)+0.1d0*XNFV(njj,nfvl)
              ENDDO !njj
            ELSE
C             ..Calculate the unit normal, distance
              XNFV(IDIST,nf)=0.d0
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                XNFV(njj,nf)=XP(1,1,nj,cnp)-XP(1,1,nj,np)
                XNFV(IDIST,nf)=XNFV(IDIST,nf)+XNFV(njj,nf)**2
              ENDDO
              XNFV(IDIST,nf)=1.d0/DSQRT(XNFV(IDIST,nf))
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                XNFV(njj,nf)=XNFV(njj,nf)*XNFV(IDIST,nf)
              ENDDO
            ENDIF !NVCNODE
            
C           ..Now have vertices and unit normal, can calculate the
C           face
C           area, and the volume contribution
            FACEAREA=0.d0
            DOTPROD=DDOT(NJ_LOC(NJL_GEOM,0,nr),VERTICES(1,1),1,XNFV(1,
     '        nf),1)
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              CENTROID(njj)=0.d0
              PROJECTION(njj)=DOTPROD*XNFV(njj,nf)
            ENDDO
            DO VERT=1,VERTTOTAL
              CALL CROSS(VERTICES(1,VERT),VERTICES(1,VERT+1),XPROD)
              DOTPROD=DDOT(NJ_LOC(NJL_GEOM,0,nr),XPROD,1,XNFV(1,nf),1)
              FACEAREA=FACEAREA+DOTPROD
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                CENTROID(njj)=CENTROID(njj)+DOTPROD*1.d0/3.d0
     '            *(VERTICES(njj,VERT)+VERTICES(njj,VERT+1) 
     '            +PROJECTION(njj))
              ENDDO
            ENDDO
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              XNFV(CENTINDX(njj),nf)=CENTROID(njj)
            ENDDO
            
            IF(DABS(FACEAREA).LT.LOOSE_TOL) THEN
c              ADD=.FALSE.
            ELSE
c              ADD=.TRUE.
              XNFV(FAREA,nf)=0.5d0*DABS(FACEAREA)
            ENDIF
            IF(FACEAREA.LT.0.d0) THEN
              VC_sum=VC_sum-1.d0/6.d0*DDOT(NJ_LOC(NJL_GEOM,0,nr),
     '          CENTROID,1,XNFV(1,nf),1)
            ELSE
              VC_sum=VC_sum+1.d0/6.d0*DDOT(NJ_LOC(NJL_GEOM,0,nr),
     '          CENTROID,1,XNFV(1,nf),1)
            ENDIF
            
          ENDDO !nf
          
          VC(nvc)=DABS(VC_sum)
C         Update the field that stores the volume information
          XP(1,1,nj_vol,np)=VC(nvc)
          VC(0)=VC(0)+VC(nvc)
          VC_field=VC_field+XP(1,1,nj_vol,np)
        ENDDO !nvc
        WRITE(OP_STRING,'('' Volume: '',F12.4)') VC_field
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
        
      CALL EXITS('UPVORO')
      RETURN
 9999 CALL ERRORS('UPVORO',ERROR)
      CALL EXITS('UPVORO')
      RETURN 1
      END


