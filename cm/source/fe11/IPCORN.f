      SUBROUTINE IPCORN(INP,NBJ,NKJE,NPF,NPNE,
     '  nr,NVHP,NVJE,NW,PG,SE,XA,XE,XP,ERROR,*)

C#### Subroutine: IPCORN
C###  Description:
C###    IPCORN defines corner nodes for boundary element problems in
C###    region nr.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER INP(NNM,NIM,NBFM),NBJ(NJM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  nr,NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3)
      REAL*8 PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,ICHAR,INFO,loop,LOOP1,LOOP2,LOOPT,nb,nc,NDUMMY,ne,
     '  ng,NG_CLOSE(10),nh,nj,nn,noelem,NOQUES,np,ns,nu,nv,nx
      REAL*8 A,DIST,DIST1,SUMXG,TN(3,2,10),XG_LOCAL(3,4),XN_LOCAL(3,10)
      LOGICAL CONTINUE,FILEIP,FOUND,INTERFACE

      CALL ENTERS('IPCORN',*9999)

      nc=1 !temporary cpb 23/11/94
      nv=1 !temporary cpb 23/11/94
      nx=1 !temporary

      FILEIP=.FALSE.
      NOQUES=0
      FOUND=.FALSE.
      ICHAR=999

      np=1
      CONTINUE=.TRUE.
      IDEFLT(1)=-1
      DO WHILE(CONTINUE)
        FORMAT='(/$,'' Enter a corner node number [Exit]: '',I2)'
        IF(IOTYPE.EQ.3) THEN
           DO WHILE((NVHP(1,np,nc).LE.2).AND.(np.LE.NPM))
             np=np+1
           ENDDO
           NDUMMY=np
           IF(np.GT.NPM)NDUMMY=IDEFLT(1)
           IDATA(1)=NDUMMY
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPT(1),
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IDATA(1).GT.IDEFLT(1)) THEN
          IF(IOTYPE.NE.3) THEN
C **        Find the number of elements that share the corner node NP
            nh=NH_LOC(1,nx) !Assumes 1 dependent variable
            np=IDATA(1)
            NVHP(nh,np,nc)=2 !Initial value
            DO ne=1,NET(nr)
              nb=NBJ(1,ne)
              DO nn=1,NNT(nb)
                IF(IDATA(1).EQ.NPNE(nn,nb,ne)) THEN
                  FOUND=.TRUE.
C cpb 23/11/94 Huh ???
C                  NCNP(NCNP(nh,0,np,nr),nh,np,nr)=NE
C                  NCNP(0,nh,np,nr)=NCNP(0,nh,np,nr)+1
                  NVHP(nh,np,nc)=ne
                  NVHP(nh,np,nc)=NVHP(nh,np,nc)+1
                ENDIF
              ENDDO
            ENDDO
            IF(FOUND) NVHP(nh,np,nc)=NVHP(nh,np,nc)-1
            FOUND=.FALSE.
            CALL ASSERT(NVHP(nh,np,nc).GT.2,
     '        '>>This is not a corner or edge node',ERROR,*9999)
C** Check whether any of the elements sharing the corner node lie on the
C** same face.  Only retain one element per face.  Only occurs when we
C** have an edge in 3d.
            IF((NVHP(nh,np,nc).GT.4).AND.
     '        (NIT(NBJ(1,NVHP(nh,np,nc))).EQ.2))THEN
              !More than 3 shared elements and 2d integrals.
              !ie probably a 3d edge
              DO loop=1,NVHP(nh,np,nc)-1
                !Calculate tangents and normals for each element
                ne=NVHP(nh,np,nc)
                nb=NBJ(1,ne)
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C ***           Find Gauss point closest to the current node
                DIST1=RMAX
                DO ng=1,NGT(nb)
                  DIST=0.0D0
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    SUMXG=0.0D0
                    DO ns=1,NST(nb)
                      SUMXG=SUMXG+PG(ns,1,ng,nb)*XE(ns,nj)
                    ENDDO
                    DIST=DIST+(SUMXG-XP(1,nv,nj,np))**2
                  ENDDO
                  IF(DIST.LE.DIST1)THEN
                    DIST1=DIST
                    NG_CLOSE(loop)=NG
                  ENDIF
                ENDDO
C ***           Find normal at this point
                ng=NG_CLOSE(loop)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  DO nu=1,NUT(nb)
                    SUMXG=0.D0
                    DO ns=1,NST(nb)
                      SUMXG=SUMXG+PG(ns,nu,ng,nb)*XE(ns,nj)
                    ENDDO
                    XG_LOCAL(nj,nu)=SUMXG
                  ENDDO
                ENDDO
                INTERFACE=.FALSE. !May need changing AJP 21-1-93
                CALL NORMAL(ne,nr,NW,XG_LOCAL,
     '            XN_LOCAL(1,ne),INTERFACE,
     '            ERROR,*9999)
                !Find unit out. normal.
C ***           Find two tangent vectors in the plane and adjust them
C ***           so that they lie in the element
                nn=1
                DO WHILE(np.NE.NPNE(nn,nb,ne))
                  nn=nn+1 !Find local node number of NP
                ENDDO
                IF(INP(nn,1,nb).GT.1)THEN !Want negative of the tangent
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,1,loop)=-XG_LOCAL(nj,2)
                  ENDDO
                ELSE
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,1,loop)=XG_LOCAL(nj,2)
                  ENDDO
                ENDIF
                IF(INP(nn,2,nb).GT.1)THEN !Want negative of the tangent
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,2,loop)=-XG_LOCAL(nj,4)
                  ENDDO
                ELSE
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,2,loop)=XG_LOCAL(nj,4)
                  ENDDO
                ENDIF
              ENDDO !End of loop
              LOOP1=1
              LOOP2=2
              LOOPT=NVHP(nh,np,nc)-1
 900          ne=NVHP(nh,np,nc)
 1000         CALL ANGLE(TN(1,1,LOOP1),TN(1,1,LOOP2),A,
     '                   XN_LOCAL(1,LOOP1),XN_LOCAL(1,LOOP2),
     '                   ERROR,*9999)
              IF(A.LT.0.1D0)THEN !assume elements on the same face
                DO noelem=LOOP2+1,NVHP(nh,np,nc)-1
                  NVHP(nh,np,nc)=
     '              NVHP(nh,np,nc)
                  !Remove LOOP2 element and move remaining
                  !elements down the list
                ENDDO
                NVHP(nh,np,nc)=NVHP(nh,np,nc)-1
                LOOP2=LOOP2+1
                IF(LOOP2.LE.LOOPT)GOTO 1000
              ENDIF
              LOOP1=LOOP1+1
              LOOP2=LOOP2+1
              IF(LOOP1.LT.LOOPT)GOTO 900
            ENDIF !End of removing edge elements on the same face
!OLD needs rewriting
C            !Set up NVJE for other corner elements
C            DO noelem=2,NVHP(nh,np,nc)
C              ne=NVHP(nh,np,nc)
C              nb=NBJ(1,ne)
C              DO nn=1,NNT(nb)
C                IF(np.EQ.NPNE(nn,nb,ne)) THEN
C                  DO nb2=1,NBFT
C                    IF(NNT(nb2).EQ.NNT(nb))THEN
C                       NVJE(nn,nb2,nh,ne)=noelem
C                     ENDIF
C                  ENDDO
C                 ENDIF
C               ENDDO
C            ENDDO
            WRITE(OP_STRING,'('' Node'',I4,'' shared by  elements '','
     '        //'4(I3,2X))') IDATA(1),(NVHP(nh,IDATA(1),nc),
     '        I=2,NVHP(nh,IDATA(1),nc))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(IOTYPE.EQ.3) THEN
             np=np+1
          ENDIF
        ELSE
          CONTINUE=.FALSE.
        ENDIF
      ENDDO
      CALL EXITS('IPCORN')
      RETURN
 9999 CALL ERRORS('IPCORN',ERROR)
      CALL EXITS('IPCORN')
      RETURN 1
      END


