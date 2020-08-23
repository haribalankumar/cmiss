      SUBROUTINE IPACTV(NPNODE,nr,nx,NYNP,XP,YP,YQ,ERROR,*)

C#### Subroutine: IPACTV
C###  Description:
C###    IPACTV inputs analytic activation time parameters for
C###    region nr and problem nx.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM),nr,nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,nc,nj,nk,nonode,NOQUES,np,nq,n_site,nv,ny
      REAL*8 RADIUS,RADIUS_MAX
      CHARACTER CHAR5*1
      LOGICAL FILEIP

      CALL ENTERS('IPACTV',*9999)

      nv=1 ! temporary

      NOQUES=0
      FILEIP=.FALSE.
      ICHAR=999

      FORMAT='(/$,'' Enter the number of initial actn sites [1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=N_SITES
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) N_SITES=IDATA(1)
      DO n_site=1,N_SITES
        WRITE(CHAR5,'(I1)') n_site
        IF(NJT.EQ.2) THEN
          FORMAT='($,'' Enter the coordinates of activation site '
     '      //CHAR5//' [0,0,0]:'',2D12.4)'
        ELSE
          FORMAT='($,'' Enter the coordinates of activation site '
     '      //CHAR5//' [0,0,0]:'',3D12.4)'
        ENDIF
        RDEFLT(1)=0.0d0
        RDEFLT(2)=0.0d0
        RDEFLT(3)=0.0d0
        IF(IOTYPE.EQ.3) THEN
          DO nj=1,NJT
            RDATA(1)=ACTVN_SITE(nj,n_site)
          ENDDO
        ENDIF
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO nj=1,NJT
            ACTVN_SITE(nj,n_site)=RDATA(nj)
          ENDDO
        ENDIF
      ENDDO !n_sites

      FORMAT='($,'' Enter the activation time interval '
     '  //'[100.0]: '',D12.4)'
      RDEFLT(1)=100.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=ACTVN_TIME
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RDELTA,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) ACTVN_TIME=RDATA(1)

      FORMAT='('' Specify whether actvn times calcd at [1]:'''//
     '  '/''   (1) Nodes'''//
     '  '/''   (2) Grid Points'''//
     '  '/''   (3) Unused'''//
     '  '/$,''   '',I1)'
      IF(IOTYPE.EQ.3)IDATA(1)=ANAL_CHOICE(nr)
      IDEFLT(1)=1
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) ANAL_CHOICE(nr)=IDATA(1)
      IF(ANAL_CHOICE(nr).EQ.2) THEN
        CALL ASSERT(CALL_GRID,'>>Define Grid first',ERROR,*9999)
      ENDIF

C     Calculate an activation sequence based on distance from the
C     sites of activation.  If a hermite description of the activation
C     sequence is required, then one should calculate activation times
C     at all grid points, and then fit these.
      IF(ANAL_CHOICE(nr).EQ.1) THEN !Calculate at nodes
        nc=1 ! variables first
        nk=1 ! value only
        RADIUS_MAX=-RMAX
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nr) !global var number
          YP(ny,7)=RMAX
          DO n_site=1,N_SITES
            IF(NJT.EQ.2) THEN
              RADIUS=DSQRT((XP(1,1,1,np)-ACTVN_SITE(1,n_site))**2+
     '          (XP(1,1,2,np)-ACTVN_SITE(2,n_site))**2)
            ELSE
              RADIUS=DSQRT((XP(1,1,1,np)-ACTVN_SITE(1,n_site))**2+
     '          (XP(1,1,2,np)-ACTVN_SITE(2,n_site))**2+
     '          (XP(1,1,3,np)-ACTVN_SITE(3,n_site))**2)
            ENDIF
            IF(RADIUS.LT.YP(ny,7)) YP(ny,7)=RADIUS
            IF(RADIUS.GT.RADIUS_MAX)RADIUS_MAX=RADIUS
          ENDDO !n_site
        ENDDO !np
C ***   Adjust the times
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nr) !global var number
          YP(ny,7)=YP(ny,7)*ACTVN_TIME/RADIUS_MAX
        ENDDO !nonode
      ELSEIF(ANAL_CHOICE(nr).EQ.2) THEN !Calculate at grid points
        DO nq=1,NQM
          YQ(nq,1,1)=RMAX
          DO n_site=1,N_SITES
            IF(NJT.EQ.2) THEN
              RADIUS=DSQRT((XP(1,nv,1,np)-ACTVN_SITE(1,n_site))**2+
     '          (XP(2,nv,1,np)-ACTVN_SITE(2,n_site))**2)
            ELSE
              RADIUS=DSQRT((XP(1,nv,1,np)-ACTVN_SITE(1,n_site))**2+
     '          (XP(2,nv,1,np)-ACTVN_SITE(2,n_site))**2+
     '          (XP(3,nv,1,np)-ACTVN_SITE(3,n_site))**2)
            ENDIF
            IF(RADIUS.LT.YQ(nq,1,1)) YQ(nq,1,1)=RADIUS
            IF(RADIUS.GT.RADIUS_MAX)RADIUS_MAX=RADIUS
          ENDDO !n_site
        ENDDO !nq
C ***   Adjust the times
        DO nq=1,NQM
          YQ(nq,1,1)=YQ(nq,1,1)*ACTVN_TIME/RADIUS_MAX
        ENDDO !nonode
      ELSE
        ERROR='>> Not implemented yet'
      ENDIF

      CALL EXITS('IPACTV')
      RETURN
 9999 CALL ERRORS('IPACTV',ERROR)
      CALL EXITS('IPACTV')
      RETURN 1
      END
CC AJPe


