      SUBROUTINE SGRESI(INDEX,ISEG,ISRESI,iw,
     '  NPNY,nr,nx,NYNO,RESID,ZP,CSEG,ERROR,*)

C#### Subroutine: SGRESI
C###  Description:
C###    SGRESI creates increment segment ISRESI.
C**** X1(nj) are the curvilinear coords of the arrow head
C**** Z1(nj) are the rect. cart. coords of the arrow head
C**** X2(nj) are the curvilinear coords of the arrow tail
C**** Z2(nj) are the rect. cart. coords of the arrow tail

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'reac00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISRESI,iw,
     '  NPNY(0:6,NYM,0:NRCM),nr,nx,ny1,ny2,NYNO(0:NYOM,NOOPM)
      REAL*8 RESID(*),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX_OLD,nc,nh1,nh2,nj,no,np,nv
      REAL*8 DZ,
     '  LENGTH,PTS(3,2),sin_THETA,cos_THETA,
     '  X1(3),X2(3),Z1(3),Z2(3),Z3(3),Z4(3),Z5(3)

      CALL ENTERS('SGRESI',*9999)
      CALL OPEN_SEGMENT(ISRESI,ISEG,iw,'RESI',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      nc=1  !temporary
c      nrc=2 !temporary AJP 1-12-94
      nv=1  !Temporary

C PJH 6Mar96 Following commented out until clarify application
c     DO nonode=1,NPNODE(0,nr)
c       np=NPNODE(nonode,nr)
c       DO nj=1,NJT
c         X1(nj)=ZP(1,nv,nj,np,nc)
c       ENDDO
!     Find the first no# corr. to np
c       NODE_RESIDUAL=.FALSE.
c       no=0
c       DO WHILE(.not.node_residual.and.no.LT.NOT(nrc,1,nr,nx))
c         no=no+1
c         IF(NPNY(4,NYNO(1,no),0).eq.np) NODE_RESIDUAL=.TRUE.
c       ENDDO
c       NODE_RESIDUAL=.TRUE.
c       IF(NODE_RESIDUAL) THEN
c         nores=no
c         DO nhx=1,NHP(np)
c           nh=NH_LOC(nhx,nx)
c           DO no=nores,NOT(nrc,1,nr,nx) !to find the first no# corr.
c                                        !to np & nh
c             IF(NPNY(4,NYNO(1,no),0).eq.np.AND.
c    '           NPNY(3,NYNO(1,no),0).eq.nh) THEN
c               nores=no
c               GOTO 20
c             ENDIF
c           ENDDO !no
c           WRITE(OP_STRING,'('' Warning!!Cannot find no in SGRESI'')')
c           CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c20         X2(nh)=ZP(1,nv,nh,np,nc)+FACTOR*RESID(nores)

      if(dop) then
        write(op_string,'('' NOT for nrc=1 is '',I4)') NOT(1,1,nr,nx)
        call writes(iodi,op_string,error,*9999)
      endif

      DO no=1,NOT(1,1,nr,nx),2 !loop over optimisation variables
        ny1=NYNO(1,no)
        ny2=NYNO(1,no+1)
        nh1=NPNY(3,ny1,1)
        nh2=NPNY(3,ny2,1)
        np =NPNY(4,ny1,1)
        if(dop) then
          write(op_string,'('' no='',I5,'' np='',I5,'
     '      //''' nh1='',I1,'' nh2='',I1)') no,np,nh1,nh2
          call writes(iodi,op_string,error,*9999)
        endif
        X1(1)=ZP(1,nv,nh1,np,nc)
        X1(2)=ZP(1,nv,nh2,np,nc)
        X2(1)=X1(1)+FACTOR*RESID(no)
        X2(2)=X1(2)+FACTOR*RESID(no+1)
          CALL XZ(ITYP10(nr),X1,Z1)
          CALL XZ(ITYP10(nr),X2,Z2)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Z1(nj)='',3E12.4)') (Z1(nj),nj=1,NJT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Z2(nj)='',3E12.4)') (Z2(nj),nj=1,NJT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C ***     Draw residual vector
          PTS(1,1)=Z1(1)
          PTS(1,2)=Z2(1)
          PTS(2,1)=Z1(2)
          PTS(2,2)=Z2(2)
          PTS(3,1)=Z1(3)
          PTS(3,2)=Z2(3)
          CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)

C ***     Draw arrow head
!         Arrow length is scaled by reaction YP(ny,iy)
!         and user defined factor.
          LENGTH=DSQRT((Z2(1)-Z1(1))**2+(Z2(2)-Z1(2))**2)
          IF(LENGTH.GT.1.D-12) THEN
            DZ=LENGTH/4.0D0 !is arrow head size
!           Z3(j) are coords of pt on vector at start of head
            Z3(1)=Z1(1)+(LENGTH-DZ*0.866d0)/LENGTH*(Z2(1)-Z1(1))
            Z3(2)=Z1(2)+(LENGTH-DZ*0.866d0)/LENGTH*(Z2(2)-Z1(2))
            sin_THETA=(Z2(2)-Z1(2))/LENGTH
            cos_THETA=(Z2(1)-Z1(1))/LENGTH
!           Z4(j) & Z5(j) are coords at end of arrow head
            Z4(1)=Z3(1)+0.5d0*DZ*sin_THETA
            Z4(2)=Z3(2)-0.5d0*DZ*cos_THETA
            Z5(1)=Z3(1)-0.5d0*DZ*sin_THETA
            Z5(2)=Z3(2)+0.5d0*DZ*cos_THETA
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Z3(nj)='',3E12.4)')
     '          (Z3(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Z4(nj)='',3E12.4)')
     '          (Z4(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Z5(nj)='',3E12.4)')
     '          (Z5(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            PTS(1,1)=Z2(1)
            PTS(1,2)=Z4(1)
            PTS(2,1)=Z2(2)
            PTS(2,2)=Z4(2)
            PTS(3,1)=Z2(3)
            PTS(3,2)=Z4(3)
            CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)
            PTS(1,1)=Z2(1)
            PTS(1,2)=Z5(1)
            PTS(2,1)=Z2(2)
            PTS(2,2)=Z5(2)
            PTS(3,1)=Z2(3)
            PTS(3,2)=Z5(3)
            CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)

          ELSE
!           WRITE(OP_STRING,'('' Arrow length is zero'')')
!           CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
      ENDDO !no
c       ENDIF !node_residual
c     ENDDO !nonode

      CALL CLOSE_SEGMENT(ISRESI,iw,ERROR,*9999)

      CALL EXITS('SGRESI')
      RETURN
 9999 CALL ERRORS('SGRESI',ERROR)
      CALL EXITS('SGRESI')
      RETURN 1
      END


