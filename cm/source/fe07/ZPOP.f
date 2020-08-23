      SUBROUTINE ZPOP(ip,NBH,nc,NEELEM,NHP,NKH,NPNODE,nr,
     '  NVHP,nx,NYNE,NYNP,RP,ZA,ZP,FIX,ERROR,*)

C#### Subroutine: ZPOP
C###  Description:
C###    ZPOP prints solutions ZP(nk,nv,nh,np,nc) and ZA(na,nh,nc,ne)
C###    and, if IWRIT3(nr,nx)=2, the residuals RP.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ip,NBH(NHM,NCM,NEM),nc,NEELEM(0:NE_R_M,0:NRM),
     '  NHP(NPM),NKH(NHM,NPM,NCM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVHP(NHM,NPM,NCM),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 RP(NYM),ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM)
!     Local Variables
      INTEGER na,nb,ne,nh,nhx,nk,NKMAX,noelem,nonode,np,nv
      CHARACTER CHAR1*1,CHAR2*2,TITLE(4)*24
      DATA TITLE/'Initial equil. solution ','Incremented soln vector ',
     '           'Intermediate soln vector','Equilibrium soln vector '/

      CALL ENTERS('ZPOP',*9999)

      WRITE(OP_STRING,'(/1X,A,'' for region '',I2)') TITLE(ip),nr
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NKMAX=1
        DO nhx=1,NHP(np)
           nh=NH_LOC(nhx,nx)
          IF(NKH(nh,np,nc).GT.NKMAX) NKMAX=NKH(nh,np,nc)
        ENDDO
        DO nhx=1,NHP(np)
          nh=NH_LOC(nhx,nx)
          DO nv=1,NVHP(nh,np,nc)
            WRITE(CHAR1,'(I1)') NKH(nh,np,nc)
            WRITE(CHAR2,'(I2)') 12*(NKMAX-NKH(nh,np,nc))+1
            IF(IWRIT3(nr,nx).EQ.1) THEN
              FORMAT='('' ZP(nk,nv='',I2,'',nh='',I1,'',np='',I5,'
     '          //''',nc='',I1,''): '','//CHAR1(1:1)//'D12.4)'
              WRITE(OP_STRING,FORMAT) nv,nh,np,nc,
     '          (ZP(nk,nv,nh,np,nc),nk=1,NKH(nh,np,nc))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(IWRIT3(nr,nx).EQ.2) THEN
              FORMAT='('' ZP(nk,nv='',I2,'',nh='',I1,'',np='',I5,'
     '          //''',nc='',I1,''): '','//CHAR1(1:1)
     '          //'D12.4,'//CHAR2(1:2)
     '          //'1X,''  RP: '','//CHAR1(1:1)//'D12.4,'' ('''
     '          //CHAR1(1:1)//'L1,'' )'')'
              WRITE(OP_STRING,FORMAT) nv,nh,np,nc,
     '          (ZP(nk,nv,nh,np,nc),nk=1,NKH(nh,np,nc)),
     '          (RP(NYNP(nk,nv,nh,np,0,nc,nr)),nk=1,NKH(nh,np,nc)),
     '          (FIX(NYNP(nk,nv,nh,np,0,nc,nr)),nk=1,NKH(nh,np,nc))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nv
        ENDDO !nh
      ENDDO !nonode (np)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh,nc,ne)
          IF(NAT(nb).GT.0) THEN
            FORMAT='(/'' ZA(na,nh='',I1,'',nc='',I1,'',ne='',I5,'
     '        //'''): '',3(D12.4,5X),/(28X,3(D12.4,5X)))'
            WRITE(OP_STRING,FORMAT)
     '        nh,nc,ne,(ZA(na,nh,nc,ne),na=1,NAT(nb))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            IF(IWRIT3(nr,nx).EQ.2) THEN
              FORMAT='(24X,''RP: '',3(D12.4,'' ('',L1,'') ''),'
     '          //'/(28X,3(D12.4,'' ('',L1,'') '')))'
              WRITE(OP_STRING,FORMAT) (RP(NYNE(na,nh,2,nc,ne)),
     '          FIX(NYNE(na,nh,2,nc,ne)),na=1,NAT(nb))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('ZPOP')
      RETURN
 9999 CALL ERRORS('ZPOP',ERROR)
      CALL EXITS('ZPOP')
      RETURN 1
      END


