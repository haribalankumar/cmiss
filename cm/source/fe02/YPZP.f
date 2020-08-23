      SUBROUTINE YPZP(iy,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,
     '  NYNE,NYNP,YP,ZA,ZP,ERROR,*)

C#### Subroutine: YPZP
C###  Description:
C###    YPZP transfers global vector YP(ny,iy) to global node parameters
C###    ZP(nk,nv,nh,np,nc) and auxiliary element parameters
C###    ZA(na,nh,nc,ne).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER iy,NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NHP(NPM),NKH(NHM,NPM,NCM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVHP(NHM,NPM,NCM),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER na,nc,ne,nh,nhx,nk,noelem,nonode,np,nrc,nv,ny
      CHARACTER CHAR*1
!     EXTERNAL FUNCTION
      INTEGER IDIGITS

      CALL ENTERS('YPZP',*9999)

C***  Subroutine as of 1/12/94 ajp,cpb.

      !KSB 04/05/05 - Adding error call to increase NIYM
      IF(NIYM.LT.iy) THEN
        WRITE(CHAR,'(I1)') IDIGITS(iy)
        WRITE(ERROR,'(''Increase NIYM to '',I'//CHAR//')') iy
        GOTO 9999
      ENDIF

      nrc=0 !only want the global variables
      DO nc=1,NCT(nr,nx)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C MPN 3Mar96 DO nhx=1,NH_LOC(0,nx)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,nc)
              DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
              ZP(nk,nv,nh,np,nc)=YP(ny,iy)
              ENDDO !nk
            ENDDO !nv
          ENDDO !nh
        ENDDO !nonode (np)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            DO na=1,NAT(NBH(nh,nc,ne))
              ny=NYNE(na,nh,nrc,nc,ne)
              ZA(na,nh,nc,ne)=YP(ny,iy)
            ENDDO !na
          ENDDO !nh
        ENDDO !noelem (ne)
      ENDDO !nc

      CALL EXITS('YPZP')
      RETURN
 9999 CALL ERRORS('YPZP',ERROR)
      CALL EXITS('YPZP')
      RETURN 1
      END


