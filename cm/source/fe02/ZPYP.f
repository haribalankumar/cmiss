      SUBROUTINE ZPYP(iy,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '  nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*)

C#### Subroutine: ZPYP
C###  Description:
C###    ZPYP transfers global node parameters ZP(nk,nv,nh,np,nc) and
C###    auxiliary element parameters ZA(na,nh,nc,ne) to global vector
C###    YP(ny,iy).

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
      INTEGER na,nc,ne,nk,nh,nhx,noelem,nonode,np,nrc,nv,ny
      CHARACTER CHAR*1
!     EXTERNAL FUNCTION
      INTEGER IDIGITS

      CALL ENTERS('ZPYP',*9999)

C***  Subroutine as of 25/11/94 ajp,cpb.

      !KSB 04/05/05 - Adding error call to increase NIYM
      IF(NIYM.LT.iy) THEN
        WRITE(CHAR,'(I1)') IDIGITS(iy)
        WRITE(ERROR,'(''Increase NIYM to '',I'//CHAR//')') iy
        GOTO 9999
      ENDIF
      
      nrc=0 !only want the global variables
C cpb 28/3/00 Might call this for problems that don't have nc=2 defined.
C      DO nc=1,2
      DO nc=1,NCT(nr,nx)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C MPN 3Mar96 DO nhx=1,NH_LOC(0,nx)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,nc)
              DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                YP(ny,iy)=ZP(nk,nv,nh,np,nc)
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
              YP(ny,iy)=ZA(na,nh,nc,ne)
            ENDDO !na
          ENDDO !nh
        ENDDO !noelem (ne)
      ENDDO !nc

      CALL EXITS('ZPYP')
      RETURN
 9999 CALL ERRORS('ZPYP',ERROR)
      CALL EXITS('ZPYP')
      RETURN 1
      END


