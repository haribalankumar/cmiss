      SUBROUTINE SGSTRA(INDEX,IBT,IDO,INP,IPOINTTYP,ISEG,ISSTRA,iw,
     '  LD,NAN,NBH,NBJ,ne,NHE,NKHE,NKJE,NPF,NPNE,nr,
     '  NVHE,NVJE,NW,nx,CURVCORRECT,
     '  PG,RGX,SE,STATIC,
     '  XA,XE,XG,XID,XIG,XP,ZA,ZE,ZG,ZP,CSEG,ERROR,*)

C#### Subroutine: SGSTRA
C###  Description:
C###    SGSTRA creates element strain segments ISSTRA(ne).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  IPOINTTYP,INP(NNM,NIM,NBFM),ISEG(*),ISSTRA(NEM,NGRSEGM),
     '  iw,LD(NDM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM),NBJ(NJM),ne,NHE,
     '  NKHE(NKM,NNM,NHM),NKJE(NKM,NNM,NJM),NPF(9),NPNE(NNM,NBFM),nr,
     '  NVHE(NNM,NBFM,NHM),NVJE(NNM,NBFM,NJM),NW,nx
      REAL*8 CURVCORRECT(2,2,NNM),PG(NSM,NUM,NGM,NBM),RGX(NGM),
     '  SE(NSM,NBFM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL STATIC
!     Local Variables
      INTEGER INDEX_OLD,na,NBFF,NBG,nc,nj,NSF,NSG,nk,nn
      REAL*8 PF1,TEMP,TEMP2
      CHARACTER CLABEL*52

      CALL ENTERS('SGSTRA',*9999)
C GBS 10-NOV-1994
      nc=1   !Temporary
      CLABEL='STRA'
      IF(STATIC) THEN
        CLABEL(43:52)='          '
      ELSE
        WRITE(CLABEL(43:52),'(E10.3)') TIME
      ENDIF
      CALL OPEN_SEGMENT(ISSTRA(ne,NTSTRA),ISEG,iw,'STRA',INDEX,
     '  INDEX_OLD,ne,1,CSEG,ERROR,*9999)

      IF(STATIC) THEN
        CALL XPXE(NBJ,NKJE,NPF,NPNE,nr,NVJE,
     '    SE,XA(1,1,ne),XE,XP,ERROR,*9999)
        CALL ZPZE(NBH,nc,NHE,NKHE,NPF,NPNE,nr,NVHE,NW,nx,
     '    CURVCORRECT,SE,ZA,ZE,ZP,ERROR,*9999)
! PJH 12Nov94 Commented this out because CG not used?
!        CALL CPCG(NW,NBH(NH_LOC(1,nx),nc),NPNE,nr,nx,CE,CG,CP,PG,
!                  ERROR,*9999)
        CALL STRAIN2(INDEX,IBT,IDO,INP,IPOINTTYP,iw,LD,NAN,NBH,NBJ,
     '    ne,NHE,nr,nx,PG,RGX,XE,XG,XID,XIG,ZE,ZG,ERROR,*9999)

      ELSE IF(.NOT.STATIC) THEN   !time dependent field stored in YP/ZP
C       following lines added 14 Dec 1990 AAY
C       Because strain2 uses the PG array it is easiest to interpolate time
C       before we call it. Then it can call zeex50 in the normal way.
C       The ZE array is interpolated at TIME and stored in ZE for
C       ns=1..NST(NBJ(nj)). Note this means NBH(nj,nc) (the Fourier basis)
C       must have the same spatial interpolation as NBJ(nj). If ZE is a
C       displacement from XE it gets added to XE and the result is stored in
C       ZE. This is passed to zeex50 as the deformed geometry. Note this
C       implies that the displacement field ze and the geometric field xe are
C       defined in the same coord system. Similarly the ZP array is
C       interpolated at time TIME_REF, added to XE and stored in XE. Zeex50
C       will take this as the undeformed geometry. Strain2 is called with XE
C       and ZE as previously EXCEPT NBJ IS SUBSTITUTED FOR NBH. This is
C       because NBH is a Fourier basis and zeex50 wants the geometric basis.
C
        CALL XPXE(NBJ,NKJE,NPF,NPNE,nr,NVJE,SE,XA(1,1,ne),XE,XP,ERROR,
     '    *9999)
        CALL ZPZE(NBH,nc,NHE,NKHE,NPF,NPNE,nr,NVHE,NW,nx,CURVCORRECT,SE,
     '    ZA,ZE,ZP,ERROR,*9999)
        DO nj=1,NJT
          NBFF=NBH(nj,nc)  !time basis
          NBG=NBJ(nj)   !corresponding spatial basis
          NSG=0
          DO nn=1,NNT(NBG)
            DO nk=1,NKT(0,NBG)
              NSG=NSG+1
              TEMP=0.0D0
              TEMP2=0.0D0
              DO na=1,IBT(2,NIT(NBFF),NBFF) !#terms in series
                NSF=(nn-1)*NKT(0,NBFF)+(na-1)*NKT(0,NBG)+nk
                TEMP=TEMP+PF1(na,1,TIME)*ZE(NSF,nj)
                TEMP2=TEMP2+PF1(na,1,TIME_REF)*ZE(NSF,nj)
              ENDDO
              ZE(NSG,nj)=TEMP+XE(NSG,nj)
              XE(NSG,nj)=TEMP2+XE(NSG,nj)
            ENDDO
          ENDDO
        ENDDO
        CALL STRAIN2(INDEX,IBT,IDO,INP,IPOINTTYP,
     '    iw,LD,NAN,NBJ,NBJ,ne,NHE,nr,nx,
     '    PG,RGX,XE,XG,XID,XIG,ZE,ZG,ERROR,*9999)

      ENDIF

      CALL CLOSE_SEGMENT(ISSTRA(ne,NTSTRA),iw,ERROR,*9999)

      CALL EXITS('SGSTRA')
      RETURN
 9999 CALL ERRORS('SGSTRA',ERROR)
      CALL EXITS('SGSTRA')
      RETURN 1
      END


