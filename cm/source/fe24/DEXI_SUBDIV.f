      SUBROUTINE DEXI_SUBDIV(
     '  DEFORM,IBT,IDO,INP,LD,NBH,NBJ,ND0,ND1,
     '  NELIST,NHE,NKHE,NKJE,
     '  NOXIPT,NPF,NPNE,NRE,NVHE,NVJE,NW,nx,CURVCORRECT,SE,SQ,
     '  SQMAX,XA,XE,XID,XP,ZA,ZD,ZP,
     '  ERROR,*)

C#### Subroutine: DEXI_SUBDIV
C###  Description:
C###    DEXI_SUBDIV uses subdivision to calculate approximate
C###    xi positions.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),SQMAX,SQ(NDM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),LD(NDM),ND0,ND1,NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NELIST(0:NEM),NHE(NEM,NXM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),
     '  NOXIPT(3),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3,NXM),nx
      LOGICAL DEFORM

!     Local Variables
      INTEGER nb,nd,ne,nj,nk,nn,np,nolist,ns,i1,i2,i3
      REAL*8 dist,DXI1,DXI2,DXI3,XI(3),PSI,
     '  SUB_PSI(NSM,0:20,0:20,0:20,NJT),XD(3),Z1(3)

      CALL ENTERS('DEXI_SUBDIV',*9999)
!     Output arguments
      IF(DOP) THEN
        WRITE(OP_STRING,'('' NOXIPT='',I5,'','',I5,'','',I5)') 
     '    NOXIPT(1),NOXIPT(2),NOXIPT(3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' SQMAX='',E11.2)') SQMAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ND=['',I5,'','',I5,'']'')') ND0,ND1
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO ne=1,NELIST(0)
          WRITE(OP_STRING,'('' NELIST['',I5,'']='',I5)') ne,NELIST(ne)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF
!     A couple of checks
      DO nj=1,3
        CALL ASSERT(NOXIPT(nj).LE.20, ' >> sub values should be <= 20',
     '    ERROR,*9999)
      ENDDO
      DO ne=1,NELIST(0)
         DO nj=1,3
           CALL ASSERT(NBJ(nj,NELIST(ne)).EQ.NBJ(nj,NELIST(1)),
     '       ' >> elements need to have same basis function',
     '       ERROR,*9999)
         ENDDO
      ENDDO

! Don't initialise to allow recalling with different elements
! this is down in DEXI using NEW/OLD flags.
!      DO nd=ND0,ND1
!        SQ(nd)=SQMAX
!        LD(nd)=0
!      ENDDO
      DXI1=1.0d0 / DBLE(NOXIPT(1)+1)
      DXI2=1.0d0 / DBLE(NOXIPT(2)+1)
      DXI3=1.0d0 / DBLE(NOXIPT(3)+1)
! creates points which are evenly spread through 0.5|1.0|1.0|1.0|0.5
      DO i3=0,NOXIPT(3)
        XI(3)=(DBLE(i3)+0.5d0)*DXI3
        DO i2=0,NOXIPT(2)
          XI(2)=(DBLE(i2)+0.5d0)*DXI2
          DO i1=0,NOXIPT(1)
            XI(1)=(DBLE(i1)+0.5d0)*DXI1
            DO nj=1,NJT
              nb=NBJ(nj,NELIST(1))
              ns=0
              DO nn=1,NNT(nb)
                DO nk=1,NKT(nn,nb)
                  ns=ns+1
                  SUB_PSI(ns,i1,i2,i3,nj) =
     '              PSI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,nn,nk,1,XI)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO nolist=1,NELIST(0)
        ne=NELIST(nolist)
        IF(DEFORM)THEN
c           Note that deformed coords --> XE
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '      NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '      NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '      XE,ZP,ERROR,*9999)
        ELSE
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '      NPF(1,1),NPNE(1,1,ne),
     '      NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA,XE,XP,ERROR,*9999)
        ENDIF
        np=0
        DO i3=0,NOXIPT(3)
          XI(3)=(DBLE(i3)+0.5d0)*DXI3
          DO i2=0,NOXIPT(2)
            XI(2)=(DBLE(i2)+0.5d0)*DXI2
            DO i1=0,NOXIPT(1)
              XI(1)=(DBLE(i1)+0.5d0)*DXI1
              np=np+1
              DO nj=1,NJT
                nb=NBJ(nj,NELIST(1))
                XD(nj)=0.0d0
                DO ns=1,NST(nb)
                  XD(nj)=XD(nj)+SUB_PSI(ns,i1,i2,i3,nj)*XE(ns,nj)
                ENDDO
              ENDDO
              CALL XZ(JTYP3,XD,Z1)

              DO nd=ND0,ND1
                 dist=0.0d0
                 DO nj=1,NJT
                   dist=dist+(ZD(nj,nd)-Z1(nj))**2
                 ENDDO
                 IF(LD(nd).EQ.0.OR.dist.LT.SQ(nd)) THEN
                   LD(nd)=ne
                   XID(1,nd)=XI(1)
                   XID(2,nd)=XI(2)
                   XID(3,nd)=XI(3)
                   SQ(nd)=dist
                 ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDDO
      CALL EXITS('DEXI_SUBDIV')
      RETURN
 9999 CALL ERRORS('DEXI_SUBDIV',ERROR)
      CALL EXITS('DEXI_SUBDIV')
      RETURN 1
      END


