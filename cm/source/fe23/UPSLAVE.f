      SUBROUTINE UPSLAVE(IBT,IDO,INP,NBJ,ne1,NEELEM,NKJE,
     '  np2,NPF,NPNE,NRE,nr1,nr2,NVJE,NVJP,SE,XIP,XP,ERROR,*)

C#### Subroutine: UPSLAVE
C###  Description:
C###    UPSLAVE updates the nodal values and their derivatives of the
C###    slave nodes after host-mesh fitting. Math related to this
C###    routine can be found at /product/cmiss/documents/pdfarchive/
C###    hostmesh_fit.pdf
C**** Written by Kumar Mithraratne, Nov. 2002.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),ne1,np2,NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  nr1,nr2,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),XIP(NIM,NPM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IPIV(3),INFO,nbf_host,nbf_slave,ne2,
     '  NENP(NPM,0:NEPM,0:NRM),nih,nis,NITB_H,NITB_S,nj,njj,nj_field,
     '  nj_max_H,nj_max_S,nk,nu,i,ndim,ndimm,ii,nv,nvt_S
      REAL*8 A(3,3),B(3),AO(3,3,3),AN(3,3,3),C(3),DEF(3,3),
     '  dxidn(3,3,NVM),d2xidn2(3,3,3),PXI,UNDEF(3,3),XA(NAM,NJM),
     '  XE(NSM,NJM)
      EXTERNAL DGETRS,DGETRF

      CALL ENTERS('UPSLAVE',*9999)

      nbf_host=NBJ(1,ne1)
      NITB_H=NIT(nbf_host)                 ! no. of xi coords in host elem.
      nj_max_H=NJ_LOC(NJL_GEOM,0,NRE(ne1)) ! no. of global coords in host elem.

      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)
      ne2=NENP(np2,1,nr2)
      nbf_slave=NBJ(1,ne2)
      NITB_S=NIT(nbf_slave)                ! no. of xi coords in slave elem.
      nj_max_S=NJ_LOC(NJL_GEOM,0,NRE(ne2)) ! no. of glob. coords in slave elem.
      nvt_S=NVJP(1,np2)  ! It is assumed that all nj s have same no. of versions

      CALL ASSERT(NITB_H.EQ.nj_max_H,'No. of local and global '//
     '  'coordinates of the host mesh should be the same',ERROR,*9999)
      CALL ASSERT(nj_max_S.EQ.nj_max_H,'No. of global coordinates '//
     '  'of the slave and host should be the same',ERROR,*9999)

C     New nodal values of slave node.
C     It is assumed that both host and slave meshes have the same
C     no. of global coordinates.

      CALL XPXE(NBJ(1,ne1),NKJE(1,1,1,ne1),NPF(1,1),NPNE(1,1,ne1),
     '  nr1,NVJE(1,1,1,ne1),SE(1,1,ne1),XA,XE,XP,ERROR,*9999)
      DO nv=1,nvt_S
C JWF 10/08/04 changed NJL_FIEL to NJL_GEOM
C The max no. of geometric variables should control loop
C and not the max. no. of field variables.  
        DO njj=1,NJ_LOC(NJL_GEOM,0,NRE(ne1))
C        DO njj=1,NJ_LOC(NJL_FIEL,0,NRE(ne1))
          nj_field=NJ_LOC(NJL_FIEL,njj,NRE(ne1))
          nj=NJ_LOC(NJL_GEOM,njj,NRE(ne1))
          XP(1,nv,nj,np2)=PXI(IBT(1,1,nbf_host),IDO(1,1,0,nbf_host),
     '      INP(1,1,nbf_host),nbf_host,1,XIP(1,np2),XE(1,nj_field))
          IF (ITYP10(nr1).EQ.4) THEN !prolate spheroidal coordinates
C ensures the theta coordinates are between zero and 2pi
            IF(nj.EQ.NJ_LOC(NJL_FIEL,3,nr2)) THEN !theta coordinate
              DO WHILE (XP(1,nv,nj,np2).GT.(2.0d0*PI))
                XP(1,nv,nj,np2)=XP(1,nv,nj,np2)-(2.0d0*PI)
              ENDDO
            ENDIF!theta coordinate
          ENDIF!prolate spheroidal coordinates
        ENDDO !njj
      ENDDO !nv

C     New first derivatives of slave node.

      DO njj=1,nj_max_H
        nj_field=NJ_LOC(NJL_FIEL,njj,NRE(ne1))
        nj=NJ_LOC(NJL_GEOM,njj,NRE(ne1))

        DO nih=1,NITB_H
          nu=INT((nih*nih+nih+2+0.1)/2.0)  ! 2,4,7
          UNDEF(njj,nih)=PXI(IBT(1,1,nbf_host),IDO(1,1,0,nbf_host),
     '      INP(1,1,nbf_host),nbf_host,nu,XIP(1,np2),XE(1,nj))
                                                          !Undeformed Jacobian
          DEF(njj,nih)=PXI(IBT(1,1,nbf_host),IDO(1,1,0,nbf_host),
     '      INP(1,1,nbf_host),nbf_host,nu,XIP(1,np2),XE(1,nj_field))
                                                          !Deformed Jacobian
        ENDDO !nih
      ENDDO !njj

      DO nis=1,NITB_S
        nk=INT((nis*nis-nis+4+0.1)/2.0)  ! 2,3,5
        DO nv=1,nvt_S
          IF(nj_max_H.EQ.3) THEN !volume host mesh
            DO nj=1,nj_max_H
              DO nih=1,NITB_H
                A(nj,nih)=UNDEF(nj,nih)
              ENDDO
              B(nj)=XP(nk,nv,nj,np2)
            ENDDO

            CALL DGETRF(NITB_H,nj_max_H,A,3,IPIV,INFO)
            CALL DGETRS('N',nj_max_H,1,A,3,IPIV,B,3,INFO)
            DO nih=1,NITB_H
              dxidn(nih,nis,nv)=B(nih)
            ENDDO
          ELSEIF(nj_max_H.EQ.2) THEN !area host mesh

          ENDIF

          DO nj=1,nj_max_S
            XP(nk,nv,nj,np2)=0.0d0
            DO nih=1,NITB_H
              XP(nk,nv,nj,np2)=XP(nk,nv,nj,np2)+
     '          dxidn(nih,nis,nv)*DEF(nj,nih)
            ENDDO
          ENDDO
        ENDDO  !nv
      ENDDO   !nis

C     New second order derivatives

      IF(NITB_S.EQ.1) THEN  ! line slave mesh(no 2nd order derivatives)
        ndimm=0
      ELSEIF (NITB_S.EQ.2) THEN !area slave mesh one 2nd order derivative)
        ndimm=1
      ELSEIF (NITB_S.EQ.3) THEN !volume slave mesh(three 2nd order derivatives)
        ndimm=3
      ENDIF

      DO ndim=1,ndimm

        nk=INT((-ndim*ndim+7*ndim+2+0.1)/2.0)     ! 4,6,7
        ii=INT((-3*ndim*ndim+11*ndim-4+0.1)/2.0)  ! 2,3,1
        DO nv=1,nvt_S
          DO njj=1,nj_max_H
            nj_field=NJ_LOC(NJL_FIEL,njj,NRE(ne1))  !deformed host mesh coord.
            nj=NJ_LOC(NJL_GEOM,njj,NRE(ne1))        !undeformed host mesh coord.

            AO(1,2,njj)=0.0d0
            AN(1,2,njj)=0.0d0
            DO i=1,NITB_H
              nu=3*i    !3,6,9
              AO(1,2,njj)=AO(1,2,njj)+PXI(IBT(1,1,nbf_host),
     '          IDO(1,1,0,nbf_host),INP(1,1,nbf_host),nbf_host,nu,
     '          XIP(1,np2),XE(1,nj))*dxidn(i,ii,nv)
              AN(1,2,njj)=AN(1,2,njj)+PXI(IBT(1,1,nbf_host),
     '          IDO(1,1,0,nbf_host),INP(1,1,nbf_host),nbf_host,nu,
     '          XIP(1,np2),XE(1,nj_field))*dxidn(i,ii,nv)
            ENDDO
            AO(2,2,njj)=0.0d0
            AN(2,2,njj)=0.0d0
            DO i=1,NITB_H
              nu=3*i*i-10*i+13  !6,5,10
              AO(2,2,njj)=AO(2,2,njj)+PXI(IBT(1,1,nbf_host),
     '          IDO(1,1,0,nbf_host),INP(1,1,nbf_host),nbf_host,nu,
     '          XIP(1,np2),XE(1,nj))*dxidn(i,ii,nv)
              AN(2,2,njj)=AN(2,2,njj)+PXI(IBT(1,1,nbf_host),
     '          IDO(1,1,0,nbf_host),INP(1,1,nbf_host),nbf_host,nu,
     '          XIP(1,np2),XE(1,nj_field))*dxidn(i,ii,nv)
            ENDDO
            AO(3,2,njj)=0.0d0
            AN(3,2,njj)=0.0d0
            DO i=1,NITB_H
              nu=INT((-3*i*i+11*i+10+0.1)/2.0)  ! 9,10,8
              AO(3,2,njj)=AO(3,2,njj)+PXI(IBT(1,1,nbf_host),
     '          IDO(1,1,0,nbf_host),INP(1,1,nbf_host),nbf_host,nu,
     '          XIP(1,np2),XE(1,nj))*dxidn(i,ii,nv)
              AN(3,2,njj)=AN(3,2,njj)+PXI(IBT(1,1,nbf_host),
     '          IDO(1,1,0,nbf_host),INP(1,1,nbf_host),nbf_host,nu,
     '          XIP(1,np2),XE(1,nj_field))*dxidn(i,ii,nv)
            ENDDO
          ENDDO !njj

          DO nj=1,nj_max_H
            DO nih=1,NITB_H
              A(nj,nih)=UNDEF(nj,nih)
            ENDDO
            B(nj)=XP(nk,nv,nj,np2)-dxidn(1,ndim,nv)*AO(1,2,nj)-
     '        dxidn(2,ndim,nv)*AO(2,2,nj)-dxidn(3,ndim,nv)*AO(3,2,nj)
          ENDDO

          CALL DGETRF(NITB_H,nj_max_H,A,3,IPIV,INFO)
          CALL DGETRS('N',nj_max_H,1,A,3,IPIV,B,3,INFO)

          DO nj=1,nj_max_S
            d2xidn2(nj,ndim,ii)=B(nj)
          ENDDO

          DO nj=1,nj_max_S
            C(nj)=0.0d0
            DO nih=1,NITB_H
              C(nj)=C(nj)+DEF(nj,nih)*d2xidn2(nih,ndim,ii)
            ENDDO
          ENDDO

          DO nj=1,nj_max_H
            XP(nk,nv,nj,np2)=C(nj)+dxidn(1,ndim,nv)*AN(1,2,nj)+
     '        dxidn(2,ndim,nv)*AN(2,2,nj)+dxidn(3,ndim,nv)*AN(3,2,nj)
          ENDDO
        ENDDO  !nv

      ENDDO  !ndim

      CALL EXITS('UPSLAVE')
      RETURN
 9999 CALL ERRORS('UPSLAVE',ERROR)
      CALL EXITS('UPSLAVE')
      RETURN 1
      END


CC AJPs - 191297 - rgb
