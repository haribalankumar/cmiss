      SUBROUTINE ASSEMBLE5_FACE(IBT,IDO,IDOXFT,INP,ISC_GK,ISR_GK,MYMS_F,
     '  NBH,NBHF,NBJ,NBJF,nf,NHE,NKB,NKHE,NKJE,NNB,NNF,NPF,NPNE,
     '  nr,NSB,NSFE,NVHE,NVJE,nx,NYNP,NYNS_E,CE,CG,CP,FS,GK,
     '  PG,RDF1,RDF2,SE,SM_F,SN_E,XDF,XP,YGF,YP,ZDF,ERROR,*)

C#### Subroutine: ASSEMBLE5_FACE
C###  Description:
C###    ASSEMBLE5_FACE calculates the contribution to global tangent
C###    stiffness matrix GK from one face.


      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDOXFT(NHM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),
     '  MYMS_F(0:NSFM,2,2,NHM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM),
     '  NBJ(NJM,NEM),NBJF(NJM),nf,NHE(NEM),
     '  NKB(2,2,2,NNM,NBFM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NPF(9),NPNE(NNM,NBFM,NEM),nr,
     '  NSB(NKM,NNM,NBFM),NSFE(NSM,2,2,NHM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),NYNS_E(0:NSM,2,2,NHM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CP(NMM,NPM),
     '  FS(NSFM,2,NSFM,2,2,NHM),GK(NZ_GK_M),PG(NSM,NUM,NGM,NBM),
     '  RDF1(NSFM,2,NHM),RDF2(NSFM,2,NHM),
     '  SE(NSM,NBFM,NEM),SM_F(NSFM,2,2,NHM),SN_E(NSM,2,2,NHM),
     '  XDF(NSFM,2,2,NJM),XP(NKM,NVM,NJM,NPM),
     '  YGF(NIYGFM,NGFM),YP(NYM,NIYM),ZDF(NSFM,2,2,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER idoxf,iFS,ijdoxf,jdoxf,JDOXFT,
     '  mb_f,mh,mhx,ms_f,MSTB,my,NEA(2),neaFS,NEAFST,neax,NEAXT,
     '  NF_EA(2),NHF,NIEF(0:2),NORMALSIGN1,ns_e,ns_f,ny,nz
      LOGICAL FLUX,INTEGRATE

      CALL ENTERS('ASSEMBLE5_FACE',*9999)

C***  Find connectivity, geometry, and material information.
      CALL FACE_INT_PREP(NEA,NEAXT,NF_EA,NHE,NHF,NIEF,
     '  NPF,nr,nx,INTEGRATE,ERROR,*9999)

      IF(INTEGRATE) THEN !false if integral is known to be zero

C***    Find necessary upwinding info.
        CALL FACE_INT_LGE(IBT,IDO,IDOXFT,INP,JDOXFT,MYMS_F,
     '    NBH,NBHF(1,1),NEA,NEAXT,NF_EA,NHF,NIEF,NKB,NKHE,NORMALSIGN1,
     '    NNB,NNF,NPNE,nr,NSB,NSFE,NVHE,NYNP(1,1,1,1,0,1),NYNS_E,nx,
     '    SE,SM_F,SN_E,ERROR,*9999)

C***    Determine face geometry nodal out-of-face derivative information.
C***    Derivatives in the face are interpolated from facial nodal values.
        CALL XPXF(IBT,IDO,NBJ,NBJF,NEA,1,NF_EA,NIEF,
     '    NKB,NKJE,NNF,NPNE,nr,NSB,NVJE,SE,XDF,XP,ERROR,*9999)

C***    Calculate cross face derivatives or deriv discontinuities.
        CALL ZPZF(IDOXFT,JDOXFT,NBHF(1,1),NEAXT,NHF,nr,NSFE,
     '    NYNS_E,nx,SN_E,YP,ZDF,ERROR,*9999)

C***    Transfer material parameters to CG.
C        DO il=1,ILT(1,nr,nx)
C          CALL ASSERT(ABS(ILP(il,1,nr,nx)).EQ.1,
C     '      '>>Material parameters must be constant',ERROR,*9999)
C        ENDDO !il
C        CALL CPCG(1,NBJF(1),NPNE,nr,nx, !is NBJF(1) always what we want
C     '    CE(1,NEA(1)),CG,CGE(1,1,NEA(1)),CP,PG,ERROR,*9999)
          CALL CPCGF(NBJF(1),NEA,NF_EA,NNF,NPNE,nr,nx,CE,CG,CP,PG,
     '      ERROR,*9999)

C***  Calculate face stiffness matrix.
        CALL ZFFS(IDOXFT,JDOXFT,NBHF(1,1),NBJF,NEAXT,NHF,NIEF,
     '    NORMALSIGN1,nr,nx,CG,FS,PG,RDF1,RDF2,XDF,YGF,ZDF,ERROR,*9999)

C!!!  need to check order of stabilizing derivs
        FLUX=ITYP15(nr,nx).NE.3.OR..TRUE.

C***    Assemble element stiffness matrix into global system.

C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP   CRITICAL(ASSEMBLE5_FACE_1)
        DO mhx=1,NHF

          IF(IWRIT4(nr,nx).GE.5) THEN
            mh=NH_LOC(mhx,nx)
            mb_f=NBHF(mh,1)
            MSTB=NST(mb_f)
            WRITE(OP_STRING,
     '        '('' Face '',I5,'' stiffness matrix FS:'')') nf
            IF(FLUX) THEN
              NEAFST=NEAXT
            ELSE
              NEAFST=1 !derivs are the same in each neighbouring element
            ENDIF
            DO neax=1,NEAFST
              DO jdoxf=1,JDOXFT
                DO idoxf=1,IDOXFT(mhx)
                  DO ms_f=1,MSTB
                    WRITE(OP_STRING,
     '                '('' FS('',I2,'','',I1,'',1..,'',I1,'','','//
     '                'I1,'','',I1,'')='','//'8D12.4:/(18X,8D12.4))')
     '                ms_f,idoxf,jdoxf,neax,mhx,
     '                (FS(ms_f,idoxf,ns_f,jdoxf,neax,mhx),ns_f=1,MSTB)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          DO neax=1,NEAXT
            IF(FLUX) THEN
              neaFS=neax
            ELSE
              neaFS=1
            ENDIF
            DO idoxf=1,IDOXFT(mhx) !cross face derivs of weight func
              DO ms_f=1,MYMS_F(0,idoxf,neax,mhx)
                my=MYMS_F(ms_f,idoxf,neax,mhx)
                DO jdoxf=1,JDOXFT !cross face derivs of interp func
                  IF(FLUX) THEN !Flux Term
                    ijdoxf=jdoxf !SN_E depends on interp func deriv
                    iFS=idoxf !FS depends on weight func deriv
                  ELSE !High order derivative discontinuity term
                    ijdoxf=idoxf !SN_E depends on weight func deriv
                    iFS=1 !FS does not depend on weight func deriv
                  ENDIF
                  DO ns_e=1,NYNS_E(0,jdoxf,neax,mhx)
                    ny=NYNS_E(ns_e,jdoxf,neax,mhx)
                    IF(ny.NE.0) THEN
                      CALL SPARSE(my,ny,NYT(1,1,nx),nz,NZ_GK_M,
     '                  NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                      ns_f=NSFE(ns_e,jdoxf,neax,mhx)
C$OMP                 ATOMIC
                      GK(nz)=GK(nz)+FS(ms_f,iFS,ns_f,jdoxf,neaFS,mhx)*
     '                  SM_F(ms_f,idoxf,neax,mhx)*
     '                  SN_E(ns_e,ijdoxf,neax,mhx)
                    ENDIF !ny!=0
                  ENDDO !ns_e
                ENDDO !jdoxf
              ENDDO !ms_f
            ENDDO !idoxf
          ENDDO !neax
        ENDDO !mhx
CC$OMP   END CRITICAL(ASSEMBLE5_FACE_1)

      ENDIF !INTEGRATE

      CALL EXITS('ASSEMBLE5_FACE')
      RETURN
 9999 CALL ERRORS('ASSEMBLE5_FACE',ERROR)
      CALL EXITS('ASSEMBLE5_FACE')
      RETURN 1
      END


