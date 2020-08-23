      SUBROUTINE ZPRP_FACE(IBT,IDO,IDOXFT,INP,MYMS_F,NBH,NBHF,NBJ,NBJF,
     '  nc,nf,NHE,NKB,NKHE,NKJE,NNB,NNF,NPF,NPNE,NPNY,nr,NSB,NSFE,
     '  NVHE,NVJE,nx,NYNP,NYNS_E,CE,CG,CP,PG,RDF,SE,SM_F,SN_E,
     '  XDF,XP,YGF,YP,ZDF,ERROR,*)

C#### Subroutine: ZPRP_FACE
C###  Description:
C###    ZPRP_FACE calculates the contribution to global residual vector
C###    YP(my,4) at current solution YP(ny,1) from one face.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDOXFT(NHM),
     '  INP(NNM,NIM,NBFM),MYMS_F(0:NSFM,2,2,NHM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM),NBJ(NJM,NEM),NBJF(NJM),nc,nf,
     '  NHE(NEM),NKB(2,2,2,NNM,NBFM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NPF(9),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),nr,
     '  NSB(NKM,NNM,NBFM),NSFE(NSM,2,2,NHM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),NYNS_E(0:NSM,2,2,NHM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CP(NMM,NPM),
     '  PG(NSM,NUM,NGM,NBM),RDF(NSFM,2,NHM),SE(NSM,NBFM,NEM),
     '  SM_F(NSFM,2,2,NHM),SN_E(NSM,2,2,NHM),
     '  XDF(NSFM,2,2,NJM),XP(NKM,NVM,NJM,NPM),YGF(NIYGFM,NGFM),
     '  YP(NYM,NIYM),ZDF(NSFM,2,2,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER idoxf,JDOXFT,ms_f,my,
     '  NEA(2),neax,NEAXT,NF_EA(2),nhx,NHF,NIEF(0:2),NORMALSIGN1,np
      LOGICAL INTEGRATE

      CALL ENTERS('ZPRP_FACE',*9999)

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(ZPRP_FACE_1)
        WRITE(OP_STRING,'(/'' =========='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' Face '',I5)') nf
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' =========='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(ZPRP_FACE_1)
      ENDIF

      IF(nc.EQ.1) THEN !upwinding is only applied on this matrix

C***    Find connectivity information.
        CALL FACE_INT_PREP(NEA,NEAXT,NF_EA,NHE,NHF,NIEF,
     '    NPF,nr,nx,INTEGRATE,ERROR,*9999)

        IF(INTEGRATE) THEN !false if integral is known to be zero

C***      Find necessary upwinding info.
          CALL FACE_INT_LGE(IBT,IDO,IDOXFT,INP,JDOXFT,MYMS_F,
     '      NBH,NBHF(1,1),NEA,NEAXT,NF_EA,NHF,NIEF,NKB,NKHE,NORMALSIGN1,
     '      NNB,NNF,NPNE,nr,NSB,NSFE,NVHE,NYNP(1,1,1,1,0,1),NYNS_E,
     '      nx,SE,SM_F,SN_E,ERROR,*9999)

C***      Determine face geometry nodal out-of-face derivative information.
C***      Derivatives in the face are interpolated from facial nodal values.
          CALL XPXF(IBT,IDO,NBJ,NBJF,NEA,1,NF_EA,NIEF,
     '      NKB,NKJE,NNF,NPNE,nr,NSB,NVJE,SE,XDF,XP,ERROR,*9999)

C***      Calculate cross face derivatives or deriv discontinuities.
          CALL ZPZF(IDOXFT,JDOXFT,NBHF(1,1),NEAXT,NHF,nr,NSFE,
     '      NYNS_E,nx,SN_E,YP,ZDF,ERROR,*9999)

C***      Transfer material parameters to CG.
C          DO il=1,ILT(1,nr,nx)
C            CALL ASSERT(ABS(ILP(il,1,nr,nx)).EQ.1,
C     '        '>>Material parameters must be constant',ERROR,*9999)
C          ENDDO !il
C          CALL CPCG(1,NBJF(1),NPNE,nr,nx, !is NBJF(1) always what we want
C     '      CE(1,NEA(1)),CG,CGE(1,1,NEA(1)),CP,PG,ERROR,*9999)
          CALL CPCGF(NBJF(1),NEA,NF_EA,NNF,NPNE,nr,nx,CE,CG,CP,PG,
     '      ERROR,*9999)

C***      Calculate modification to facial residual vector.
          CALL ZFRF30(IDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,NORMALSIGN1,nr,nx,
     '      CG,PG,RDF,XDF,YGF,ZDF,ERROR,*9999)

C***      Assemble facial residual modifications RDF
C***      into global residual YP(*,4).

          neax=1
C KAT 20Mar01: Can't branch out of critical section.
C              Atomic directive is probably better anyway.
CC$OMP     CRITICAL(ZPRP_FACE_2)
          DO nhx=1,NHF
            DO idoxf=1,IDOXFT(nhx) !cross face derivatives
              DO ms_f=1,MYMS_F(0,idoxf,neax,nhx)
                my=MYMS_F(ms_f,idoxf,neax,nhx)
                IF(NPNY(0,my,0).EQ.1) THEN
                  np=NPNY(4,my,0)
C                 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF
C$OMP           ATOMIC
                YP(my,4)=YP(my,4)+
     '            RDF(ms_f,idoxf,nhx)*SM_F(ms_f,idoxf,neax,nhx)
              ENDDO !ms_f
            ENDDO !idoxf
          ENDDO !nhx
CC$OMP     END CRITICAL(ZPRP_FACE_2)
        ENDIF !INTEGRATE
      ENDIF !nc=1

      CALL EXITS('ZPRP_FACE')
      RETURN
 9999 CALL ERRORS('ZPRP_FACE',ERROR)
      CALL EXITS('ZPRP_FACE')
      RETURN 1
      END



