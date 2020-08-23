      SUBROUTINE OPC50(IBT,IDO,INP,NAN,NBH,NBJ,ne,NHE,NPNE,nr,nx,
     '  CE,CG,CP,FEXT,PG,RGX,XE,XG,XIG,YG,ZE,ZG,ERROR,*)

C#### Subroutine: OPC50
C###  Description:
C###    OPC50 outputs element ne reference and fibre Gauss point
C###    stress and strain fields for current solution.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ne,NHE,NPNE(NNM,NBFM),
     '  nr,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),RGX(NGM),XE(NSM,NJM),
     '  XG(NJM,NUM),XIG(NIM,NGM),YG(NIYGM,NGM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ng,IPOINTTYP
      REAL*8 DZDX(3,3),EG(3,3),PHI(3),PST(3),
     '  R(3,3),RGX2D,RGZ,RGZ2D,RI1,RI2,RI3,RM(3,3),
     '  TC(3,3),TG(3,3),TN(3,3),TNA,U(3,3)
      CHARACTER STRESSTYPE*17

      PARAMETER(IPOINTTYP=1) !Gauss point strains/stresses

      DATA STRESSTYPE/'Total'/

      CALL ENTERS('OPC50',*9999)
      WRITE(OP_STRING,
     '        '(/'' Element'',I3,'' Gauss point solutions:'')') NE
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO ng=1,NGT(NBH(NH_LOC(1,nx)))
        CALL ZEEX50('Reference',IBT,IDO,INP,NAN,NBH,NBJ,
     '    ng,NHE,NPNE,nr,nx,
     '    DZDX,CE,CG,CP,EG,PG,PHI,PST,
     '    R,RGX(ng),RI1,RI2,RI3,RM,U,
     '    XE,XG,XIG(1,ng),ZE,ZG,ERROR,*9999)
        CALL OPEG50('Reference',IPOINTTYP,NBH,ng,nr,nx,DZDX,EG,
     '    PHI,PST,R,RI1,RI2,RI3,RM,U,XG,XIG(1,ng),ZG,.FALSE.,
     '    .FALSE.,ERROR,*9999)
        CALL ZEEX50('Fibre',IBT,IDO,INP,NAN,NBH,NBJ,
     '    ng,NHE,NPNE,nr,nx,
     '    DZDX,CE,CG,CP,EG,PG,PHI,PST,
     '    R,RGX(ng),RI1,RI2,RI3,RM,U,
     '    XE,XG,XIG(1,ng),ZE,ZG,ERROR,*9999)
        CALL OPEG50('Fibre',IPOINTTYP,NBH,ng,nr,nx,DZDX,EG,
     '    PHI,PST,R,RI1,RI2,RI3,RM,U,XG,XIG(1,ng),ZG,.FALSE.,
     '    .FALSE.,ERROR,*9999)
        CALL ZETX50('Reference','Cauchy',STRESSTYPE,IBT,IDO,INP,
     '    NAN,NBH,NBJ,ng,
     '    NHE,NPNE,nr,ne,nx,CE,CG,CP,FEXT(1,ng),PG,PHI,PST,
     '    RGX(ng),RGX2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,
     '    XE,XG,XIG(1,ng),YG(1,ng),ZE,ZG,ERROR,*9999)
        CALL OPTG50('Reference',IPOINTTYP,'Total',NBH,ng,nr,nx,
     '    PHI,PST,RM,TC,TG,TN,TNA,XG,XIG(1,ng),ZG,.TRUE.,
     '    ERROR,*9999)
        CALL ZETX50('Fibre','Cauchy',STRESSTYPE,IBT,IDO,INP,
     '    NAN,NBH,NBJ,ng,
     '    NHE,NPNE,nr,ne,nx,CE,CG,CP,FEXT(1,ng),PG,PHI,PST,
     '    RGX(ng),RGX2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,
     '    XE,XG,XIG(1,ng),YG(1,ng),ZE,ZG,ERROR,*9999)
        CALL OPTG50('Fibre',IPOINTTYP,'Total',NBH,ng,nr,nx,
     '    PHI,PST,RM,TC,TG,TN,TNA,XG,XIG(1,ng),ZG,.TRUE.,
     '    ERROR,*9999)
      ENDDO

      CALL EXITS('OPC50')
      RETURN
 9999 CALL ERRORS('OPC50',ERROR)
      CALL EXITS('OPC50')
      RETURN 1
      END


