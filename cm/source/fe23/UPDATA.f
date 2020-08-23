      SUBROUTINE UPDATA(FD,IBT,IDO,INP,ISIZE_MFI,LD,LD_NP,NAN,NBH,NBJ,
     &  NBJF,NDDATA,NEELEM,NELIST,NFF,NHE,NKEF,NKHE,NKJE,NLL,NNF,NPF,
     &  NPNE,NPNF,NVHE,NVJE,NVJF,NRE,NRLIST,NW,Z_CONT_LIST,CURVCORRECT,
     &  MFI,PG,RG,SE,SF,SQ,WD,WG,XA,XAB,XE,XG,XID,XP,YQS,ZA,Z_CONT,ZD,
     &  ZE,ZG,ZP,STRING,ERROR,*)

C#### Subroutine: UPDATA
C###  Description:
C###    UPDATA updates data weights or fields.
C###  See-Also: COMPDAT

C LKC  2-JUL-1999 Adding updata from fields evaluated from XP
C DPN  4-DEC-2003 Adding updata from grid fields

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nonl00.cmn'      
      INCLUDE 'trsf00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER FD(NDM),IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     &  INP(NNM,NIM,NBFM),ISIZE_MFI(3,NSSM),LD(NDM),LD_NP(NDM),
     &  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     &  NDDATA(0:NDM,0:NRM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NFF(6,NEM),NHE(NEM,NXM),NKEF(0:4,16,6,NBFM),
     &  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     &  NNF(0:17,6,NBFM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     &  NRE(NEM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NW(NEM,3,NXM),
     &  Z_CONT_LIST(NDM,2,7),NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),MFI(NDM,NTSM,3,NSSM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     &  SQ(NDM),WD(NJM,NDM),WG(NGM,NBM),XA(NAM,NJM,NEM),XAB(NORM,NEM),
     &  XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     &  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),ZD(NJM,NDM),
     &  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER ERR,I,IBEG,ICNT,IEND,II,IP,ISTRIPE_END(2),J,N3CO,nb,
     '  ncount,nd,ne,nef,nf,ni,nj,njj,NJ_FROM,NJ_TO,np,nr,nt,NT_FIELD,
     &  nss,nx
      REAL*8 A(2,2),AVANG,COSAVANG,det_A,det_m,DIST,DOT,GGQ,
     &  inv_A(2,2),inv_m(2,2),m(2,2),MAXDIST,NORMAL(3),
     '  POS(3),projected(3),SUMANG,SMINDIST,SMINOLD,TAN_DERI(3,2,2),
     '  TANGENT(3,2),TIME,THRESH,VEC_MAG,WCON1,WCON2,WCON3,WCON4,XFF1,
     '  XFF2,XI(3)
      CHARACTER TYPE*10,TYPE2*16
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,DEFORMED,INWARD,OUTWARD

!     Functions
      INTEGER CALC_SAMPLE_FROM_TIME,IFROMC
      REAL*8 PXI,RFROMC

      CALL ENTERS('UPDATA',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update data
C###  Parameter:      <(weights)[weights]>
C###  Description: sets array WD according to whether data
C###      obtained from polylines is oriented primarily vertically
C###      (at an angle of [80,100] degrees) or horizontally (at
C###      an angle of [-10,10] degrees) or neither
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers which data belongs to.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<(weights)[weights]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update data field
C###  Description:
C###    Updates the fields associated with data points
C###  Parameter:      <from (projection/gap/MFI/YQS)[projection]>
C###    PROJECTION:
C###    Updates the data with field values which are the orthogonal
C###    distance of the data to the surface. If only distances
C###    which are inside (inward) or outside (outward) [as determined
C###    by the normal to the surface] are required in the field
C###    position of the data, then the appropriate option should
C###    be specified.
C###    GAP:
C###    Updates the data with field values which are the gap function
C###    for contact mechanics. The projection distances are +ve inside
C###    the volume and -ve outside [as determined by the normal to the
C###    surface].
C###    MFI:
C###    Updates the data with field values from the MFI array.
C###    YQS:
C###    Updates the data with field values from the YQS array.
C###  Parameter:      <to slave/master>
C###    SLAVE:
C###    Updates Z_CONT,Z_CONT_LIST with the contact points information
C###    for the slave surface, the gap values contained in ZD, contact elements
C###    and contact XI locations are updated.
C###    MASTER:
C###    Updates Z_CONT,Z_CONT_LIST with the contact points information
C###    for the master surface, the gap values contained in ZD, contact elements
C###    and contact XI locations are updated.
C###  Parameter:      <inward/outward>
C###    For PROJECTION:
C###    Only updates data if on the inside or outside of a mesh.
C###  Parameter:      <nss #[1]>
C###    For MFI: The MFI index to update data from.
C###    For YQS: The YQS index to update data from. Currently only one
C###             index can be given, so things won't work too well if
C###             you want to use different indices for different
C###             variants.
C###  Parameter:      <nj #[1]>
C###    For YQS: The index of the data field to be updated
C###  Parameter:      <time #[0.0]>
C###    For MFI: The time to update into the data fields.
C###  Parameter:      <frequency #[IPTRSF]>
C###    For MFI: The frequency to map between time and time-step
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers which data belongs to.

        OP_STRING(1)=STRING(1:IEND)//' field'
        OP_STRING(2)=BLANK(1:15)//
     '    '<from (projection/gap/MFI/YQS) [projection]>'
        OP_STRING(3)=BLANK(1:15)//'<to slave/master>'
        OP_STRING(4)=BLANK(1:15)//'<inward/outward>'
        OP_STRING(5)=BLANK(1:15)//'<nss #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<nj #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<time #[0.0]>'
        OP_STRING(8)=BLANK(1:15)//'<frequency #[IPTRSF]>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update data field
C###  Description:
C###    Updates the fields associated with data points
C###  Parameter:      <volume/coupled_nodes>
C###    Updates the field with volumes calculated from undeformed or
C###    deformed geometry, or with a field that is stored at a set of
C###    coupled nodes.        
C###  Parameter:      <undeformed/deformed [deformed]>
C###    UNDEFORMED:
C###    Updates the data point field with the mean volume = total
C###    volume of the region / number of data points (for 'volume'
C###    qualifier'); updates the data point field with undeformed field
C###    value in XP (for 'coupled_nodes' qualifier).
C###    DEFORMED: Updates the data point field with deformed volumes
C###    calculated using the data point xi locations as sample points
C###    (for 'volume' qualifier) .
C###  Parameter:      <from_field_num (#)[1]>
C###    Field that stores the undeformed volumes.
C###  Parameter:      <to_field_num (#)[1]>
C###    Field that stores the deformed volumes (will overwrite any
C###    current values).
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers which data belongs to.

        OP_STRING(1)=STRING(1:IEND)//' field'
        OP_STRING(2)=BLANK(1:15)//'<volume>'
        OP_STRING(3)=BLANK(1:15)//'<undeformed/deformed [deformed]>'
        OP_STRING(4)=BLANK(1:15)//'<from_field_num (#)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<to_field_num (#)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update data datafield
C###  Description:
C###    Updates the data with field values stored in XP.
C###    Assumes that the xi positions have already been calculated.
C###    This is designed from comparing field values before and
C###    after a field fit. An additional field must be defined to
C###    store the new data to avoid overwriting
C###    the field data.
C###  Parameter:      <from_field_num (#)[1]>
C###    The field to update
C###  Parameter:      <to_field_num (#)[1]>
C###    Which field number to updata into. The default
C###    overwrites the old data.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers which field data belongs to..

        OP_STRING(1)=STRING(1:IEND)//' datafield'
        OP_STRING(2)=BLANK(1:15)//'<from_field_num (#)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<to_field_num (#)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPDATA',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'WEIGHTS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='WEIGHTS'
        ELSE IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN

C LKC 5-FEB-2002 This should not be essential
C          CALL ASSERT(CALL_ELEM,'>> no elements -define elements first',
C         '      ERROR,*9999)
          
          CALL ASSERT(CALL_DATA,'>> no data or xi positions',
     '      ERROR,*9999)
          TYPE='FIELD'
          IF(CBBREV(CO,'VOLUME',3,noco+1,NTCO,N3CO)) THEN
            TYPE2='VOLUME'
          ELSE IF(CBBREV(CO,'COUPLED_NODES',3,noco+1,NTCO,N3CO)) THEN
            TYPE2='CPNODES'
          ELSE IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'PROJECTION',2)) THEN
              TYPE2='PROJECTION'
            ELSE IF(ABBREV(CO(N3CO+1),'GAP',2)) THEN
              TYPE2='GAP'
            ELSE IF(ABBREV(CO(N3CO+1),'MFI',3)) THEN
              TYPE2='MFI'
            ELSE IF(ABBREV(CO(N3CO+1),'YQS',3)) THEN
              TYPE2='YQS'
            ENDIF
            OUTWARD=.FALSE.
            INWARD=.FALSE.
            IF(CBBREV(CO,'OUTWARD',2,noco+1,NTCO,N3CO)) THEN
              OUTWARD=.TRUE.
            ELSE IF(CBBREV(CO,'INWARD',3,noco+1,NTCO,N3CO)) THEN
              INWARD=.TRUE.
            ENDIF
          ELSE IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'SLAVE',3)) THEN
              TYPE2='SLAVE'
            ElSE IF(ABBREV(CO(N3CO+1),'MASTER',3)) THEN
              TYPE2='MASTER'
            ENDIF
          ELSE
            TYPE2='DEFAULT'
          ENDIF
        ELSEIF(CBBREV(CO,'DATAFIELD',5,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
          TYPE2='DATAFIELD'
        ELSE
          ERROR='Unknown Type'
          GOTO 9999
        ENDIF

        IF(TYPE2(1:6).EQ.'VOLUME'.OR.TYPE2(1:7).EQ.'CPNODES')THEN
          IF(CBBREV(CO,'FROM_FIELD_NUM',5,noco+1,NTCO,N3CO)) THEN
            NJ_FROM=IFROMC(CO(N3CO+1))
          ELSE
            NJ_FROM=1
          ENDIF
          IF(CBBREV(CO,'TO_FIELD_NUM',2,noco+1,NTCO,N3CO)) THEN
            NJ_TO=IFROMC(CO(N3CO+1))
          ELSE
            NJ_TO=1
          ENDIF
          IF(CBBREV(CO,'UNDEFORMED',5,noco+1,NTCO,N3CO)) THEN
            DEFORMED=.FALSE.
          ELSE
            DEFORMED=.TRUE. !default
          ENDIF
          
        ENDIF

        IF(TYPE2(1:9).EQ.'DATAFIELD') THEN
          IF(CBBREV(CO,'FROM_FIELD_NUM',5,noco+1,NTCO,N3CO)) THEN
            NJ_FROM=IFROMC(CO(N3CO+1))
          ELSE
            NJ_FROM=1
          ENDIF

          IF(CBBREV(CO,'TO_FIELD_NUM',3,noco+1,NTCO,N3CO)) THEN
            NJ_TO=IFROMC(CO(N3CO+1))
          ELSE
            NJ_TO=1
          ENDIF

        ELSEIF(TYPE2(1:3).EQ.'MFI') THEN

          IF(CBBREV(CO,'TIME',3,noco+1,NTCO,N3CO)) THEN
            TIME=RFROMC(CO(N3CO+1))
          ELSE
            TIME=0.D0
          ENDIF

C Over-ride the frequency possibly set in the iptrsf file
          IF(CBBREV(CO,'FREQUENCY',3,noco+1,NTCO,N3CO)) THEN
            TRSF_FREQUENCY=RFROMC(CO(N3CO+1))
          ENDIF

          IF(CBBREV(CO,'NSS',3,noco+1,NTCO,N3CO)) THEN
            nss=IFROMC(CO(N3CO+1))
          ELSE
            nss=1
          ENDIF

        ELSEIF(TYPE2(1:3).EQ.'YQS') THEN

          CALL ASSERT(NQT.LE.NDM,'>> Increase NDM to NQM',
     &      ERROR,*9999)
          
          IF(CBBREV(CO,'NSS',3,noco+1,NTCO,N3CO)) THEN
            nss=IFROMC(CO(N3CO+1))
          ELSE
            nss=1
          ENDIF

          CALL ASSERT(nss.LE.NIQST,'>> Invalid nss',
     &      ERROR,*9999)

          IF(CBBREV(CO,'NJ',2,noco+1,NTCO,N3CO)) THEN
            NJ_TO=IFROMC(CO(N3CO+1))
          ELSE
            NJ_TO=1
          ENDIF

          CALL ASSERT(NJ_TO.LE.NJM,'>> Invalid nj',
     &      ERROR,*9999)

        ENDIF !TYPE2

C LKC f-JUL-1999 Generalise for regions
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL ASSERT(NRLIST(0).EQ.1,
     '    '>> Only implemented for 1 region',ERROR,*9999)
        nr=NRLIST(1)


        IF(TYPE(1:5).EQ.'FIELD') THEN
C news AJP 23/7/97
          IF(TYPE2(1:10).EQ.'PROJECTION') THEN
            IF(NJ_LOC(NJL_FIEL,0,1).GT.0) THEN
              NT_FIELD=NJ_LOC(NJL_FIEL,0,1)
            ELSE
              NT_FIELD=1
            ENDIF

C LKC 4-JUL-1999 Generalise for regions
C            nr=1 !temporary
            DO nd=1,NDT
              !Find whether data is inside or outside.  Firstly find
              !normal at xid
              ne=LD(nd) !element #
              nef=FD(nd) ! local face #
              IF(ne.GT.0) THEN
                XI(1)=XID(1,nd)
                XI(2)=XID(2,nd)
                CALL GET_TNVECTOR2(IBT,IDO,INP,NBJ,NBJF,ne,nef,NFF,
     '            NKEF,NKJE,NNF,NPF,NPNE,NPNF,nr,NVJE,NVJF,NORMAL,POS,
     '            SE,SF,TAN_DERI,TANGENT,XA,XE,XI,XP,ERROR,*9999)
                DOT=0.0d0
                DO nj=1,NJT
                  DOT=DOT+NORMAL(nj)*(ZD(nj,nd)-POS(nj))
                ENDDO !dot prod. of normal vector and data proj. vector
                DO nj=NJT+1,NJT+NT_FIELD
                  IF (DOT.GT.0) THEN
                    ZD(nj,nd)=DSQRT(SQ(nd))
                    IF(.NOT.INWARD) THEN
                      WD(nj,nd)=1.0d0
                    ELSE
!                      WD(nj,nd)=0.0d0
                      WD(nj,nd)=1.0d0 !test AJP to get better
                      ZD(nj,nd)=0.0d0 !fit
                    ENDIF
                  ELSE
                    ZD(nj,nd)=-DSQRT(SQ(nd))
                    IF(.NOT.OUTWARD) THEN
                      WD(nj,nd)=1.0d0
                    ELSE
!                      WD(nj,nd)=0.0d0
                      WD(nj,nd)=1.0d0 !test AJP to get better
                      ZD(nj,nd)=0.0d0 !fit
                    ENDIF
                  ENDIF
                ENDDO !nj
              ELSE !ne is 0
                DO nj=NJT+1,NJT+NT_FIELD
                  ZD(nj,nd)=0.0d0
                  WD(nj,nd)=0.0d0
                ENDDO !nj
              ENDIF !ne
            ENDDO !nd

C JWF added 14/12/01
C Updates data field with gap function for contact mechanics.
C +ve inside volume, -ve outside volume

          ELSEIF(TYPE2(1:3).EQ.'GAP') THEN

C*** 21/02/08 JHC Add check on size of NJM
            CALL ASSERT(NJM.GE.67,
     '      '>>Increase NJM to 67 for general CONTACT problems',
     '      ERROR,*9999)

            IF(NJ_LOC(NJL_FIEL,0,1).GT.0) THEN
              NT_FIELD=NJ_LOC(NJL_FIEL,0,1)
            ELSE
              NT_FIELD=1
            ENDIF

            DO nd=1,NDT

              ne=LD(nd) ! element #
              nef=FD(nd) ! local face #
              nf=NFF(nef,ne) ! global face #
              nr=NRE(ne) ! region #
              XI(1)=XID(NPF(1,nf),nd)
              XI(2)=XID(NPF(3,nf),nd)

              CALL GET_TNVECTOR2(IBT,IDO,INP,NBJ,NBJF,ne,nef,NFF,
     '          NKEF,NKJE,NNF,NPF,NPNE,NPNF,nr,NVJE,NVJF,NORMAL,POS,
     '          SE,SF,TAN_DERI,TANGENT,XA,XE,XI,XP,ERROR,*9999)
              DO nj=1,NJT
                projected(nj)=(POS(nj)-ZD(nj,nd))
              ENDDO !nj

C             Ensure that normal comes out of volume
              nb=NBJ(1,ne) !basis # of element
              IF (XID(NNF(1,nef,nb),nd).EQ.0.0d0) THEN ! Xi direction normal to face
                DO nj=1,NJT
                  NORMAL(nj)=-1.0d0*NORMAL(nj) !normal is inwards so reverse
                ENDDO !nj
              ENDIF

              DOT=0.0d0
              DO nj=1,NJT
                DOT=DOT+projected(nj)*NORMAL(nj)
              ENDDO !dot prod. of normal vector and data proj. vector
              ZD(13,nd)=DOT !norm_gap
              WD(13,nd)=1.0d0

C*** 29/09/08 JHC Removed the if statement below. Tied contact now follows the same eq. as frictional contact
CC*** 21/02/08 JHC Added if statements for tied contact which uses normalised tangent vectors                             
C              IF (Z_CONT_LIST(nd,1,4).EQ.2) THEN ! tied contact
CC               Normalise tangent_1 vector as it is not normalised in GET_TNVECTOR2.f
C                VEC_MAG=0.0d0
C                VEC_MAG=SQRT((TANGENT(1,1)**2)+(TANGENT(2,1)**2)+
C     '            (TANGENT(3,1)**2))
C                TANGENT(1,1)=TANGENT(1,1)/VEC_MAG
C                TANGENT(2,1)=TANGENT(2,1)/VEC_MAG      
C                TANGENT(3,1)=TANGENT(3,1)/VEC_MAG     
C                DOT=0.0d0
C                DO nj=1,NJT
C                  DOT=DOT+projected(nj)*TANGENT(nj,1)
C                ENDDO !dot prod. of tangent(1) vector and data proj. vector
C                ZD(14,nd)=DOT !tang_gap1
C                WD(14,nd)=1.0d0
CC               To create orthogonal coordinate system must use normal(nj) X tangent(nj,1)
CC               instead of tangent(nj,2). Therefore overwrite old TANGENT(nj,2).
C               CALL cross(NORMAL,TANGENT(1,1),TANGENT(1,2))
C               Normalise new tangent_2 vector
C                VEC_MAG=0.0d0
C                VEC_MAG=SQRT((TANGENT(1,2)**2)+(TANGENT(2,2)**2)+
C     '            (TANGENT(3,2)**2))
C                TANGENT(1,2)=TANGENT(1,2)/VEC_MAG
C                TANGENT(2,2)=TANGENT(2,2)/VEC_MAG      
C                TANGENT(3,2)=TANGENT(3,2)/VEC_MAG      
C                DOT=0.0d0
C                DO nj=1,NJT
C                  DOT=DOT+projected(nj)*TANGENT(nj,2)
C                ENDDO !dot prod. of tangent(2) vector and data proj. vector
C                ZD(15,nd)=DOT !tang_gap2
C                WD(15,nd)=1.0d0
C              ELSE ! for frictionless or frictional contact, ZD(14,nd) and ZD(15,nd) are not used

              ZD(14,nd)=0.0d0 !tang_gap1
              WD(14,nd)=0.0d0

              ZD(15,nd)=0.0d0 !tang_gap2
              WD(15,nd)=0.0d0

C              ENDIF
C NEWE JHC                                                                                                                  
              IF(ne.GT.0) THEN 
                DO nj=1,NJT
                  ZD(nj+3,nd)=NORMAL(nj)
                  ZD(nj+6,nd)=TANGENT(nj,1)
                  ZD(nj+9,nd)=TANGENT(nj,2)
C*** 21/02/08 JHC Add tengent derivative vectors to ZD
                  ZD(nj+30,nd)=TAN_DERI(nj,1,1)
                  ZD(nj+33,nd)=TAN_DERI(nj,1,2)
                  ZD(nj+36,nd)=TAN_DERI(nj,2,1)
                  ZD(nj+39,nd)=TAN_DERI(nj,2,2)

                  WD(nj+3,nd)=1.0d0
                  WD(nj+6,nd)=1.0d0
                  WD(nj+9,nd)=1.0d0

                  WD(nj+30,nd)=1.0d0
                  WD(nj+33,nd)=1.0d0
                  WD(nj+36,nd)=1.0d0
                  WD(nj+39,nd)=1.0d0
                ENDDO !nj

C*** 21/02/08 JHC Compute m, the covariant metic tensor of tangent vectors, 
C                 inverse of m, the contraviant metric tensor of tangent vectors
C                 and inverse of A. Refer to Computational Contact and Impact mechancis by
C                 T.Laursen for what they represent
                           
                DO i=1,2
                  DO j=1,2
                    m(i,j)=0.0d0
                    A(i,j)=0.0d0
                    inv_m(i,j)=0.0d0
                    inv_A(i,j)=0.0d0
                    DOT=0.0d0
                    DO nj=1,NJT
                      m(i,j)=m(i,j)+TANGENT(nj,i)*TANGENT(nj,j)
                      DOT=DOT+TAN_DERI(nj,i,j)*NORMAL(nj)
                    ENDDO !nj
                    A(i,j)=m(i,j)+ZD(13,nd)*DOT ! ND(13,nd) is the normal gap
                  ENDDO
                ENDDO

C*** 21/02/08 JHC Compute contravariant tensor inv_m of metric tensor, m          

                det_m=m(1,1)*m(2,2)-m(1,2)*m(2,1)
                IF (DABS(det_m).LE.ZERO_TOL) THEN
                  WRITE(OP_STRING,'('' >>Warning: determinant of '
     &              //'covariant tensor m at nd='',I5,'
     &              //''' is less than zero => '
     &              //'Division By Zero'')') nd
                ELSE
                  inv_m(1,1)=(1/det_m)*m(2,2)
                  inv_m(1,2)=-(1/det_m)*m(1,2)
                  inv_m(2,1)=-(1/det_m)*m(2,1)
                  inv_m(2,2)=(1/det_m)*m(1,1)
                ENDIF

C*** 21/02/08 JHC Compute inverse of A       
                det_A=A(1,1)*A(2,2)-A(1,2)*A(2,1)
                IF (DABS(det_A).LE.ZERO_TOL) THEN
                  WRITE(OP_STRING,'('' >>Warning: determinant of '
     &              //'A at nd='',I5,'
     &              //''' is less than zero => '
     &              //'Division By Zero'')') nd
                ELSE
                  inv_A(1,1)=(1/det_A)*A(2,2)
                  inv_A(1,2)=-(1/det_A)*A(1,2)
                  inv_A(2,1)=-(1/det_A)*A(2,1)
                  inv_A(2,2)=(1/det_A)*A(1,1)
                ENDIF
C*** 21/02/08 JHC Store m, inverse of m and inverse of A into ZD
                ZD(47,nd)=m(1,1)
                ZD(48,nd)=m(1,2)
                ZD(49,nd)=m(2,1)
                ZD(50,nd)=m(2,2)

                WD(47,nd)=1.0d0
                WD(48,nd)=1.0d0
                WD(49,nd)=1.0d0
                WD(50,nd)=1.0d0
      
                ZD(51,nd)=inv_m(1,1)
                ZD(52,nd)=inv_m(1,2)
                ZD(53,nd)=inv_m(2,1)
                ZD(54,nd)=inv_m(2,2)

                WD(51,nd)=1.0d0
                WD(52,nd)=1.0d0
                WD(53,nd)=1.0d0
                WD(54,nd)=1.0d0

                ZD(55,nd)=inv_A(1,1)
                ZD(56,nd)=inv_A(1,2)
                ZD(57,nd)=inv_A(2,1)
                ZD(58,nd)=inv_A(2,2)

                WD(55,nd)=1.0d0
                WD(56,nd)=1.0d0
                WD(57,nd)=1.0d0
                WD(58,nd)=1.0d0

              ELSE !ne is 0
                DO nj=1,NJT
                  ZD(nj+3,nd)=0.0d0
                  ZD(nj+6,nd)=0.0d0
                  ZD(nj+9,nd)=0.0d0

                  ZD(nj+30,nd)=0.0d0
                  ZD(nj+33,nd)=0.0d0
                  ZD(nj+36,nd)=0.0d0
                  ZD(nj+39,nd)=0.0d0

                  WD(nj+3,nd)=0.0d0
                  WD(nj+6,nd)=0.0d0
                  WD(nj+9,nd)=0.0d0

                  WD(nj+30,nd)=0.0d0
                  WD(nj+33,nd)=0.0d0
                  WD(nj+36,nd)=0.0d0
                  WD(nj+39,nd)=0.0d0
                ENDDO !nj

                ZD(47,nd)=0.0d0
                ZD(48,nd)=0.0d0
                ZD(49,nd)=0.0d0
                ZD(50,nd)=0.0d0

                WD(47,nd)=0.0d0
                WD(48,nd)=0.0d0
                WD(49,nd)=0.0d0
                WD(50,nd)=0.0d0
      
                ZD(51,nd)=0.0d0
                ZD(52,nd)=0.0d0
                ZD(53,nd)=0.0d0
                ZD(54,nd)=0.0d0

                WD(51,nd)=0.0d0
                WD(52,nd)=0.0d0
                WD(53,nd)=0.0d0
                WD(54,nd)=1.0d0

                ZD(55,nd)=0.0d0
                ZD(56,nd)=0.0d0
                ZD(57,nd)=0.0d0
                ZD(58,nd)=0.0d0

                WD(55,nd)=0.0d0
                WD(56,nd)=0.0d0
                WD(57,nd)=0.0d0
                WD(58,nd)=0.0d0
              ENDIF !ne

C             Modification for friction/tied contact
C*** 29/09/08 JHC Tied contact uses the same equation as the frictional contact
C             Remember master surface Xi coordinates at beginning of load step for 
C             calculation of distance slipped in friction/tied problems
              IF ((Z_CONT_LIST(nd,1,4).EQ.3).OR. !friction problem 
     &          (Z_CONT_LIST(nd,1,4).EQ.2)) THEN !tied contact problem
C*** 21/02/08 JHC Removed below in order to implement Coulomb's frictional law                                 
C                IF (LOAD_IT.EQ.1) THEN !beginning of load step                                      
C                  Z_CONT(nd,1,23)=XI(1)
C                  Z_CONT(nd,1,24)=XI(2)       
C                ENDIF
C                XI(1)=Z_CONT(nd,1,23)
C                XI(2)=Z_CONT(nd,1,24)
C                CALL GET_TNVECTOR2(IBT,IDO,INP,NBJ,NBJF,ne,nef,NFF,
C     '            NKEF,NKJE,NNF,NPF,NPNE,NPNF,nr,NVJE,NVJF,NORMAL,POS,
C     '            SE,SF,TAN_DERI,TANGENT,XA,XE,XI,XP,ERROR,*9999)
C                DO nj=1,NJT
C                  projected(nj)=(POS(nj)-ZD(nj,nd))
C                ENDDO !nj 
C                                                                        
CC               Recalculate tangential distance slipped since beginning of load step.
C                DOT=0.0d0
C                DO nj=1,NJT
C                  DOT=DOT+projected(nj)*ZD(nj+6,nd)
C                ENDDO !dot prod. of tangent(1) vector and data proj. vector
C                ZD(14,nd)=DOT
C                WD(14,nd)=1.0d0
C                                               
C                DOT=0.0d0
C                DO nj=1,NJT
C                  DOT=DOT+projected(nj)*ZD(nj+9,nd)
C                ENDDO !dot prod. of tangent(2) vector and data proj. vector
C                ZD(15,nd)=DOT
C                WD(15,nd)=1.0d0

C*** 21/02/08 JHC CONV_IT should be used instead of LOAD_IT as the slip distance for each augmentation needs to be 
C                 recalculated
C                IF (LOAD_IT.EQ.1) THEN
                IF (CONV_IT.EQ.1) THEN
C                 At the beginning of the each load step, the previous xi location is
C                 set to be the same as the starting point                                  
                  ZD(59,nd)=XI(1)
                  ZD(60,nd)=XI(2)   
                  ZD(43,nd)=0.0d0
                  ZD(44,nd)=0.0d0

                ELSE

C                 evaluate relative movement of xi locations 
                  ZD(43,nd)=XI(1)-ZD(59,nd) 
                  ZD(44,nd)=XI(2)-ZD(60,nd) 

                ENDIF
      
              ENDIF !friction
                                                                                                                                                         
            ENDDO !nd                                                                                                  

C JWF added 14/06/02
C Initialise and update Z_CONT,Z_CONT_LIST with slave info

          ELSEIF(TYPE2(1:5).EQ.'SLAVE') THEN

            DO nd=1,NDT

              Z_CONT_LIST(nd,1,1)=LD(nd)   ! element of projection
              Z_CONT_LIST(nd,1,2)=FD(nd)   ! local face #

              DO i=1,3
                Z_CONT(nd,1,i)=XID(i,nd) ! XI locations
              ENDDO

              Z_CONT(nd,1,16)=1.0d0 ! Jacobian

            ENDDO !nd

C JWF added 14/06/02
C Initialise and update Z_CONT,Z_CONT_LIST with master info

          ELSEIF(TYPE2(1:6).EQ.'MASTER') THEN

            DO nd=1,NDT

              Z_CONT_LIST(nd,2,1)=LD(nd)   ! element of projection
              Z_CONT_LIST(nd,2,2)=FD(nd)   ! local face #

              DO i=1,3
                Z_CONT(nd,2,i)=XID(i,nd) ! XI locations
              ENDDO

              DO i=4,6  ! normal components
                Z_CONT(nd,1,i)=ZD(i,nd) ! towards slave
                Z_CONT(nd,2,i)=-ZD(i,nd) ! towards master 
              ENDDO !i

              DO i=7,9  ! tangent(1) components
                Z_CONT(nd,1,i)=ZD(i,nd) ! towards slave
                Z_CONT(nd,2,i)=-ZD(i,nd) ! towards master 
              ENDDO !i

              DO i=10,12  ! tangent(2) components
                Z_CONT(nd,1,i)=ZD(i,nd) ! towards slave
                Z_CONT(nd,2,i)=-ZD(i,nd)! towards master 
              ENDDO !i

              DO i=13,15
                Z_CONT(nd,1,i)=ZD(i,nd) ! slave
                Z_CONT(nd,2,i)=ZD(i,nd) ! master
              ENDDO !i signed gap magnitudes normal/tangential

              Z_CONT(nd,2,16)=1.0d0 ! Jacobian

C*** 21/02/08 JHC Adding tangent derivative vectors, m, inverse of m and A to Z_CONT 
              DO i=31,33  ! derivative of tangent(1) components wrt xi1
                Z_CONT(nd,1,i)=ZD(i,nd) ! towards slave
                Z_CONT(nd,2,i)=ZD(i,nd) ! towards master 
              ENDDO !i
              
              DO i=34,36  ! derivative of tangent(1) components wrt xi2
                Z_CONT(nd,1,i)=ZD(i,nd) ! towards slave
                Z_CONT(nd,2,i)=ZD(i,nd) ! towards master 
              ENDDO !i

              DO i=37,39  ! derivative of tangent(2) components wrt xi1
                Z_CONT(nd,1,i)=ZD(i,nd) ! towards slave
                Z_CONT(nd,2,i)=ZD(i,nd) ! towards master 
              ENDDO !i

              DO i=40,42  ! derivative of tangent(2) components wrt xi2
                Z_CONT(nd,1,i)=ZD(i,nd) ! towards slave
                Z_CONT(nd,2,i)=ZD(i,nd) ! towards master 
              ENDDO !i

              DO i=43,44  ! relative movements of xi locations between each Newton step
                Z_CONT(nd,1,i)=ZD(i,nd) 
                Z_CONT(nd,2,i)=ZD(i,nd) 
              ENDDO !i

C             Z_CONT(nd,1,45),Z_CONT(nd,2,45)
C             Z_CONT(nd,1,46),Z_CONT(nd,2,46) are used to store Phi_trial and contact pressure for augmentation
C                                             See fe50/CONTACT_RESIDUAL.f

              DO i=47,50  ! covariant mectic tensor of tangent vectors in the current coordinates (m)
                Z_CONT(nd,1,i)=ZD(i,nd)
                Z_CONT(nd,2,i)=ZD(i,nd)
              ENDDO !i

              DO i=51,54  ! contravariant mectic tensor of tangent vectors in the deformed coordinates (inv_m)
                Z_CONT(nd,1,i)=ZD(i,nd)
                Z_CONT(nd,2,i)=ZD(i,nd)
              ENDDO !i

              DO i=55,58  ! inverse of A (inv_A)
                Z_CONT(nd,1,i)=ZD(i,nd)
                Z_CONT(nd,2,i)=ZD(i,nd)
              ENDDO !i

            ENDDO !nd

          ELSEIF(TYPE2(1:3).EQ.'MFI') THEN

            CALL ASSERT(NJM.GE.NJT+NJT,
     '        '>> Increase NJM to NJT*2',ERROR,*9999)

C LKC 17-FEB-2002 Check that the number of data points present is
C           correct.
            
            CALL ASSERT(NDDATA(0,nr).EQ.ISIZE_MFI(1,nss),
     '        '>> Inconsistent number of data points and MFI values',
     '        ERROR,*9999)

            
C 16-APR-2002 Put in correct mapping between time step and real time
C            nt=INT(TIME)+1
            nt=CALC_SAMPLE_FROM_TIME(TIME,ERR,ERROR)

            WRITE(OP_STRING,
     &        '('' Updating data at MFI index '',I5,'' & time '',F8.2)')
     &        nt,TIME            
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

            IF(nt.GT.ISIZE_MFI(2,nss)) THEN ! error check for array bounds
              IEND=0
              CALL APPENDC(IEND,
     &          ' Invalid time sample of nt = ',ERROR)              
              CALL APPENDI(IEND,nt,ERROR)              
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              GOTO 9999
            ENDIF

            
            DO nd=1,NDDATA(0,nr)
              DO nj=1,NJT
                ZD(NJT+nj,nd)=MFI(nd,nt,nj,nss)
              ENDDO
            ENDDO

          ELSEIF(TYPE2(1:3).EQ.'YQS') THEN

            NJ_TO=NJ_LOC(NJL_FIEL,NJ_TO,nr)
            DO nd=1,NQT
              ZD(NJ_TO,nd)=YQS(nss,nd)
            ENDDO

          ELSEIF(TYPE2(1:9).EQ.'DATAFIELD') THEN

C LKC 2-JUL-1999 adding updata data from a field stored in XP
C*** Update data field (not projections)

C*** Check fields are appropriate
            CALL ASSERT(NJ_FROM.LE.NJ_LOC(NJL_FIEL,0,nr),
     '        '>> Invalid field FROM',ERROR,*9999)
            CALL ASSERT(NJ_TO.LE.NJ_LOC(NJL_FIEL,0,nr),
     '        '>> Invalid field TO',ERROR,*9999)
            CALL ASSERT(NJ_FROM.GT.0,'>>Field FROM is -ve',ERROR,*9999)
            CALL ASSERT(NJ_TO.GT.0,'>> Field TO is -ve',ERROR,*9999)


            IEND=0
            CALL APPENDC(IEND,' Updating field #',OP_STRING(1))
            CALL APPENDI(IEND,NJ_FROM,OP_STRING(1))
            CALL APPENDC(IEND,' to field #',OP_STRING(1))
            CALL APPENDI(IEND,NJ_TO,OP_STRING(1))
            CALL APPENDC(IEND,' for region ',OP_STRING(1))
            CALL APPENDI(IEND,nr,OP_STRING(1))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)


C*** Set the actual field numbers
            NJ_FROM=NJ_LOC(NJL_FIEL,NJ_FROM,nr)
            NJ_TO=NJ_LOC(NJL_FIEL,NJ_TO,nr)

            DO nd=1,NDDATA(0,nr)
              CALL XPXE(NBJ(1,LD(nd)),NKJE(1,1,1,LD(nd)),
     '          NPF(1,1),NPNE(1,1,LD(nd)),
     '          NRE(LD(nd)),NVJE(1,1,1,LD(nd)),
     '          SE(1,1,LD(nd)),XA(1,1,LD(nd)),XE,XP,ERROR,*9999)

              DO ni=1,NIT(NBJ(NJ_TO,LD(nd)))
                XI(ni)=XID(ni,nd)
              ENDDO

C*** Interpolate to obtain the field at the data point -
C***  the resulting field is added as a new field in ZD
              nb=NBJ(NJ_TO,LD(nd))
              ZD(NJ_TO,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '          INP(1,1,nb),nb,1,XI,XE(1,NJ_FROM))
              WD(NJ_TO,nd)=1.D0
            ENDDO !nd
            OP_STRING(1)=' '
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          ELSE IF(TYPE2(1:7).EQ.'CPNODES')THEN
C*** Update the data field (to_field) from nodal field (from_field)
C***  Check validity and set field numbers
            CALL ASSERT(NJ_FROM.LE.NJ_LOC(NJL_FIEL,0,nr),
     &        '>> FROM field not defined',ERROR,*9999)
            CALL ASSERT(NJ_TO.LE.NJ_LOC(NJL_FIEL,0,nr),
     &        '>> TO field not defined',ERROR,*9999)
            NJ_FROM=NJ_LOC(NJL_FIEL,NJ_FROM,nr)
            NJ_TO=NJ_LOC(NJL_FIEL,NJ_TO,nr)
            DO nd=1,NDT
              np=LD_NP(nd) !coupled node number
              ZD(nj_to,nd)=XP(1,1,nj_from,np)
              WD(nj_to,nd)=1.d0
            ENDDO !nd
          ELSE IF(TYPE2(1:6).EQ.'VOLUME')THEN
c*** Calculate undeformed or deformed volumes
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &        ERROR,*9999)
            nx=1 !temporary
            CALL UPDATA_VOLUME(IBT,IDO,INP,LD,LD_NP,NAN,NBH,NBJ,
     &        NELIST,NFF,NHE(1,nx),NJ_FROM,NJ_TO,NKHE,NKJE,NLL,NPF,
     &        NPNE,nr,NRE,NVHE,NVJE,NW(1,1,nx),nx,CURVCORRECT,PG,RG,
     &        SE,WG,XA,XAB,XE,XG,XID,XP,ZA,ZD,ZE,ZG,ZP,DEFORMED,ERROR,
     &        *9999)
          ELSE IF(TYPE2(1:7).EQ.'DEFAULT')THEN
C MHT 7-Sept-2004 update the number of fields for data points
C set the weightings to 1.0            
            DO nd=1,NDT
              DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
                nj=NJ_LOC(NJL_FIEL,njj,nr)
                WD(nj,nd)=1.d0
              ENDDO !njj
            ENDDO !nd
          ELSE
            ERROR='TYPE2 unknown - update code'
            GOTO 9999
          ENDIF ! TYPE2
          CALL_DATA_FIELD=.TRUE.


        ELSE IF(TYPE(1:7).EQ.'WEIGHTS') THEN
          IF(CBBREV(CO,'PROJECT',3,noco+1,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'MAXDIST',3,noco+1,NTCO,N3CO)) THEN
              MAXDIST=RFROMC(CO(N3CO+1))
            ELSE
              MAXDIST=1.d8 !i.e. no limit on the distance
            ENDIF
            ncount=0
            DO nd=1,NDT
              IF(SQ(nd).GT.MAXDIST**2)THEN
                DO nj=1,NJT
                  WD(nj,nd)=0.d0
                ENDDO
                SQ(nd)=0.d0
                ncount=ncount+1
              ENDIF
            ENDDO
            WRITE(OP_STRING,'(''Set weighting on '',I6,'
     &        //''' data to 0.0'')') ncount
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

           ELSE
C*****This section sets array WD according to whether data
C         obtained from polylines is oriented primarily vertically
C         (at an angle of [80,100] degrees) or horizontally (at
C         an angle of [-10,10] degrees) or neither

          SMINOLD=999999
          DO I=2,NDT
            SMINDIST=(ZD(1,I)-ZD(1,I-1))**2.d0+(ZD(2,I)-ZD(2,I-1))**2.d0
            IF(SMINDIST.LT.SMINOLD.AND.SMINDIST.GE.1.0D-6)
     '        SMINOLD=SMINDIST
          END DO  !SMINOLD is the smallest squared dist. between pts
                  !with the exception of cases where we purposely made 2 pts
                  !lie on exactly the same position

C ... Table ISTRIPE_END(*) indicates if point ZD is last in a polyline
          DO I=2,NDT
            DIST=(ZD(1,I)-ZD(1,I-1))**2.d0+(ZD(2,I)-ZD(2,I-1))**2.d0
            THRESH=25.d0*SMINOLD
            IF(DIST.GT.THRESH) THEN
              ISTRIPE_END(I-1)=1
            ELSE
              ISTRIPE_END(I-1)=0
            ENDIF
          END DO
          ISTRIPE_END(NDT)=1

          WCON1=DCOS(80.d0/180.d0*PI)      !COS(80 DEG)
C !!!! AJP 7/10/95 SHOULD BE LINE BELOW ??? WCON2=-WCON2
          WCON2=-WCON1
          WCON3=DCOS(10.d0/180.d0*PI)      !COS(10 DEG)
          WCON4=-WCON3

          ICNT=1
443       II=1           !counter for use in setting weights
          IP=1           !counter for use in dividing SUMANG
          SUMANG=0.d0
445       ICNT=ICNT+1
          IF (ICNT.GT.NDT) GO TO 446
          GGQ=((ZD(1,ICNT)-ZD(1,ICNT-1))**2.d0
     '        +(ZD(2,ICNT)-ZD(2,ICNT-1))**2.d0)**0.5d0
          IF (GGQ.NE.0.d0) THEN
            XFF1=(ZD(1,ICNT)-ZD(1,ICNT-1))/GGQ
            XFF2=(ZD(2,ICNT)-ZD(2,ICNT-1))/GGQ
            SUMANG=SUMANG+(DATAN2(XFF2,XFF1))*180.d0/PI
            IP=IP+1
          ENDIF
          II=II+1
          IF (ISTRIPE_END(ICNT).EQ.1) THEN
            AVANG=SUMANG/DBLE(IP-1)
            COSAVANG=DCOS(AVANG*PI/180.d0)
            IF (COSAVANG.GT.WCON2.AND.COSAVANG.LT.WCON1) THEN  !verticle
              DO J=ICNT-II+1,ICNT
                WD(1,J)=1.d0
                WD(2,J)=0.d0
              ENDDO
            ELSE IF (COSAVANG.GT.WCON3.OR.COSAVANG.LT.WCON4) THEN   !horizontal
              DO J=ICNT-II+1,ICNT
                WD(1,J)=0.d0
                WD(2,J)=1.d0
              ENDDO
            ELSE   !warning
              PRINT*,'Stripe is neither horizontal or verticle:'
              PRINT*,'Weights for this stripe have been left unchanged.'
            ENDIF
            ICNT=ICNT+1
            GO TO 443
          ELSE
            GO TO 445
          ENDIF
446       I=I      !End of section
        ENDIF
       ENDIF
      ENDIF

      CALL EXITS('UPDATA')
      RETURN
 9999 CALL ERRORS('UPDATA',ERROR)
      CALL EXITS('UPDATA')
      RETURN 1
      END


