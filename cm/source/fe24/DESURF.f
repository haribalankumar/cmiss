      SUBROUTINE DESURF(IBT,IDO,INP,NBJ,NELIST,NEELEM,NLS_SURF_PSI,
     '  NLS_SURF_XI,NRLIST,STRING,ERROR,*)

C#### Subroutine: DESURF
C###  Description:
C###    DESURF calculates surface and tag subdivisions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NRLIST(0:NRM)
      REAL*8 NLS_SURF_XI(3,26,26,4),NLS_SURF_PSI(16,26,26,3,4)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER i1,i2,ibeg,iend,INDEX_SURF,INT_TEMP(1),N3CO,nb,ni,nj,nk,
     '  nn, !ISURF_ELEM(1,10), !ALISTAIR:check this: PJH 2Sept95
     '  NOXIPT(2),NT,NTXIPT
      REAL*8 DXI1,DXI2,PSI1,XI3,XI(3)
      CHARACTER VARIABLE_TYPE*8,STATUS*3
      LOGICAL ABBREV,CALCU,CBBREV,FILIO,GENER,MOUSE
C MHT 04-05-01 Declared and not referenced
C     INTEGER np,ICON_X(200),nolist
C     REAL*8 NLS_CON_PSI(16,2,1000,3,100),NLS_CON_XI(3,2,1000,100),
      CALL ENTERS('DESURF',*9999)
 1    IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define surface;c geometry
C###  Description:
C###    Calculates surface subdivisions.
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.
C###  Parameter:      <at XI_3#[0.0]>
C###  Parameter:      <surface_index #[1]>
C###  Parameter:      <by NO_PTS#[10,10]>

        OP_STRING(1)=STRING(1:IEND)//';c geometry'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<at XI_3#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<surface_index #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<by NO_PTS#[10,10]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C KAT 2001-04-10: This code is not used/working
CC#### Command: FEM define surface;c tag2d
CC###  Description:
CC###    Calculates tag subdivisions.
CC###  Parameter:      <surface_index #[1]>

C        OP_STRING(1)=STRING(1:IEND)//';c tag2d'
C        OP_STRING(2)=BLANK(1:15)//'<surface_index #[1]>'
C        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('FE24','DOC','DESURF',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS('C',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'AT',2,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSRE(CO(N3CO+1),XI3,ERROR,*9999)
        ELSE
          XI3=0.d0
        ENDIF
        IF(CBBREV(CO,'BY',2,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),2,NTXIPT,NOXIPT,ERROR,*9999)
        ELSE
          NOXIPT(1)=10
          NOXIPT(2)=10
        ENDIF
        IF(CBBREV(CO,'SURFACE_INDEX',3,NOCO+1,NTCO,N3CO)) THEN
C LKC 24-APR-1998 Need to pass in array
C          CALL PARSIL(CO(N3CO+1),1,NT,INDEX_SURF,ERROR,*9999)
          CALL PARSIL(CO(N3CO+1),1,NT,INT_TEMP,ERROR,*9999)
          INDEX_SURF=INT_TEMP(1)
        ELSE IF(ABBREV(COQU(NOCO,1),'C',1)) THEN
          INDEX_SURF=1
        ENDIF

        IF(CALCU) THEN
          IF(CBBREV(CO,'GEOMETRY',3,NOCO+1,NTCO,N3CO)) THEN
            VARIABLE_TYPE='GEOMETRY'
C KAT 2001-04-10: This code is not used/working
C          ELSE IF(CBBREV(CO,'TAG2D',5,NOCO+1,NTCO,N3CO)) THEN
C            VARIABLE_TYPE='TAG2D'
          ELSE
            VARIABLE_TYPE='GEOMETRY'
          ENDIF

          IF(VARIABLE_TYPE(1:8).EQ.'GEOMETRY') THEN
            !save number of elements in this surface

c!!! cpb 3/9/95 Alistair check this line

C            ISURF_ELEM(0,INDEX_SURF)=NELIST(0)


!            SURF_SUBDIV(1,INDEX_SURF)=NOXIPT(1)
!            SURF_SUBDIV(2,INDEX_SURF)=NOXIPT(2)
            CALL ASSERT(NELIST(0).LE.100,
     '        '>>Error: increase number of elements in /SURF04/',
     '        ERROR,*9999)
            CALL ASSERT(NOXIPT(1).LE.25,
     '        '>>Error: subdivision must be <= 25x25 /SURF04/',
     '        ERROR,*9999)
            CALL ASSERT(NOXIPT(2).LE.25,
     '        '>>Error: subdivision must be <= 25x25 /SURF04/',
     '        ERROR,*9999)

C            DO nolist=1,NELIST(0)
C              ISURF_ELEM(nolist,INDEX_SURF)=NELIST(nolist)
C            ENDDO
            XI(3)=XI3
            DXI1=1.0d0/DBLE(NOXIPT(1))
            DXI2=1.0d0/DBLE(NOXIPT(2))
            !subdivide
            DO i2=1,NOXIPT(2)+1
              XI(2)=DBLE(i2-1)*DXI2
              DO i1=1,NOXIPT(1)+1
                XI(1)=DBLE(I1-1)*DXI1
                DO ni=1,3
                  NLS_SURF_XI(ni,i1,i2,INDEX_SURF)=XI(ni)
                ENDDO
                !precalculate basis functions
                DO nj=1,NJT
                  nb=NBJ(nj,NELIST(1))
                  DO nk=1,NKT(0,nb)
                    DO nn=1,NNT(nb)
                      NLS_SURF_PSI(nk+(nn-1)*NKT(0,nb),I1,I2,nj,
     '                  INDEX_SURF)=PSI1(IBT(1,1,nb),IDO(1,1,nb),
     '                  INP(1,1,nb),nb,1,nk,nn,XI)
                    ENDDO !nn
                  ENDDO !nk
                ENDDO !nj
              ENDDO !i1
            ENDDO !i2

C KAT 2001-04-10: This code is not used/working
C          ELSE IF(VARIABLE_TYPE(1:5).EQ.'TAG2D') THEN
C            !for all tag planes
C            DO NT=1,NTAG(0)
C              !intersect surface with plane
Cc             CALL SURF_PLANE_X(NBJ,NKE,NPNE,NPF,SE,XA,XA(1,1,ne),XE,XP,
Cc    '          INDEX_SURF,NT,TAG_NORMAL(1,NT),TAG_POS(1,NT),
Cc    '          ERROR,*9999)
CC!!! ICON_X is used here before it is set
CC             precalculate basis functions
C              DO np=1,ICON_X(NT)
C                DO nj=1,NJT
C                  nb=NBJ(nj,ISURF_ELEM(1,INDEX_SURF))
C                  DO ni=1,2
C                    DO nk=1,NKT(0,nb)
C                      DO nn=1,NNT(nb)
C                        NLS_CON_PSI(nk+(nn-1)*NKT(0,nb),ni,np,nj,NT)=
C     '                    PSI1(IBT(1,1,nb),IDO(1,1,nb),INP(1,1,nb),nb,
C     '                    1,nk,nn,NLS_CON_XI(1,ni,np,NT))
C                      ENDDO !nn
C                    ENDDO !nk
C                  ENDDO !ni
C                ENDDO !nj
C              ENDDO !np
C            ENDDO !nt
          ENDIF !variable_type
        ENDIF !calcu

      ENDIF !the big if

      CALL EXITS('DESURF')
      RETURN
 9999 CALL ERRORS('DESURF',ERROR)
      CALL EXITS('DESURF')
      RETURN 1
      END


