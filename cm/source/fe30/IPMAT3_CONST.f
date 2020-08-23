      SUBROUTINE IPMAT3_CONST(nr,nx,ALL_REGIONS,ERROR,*)

C#### Subroutine: IPMAT3_CONST
C###  Description:
C###    IPMAT3_CONST inputs material parameters where the properties
C###    do not vary over the mesh.  The properties are stored in CEMAT
C###    instead of CE (see IPMAT3).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b32.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iltot00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'titl30.cmn'
!     Parameter List
      INTEGER nr,nx
      CHARACTER ERROR*(*)
      LOGICAL ALL_REGIONS
!     Local Variables
      INTEGER IB,IBEG,ICHAR,IE,IEND,il,NOQUES
      CHARACTER CHAR*100,CHAR11*11
      LOGICAL FILEIP

      CALL ENTERS('IPMAT3_CONST',*9999)

      CALL ASSERT(ITYP2(nr,nx).NE.0,
     '  '>>Solution type has not been defined',ERROR,*9999)
      CALL ASSERT(ITYP2(nr,nx).GT.2,'>>Solution type is incorrect',
     '  ERROR,*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      ILT(1,nr,nx)=ILTOT2(NJT,ITYP2(nr,nx),ITYP3(nr,nx))
      ILT(0,nr,nx)=1
      DO il=ILT(0,nr,nx),ILT(1,nr,nx)
        CHAR=TITL32(il,ITYP2(nr,nx),ITYP3(nr,nx))
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),ITYP2(nr,nx),ITYP3(nr,nx))
        WRITE(CHAR11,'(E11.4)') RDEFLT(1)
        CALL STRING_TRIM(CHAR11,IB,IE)
        FORMAT='($,'' The value of the '//CHAR(IBEG:IEND)//' is '//
     '    '['//CHAR11(IB:IE)//']: '',E11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=PULMAT(il)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,*9999)
        IF(IOTYPE.NE.3) PULMAT(il)=RDATA(1)
      ENDDO
 
C     Gas exchange - save parameters in common block variables
      IF(ITYP3(nr,nx).EQ.7)THEN
        P_ATM=PULMAT(1)                !atmospheric total pressure, mmHg 
        LUNG_ELASTANCE=PULMAT(2)       !elastance, mmHg/l       
        HCT_INITIAL=PULMAT(3)          !capillary blood hematocrit
        CAP_BLOOD_VOL=PULMAT(4)        !capillary blood volume, litre
        ALV_AIR_VOL=PULMAT(5)          !alveolar air volume, litre
        AIR_BLOOD_SURFACE=PULMAT(6)    !air-blood surface area, m2  
        BARRIER_THICKNESS=PULMAT(7)    !air-blood barrier thickness, m
      ENDIF

      IF(FILEIP.AND..NOT.ALL_REGIONS) CALL CLOSEF(IFILE,ERROR,*9999)
      CALL EXITS('IPMAT3_CONST')
      RETURN
 9999 CALL ERRORS('IPMAT3_CONST',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPMAT3_CONST')
      RETURN 1
      END



