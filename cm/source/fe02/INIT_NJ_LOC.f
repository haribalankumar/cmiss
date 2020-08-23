      SUBROUTINE INIT_NJ_LOC(NUMBER,TYPE,nr,ERROR,*)

C#### Subroutine: INIT_NJ_LOC
C###  Description:
C###    INIT_NJ_LOC initialises NJ_LOC in region nr.

C#### Variable: NJ_LOC(0:nj_type,0:njj,0:nr)
C###  Type: INTEGER
C###  Set_up: INIT_NJ_LOC,DECOOR,IPFIBR,IPSHEE,IPFIEL
C###  Description:
C###    <HTML><PRE>
C###    NJ_LOC(xxxx,0,nr) is the number of variables of type xxxx
C###    defined in region nr, where xxxx is either NJL_GEOM, NJL_FIBR
C###    or NJL_FIEL, constants defined in the common block file
C###    loc00.cmn.
C###    NJ_LOC(NJL_GEOM,1..,nr) are nj indices in XP etc for geom. vars
C###    NJ_LOC(NJL_FIBR,1..,nr)  "  "     "    "  "   "   "  fibre  "
C###    NJ_LOC(NJL_FIEL,1..,nr)  "  "     "    "  "   "   "  field  "
C###    NJ_LOC(xxxx,0,0) is total number of nj indicies of type xxxx
C###    defined over all regions. NJ_LOC(0,0,nr) is the total number
C###    of nj locations defined in region nr and NJ_LOC(0,0,0) is the
C###    total number of nj locations defined in all regions.
C###    </PRE></HTML>
C###  See-Also: loc00.cmn

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
!     Parameter List
      INTEGER NUMBER,TYPE,nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,njj1,njj2,nrr

      CALL ENTERS('INIT_NJ_LOC',*9999)

      CALL ASSERT(NUMBER.LE.NJ_LOC_MX,
     '  '>>Increase NJ_LOC_MX in loc00.cmn',ERROR,*9999)

      DO nj=1,NUMBER
        NJ_LOC(TYPE,nj,nr)=NJ
        NJ_TYPE(nj,1)=TYPE
        NJ_TYPE(nj,2)=NJ
      ENDDO
      NJ_LOC(TYPE,0,nr)=NUMBER
      IF(NUMBER.GT.NJ_LOC(TYPE,0,0)) NJ_LOC(TYPE,0,0)=NUMBER
      NJ_LOC(0,0,nr)=0
      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=NJ
        ENDDO
      ENDDO
      DO nrr=1,NRT
        IF(NJ_LOC(0,0,nrr).GT.NJ_LOC(0,0,0))
     '    NJ_LOC(0,0,0)=NJ_LOC(0,0,nrr)
      ENDDO !nrr

      CALL EXITS('INIT_NJ_LOC')
      RETURN
 9999 CALL ERRORS('INIT_NJ_LOC',ERROR)
      CALL EXITS('INIT_NJ_LOC')
      RETURN 1
      END


