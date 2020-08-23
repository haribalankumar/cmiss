      SUBROUTINE IPREGI(ERROR,*)

C#### Subroutine: IPREGI
C###  Description:
C###    IPREGI inputs total number of regions.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,nj,njj1,njj2,nr,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPREGI',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='($,'' The total number of regions is [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NRT
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NRM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NRT=IDATA(1)

      DO nr=1,NRT
        DO nj=1,NJT
          NJ_LOC(NJL_GEOM,nj,nr)=NJ
          NJ_TYPE(nj,1)=NJL_GEOM
          NJ_TYPE(nj,2)=NJ
        ENDDO
        NJ_LOC(NJL_GEOM,0,nr)=NJT
        IF(NJT.GT.NJ_LOC(NJL_GEOM,0,0)) NJ_LOC(NJL_GEOM,0,0)=NJT
        NJ_LOC(0,0,nr)=0
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=NJ
          ENDDO
        ENDDO
      ENDDO
      DO nr=1,NRT
        IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0))
     '    NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
      ENDDO !nrr

      CALL EXITS('IPREGI')
      RETURN
 9999 CALL ERRORS('IPREGI',ERROR)
      CALL EXITS('IPREGI')
      RETURN 1
      END


