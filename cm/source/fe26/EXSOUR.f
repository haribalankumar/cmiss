      SUBROUTINE EXSOUR(NDIPOLES,DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '  NRLIST,NXLIST,
     '  DIPOLE_DIR,DIPOLE_CEN,STRING,ERROR,*)

C#### Subroutine: EXSOUR
C###  Description:
C###    EXSOUR exports dipole sources
C Created 5-APR-2002 by LKC

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NDIPOLES(NRM),DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),
     '  NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,N3CO,ndipole,nj,
     '  nr,nt,NCENTER,NORIENT,nx,nxc,OFFSET
      REAL*8 TIME,CENTRE(4),DIRECTION(4)
      CHARACTER FILE*100,GROUPNAME*50,LINE*132
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,TIME_SET

!     Functions
      INTEGER IFROMC
      REAL*8 RFROMC

      CALL ENTERS('EXSOUR',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM export source<;FILENAME[default]>
C###  Parameter:      <as NAME[dipole]>
C###    Specify a label for the group
C###  Parameter:      <offset OFFSET[10000]>
C###    Specify the offset of the labels
C###  Parameter:      <region #[all]>
C###    Specify the region number.
C###  Parameter:      <using (fit/solve)[solve]>
C###    Specify the nx type to use.
C###  Parameter:      <time #>
C###    Specify a specific time to export. If not specified export all
C###    time points. Will find the integer time point closest to, but
C###    not greater than, the specified time.
C###  Parameter:      <class #[1]>
C###    Specify the class number.
C###  Description:
C###    Add OFFSET to data points reference numbers when

        OP_STRING(1)=STRING(1:IEND)//
     '    '<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<as NAME[dipole]>'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<offset #[10000]>'
        OP_STRING(5)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(6)=BLANK(1:15)//'<time #>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXSIGN',ERROR,*9999)
      ELSE

        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)

        IF(CBBREV(CO,'AS',2,noco,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          GROUPNAME=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          GROUPNAME='DIPOLE'
        ENDIF

        IF(CBBREV(CO,'OFFSET',2,noco+1,NTCO,N3CO)) THEN
          OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          OFFSET=10000
        ENDIF

        IF(CBBREV(CO,'TIME',4,noco+1,NTCO,N3CO)) THEN
          TIME=RFROMC(CO(N3CO+1))
          TIME_SET=.TRUE.
        ELSE
          TIME_SET=.FALSE.
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '        ERROR,*9999)
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ENDIF
        ELSE
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ENDIF

        CALL ASSERT(NDIPOLES(nr).GE.1,
     '    '>> No dipoles to export in this region',ERROR,*9999)

        CALL STRING_TRIM(FILE,IBEG,IEND)
        CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.exdata','NEW',
     '    'SEQUEN','FORMATTED',132,ERROR,*9999)

        CALL STRING_TRIM(GROUPNAME,IBEG,IEND)
        WRITE(IFILE,'( '' Group name: '',A)') GROUPNAME(IBEG:IEND)


        ndipole=1 !assume all the dipoles are structured the same

        IF(DIPOLE_CEN_NTIME(ndipole,nr,nx).GT.1) THEN
          NCENTER=DIPOLE_CEN_NTIME(ndipole,nr,nx)
        ELSE
          NCENTER=1
        ENDIF

        IF(DIPOLE_DIR_NTIME(ndipole,nr,nx).GT.1) THEN
          NORIENT=DIPOLE_DIR_NTIME(ndipole,nr,nx)
        ELSE
          NORIENT=1
        ENDIF

        IF(TIME_SET) THEN
          NCENTER=1
          NORIENT=1
        ENDIF
        
        IEND=0
        CALL APPENDC(IEND,' #Fields=',LINE)
        CALL APPENDI(IEND,NCENTER+NORIENT,LINE)
        WRITE(IFILE,'(A)') LINE(1:IEND)

C*** Write the center headers
        DO nt=1,NCENTER
          IEND=0
          CALL APPENDC(IEND,' ',LINE)
          CALL APPENDI(IEND,nt,LINE)
          CALL APPENDC(IEND,') center.',LINE)
          CALL APPENDI(IEND,nt,LINE)
          IF(NJ_LOC(NJL_GEOM,0,nr).EQ.1) THEN
            CALL APPENDC(IEND,', coordinate, '//
     '        'rectangular cartesian, #Components=1',LINE)
            WRITE(IFILE,'(A)') LINE(1:IEND)
            WRITE(IFILE,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
          ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
            CALL APPENDC(IEND,', coordinate, rectangular cartesian'//
     '        ', #Components=2',LINE)
            WRITE(IFILE,'(A)') LINE(1:IEND)
            WRITE(IFILE,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
            WRITE(IFILE,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
          ELSE
            CALL APPENDC(IEND,', coordinate, rectangular cartesian'//
     '        ', #Components=3',LINE)
            WRITE(IFILE,'(A)') LINE(1:IEND)
            WRITE(IFILE,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
            WRITE(IFILE,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
            WRITE(IFILE,'(1X,''  z.  Value index= 3, #Derivatives=0'')')
          ENDIF
        ENDDO

C*** Write the orientation headers
        DO nt=1,NORIENT
          IEND=0
          CALL APPENDC(IEND,' ',LINE)
          CALL APPENDI(IEND,nt+NCENTER,LINE)
          CALL APPENDC(IEND,') orient.',LINE)
          CALL APPENDI(IEND,nt,LINE)

          IF(NJ_LOC(NJL_GEOM,0,nr).EQ.1) THEN
            CALL APPENDC(IEND,
     '        ', coordinate, rectangular cartesian, #Components=1',LINE)
            WRITE(IFILE,'(A)') LINE(1:IEND)
            WRITE(IFILE,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
          ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
            CALL APPENDC(IEND,
     '        ', coordinate, rectangular cartesian, #Components=2',LINE)
            WRITE(IFILE,'(A)') LINE(1:IEND)
            WRITE(IFILE,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
            WRITE(IFILE,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
          ELSE
            CALL APPENDC(IEND,
     '        ', coordinate, rectangular cartesian, #Components=3',LINE)
            WRITE(IFILE,'(A)') LINE(1:IEND)
            WRITE(IFILE,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
            WRITE(IFILE,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
            WRITE(IFILE,'(1X,''  z.  Value index= 3, #Derivatives=0'')')
          ENDIF
        ENDDO

        DO ndipole=1,NDIPOLES(nr)
          WRITE(IFILE,'(1X,''Node: '',I9)') ndipole+OFFSET

          IF(TIME_SET) CALL GETDIPOLE(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     &      ndipole,nr,CENTRE,DIPOLE_CEN,DIPOLE_DIR,DIRECTION,TIME,
     &      ERROR,*9999)
          
C***      Center
          IF(NCENTER.EQ.1) THEN !not time dependent
            IF(TIME_SET) THEN
              WRITE(IFILE,'(1X,3E13.5)') (CENTRE(nj),nj=1,
     &          NJ_LOC(NJL_GEOM,0,nr))
            ELSE
              WRITE(IFILE,'(1X,3E13.5)')
     '          (DIPOLE_CEN(nj,0,ndipole,nr,nx),
     '          nj=1,NJ_LOC(NJL_GEOM,0,nr))
            ENDIF
          ELSE
            DO nt=1,NCENTER
              WRITE(IFILE,'(1X,3E13.5)')
     '          (DIPOLE_CEN(nj,nt,ndipole,nr,nx),
     '          nj=1,NJ_LOC(NJL_GEOM,0,nr))
            ENDDO
          ENDIF

C***      Orientation
          IF(NORIENT.EQ.1) THEN !not time dependent
            IF(TIME_SET) THEN
              WRITE(IFILE,'(1X,3E13.5)') 
     '          (DIRECTION(nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
            ELSE
              WRITE(IFILE,'(1X,3E13.5)') 
     '          (DIPOLE_DIR(nj,0,ndipole,nr,nx),
     '          nj=1,NJ_LOC(NJL_GEOM,0,nr))
            ENDIF
          ELSE
            DO nt=1,NORIENT
              WRITE(IFILE,'(1X,3E13.5)') 
     '          (DIPOLE_DIR(nj,nt,ndipole,nr,nx),
     '          nj=1,NJ_LOC(NJL_GEOM,0,nr))
            ENDDO
          ENDIF
        ENDDO

        CALL CLOSEF(IFILE,ERROR,*9999)

      ENDIF
      CALL EXITS('EXSOUR')
      RETURN

 9999 CALL ERRORS('EXSOUR',ERROR)
      CALL EXITS('EXSOUR')
      RETURN 1
      END




