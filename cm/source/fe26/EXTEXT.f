      SUBROUTINE EXTEXT(NEELEM,NELIST,NQNE,NQS,NQXI,NRLIST,
     '  YQ,STRING,ERROR,*)

C#### Subroutine: EXTEXT
C###  Description:
C###    EXTEXT exports textures from the grid pt data base.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'nqloc00.inc'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM)
      REAL*8 YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IFROMC,
     '  Material,n,na,ne,nee,niq,N3CO,no_nelist,no_nrlist,
     '  N_cells,N_cells_Xi1,N_cells_Xi2,N_cells_Xi3,
     '  n1,n2,n3,N_points,N_pts_Xi1,N_pts_Xi2,N_pts_Xi3,
     '  N_graphical_materials,nq,nr,nx,SCHEME,
     '  Xi1_min,Xi2_min,Xi3_min,
     '  Xi1_max,Xi2_max,Xi3_max
      REAL*8 Isovalue,SCALAR(20)
      CHARACTER CHAR4*4,CHAR6*6,
     '  FILE*100,FILENAME*100,TEXTURE_NAME*50
      LOGICAL ALL_REGIONS,CBBREV,FOUND,NAME_ASSIGNED

      CALL ENTERS('EXTEXT',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM export textures<;FILENAME[default]>
C###  Parameter:      <region #[1]>
C###  Parameter:      <region (#s/all)[1]>
C###    Limit to the specified regions.
C###  Parameter:      <element (GROUP/#s/all)[all]>
C###    Specify the elements.
C###  Parameter:      <index #[1]>
C###    Specify the niq index in the YQ array for the desired grid
C###    variable.
C###  Parameter:      <level #[1]>
C###    Specify the grid level.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <as NAME[0..0]>
C###    Label the file with a character name.
C###  Description:
C###    Export textures from the grid point data base.

        OP_STRING(1)=STRING(1:IEND)
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<element (GROUP/#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<index #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<level #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<as NAME[element#]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXTEXT',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'INDEX',1,noco+1,NTCO,N3CO)) THEN
          niq=IFROMC(CO(N3CO+1))
        ELSE
          CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq,NIQ_V,ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'LEVEL',1,noco+1,NTCO,N3CO)) THEN
          na=IFROMC(CO(N3CO+1))
        ELSE
          na=1
        ENDIF

        IF(CBBREV(CO,'CLASS',1,noco+1,NTCO,N3CO)) THEN
          nx=IFROMC(CO(N3CO+1))
        ELSE
          nx=1
        ENDIF

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          TEXTURE_NAME=CO(N3CO+1)(IBEG1:IEND1)
          NAME_ASSIGNED=.TRUE.
        ELSE
          NAME_ASSIGNED=.FALSE.
        ENDIF
C GMH 26/12/96 NAME_ASSIGNED is not used
        NAME_ASSIGNED=NAME_ASSIGNED
C GMH 26/12/96 TEXTURE_NAME is not used
        TEXTURE_NAME=TEXTURE_NAME

        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)

          Xi1_min=0
          Xi2_min=0
          Xi3_min=0
          Xi1_max=1
          Xi2_max=1
          Xi3_max=1

          DO nee=1,NEELEM(0,nr)
            ne=NEELEM(nee,nr)
            SCHEME=NQS(ne)

            N_cells_Xi1=max(NQXI(1,SCHEME)-1,1) !#cells in Xi1 dir.n
            N_cells_Xi2=max(NQXI(2,SCHEME)-1,1) !#cells in Xi2 dir.n
            N_cells_Xi3=max(NQXI(3,SCHEME)-1,1) !#cells in Xi3 dir.n
            N_cells=N_cells_Xi1*N_cells_Xi2*N_cells_Xi3
            CALL ASSERT(N_cells_Xi1.LE.20,
     '        'Increase dimension of SCALAR',ERROR,*9999)

            N_pts_Xi1=N_cells_Xi1+1         !#points in Xi1 dir.n
            N_pts_Xi2=N_cells_Xi2+1         !#points in Xi2 dir.n
            N_pts_Xi3=N_cells_Xi3+1         !#points in Xi3 dir.n
            N_points=N_pts_Xi1*N_pts_Xi2*N_pts_Xi3

            N_graphical_materials=1
            Material=1

            FOUND=.FALSE.
            DO no_nelist=1,NELIST(0)
              IF(NELIST(no_nelist).EQ.ne) FOUND=.TRUE.
            ENDDO
            IF(FOUND) THEN
C **          open file with element# appended
              WRITE(CHAR6,'(I6)') ne
              CALL STRING_TRIM(CHAR6,IBEG1,IEND1)
              TEXTURE_NAME='element '//CHAR6(IBEG1:IEND1)
              CALL STRING_TRIM(FILE,IBEG,IEND)
              WRITE(CHAR4,'(I4)') 1000+ne
              FILENAME=FILE(IBEG:IEND)//'_'//CHAR4//'.extext'
              WRITE(OP_STRING,'('' Creating '',A)') FILENAME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL OPENF(IFILE,'DISK',FILENAME,'NEW','SEQUEN',
     '          'FORMATTED',132,ERROR,*9999)

C **          write the group name
c             CALL STRING_TRIM(TEXTURE_NAME,IBEG,IEND)
c             WRITE(IFILE,'('' Texture name: '',A)')
c    '          TEXTURE_NAME(IBEG:IEND)
c             WRITE(IFILE,'(1X,''Element: '',I5)') ne

              WRITE(IFILE,'(1X,3I3)') Xi1_min,Xi2_min,Xi3_min
              WRITE(IFILE,'(1X,3I3)') Xi1_max,Xi2_max,Xi3_max
              WRITE(IFILE,'(1X,3I3)') N_cells_Xi1,N_cells_Xi2,
     '          N_cells_Xi3
              WRITE(IFILE,'(1X, I3)') N_graphical_materials
              WRITE(IFILE,'(1X,'' default'')')
              WRITE(IFILE,'((1X,8I3))') (Material,n=1,N_cells)

              DO n3=1,N_cells_Xi3
                DO n2=1,N_cells_Xi2
                  DO n1=1,N_cells_Xi1 !calc cell values for 8 cells

C                    ng1=n1+(n2-1)*N_pts_Xi1 +(n3-1)*N_pts_Xi1*N_pts_Xi2
C                    ng2=ng1+1
C                    ng3=ng1+N_pts_Xi1
C                    ng4=ng3+1
C                    ng5=n1+(n2-1)*N_pts_Xi1 + n3*N_pts_Xi1*N_pts_Xi2
C                    ng6=ng5+1
C                    ng7=ng5+N_pts_Xi1
C                    ng8=ng7+1
C                    IF(DOP) THEN
C                      WRITE(*,'('' ng1..8: '',8I3)')
C     '                  ng1,ng2,ng3,ng4,ng5,ng6,ng7,ng8
C                    ENDIF
C                    SCALAR(n1)=0.125d0*(YQ(NQGE(ng1,ne,nb),niq,na,nx)
C     '                                 +YQ(NQGE(ng2,ne,nb),niq,na,nx)
C     '                                 +YQ(NQGE(ng3,ne,nb),niq,na,nx)
C     '                                 +YQ(NQGE(ng4,ne,nb),niq,na,nx)
C     '                                 +YQ(NQGE(ng5,ne,nb),niq,na,nx)
C     '                                 +YQ(NQGE(ng6,ne,nb),niq,na,nx)
C     '                                 +YQ(NQGE(ng7,ne,nb),niq,na,nx)
C     '                                 +YQ(NQGE(ng8,ne,nb),niq,na,nx))

                    SCALAR(n1)=1.0d0 !temp for now
                  ENDDO !n1
                  WRITE(IFILE,'((1X,9E12.4))')
     '              (SCALAR(n1),n1=1,N_cells_Xi1)
                ENDDO !n2
              ENDDO !n3

              WRITE(IFILE,'((1X,9I3))') (Material,n=1,N_points)
              Isovalue=1.d0
              WRITE(IFILE,'((1X,E12.4))') Isovalue
              WRITE(IFILE,'('' nodal_values'')')
              WRITE(IFILE,'((1X,9E12.4))')
     '          (YQ(NQNE(ne,nq),niq,na,nx),nq=1,N_points)

              CALL CLOSEF(IFILE,ERROR,*9999)
            ENDIF !found
          ENDDO !no_nelist
        ENDDO !no_nrlist

      ENDIF

      CALL EXITS('EXTEXT')
      RETURN
 9999 CALL ERRORS('EXTEXT',ERROR)
      CALL EXITS('EXTEXT')
      RETURN 1
      END


