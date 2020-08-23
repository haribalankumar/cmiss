      SUBROUTINE UPFLUX(FIX,NPNODE,NRLIST,NYNP,STRING,XP,YP,ERROR,*)

C#### Subroutine: UPFLUX
C###  Description:
C###    Updates boundary condition fluxes from a constant, field, solution,
C###      or material array. The flux can be substituted, added, subtracted,
C###      multiplied or divided by the array. The array can be multiplied by
C###      a scale factor before it operates on the fluxes.


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER FLUXBC,IBEG,IEND,IFROMC,
     '  N3CO,nh1,nhc1,nj1,njf2,nonode,np,nr,nrc1,
     '  ny_fluxbc,nx,nxc1,PART2
      REAL*8 CONST,RFROMC,SCALE
      LOGICAL ALL_REGIONS,CBBREV
      CHARACTER OPERATION*16,TYPE*16,WITH*8

      CALL ENTERS('UPFLUX',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update material
C###  Parameter:      <parameter #[1]>
C###    Specify the material parameter number to update.
C###  Parameter:      <depvar #[1]>
C###    Specify the dependent variable number to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number to update.
C###  Parameter:      <region #[1]>
C###    Specify the region number to update.
C###  Parameter:      <at (nodes/*elements/*gauss_pts/*grid_pts[nodes]>
C###    Specify how to update material parameters
C###  Parameter:      substitute/add/substract/multiply/divide
C###    Specify how to update the material parameters
C###  Parameter:          <constant VALUE#>
C###    Specify the constant to update the material parameters with
C###  Parameter:          <field FIELD#>
C###    Specify the field to update the material parameters with
C###  Parameter:          <*solution>
C###    Specify the solution to update the material parameters with
C###  Parameter:              <yp_index #[1]>
C###    Specify the index in YP to update the material parameters
C###  Parameter:              <depvar #[1]>
C###    Specify the dependent variable number to update the material parameters
C###  Parameter:              <class #[1]>
C###    Specify the class number to update the material parameters
C###  Parameter:              <region #[1]>
C###    Specify the region number to update the material parameters
C###  Parameter:          <*material>
C###    Specify the solution to update the material parameters with
C###  Parameter:              <parameter #[1]>
C###    Specify the material parameter number to update the material parameters
C###  Parameter:              <depvar #[1]>
C###    Specify the dependent variable number to update the material parameters
C###  Parameter:              <class #[1]>
C###    Specify the class number to update the material parameters
C###  Parameter:              <region #[1]>
C###    Specify the region number to update the material parameters
C###  Parameter:      <scale_factor VALUE[1.0]>
C###    Specify how to update material parameters
C###  Description: Updates the material parameters with
C###    a constant, field, solution (YP array), field (XP array)
C###    or material (CP/CE/CG array). These array can substitute, add,
C###    subtract, multiply or divide the material parameters. The arrays
C###    can be scaled before the material parameters are updated.
C###    The update can occur on the node, element, or gauss material
C###    parameter array.
C###    Use 'fem list material' the list the material parameters
C###    to operate on.
C###    NOTE: A number of these features have not been implemented. These
C###    features will be implements as they are needed.The '*' indicate
C###    the features which have not been implemented.
C###

        OP_STRING(1)=BLANK(1:IEND)//'<depvar #[1]>'
        OP_STRING(2)=BLANK(1:IEND)//'<class #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<region #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'substitute/add/substract/'
     '    //'multiply/divide'
        OP_STRING(5)=BLANK(1:IEND)//'    <constant VALUE#>'
        OP_STRING(6)=BLANK(1:IEND)//'    <field FIELD#>'
        OP_STRING(7)=BLANK(1:IEND)//'    <*solution>'
        OP_STRING(8)=BLANK(1:IEND)//'         <yp_index #[1]>'
        OP_STRING(9)=BLANK(1:IEND)//'         <depvar #[1]>'
        OP_STRING(10)=BLANK(1:IEND)//'         <class #[1]>'
        OP_STRING(11)=BLANK(1:IEND)//'         <region #[1]>'
        OP_STRING(12)=BLANK(1:IEND)//'    <*material>'
        OP_STRING(13)=BLANK(1:IEND)//'         <parameter #[1]>'
        OP_STRING(14)=BLANK(1:IEND)//'         <depvar #[1]>'
        OP_STRING(15)=BLANK(1:IEND)//'         <class #[1]>'
        OP_STRING(16)=BLANK(1:IEND)//'         <region #[1]>'
        OP_STRING(17)=BLANK(1:IEND)//'<scale_factor VALUE#[1.0]>'
        OP_STRING(18)=BLANK(1:IEND)//'* - not implemented'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPFLUX',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'SUBSTITUTE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='SUBSTITUTE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'ADD',1,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='ADD'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'SUBTRACT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='SUBTRACT'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'MULTIPLY',1,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='MULTIPLY'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'DIVIDE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='DIVIDE'
          PART2=N3CO
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(TYPE(1:9).EQ.'OPERATION') THEN

          IF(CBBREV(CO,'DEPVAR',2,noco+1,PART2,N3CO)) THEN
            nhc1=IFROMC(CO(N3CO+1))
          ELSE
            nhc1=1
          ENDIF

          IF(CBBREV(CO,'CLASS',2,noco+1,PART2,N3CO)) THEN
            nxc1=IFROMC(CO(N3CO+1))
          ELSE
            nxc1=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc1,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'Invalid class',ERROR,*9999)

          IF(CBBREV(CO,'REGION',2,noco+1,PART2,N3CO)) THEN
            nrc1=IFROMC(CO(N3CO+1))
          ELSE
            nrc1=1
          ENDIF

          IF(CBBREV(CO,'CONSTANT',2,PART2,NTCO,N3CO))THEN
            WITH='CONSTANT'
            CONST=RFROMC(CO(N3CO+1))
          ELSEIF(CBBREV(CO,'FIELD',2,PART2,NTCO,N3CO))THEN
            WITH='FIELD'
            njf2=IFROMC(CO(N3CO+1))
          ELSEIF(CBBREV(CO,'SOLUTION',2,PART2,NTCO,N3CO))THEN
            WITH='SOLUTION'
            ERROR='>> ''SOLUTION'' feature not implemented'
            GOTO 9999
          ELSEIF(CBBREV(CO,'MATERIAL',2,PART2,NTCO,N3CO))THEN
            WITH='MATERIAL'
            ERROR='>> ''MATERIAL'' feature not implemented'
            GOTO 9999
          ELSE
            ERROR='>> Not implemented'
            GOTO 9999
          ENDIF
          IF(CBBREV(CO,'SCALE_FACTOR',2,PART2,NTCO,N3CO))THEN
            SCALE=RFROMC(CO(N3CO+1))
          ELSE
            SCALE=1.0d0
          ENDIF

        ELSE
          ERROR='Code needs updating'
          GOTO 9999
        ENDIF !type



        IF(TYPE(1:9).EQ.'OPERATION')THEN

          FLUXBC=2 ! flux boundary conditions
          IF(WITH(1:8).EQ.'CONSTANT')THEN
            nr=nrc1
            nh1=NH_LOC(nhc1,nxc1)
            IF(OPERATION(1:10).EQ.'SUBSTITUTE')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=SCALE*CONST
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:3).EQ.'ADD')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)+SCALE*CONST
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:8).EQ.'SUBTRACT')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)-SCALE*CONST
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:8).EQ.'MULTIPLY')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)*SCALE*CONST
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:8).EQ.'DIVIDE')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)/(SCALE*CONST)
                ENDIF
              ENDDO
            ENDIF
          ELSEIF(WITH(1:5).EQ.'FIELD')THEN
            nr=nrc1
            nh1=NH_LOC(nhc1,nxc1)
            nj1=NJ_LOC(NJL_FIEL,njf2,nr)
            IF(OPERATION(1:10).EQ.'SUBSTITUTE')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=SCALE*XP(1,1,nj1,np)
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:3).EQ.'ADD')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)+
     '              SCALE*XP(1,1,nj1,np)
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:8).EQ.'SUBTRACT')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)-
     '              SCALE*XP(1,1,nj1,np)
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:8).EQ.'MULTIPLY')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)*
     '              SCALE*XP(1,1,nj1,np)
                ENDIF
              ENDDO
            ELSEIF(OPERATION(1:8).EQ.'DIVIDE')THEN
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny_fluxbc=NYNP(1,1,nh1,np,0,FLUXBC,nr)
                IF(FIX(ny_fluxbc,1,nx).AND.
     '            (.NOT.FIX(ny_fluxbc,2,nx)))THEN
                  YP(ny_fluxbc,1,nx)=YP(ny_fluxbc,1,nx)/
     '              (SCALE*XP(1,1,nj1,np))
                ENDIF
              ENDDO
            ENDIF
          ELSE
            ERROR='>> Not implemented'
            GOTO 9999
          ENDIF
        ENDIF
      ENDIF


      CALL EXITS('UPFLUX')
      RETURN
 9999 CALL ERRORS('UPFLUX',ERROR)
      CALL EXITS('UPFLUX')
      RETURN 1
      END


