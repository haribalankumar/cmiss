      SUBROUTINE EVINTE(NLLINE,NLLIST,NPL,NRLIST,
     &          XP,STRING,ERROR,*)

C#### Subroutine: EVINVT
C###  Description:
C###    EVINVT evaluates the integral of a given field for a
C###      list of point, line, face, or volume.
C**** Written by Duane Malcolm, 29 January 2003

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),
     &  NPL(5,0:3,NLM),NRLIST(0:NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER STRING*(MXCH),ERROR*(*)

!     Local Variables
      INTEGER ERR,IBEG,IEND,IFROMC,line,N3CO,njf,nl,nr
      REAL*8 INTEGRAL,SUM,CHXE(NSM,4)
      LOGICAL ALL_REGIONS,CBBREV,VERBOSE,TERMINAL,VARIABLE
      CHARACTER VARNAME*255

      CALL ENTERS('EVINTE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate integral
C###  Parameter:        <field #[1]>
C###    Specify the field to integrate.
C###  Parameter:        <points|lines|faces|volumes #s>
C###    Specify the type of elements to integrate over.
C###  Parameter:        <verbose>
C###    To output element integrals.
C###  Parameter:        <terminal>
C###    Outputs the total integral to the terminal.
C###  Parameter:        <variable VARIABLE_NAME>
C###    Outputs total integral to a interpreter variable.
C###  Parameter:        <region #[1]>
C###    Specify the region.
C###  Description:
C###    Evaluates the integral of a 1D field and saves the result in
C###    another field. A reference node and value can be specified.
C###    The is no 2D and 3D functionality yet.

C DMAL 29 OCT 2003
C This is a new command that will sum the integrals of 
C point, lines, surfaces, and volumes.
C This has been implemented for lines but the structure
C is there to expand it to cope with the other element types.

C This command should be used with the grouping commands to allow the user
C to specify the region to integrate. This allows flexibility and reduces
C the replication of code.
C For example,
C >>fem group lines between nodes 1,10 in xi1-direction region 1 as BOUNDARY
C >>fem evaluate integral field 3 lines BOUNDARY region 1

C The output should be extended to pass the result to a intepreter variable.
C See "fem update variable" command for help.

C Currently, this command only copes with cubic-Hermite lines. The gauss points
C and weights and the basis functions are calculated locally instead of using 
C the cmiss basis functions. This was so the user doesn't have to setup every 
C permutation of basis that exists in the domamin. For instance, if a 3d
C domain is defined, the user is requires to setup a 2d face basis for
C every permutation of basis in the 3d domain. Either cmiss should
C automaically setup every permutation for user defined domains or there should
C be simple functions to create the basis functions on the fly (this might 
C be supplied by a library) or both.


        OP_STRING(1)=STRING(1:IEND)//' <field #[1]>'
        OP_STRING(2)=BLANK(1:15)//'<[points|lines|faces|volumes] #s>'
        OP_STRING(3)=BLANK(1:15)//'<verbose>'
        OP_STRING(4)=BLANK(1:15)//'<terminal>'
        OP_STRING(5)=BLANK(1:15)//'<variable VARIABEL_NAME>'
        OP_STRING(6)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
      
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVINTE',ERROR,*9999)
      ELSE
        
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF (NRLIST(0).EQ.1) THEN
          nr=NRLIST(1)
        ELSE
          CALL ASSERT(.FALSE.,'>> Mutliple regions not implemented',
     &     ERROR,*9999)
        ENDIF
        
        
C Parse field to integrate
        IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
          njf=IFROMC(CO(N3CO+1))
        ELSE
          njf=1
        ENDIF

C Parse output
        IF(CBBREV(CO,'VERBOSE',3,noco+1,NTCO,N3CO)) THEN
          VERBOSE=.TRUE.
        ELSE
          VERBOSE=.FALSE.
        ENDIF
        
        IF(CBBREV(CO,'TERMINAL',3,noco+1,NTCO,N3CO)) THEN
          TERMINAL=.TRUE.
        ELSE
          TERMINAL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'VARIABLE',3,noco+1,NTCO,N3CO)) THEN
          VARIABLE=.TRUE.
          VARNAME=CO(N3CO+1)
        ELSE
          VARIABLE=.FALSE.
        ENDIF

C Parse node, line, face or volume numbers
        IF(CBBREV(CO,'NODES',2,noco+1,NTCO,N3CO))THEN
          CALL ASSERT(.FALSE.,'>> Nodes not implemented',
     &     ERROR,*9999)
        ELSEIF(CBBREV(CO,'LINES',2,noco+1,NTCO,N3CO))THEN
          CALL PARSE_LINES(NLLINE,NLLIST,noco-1,NRLIST,NTCO,CO,
     &      ERROR,*9999)
          
          SUM=0.0D0
          DO line=1,NLLIST(0)
            nl=NLLIST(line)
            
            ! Convert basis to cubic-Hermite
            CALL XPXE_CHLINE(CHXE,njf,nl,NPL,nr,XP,ERROR,*9999)
            ! Calculate line integral
            CALL INTEGRATE_LINE(CHXE,INTEGRAL,ERROR,*9999)
            
            SUM=SUM+INTEGRAL
            
            IF(VERBOSE) WRITE(IOOP,'('' Integral ('',I5,'
     &        //''') = '',E12.5)') nl,INTEGRAL
            
          ENDDO
          
          IF(TERMINAL)THEN
            WRITE(IOOP,'('' Total Integral='',E12.5)') SUM
          ENDIF

          IF(VARIABLE)THEN
            CALL STRING_TRIM(VARNAME,IBEG,IEND)
            CALL SET_USER_DOUBLE(VARNAME(IBEG:IEND),SUM,ERR)
          ENDIF
          
        ELSEIF(CBBREV(CO,'FACES',2,noco+1,NTCO,N3CO))THEN
          CALL ASSERT(.FALSE.,'>> Faces not implemented',
     &      ERROR,*9999)
        ELSEIF(CBBREV(CO,'VOLUMES',2,noco+1,NTCO,N3CO))THEN
          CALL ASSERT(.FALSE.,'>> Volumes not implemented',
     &      ERROR,*9999)
        ENDIF
          
C Parse output (Terminal)

      ENDIF

      CALL EXITS('EVINTE')
      RETURN
 9999 CALL ERRORS('EVINTE',ERROR)
      CALL EXITS('EVINTE')
      RETURN 1
      END


