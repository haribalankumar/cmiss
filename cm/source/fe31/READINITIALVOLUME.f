      SUBROUTINE READINITIALVOLUME(MECHANICS_FILETYPE,NBJ,NDP,NEELEM,
     &  NENP,NPNE,NXI,BBM,CE,sumvolume,undef,XP,ZD,FIRST_READ,
     &  ERROR,*)

C#### Subroutine: ReadInitialVolume
C###  Description:
C###    ReadInitialVolume calculates the initial (FRC) volume of each
C###    acinar unit. Requires that ratio of deformed to undeformed volume
C###    has been read in prior to fem evaluate flow command. Volume ratio
C###    must be in field 5 of ipdata file (generally, fields in this file
C###    are: x,y,z,Pe,vol_ratio,elasticity

C***  BBM(1,ne) = volume of acinar unit subtended by terminal bronchiole ne
C***  BBM(2,ne) = concentration in acinar unit subtended by terminal bronchiole ne
C***  ZD(nj_Vratio,nd) = ratio of deformed to undeformed volume (read in from ipdata file, where nj_Vratio is typically 5)
C***        OR...
C***  XP(1,1,nj_Vratio,np) = ratio of deformed to undeformed volume (read in from ipfiel file)

      
      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn' 
!     Parameter List
      INTEGER MECHANICS_FILETYPE,NBJ(NJM,NEM),NDP(NDM),NEELEM(0:NE_R_M),
     &  NENP(NPM,0:NEPM),NPNE(NNM,NBFM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),sumvolume,undef,
     &  XP(NKM,NVM,NJM,NPM),
     &  ZD(NJM,NDM)
      CHARACTER ERROR*(*)
      LOGICAl FIRST_READ
!     Local variables
      INTEGER nb,nd,ne,nj,noelem,np,Nterms
      REAL*8 mean_ratio,Vratio 
      CHARACTER STRING*50

      sumvolume=0.d0
      mean_ratio=0.d0
      Nterms=0

C    AJS Add option to read mechanics data from field file
      IF(MECHANICS_FILETYPE.EQ.1)THEN!'IPDATA')THEN !read mechanics results from ipdata file
        DO nd=1,NDT
          np=nd !have to assume sequential numbering else not working!!!
c        np=NDP(nd)
c        ne=NENP(np,1) !the first element that the node is in
          ne=np-1
          IF(ne.GT.0)THEN
            IF(NXI(1,0,ne).EQ.0)THEN !terminal
              Vratio=ZD(nj_Vratio,nd)
              IF(Vratio.LE.0.0d0)THEN !check volume ratio is >0
                WRITE(STRING,'('' Volume ratio <=0 for node '',I5)') np
                CALL ASSERT(.FALSE.,STRING,ERROR,*9999)
              ENDIF
              BBM(1,ne)=Vratio*undef !acinar volume mm3
              mean_ratio=mean_ratio+Vratio
              IF(FIRST_READ) CE(nm_vinit,ne)=BBM(1,ne)
              Nterms=Nterms+1
              BBM(2,ne)=0.d0 !acinar concentration
              sumvolume=sumvolume+BBM(1,ne) !sum total volume of alveolar air (respiratory volume)
            ENDIF
          ENDIF !ne
        ENDDO

      ELSEIF(MECHANICS_FILETYPE.EQ.2)THEN!'IPFIEL')THEN !read mechanics results from ipfiel file
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          IF(NXI(1,0,ne).EQ.0)THEN !terminal
            nb=NBJ(1,ne)
            np=NPNE(2,nb,ne) !end node
            Vratio=XP(1,1,nj_Vratio,np)
            IF(Vratio.LE.0.0d0)THEN !check volume ratio is >0
              WRITE(STRING,'('' Volume ratio <=0 for node '',I5)') np
              CALL ASSERT(.FALSE.,STRING,ERROR,*9999)
            ENDIF
            BBM(1,ne)=Vratio*undef !acinar volume mm3
            mean_ratio=mean_ratio+Vratio
            IF(FIRST_READ) CE(nm_vinit,ne)=BBM(1,ne)
            Nterms=Nterms+1
            BBM(2,ne)=0.d0 !acinar concentration
            sumvolume=sumvolume+BBM(1,ne) !sum total volume of alveolar air (respiratory volume)
          ENDIF
        ENDDO !noelem
      ENDIF !MECHANICS_FILETYPE

      RETURN
 9999 CALL ERRORS('ReadInitialVolume',ERROR)
      RETURN
      END

