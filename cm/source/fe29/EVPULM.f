      SUBROUTINE EVPULM(NBJ,NEELEM,NELIST,NPNE,NRLIST,NVJE,NXI,NXLIST,
     &  NYNP,CE,CP,XP,YP,STRING,ERROR,*)

C#### Subroutine: EVPULM
C###  Description:
C###    EVPULM calculates and lists parameters from pulmonary 
C###    simulations.  e.g. the normalised slopes over multiple gas
C###    mixing breaths.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn' 
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 CE(NMM,NEM,NXM),CP(NMM,NPM,NXM),XP(NKM,NVM,NJM,NPM),
     &  YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,nb,NBREATH,ne,ne_CC,ne_front,ne0,nh,
     &  nobreath,noelem,norec,np,np1,np2,nx,ny,ny1,ny2
      REAL*8 acinus_Pe,CONC,CONCFRONT,CUM_EV,distance_acinus,
     &  distance_acinus0,distance_trachea,EV,EV_END,EV_PREV,EV_START,
     &  gradient,gradient0,LconvLdiff,length,locationDF,
     &  locationCC,norm_slope,Pe,pxi,radius,Sacin,
     &  Scond,slope,Sn0,Sn1,Sn2,SUMX,SUMXSQ,SUMXY,SUMY,SUM_N2,TO,TO0,
     &  TO1,TO2,VTIDAL
      CHARACTER FILE*(MXCH)
      LOGICAL ALL_REGIONS,CBBREV,FOUND0,FOUND1,FOUND2,FRONT,NEXT,OPFILE,
     &  SN,START,TRANSPORT
!     Functions
      INTEGER IFROMC
      REAL*8 LENGTH_1D,RADIUS_1D,RFROMC

      CALL ENTERS('EVPULM',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate pulmonary<;FILENAME>
C###  Parameter:        <sn>
C###    Specifies that the normalised slopes will be evaluated
C###  Parameter:        <for (#breaths/1)[1]>
C###    Slopes will be evaluated for a default of 1 breath, or for the
C###    specified # of breaths
C###  Parameter:        <ev_start [0.75]>
C###    Specifies the start expired volume for Sn calculation
C###  Parameter:        <ev_end [0.95]>
C###    Specifies the end expired volume for Sn calculation
C###  Parameter:        <vdf>
C###    Specifies that the Fowler dead space will be evaluated
C###  Parameter:        <vdb>
C###    Specifies that the Bohr dead space will be evaluated
C###  Description:
C###    EVPULM evaluates normalised slopes from .history file output

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> '
        OP_STRING(2)=BLANK(1:15)//'<sn>'
        OP_STRING(3)=BLANK(1:15)//'<for (#breaths/1)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<ev_start [0.75]>'
        OP_STRING(5)=BLANK(1:15)//'<ev_end [0.95]>'
        OP_STRING(6)=BLANK(1:15)//'<vt (tidal vol/vt)[vt]>'
        OP_STRING(7)=BLANK(1:15)//'<vdb>'
        OP_STRING(8)=BLANK(1:15)//'<vdf>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVPULM',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          CALL OPENF(IOFILE4,'DISK',FILE(IBEG:IEND)//'.oppulm','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL ASSERT(OPFILE,'>>Files were not opened',ERROR,*9999)
        SN=.FALSE.
        FRONT=.FALSE.
        TRANSPORT=.FALSE.
        VTIDAL=1.d0
        IF(CBBREV(CO,'SN',2,noco+1,NTCO,N3CO)) THEN
          CALL OPENF(IOFILE5,'DISK',FILE(IBEG:IEND)//'.history','OLD',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)

          SN=.TRUE.
          IF(CBBREV(CO,'FOR',3,noco+1,NTCO,N3CO))THEN
            NBREATH=IFROMC(CO(N3CO+1))
          ELSE
            NBREATH=1
          ENDIF
          IF(CBBREV(CO,'EV_START',5,noco+1,NTCO,N3CO))THEN
            EV_START=RFROMC(CO(N3CO+1))
          ELSE
            EV_START=0.75d0
          ENDIF
          IF(CBBREV(CO,'EV_END',5,noco+1,NTCO,N3CO))THEN
            EV_END=RFROMC(CO(N3CO+1))
          ELSE
            EV_END=0.95d0
          ENDIF
          IF(CBBREV(CO,'VT',2,noco+1,NTCO,N3CO)) THEN
            VTIDAL=RFROMC(CO(N3CO+1))
          ENDIF

        ELSE IF(CBBREV(CO,'FRONT',3,noco+1,NTCO,N3CO)) THEN
          FRONT=.TRUE.
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     &         *9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nx=NXLIST(1)
          IF(CBBREV(CO,'ELEMENTS',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,ERROR,
     &        *9999)
          ENDIF
        ELSE IF(CBBREV(CO,'TRANSPORT',3,noco+1,NTCO,N3CO)) THEN
          TRANSPORT=.TRUE.
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     &         *9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nx=NXLIST(1)
          IF(CBBREV(CO,'ELEMENTS',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,ERROR,
     &        *9999)
          ENDIF
        ELSE IF(CBBREV(CO,'THERMAL',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     &      *9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nx=NXLIST(1)
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,ERROR,
     &      *9999)
        ENDIF

        IF(SN)THEN
          FOUND0=.FALSE.
          FOUND1=.FALSE.
          FOUND2=.FALSE.
          CUM_EV=0.d0
          EV_PREV=0.d0
          DO nobreath=1,NBREATH
            SUMX=0.d0
            SUMXY=0.d0
            SUMXSQ=0.d0
            SUMY=0.d0
            SUM_N2=0.d0
            norec=0
            START=.FALSE.
            DO WHILE(.NOT.START) !find first calculation point
              READ(IOFILE5,'(2(F15.5))') EV,CONC
              IF(EV.GE.(EV_START*VTIDAL+CUM_EV))THEN
                START=.TRUE.
              ELSE
                SUM_N2=SUM_N2+(EV-EV_PREV)*CONC
                EV_PREV=EV
              ENDIF !EV
            ENDDO !START
            NEXT=.FALSE.
            DO WHILE(.NOT.NEXT)
              SUMX=SUMX+EV
              SUMY=SUMY+CONC
              SUMXSQ=SUMXSQ+EV*EV
              SUMXY=SUMXY+EV*CONC
              SUM_N2=SUM_N2+(EV-EV_PREV)*CONC
              norec=norec+1
              EV_PREV=EV
              READ(IOFILE5,'(2(F15.5))') EV,CONC
              IF(EV.GT.(EV_END*VTIDAL+CUM_EV)) NEXT=.TRUE.
            ENDDO !NEXT
            DO WHILE(EV.LT.VTIDAL+CUM_EV)
              READ(IOFILE5,'(2(F15.5))') EV,CONC
              SUM_N2=SUM_N2+(EV-EV_PREV)*CONC
              EV_PREV=EV
            ENDDO !EV
            slope=-(norec*SUMXY-SUMX*SUMY)/(norec*SUMXSQ-SUMX*SUMX)
            SUM_N2=1.d0-SUM_N2/VTIDAL !g/mm^3
            norm_slope=slope/SUM_N2

            TO=nobreath*VTIDAL/(CURRENT_VOLUME(nx)/1.d6)
            IF(.NOT.FOUND0)THEN
              Sn0=norm_slope
              TO0=TO
              FOUND0=.TRUE.
            ELSE IF(.NOT.FOUND1.AND.TO.GE.1.5d0)THEN
              TO1=TO
              Sn1=norm_slope
              FOUND1=.TRUE.
            ENDIF
            IF(FOUND1.AND.TO.GE.6.d0)THEN
              IF(.NOT.FOUND2)THEN
                TO2=TO
                Sn2=norm_slope
                FOUND2=.TRUE.
              ENDIF
            ENDIF

            WRITE(IOFILE4,'(I4,4(D13.5))') nobreath,TO,norm_slope,
     '        slope,SUM_N2
            WRITE(OP_STRING,'('' breath'',I3,''  TO='',F7.2,
     '        ''  Sn='',F10.5,
     '        ''  slope='',D10.3,''  N2='',F10.5)')
     '        nobreath,TO,norm_slope,slope,sum_n2
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CUM_EV=CUM_EV+VTIDAL
          ENDDO !nobreath

          IF(FOUND1.AND.FOUND2)THEN
            Scond=(Sn2-Sn1)/(TO2-TO1)
            Sacin=Sn0-Scond*TO0
            WRITE(OP_STRING,'('' Sacin = '',F10.4,
     '        '' L.-1, Scond = '',F10.4)') Sacin,Scond
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Calculated using TO = '',F10.4,
     '        '' and '',F10.4)') TO1,TO2
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(IOFILE4,'('' Sacin = '',F10.4,
     '        '' L.-1, Scond = '',F10.4)') Sacin,Scond
          ENDIF
          IF(TO2.LT.6.d0)THEN
            WRITE(OP_STRING,
     &        '('' Warning: TO did not reach 6.0, run longer'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE IF(FRONT.OR.TRANSPORT)THEN !write out the diffusion front information
          CALL ASSERT(NMM.GE.6,'>> Set NMM >= 6',ERROR,*9999)

c          IF(FRONT)THEN
c            WRITE(IOFILE4,'(''  TBelem''
c     &        ''   TrachDist    z height     ConcTB   ''
c     &        '' %FlowTB      PeTB      DFelem    distAcn  ''
c     &        '' PeDF   acinusPE LconvLdiff ne_CC locCC  PECC'')')
c          ELSE IF (TRANSPORT)THEN
c            WRITE(IOFILE4,'(''  Elem ''
c     &        ''   TrachDist    TBDist     ConcGrad    %Flow  ''
c     &        '' ConcMean      Pe '')')
c          ENDIF
          DO noelem=1,NEELEM(0,1)
            ne=NEELEM(noelem,1)
            nb=NBJ(1,ne)
            np=NPNE(2,nb,ne)
            ny1=NYNP(1,1,1,NPNE(1,nb,ne))
            ny2=NYNP(1,1,1,NPNE(2,nb,ne))
C... Calculate the gradient
            gradient = (YP(ny2,1,nx)-YP(ny1,1,nx))/
     &        LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
            CE(5,ne,nx)=gradient
C... Calculate the distance from the entry
            distance_trachea=0.d0
            ne0=NXI(-1,1,ne)
            DO WHILE(ne0.NE.0)
                length=LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP) !length of this element
                distance_trachea=distance_trachea+length
                ne0=NXI(-1,1,ne0)
            ENDDO              !while not at stem
            CE(6,ne,nx)=distance_trachea
          ENDDO

          DO noelem=1,NELIST(0)
             ne=NELIST(noelem)
             nb=NBJ(1,ne) !geometry basis function
             np=NPNE(2,nb,ne)
             ny1=NYNP(1,1,1,NPNE(1,nb,ne))
             ny2=NYNP(1,1,1,NPNE(2,nb,ne))
C... Calculate Peclet number
             radius=RADIUS_1D(NBJ,ne,NPNE,NVJE,XP)
     &            *DSQRT(XP(1,1,nj_alveoli,np))
             Pe=2.d0*radius*XP(1,1,nj_flow,np)
     &            *INLET_FLOW(1)/(PI*radius**2)/PULMAT(1)
             XP(1,1,9,np)=Pe
             IF(TRANSPORT)THEN
                WRITE(IOFILE4,'(I6,3(F12.4),F16.8,2(F12.4))') ne,
     &               CE(6,ne,nx),0.d0,CE(5,ne,nx),XP(1,1,nj_flow,np)
     &               *100.d0,0.5d0*(YP(ny1,1,nx)
     &               +YP(ny2,1,nx)),Pe
             ENDIF
C search backwards to write out all results
             ne0=NXI(-1,1,ne)
             DO WHILE(ne0.NE.0)
                ny=NYNP(1,1,1,NPNE(1,nb,ne0))
                np=NPNE(2,nb,ne0)
                ny1=NYNP(1,1,1,NPNE(1,nb,ne0))
                ny2=NYNP(1,1,1,NPNE(2,nb,ne0))
C... Calculate Peclet number
                radius=RADIUS_1D(NBJ,ne0,NPNE,NVJE,XP)
     &               *DSQRT(XP(1,1,nj_alveoli,np))
                Pe=2.d0*radius*XP(1,1,nj_flow,np)
     &               *INLET_FLOW(1)/(PI*radius**2)/PULMAT(1)
                XP(1,1,9,np)=Pe

                IF(TRANSPORT)THEN
                   WRITE(IOFILE4,'(I6,3(F12.4),F16.8,2(F12.4))') 
     &               ne0,CE(6,ne0,nx),CE(6,ne0,nx)-CE(6,ne,nx),
     &               CE(5,ne0,nx),XP(1,1,nj_flow,np)*100.d0,
     &               0.5d0*(YP(ny1,1,nx)+YP(ny2,1,nx)),Pe
                ENDIF               
                ne0=NXI(-1,1,ne0)
             ENDDO 
! distance_trachea == distance to end of terminal bronchiole
C search forwards to write out all results
             ne0=NXI(1,1,ne)
             DO WHILE(ne0.NE.0)
                ny=NYNP(1,1,1,NPNE(1,nb,ne0))
                ny1=NYNP(1,1,1,NPNE(1,nb,ne0))
                ny2=NYNP(1,1,1,NPNE(2,nb,ne0))
                np=NPNE(2,nb,ne0)
C... Calculate Peclet number
                radius=RADIUS_1D(NBJ,ne0,NPNE,NVJE,XP)
     &               *DSQRT(XP(1,1,nj_alveoli,np))
                Pe=2.d0*radius*XP(1,1,nj_flow,np)
     &               *INLET_FLOW(1)/(PI*radius**2)/PULMAT(1)
                XP(1,1,9,np)=Pe

                IF(TRANSPORT)THEN
                  WRITE(IOFILE4,'(I6,3(F12.4),F16.8,2(F12.4))') ne0,
     &               CE(6,ne0,nx),CE(6,ne0,nx)-CE(6,ne,nx),CE(5,ne0,nx),
     &               XP(1,1,nj_flow,np)*100.d0,0.5d0*(YP(ny1,1,nx)+
     &               YP(ny2,1,nx)),Pe   
                ENDIF           
                ne0=NXI(1,1,ne0)
             ENDDO 
          ENDDO

C.........Note that this element list should only be the terminal
C.........bronchioles!!!
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            nb=NBJ(1,ne)
            np=NPNE(2,nb,ne) !the end node, will use version 1
            IF(NXI(1,0,ne).GT.0)THEN !not appending a LPM so write out
C.............Calculate the location of the diffusion front
C....NOTE THAT THIS IS FOR A SYMMETRIC ACINUS ONLY
   ! find the end of the acinus
              ne0=ne !set equal to the terminal bronchiole
              distance_acinus=0.d0
              DO WHILE(NXI(1,0,ne0).NE.0)
                 ne0=NXI(1,1,ne0)
                 distance_acinus = distance_acinus+
     &             LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP)
              ENDDO
              distance_acinus0=distance_acinus
   ! ne0 == the element number of the terminal alveolated branch
              gradient = CE(5,ne0,nx)
              FOUND0=.FALSE.
              DO WHILE(ne0.NE.ne) !i.e. while not at the TB
                 distance_acinus=distance_acinus-
     &             LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP)
                 ny1=NYNP(1,1,1,NPNE(1,nb,ne0))
                 ny2=NYNP(1,1,1,NPNE(2,nb,ne0))

c                 IF(DABS(CE(5,ne0,nx)).GT.DABS(gradient))THEN
c                    gradient = CE(5,ne0,nx)
c                    ne_front = ne0
c                    locationDF = distance_acinus + 0.5d0*
c     &                LENGTH_1D(NBJ,ne_front,NPNE,NVJE,XP)
c                 ENDIF
C... Determine the DF as found when the concentration is less than 99% inspired
                 IF(.NOT.FOUND0)THEN
                    IF(YP(ny1,1,nx).GE.0.99d0)THEN
                       locationCC = distance_acinus + 0.5d0*
     &                      LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP)
                       ne_CC = ne0
                       FOUND0=.TRUE.
                    ENDIF
                 ENDIF
                 ne0=NXI(-1,1,ne0) !move one element proximal
              ENDDO

              distance_acinus = 0.d0
              IF(.NOT.FOUND0)THEN !check whether at acinus entry
                 ny1=NYNP(1,1,1,NPNE(1,nb,ne)) !start of TB
                 IF(YP(ny1,1,nx).GE.0.99d0)THEN
                    locationCC = 0.d0
                    ne_CC = ne
                    FOUND0=.TRUE.
                 ENDIF
              ENDIF

   ! also check that the front is not proximal to the acinus
              DO WHILE(.NOT.FOUND0)
                 distance_acinus=distance_acinus-
     &             LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP)
c                 IF(DABS(CE(5,ne0,nx)).GT.DABS(gradient))THEN
c                    gradient = CE(5,ne0,nx)
c                    ne_front = ne0
c                    locationDF = distance_acinus + 0.5d0*
c     &                   LENGTH_1D(NBJ,ne_front,NPNE,NVJE,XP)
c                 ENDIF
                 IF(.NOT.FOUND0)THEN
                    ny1=NYNP(1,1,1,NPNE(1,nb,ne0))
                    IF(YP(ny1,1,nx).GE.0.99d0)THEN
                       locationCC = distance_acinus + 0.5d0*
     &                      LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP)
                       ne_CC = ne0
                       FOUND0=.TRUE.
                    ENDIF
                 ENDIF
                 ne0=NXI(-1,1,ne0) !move one element proximal
              ENDDO

C              np2=NPNE(2,nb,ne_front)
              np2=NPNE(2,nb,ne_CC) ! calculating for front based on concentration
              
c              ne=NXI(1,1,ne) !set equal to the transitional bronchiole
              np1=NPNE(1,nb,ne)
              np=NPNE(2,nb,ne)

C              IF(ne_front.NE.0)THEN
              IF(ne_CC.NE.0)THEN

                ny=NYNP(1,1,1,np1) !concentration at acinar entrance

C... Calculate the acinus Peclet number at np2 (Sapoval)
                acinus_Pe = XP(1,1,nj_flow,np2)*INLET_FLOW(nx)/
     &            (PI*radius**2)*(distance_acinus0-locationCC)/
     &            PULMAT(1)
c                acinus_Pe = XP(1,1,nj_flow,np2)*INLET_FLOW(nx)/
c     &            (PI*radius**2)*(distance_acinus0-locationDF)/
c     &            PULMAT(1)

C... Calculate Lconv/Ldiff (Snitzman) for 4 second breath
                LconvLdiff = XP(1,1,nj_flow,np2)*INLET_FLOW(nx)/
     &            (PI*radius**2)*DSQRT(4.d0/(2.d0*PULMAT(1)))
                
                distance_trachea = CE(6,ne,nx)
                IF(FRONT)THEN
c                  WRITE(IOFILE4,
c     &              '(I7,5(F10.3),I7,F12.4,3(D12.4),I8,2(D12.4))')
c     &              ne,distance_trachea,XP(1,1,3,np1),YP(ny,1,nx),
c     &              XP(1,2,nj_flow,np)*INLET_FLOW(nx),XP(1,1,9,np),
c     &              ne_front,locationDF,XP(1,1,9,np2),
c     &              acinus_Pe,LconvLdiff,ne_CC,locationCC,
c     &              XP(1,1,9,NPNE(2,nb,ne_CC))
                  WRITE(IOFILE4,
     &              '(I7,5(F10.3),I7,F12.4,3(D12.4),I8,2(D12.4))')
     &              ne,distance_trachea,XP(1,1,3,np1),YP(ny,1,nx),
     &              XP(1,2,nj_flow,np)*INLET_FLOW(nx),XP(1,1,9,np),
     &              ne_CC,locationCC,
     &              XP(1,1,9,NPNE(2,nb,ne_CC)),acinus_Pe,LconvLdiff
                ENDIF
              ENDIF
            ENDIF !NXI not terminal
          ENDDO !noelem

        ELSE
          ne=NELIST(1)
          nb=NBJ(1,ne)
          np=NPNE(1,nb,ne) !the end node, will use version 1
          distance_trachea=0.d0
          nh=NH_LOC(1,nx)
          ny1=NYNP(1,1,nh,np) !for temperature in K
          ny2=NYNP(1,1,nh+1,np) !for humidity in g/mm^3
          WRITE(IOFILE4,'(I8,4(F12.4))')
     &      ne,distance_trachea,YP(ny1,1,nx)-273.15,CP(4,np,1)-273.15,
     &      YP(ny2,1,nx)
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            nb=NBJ(1,ne)
            np=NPNE(2,nb,ne) !the end node, will use version 1
            distance_trachea=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP) !length of this element
            ne0=NXI(-1,1,ne)
            DO WHILE(ne0.NE.0)
              length=LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP) !length of this element
              distance_trachea=distance_trachea+length
              ne0=NXI(-1,1,ne0)
            ENDDO !while not at stem
            nh=NH_LOC(1,nx)
            ny1=NYNP(1,1,nh,np) !for temperature in K
            ny2=NYNP(1,1,nh+1,np) !for humidity in g/mm^3
            WRITE(IOFILE4,'(I8,4(F12.4))')
     &        np,distance_trachea,YP(ny1,1,nx)-273.15,CP(4,np,1)-273.15,
     &        YP(ny2,1,nx)
          ENDDO !noelem
        ENDIF !SN
        
        IF(OPFILE) THEN
          CALL CLOSEF(IOFILE5,ERROR,*9999)
          CALL CLOSEF(IOFILE4,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVPULM')
      RETURN
 9999 CALL ERRORS('EVPULM',ERROR)
      CALL EXITS('EVPULM')
      RETURN 1
      END



