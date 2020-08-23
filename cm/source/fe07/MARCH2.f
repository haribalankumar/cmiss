      SUBROUTINE MARCH2(NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,
     '  NYNE,NYNP,NYNR,XP,YP,ZA,ZP,ERROR,*)

C#### Subroutine: MARCH2
C###  Description:
C###    MARCH2 solves PDE of the form    A Uxx + B Uxt + C Utt = F
C###    using the method of characteristics.

C**** Calls Subroutine STARTCON to set initial conditions.
C**** Calls Subroutine CHARAC to calculate solution at point R
C**** from known values at points P & Q on previous solution front.
C**** At each soln point nine params are stored in array SURF(IT,IX,IS)
C**** These params are in order (IS=1-9) A,B,C,F,P,Q,X,T,U (P=Ux, Q=Ut)
C**** SURF(IT,nonode,7) is x-position of node nonode at iteration IT
C**** SURF(IT,nonode,8) is t-position of node nonode at iteration IT
C**** SURF(IT,nonode,9) is value of node nonode at iteration IT
C**** YP(ny,1) is solution vector at time T+DT
C**** YP(ny,3) is incremental boundary conditions
C****       4  "  solution vector at time T
C****       5  "  reaction vector at time T+DT

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'load00.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVHP(NHM,NPM,NCM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,it,IXSTEP,N,nonode,no_nynr
      REAL*8 APP,AQQ,ARR,BPP,BQQ,BRR,CPP,CQQ,CRR,FPP,FQQ,FRR,PPP,PQQ,
     '  PRR,QPP,QQQ,QRR,TPP,TQQ,TRR,UPP,UQQ,URR,XPP,XQQ,XRR,
     '  SURF(100,11,9)

      CALL ENTERS('MARCH2',*9999)

C**** Set initial conditions
      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '  nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
      CALL STARTCON(NPNODE,SURF,XP,ZP,ERROR,*9999)

      CALL STRING_TRIM(FILE02,IBEG,IEND)
C!!! same unit as IOOUT for `set output'!
      CALL OPENF(9,'DISK',FILE02(IBEG:IEND)//'.history','NEW',
     '  'SEQUEN','FORMATTED',132,ERROR,*9999)
C!!!!!! T not defined
C      WRITE(9,'('' YP(ny,1,nx) at t='',D11.4,'' :'')') T
      WRITE(9,'(10D13.5)') (YP(NYNR(no_nynr,0,1),1),
     '  no_nynr=1,NYNR(0,0,1))
C!!!!!! T not defined
C      WRITE(9,'('' YP(ny,2,nx) at t='',D11.4,'' :'')') T
      WRITE(9,'(10D13.5)') (YP(NYNR(no_nynr,0,1),2),
     '  no_nynr=1,NYNR(0,0,1))

      DO nonode=1,NPNODE(0,nr)
C****   Write (*,*)'A=', SURF(1,nonode,1)
C****   Write (*,*)'B=', SURF(1,nonode,2)
C****   Write (*,*)'C=', SURF(1,nonode,3)
C****   Write (*,*)'F=', SURF(1,nonode,4)
        WRITE(OP_STRING,*)'X=', SURF(1,nonode,7)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)'T=', SURF(1,nonode,8)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C****   Write (*,*)'P=', SURF(1,nonode,5)
C****   Write (*,*)'Q=', SURF(1,nonode,6)
        WRITE(OP_STRING,*)'U=', SURF(1,nonode,9)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDDO

C**** Outer Loop. Each new solution front designated by new value of IT.
      DO it=1,NTLOAD
        N=NPNODE(0,1)-it

C**** Step along current solution front nonode = 1 to N.
C**** Read values for current points P,Q from array SURF.

        DO nonode =1,N
          APP=SURF(it,nonode,1)
          BPP=SURF(it,nonode,2)
          CPP=SURF(it,nonode,3)
          FPP=SURF(it,nonode,4)
          PPP=SURF(it,nonode,5)
          QPP=SURF(it,nonode,6)
          XPP=SURF(it,nonode,7)
          TPP=SURF(it,nonode,8)
          UPP=SURF(it,nonode,9)

          AQQ=SURF(it,nonode+1,1)
          BQQ=SURF(it,nonode+1,2)
          CQQ=SURF(it,nonode+1,3)
          FQQ=SURF(it,nonode+1,4)
          PQQ=SURF(it,nonode+1,5)
          QQQ=SURF(it,nonode+1,6)
          XQQ=SURF(it,nonode+1,7)
          TQQ=SURF(it,nonode+1,8)
          UQQ=SURF(it,nonode+1,9)


C****     Calculate solution at new point R on characteristics
C****     through current P & Q
          CALL CHARAC(APP,AQQ,ARR,BPP,BQQ,BRR,CPP,CQQ,CRR,FPP,FQQ,FRR,
     ,                XPP,XQQ,XRR,TPP,TQQ,TRR,
     ,                PPP,PQQ,PRR,QPP,QQQ,QRR,UPP,UQQ,URR)

C****     Store solution at R and reset P,Q variables for next space step.
C****     IXSTEP = 0 or 1 alternately?

          IXSTEP=0

          SURF(it+1,nonode+IXSTEP,1)=ARR
          SURF(it+1,nonode+IXSTEP,2)=BRR
          SURF(it+1,nonode+IXSTEP,3)=CRR
          SURF(it+1,nonode+IXSTEP,4)=FRR
          SURF(it+1,nonode+IXSTEP,5)=PRR
          SURF(it+1,nonode+IXSTEP,6)=QRR
          SURF(it+1,nonode+IXSTEP,7)=XRR
          SURF(it+1,nonode+IXSTEP,8)=TRR
          SURF(it+1,nonode+IXSTEP,9)=URR
        ENDDO

        DO nonode=1,N
          WRITE (*,*) nonode
C****     Write (*,*)'A=', SURF(it+1,nonode,1)
C****     Write (*,*)'B=', SURF(it+1,nonode,2)
C****     Write (*,*)'C=', SURF(it+1,nonode,3)
C****     Write (*,*)'F=', SURF(it+1,nonode,4)
          WRITE(OP_STRING,*)'X=', SURF(it+1,nonode,7)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)'T=', SURF(it+1,nonode,8)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C****     Write (*,*)'P=', SURF(it+1,nonode,5)
C****     Write (*,*)'Q=', SURF(it+1,nonode,6)
          WRITE(OP_STRING,*)'U=', SURF(it+1,nonode,9)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ZP(1,1,1,nonode,1)=SURF(it+1,nonode,9)
        ENDDO
C!!!!!! T not defined
C       WRITE(9,'('' YP(ny,1,nx) at t='',D11.4,'' :'')') T+DT
        WRITE(9,'(10D13.5)') (YP(NYNR(no_nynr,0,1),1),
     '    no_nynr=1,NYNR(0,0,1))
C!!!!!! T not defined
C       WRITE(9,'('' YP(ny,2,nx) at t='',D11.4,'' :'')') T+DT
        WRITE(9,'(10D13.5)') (YP(NYNR(no_nynr,0,1),2),
     '    no_nynr=1,NYNR(0,0,1))
      ENDDO

      CALL EXITS('MARCH2')
      RETURN
 9999 CALL ERRORS('MARCH2',ERROR)
      CALL EXITS('MARCH2')
      RETURN 1
      END



