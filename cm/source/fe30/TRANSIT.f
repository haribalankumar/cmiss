      SUBROUTINE TRANSIT(NBJ,NEELEM,NENP,NENQ,NLL,NPNE,NQET,NQNE,
     '  NQS,nr,NYNQ,DL,XQ,YQ,ERROR,*)

C#### Subroutine: TRANSIT
C###  Description:
C###    TRANSIT calculates blood transit times through a blood vessel
C###    network.

C***  Created by KSB 17th July, 2003.
C***  Created for use in the pulmonary arterial and venous networks.

C***  PATH_END(path#) stores last element number of pathway
      
C***  PATHWAY(path#,2) stores sum of RBC transit time
C***  PATHWAY(path#,3) stores sum of the length of the pathway
C***  PATHWAY(path#,4) stores the number of elements in a pathway
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NLL(12,NEM),
     '  NPNE(NNM,NBFM,NEM),NQET(NQSCM),NQNE(NEQM,NQEM),NQS(NEQM),nr,
     '  NYNQ(NHM,NQM,0:NRCM)
      REAL*8 DL(3,NLM),XQ(NJM,NQM),YQ(NYQM,NIQM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,j,nb,ne,ne2,nl,noelem,no_nq,np2,np3,nq,ny_v,
     &  path,PATH_END(0:MAX_PATH),path_number,TOTAL_PATHS
      REAL*8 flow,LENGTH,PATHWAY(0:MAX_PATH,FACTORS),
     &  TIME,TIME_VEIN,TRANS_TIME(NEM),
     &  TRANS_TIME_VEIN(NEM),TOL
      LOGICAL CONTINU,NEW
      
      CALL ENTERS('TRANSIT',*9999)

      TOL=1.d-6 !tolerance
      DO i=0,MAX_PATH !initialise
        DO j=1,FACTORS !stores required factors for a given pathway
          PATHWAY(i,j)=0.d0
        ENDDO
        PATH_END(i)=0 !initialise integer array
      ENDDO
      DO noelem=1,NEM !initialise
        TRANS_TIME(noelem)=0.d0
        TRANS_TIME_VEIN(noelem)=0.d0
      ENDDO
      nb=NBJ(1,NENQ(1,NQ_START(nr)))
C... First calculate the transit time through each element      
C... Sum the transit time over an elem (ne) at each grid point (nq)
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        TIME=0.d0
        TIME_VEIN=0.d0 
        nl=NLL(1,ne)
        LENGTH=DL(3,nl)/NQET(NQS(ne))
        DO no_nq=1,NQET(NQS(ne))
          nq=NQNE(ne,no_nq)
          ny_v=NYNQ(3,nq,0) !ny # for velocity solution
          IF(DABS(YQ(ny_v,1)).GT.TOL) !If velocity is greater than 0
          !sums time for each grid point in an element
     &      TIME=TIME+LENGTH/YQ(ny_v,1) !time=length/velocity
          IF(VENOUS_NETWORK.EQ.'Y'.AND.N_VENOUS_GEOM.EQ.1) THEN
           !identical venous tree
            ny_v=NYNQ(6,nq,0) !ny # for venous velocity
            IF(DABS(YQ(ny_v,1)).GT.TOL)
     &        TIME_VEIN=TIME_VEIN+LENGTH/YQ(ny_v,1) !time for identical venous tree
          ENDIF
        ENDDO
        TRANS_TIME(ne)=TIME
        TRANS_TIME_VEIN(ne)=TIME_VEIN !identical venous tree transit time
      ENDDO !neelem
C... Calculate transit time for each flow pathway      
c      count=1
      TOTAL_PATHS=1       
c     PATHWAY(TOTAL_PATHS,1)=DBLE(NENQ(1,NQ_START(nr))) !start element #
      PATH_END(TOTAL_PATHS)=NENQ(1,NQ_START(nr)) !start element #
      PATHWAY(TOTAL_PATHS,2)=TRANS_TIME(NENQ(1,NQ_START(nr)))
      PATHWAY(TOTAL_PATHS,3)=DL(3,NLL(1,NENQ(1,NQ_START(nr))))
      PATHWAY(TOTAL_PATHS,4)=1.d0
      CONTINU=.TRUE.
C... Currently, array PATHWAY stores all possible flow pathways     
      DO WHILE(CONTINU.AND.TOTAL_PATHS.LT.MAX_PATH)
C        count=count+1                         
        path=0
        DO WHILE(path.LT.TOTAL_PATHS)
          path=path+1
          NEW=.FALSE.
          ne=PATH_END(path) !element # at end of previous path
c          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          IF(NENP(np2,0,nr).GT.1) THEN !more possible flow pathways
            DO noelem=1,NENP(np2,0,nr) 
              ne2=NENP(np2,noelem,nr) !global element #
              IF(ne.NE.ne2) THEN !makes sure doesn't use inlet ne
                np3=NPNE(2,nb,ne2) !other node of possible pathway
                IF(np3.EQ.np2) np3=NPNE(1,nb,ne2) !element, ne2
                IF(DABS(YQ(NYNQ(3,NQNE(ne2,1),0),1)).GT.TOL) THEN 
C... velocity at 1st grid point of ne2
                  IF(NEW.AND.TOTAL_PATHS.LT.MAX_PATH) THEN
                    TOTAL_PATHS=TOTAL_PATHS+1
                    DO j=1,FACTORS !copies previous entries into new
                      !NB/ CURRENLTY ONLY STORING 4 FACTORS (FACTORS=7)
                      PATHWAY(TOTAL_PATHS,j)=PATHWAY(path,j) !pathway
                    ENDDO !j
C... have to take away the values of the ne just done
                    PATHWAY(TOTAL_PATHS,2)=PATHWAY(path,2)-
     '                TRANS_TIME(PATH_END(path))
                    PATHWAY(TOTAL_PATHS,3)=PATHWAY(path,3)-
     '                DL(3,NLL(1,PATH_END(path)))
                    PATHWAY(TOTAL_PATHS,4)=PATHWAY(path,4)-1.d0
                    path_number=TOTAL_PATHS
                  ELSE
                    path_number=path
                  ENDIF !NEW
C                  PATHWAY(0,path_number)=PATHWAY(0,path_number)+1
c                 PATHWAY(path_number,1)=DBLE(ne2) !stores element #
                  PATH_END(path_number)=ne2 !stores element #
                  PATHWAY(path_number,2)=PATHWAY(path_number,2)+
     '              TRANS_TIME(ne2) 
                  PATHWAY(path_number,3)=PATHWAY(path_number,3)+
     '               DL(3,NLL(1,ne2))
                  PATHWAY(path_number,4)=PATHWAY(path_number,4)+1.d0
                  NEW=.TRUE.
                ELSE !if velocity=0 remove pathway (incomplete pathway)
                  
                ENDIF !velocity.GT.0
              ENDIF !ne.NE.ne2
            ENDDO !noelem
          ENDIF !NENP(np2,0,nr).GT.1
        ENDDO !WHILE path
        CONTINU=.FALSE. !each pathway has reached outlet
C        PATHWAY(0,4)=0.d0 
        DO path=1,TOTAL_PATHS
          IF(NENP(NPNE(2,nb,PATH_END(path)),0,nr).GT.1) THEN
            DO noelem=1,NENP(NPNE(2,nb,PATH_END(path)),0,nr)
              ne=NENP(NPNE(2,nb,PATH_END(path)),noelem,nr)
              IF(ne.NE.PATH_END(path)) THEN
                IF(DABS(YQ(NYNQ(3,NQNE(ne2,1),0),1)).GT.TOL)
     '            CONTINU=.TRUE. !unfinished pathways remain
              ENDIF
            ENDDO
          ENDIF
        ENDDO !path
      ENDDO !WHILE CONTINU
      
      IF(TOTAL_PATHS.EQ.MAX_PATH.AND.CONTINU) THEN
        WRITE(OP_STRING,'('' INCREASE MAX # PATHWAYS ALLOWED '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C.. Open file to print pathway results to:
      CALL OPENF(IOFILE2,'DISK','transit.out','NEW',
     '  'SEQUEN','FORMATTED',132,ERROR,*9999)
      WRITE(OP_STRING,
     '  '('' Path | Grid pt |   x   |   y   |   z  |Path length|'//
     '  'T-time| Pressure | Radius | Velocity | Flow |'//
     &  ' # elems '')')
      CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)
      
      DO path=1,TOTAL_PATHS
        np2=NPNE(2,nb,PATH_END(path))
        ne=PATH_END(path)
        nq=NQNE(ne,NQET(NQS(ne)))
        IF(NENP(np2,0,nr).EQ.1) THEN
          !only do this if pathway complete
          PATHWAY(0,2)=PATHWAY(0,2)+PATHWAY(path,2) !sum transit time
          PATHWAY(0,3)=PATHWAY(0,3)+PATHWAY(path,3) !sum path length
          PATHWAY(0,4)=PATHWAY(0,4)+PATHWAY(path,4) !total ne
          flow=YQ(NYNQ(3,nq,0),1)*PI*YQ(NYNQ(2,nq,0),1)**2.d0 !Q=V*pi*R**2
          WRITE(IOFILE2,'(I5,1X,I6,3(2X,F7.2),5(1X,F8.3),1X,F8.6,1X,
     &      I5)') path,nq,XQ(1,nq),XQ(2,nq),XQ(3,nq),PATHWAY(path,3),
     &      PATHWAY(path,2),YQ(NYNQ(1,nq,0),1),YQ(NYNQ(2,nq,0),1),
     &      YQ(NYNQ(3,nq,0),1),flow,INT(PATHWAY(path,4))
        ENDIF
      ENDDO !path
C... close output file      
      CALL CLOSEF(IOFILE2,ERROR,*9999)
      
      PATHWAY(0,2)=PATHWAY(0,2)/TOTAL_PATHS !avge transit time
      PATHWAY(0,3)=PATHWAY(0,3)/TOTAL_PATHS !avge path length
      PATHWAY(0,4)=PATHWAY(0,4)/TOTAL_PATHS !mean ne's in a pathway
C      IF(DABS(PATHWAY(0,2)).GT.0.d0) THEN
        WRITE(OP_STRING,'('' Average transit time (s) = '',F15.5)')
     '    PATHWAY(0,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average path length (mm) = '',F15.5)')
     '    PATHWAY(0,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average # elements in path = '',F12.2)')
     '    PATHWAY(0,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDIF
        
      IF(VENOUS_NETWORK.EQ.'Y'.AND.N_VENOUS_GEOM.EQ.1) THEN
C... Calculate transit time for identical venous geometry
c        count=1
        TOTAL_PATHS=1       
c        PATHWAY(TOTAL_PATHS,1)=DBLE(NENQ(1,NQ_START(nr))) !start element #
        PATH_END(TOTAL_PATHS)=NENQ(1,NQ_START(nr)) !start element #
        PATHWAY(TOTAL_PATHS,2)=TRANS_TIME_VEIN(NENQ(1,NQ_START(nr)))
        PATHWAY(TOTAL_PATHS,3)=DL(3,NLL(1,NENQ(1,NQ_START(nr))))
        PATHWAY(TOTAL_PATHS,4)=1.d0
        CONTINU=.TRUE.
C... Currently, array PATHWAY stores all possible flow pathways     
        DO WHILE(CONTINU.AND.TOTAL_PATHS.LT.MAX_PATH)
          path=0
          DO WHILE(path.LT.TOTAL_PATHS)
            path=path+1
            NEW=.FALSE.
            ne=PATH_END(path) !element # at end of previous path
c            np1=NPNE(1,nb,ne)
            np2=NPNE(2,nb,ne)
            IF(NENP(np2,0,nr).GT.1) THEN !more possible flow pathways
              DO noelem=1,NENP(np2,0,nr) 
                ne2=NENP(np2,noelem,nr) !global element #
                IF(ne.NE.ne2) THEN !makes sure doesn't use inlet ne
                  np3=NPNE(2,nb,ne2) !other node of possible pathway
                  IF(np3.EQ.np2) np3=NPNE(1,nb,ne2) !element, ne2
                  IF(DABS(YQ(NYNQ(3,NQNE(ne2,1),0),1)).GT.TOL) THEN 
                    IF(NEW.AND.TOTAL_PATHS.LT.MAX_PATH) THEN
                      TOTAL_PATHS=TOTAL_PATHS+1
                      DO j=1,FACTORS !copies previous entries into new
                        PATHWAY(TOTAL_PATHS,j)=PATHWAY(path,j) !pathway
                      ENDDO !j
                      PATHWAY(TOTAL_PATHS,2)=PATHWAY(path,2)-
     '                  TRANS_TIME_VEIN(PATH_END(path))
                      PATHWAY(TOTAL_PATHS,3)=PATHWAY(path,3)-
     '                  DL(3,NLL(1,PATH_END(path)))
                      PATHWAY(TOTAL_PATHS,4)=PATHWAY(path,4)-1.d0
                      path_number=TOTAL_PATHS
                    ELSE
                      path_number=path
                    ENDIF !NEW
c                    PATHWAY(path_number,1)=DBLE(ne2) !stores element #
                    PATH_END(path_number)=ne2 !stores element #
                    PATHWAY(path_number,2)=PATHWAY(path_number,2)+
     '                TRANS_TIME_VEIN(ne2) 
                    PATHWAY(path_number,3)=PATHWAY(path_number,3)+
     '                DL(3,NLL(1,ne2))
                    PATHWAY(path_number,4)=PATHWAY(path_number,4)+1.d0
                    NEW=.TRUE.
C                  ELSE !if velocity=0 remove pathway (incomplete pathway)
                  ENDIF !velocity.GT.0
                ENDIF !ne.NE.ne2
              ENDDO !noelem
            ENDIF !NENP(np2,0,nr).GT.1
          ENDDO !WHILE path
          CONTINU=.FALSE. !each pathway has reached outlet
C         PATHWAY(0,4)=0.d0 
          DO path=1,TOTAL_PATHS
            IF(NENP(NPNE(2,nb,PATH_END(path)),0,nr).GT.1) THEN
              DO noelem=1,NENP(NPNE(2,nb,PATH_END(path)),0,nr)
                ne=NENP(NPNE(2,nb,PATH_END(path)),noelem,nr)
                IF(ne.NE.PATH_END(path)) THEN
                  IF(DABS(YQ(NYNQ(3,NQNE(ne2,1),0),1)).GT.TOL)
     '              CONTINU=.TRUE. !unfinished pathways remain
                ENDIF
              ENDDO
            ENDIF
          ENDDO !path
        ENDDO !WHILE CONTINU
        IF(TOTAL_PATHS.EQ.MAX_PATH.AND.CONTINU) THEN
          WRITE(OP_STRING,'('' INCREASE MAX # PATHWAYS ALLOWED '')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO path=1,TOTAL_PATHS
          IF(NENP(NPNE(2,nb,PATH_END(path)),0,nr).EQ.1) THEN
 !only do this if pathway complete
            PATHWAY(0,2)=PATHWAY(0,2)+PATHWAY(path,2) !sum transit time
            PATHWAY(0,3)=PATHWAY(0,3)+PATHWAY(path,3) !sum path length
            PATHWAY(0,4)=PATHWAY(0,4)+PATHWAY(path,4) !total ne
          ENDIF
        ENDDO !path
        PATHWAY(0,2)=PATHWAY(0,2)/TOTAL_PATHS !avge transit time
        PATHWAY(0,3)=PATHWAY(0,3)/TOTAL_PATHS !avge path length
        PATHWAY(0,4)=PATHWAY(0,4)/TOTAL_PATHS !mean ne's in a pathway
        WRITE(OP_STRING,
     &    '('' Average venous transit time (s) = '',F15.5)')
     &    PATHWAY(0,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     &    '('' Average venous path length (mm) = '',F15.5)')
     &    PATHWAY(0,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     &    '('' Average # elements in venous path = '',F12.2)')
     &    PATHWAY(0,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF !identical venous geometry solution
      
      
      CALL EXITS('TRANSIT')
      RETURN
 9999 CALL ERRORS('TRANSIT',ERROR)
      CALL EXITS('TRANSIT')
      RETURN 1
      END


