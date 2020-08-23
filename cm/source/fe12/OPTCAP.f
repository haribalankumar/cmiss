      SUBROUTINE OPTCAP(INLET,MAX_PATH_I,TOTAL_PATHS,AVGE_WBC,
     '  NO_BLOCK,PATHWAY,ERROR,*)

C#### Subroutine: OPTCAP
C###  Description:
C###    OPTCAP outputs red blood cell (RBC) and white blood
C###    cell (WBC) transit time results through the pulmonary
C###    capillaries, to screen and to file.

C***  Created by KSB, 29 November 2002.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER INLET,MAX_PATH_I,TOTAL_PATHS
      REAL*8 AVGE_WBC,NO_BLOCK,PATHWAY(0:MAX_PATH,FACTORS)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,IBEG,IEND,path,NO_INT,tint
      REAL*8 high,length_avge,low,max_time,max_time_wbc,min_time,
     '  min_time_wbc,NO_STUCK_PATH,PLOT_DATA(300,FACTORS),SIZE_INT,
     '  stuck,WBC_TIME
      CHARACTER CHAR1*1,STRING*(50)

      CALL ENTERS('OPTCAP',*9999)

      length_avge=0.d0
      NO_STUCK_PATH=0.d0
      WBC_TIME=0.d0
      SIZE_INT=0.10d0
      low=0.d0
      high=30.d0
      NO_INT=NINT((high-low)/SIZE_INT)
      DO tint=1,NO_INT
        DO i=1,FACTORS
          PLOT_DATA(tint,i)=0.d0 !initialise
        ENDDO
      ENDDO
      DO path=(TOTAL_PATHS-MAX_PATH_I+1),TOTAL_PATHS
        PATHWAY(0,3)=PATHWAY(0,3)+PATHWAY(path,3) !sums RBC fraction
      ENDDO !path
      DO path=(TOTAL_PATHS-MAX_PATH_I+1),TOTAL_PATHS !calculating averages
        PATHWAY(0,1)=PATHWAY(0,1)+PATHWAY(path,2) !average,not weighted
        PATHWAY(0,2)=PATHWAY(0,2)+(PATHWAY(path,3)/PATHWAY(0,3))
     '    *PATHWAY(path,2) !weighted average
        PATHWAY(0,4)=PATHWAY(0,4)+PATHWAY(path,4) !sums path length
        length_avge=length_avge+(PATHWAY(path,3)/PATHWAY(0,3))
     '    *PATHWAY(path,4) !weighted average
        PATHWAY(0,5)=PATHWAY(0,5)+PATHWAY(path,5) !sums # segments
        PATHWAY(0,7)=PATHWAY(0,7)+PATHWAY(path,7)
        IF(PATHWAY(path,7).GT.0.d0) NO_STUCK_PATH=NO_STUCK_PATH+1.d0
        WBC_TIME=WBC_TIME+PATHWAY(path,6) !sums WBC time
        PATHWAY(0,6)=PATHWAY(0,6)+(PATHWAY(path,3)/PATHWAY(0,3))*
     '    PATHWAY(path,6) !flow weighted average WBC transit time
C... putting transit time data into intervals for plotting distribution
        DO tint=1,NO_INT
          IF(PATHWAY(path,2).GE.(low+(tint-1)*SIZE_INT).AND.
     '      PATHWAY(path,2).LT.(low+tint*SIZE_INT)) THEN !inside interval
            PLOT_DATA(tint,1)=PLOT_DATA(tint,1)+
     '        (PATHWAY(path,3)/PATHWAY(0,3))*100.d0 !% Q
            PLOT_DATA(tint,2)=PLOT_DATA(tint,2)+(1.d0/TOTAL_PATHS)
     '        *100.d0 !# paths in inter
          ENDIF
          IF(PATHWAY(path,6).GE.(low+(tint-1)*SIZE_INT).AND.
     '      PATHWAY(path,6).LT.(low+tint*SIZE_INT)) THEN
            PLOT_DATA(tint,4)=PLOT_DATA(tint,4)+
     '        (PATHWAY(path,3)/PATHWAY(0,3))*100.d0 !% Q
            PLOT_DATA(tint,5)=PLOT_DATA(tint,5)+(1.d0/TOTAL_PATHS)
     '        *100.d0 !# paths in inter
          ENDIF
        ENDDO !tint
      ENDDO !path
C... write results to file
      STRING='time'
      WRITE(CHAR1,'(I1)') INLET
      CALL STRING_TRIM(STRING,IBEG,IEND)
      CALL APPENDC(IEND,CHAR1,STRING)
      CALL OPENF(IOFILE2,'DISK',STRING(IBEG:IEND)//'.out','NEW',
     '  'SEQUEN','FORMATTED',132,ERROR,*9999)
      DO tint=1,NO_INT
        WRITE(IOFILE2,'(3X,5(F12.6))') (tint-1+tint)*SIZE_INT*0.5d0,
     '    PLOT_DATA(tint,1),PLOT_DATA(tint,2),PLOT_DATA(tint,4),
     '    PLOT_DATA(tint,5) !x @ center of interval
      ENDDO
      CALL CLOSEF(IOFILE2,ERROR,*9999)
      PATHWAY(0,1)=PATHWAY(0,1)/MAX_PATH_I !non-weighted avge time
      PATHWAY(0,4)=PATHWAY(0,4)/MAX_PATH_I
      PATHWAY(0,5)=PATHWAY(0,5)/MAX_PATH_I
      WBC_TIME=WBC_TIME/MAX_PATH_I !avge WBC transit time
      PATHWAY(0,7)=PATHWAY(0,7)/MAX_PATH_I
C... write out results to screen
      WRITE(OP_STRING,'('' # pathways for inlet = '',I6)')MAX_PATH_I
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' % of total blood flow = '',F12.4)')
     '  PATHWAY(0,3)*100.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Average # segments in path = '',F12.4)')
     '  PATHWAY(0,5)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Average RBC transit time (s)= '',F12.4)')
     '  PATHWAY(0,2)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Non-weighted RBC transit time(s)= '',F12.4)')
     '  PATHWAY(0,1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Avge weighted path length (mm)= '',F12.4)')
     '  length_avge
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Average path length (mm)= '',F12.4)')
     '  PATHWAY(0,4)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Average RBC velocity (mm/s)= '',F12.4)')
     '  PATHWAY(0,4)/PATHWAY(0,2)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Average WBC transit time (s)= '',F14.4)')
     '  PATHWAY(0,6)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      NO_BLOCK=NO_BLOCK/MAX_PATH_I*100.d0
      WRITE(OP_STRING,'('' % of WBCs stuck in network = '',F12.4)')
     '  NO_BLOCK
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      stuck=(NO_STUCK_PATH/MAX_PATH_I)*100.d0
      AVGE_WBC=2.d0*AVGE_WBC/MAX_PATH_I !average diameter of WBC
      WRITE(OP_STRING,'('' Average WBC diameter (mm) = '',F12.6)')
     '  AVGE_WBC
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' % pathways WBCs get stopped in = '',F12.4)')
     '  stuck
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Average # times WBC stop in path= '',F12.4)')
     '  PATHWAY(0,7)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C... calculating range of transit times
      min_time=1.d6 !inititalise
      max_time=0.d0
      min_time_wbc=1.d6
      max_time_wbc=0.d0
      DO path=(TOTAL_PATHS-MAX_PATH_I+1),TOTAL_PATHS
        min_time=MIN(min_time,PATHWAY(path,2)) !min transit time
        max_time=MAX(max_time,PATHWAY(path,2)) !max transit time
        min_time_wbc=MIN(min_time_wbc,PATHWAY(path,6))
        max_time_wbc=MAX(max_time_wbc,PATHWAY(path,6))
      ENDDO
      WRITE(OP_STRING,'('' RBC transit time range (s)= '',F12.4,'//
     '  ''','',F12.4)') min_time,max_time
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' WBC transit time range (s)= '',F14.4,'//
     '  ''','',F14.4)')  min_time_wbc,max_time_wbc
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)


      CALL EXITS('OPTCAP')
      RETURN
 9999 CALL ERRORS('OPTCAP',ERROR)
      CALL EXITS('OPTCAP')
      RETURN 1
      END

      
