      INTEGER FUNCTION BMG_IJK_MAP_SIZE( Nd, Nx, Ny, Nz )

C =======================================================================
C 
C   BMG_IJK_MAP_SIZE
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_IJK_MAP_SIZE returns size of BMG_IJK_MAP
C
C   Written:    2005/05/30 (TMA)
C
C =======================================================================
C   INPUT:
C ========================
C
C   ----------------------
C   Dimensions:
C   ---------------------- 
C
C   Nd       Dimension of problem
C   Nx       Number of points in x-direction (excluding ghost points)
C   Ny       Number of points in y-direction (excluding ghost points)
C   Nz       Number of points in z-direction (excluding ghost points)
C
C =======================================================================
C   LOCAL:
C ========================
C
C ==========================================================================
      

      IMPLICIT NONE
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'
      
C ========================================================================== 

C ---------------- 
C     Arguments  
C
      INTEGER  Nd, Nx, Ny, Nz

      INTEGER  kg, Nxc, Nyc, Nzc, size

C ========================================================================== 


      IF(Nd.EQ.2) THEN

         kg=1

         size = (Nx+2)*(Ny+2)*2

 10      CONTINUE
      
            Nxc = (Nx-1)/2**kg+1
            Nyc = (Ny-1)/2**kg+1

            size = size + (Nxc+2)*(Nyc+2)*2

            IF (MIN(Nxc,Nyc).GE.2) THEN   
               kg=kg+1
               GOTO 10
            ENDIF
        
 20      CONTINUE
 
         BMG_IJK_MAP_SIZE = size

      ELSEIF(Nd.EQ.3) THEN

         kg=1

         size = (Nx+2)*(Ny+2)*(Nz+2)*3

 30      CONTINUE
      
            Nxc = (Nx-1)/2**kg+1
            Nyc = (Ny-1)/2**kg+1
            Nzc = (Nz-1)/2**kg+1

            size = size + (Nxc+2)*(Nyc+2)*(Nzc+2)*3

            IF (MIN(Nxc,Nyc,Nzc).GE.2) THEN   
               kg=kg+1
               GOTO 30
            ENDIF
        
 40      CONTINUE
 
         BMG_IJK_MAP_SIZE = size

      ENDIF

      RETURN
      END
