!
! Parameters.Inc
!
! Contains all parameter declarations for CMISS.  If this file
! is in a local directory, then change the logical cmiss$parameters
! to point to that directory.
!

C#### NAMX is maximum number of auxiliary parameters
C###  NBMX is    "      "     " basis functions
C###  NCMX is    "      "     " #var types for a given dependent var
C###  NDMX is    "      "     " data points
C###  NEGMX is the maximum size of a Gauss point group array
C###  NEMX is    "      "     " elements
C###  NE_R_MX is "      "     " elements in a region
C###  NEIMX is   "      "     " adjacent elements in an Xi direction
C###  NELMX is   "      "     " elements along a line
C###  NFMX is    "      "     " global face segments
C###  NFVCMX is  "      "     " max number of voronoi cell faces
C###  NFVMX is   "      "     " max number of voronoi cell connections
C###  NF_R_MX is "      "     " face segments in a region
C###  NGMX is    "      "     " Gauss points per element
C###  NGRSEGMX is       "     " number of graphical segments
C###  NHMX is    "      "     " dependent variables
C###  NIMX is    "      "     " local Xi-coordinates
C###  NISR_G?MX is      "     " Amount of sparse row storage
C###  NISC_G?MX is      "     " Amount of sparse column storage
C###  NJMX is    "      "     " global reference coordinates
C###  NKMX is    "      "     " derivatives per variable
C###  NLMX is    "      "     " global line segments
C###  NL_R_MX is "      "     " lines in a region
C###  NMMX is    "      "     " material parameters
C###  NNMX is    "      "     " element nodes
C###  NOMX is    "      "     " degrees of freedom
C###  NPMX is    "      "     " global nodes
C###  NP_R_MX is "      "     " global nodes in a region
C###  NQMX is    "      "     " global grid points
C###  NYQMX is    "      "    " global grid degrees dof's
C###  NRMX is    "      "     " regions & contours per elem & vari segs
C###  NSMX is    "      "     " element dofs per variable (<=NNMX*NKMX)
C###  NSSMX is   "      "     " signal signal sets
C###  NSFMX is   "      "     " face dofs per variable (<=NNMX*NKMX)
C###  NTMX is    "      "     " eigenvalues
C###  NTSMX is   "      "     " time samples
C###  NUMX is    "      "     " derivative terms up to 2nd order
C###  NVMX is    "      "     " # of versions of a particular value.
C###  NVCBMX is  "      "     " max number of voronoi boundary nodes
C###  NVCMX is   "      "     " max number of voronoi cells
C###  NWMX is    "      "     " workstations for segment arrays (# graphic wins)
C###  NXMX is    "      "     " the number of problem types
C###  NYMX is    "      "     " mesh dofs
C###  NY_R_MX is "      "     " mesh dofs in a region
C###  NZMX is    "      "     " coeff.s in 1D global stiffness vector
C###  NZ_arrayname_MX   "     " nz dimension of arrayname
C###  NZ_MINOS   "      "     " size of the minos sparse arrays
C###  NBFMX is   "      "     " basis function families
C###  NCOMX is   "      "     " num of non-lin constraint for opti probs
C###  NDEMX is   "      "     " data points in one element
C###  NDIPOLEMX is "    "     " dipoles
C###  NDIPTIMMX is "    "     " dipole time points
C###  NIQMX is   "      "     " variables per grid point
C###  NIQSMX is  "      "     " cell variables (state+derived)
C###  NIFEXTMX is"      "     " number of storage locations for FEXT
C###  NIYMX is   "      "     " variables per mesh dof
C###  NIYFIXMX is"      "     " number of storage locations for FIX
C###  NIYGMX is  "      "     " variables stored for gauss points
C###  NIYGFMX is "      "     " variables stored for face gauss points
C###  NLCMX is   "      "     " num of lin constraints for opti probs
C###  NMAQMX is  "      "     " auxiliary grid parameters
C###  NMQMX is   "      "     " cell material parameters
C###  NOPMX is   "      "     " optimisation variables
C###  NORMX is   "      "     " size of fractal tree order arrays
C###  NOYMX is   "      "     " solution dofs attached to one mesh dof
C###  NOOPMX is  "      "     " optimisation or solution variables
C###  NPDMX is   "      "     " domain nodes for BE problems
C###  NQEMX is   "      "     " grid points in any element
C###  NQGMX is   "      "     " non-zeros in grid matrix row
C###  NQIMX is   "      "     " cell integer parameters
C###  NQRMX is   "      "     " cell real parameters
C###  NQISVMX is "      "     " spatially varying cell int. parameters
C###  NQRSVMX is "      "     " spatially varying cell real parameters
C###  NQSCMX is  "      "     " grid point schemes
C###  NQVMX is   "      "     " cell variants
C###  NRCMX is   "      "     " number of rows and columns (s.b. 2)
C###  NREMX is   "      "     " residuals for optimization problems
C###  NTIMEPOINTSMX is ""     " time points for a time variable
C###  NTIMEVARSMX is "  "     " time variables
C###  NYOMX is   "      "     " mesh dofs attached to one solution dof
C###  MROWMX is  "      "     " rows in a problem
C###  NIMAGEMX is       "     " image cell array dimension
C###  NY_TRANSFER_MX is maximum size of the transfer matrix
C###  NYYMX is   "      "     " mesh dofs mapped to one mesh dof
C###
C###  USE_BEMX       is 1 if bem is to be used, 0 if not
C###  USE_CELLX      is 1 if cell models are to be used, 0 if not
C###  USE_DATAX      is 1 if data is to be used, 0 if not
C###  USE_DIPOLEX    is 1 if dipoles are to be used, 0 if not
C###  USE_GAUSS_PT_MATERIALSX  is 1 if any material propterties are to be defined by Gauss points, 0 if not
C###  USE_GRAPHICSX  is 1 if graphics is to be used, 0 if not
C###  USE_GRIDX      is 1 if collocation grids are to be used, 0 if not
C###  USE_LUNGX      is 1 if lung problems are to be used, 0 if not
C###  USE_MINOSX     is 1 if minos is to be used, 0 if not
C###  USE_NLSTRIPEX  is 1 if non-lin stripe fit is to be used, 0 if not
C###  USE_NONLINX    is 1 if non-linear problems are to be used, 0 if not
C###  USE_NPSOLX     is 1 if npsol is to be used, 0 if not
C###  USE_SPARSEX    is 1 if sparse arrays are to be used, 0 if not
C###  USE_TIMEX      is 1 if time-dep arrays are used for storage
C###  USE_TRANSFERX  is 1 if a transfer matrix is to be used/calculated
C###  USE_VORONOIX   is 1 if if Voronoi cells are to be used
C****
C****
C###  Note: NHMX must be equal to NJMX for problems where the dependent
C###        variable array carries deformed coordinates (eg FE50 problems).
C###  Note: NZMX must be > NOMX^2 and NYMX^2 for some problems. It no longer
C###        longer needs to be large enough to take the fill in during sparse
C###        matrix solution since Harwell/Umfpack-2.0 are not used.
C###  Note: NCOMX used to be 20. Currently set to 0 since we don't have 
C###        nonlinear constraints with MinLSSQP.

      INTEGER NAMX,NBFMX,NBMX,NCMX,NCOMX,NCVMX,NDMX,NDEMX,
     '  NDIPOLEMX,
     '  NDIPTIMMX,NEGMX,NEMX,NEIMX,NELMX,NEPMX,
     '  NNEPMX,NEQMX,NE_R_MX,NFVCMX,NFVMX,
     '  NFMX,NF_R_MX,NGMX,NGRSEGMX,NHMX,NIMX,NIQMX,NIQSMX,NISC_GDMX,
     '  NISC_GKMX,NISC_GKKMX,NISC_GMMX,NISC_GMMMX,NISC_GQMX,NISR_GDMX,
     '  NISR_GKMX,NISR_GKKMX,NISR_GMMX,NISR_GMMMX,NISR_GQMX,
     '  NIFEXTMX,NIYMX,NIYFIXMX,NIYGMX,NIYGFMX,NJMX,NKMX,NLMX,NLISTMX,
     '  NL_R_MX,NMAQMX,NMMX,NMQMX,NNMX,NOMX,NOPMX,NOOPMX,
     '  NORMX,NOYMX,NPMX,NPDMX,NP_R_MX,NQEMX,NQGMX,
     '  NQMX,NQIMX,NQISVMX,NQSCMX,
     '  NQRMX,NQRSVMX,NQVMX,NRMX,NRCMX,NREMX,
     '  NSFMX,NSSMX,NSMX,NTIMEPOINTSMX,NTIMEVARSMX,NTMX,NTSMX,NUMX,NVMX,
     '  NVCBMX,NVCMX,NWMX,
     '  NXMX,NYMX,NYOMX,
     '  NYQMX,NY_R_MX,NYROWMX,NZMX,NZ_FNY_MX,NZ_GD_MX,
     '  NZ_GK_MX,NZ_GKK_MX,
     '  NZ_GM_MX,NZ_GMM_MX,NZ_GQ_MX,
     '  NLCMX,NZ_MINOS,NIMAGEMX,NY_TRANSFER_MX,NYYMX,
     '  USE_BEMX,USE_CELLX,USE_DATAX,USE_DIPOLEX,
     '  USE_GAUSS_PT_MATERIALSX,USE_GRAPHICSX,USE_GRIDX,USE_LUNGX,
     '  USE_MAGNETICX,
     '  USE_MINOSX,USE_NLSTRIPEX,USE_NONLINX,USE_NPSOLX,USE_SPARSEX,
     '  USE_TIMEX,USE_TRANSFERX,USE_VORONOIX
      PARAMETER
     ' (NAMX=5,          NBMX=42,      NBFMX=13,     NCMX=2,
     '  NCOMX=0,         NCVMX=3,      NDMX=1700,    NDEMX=500,
     '  NDIPOLEMX=20,    NDIPTIMMX=500,NEMX=120,     NEIMX=4,
     '  NEPMX=20,        NE_R_MX=NEMX, NFMX=220,     NF_R_MX=NFMX,
     '  NGMX=81,         NHMX=6,       NIMX=3,       NIMAGEMX=1,
     '  NIQMX=6,         NIQSMX=50,    NIFEXTMX=8,   NIYMX=16,
     '  NIYFIXMX=5,      NIYGMX=6,     NIYGFMX=0,    NFVCMX=6,
     '  NFVMX=300,       NJMX=6,       NKMX=8,       NLMX=259,
     '  NLCMX=1,         NL_R_MX=NLMX, NMMX=35,      NMQMX=99,
     '  NNMX=12,         NMAQMX=6,     NOMX=450,     NOPMX=125,
     '  NORMX=20,        NOYMX=1,      NPMX=200,     NPDMX=1,
     '  NP_R_MX=NPMX,    NQEMX=81,     NQGMX=22,     NQMX=4225,
     '  NQIMX=1,
     '  NQISVMX=1,       NQRMX=1,      NQRSVMX=1,    NQSCMX=9,
     '  NQVMX=1,         NRMX=6,       NRCMX=2,      NVCBMX=100,
     '  NVCMX=NPMX,      NREMX=400,    NSFMX=16,     NSMX=64,
     '  NSSMX=1,
     '  NTIMEPOINTSMX=1, NTIMEVARSMX=1,NTMX=20,      NTSMX=200,
     '  NUMX=11,         NVMX=16,      NWMX=3,       NXMX=3,
     '  NYMX=450,        NYOMX=4,      NYQMX=1,      NY_R_MX=NYMX, !only one region
     '  NYROWMX=NY_R_MX, NY_TRANSFER_MX=200,NYYMX=0,
     '  NZMX=NYROWMX*NYROWMX,
     '  NZ_FNY_MX =1,
     '  NZ_GD_MX  =NZMX,
     '  NZ_GK_MX  =NZMX,
     '  NZ_GKK_MX =NOMX*NOMX,
     '  NZ_GM_MX  =NZMX,
     '  NZ_GMM_MX =NOMX*NOMX,
     '  NZ_GQ_MX  =NZMX,
!       Setting sparsity array sizes
!       See above sparsity/solver comment before changing sizes
     '  NISC_GDMX  =NZ_GD_MX,           !compressed row
     '  NISR_GDMX  =NYROWMX+1,          !compressed row
     '  NISC_GKMX  =NZ_GK_MX,           !compressed row
     '  NISR_GKMX  =NYROWMX+1,          !compressed row
     '  NISC_GKKMX =NZ_GKK_MX,          !compressed row
     '  NISR_GKKMX =NOMX+1,             !compressed row
     '  NISC_GMMX  =NZ_GM_MX,           !compressed row
     '  NISR_GMMX  =NYROWMX+1,          !compressed row
     '  NISC_GMMMX =NZ_GMM_MX,          !compressed row
     '  NISR_GMMMX =NOMX+1,             !compressed row
     '  NISC_GQMX  =NZ_GQ_MX,           !compressed row
     '  NISR_GQMX  =NYROWMX+1,          !compressed row
     '  NZ_MINOS=NOPMX*(NCOMX+NLCMX),
     '  NELMX=16,
! CS 2/7/98 The MAX function should not be used here, instead
! initialise parameters with the correct values directly
!       NNEPMX=MAX(NEMX,NPMX)
!       NEQMX=MAX(NEMX,NQMX)
!       NLISTMX=MAX(NEMX,NPMX,NLMX,NDMX)
!       NOOPMX=MAX(NOMX,NOPMX)
     '  NNEPMX=NPMX,
     '  NEGMX=1+NEMX*(NGMX+2),
     '  NEQMX=NQMX,
     '  NGRSEGMX=6,
     '  NLISTMX=NDMX,
     '  NOOPMX=NOMX,
     '  USE_BEMX=1,          USE_CELLX=1,         USE_DATAX=1,
     '  USE_DIPOLEX=1,       USE_GAUSS_PT_MATERIALSX=0,USE_GRAPHICSX=1,
     '  USE_GRIDX=1,         USE_LUNGX=1,         USE_MAGNETICX=0,
     '  USE_MINOSX=0,        USE_NLSTRIPEX=1,     USE_NONLINX=1,
     '  USE_NPSOLX=1,        USE_SPARSEX=1,       USE_TIMEX=1,
     '  USE_TRANSFERX=1,     USE_VORONOIX=0)


