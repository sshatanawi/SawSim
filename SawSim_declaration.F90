module SawSim_declaration
! *******************************************************************
! Contains all data decalarations(domain, arrays, BC information)
! Parameters / Domain
! *******************************************************************
    use petsc
    implicit none

    ! ---------------------------------------------------------------------
    ! Domain geometry From Noah MP
    ! ---------------------------------------------------------------------
    integer         :: NCOL         ! Number of columns
    integer         :: NROW         ! Number of rows
    real(kind=8)    :: T_SOIL_THK   ! Total thickness of soil layer in (NOAH MP)
    real(kind=8)    :: DELTAC       ! Along column, the width of row i (in the y direction)
    real(kind=8)    :: DELTAR       ! Along row, the width of column j (in the x direction)
    real(kind=8), allocatable,dimension (:,:)       ::  TERRAIN ! Ground Surface Elevations

    ! ---------------------------------------------------------------------
    !  Groundwater geometry and properties (2 layers: top unconfined, bottom confined)
    ! ---------------------------------------------------------------------

    integer         :: NLAY         ! Number of groundwater layers (2 layers only)
    real(kind=8)    :: VAD_THK      ! Vadose layer thickness (m) (the 1/2 of Soil Thickness, MMF scheme)
    real(kind=8)    :: BED_THK      ! Confining bed thickness (m)
    real(kind=8)    :: BED_K        ! Confining bed vertical hydraulic conductivity (m/s)
    real(kind=8)    :: UNCONF_THK   ! Confined aquifer thickness (m) if the WT at bottom of Vadose layer
    real(kind=8)    :: CONF_AQ_THK  ! Confined aquifer thickness (m)
    real(kind=8), allocatable,dimension (:,:,:)     :: KR !Hydraulic conductivity along row (assoicated with I) in x direction
    real(kind=8), allocatable,dimension (:,:,:)     :: KC !Hydraulic conductivity along column(assoicated of J) in y direction
    real(kind=8), allocatable,dimension (:,:,:)     :: KV !Vertival hydraulic conductivity
    real(kind=8), allocatable,dimension (:,:)       :: SY !Specific Yield for the unconfined aquifer
    real(kind=8), allocatable,dimension (:,:)       :: SS !Specific Storage for the confined aquifer
    ! ---------------------------------------------------------------------
    ! Global size
    ! ---------------------------------------------------------------------
    integer         :: N            ! Number of groundwater total grid cells NCOL x NROW X NLAY 
    
    ! ---------------------------------------------------------------------
    ! Time stepping
    ! ---------------------------------------------------------------------
    integer         :: NTSP         ! Number of time steps
    real(kind=8)    :: DELTAT       ! time step size. For the groundwater calcuation it should be covert to seconds
    integer         :: MAX_ITER     ! Number of iterations for picard iteration
    real(kind=8)    :: picard_tol    ! Picard iteration tolerance 
    real(kind=8)    :: Time_start   ! The begining of time series
    
    ! ---------------------------------------------------------------------
    ! Boundary condition types and values
    ! ---------------------------------------------------------------------
    
    integer, parameter  :: BC_NONE   =   0 ! Not a boundary condition (interior)
    integer, parameter  :: BC_FIXED  =   1 ! Constant head boundary condition
    integer, parameter  :: BC_NOFLOW =   2 ! No flow boundary condition
    integer, allocatable,dimension(:,:,:)       :: BC_type     
    real(kind=8), allocatable,dimension(:,:,:)  :: BC_value 

    ! ---------------------------------------------------------------------
    ! Source types
    ! ---------------------------------------------------------------------
    integer, parameter  :: SRC_NONE  = 0  ! Cell has NO external source
    integer, parameter  :: SRC_WELL  = 1  ! Cell contains a well
    integer, parameter  :: SRC_RIVER = 2  ! Cell contains a river
    integer, allocatable,dimension(:,:,:)      :: source_type   
    real(kind=8), allocatable,dimension(:,:,:) :: source_value
    
    ! ---------------------------------------------------------------------
    ! Initial condition ( Input or Noah MP's spin up)
    ! ---------------------------------------------------------------------
    real(kind=8), allocatable, dimension (:,:,:)      ::  H_initial ! initial head
       
    ! ---------------------------------------------------------------------
    ! External source/sink : Recharge, pumping, River Exchange
    ! ---------------------------------------------------------------------
        ! recharge rates(for testing)
    real(kind=8), allocatable,dimension (:,:)       :: recharge_rate    
        ! recharge from Noah MP
    real(kind=8), allocatable, dimension (:,:)      :: WTD              ! Water Table depth (m)
    !real(kind=8),allocatable,dimension(:,:)        :: RECH
    !real(kind=8),allocatable,dimension(:,:)       :: DEEPRECH
    
        ! pumping rate
    real(kind=8), allocatable,dimension (:,:,:,:)     :: pumping_rate     
    ! ---------------------------------------------------------------------
    ! River Propertiers
    ! ---------------------------------------------------------------------
    real(kind=8)    :: K_RIVERBED   ! River bed hydraulic conductivty
    real(kind=8)    :: RIVERBED_THK   ! River bed thickness
    real(kind=8), allocatable, dimension (:,:,:)    :: W_RIVER ! River width
    real(kind=8), allocatable, dimension (:,:,:)    :: L_RIVER ! River length
    real(kind=8), allocatable, dimension (:,:,:)    :: H_RIVER ! River water level
    real(kind=8), allocatable, dimension (:,:,:)    :: RIVBED_BOT ! Riverbed bottom elevation  
    real(kind=8), allocatable, dimension (:,:,:)    :: C_RIVBED ! Riverbed conductance  
    

    ! ---------------------------------------------------------------------
    ! Arrays for variables
    ! --------------------------------------------------------------------- 
    real(kind=8), allocatable,dimension (:)         :: TIME
    real(kind=8), allocatable,dimension (:,:,:,:)   :: GW_HEAD
    real(kind=8), allocatable, dimension (:,:)      :: DELTAV       ! thickness of saturated layer in the unconfined aquifer only 
    real(kind=8), allocatable, dimension (:,:,:)    ::  HNEW        ! The current timestep for peizometric head
    real(kind=8), allocatable, dimension (:,:,:)    ::  HOLD        ! The previous timstep for peizometric head
    real(kind=8), allocatable, dimension (:,:,:)    :: GWTOP        ! The elevation of the top of GW layer
    real(kind=8), allocatable, dimension (:,:,:)    :: GWBOT        ! The elevation of the bottom of GW layer
    real(kind=8), allocatable, dimension (:,:,:)    :: CC           ! Conducatnce along column
    real(kind=8), allocatable, dimension (:,:,:)    :: CR           ! Conducatnce along row
    real(kind=8), allocatable, dimension (:,:,:)    :: CV           ! Vertical Conducatnce
    real(kind=8), allocatable, dimension (:,:,:)    :: B_RIVER      ! The constant coefficient from GW_SW exchange
    real(kind=8), allocatable, dimension (:,:,:)    :: P_RIVER      ! The head dependent coefficient from GW_SW exchange

    ! ---------------------------------------------------------------------
    ! Input files paths
    ! ---------------------------------------------------------------------    
    character(len=200)  ::  K_FILE              ! Hydraulic conductivity
    character(len=200)  ::  STORAGE_FILE        ! Aquifer storage prop.
    character(len=200)  ::  SOURCE_FILE         ! Source file
    character(len=200)  ::  BC_FILE             ! Boundary conditions file
    character(len=200)  ::  HEAD0_FILE           ! initial conditons
    character(len=200)  ::  RIVER_FILE          ! River geomotry file
    character(len=200)  ::  WELL_FILE           ! WELL file

    ! --------------------------------------------------------------------
    ! WATER BUDGET variables 
    ! Arrays for water fluxes (per grid cell), sotage changes, and net flux
    ! --------------------------------------------------------------------
    !real(kind=8), allocatable :: Q_recharge(:,:,:) ! Water entering from recharge (calculated in your GW_recharge subroutine)
    !real(kind=8), allocatable :: Q_pumping(:,:,:) ! Water removed or added by wells (from your pumping data)
    !real(kind=8), allocatable :: Q_river(:,:,:) ! Exchange flux calculated in GW_RIVER_interaction
    !real(kind=8), allocatable :: Q_lat_xri(:,:,:) ! Lateral flux, net flow in the x-direction computed using the conductance array CR and head differences 
    !real(kind=8), allocatable :: Q_lat_ycj(:,:,:) ! Lateral flux, net flow in the y-direction computed using the conductance array CC and head differences 
    !real(kind=8), allocatable :: Q_vert(:,:,:) ! Flow between layers using CV and vertical head differences 
    !real(kind=8), allocatable :: dStorage_m3(:,:,:) ! storage change (ΔStorage)
    !real(kind=8), allocatable :: Net_flow(:,:,:)    ! net flux (total flux) per grid cell

    real(kind=8), allocatable :: recharge_flux(:,:,:) ! Water entering from recharge (calculated in your GW_recharge subroutine) unit m
    real(kind=8), allocatable :: pumping_flux(:,:,:) ! Water removed or added by wells (from your pumping data)
    real(kind=8), allocatable :: river_flux(:,:,:) ! Exchange flux calculated in GW_RIVER_interaction
    real(kind=8), allocatable :: lat_xri_flux(:,:,:) ! Lateral flux, net flow in the x-direction computed using the conductance array CR and head differences 
    real(kind=8), allocatable :: lat_ycj_flux(:,:,:) ! Lateral flux, net flow in the y-direction computed using the conductance array CC and head differences 
    real(kind=8), allocatable :: vertical_flux(:,:,:) ! Flow between layers using CV and vertical head differences 
    real(kind=8), allocatable :: dStorage(:,:,:) ! storage change (ΔStorage)
    real(kind=8), allocatable :: Net_flux(:,:,:)    ! net flux (total flux) per grid cell
end module SawSim_declaration
