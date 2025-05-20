module SawSim_subroutines
!************************************************************************
! Subroutines for:
!  - read_config (GW_config.txt) => scalar domain/time
!  - read_{K, stor, BC, Source} from separate files
!  - build harmonic conduction arrays (Cx,Cy,Cz)
!  - a small function harmonic_mean
!************************************************************************
   use SawSim_declaration
   use petsc
   implicit none

   contains

!**********************************************************************
! SUBROUTINE GW_initialize => 
!   1) read_config => sets NCOL,NROW,NLAY
!   2) allocate arrays
!   3) read K, stor, BC, source from separate files
!   4) build conduction arrays from harmonic means
!**********************************************************************
   subroutine GW_initialize()
   
      implicit none
      
      integer :: i,j,k
      logical :: foundriver ! to read river data

      foundriver = .false.
      
      call read_config()

      ! Calculate N The total size of grid cells  
      N=NCOL*NROW*NLAY
      ! Allocate arrays and initialize to defaults
      if (.not. allocated(KR))then
         allocate(KR(NCOL,NROW,NLAY))
         KR = 0.0
      end if 
      if (.not. allocated(KC))then
         allocate(KC(NCOL,NROW,NLAY))
         KC = 0.0
      end if
      if (.not. allocated(KV))then
         allocate(KV(NCOL,NROW,NLAY))
         KV = 0.0
      end if
      if (.not. allocated(SY))then
         allocate(SY(NCOL,NROW))
         SY = 0.0
      end if
      if (.not. allocated(SS))then
         allocate(SS(NCOL,NROW))
         SS = 0.0
      end if
      if (.not. allocated(BC_type))then
         allocate(BC_type(NCOL,NROW,NLAY))
         BC_type = BC_NONE
      end if
      if (.not. allocated(BC_value))then
         allocate(BC_value(NCOL,NROW,NLAY))
         BC_value = 0.0
      end if
      if (.not. allocated(source_type))then
         allocate(source_type(NCOL,NROW,NLAY))
         source_type = SRC_NONE
      end if
      if (.not. allocated(source_value))then
         allocate(source_value(NCOL,NROW,NLAY))
         source_value = 0.0
      end if
      if (.not. allocated(CR))then
         allocate(CR(NCOL,NROW,NLAY))
         CR = 0.0
      end if
      if (.not. allocated(CC))then
         allocate(CC(NCOL,NROW,NLAY))
         CC = 0.0
      end if
      if (.not. allocated(CV))then
         allocate(CV(NCOL,NROW,NLAY))
         CV = 0.0
      end if
      if (.not. allocated(DELTAV))then
         allocate(DELTAV(NCOL,NROW))
         DELTAV = 0.0
      end if
      if (.not. allocated(HNEW))then
         allocate(HNEW(NCOL,NROW,NLAY))
         HNEW = 0.0
      end if
      if (.not. allocated(HOLD))then
         allocate(HOLD(NCOL,NROW,NLAY))
         HOLD = 0.0
      end if
      if (.not. allocated(WTD))then
         allocate(WTD(NCOL,NROW))
         WTD = 0.0
      end if
      if (.not. allocated(TERRAIN)) then
         allocate(TERRAIN(NCOL,NROW))
         TERRAIN = 0.0
      end if
      if (.not. allocated(recharge_rate)) then
         allocate(recharge_rate(NCOL,NROW))
         recharge_rate = 0.0
      end if
      if (.not. allocated(pumping_rate)) then
         allocate(pumping_rate(NCOL,NROW,NLAY,NTSP))
         pumping_rate = 0.0
      end if
      if (.not. allocated(H_initial)) then
         allocate(H_initial(NCOL,NROW,NLAY))
         H_initial = 0.0
      end if
      if (.not. allocated(GWTOP))then
         allocate (GWTOP(NCOL,NROW,NLAY))
         GWTOP = 0.0
      end if
      if (.not. allocated(GWBOT))then
        allocate(GWBOT(NCOL,NROW,NLAY))
        GWBOT = 0.0
      end if
            if(.not. allocated(W_RIVER)) then
        allocate(W_RIVER(NCOL,NROW,NLAY))
        W_RIVER=0.0
      end if

      if(.not. allocated(L_RIVER)) then
        allocate(L_RIVER(NCOL,NROW,NLAY))
        L_RIVER=0.0
      end if

      if(.not. allocated(H_RIVER)) then
        allocate(H_RIVER(NCOL,NROW,NLAY))
        H_RIVER=0.0
      end if

      if(.not. allocated(RIVBED_BOT)) then
         allocate(RIVBED_BOT(NCOL,NROW,NLAY))
         RIVBED_BOT=0.0
      end if

      if(.not. allocated(C_RIVBED)) then
         allocate(C_RIVBED(NCOL,NROW,NLAY))
         C_RIVBED=0.0
      end if
    
      if(.not. allocated(TIME)) then
        allocate(TIME(NTSP))
        TIME=0.0
      end if

      if(.not. allocated(GW_HEAD)) then
         allocate(GW_HEAD(NCOL,NROW,NLAY,NTSP))
         GW_HEAD=0.0
      end if

 if(.not. allocated(recharge_flux)) then
         allocate(recharge_flux(NCOL,NROW,NLAY))
         recharge_flux=0.0
      end if

      if(.not. allocated(pumping_flux)) then
         allocate(pumping_flux(NCOL,NROW,NLAY))
         pumping_flux=0.0
      end if

      if(.not. allocated(river_flux)) then
         allocate(river_flux(NCOL,NROW,NLAY))
         river_flux=0.0
      end if


      if(.not. allocated(lat_xri_flux)) then
         allocate(lat_xri_flux(NCOL,NROW,NLAY))
         lat_xri_flux=0.0
      end if


      if(.not. allocated(lat_ycj_flux)) then
         allocate(lat_ycj_flux(NCOL,NROW,NLAY))
         lat_ycj_flux=0.0
      end if


      if(.not. allocated(vertical_flux)) then
         allocate(vertical_flux(NCOL,NROW,NLAY))
         vertical_flux=0.0
      end if

      if(.not. allocated(dStorage)) then
         allocate(dStorage(NCOL,NROW,NLAY))
         dStorage=0.0
      end if

      if(.not. allocated(Net_flux)) then
         allocate(Net_flux(NCOL,NROW,NLAY))
         Net_flux=0.0
      end if


      ! Now read each from a separate file:
      call read_K_file(K_FILE)
      call read_stor_file(STORAGE_FILE)
      call read_BC_file(BC_FILE)
      call read_source_file(SOURCE_FILE)
      call read_initial_head_file(HEAD0_FILE)
      
      ! Check if_there is a river?
      ! Call subroutine to read RIVER information 
      ! decide if there is a RIVER to read the information 

      do k = 1,NLAY
         do j=1, NROW
            do i=1, NCOL
               if (source_type(i,j,k)==SRC_RIVER)then
                  foundriver = .true.
                  exit
               end if
            end do
            if (foundriver) exit
         end do
         if (foundriver) exit
      end do

      if (foundriver)then
         call read_RIVER_file(RIVER_FILE)
      end if

      ! -----------------------------------------------------------------
      ! Call Surface elevation from Noah MP
      ! -----------------------------------------------------------------

      ! ----------------------------------------------------------------
      ! For testing : Flat surface only
      ! ---------------------------------------------------------------- 
      do j= 1, NROW
         do i= 1, NCOL
            TERRAIN (i,j)= 153.0
         end do
      end do
      
      ! --------------------------------------------------------------------
      ! Set HOLD before time loop to initial head
      ! ---------------------------------------------------------------------
      do k = 1, NLAY
         do j= 1, NROW
            do i= 1, NCOL
               HOLD(i,j,k)= H_initial (i,j,k)
            end do
         end do
      end do

      ! ----------------------------------------------------------------
      ! Speical case: if one of the two layers is inactive
      ! ---------------------------------------------------------------
      ! Check if the upper (unconfined) layer is inactive.
      if (UNCONF_THK == 0.0) then
         print*, "Upper layer inactive: constraining layer 1 to fixed initial heads."
         do j = 1, NROW
             do i = 1, NCOL
                ! Set the boundary condition for layer 1 (k=1) to fixed head.
                BC_type(i,j,1) = BC_FIXED
                BC_value(i,j,1) = H_initial(i,j,1)
             end do
         end do
      endif
       
      ! Check if the lower (confined) layer is inactive.
      if (CONF_AQ_THK == 0.0) then
          print*, "Lower layer inactive: constraining layer 2 to fixed initial heads."
          do j = 1, NROW
             do i = 1, NCOL
               ! Set the boundary condition for layer 2 (k=2) to fixed head.
               BC_type(i,j,2) = BC_FIXED
               BC_value(i,j,2) = H_initial(i,j,2)
               end do   
          end do
     end if

      print*, "GW_initialize done => N=", N
   end subroutine GW_initialize

! **********************************************************************  
! SUBROUTINE read_config
!  - read_config (GW_config.txt) => scalar domain/time
! **********************************************************************
   subroutine read_config()

      implicit none

      integer :: unit=99, ios
      character(len=256) :: line
      integer :: eqPos
      character(len=128) :: param, val

      open(unit=unit, file="GW_config.txt", status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error: Could not open GW_config.txt"
         stop
      end if

      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading line from config"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         eqPos=index(line,"=")
         if (eqPos==0) cycle

         param=trim(line(1:eqPos-1))
         val=trim(line(eqPos+1:))

         select case(param)
            case("NCOL");           read(val,*) NCOL
            case("NROW");           read(val,*) NROW
            case("NLAY");           read(val,*) NLAY
            case("DELTAC");         read(val,*) DELTAC
            case("DELTAR");         read(val,*) DELTAR
            case("NTSP");           read(val,*) NTSP
            case("DELTAT");         read(val,*) DELTAT
            case("MAX_ITER");       read(val,*) MAX_ITER
            case("Time_start");     read(val,*) Time_start
            case("picard_tol");     read(val,*) picard_tol            
            case("VAD_THK");        read(val,*) VAD_THK
            case("UNCONF_THK");     read(val,*) UNCONF_THK
            case("CONF_AQ_THK");    read(val,*) CONF_AQ_THK
            case("BED_THK");        read(val,*) BED_THK
            case("BED_K");          read(val,*) BED_K
            case("T_SOIL_THK");     read(val,*) T_SOIL_THK
            case("K_RIVERBED");     read(val,*) K_RIVERBED
            case("RIVERBED_THK");   read(val,*) RIVERBED_THK
            ! ---------------------------
            !   FILE NAMES
            ! ---------------------------
            case("K_FILE")
               K_FILE = val
            case("STORAGE_FILE")
               STORAGE_FILE = val
            case("BC_FILE")
               BC_FILE = val
            case("SOURCE_FILE")
               SOURCE_FILE = val
            case("HEAD0_FILE")
               HEAD0_FILE = val
            case("RIVER_FILE")
               RIVER_FILE = val
            case("WELL_FILE")
               WELL_FILE = val
            case default
               print*, "Ignoring unknown param =>", param
         end select
      end do

      close(unit)
   
      ! calculate N, the number of total grid cells
      N=NCOL*NROW*NLAY

      print*, "Read config => NCOL=",NCOL," NROW=",NROW," NLAY=",NLAY, "N", N
      print*, "DELTAC=",DELTAC," DELTAR=",DELTAR," NTSP=",NTSP," dt=",DELTAT
      print*, "MAX_ITER=",MAX_ITER," picard_tol=",picard_tol
      print*, "VAD_THK=",VAD_THK," UNCONF_THK=",UNCONF_THK, "CONF_AQ_THK", CONF_AQ_THK
      print*, "BED_THK=",BED_THK," BED_K=",BED_K, "T_SOIL_THK", T_SOIL_THK
      print*, "Riverbed_K", K_RIVERBED, "Riverbed thickness", RIVERBED_THK
   end subroutine read_config

!**********************************************************************
! SUBROUTINE read_K_file(filename)
!  file lines: i j k  KxVal KyVal KzVal
!  i in [1..NCOL], j in [1..NROW], k in [1..NLAY]
!**********************************************************************
   subroutine read_K_file(filename)

      implicit none

      character(len=*),intent(in):: filename
      integer :: unit=20, ios
      character(len=256):: line
      integer :: i,j,k
      real(kind=8) :: kxVal, kyVal, kzVal
      integer :: count

      print*, "Reading K from => ", filename
      open(unit=unit, file=filename, status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error open K file =>", filename
         stop
      end if

      count=0
      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading line in K file"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         read(line,*, iostat=ios) i, j, k, kxVal, kyVal, kzVal
         if (ios==0) then
            if (i>=1 .and. i<=NCOL .and. j>=1 .and. j<=NROW .and. k>=1 .and. k<=NLAY) then
               KR(i,j,k)= kxVal
               KC(i,j,k)= kyVal
               KV(i,j,k)= kzVal
               count=count+1
            else
               print*, "Out-of-range => (i,j,k)=",i,j,k
            end if
         else
            ios=0
            print*, "invalid line => ", line
         end if
      end do
      close(unit)
      print*, "Done reading K => assigned ", count," cells"
   end subroutine read_K_file

!**********************************************************************
! SUBROUTINE read_stor_file => lines: i j SyVal SsVal
!  but let's do Sy(i,j) to keep i= x, j= y. Actually let's keep the array
!  definition consistent with Kx which is Kx(i,j,k).
!  We'll store Sy(i,j), Ss(i,j).
!**********************************************************************
   subroutine read_stor_file(filename)

      implicit none

      character(len=*),intent(in):: filename
      integer :: unit=21, ios
      character(len=256):: line
      integer :: i,j
      real(kind=8) :: syVal, ssVal
      integer :: count

      print*, "Reading Storage from => ", filename
      open(unit=unit, file=filename, status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error open stor file =>", filename
         stop
      end if

      count=0
      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading line in stor file"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         read(line,*, iostat=ios) i, j, syVal, ssVal
         if (ios==0) then
            if (i>=1 .and. i<=NCOL .and. j>=1 .and. j<=NROW) then
               SY(i,j)= syVal
               SS(i,j)= ssVal
               count=count+1
            else
               print*, "Out-of-range stor => i=",i,"j=",j 
            end if
         else
            ios=0
            print*, "invalid stor line =>", line
         end if
      end do
      close(unit)
      print*, "Done reading stor => assigned ", count," cells"
   end subroutine read_stor_file

!**********************************************************************
! SUBROUTINE read_BC_file => i j k BCtype BCval
!**********************************************************************
   subroutine read_BC_file(filename)

      implicit none

      character(len=*),intent(in):: filename
      integer :: unit=22, ios
      character(len=256):: line
      integer :: i,j,k, bcT
      real(kind=8) :: bcV
      integer :: count

      print*, "Reading BC => ",filename
      open(unit=unit, file=filename, status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error open BC file =>",filename
         stop
      end if

      count=0
      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading BC line"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         read(line,*, iostat=ios) i,j,k, bcT, bcV
         if (ios==0) then
            if (i>=1 .and. i<=NCOL .and. j>=1 .and. j<=NROW .and. k>=1 .and. k<=NLAY) then
               BC_type(i,j,k)= bcT
               BC_value(i,j,k)= bcV
               count=count+1
            else
               print*, "Out-of-range BC =>", i,j,k
            end if
         else
            ios=0
            print*, "invalid BC line =>", line
         end if
      end do
      close(unit)
      print*, "Done reading BC => assigned ", count," lines."
   end subroutine read_BC_file

!**********************************************************************
! SUBROUTINE read_source_file => i j k srcType
!**********************************************************************
   subroutine read_source_file(filename)

      implicit none

      character(len=*),intent(in):: filename
      integer :: unit=23, ios
      character(len=256):: line
      integer :: i,j,k, sT
      !real(kind=8) :: sV
      integer :: count

      print*, "Reading Source => ",filename
      open(unit=unit, file=filename, status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error open source file =>",filename
         stop
      end if

      count=0
      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading source line"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         read(line,*, iostat=ios) i,j,k, sT !, sV
         if (ios==0) then
            if (i>=1 .and. i<=NCOL .and. j>=1 .and. j<=NROW .and. k>=1 .and. k<=NLAY) then
               source_type(i,j,k)= sT
               !source_value(i,j,k)= sV
               count=count+1
            else
               print*, "Out-of-range source =>", i,j,k
            end if
         else
            ios=0
            print*, "invalid source line =>", line
         end if
      end do
      close(unit)
      print*, "Done reading source => assigned ", count," lines."
   end subroutine read_source_file

!**********************************************************************
! SUBROUTINE read_RIVER_file => i j k  width length height bottomelev 
!**********************************************************************
   subroutine read_RIVER_file(filename)
   
      implicit none
           
      character(len=*),intent(in):: filename
      integer :: unit=26, ios
      character(len=256):: line
      integer :: i,j,k
      real(kind=8) :: width, length, height, bottomelev
      integer :: count



      print*, "Reading RIVER => ",filename
      open(unit=unit, file=filename, status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error open RIVER file =>",filename
         stop
      end if

      count=0
      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading RIVER line"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         ! Expect format: i  j  k  width  length  stage  botElev
         read(line,*, iostat=ios) i,j,k, width,length,height,bottomelev
         if (ios==0) then
            if (i>=1 .and. i<=NCOL .and. j>=1 .and. j<=NROW .and. k>=1 .and. k<=NLAY) then
               W_RIVER(i,j,k)= width
               L_RIVER(i,j,k)= length
               H_RIVER(i,j,k)= height
               RIVBED_BOT(i,j,k)= bottomelev
               count=count+1
            else
               print*, "Out-of-range RIVER =>", i,j,k
            end if
         else
            ios=0
            print*, "invalid RIVER line =>", line
         end if
      end do
      close(unit)
      print*, "Done reading RIVER  => assigned ", count," lines."
   end subroutine read_RIVER_file



!**********************************************************************
! SUBROUTINE read_initial_head_file => i j k h0_value
!**********************************************************************
   subroutine read_initial_head_file(fileName)

      implicit none

      character(len=*),intent(in):: filename
      integer :: unit=24, ios
      character(len=256):: line
      integer :: i,j,k
      real(kind=8) :: h0_value
      integer :: count

      print*, "Reading initial head => ",filename
      open(unit=unit, file=filename, status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error open initial head file =>",filename
         stop
      end if

      count=0
      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading initial head line"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         read(line,*, iostat=ios) i,j,k, h0_value
         if (ios==0) then
            if (i>=1 .and. i<=NCOL .and. j>=1 .and. j<=NROW .and. k>=1 .and. k<=NLAY) then
               
               H_initial(i,j,k)= h0_value
               count=count+1
            else
               print*, "Out-of-range initial head =>", i,j,k
            end if
         else
            ios=0
            print*, "invalid initial head line =>", line
         end if
      end do
      close(unit)
      print*, "Done reading initial head => assigned ", count," lines."      
   end subroutine read_initial_head_file

!**********************************************************************
! subroutine SawSim_solve
!   - create local PETSc objects
!   - do time loop (NTSP steps) with PCG + Picard
!**********************************************************************
   subroutine GW_solver()

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>   
      use petsc
      use petscvec
      use petscmat
      use petscksp
      use petscsys

      implicit none

      Mat               A_Matrix
      Vec               RHS_vector  
      Vec               h_old_vector
      Vec               h_new_vector
      KSP               ksp_GW
      PC                pc_GW
      PetscErrorCode    ierr  
      PetscInt          step, Iter
      PetscReal         diffNorm

      real(kind=8), allocatable :: HOLD_previous_timestep(:,:,:)

      ! Create matrix & vectors
      PetscCall(MatCreate(PETSC_COMM_WORLD, A_Matrix, ierr))
      PetscCall(MatSetSizes(A_Matrix, PETSC_DECIDE, PETSC_DECIDE, N, N, ierr))
      PetscCall(MatSetFromOptions(A_Matrix, ierr))
      PetscCall(MatSetUp(A_Matrix, ierr))

      PetscCall(VecCreate(PETSC_COMM_WORLD, RHS_vector, ierr))
      PetscCall(VecSetSizes(RHS_vector, PETSC_DECIDE, N, ierr))
      PetscCall(VecSetFromOptions(RHS_vector, ierr))

      PetscCall(VecDuplicate(RHS_vector, h_new_vector, ierr))
      PetscCall(VecDuplicate(RHS_vector, h_old_vector, ierr))

      ! Initialize the matrix and vector
      PetscCall(MatZeroEntries(A_Matrix,ierr))
      PetscCall(VecSet(RHS_vector,0.d0,ierr))
      PetscCall(VecSet(h_new_vector, 0.d0, ierr))
      PetscCall(VecSet(h_old_vector, 0.d0, ierr))
      
      ! Set up the linear , Krylov solver
      PetscCall(KSPCreate(PETSC_COMM_WORLD, ksp_GW, ierr))
      PetscCall(KSPSetOperators(ksp_GW, A_Matrix, A_Matrix, ierr))
      PetscCall(KSPGetPC(ksp_GW, pc_GW, ierr))
      ! Preconditionar type: Jacobi
      PetscCall(PCSetType(pc_GW, PCJACOBI, ierr))
      ! Krylov solver: GMRES, Note: CG does not work here 
      PetscCall(KSPSetType(ksp_GW, KSPGMRES, ierr))
      PetscCall(KSPSetFromOptions(ksp_GW, ierr))
      ! After creating the KSP object (ksp_GW) and setting the operators:
      PetscCall(KSPSetTolerances(ksp_GW, 1.0d-6, 1.0d-8, 1.0d4, 1000, ierr))

      !---------------------------------------------------
      ! Initialize heads: h_old_vector and h_new_vector
      !---------------------------------------------------
      
      call set_initial_heads(h_old_vector)
      call set_initial_heads(h_new_vector)
      
      ! allocate HOLD_previous_timestep
      allocate(HOLD_previous_timestep(NCOL,NROW,NLAY))
      !---------------------------------------------------
      ! Time loop with Picard iteration
      !---------------------------------------------------
      do step=1, NTSP
         print*, "Time step ", step, " of ", NTSP
         
         TIME(step) = Time_start + (step - 1) * DELTAT  ! in seconds
         
     

         ! Update the fortran HOLD
         call update_HOLD(h_old_vector)
         HOLD_previous_timestep(:,:,:)= HOLD(:,:,:)
         
         ! call head independednt subroutines
         call GW_recharge()
         
         !See if there WELL in the domain to call it
         call is_there_well()  

         ! build confined conductances because it doesnot depend of head
         call build_confined_conductances()

         ! build riverbed conductance based on river geomerty      
         if (any(source_type==SRC_RIVER)) then
            call build_river_conductance 
         end if
         
         do Iter=1, MAX_ITER

            ! get the latest “HOLD” from 
            call update_HOLD(h_old_vector)     
            ! Compute subsurface elevations only once per time step
            call calculate_subsurface_elevation()

            ! Compute conductances for unconfined only
            call build_unconfined_conductances()

            if (any(source_type==SRC_RIVER)) then
               call GW_RIVER_interaction()           ! head‐dependent P_RIVER / B_RIVER
            end if

            ! Build the matrix and RHS for the current time step
            call GW_AMatrix_RHS(step,h_old_vector, h_new_vector, A_Matrix, RHS_vector, ierr)

            ! View h_old_vector
            PetscCall(VecView(h_old_vector, PETSC_VIEWER_STDOUT_WORLD, ierr))
            
            ! Solve
            PetscCall(KSPSolve(ksp_GW, RHS_vector, h_new_vector, ierr))

            ! View h_new_vector
            PetscCall(VecView(h_new_vector, PETSC_VIEWER_STDOUT_WORLD, ierr))
            
            ! Compare
            call compare_heads(h_old_vector, h_new_vector, diffNorm, ierr)
            
            if (diffNorm< picard_tol) exit
            
            ! prepare for the next picard iteration
            PetscCall(VecCopy(h_new_vector, h_old_vector, ierr))

         end do
         
         
         ! write out h_new_vector to HNEW(:,:,:)
         call update_solution(h_new_vector, ierr)
      
         HOLD (:,:,:) = HOLD_previous_timestep(:,:,:)
         
         ! Save the resulrs to a CSV file
         call save_csv_results( 'GW_results.csv', step)
         
         ! Compute the water fluxes
         call compute_water_fluxes (step)
         
         ! Save the fluxes

         call save_fluxes_results ('flux_results.csv', step)

         PetscCall(VecCopy(h_new_vector, h_old_vector, ierr))
      end do

      ! print final
      !call update_solution(h_new_vector, ierr)
      deallocate(HOLD_previous_timestep)
      ! cleanup local
      PetscCall(KSPDestroy(ksp_GW, ierr))
      PetscCall(MatDestroy(A_Matrix, ierr))
      PetscCall(VecDestroy(RHS_vector, ierr))
      PetscCall(VecDestroy(h_new_vector, ierr))
      PetscCall(VecDestroy(h_old_vector, ierr))
   end subroutine GW_solver

!**********************************************************************
!  SUBROUTINE: compute_idx(i,j,k,iglobal,iPETSc)
!   to compute the index for the A_Matrix
!   Returns 1-based global index for default fortran
!   Returns 0-based global index for PETSc.
!**********************************************************************
   subroutine compute_idx(i, j, k,ipetsc)
   
      implicit none
   
      integer, intent(in)  :: i,j,k
      integer, intent(out) :: ipetsc
      integer              :: iglobal
      iglobal = i + NCOL*((j-1)+NROW*(k - 1))
      !ipetsc = (i-1) + NCOL*((j- 1)+NROW*(k - 1))
      ipetsc = iglobal-1
   end subroutine compute_idx

!**********************************************************************
!  SUBROUTINE: set_initial heads(h_old_vector)
!   to assign the initial heads to vector
!**********************************************************************
   subroutine set_initial_heads(head_vec)

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
      use petsc
      use petscvec
      use petscsys

      implicit none
      
      PetscErrorCode          ierr
      PetscInt                i,j,k,row
      PetscReal            :: val
   
      Vec, intent(inout)   :: head_vec
   
     ! PetscCall(VecView(head_vec, PETSC_VIEWER_STDOUT_WORLD, ierr))

      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               call compute_idx(i, j, k,row)
               !row=(i-1)+NCOL*((j-1)+NROW*(k-1))
               val= H_initial(i,j,k)
               PetscCall(VecSetValue(head_vec, row, val, INSERT_VALUES, ierr))
            end do
         end do
      end do
   
      PetscCall(VecAssemblyBegin(head_vec, ierr))
      PetscCall(VecAssemblyEnd(head_vec,   ierr))
    !  PetscCall(VecView(head_vec, PETSC_VIEWER_STDOUT_WORLD, ierr))
   end subroutine set_initial_heads

! **********************************************************************
!  SUBROUTINE: GW_is_there_well
!  Call subroutine to read pumping information 
! **********************************************************************
   subroutine is_there_well()
   
      implicit none
      !integer, intent(in) :: t_s
      integer :: i, j, k 
      logical :: foundwell

      foundwell = .false.
   
      ! decide if there is a WELL to read the information 
      do k = 1,NLAY
         do j=1, NROW
            do i=1, NCOL
               if (source_type(i,j,k)==SRC_WELL)then
                  foundwell = .true.
                  exit
               end if
            end do
            if (foundwell) exit
         end do
         if (foundwell) exit
      end do
      if (foundwell) then
         call read_WELL_file(WELL_FILE)
      end if

   end subroutine is_there_well


!**********************************************************************
! SUBROUTINE read_WELL_file => i j k pumpVal
!**********************************************************************
   subroutine read_WELL_file(filename)
      
      implicit none
   
      character(len=*),intent(in):: filename
      integer :: unit=25, ios
      character(len=256):: line
      integer :: t,i,j,k
      real(kind=8) :: pumpVal
      integer :: count

      print*, "Reading WELL => ",filename
      open(unit=unit, file=filename, status="old", action="read", iostat=ios)
      if (ios/=0) then
         print*, "Error open WELL file =>",filename
         stop
      end if

      count=0
      do
         read(unit,'(A)', iostat=ios) line
         if (ios<0) exit
         if (ios>0) then
            print*, "Error reading WELL line"
            stop
         end if
         line=trim(adjustl(line))
         if (len(line)==0 .or. line(1:1)=="#") cycle

         ! Expect format: t  i  j  k  pumpVal

         read(line,*, iostat=ios) t,i,j,k,pumpVal
         if (ios==0) then
            if (i>=1 .and. i<=NCOL .and. j>=1 .and. j<=NROW .and. k>=1 .and. k<=NLAY .and. t>=1 .and. t<=NTSP) then
               pumping_rate(i,j,k,t)= pumpVal /DELTAT ! devide my deltaT to make it m3/s 
               count=count+1
            else
               print*, "Out-of-range WELL =>", i,j,k, "time index=", t
            end if
         else
            ios=0
            print*, "invalid WELL line =>", line
         end if
      end do
      close(unit)
      print*, "Done reading WELL => assigned ", count," lines."
   end subroutine read_WELL_file



!**********************************************************************
! SUBROUTINE: update_HOLD
!  Here you update the HOLD in Fortran from h_old_vector in PETSc
! to use it in other subroutines
!**********************************************************************
   subroutine update_HOLD(hVec)

#include <petsc/finclude/petscvec.h>
       use petsc
       use petscvec

      implicit none

      Vec, intent(inout)  :: hVec
      PetscReal, pointer  :: arr(:)
      PetscErrorCode         ierr
      PetscInt               i,j,k, row

      ! Reset HOLD to 0.0 before updating it
      HOLD(:,:,:) = 0.0

      PetscCall(VecGetArrayF90(hVec, arr, ierr))
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               call compute_idx(i,j,k, row)
              if  (BC_type(i,j,k)==BC_NOFLOW .or. BC_type(i,j,k)==BC_FIXED)  then
                   HOLD(i,j,k) = BC_value(i,j,k)
               else
                   HOLD(i,j,k)=arr(row+1)
               end if
            end do
         end do
      end do
      PetscCall(VecRestoreArrayF90(hVec, arr, ierr))
   end subroutine update_HOLD


!**********************************************************************
! SUBROUTINE calculate_subsurface_elevation => fill GWTOP, GWBOT
!   with elevation based of surface elevation and soil thickness
!**********************************************************************

   subroutine calculate_subsurface_elevation () ! for testing (before integration with Noah MP)
   !subroutine calculate_subsurface_elevation (TERRAIN) for integration with Noah MP
      
      implicit none
      
      integer :: i,j,k
      if (.not. allocated(GWTOP))then    
         allocate (GWTOP(NCOL,NROW,NLAY))
         GWTOP = 0.0
      end if 
      if (.not. allocated(GWBOT))then
         allocate(GWBOT(NCOL,NROW,NLAY))       
         GWBOT = 0.0
      end if

      GWTOP(:,:,:) = 0.0
      GWBOT(:,:,:) = 0.0

      ! ----------------------------------------------------------------
      ! Calculate the elevation of the top and Bottom of each of GW layers 
      ! ---------------------------------------------------------------- 
      do k = 1, NLAY
         do j= 1, NROW
            do i= 1, NCOL
               ! Unconfined aquifer
               if (k==1)then
                  GWTOP(i,j,k) = HOLD(i,j,1)
                  GWBOT (i,j,k) = TERRAIN(i,j) - T_SOIL_THK - VAD_THK - UNCONF_THK
               else ! Confined aquifer
                  GWTOP (i,j,k) = GWBOT(i,j,k-1) - BED_THK
                  GWBOT (i,j,k) = GWTOP(i,j,k)- CONF_AQ_THK
               end if
            end do
         end do
      end do

   end subroutine calculate_subsurface_elevation

!**********************************************************************
! SUBROUTINE build_confined_conductances => fill TR, TC, DELTAV
!   with harmonic means of KR, CC, CV
!**********************************************************************
   subroutine build_confined_conductances()
      
      implicit none
      
      integer :: i,j,k
      real(kind=8),allocatable,dimension(:,:,:)  :: TR  ! Transmissivity along row x-direction (I)
      real(kind=8),allocatable,dimension(:,:,:) :: TC ! Transmissivity along column y direction (J)
     
      if (.not. allocated(TR)) then
         allocate(TR(NCOL,NROW,NLAY))
         TR = 0.0
      end if 
      if (.not. allocated(TC))then
         allocate(TC(NCOL,NROW,NLAY))
         TC =0.0
      end if

      TR(:,:,:) = 0.0
      TC(:,:,:) = 0.0
      DELTAV(:,:) = 0.0
      CC(:,:,:) = 0.0
      CR(:,:,:) = 0.0
      CV(:,:,:) = 0.0

      ! Loop to calculate transmissivity
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
                  !if (k==1)then ! Unconfined layer
                  !   ! calcuate the thickness of layer 1 (unconfined layer)
                  !   DELTAV(i,j)= HOLD(i,j,1)-GWBOT(i,j,1)
                  !   if (DELTAV(i,j) < 0.0)then
                  !       DELTAV(i,j)=0.0
                  !   end if
                  !   ! calculate transmissivity
                  !   TR(i,j,k) = KR(i,j,k)*DELTAV(i,j)
                  !   TC(i,j,k) = KC(i,j,k)*DELTAV(i,j)
                  if (k==2) then ! Confined layer
                     TR(i,j,k) = KR(i,j,k)*CONF_AQ_THK
                     TC(i,j,k) = KC(i,j,k)*CONF_AQ_THK
                  endif
            end do
         end do
      end do


      ! x direction => between i->i+1
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               if (i < NCOL) then
                  ! If cell i is at the left boundary and marked no-flow, then its right interface is blocked.
                  if (i == 1 .and. BC_type(i,j,k)==BC_NOFLOW) then
                     CR(i,j,k) = 0.0
                  ! If the cell i+1 is at the right boundary and marked no-flow, then its left interface is blocked.
                  else if (i+1 == NCOL .and. BC_type(i+1,j,k)==BC_NOFLOW) then
                     CR(i,j,k) = 0.0
                  else
                     ! calculate harmonic condutance
                     if (TR(i,j,k)==0.0)then
                        CR(i,j,k)=0.0
                     else    
                        CR(i,j,k)=2.0*DELTAC*(TR(i,j,k)*TR(i+1,j,k)/(TR(i,j,k)*DELTAR+TR(i+1,j,k)*DELTAR))
                     end if
                  end if
               end if
            end do
         end do
      end do


      ! y direction => between j->j+1
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               if (j < NROW) then
                  ! If cell j is at the upper boundary and marked no-flow, then its lower interface is blocked.
                  if (j == 1 .and. BC_type(i,j,k)==BC_NOFLOW) then
                     CC(i,j,k) = 0.0
                  ! If the cell j+1 is at the lower boundary and marked no-flow, then its upper interface is blocked.
                  else if (j+1 == NROW .and. BC_type(i,j+1,k)==BC_NOFLOW) then
                     CC(i,j,k) = 0.0
                  else
                     ! calculate harmonic condutance
                     if (TC(i,j,k)==0.0)then
                        CC(i,j,k)=0.0
                     else
                        CC(i,j,k)=2.0*DELTAR*(TC(i,j,k)*TC(i,j+1,k)/(TC(i,j,k)*DELTAC+TC(i,j+1,k)*DELTAC))
                     end if 
                 end if
               end if
            end do
         end do
      end do

      ! z direction => between k->k+1
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               if (k<NLAY) then
                  if (DELTAV(i,j)==0.0)then
                        CV(i,j,k)= 0.0
                  else
                       CV(i,j,k) =DELTAR*DELTAC/(0.5*(DELTAV(i,j)/KV(i,j,k))+(BED_THK/BED_K)+(0.5*CONF_AQ_THK/KV(i,j,k+1)))
                  end if
               end if
            end do
         end do
      end do

      deallocate (TC)
      deallocate (TR)

    !  print*, "confined conduction arrays done."
   end subroutine build_confined_conductances


!**********************************************************************
! SUBROUTINE build_unconfined_conductances => fill TR, TC, DELTAV
!   with harmonic means of KR, CC, CV
!**********************************************************************
   subroutine build_unconfined_conductances()
      
      implicit none
      
      integer :: i,j,k
      real(kind=8),allocatable,dimension(:,:,:)  :: TR  ! Transmissivity along row x-direction (I)
      real(kind=8),allocatable,dimension(:,:,:) :: TC ! Transmissivity along column y direction (J)
     
      if (.not. allocated(TR)) then
         allocate(TR(NCOL,NROW,NLAY))
         TR = 0.0
      end if 
      if (.not. allocated(TC))then
         allocate(TC(NCOL,NROW,NLAY))
         TC =0.0
      end if

      TR(:,:,:) = 0.0
      TC(:,:,:) = 0.0
      DELTAV(:,:) = 0.0
      CC(:,:,:) = 0.0
      CR(:,:,:) = 0.0
      CV(:,:,:) = 0.0

      ! Loop to calculate transmissivity
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
                  if (k==1)then ! Unconfined layer
                     ! calcuate the thickness of layer 1 (unconfined layer)
                     DELTAV(i,j)= HOLD(i,j,1)-GWBOT(i,j,1)
                     if (DELTAV(i,j) < 0.0)then
                         DELTAV(i,j)=0.0
                     end if
                     ! calculate transmissivity
                     TR(i,j,k) = KR(i,j,k)*DELTAV(i,j)
                     TC(i,j,k) = KC(i,j,k)*DELTAV(i,j)
                  !else ! Confined layer
                  !   TR(i,j,k) = KR(i,j,k)*CONF_AQ_THK
                  !   TC(i,j,k) = KC(i,j,k)*CONF_AQ_THK
                  endif
            end do
         end do
      end do


      ! x direction => between i->i+1
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               if (i < NCOL) then
                  ! If cell i is at the left boundary and marked no-flow, then its right interface is blocked.
                  if (i == 1 .and. BC_type(i,j,k)==BC_NOFLOW) then
                     CR(i,j,k) = 0.0
                  ! If the cell i+1 is at the right boundary and marked no-flow, then its left interface is blocked.
                  else if (i+1 == NCOL .and. BC_type(i+1,j,k)==BC_NOFLOW) then
                     CR(i,j,k) = 0.0
                  else
                     ! calculate harmonic condutance
                     if (TR(i,j,k)==0.0)then
                        CR(i,j,k)=0.0
                     else    
                        CR(i,j,k)=2.0*DELTAC*(TR(i,j,k)*TR(i+1,j,k)/(TR(i,j,k)*DELTAR+TR(i+1,j,k)*DELTAR))
                     end if
                  end if
               end if
            end do
         end do
      end do


      ! y direction => between j->j+1
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               if (j < NROW) then
                  ! If cell j is at the upper boundary and marked no-flow, then its lower interface is blocked.
                  if (j == 1 .and. BC_type(i,j,k)==BC_NOFLOW) then
                     CC(i,j,k) = 0.0
                  ! If the cell j+1 is at the lower boundary and marked no-flow, then its upper interface is blocked.
                  else if (j+1 == NROW .and. BC_type(i,j+1,k)==BC_NOFLOW) then
                     CC(i,j,k) = 0.0
                  else
                     ! calculate harmonic condutance
                     if (TC(i,j,k)==0.0)then
                        CC(i,j,k)=0.0
                     else
                        CC(i,j,k)=2.0*DELTAR*(TC(i,j,k)*TC(i,j+1,k)/(TC(i,j,k)*DELTAC+TC(i,j+1,k)*DELTAC))
                     end if 
                 end if
               end if
            end do
         end do
      end do

      ! z direction => between k->k+1
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               if (k<NLAY) then
                  if (DELTAV(i,j)==0.0)then
                        CV(i,j,k)= 0.0
                  else
                       CV(i,j,k) =DELTAR*DELTAC/(0.5*(DELTAV(i,j)/KV(i,j,k))+(BED_THK/BED_K)+(0.5*CONF_AQ_THK/KV(i,j,k+1)))
                  end if
               end if
            end do
         end do
      end do

      deallocate (TC)
      deallocate (TR)

    !  print*, "Unconfind conduction arrays done."
   end subroutine build_unconfined_conductances


!**********************************************************************
!  SUBROUTINE: assemble_matrix_and_rhs
!   - Builds matrix A and vector b(RHS) for the current time step and iteration
!   - Storage term: (Storage / dt)
!   - Conductances: CR, CC, CV
!    - Uses BC arrays (BC_type, BC_value) for Dirichlet/no-flow
!   - Uses storage approach: if k=1 => Sy, if k=2 => Ss, for transient
!   - Also includes any source (wells, rivers) from source arrays
!**********************************************************************
   subroutine GW_AMatrix_RHS(t_s,h_old_vector, h_new_vector,A_Matrix,RHS_vector,ierr)
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
      use petsc
      use petscvec
      use petscmat
   
      implicit none
      integer, intent(in)         :: t_s
      Mat, intent(inout)          :: A_Matrix
      Vec, intent(inout)          :: RHS_vector
      Vec, intent(inout)          :: h_old_vector
      Vec, intent(inout)          :: h_new_vector 
      PetscErrorCode              ierr
      PetscInt                    i,j,k, row, col
      PetscReal, pointer          :: arr_old(:), arr_new(:)
      PetscScalar                 diagVal, off_diaVal,stor_Val, oldHead 
   
      real(kind=8),allocatable,dimension(:,:,:) :: HCOF ! Head Coefficent
      real(kind=8),allocatable,dimension(:,:,:) :: B !The constant cofficient from all exteral sources
      real(kind=8),allocatable,dimension(:,:,:) :: P !Head dependent cofficient from all external sources 
   
      if (.not. allocated (HCOF))then
         allocate(HCOF(NCOL,NROW,NLAY))
         HCOF = 0.0
      end if
      if (.not. allocated(B))then
         allocate (B(NCOL,NROW,NLAY))
         B = 0.0
      end if
      if (.not. allocated(P))then
         allocate(P(NCOL,NROW,NLAY))
         P = 0.0
      end if

      HCOF(:,:,:) = 0.0
      B(:,:,:) = 0.0
      P(:,:,:) = 0.0


      
      ! Access the Fortran arrays of the PETSc vector
      PetscCall(VecGetArrayF90(h_old_vector, arr_old, ierr))
      PetscCall(VecGetArrayF90(h_new_vector, arr_new, ierr))
   
      ! Clear the matrix and vector
      PetscCall(MatZeroEntries(A_Matrix,ierr))
      PetscCall(VecSet(RHS_vector,0.d0,ierr))

      ! ------------------------------------------------------------------------------------
      ! Calculate P(i,j,k) and B(i,j,k)
      !   There are 3 external sources 
      !    1) Recharge: Always there and it does not depend on the current GW head
      !    2) Pumpimg(Well): (-ive) for extraction and (+ive) for injection
      !       it does not depends on GW head
      !    3) River: it depends on the current GW head , it is interaction relationship 
      !       that depends on the water head difference.          
      !       river-GW interconnection (Q_river)
      !       Q_river= C_river * (H_river- H_GW)
      !       Q_river=C_river*H_river-C_river*H_GW   
      !       Q_river=B+P*H_GW, where B= C_river*H_river; and P = -C_river  
      ! -----------------------------------------------------------------------------------
      
!      ! Call recharge from NoahMP
!      call GW_recharge()

!      ! Call GW_river_exchange if there is a river 
!      do k = 1, NLAY
!         do j = 1, NROW
!             do i = 1, NCOL
!               if (source_type(i,j,k) ==SRC_RIVER)then
!                  call GW_RIVER_interaction()
!               endif
!            end do
!         end do
!      end do

      do k = 1, NLAY
         do j = 1, NROW
             do i = 1, NCOL
               if (source_type(i,j,k) ==SRC_RIVER)then
                  P(i, j, k) = P_RIVER (i,j,k)
                  B(i, j, k) = recharge_rate(i,j)+B_RIVER(i,j,k)
               elseif (source_type(i,j,k) ==SRC_WELL)then
                  P(i, j, k) = 0.0
                  B(i, j, k) = recharge_rate(i,j)+ pumping_rate(i,j,k,t_s)
               else ! source_type(i,j,k)==SCR_NONE (remember there is always recharge)
                  P(i, j, k) = 0.0
                  B(i, j, k) = recharge_rate(i,j)
               endif
            end do
         end do
      end do

      do k = 1, NLAY
         do j = 1, NROW
             do i = 1, NCOL
               ! Calculate HCOF(i,j,k)
               if (k == 1) then
                  HCOF(i, j, k) = P(i, j, k) - (SY(i, j) * DELTAC *DELTAR / DELTAT)
               else
                  HCOF(i, j, k) = P(i, j, k) - (SS(i, j) * DELTAC *DELTAR * CONF_AQ_THK / DELTAT)
               end if
            end do
         end do
      end do
   
      ! assemble the matrix   
      do k=1, NLAY
         do j=1, NROW
            do i=1, NCOL
               
               !Set the intital values to zero
               diagVal= 0.0
               off_diaVal = 0.0
               stor_Val = 0.0

               call compute_idx(i,j,k,row)
            
               ! Dirichlet
               if (BC_type(i,j,k) == BC_FIXED) then
                  PetscCall(MatSetValue(A_Matrix, row, row, 1.d0, INSERT_VALUES, ierr))
                  PetscCall(VecSetValue(RHS_vector, row, BC_value(i,j,k), INSERT_VALUES, ierr))
                  cycle
               else if (BC_type(i,j,k) == BC_NOFLOW) then
                  ! This means treat no flow from the outside boundary
                  ! We'll handle that by not connecting out-of-domain neighbors
                  ! (just normal logic: if i>1 => connect left, etc.)
                  ! so do nothing special except handle interior flows below
               end if

               
               !-------------------------------------------------------------
               !   Add_flow_x
               !   Adds +/- x conduction from cell (i,j,k) to matrix
               !---------------------------------------------------------------
            
               ! +x neighbor
               if (i < NCOL) then
                  off_diaVal = CR(i,j,k)
                  diagVal= diagVal+ off_diaVal
                  call compute_idx(i+1,j,k, col)
                  PetscCall(MatSetValue(A_Matrix, row, col, -off_diaVal, INSERT_VALUES, ierr))
               end if
            
               ! -x neighbor
               if (i > 1) then
                  off_diaVal = CR(i-1,j,k)
                  diagVal = diagVal+ off_diaVal
                  call compute_idx(i-1,j,k, col)
                  PetscCall(MatSetValue(A_Matrix, row, col, -off_diaVal, INSERT_VALUES, ierr))
               end if
               
               !----------------------------------------------------------------------
               ! Add_flow_y
               !   Adds +/- y conduction from cell (i,j,k) to matrix
               !----------------------------------------------------------------------
               ! +y
               if (j < NROW) then
                  off_diaVal = CC(i,j,k)
                  diagVal = diagVal+ off_diaVal
                  call compute_idx(i, j+1, k, col)
                  PetscCall (MatSetValue(A_Matrix, row, col, -off_diaVal, INSERT_VALUES, ierr))
               end if
            
               ! -y
               if (j > 1) then
                  off_diaVal = CC(i, j-1, k)
                  diagVal = diagVal+ off_diaVal
                  call compute_idx(i, j-1, k, col)
                  PetscCall(MatSetValue(A_Matrix, row, col, -off_diaVal, INSERT_VALUES, ierr))
               end if

               !----------------------------------------------------------------------
               ! Add_flow_z
               !   Adds +/- z conduction from cell (i,j,k) to matrix
               !----------------------------------------------------------------------
            
               ! +z
               if (k < NLAY) then
                  off_diaVal = CV(i,j,k)
                  diagVal = diagVal+ off_diaVal
                  call compute_idx(i, j, k+1, col)
                  PetscCall(MatSetValue(A_Matrix, row, col, -off_diaVal, INSERT_VALUES, ierr))
               end if   
               ! -z
               if (k > 1) then
                  off_diaVal = CV(i,j,k-1)
                  diagVal = diagVal+ off_diaVal
                  call compute_idx(i, j, k-1, col)
                  PetscCall(MatSetValue(A_Matrix, row, col, -off_diaVal, INSERT_VALUES, ierr))
               end if

               ! finalize diagonal
               PetscCall(MatSetValue(A_Matrix, row, row, diagVal-HCOF(i,j,k), INSERT_VALUES, ierr))
                

               ! Add storage term to the RHS
               ! Storage: if k=1 => use Sy(i,j), if k=2 => Ss(i,j)
               if (k == 1) then
                  stor_Val = (SY(i, j)*DELTAC*DELTAR/DELTAT)
               else
                  stor_Val = (SS(i,j)*DELTAC*DELTAR*CONF_AQ_THK/DELTAT)
               end if
               oldHead = 0.0 
               oldHead = arr_old(row)
               !call VecSetValue(RHS_vector, row, storVal*oldHead, ADD_VALUES, iErr)
               PetscCall(VecSetValue(RHS_vector, row, B(i,j,k)+ stor_Val*oldHead, INSERT_VALUES, ierr))                           
            end do
         end do
      end do

      PetscCall(MatAssemblyBegin(A_Matrix, MAT_FINAL_ASSEMBLY, ierr))
      PetscCall(MatAssemblyEnd(A_Matrix,   MAT_FINAL_ASSEMBLY, ierr))
      PetscCall(VecAssemblyBegin(RHS_vector, ierr))
      PetscCall(VecAssemblyEnd(RHS_vector,   ierr))

      !PetscCall(MatView(A_Matrix, PETSC_VIEWER_STDOUT_WORLD, ierr))
      !PetscCall(VecView(RHS_vector, PETSC_VIEWER_STDOUT_WORLD, ierr))
      PetscCall(VecRestoreArrayF90(h_old_vector, arr_old, ierr))
      PetscCall(VecRestoreArrayF90(h_new_vector, arr_new, ierr))

      ! deallocate 
      deallocate (HCOF)
      deallocate (B)
      deallocate (P)
   end subroutine GW_AMatrix_RHS




!**********************************************************************
!  SUBROUTINE: GW_recharge (RECH,DEEPRECH)
!  Recharge to GW based on the unconfined GW head
!  RECH and DEEPRECH have unit of : Length/Time , needs to convert to L3/m (*DELTAC*DELTAR) 
!**********************************************************************

   !subroutine GW_recharge (RECH,DEEPRECH,WTD) ! When coupling with NoahMP
   subroutine GW_recharge()
   
      implicit none
      
      integer :: i, j, k 
      
      !real(kind=8), intent (in),dimension(:,:)         :: RECH
      !real(kind=8), intent (in),dimension(:,:)         :: DEEPRECH
      !real(kind=8), intent (inout), dimension(:,:)       :: WTD
      if (.not. allocated(WTD))then
        allocate (WTD(NCOL,NROW))
        WTD = 0.0
      end if

      if (.not. allocated(recharge_rate))then
        allocate(recharge_rate(NCOL,NROW))
        recharge_rate = 0.0
      end if

      WTD(:,:) = 0.0
   
   ! --------------------When Coupling with Noah MP--------------------------
   !   do j = 1, NROW
   !      do i = 1, NCOL
   !         WTD (i,j) = TERRAIN(i,j)-HOLD(i,j,1)
   !         if (WTD (i,j) < T_SOIL_THK) then!
   !            ! the WT is withon soil layers 
   !            call SHALLOWWATERTABLE(WTD) ! From NOAH MP
   !            recharge_rate(i,j) = RECH(i,j)*DELTAC*DELTAR ! RECH from SHALLOWWATERTABLE subroutine Noah MP (SS: double check)
   !         else 
   !            ! the WT is below soil layers 
   !            call WTABLE_mmf_noahmp(WTD) ! From Noah MP
   !            recharge_rate (i,j) = DEEPRECH(i,j)*DELTAC*DELTAR
   !         end if 
   !     end do
   ! end do
   ! ------------------------ For Testing purposes ----------------------------
   
      do j = 1, NROW
         do i = 1, NCOL
            if (HOLD(i,j,1)==0.0)then
                WTD(i,j)=0.0
            else
                WTD (i,j) = TERRAIN(i,j)-HOLD(i,j,1)
            end if
            recharge_rate(i,j) = 8.0e-9 *DELTAC*DELTAR
         end do
      end do
   end subroutine GW_recharge


! **********************************************************************
!  SUBROUTINE: build river conductance
!  - Declare the River gemotry,head, and riverbed_elevation
!  - Read River input file
!  - Calculate Riverbed conductance 
! **********************************************************************
   subroutine build_river_conductance()
      
      implicit none
          
      integer :: i,j,k
      
      if (.not. allocated(C_RIVBED))then
         allocate(C_RIVBED(NCOL,NROW,NLAY))
         C_RIVBED = 0.0
      end if
      ! Calculate River bed conductance 
      do k = 1,NLAY
         do j=1, NROW
            do i=1, NCOL
            C_RIVBED(i,j,k) = K_RIVERBED *L_RIVER(i,j,k)*W_RIVER(i,j,k)/RIVERBED_THK
            end do
         end do
      end do      
   end subroutine build_river_conductance

! **********************************************************************
!  SUBROUTINE: GW_RIVER_interaction
!  - Declare the River gemotry,head, and riverbed_elevation
!  - Read River input file
!  - Calculate Riverbed conductance 
! **********************************************************************
   subroutine GW_RIVER_interaction()
      
      implicit none
          
      
      integer :: i,j,k
      
      if (.not. allocated(C_RIVBED))then
         allocate(C_RIVBED(NCOL,NROW,NLAY))
         C_RIVBED = 0.0
      end if
      
      if (.not. allocated(P_RIVER))then
         allocate(P_RIVER(NCOL,NROW,NLAY))
         P_RIVER = 0.0
      end if
      
      if (.not. allocated(B_RIVER))then
         allocate(B_RIVER(NCOL,NROW,NLAY))
         B_RIVER = 0.0
      end if
      
      ! Does the river recieve water from GW OR leak water to GW
      do k = 1,NLAY
         do j=1, NROW
            do i=1, NCOL
               if (HOLD(i,j,k) > RIVBED_BOT(i,j,k)) then
                  P_RIVER(i,j,k) = -C_RIVBED(i,j,k)
                  B_RIVER(i,j,k) = C_RIVBED(i,j,k)*H_RIVER(i,j,k)
               else
                  P_RIVER(i,j,k) = 0.0     
                  B_RIVER (i,j,k) = C_RIVBED(i,j,k)*(H_RIVER(i,j,k)-RIVBED_BOT(i,j,k))
               end if
            end do
         end do
      end do
      
   end subroutine GW_RIVER_interaction


!**********************************************************************
!  SUBROUTINE: compare_heads
!   to compute L-infinity norm of difference (h_new - h_old)
!**********************************************************************
   subroutine compare_heads(h_old, h_new, diffNorm,ierr)

#include <petsc/finclude/petscvec.h>
      use petscvec

      implicit none 

      Vec, intent(in)         ::  h_old, h_new
      PetscReal, intent(out)  ::  diffNorm
      Vec                         diffVec
      PetscErrorCode              ierr
      
      diffNorm = 0.0

      PetscCall(VecDuplicate(h_new, diffVec, ierr))
      PetscCall(VecCopy(h_new, diffVec, ierr))              ! diffVec = h_new
      PetscCall(VecAXPY(diffVec, -1.d0, h_old,ierr))    ! diffVec = h_new - h_old
      PetscCall(VecNorm(diffVec, NORM_INFINITY, diffNorm, ierr))
      PetscCall(VecDestroy(diffVec, ierr))
   end subroutine compare_heads


!**********************************************************************
!  SUBROUTINE: compute water fluxes
!  If the neighboring head is higher, Q_>0 → positive inflow into cell.
!  If it’s lower, Q<0 → negative → outflow from cell.
!  i increases from west→east
!  j increases (moves) from north→south
!  k increases (moves) from up→down
!  head_diff = HNEW(neighbor) - HNEW(i,j,k)
!  Net_flow = sum_of_all_Q_face + external_sources
!**********************************************************************
   subroutine compute_water_fluxes (t_s)
      implicit none

      integer, intent(in)         :: t_s
      integer        ::    i,j,k
      real(kind=8)   ::    head_diff_x, head_diff_y, head_diff_z
      real(kind=8),allocatable,dimension(:,:,:)  :: Q_lat_xri_west  ! ! flow across the west face, from (i−1,j,k) → (i,j,k) in m3/s
      real(kind=8),allocatable,dimension(:,:,:) :: Q_lat_xri_east ! ! flow across the east face, from (i,j,k) → (i+1,j,k)   in m3/s
      real(kind=8),allocatable,dimension(:,:,:) :: Q_lat_ycj_north ! ! flow across the north face, from (i,j-1,k) → (i,j,k)   in m3/s
      real(kind=8),allocatable,dimension(:,:,:) :: Q_lat_ycj_south! ! flow across the south face, from (i,j,k) → (i,j+1,k)   in m3/s
      real(kind=8),allocatable,dimension(:,:,:) :: Q_vert_up ! ! flow across the up face, from (i,j,k-1) → (i,j,k)   in m3/s
      real(kind=8),allocatable,dimension(:,:,:) :: Q_vert_down ! ! flow across the down face, from (i,j,k) → (i,j,k+1)   in m3/s
      
      real(kind=8) ::imbalance

      if (.not. allocated(Q_lat_xri_east)) then
         allocate(Q_lat_xri_east(NCOL,NROW,NLAY))
         Q_lat_xri_east = 0.0
      end if 
      
      if (.not. allocated(Q_lat_xri_west))then
         allocate(Q_lat_xri_west(NCOL,NROW,NLAY))
         Q_lat_xri_west =0.0
      end if

      if (.not. allocated(Q_lat_ycj_north)) then
         allocate(Q_lat_ycj_north(NCOL,NROW,NLAY))
         Q_lat_ycj_north = 0.0
      end if 

      if (.not. allocated(Q_lat_ycj_south))then
         allocate(Q_lat_ycj_south(NCOL,NROW,NLAY))
         Q_lat_ycj_south =0.0
      end if

      if (.not. allocated(Q_vert_down)) then
         allocate(Q_vert_down(NCOL,NROW,NLAY))
         Q_vert_down = 0.0
      end if

      if (.not. allocated(Q_vert_up))then
         allocate(Q_vert_up(NCOL,NROW,NLAY))
         Q_vert_up =0.0
      end if


      ! zero out your module arrays (allocated in GW_initialize)
      recharge_flux = 0.0
      pumping_flux  = 0.0
      river_flux    = 0.0
      lat_xri_flux  = 0.0
      lat_ycj_flux  = 0.0
      vertical_flux = 0.0
      dStorage      = 0.0
      Net_flux      = 0.0

      !VERY IMP: (delta T in second to covert to the unit of timestep)    
      do k = 1, NLAY
         do j = 1, NROW
            do i = 1, NCOL

               ! ---- External fluxes -----
               ! Q_recharge comes from GW_recharge already and it applies on Layer 1 only 
               if (k == 1) then
                  !Q_recharge(i,j,k) = recharge_rate(i,j) ! in m3/s
                  recharge_flux (i,j,k) = recharge_rate (i,j)/(DELTAC *DELTAR)*DELTAT ! in m
               else
                  !Q_recharge(i,j,k) = 0.0
                  recharge_flux (i,j,k)= 0.0
               end if

               ! Check Units here, The unit of pumping is important
               ! Q_pumping comes from well data
               if (source_type(i,j,k) == SRC_WELL) then
                  !Q_pumping(i,j,k) = pumping_rate(i,j,k, t_s)  ! NEED TO MAKE SURE THAT THE UNIT in m3/s (per second)
                  pumping_flux(i,j,k) = pumping_rate(i,j,k,t_s)/ (DELTAC * DELTAR) *DELTAT ! NEED TO MAKE SURE THAT THE UNIT in m/timestep (The unit of time step)
               else
                  !Q_pumping(i,j,k) = 0.0
                  pumping_flux(i,j,k) = 0.0 
               end if
               
               ! Q_river comes from GW_river_interaction
               if (source_type(i,j,k) == SRC_RIVER) then
                  !Q_river(i,j,k) = B_RIVER(i,j,k) + P_RIVER(i,j,k) * HNEW(i,j,k)  ! in m3/s
                  river_flux(i,j,k) = (B_RIVER(i,j,k) + P_RIVER(i,j,k) * HNEW(i,j,k)) / (DELTAC * DELTAR)*DELTAT  ! in m
               else
                  !Q_river(i,j,k) = 0.0
                  river_flux(i,j,k) = 0.0 
               end if

               ! ---- Lateral fluxes -----

               ! X-direction 
               ! flow across the west face, from (i−1,j) → (i,j)
               !––––– west face (neighbor i-1,j) –––––
               if (i>1) then
                  head_diff_x= HNEW(i,j,k) - HNEW(i-1,j,k)
                  Q_lat_xri_west(i,j,k) = - CR(i-1,j,k)*head_diff_x
               else
                  Q_lat_xri_west(i,j,k) = 0.0
               end if

               ! flow across the east face, from (i,j) → (i+1,j)
               !––––– east face (neighbor i+1,j) –––––
               if (i < NCOL) then
                  head_diff_x = HNEW(i+1,j,k) - HNEW(i,j,k) ! or HOLD ?
                  Q_lat_xri_east(i,j,k) = - CR(i,j,k) * head_diff_x
               else
                  Q_lat_xri_east(i,j,k) = 0.0
               end if



               ! Y-direction 
               ! flow across the north face, from (i,j-1) → (i,j)
               !––––– north face (neighbor i,j-1) –––––
               if (j>1) then
                  head_diff_y = HNEW(i,j,k)-HNEW(i,j-1,k)
                  Q_lat_ycj_north(i,j,k) = - CC(i,j-1,k)*head_diff_y
               else
                  Q_lat_ycj_north(i,j,k) = 0.0
               end if

               ! flow across the south face, from (i,j) → (i,j+1)
               !––––– south face (neighbor i,j+1) –––––
               if (j < NROW) then
                  head_diff_y = HNEW(i,j+1,k) - HNEW(i,j,k) ! or HOLD ?
                  Q_lat_ycj_south(i,j,k) =  - CC(i,j,k) * head_diff_y
               else
                  Q_lat_ycj_south(i,j,k) = 0.0
               end if

               ! ---- Vertical Flux  Z- direction ----
               ! flow across the up face, from (i,j,k-1) → (i,j,k)
               !––––– up face (neighbor i,j,k-1) –––––
               if (k > 1) then
                  head_diff_z = HNEW(i,j,k) - HNEW(i,j,k-1) ! ! or HOLD ?
                  Q_vert_up(i,j,k) =  - CV(i,j,k-1) * head_diff_z
               else
                  Q_vert_up(i,j,k) = 0.0
               end if

                ! flow across the down face, from (i,j,k) → (i,j,k+1)  
               !––––– down face (neighbor i,j,k+1) –––––
               if (k < NLAY) then
                  head_diff_z = HNEW(i,j,k+1) - HNEW(i,j,k) ! ! or HOLD ?
                  Q_vert_down(i,j,k) = - CV(i,j,k) * head_diff_z
               else
                  Q_vert_down(i,j,k) = 0.0
               end if


               ! ---- Storage change ----
               ! dStorage = storage_coefficient * (H_new - H_old)
               if (k == 1) then
                  !dStorage_m3(i,j,k) = (HNEW(i,j,k) - HOLD(i,j,k)) * (SY(i,j) * DELTAC*DELTAR) ! m3
                  dStorage(i,j,k)    = (HNEW(i,j,k) - HOLD(i,j,k)) * SY(i,j) ! m
               else
                  !dStorage_m3(i,j,k) = (HNEW(i,j,k) - HOLD(i,j,k)) * (SS(i,j) * DELTAC*DELTAR * CONF_AQ_THK) !in m3
                  dStorage(i,j,k)    = (HNEW(i,j,k) - HOLD(i,j,k)) * (SS(i,j) * CONF_AQ_THK) !in m
               end if

               ! ---- Net flux ----
               ! Inflow - Outflow 
               ! west (>) In
               ! east (>) out 
               !Q_lat_xri (i,j,k) = Q_lat_xri_west (i,j,k) - Q_lat_xri_east (i,j,k) ! in m3/s
               !Q_lat_ycj (i,j,k) = Q_lat_ycj_north (i,j,k) - Q_lat_ycj_south (i,j,k) ! in m3/s
               !Q_vert (i,j,k) = Q_vert_up (i,j,k) - Q_vert_down (i,j,k) ! in m3/s

               !VERY IMP: (delta T in second to covert to the unit of timestep)    
               lat_xri_flux (i,j,k) = (Q_lat_xri_west (i,j,k) - Q_lat_xri_east (i,j,k)) / (DELTAC * DELTAR)*DELTAT  ! in m
               lat_ycj_flux (i,j,k) = (Q_lat_ycj_north (i,j,k) - Q_lat_ycj_south (i,j,k)) / (DELTAC * DELTAR)*DELTAT  ! in m
               vertical_flux (i,j,k) = (Q_vert_up (i,j,k) - Q_vert_down (i,j,k)) / (DELTAC * DELTAR)*DELTAT  ! in m

               !Net_flow(i,j,k) = Q_recharge(i,j,k) + Q_pumping(i,j,k) + Q_river(i,j,k) + Q_lat_xri(i,j,k) + Q_lat_ycj(i,j,k)+ Q_vert(i,j,k) ! in m3/s
               Net_flux (i,j,k)= recharge_flux(i,j,k) + pumping_flux(i,j,k)+river_flux(i,j,k) +lat_xri_flux(i,j,k)+lat_ycj_flux(i,j,k)+vertical_flux(i,j,k) ! in m
               ! The water budget in each cell should satisfy:
               ! Net_flow(i,j,k)*DELTAT = dStorage_m3(i,j,k)   (within numerical error) m3

               imbalance= Net_flux(i,j,k)-dStorage(i,j,k)

               if (i==10 .and. j==9 .and. k==2) then
                  write(*,'(A)') '--- BUDGET BREAKDOWN for cell (10,9,2) ---'
                  write(*,'(A,F12.6)') '  recharge_flux = ', recharge_flux(i,j,k)
                  write(*,'(A,F12.6)') '  pumping_flux = ',  pumping_flux(i,j,k)
                  write(*,'(A,F12.6)') '  river_flux = ',    river_flux(i,j,k)
                  write(*,'(A,F12.6)') '  x-flux west  = ', Q_lat_xri_west(i,j,k)/(DELTAC*DELTAR)*DELTAT
                  write(*,'(A,F12.6)') '  x-flux east  = ', Q_lat_xri_east(i,j,k)/(DELTAC*DELTAR)*DELTAT
                  write(*,'(A,F12.6)') '  y-flux north = ', Q_lat_ycj_north(i,j,k)/(DELTAC*DELTAR)*DELTAT
                  write(*,'(A,F12.6)') '  y-flux south = ', Q_lat_ycj_south(i,j,k)/(DELTAC*DELTAR)*DELTAT
                  write(*,'(A,F12.6)') '  z-flux up    = ', Q_vert_up(i,j,k)/(DELTAC*DELTAR)*DELTAT
                  write(*,'(A,F12.6)') '  z-flux down  = ', Q_vert_down(i,j,k)/(DELTAC*DELTAR)*DELTAT
                  write(*,'(A,F12.6)') '  Net_flux     = ', Net_flux(i,j,k)
                  write(*,'(A,F12.6)') '  dStorage     = ', dStorage(i,j,k)
                  write(*,'(A,F12.6)') '  imbalance    = ', imbalance
                  write(*,'(A)') '-------------------------------------------'
               end if

            end do
         end do
      end do

      ! from flow m3/s to flux m/timstep

 
      deallocate (Q_lat_xri_east)
      deallocate (Q_lat_xri_west)
      deallocate (Q_lat_ycj_north)
      deallocate (Q_lat_ycj_south)
      deallocate (Q_vert_down)
      deallocate (Q_vert_up)

   end subroutine compute_water_fluxes

!**********************************************************************
! SUBROUTINE: update_solution
!   Utility to update HNEW
!**********************************************************************
   subroutine update_solution(hVec,ierr)

#include <petsc/finclude/petscvec.h>
      use petscvec
      use petsc
      implicit none

      Vec, intent(in)         :: hVec
      PetscReal, pointer      :: arr(:)
      PetscErrorCode            ierr
      PetscInt                  i,j,k, row
      !real(kind=8)            :: h_final


      PetscCall(VecGetArrayF90(hVec, arr, ierr))
      do k=1,NLAY
         do j=1,NROW
            do i=1,NCOL
               call compute_idx(i,j,k, row)
               HNEW(i,j,k)= arr(row+1)
               if (BC_type(i,j,k)==BC_NOFLOW .or. BC_type(i,j,k)==BC_FIXED) then
                  HNEW(i,j,k)= BC_value(i,j,k)
               end if                
               write(*,'("(i,j,k)=(",3i2,") => h=",f8.3)') i,j,k, HNEW(i,j,k)
            end do
         end do
      end do
      PetscCall(VecRestoreArrayF90(hVec, arr, ierr))
   end subroutine update_solution


!**********************************************************************
!  SUBROUTINE: save_csv_results
!  Save the results as CSV file       
!**********************************************************************

   subroutine save_csv_results(filename, time_step)
      use petsc
      use petscvec
      implicit none
    
      ! Input arguments

      character(len=*), intent(in) :: filename
      integer, intent(in) :: time_step

      PetscErrorCode :: ierr
      PetscInt :: i, j, k, row
      integer :: unit, ios
      logical :: file_exists
      real(kind=8)            :: h_final

      ! Check if the file already exists. If not, create and write header.
      inquire(file=filename, exist=file_exists)
      if (.not. file_exists) then
          unit = 99
          open(unit=unit, file=filename, status="replace", action="write", form="formatted", iostat=ios)
          if (ios /= 0) then
             print*, "Error opening file ", filename
             stop
          end if
          ! Write a header line (time, i, j, k, head)
          write(unit,'("time,i,j,k,head")')
      else
          unit = 99
          open(unit=unit, file=filename, status="old", action="write", &
               position="append", form="formatted", iostat=ios)
          if (ios /= 0) then
             print*, "Error opening file ", filename
             stop
          end if
      end if    
      do k = 1, NLAY
         do j = 1, NROW
            do i = 1, NCOL
               h_final= HNEW(i,j,k)
               if (BC_type(i,j,k)==BC_NOFLOW .or. BC_type(i,j,k)==BC_FIXED) then
                   h_final = BC_value(i,j,k)
               end if     
               ! Write a line: time_step, i, j, k, head
               write(unit,'(F12.6,1X,I3,1X,I3,1X,I3,1X,F12.6)') time_step, i, j, k, h_final
            end do
         end do
      end do
      close(unit)
   end subroutine save_csv_results

!**********************************************************************
!  SUBROUTINE: save_fluxes
!***********************************************************************
   subroutine save_fluxes_results(filename, time_step)
      
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: time_step
      integer :: i, j, k, unit, ios
      logical :: file_exists

      ! Check if file exists; if not, write header.
      inquire(file=filename, exist=file_exists)
      if (.not. file_exists) then
         unit = 101
         open(unit=unit, file=filename, status="replace", action="write", form="formatted", iostat=ios)
         write(unit,*) 'time,i,j,k,recharge_flux,pumping_flux,river_flux,lat_xri_flux,lat_ycj_flux,vertical_flux,dStorage,Net_flux'
      else
         unit = 101
         open(unit=unit, file=filename, status="old", action="write", position="append", form="formatted", iostat=ios)
      end if

      do k = 1, NLAY
         do j = 1, NROW
            do i = 1, NCOL
               write(unit,*) time_step, i, j, k, recharge_flux(i,j,k), pumping_flux(i,j,k),river_flux(i,j,k), &
                             lat_xri_flux(i,j,k), lat_ycj_flux(i,j,k), vertical_flux(i,j,k), dStorage(i,j,k), Net_flux(i,j,k)
            end do
         end do
      end do
      close(unit)
   end subroutine save_fluxes_results


!**********************************************************************
!  SUBROUTINE: gw_cleanup
!  Deallocate arrays
!***********************************************************************
   subroutine GW_cleanup()

      implicit none

      if (allocated(BC_type))then
          deallocate(BC_type) 
      end if
      
      if (allocated(BC_value))then
         deallocate(BC_value)
      end if
      
      if (allocated(source_type))then
         deallocate(source_type)
      end if
      
      if (allocated(source_value))then
         deallocate(source_value)
      end if
      
      if (allocated(KC))then
         deallocate(KC)
      end if
      
      if (allocated(KR))then
        deallocate(KR)
      end if
      
      if(allocated(KV))then
        deallocate(KV)
      end if
      
      if (allocated(CC))then
        deallocate(CC)
      end if
      
      if (allocated(CR))then
        deallocate(CR)
      end if
      
      if (allocated(CV))then
        deallocate(CV)
      end if
      
      if (allocated(SY))then
        deallocate(SY)
      end if
      
      if (allocated(SS))then
        deallocate(SS)
      end if
      
      if (allocated(TERRAIN))then
        deallocate(TERRAIN)
      end if
      
      if (allocated(H_initial))then
        deallocate(H_initial)
      end if
      
      if (allocated(recharge_rate))then
        deallocate(recharge_rate)
      end if
      
      if (allocated(WTD))then
        deallocate(WTD)
      end if
      
      if (allocated(pumping_rate))then
        deallocate(pumping_rate)
      end if
      
      if (allocated(TIME))then
        deallocate(TIME)
      end if
      
      if (allocated(GW_HEAD))then
        deallocate(GW_HEAD)
      end if
      
      if ( allocated(DELTAV))then
        deallocate(DELTAV)
      end if
      
      if (allocated(HNEW))then
        deallocate(HNEW)
      end if
      
      if (allocated(HOLD))then
        deallocate(HOLD)
      end if
      
      if (allocated(GWTOP))then
        deallocate(GWTOP)
      end if
      
      if (allocated(GWBOT))then
        deallocate(GWBOT)
      end if
      
      if (allocated(B_RIVER))then
        deallocate(B_RIVER)
      end if
      
      if (allocated(P_RIVER))then
        deallocate(P_RIVER)
      end if        
      
      if (allocated(H_RIVER))then 
         deallocate(H_RIVER)
      end if
      
      if(allocated(W_RIVER))then
         deallocate(W_RIVER)
      end if

      if(allocated(L_RIVER))then
         deallocate(L_RIVER)
      end if

      if (allocated(RIVBED_BOT))then
        deallocate(RIVBED_BOT)
      end if

      if (allocated(C_RIVBED))then
        deallocate(C_RIVBED)
      end if

      if (allocated(recharge_flux))then
        deallocate(recharge_flux)
      end if

      if (allocated(pumping_flux))then
        deallocate(pumping_flux)
      end if

      if (allocated(river_flux))then
        deallocate(river_flux)
      end if

      if (allocated(lat_xri_flux))then
        deallocate(lat_xri_flux)
      end if

      if (allocated(lat_ycj_flux))then
        deallocate(lat_ycj_flux)
      end if

      if (allocated(vertical_flux))then
        deallocate(vertical_flux)
      end if

      if (allocated(dStorage))then
        deallocate(dStorage)
      end if

      if (allocated(Net_flux))then
        deallocate(Net_flux)
      end if

   end subroutine GW_cleanup

end module SawSim_subroutines
