module m_thermalplasma

 !--- DEPENDENCIES  -----------------------------------------!
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_particlelogistics
  implicit none

  
  !--- PRIVATE variables -----------------------------------------!
  real, private :: t_crit = 1, t_nonrel = 0.01      ! Should be (t_crit > t_nonrel)  --> chooses method of particle injection
  !...............................................................!


  !--- PRIVATE functions -----------------------------------------!
  private :: tabulateMaxwellian
  private :: deallocateMaxwellian
  private :: sample_u_from_table    ! Added -- not in Vanilla version.
  !...............................................................!

  
  
  
  
  ! ----------------------------------------------------------------------------
  !         type      maxwellian
  ! 
  ! Defines a Maxwellian PDF that can then be filled later. The filling method 
  ! depends on the temperature kT/mc^2. 
  ! ----------------------------------------------------------------------------
  type :: maxwellian
    ! tabulated maxwellian is used 1d/2d for any T, and 3d for T < t_crit
    ! ... in all cases if T < t_nonrel -- we use non-relativistic maxwellian
    ! ... and `u_table` contains `beta` instead of `4-velocity`
    real :: temperature, shift_gamma
    real :: n1 = 0.0, n2 = 0.0, n3 = 1.0            ! unit vector for arbitrary drift (used if abs(shift_dir)==9) 
    real, allocatable, dimension(:) :: DF_table     ! Table of probabilites
    real, allocatable, dimension(:) :: u_table      ! Velocities corresponding to the probs
    logical :: generated = .false., shift_flag = .false.
    integer :: npoints, shift_dir
    integer :: dimension = 3
  end type maxwellian
  !.............................................................................
  
  
  
  
  
  
  
  
!-------------------------------------------------------------------------------
!                       DEFINE FUNCTIONS BELOW
!-------------------------------------------------------------------------------
  
  

contains
  ! See more details in Zenitani 2015
  !   arXiv:1504.03910v1
  
  
  
  !-------------------------------------------------------------------------------
  !           subroutine   deallocateMaxwellian() 
  !
  ! Deallocate a Maxwellian if need be to save space
  !-------------------------------------------------------------------------------
  subroutine deallocateMaxwellian(maxw)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    if (allocated(maxw % DF_table)) deallocate (maxw % DF_table)
    if (allocated(maxw % u_table)) deallocate (maxw % u_table)
  end subroutine deallocateMaxwellian
  
  
  
  
  !-------------------------------------------------------------------------------
  !           subroutine   tabulateMaxwellian() 
  !
  ! Generate a Maxwellian CDF for a given temperature that can later be interpolated.
  ! Currently, only used if t_nonrel < t < t_crit.
  ! More details can be found in Sect. 2.A. of Zenitani (2015)
  !
  !   NOTE: Since this is for one temperature, the CDF is normalized s.t. there is
  !         no need for the 2nd-order Bessel function.
  !   NOTE: By default, this is tabulated with 1,000,000 points and a maximum velocity
  !         of MAX( 20*T , 1)  (T=temperature).
  !-------------------------------------------------------------------------------
  subroutine tabulateMaxwellian(maxw)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    integer :: iter
    real :: u_1, u_2, u_max, df = 0.0, temp, norm, mn, mx
    real :: g1, g2, e1, e2
    !.........  ..............  .............


    ! Check table existance.
    if (maxw % generated) then
      call throwError('ERROR: maxwell table already generated.')
    end if
    
    
    ! Create velocity grid based on temperature and allocate tables 
    temp = maxw % temperature
    maxw % npoints = 1000000
    u_max = MAX( 20.*temp , 1.)       ! Max velocity for velocity grid

    allocate (maxw % DF_table(maxw % npoints))
    allocate (maxw % u_table(maxw % npoints))


    ! Loop over vel. grid and fill in relative probs based on dim of CDF
    do iter = 1, maxw % npoints
    
      ! Update prev vel (u_1) and vel to find CDF value (u_2) 
      u_1 = u_max * REAL(iter - 1) / REAL(maxw % npoints)   
      u_2 = u_max * REAL(iter) / REAL(maxw % npoints)
      
      ! Re-normalize by factoring out the exponent ---> **NEEDED TO AVOID UNDERFLOW IN NORM**
      !     exp(-γ/T) = exp(-1/T)exp(-(γ-1)/T) --> would underflow with (γ,T)=(1.0,0.1) --> exp(-1.0/0.01)~3e-44
      !                                             Due to norm later, factor of exp(-1/T) does not matter
      g1 = sqrt(1.0 + u_1*u_1)
      g2 = sqrt(1.0 + u_2*u_2)
    
      e1 = exp( -(g1 - 1.)/temp )
      e2 = exp( -(g2 - 1.)/temp )
    
      if (maxw % dimension .eq. 1) then
        df = (u_2 - u_1) * 0.5 * (e1 + e2)
      else if (maxw % dimension .eq. 2) then
        df =  (u_2 - u_1) * 0.5 * (e1*u_1 + e2*u_2)
      else if (maxw % dimension .eq. 3) then
        df =  (u_2 - u_1) * 0.5 * (e1*u_1*u_1 + e2*u_2*u_2)
      end if
      
      ! Allocate to table
      maxw % u_table(iter) = u_2
      
      if (iter .eq. 1) then
        maxw % DF_table(iter) = df
      else
        maxw % DF_table(iter) = maxw%DF_table(iter-1) + df
      end if
      
    end do
    
    
    
    ! Normalize the table to the last probability (since CDF)
    norm = maxw % DF_table(maxw % npoints)
    
    if (norm /= norm) then
      call throwError('ERROR: norm is NaN.')
    end if
    
    if (abs(norm) > huge(norm)) then
      call throwError('ERROR: norm is Inf.')
    end if


    
    if (norm <= tiny(norm)) then
      call throwError('ERROR: Maxwellian CDF underflow: DF_table(end) ~ 0. Reduce u_max and/or use higher precision.')
    end if
    
    do iter = 1, maxw % npoints
      maxw % DF_table(iter) = maxw % DF_table(iter) / norm
    end do
    
    
    
    ! CDF should start near 0 and end at 1
    mn = minval(maxw % DF_table)
    mx = maxval(maxw % DF_table)
    
    if (mn /= mn) then
      call throwError('ERROR: mn is NaN.')
    end if
    
    if (mx /= mx) then
      call throwError('ERROR: mn is NaN.')
    end if
    
    if (abs(maxw % DF_table(maxw % npoints) - 1.0) > 1e-6) then
      call throwError('ERROR: Maxwellian CDF does not end at 1 after normalization.')
    end if
    
    ! if (mn < -1e-6 .or. mx > 1.0 + 1e-6) then
    !   call throwError('ERROR: Maxwellian CDF out of [0,1] bounds after normalization.')
    ! end if
    
    ! Monotonic check (binary search requires this)
    do iter = 2, maxw % npoints
      if (maxw % DF_table(iter) < maxw % DF_table(iter-1)) then
        call throwError('ERROR: Maxwellian CDF is not monotonic increasing.')
      end if
    end do
    
    ! Check allocation final time and specify generated
    if (.not. allocated(maxw % DF_table) .or. .not. allocated(maxw % u_table)) then
      call throwError('ERROR: Allocation failed within tabulateMaxwellian().')
    else
      maxw % generated = .true.
    end if
        
  end subroutine tabulateMaxwellian
  !.....................................................................................
  
  
  
  
  


  !-------------------------------------------------------------------------------
  !           function    sample_u_from_table()
  !
  ! Draw a random velocity from a Mawellian CDF. Uses a CDF generated from 
  ! tabulateMaxwellian(). 
  !
  ! Does so by drawing a random number from a uniform distribution, then using a binary
  ! search to place on CDF and invert to a velocity. The velocity is found by interpolating
  ! the CDF.
  !-------------------------------------------------------------------------------
  function sample_u_from_table(maxw)
    implicit none
    real :: sample_u_from_table
    type(maxwellian), intent(in) :: maxw
    integer :: lo, hi, mid
    real :: x, f1, f2, u1, u2, t
    !......   .......  ............
    
    
    if (.not. allocated(maxw % DF_table) .or. .not. allocated(maxw % u_table)) then
      call throwError('ERROR: sample_u_from_table called before CDF table allocated on this rank.')
    end if


    ! Draw from Uniform Dist
    x = random(dseed)


    ! If outside the probability range of the CDF (shouldn't happen but just to be safe)
    if (x <= maxw%DF_table(1)) then
      sample_u_from_table = maxw%u_table(1)    ! Below min prob of CDF
      return
    end if
    if (x >= maxw%DF_table(maxw%npoints)) then
      sample_u_from_table = maxw%u_table(maxw%npoints)  ! Above max prob of CDF
      return
    end if


    ! binary search: find lo,hi s.t. DF(lo) <= x < DF(hi)
    lo = 1
    hi = maxw%npoints
    do while (hi - lo > 1)
      mid = (lo + hi) / 2       ! Integer division
      if (maxw%DF_table(mid) < x) then
        lo = mid
      else
        hi = mid
      end if
    end do
    
    ! Store for interpolation
    f1 = maxw%DF_table(lo);  f2 = maxw%DF_table(hi)
    u1 = maxw%u_table(lo);   u2 = maxw%u_table(hi)

    if (f2 <= f1) then
      ! should not happen if DF_table is monotonic, but guard anyway
      sample_u_from_table = u1
      return
    end if
    
    ! Interpolate the velocity based on the CDF
    t = (x - f1) / (f2 - f1)
    sample_u_from_table = u1 + t*(u2 - u1)
    return
  end function sample_u_from_table
  !...............................................................................
  
  
  
  
  
  
  
  
  
  

  !-------------------------------------------------------------------------------
  !           subroutine    generateFromMaxwellian()    
  !
  ! Draw a velocity from a Maxwellian distribution with a method based on the temp.
  !     
  !                t < t_nonrel     --> Draw from normal distribution via Leva (1992)
  !         t_nonrel <= t < t_crit   --> Draw from pre-tabulated MJ CDF
  !           t_crit <= t            --> Draw from MJ via the Sobol method
  !     
  !     NOTE: For middle case, if table does not exist for the inputted dist, then
  !           tabulateMaxwellian() is called to generate the CDF.
  !-------------------------------------------------------------------------------
  subroutine generateFromMaxwellian(maxw, u_, v_, w_)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    real, intent(out) :: u_, v_, w_
    real :: U = 0.0, ETA = 0.0, X1, X2, X3, X4, X5, X6, X7, X8, dx1, dx2, BETA, gamma
    real :: n1, n2, n3, beta0, gamma0, up, uperp1, uperp2, uperp3, b_shift, g_shift
    logical :: flag
    integer :: iter, dim_
    !.......   .........   ............    ..........   .........   ...........
    
    
    ! Velocity components to return
    u_ = 0.0; v_ = 0.0; w_ = 0.0
    
    
    ! --------------------------
    ! Generate vel. w/ method based on temp.
    ! --------------------------

    ! NR --> NORMAL DISTRIBUTION VIA REJECTION
    if (maxw % temperature .lt. t_nonrel) then
      
      ! Begin Leva Rejection for each dimension
      do dim_ = 1, maxw % dimension
        X8 = 1.0; X1 = 1.0; X2 = 1.0
        do while ((X8 .gt. 0.27597) .and. &
                  ((X8 .gt. 0.27846) .or. (X1 .lt. 1e-16) .or. (X2**2 .gt. -4 * log(X1) * X1**2)))
          X1 = random(dseed)
          X2 = 1.7156 * (random(dseed) - 0.5)
          X3 = X1 - 0.449871
          X4 = abs(X2) + 0.386595
          X8 = X3**2 + X4 * (0.196 * X4 - 0.25472 * X3)
        end do
      
        ! Assign velocity component based on dimension being considered
            ! vel. is rescaled based on sqrt(temperature)
        if (dim_ .eq. 1) then
          u_ = X2 / X1 * sqrt(maxw % temperature)
        else if (dim_ .eq. 2) then
          v_ = X2 / X1 * sqrt(maxw % temperature)
        else if (dim_ .eq. 3) then
          w_ = X2 / X1 * sqrt(maxw % temperature)
        end if
      end do
      
      
    ! MILDLY REL --> DRAW FROM TABULATED MJ CDF  
    else if ( (maxw % temperature .ge. t_nonrel) .and. (maxw % temperature .lt. t_crit) ) then
      
      ! Ensure table exists on this rank
      if (.not. allocated(maxw%DF_table) .or. .not. allocated(maxw%u_table)) then
        write(*,'(A,1X,I0,1X,A,1X,L1,1X,L1,1X,A,1X,L1)') &
          "RANK", mpi_rank, "alloc(DF,u)=", allocated(maxw%DF_table), allocated(maxw%u_table), &
          " generated=", maxw%generated
        call tabulateMaxwellian(maxw)
      else
        ! If someone set generated=.true. but didn’t allocate (bug), catch it early
        if (.not. allocated(maxw%DF_table) .or. .not. allocated(maxw%u_table)) then
          call throwError("ERROR: maxw%generated true but CDF table not allocated on this rank.")
        end if
      end if
      write(*,'(A,1X,I0,1X,A,1X,L1,1X,L1,1X,A,1X,L1)') &
          "RANK", mpi_rank, "alloc(DF,u)=", allocated(maxw%DF_table), allocated(maxw%u_table), &
          " generated=", maxw%generated

      U = sample_u_from_table(maxw)     ! Draw Vel. magnitude
      
    
    ! HIGHLY REL --> MJ VIA REJECTION
    else if (maxw % temperature .ge. t_crit) then
    
      ! Begin Sobol Rejection
      ETA = 0.0; U = 0.0
      do while (ETA**2 - U**2 .le. 1)
        X4 = random(dseed); X5 = random(dseed)
        X6 = random(dseed); X7 = random(dseed)
        X8 = X4 * X5 * X6 * X7
        if (X8 .lt. 1e-16) cycle
        U = -maxw % temperature * log(X8 / X7)
        ETA = -maxw % temperature * log(X8)
      end do
    
    end if
    
    
    ! -------------------
    ! Proj. vel. mag. (already done in NR case)
    ! -------------------
    
    if ( maxw % temperature .ge. t_nonrel ) then
      
      ! Project U mag from spherical coords (see Zenitani 2015)
      
      if (maxw % dimension .eq. 1) then
        u_ = U * SIGN(1.0, 0.5 - random(dseed))
        
      else if (maxw % dimension .eq. 2) then
        X1 = 2.0 * M_PI * random(dseed)
        u_ = U * cos(X1)
        v_ = U * sin(X1)
        
      else if (maxw % dimension .eq. 3) then
        X1 = 1.0 - 2.0 * random(dseed)
        X2 = 2.0 * M_PI * random(dseed)
        w_ = U * X1
        X1 = sqrt(1.0 - X1**2)
        u_ = U * X1 * cos(X2)
        v_ = U * X1 * sin(X2)
      end if
    end if
      
      
    ! ----------------------
    ! Shift dist. if there's a bulk flow
    !       --> also corrects volume (see Zenitani Sect 3B)
    ! ----------------------
    
    if (maxw % shift_flag) then
        
      gamma = sqrt(1.0 + u_**2 + v_**2 + w_**2)         ! Lorentz-factor of particle
      g_shift = maxw % shift_gamma                      ! Lorentz factor of bulk flow
      b_shift = sqrt( 1.0 - g_shift**(-2) ) ! β = v/c = u/γ of bulk flow
      
      X8 = random(dseed)    ! Draw from uniform dist to help with volume correction a la Zenitani (2015)
      
      
      ! SHIFT BASED ON SPECIFIED DIRECTION
      !     --> In each case, vol. correction is done using rejection algorithm with random variable
      select case (maxw % shift_dir)
      
      case (+1) ! +x
        if (-b_shift * u_ / gamma .gt. X8) u_ = -u_
        u_ = g_shift * (u_ + b_shift * gamma)
      case (-1) ! -x
        b_shift = -b_shift
        if (-b_shift * u_ / gamma .gt. X8) u_ = -u_
        u_ = g_shift * (u_ + b_shift * gamma)
      case (+2) ! +y
        if (-b_shift * v_ / gamma .gt. X8) v_ = -v_
        v_ = g_shift * (v_ + b_shift * gamma)
      case (-2) ! -y
        b_shift = -b_shift
        if (-b_shift * v_ / gamma .gt. X8) v_ = -v_
        v_ = g_shift * (v_ + b_shift * gamma)
      case (+3) ! +z
        if (-b_shift * w_ / gamma .gt. X8) w_ = -w_
        w_ = g_shift * (w_ + b_shift * gamma)
      case (-3) ! -z
        b_shift = -b_shift
        if (-b_shift * w_ / gamma .gt. X8) w_ = -w_
        w_ = g_shift * (w_ + b_shift * gamma)
        
        
      case default
        if ( abs(maxw % shift_dir) == 9 ) then
        
          n1=maxw%n1;   n2=maxw%n2;    n3=maxw%n3     ! Components of direction unit vector of bulk flow
          up = u_*n1 + v_*n2 + w_*n3    ! Vel. component along shift direction
          
          ! Flip sign of shift if shift_dir == -9
          if ( maxw % shift_dir < 0 ) b_shift  = -b_shift   
          
          ! Reflect across plane perpendicular to n (if needed)
          !     i.e. make up --> -up       (equiv. to u_ --> -u_ for cases +/- 1)
          if ( -b_shift * up / gamma  .gt.  X8) then
            u_ = u_ - 2*up*n1
            v_ = v_ - 2*up*n2
            w_ = w_ - 2*up*n3
            
            up = -up
          end if
          
          ! Decompose into perpendicular components
          uperp1 = u_ - up*n1
          uperp2 = v_ - up*n2
          uperp3 = w_ - up*n3

          ! Lorentz boost of four-velocity along n:
          ! u_parallel' = γ_s * (u_parallel + β_s * γ)
          up = g_shift * (up + b_shift * gamma)

          ! perpendicular component unchanged; recompose
          u_ = uperp1 + up * n1
          v_ = uperp2 + up * n2
          w_ = uperp3 + up * n3
          
        end if
        
      end select
    end if
    
  end subroutine generateFromMaxwellian
  !.................................................................................
  
  
  
  
  
  
  
 
  !-------------------------------------------------------------------------------
  !           subroutine    fillingRegionWithThermalPlasma()    
  !
  ! Given a specified region of the simulation, fills the region with a plasma 
  ! 
  !     NOTE: The velocity distribution used for the plasma is based on the temperature of the plasma
  !           and is decided within generateFromMaxwellian()
  !     NOTE: This assumes that the total charge density of the injected species is 0 to maintain
  !           charge neutrality since the particles of each species are injected at the same location
  !-------------------------------------------------------------------------------
  subroutine fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp, &
                                         temperature, shift_gamma, shift_dir, zero_current, &
                                         dimension, weights, spat_distr_ptr, &
                                         dummy1, dummy2, dummy3)
    implicit none
    ! assuming that the charges of all species given in `fill_species` add up to `0`
    type(region), intent(in) :: fill_region
    integer, intent(in) :: num_species
    integer, intent(in) :: fill_species(num_species)
    real, intent(in) :: ndens_sp, temperature
    real, optional, intent(in) :: shift_gamma
    integer, optional, intent(in) :: shift_dir, dimension
    type(maxwellian) :: fill_maxwellian
    integer :: num_part, n, s, spec_, dimension_
    integer(kind=2) :: xi_, yi_, zi_
    real :: fill_xmin, fill_xmax, &
            fill_ymin, fill_ymax, &
            fill_zmin, fill_zmax
    real :: u_, v_, w_, dx_, dy_, dz_
    real :: x_, y_, z_, rnd, num_part_r
    real :: x_glob, y_glob, z_glob

    real, intent(in), optional :: weights
    real :: weights_, rnd_num
    logical, intent(in), optional :: zero_current
    logical :: zero_current_

    procedure(spatialDistribution), pointer, intent(in), optional :: spat_distr_ptr
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real :: dummy1_, dummy2_, dummy3_
    real :: nrm, beta0
    
    
    ! ------------ CHECK EXISTANCE OF VARIABLES ---------------------

    if (.not. present(dimension)) then  ! Check Dimension
      dimension_ = 3
    else
      dimension_ = dimension
    end if

    if (.not. present(zero_current)) then
      zero_current_ = .false.
    else
      zero_current_ = zero_current
    end if

    if (.not. present(weights)) then
      weights_ = 1.0
    else
      weights_ = weights
    end if

    if (present(dummy1)) then
      dummy1_ = dummy1
    else
      dummy1_ = 0.0
    end if
    if (present(dummy2)) then
      dummy2_ = dummy2
    else
      dummy2_ = 0.0
    end if
    if (present(dummy3)) then
      dummy3_ = dummy3
    else
      dummy3_ = 0.0
    end if
    ! ===================================================================
    
    
    ! ----------------- BEGIN FILLING MAXWELLIAN -------------------------------

    fill_maxwellian % dimension = dimension_    ! Set dimension of Maxwellian
    fill_maxwellian % generated = .false.       ! Determine that Maxwellian has not yet be generated (will set when tabulated)
    if (present(shift_gamma)) then
      fill_maxwellian % shift_gamma = abs(shift_gamma)
      fill_maxwellian % shift_flag = .true.
    else
      fill_maxwellian % shift_flag = .false.
    end if
    
    if (present(shift_dir)) then
      if (abs(shift_dir) == 9) then
        nrm = sqrt(dummy1_**2 + dummy2_**2 + dummy3_**2)
        if (nrm > 0.0) then
          fill_maxwellian % n1 = dummy1_ / nrm
          fill_maxwellian % n2 = dummy2_ / nrm
          fill_maxwellian % n3 = dummy3_ / nrm
        else
          fill_maxwellian % n1 = 0.0; fill_maxwellian % n2 = 0.0; fill_maxwellian % n3 = 1.0
        end if
        
        if (zero_current_) then
          fill_maxwellian % shift_dir = 9
        else
          fill_maxwellian % shift_dir = 9 * INT(SIGN(1.0, species(spec_) % ch_sp))
        end if
        
      else
        if (zero_current_) then
          fill_maxwellian % shift_dir = shift_dir
        else
          fill_maxwellian % shift_dir = INT(SIGN(1.0, species(spec_) % ch_sp)) * shift_dir
        end if
      end if
    end if
    
    

    ! global to local coordinates
#ifdef oneD
    call globalToLocalCoords(fill_region % x_min, 0.0, 0.0, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, 0.0, 0.0, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin)
#elif defined(twoD)
    call globalToLocalCoords(fill_region % x_min, fill_region % y_min, 0.0, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, fill_region % y_max, 0.0, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin) &
                 * (fill_ymax - fill_ymin)
#elif defined(threeD)
    call globalToLocalCoords(fill_region % x_min, fill_region % y_min, fill_region % z_min, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, fill_region % y_max, fill_region % z_max, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin) &
                 * (fill_ymax - fill_ymin) &
                 * (fill_zmax - fill_zmin)
#endif



    if (num_part_r .lt. 10.0) then
      if (num_part_r .ne. 0.0) then
        num_part_r = poisson(num_part_r)
      else
        num_part_r = 0.0
      end if
    else
      num_part_r = CEILING(num_part_r)
    end if
    num_part = INT(num_part_r)

    n = 0
    do while (n .lt. num_part)
      ! generate coords for all species
      call generateCoordInRegion(fill_xmin, fill_xmax, fill_ymin, fill_ymax, fill_zmin, fill_zmax, &
                                 x_, y_, z_, xi_, yi_, zi_, dx_, dy_, dz_)

      ! if spatial distribution function is present, compute it
      !   otherwise use uniform distribution
      if (present(spat_distr_ptr)) then
        x_glob = REAL(this_meshblock % ptr % x0) + x_
        y_glob = REAL(this_meshblock % ptr % y0) + y_
        z_glob = REAL(this_meshblock % ptr % z0) + z_
        rnd = spat_distr_ptr(x_glob=x_glob, y_glob=y_glob, z_glob=z_glob, &
                             dummy1=dummy1_, dummy2=dummy2_, dummy3=dummy3_)
      else
        rnd = 1.0
      end if
      rnd_num = random(dseed)
      if ((.not. present(spat_distr_ptr)) .or. (rnd_num .lt. rnd)) then
        do s = 1, num_species
          ! generate momenta for every species individually
          spec_ = fill_species(s)
          if ((spec_ .le. 0) .or. (spec_ .gt. nspec)) then
            call throwError('Wrong species specified in fillRegionWithThermalPlasma.')
          end if
          
                

          if (temperature .gt. 0) then
            fill_maxwellian % temperature = temperature / species(spec_) % m_sp             ! **********************************************************
            call generateFromMaxwellian(fill_maxwellian, u_, v_, w_)
            
          else
            u_ = 0.0; v_ = 0.0; w_ = 0.0
            if (fill_maxwellian % shift_flag) then
              if (abs(fill_maxwellian % shift_dir) .eq. 1) then
                u_ = SIGN(1, fill_maxwellian % shift_dir) * fill_maxwellian % shift_gamma * &
                     sqrt(1.0 - fill_maxwellian % shift_gamma**(-2))
              else if (abs(fill_maxwellian % shift_dir) .eq. 2) then
                v_ = SIGN(1, fill_maxwellian % shift_dir) * fill_maxwellian % shift_gamma * &
                     sqrt(1.0 - fill_maxwellian % shift_gamma**(-2))
              else if (abs(fill_maxwellian % shift_dir) .eq. 3) then
                w_ = SIGN(1, fill_maxwellian % shift_dir) * fill_maxwellian % shift_gamma * &
                     sqrt(1.0 - fill_maxwellian % shift_gamma**(-2))
              else if (abs(fill_maxwellian % shift_dir) .eq. 9) then
                beta0 = sqrt(max(0.0, 1.0 - fill_maxwellian % shift_gamma**(-2)))
                if (fill_maxwellian % shift_dir < 0) beta0 = -beta0
                u_ = fill_maxwellian % shift_gamma * beta0 * fill_maxwellian % n1
                v_ = fill_maxwellian % shift_gamma * beta0 * fill_maxwellian % n2
                w_ = fill_maxwellian % shift_gamma * beta0 * fill_maxwellian % n3

              end if
            end if
          end if
          call createParticle(spec_, xi_, yi_, zi_, dx_, dy_, dz_, u_, v_, w_, &
                              weight=weights_)
        end do
      end if
      n = n + 1
    end do
    call deallocateMaxwellian(fill_maxwellian)
  end subroutine fillRegionWithThermalPlasma
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
end module m_thermalplasma
