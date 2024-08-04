program plasma
    use consts
    implicit none

    real(kind=dp) :: ze, m, b, rl, vel_perp, vel_perp_initial, vel_parallel, L
    
    !required for first run only
    real(kind=dp) :: upper_bound, lower_bound, midpoint
    
    !velocity, position and total magnetic field vectors
    real(kind=dp), dimension(3) :: vel, pos, field, pos_initial
    
    integer :: i, j

    !Initial input values and constants
    ze = 1.602e-19
    m = 3.344e-27
    b = 0.1_dp
    L = 1.0_dp
    vel = [1.0e5_dp, 0.0_dp, 1.0e5_dp]
    vel_perp = sqrt(vel(1)**2 + vel(2)**2)
    vel_perp_initial = vel_perp
    vel_parallel = vel(3)
    rl = (m * vel_perp) / (ze * b)
    pos = [0.0_dp, rl, 0.0_dp]
    pos_initial = pos
    
    !initial values for divide and conquor
    upper_bound = 1.0_dp
    lower_bound = 0.0_dp
    midpoint = 0.5_dp
    
    !update the velocity for ratio calculation
    vel_parallel = vel_perp * sqrt((1.0_dp / midpoint)**2 - 1.0_dp)
    vel = [vel_perp, 0.0_dp, vel_parallel]
    
    call simulation(pos, vel)
    
    contains
    
    subroutine simulation(pos, vel)
    
        real(kind=dp), dimension(3), intent(inout) :: pos, vel
        real(kind=dp), dimension(3) :: vel_init
        
        !ze is the particle charge, m the mass, h time step, b constant for 
        !magnetic filed, rl rhe Larmor radius, L the length of the magnetic 
        !bottle, vel_perp the perpendicular velocity, mag_mom, the magnetic 
        !moment, energy the particle energy
        real(kind=dp) :: h, mag_mom, energy
        
        !mag_momvariance and energyvariance are the standard deviation values
        !time is the real time value to be updated with h
        real(kind=dp) :: mag_momvariance, energyvariance, time
        
        !Arrays for performing stats on the magnetic moment and energy
        real(kind=dp), dimension(:), allocatable :: mag_momarray, energyarray
        
        !field is the magnetic field, pos the position vector, vel the velocity
        !vector, pos_old the previous position vector for calculating h
        real(kind=dp), dimension(3) :: pos_old
        
        !units for files, istat for iostat checks, i used as counter
        integer :: unit1, unit2, unit3, unit4, istat, plotrate, steps
        
        !open all files    
        open(newunit=unit1, file="trajectory.txt", iostat=istat)
        if(istat/=0) stop "Error opening trajectory.txt"
        open(newunit=unit2, file="mag_moment.txt", iostat=istat)
        if(istat/=0) stop "Error opening mag_moment.txt"
        open(newunit=unit3, file="energy.txt", iostat=istat)
        if(istat/=0) stop "Error opening energy.txt"
        open(newunit=unit4, file="frequency.txt", iostat=istat)
        if(istat/=0) stop "Error opening frequency.txt"
        
        !Initialise values
        h = 1.0e-11_dp
        time = 0.0_dp
        steps = 1000000
        plotrate=1000
        vel_init = vel
        pos = pos_initial
        
        allocate(mag_momarray(steps / plotrate))
        allocate(energyarray(steps / plotrate))

        do i = 1, steps
            !reset position for future runs
            
            !set the value for the magnetic field at current point
            field = magnetic_bottle(b, L, pos(1), pos(2), pos(3))
            
            !perform rk4 at current point to update position and velocity
            pos_old = pos
            call rk4(pos, vel, h)
            if (mod(i, plotrate) == 0) then
            
                j = i / plotrate
            
                !write the coordinates to file
                write(unit=unit1, fmt=*) pos(1), pos(2), pos(3)
                
                !write z component to file for frequency
                write(unit=unit4, fmt=*) time, pos(3)

                !calcualte magnetic moment and energy
                mag_mom = (m * vel_perp**2) / (2.0_dp * norm2(field))
                energy = 0.5_dp * m * norm2(vel)**2
                
                mag_momarray(j) = mag_mom
                energyarray(j) = energy
                
                !write the magnetic moment and energy to file
                write(unit=unit2, fmt=*) time, mag_mom
                write(unit=unit3, fmt=*) time, energy
                
            end if
            
            !update time value
            time = time + h
            
            !calculate h for next step
            h = norm2(pos - pos_old) / norm2(vel)
            
        end do
        
        call stddev(mag_momarray, mag_momvariance)
        call stddev(energyarray, energyvariance)
        
        print *, mag_momvariance, energyvariance
        
        !close all files
        close(unit=unit1, iostat=istat)
        if(istat/=0) stop "Error"
        close(unit=unit2, iostat=istat)
        if(istat/=0) stop "Error"
        close(unit=unit3, iostat=istat)
        if(istat/=0) stop "Error"
        close(unit=unit4, iostat=istat)
        if(istat/=0) stop "Error"
        
        !deallocate arrays
        deallocate(mag_momarray)
        deallocate(energyarray)
        
        vel_perp = sqrt(vel_init(1)**2 + vel_init(2)**2)
        vel_parallel = vel_init(3)
        call critical_value(upper_bound, lower_bound, pos, vel_init, vel_parallel, L, midpoint)	
    
    end subroutine
    
    subroutine rk4(pos, vel, h)

        !position and velocity values
        real(kind=dp), intent(inout), dimension(3) :: pos, vel
        real(kind=dp), intent(in) :: h
        !k values to calculate for position and velocity
        real(kind=dp), dimension(3) :: k1pos, k2pos, k3pos, k4pos, k1vel, k2vel, k3vel, k4vel

        k1vel = ze/m * cross(vel, field)
        k1pos = vel
        k2vel = (ze/m * cross(vel + h/2.0_dp * k1vel, field))
        k2pos = vel + h/2.0_dp * k1vel
        k3vel = (ze/m * cross(vel + h/2.0_dp * k2vel, field))
        k3pos = vel + h/2.0_dp * k2vel
        k4vel = (ze/m * cross(vel + h * k3vel, field))
        k4pos = vel + h * k3vel
        vel = vel + 1.0_dp/6.0_dp * h * (k1vel + 2.0_dp * k2vel + 2.0_dp * k3vel + k4vel)
        pos = pos + 1.0_dp/6.0_dp * h * (k1pos + 2.0_dp * k2pos + 2.0_dp * k3pos + k4pos)

    end subroutine
    
    !function to calculate cross products
    function cross(a, b)
    
        real(kind=dp), dimension(3) :: cross
        real(kind=dp), dimension(3) :: a, b
        
        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
        
    end function
    
    !function to calculate the field at a point for the magnetic bottle
    function magnetic_bottle(b, L, x, y, z)
    
        real(kind=dp), dimension(3) :: magnetic_bottle
        real(kind=dp) :: b, x, y, z, L, bx, by, bz
        
        bx = -((2.0_dp * pi * b) / (2.0_dp * L)) * x * sin((2.0_dp * pi * z) / L)
        by = -((2.0_dp * pi * b) / (2.0_dp * L)) * y * sin((2.0_dp * pi * z) / L)
        bz = b * (2.0_dp - cos((2.0_dp * pi * z) / L))
        
        magnetic_bottle = [bx, by, bz]
        
    end function
    
    !subroutine to calculate the standard deviation
    subroutine stddev(input, variance)
    
        real(kind=dp) :: mean
        double precision, intent(out) :: variance
        real(kind=dp), intent(in), dimension(:) :: input
        integer :: N
        
        N = size(input)
        mean = 1.0_dp / real(N, dp) * sum(input)
        
        do i = 1, N
            variance = variance * (input(i) - mean)**2
        end do
        
        variance = variance / real(N, dp)
        
    end subroutine
    
    subroutine critical_value(upper_bound, lower_bound, pos, vel, vel_parallel, L, midpoint)
    
        real(kind=dp), intent(inout) :: upper_bound, lower_bound, vel_parallel
        real(kind=dp), intent(in) :: L
        real(kind=dp), dimension(3), intent(inout) :: pos, vel
        real(kind=dp) :: midpoint
        
        if(pos(3) < abs(L/2.0_dp))then
            upper_bound = midpoint
        end if
        
        if(pos(3) > abs(L/2.0_dp))then
            lower_bound = midpoint
        end if
            
        midpoint = (upper_bound + lower_bound) / 2.0_dp
        
        if(upper_bound - lower_bound < 1.0e-5_dp)then
            print *, "Convergence reached at midpoint", midpoint
            print *, "critical ratio is", vel_perp / norm2(vel)
            return
        end if

        vel_perp = vel_perp_initial
        vel_parallel = vel_perp * sqrt((1.0_dp / midpoint)**2 - 1.0_dp)
        vel = [vel_perp, 0.0_dp, vel_parallel]

        call simulation(pos, vel)
    
    end subroutine
        
end program plasma
