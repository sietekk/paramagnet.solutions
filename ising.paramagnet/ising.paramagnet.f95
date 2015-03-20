!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION:                                                       !
! MCMA simulation of the 1D Ising Model. Computes <E>, <M>, C (heat  !
! capacity), and X (magnetic susceptibility). Note: b=1, and Kb=1    !
!                                                                    !
! AUTHOR INFO:                                                       !
! Author: Michael Conroy                                             !
! Email: sietekk@gmail.com                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program isingparamagnet
implicit none

!DECLARE VARIABLES
integer, dimension(12) :: seed !Seed variable for the random number generator
integer :: s
integer :: N !Number of sites/spins
integer :: Nequil !Number of equlibrium phase steps
integer :: Nmeasure !Number of measurement phase steps
INTEGER :: Nskip !Number of skip sweeps between measurements
integer :: tempstep !Number of iterations to run MCMA
integer :: anglestep !Number of divisions of angle
real :: Tmax, Tmin !Temp values
real :: gamma !Angle changein rads; adjustable to so the Metropolis acceptance rate is about 50%
INTEGER :: sweep, Nskipcount, Nmeasurecount, tempcount, Nequilcount
integer :: t !Counter
real :: dT !Interval change in temperature
real :: beta !Beta value; 1/T
real :: d_theta !Interval change in angle
integer :: thetacount !Counter
real, parameter :: pi = 3.141592653589793239
real :: new_spin_angle, current_spin_angle !Spin values to calculate change in energy
INTEGER :: currentspin !index to identify spin in lattice(:)
real :: energy, magnetization
real :: deltaH
real :: sinF
real :: accept = 0., reject = 0.
integer :: runs = 0 !Acceptance ratio
integer :: ecount, mcount
real, allocatable :: Tarray(:) !Temp values across 'tempstep' intervals
real, allocatable :: lattice(:), coslattice(:) !Spin configuration and cos value array
real, allocatable :: EnergyArray(:,:,:) !For each T(n), holds Energies (n,1,:) and Energies^2 (n,2,:)
real, allocatable :: MagArray(:,:,:) !For each T(n), holds Magnetizations (n,1,:) and Magnetizations^2 (n,2,:)
real, allocatable :: DataArray(:,:) !T (n,1), <E>/N (n,2), <M>/N (n,3), C(n,4), X (n,5)
real, allocatable :: AngleArray(:) !Angle values across 'anglestep' intervals
real, allocatable :: AcceptArray(:) !Holds acceptance ratio of each run

!I/O DATA ENTRY
!Pulls input data from external text file in same directory as executable
open (unit=1, file='input.dat', status='old', action='read')
read(1, *), N, Nequil, Nmeasure, Nskip, Tmin, Tmax, tempstep, anglestep, gamma
close(1)

ALLOCATE(Tarray(tempstep),lattice(N),coslattice(N),EnergyArray(tempstep,2,Nmeasure),&
        MagArray(tempstep,2,Nmeasure),DataArray(tempstep,11), &
        AngleArray(0:anglestep),AcceptArray(Nskip*Nmeasure*N))

!SEED RANDOM NUMBER GENERATOR
!Pull from same pseudo-random number pool each time program is run
seed = (/(12345678, s=1,12)/) !Seed var must be 1x12; given 12345678 for each value
call random_seed(put = seed) !Seed the random number generator

!CONSTRUCT TEMPERATURE ARRAY
dT = (Tmax - Tmin)/(tempstep)
Tarray(:) = (/( Tmax-(t*dT), t=1,tempstep )/) !Populate temperature array

!CONSTRUCT ANGLE ARRAY
d_theta = pi/anglestep
AngleArray(:) = (/( thetacount*d_theta, thetacount=1,anglestep )/)!Populate angle array

!SPIN CONFIGURATION
CALL spin_configuration

!MC ALGORITHM
do tempcount=1,tempstep !Loop through each temperature in Tarray, which also is total # runs for MCMA
    !CALCULATE BETA
    beta = 1. / Tarray(tempcount)

    !EQUILIBRATION PHASE
    DO Nequilcount=1,Nequil
        DO sweep=1,N
            CALL spin_update
            !print *,'Equil phase:',Nequilcount,sweep, runs
        END DO
    END DO

    !MEASUREMENT PHASE
    DO Nmeasurecount=1,Nmeasure
        runs = runs + 1
        DO Nskipcount=1,Nskip
            runs = runs + 1
            DO sweep=1,N
                CALL spin_update
                runs = runs + 1
                !print*,'runs',runs
                !Calculate and Print Acceptance Ratio
                !AcceptArray(runs) = accept/(accept+reject)

                !print*,'Acceptances: ', accept,'Rejections: ',reject
                !print *,'Meas phase:',Nmeasurecount,Nskipcount,sweep, runs
            END DO
        END DO

        !MEASUREMENT
        DO ecount=1,N
             coslattice(ecount) = COS(lattice(ecount))
        END DO
        energy = -(sum(coslattice(:)))
        !print *,'Energy: ',energy
        EnergyArray(tempcount,1,Nmeasurecount) = energy !Stores MCMA energy for a temperature
        EnergyArray(tempcount,2,Nmeasurecount) = energy**2 !Stores the square of that energy (for calculating Cv)
        magnetization = abs(sum(coslattice(:)))
        !print *,'magnetization: ',magnetization
        MagArray(tempcount,1,Nmeasurecount) = magnetization
        MagArray(tempcount,2,Nmeasurecount) = magnetization**2

    END DO

    !print*,'Run: ',tempcount

    !AVERAGES OF OBSERVABLES
    !Temperature
    DataArray(tempcount,1) = Tarray(tempcount)
    !print*,DataArray(tempcount,1),Tarray(tempcount)
    !<E>
    DataArray(tempcount,2) = sum( EnergyArray(tempcount,1,:) ) / Nmeasure
    !print*,DataArray(tempcount,2), sum( EnergyArray(tempcount,1,:) ) / Nmeasure
    !<E^2>
    DataArray(tempcount,3) = sum( EnergyArray(tempcount,2,:) ) / Nmeasure
    !print*,DataArray(tempcount,3), sum( EnergyArray(tempcount,2,:) ) / Nmeasure
    !<E>^2
    DataArray(tempcount,4) = (DataArray(tempcount,2))**2
    !print*, DataArray(tempcount,4), (DataArray(tempcount,1))**2
    !<M>
    DataArray(tempcount,5) = sum( MagArray(tempcount,1,:) ) / Nmeasure
    !print*,DataArray(tempcount,5),sum( MagArray(tempcount,1,:) ) / Nmeasure
    !<M^2>
    DataArray(tempcount,6) = sum( MagArray(tempcount,2,:) ) / Nmeasure
    !print*,DataArray(tempcount,6),sum( MagArray(tempcount,2,:) ) / Nmeasure
    !<M>^2
    DataArray(tempcount,7) = (DataArray(tempcount,5))**2
    !print*,DataArray(tempcount,7),(DataArray(tempcount,5))**2
    !C
    DataArray(tempcount,8) = (DataArray(tempcount,3) - DataArray(tempcount,4)) /(N * Tarray(tempcount)**2)
    !print*,DataArray(tempcount,8), (beta**2 / N) * (DataArray(tempcount,3) - DataArray(tempcount,4))
    !X
    DataArray(tempcount,9) = (beta / N) * (DataArray(tempcount,6) - DataArray(tempcount,7))
    !print*,DataArray(tempcount,9), (beta / N) * (DataArray(tempcount,6) - DataArray(tempcount,7))
    !<E>/N
    DataArray(tempcount,10) = DataArray(tempcount,2) / N
    !print*,DataArray(tempcount,10),DataArray(tempcount,2) / N
    !<M>/N
    DataArray(tempcount,11) = DataArray(tempcount,5) / N
    !print*,DataArray(tempcount,11), DataArray(tempcount,5) / N

!END MC ALGORITHM
END DO

CALL save_data
print*,'runs ', runs
!print*,'accept ratio ', sum(AcceptArray(:))/size(AngleArray)
print*,'Acceptance Ratio: ',accept/(accept+reject)

CONTAINS
!------------------------------------------------------------------------!
    SUBROUTINE spin_configuration
    INTEGER :: count, spincount
    REAL :: randomangle, rand, scaled
    DO spincount=1,N !Loop through lattice
        CALL random_number(rand)
        scaled = rand * anglestep !Scales to number of sites
        DO count=1,(anglestep+1) !Picks spin corresponding to scaled random number
            IF (scaled>=(count-1) .AND. scaled<count) THEN
                randomangle = AngleArray((count-1))
            !ELSE
            !	print *,'scaled, count, randomangle: ',scaled, count, randomangle
            !    STOP 'ERROR: spin_configuration randomangle generation.'
            END IF
        END DO
        lattice(spincount) = randomangle
    END DO
    END SUBROUTINE spin_configuration
!------------------------------------------------------------------------!
    SUBROUTINE spin_update
    REAL :: rand, compare, new_theta
    !Current spin identification and angle
    currentspin = randomspin(N)
    current_spin_angle = lattice(currentspin)

    !Propose change in spin angle
    CALL random_number(rand)
    new_theta = (rand-0.5)*gamma + current_spin_angle
    !print*,'new_theta: ',new_theta

    !Test if theta is in bounds: 0<=theta<=pi
     IF (new_theta <= pi .AND. new_theta >= 0.) THEN !Add temp, if addition is less than pi
        new_spin_angle = new_theta
    	!print*,'< pi: ', new_spin_angle
    ELSE IF (new_theta > pi .AND. new_theta < 2.0*pi) THEN !Subtract temp, if adding is greater than pi
         new_spin_angle = 2.0*pi - new_theta
    	!print*,'> pi: ', new_spin_angle
    ELSE IF (new_theta < 0. .AND. new_theta > -pi) THEN
        new_spin_angle = ABS(new_theta)
    ELSE
        STOP 'Theta bounding test failed!'
    END IF

    !Change in energy
    deltaH =  -(COS(new_spin_angle) - COS(current_spin_angle))

    !Spin angle acceptance/rejections
    CALL random_number(rand)
    sinF = SIN(new_spin_angle)/SIN(current_spin_angle)
    compare = sinF * exp(-(beta * deltaH))
    IF (rand <= compare) THEN !Accept
         lattice(currentspin) = new_spin_angle
         accept = accept + 1.
    ELSE
        reject = reject + 1.
    END IF
    END SUBROUTINE spin_update
!------------------------------------------------------------------------!
    SUBROUTINE save_data
    INTEGER :: i,j,k,l,m,n
    !<E>-------------------------------------------------------------
    OPEN(unit=99,file='e_average.dat',status='replace')
        DO i=1,tempstep
            WRITE(99,*) DataArray(i,1), ' ', DataArray(i,2)
        END DO
    CLOSE(unit=99)
    PRINT *, '<E> data has been saved!'
    !<M>-------------------------------------------------------------
    OPEN(unit=99,file='m_average.dat',status='replace')
        DO j=1,tempstep
            WRITE(99,*) DataArray(j,1), ' ', DataArray(j,5)
        END DO
    CLOSE(unit=99)
    PRINT *, '<M> data has been saved!'
    !Specific Heat-----------------------------------------------------
    OPEN(unit=99,file='specificheat.dat',status='replace')
        DO k=1,tempstep
            WRITE(99,*) DataArray(k,1), ' ', DataArray(k,8)
        END DO
    CLOSE(unit=99)
    PRINT *, 'Specific heat data has been saved!'
    !Magnetization-----------------------------------------------------
    OPEN(unit=99,file='magnetization.dat',status='replace')
        DO l=1,tempstep
            WRITE(99,*) DataArray(l,1), ' ', DataArray(l,9)
        END DO
    CLOSE(unit=99)
    PRINT *, 'Magnetization data has been saved!'
    !<E>/N-------------------------------------------------------------
    OPEN(unit=99,file='e_avg_per_site.dat',status='replace')
        DO m=1,tempstep
            WRITE(99,*) DataArray(m,1), ' ', DataArray(m,10)
        END DO
    CLOSE(unit=99)
    PRINT *, '<E>/N data has been saved!'
    !<M>/N-------------------------------------------------------------
    OPEN(unit=99,file='m_avg_per_site.dat',status='replace')
        DO n=1,tempstep
            WRITE(99,*) DataArray(n,1), ' ', DataArray(n,11)
        END DO
    CLOSE(unit=99)
    PRINT *, '<M>/N data has been saved!'
    END SUBROUTINE save_data
!------------------------------------------------------------------------!
    FUNCTION randomspin(N)
    INTEGER :: N, count, randomspin
    REAL :: rand, scaled
    CALL random_number(rand)
    scaled = rand * N !Scales to number of sites
    DO count=1,(N+1) !Picks spin corresponding to scaled random number
        IF (scaled>=(count-1) .AND. scaled<count) randomspin = count
    END DO
    END FUNCTION randomspin
!------------------------------------------------------------------------!

!DEALLOCATE(Tarray)
!DEALLOCATE(lattice)
!DEALLOCATE(EnergyArray)
!DEALLOCATE(MagArray)
!DEALLOCATE(DataArray)

!END OF PROGRAM
end program isingparamagnet
