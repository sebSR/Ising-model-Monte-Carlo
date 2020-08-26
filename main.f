      implicit none

      real reducedTemp
      integer MonteCarloSteps, latticeSize
      parameter(MonteCarloSteps=500000, latticeSize=14, reducedTemp=2.3)

      integer lattice(latticeSize, latticeSize)
      integer next(latticeSize), previous(latticeSize)
      integer currentStep, i, j, i2, j2
      real energyDifference, random, magnetisation, sum, w

      open(1,file = 'dat1.dat')
      open(2,file = 'dat2.dat')

      !=========================================================================

      ! tables of neighbourhoods
      do i = 1, latticeSize
         next(i) = i + 1
         previous(i) = i - 1
      enddo
      next(latticeSize) = 1
      previous(1) = latticeSize

      !=========================================================================

      ! initial random configuration
      do i = 1, latticeSize
        do j = 1, latticeSize
        call random_number(random)
          if (random < 0.5) then
            lattice(i,j) = 1
          else
            lattice(i,j) = -1
          end if
        enddo
      enddo

      !=========================================================================

      do currentStep = 1, MonteCarloSteps
        do i = 1, latticeSize
          do j = 1, latticeSize
          sum = lattice(next(i),j) + lattice(previous(i),j)
     &          + lattice(i,next(j)) + lattice(i,previous(j))
          energyDifference = 2.0*lattice(i,j)*sum
          if (energyDifference < 0) then
            lattice(i,j) = -lattice(i,j) ! accept new configuration
          else
            call random_number(random)
            w = EXP(-(energyDifference)/reducedTemp)
            if (random <= w) then
              lattice(i,j) = -lattice(i,j) ! accept new configuration
            end if
          end if
          enddo
        enddo

        ! thermalization of system
        if(currentStep >= 30000) then
          ! uncorrelated configurations
          if(mod(currentStep, 100) == 0) then
            magnetisation = 0
            do i2 = 1, latticeSize
              do j2 = 1, latticeSize
                magnetisation = magnetisation + lattice(i2,j2)
              enddo
            enddo
            magnetisation = magnetisation/float(latticeSize*latticeSize)
            write(1,*) currentStep, magnetisation
          endif
        endif
      enddo

      !=========================================================================

1000  format(1x,300I4)
      do j = 1,latticeSize
        write(2, 1000)(lattice(i,j), i = 1, latticeSize)
      enddo

      !=========================================================================

      end
