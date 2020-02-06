      implicit none
      integer MC        ! number of Monte Carlo steps
      integer L         ! layout size
      real TR           ! reduced temperature
      integer i,j,i2,j2,d,K
      real Kb           ! reduced Boltzman constant *10^(23)
      parameter(MC = 20 000, L = 300, Kb = 1.38)
      integer s(L,L)    ! node table, s(i,j) = +-1
      integer next(L),previous(L),meter
      real dU,w,random,m,ms,ms2,ms4,ER,ERS,ERS2,KSI,U,C
      
      open(1,file = 'dat1.dat')
      open(2,file = 'dat2.dat')
      open(3,file = 'dat3.dat')
      open(4,file = 'dat4.dat')
      open(5,file = 'dat5.dat')
      open(6,file = 'dat6.dat')
      
      do i = 1, L    ! tables of neighbourhood nodes
         next(i) = i + 1
         previous(i) = i - 1
      enddo
      next(L) = 1
      previous(1) = L
      
      do i = 1, L    ! a initial random configuration
        do j = 1, L
        call random_number(random)
      if (random < 0.5) then
         s(i,j) = 1
      else
         s(i,j) = -1
      end if
      enddo
      enddo
      
      TR = 2.0
      !do TR = 1.0, 3.0, 0.5     ! loop after temperatures
      ms = 0        ! magnetization < m >, which is the average value of the spin
      ms2 = 0       ! < m^2 >
      ms4 = 0       ! < m^4 >
      ERS = 0       ! energy < E > (for one node)
      ERS2 = 0      ! < E^2 >  (for one node)
      do K = 1, MC  ! loop after Monte Carlo steps
       meter = 0    ! counter
       do i = 1, L  ! loop after every node
        do j = 1, L
        dU = (-2.0)*(-s(i,j))*(s(next(i),j) + s(previous(i),j)
     &       + s(i,next(j))+ s(i,previous(j))) ! energy difference for new configuration
            if (dU < 0) then
               s(i,j) = -s(i,j) ! accept new configuration
               ! we always accept the new configuration if its energy is less or
               ! equal energy from the previous configuration
            else
               call random_number(random)
               w = EXP(-(dU)/TR)
               if (random <= w) then
                  s(i,j) = -s(i,j) ! accept new configuration
                  ! new configuration that increases the energy of the system,
                  ! we only accept with some probability
               end if
            end if
      enddo
      enddo
      if( K >= 30 000 ) then ! thermalization of system
          if( mod(K, 1 00) == 0 ) then ! uncorrelated configurations for calculations
           m = 0
           ER = 0
           meter = meter + 1  ! counter for average values
           do i2 = 1, L      ! magnetization and energy for all system
               do j2 = 1, L
                  m = m + s(i2,j2)  ! magnetization all system
                  ER = ER + (-1.0)*(s(i2,j2))*s(next(i2),j2)  ! energy all system
               enddo
           enddo
           m = m/float(L*L)
           ms = ms+abs(m)
           ms2 = ms2+m*m
           ms4 = ms4+m**4.0
           ERS = ERS+ER
           ERS2 = ERS2+ER*ER
           write(1,*) K, m  ! magnetization in each Monte Carlo step
          endif
      endif
      enddo
      ms = ms/float(meter)
      ms2 = ms2/float(meter)
      ms4 = ms4/float(meter)
      ERS = ERS/float(meter)
      ERS2 = ERS2/float(meter)
      KSI = (L*L/TR) * (ms2 - ms**2.0)
      C = Kb / (L*L*TR**2.0) * (ERS2 - ERS**2.0) ! C *10^(-23)
      U = 1 - (ms4)/(3.0*ms2*ms2)
      write(2,*) TR, m
      write(3,*) TR, KSI
      write(4,*) TR, C
      write(5,*) LOG(abs(1-(TR/2.269))*L),LOG(ms*L**(0.125))
      !enddo
      
1000  format(1x,300I4)   ! write map all of nodes
      do j = 1,L
         write(6,1000)(s(i,j), i = 1, 300)  ! L -> 300
      enddo
      !  GNUPLOT
      !  set palette defined ( 0 0 0 0, 1 1 1 1 )
      !  plot '...' matrix with image
      
      end
