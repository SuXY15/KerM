program mixing_test
    use emst_subs, only: emst_rnu

    implicit none
	character(len=64)::filename
    character(len=64)::casename = "PoF_1996_Fig9b"
    character(len=64)::tmp_str
    character(len=64)::model_name
	integer:: npd                             ! dimension for (maximum number of) particles 
	integer:: ncomp                           ! dimension for compositions
    integer:: np, info, i, j, lu=7, nvar, ivar
    real(kind(1.e0)):: omega, dt, omdt, omt, omt_end
    real(kind(1.e0)):: vf0, vf, u(1), cvars(6), tmp_input
    real(kind(1.e0)):: t0, t1, timecost
    real(kind(1.e0)), allocatable:: vars(:)
    real(kind(1.e0)), allocatable:: wt(:)     ! (npd), weight of each particle
    real(kind(1.e0)), allocatable:: state(:)  ! (npd), state (age) of each particle
	real(kind(1.e0)), allocatable:: f(:,:)    ! (npd,ncomp) scalar values of particle compositions
    real(kind(1.e0)), allocatable:: fscale(:) ! (ncomp) scale in each composition dimension
                                              ! in this test, scales are set to be uniform (1)
                                              ! in typical combustion cases:
                                              !     scale for species = 0.01
                                              !     scale for enthaply = 1e16
    
    model_name = "KerM"

    ! ==============================
    ! settings for number of particles and compositions
    npd = 50000
    ncomp = 1
    np = npd
    
    ! settings for mixing times
    omega = 2.0
    dt = 4e-3
    omt_end = 10.0 	    ! duration of whole mixing simulation
	omdt = omega*dt     ! time step in mixing process

    ! ==============================
    ! preparation
    allocate( f(npd,ncomp), state(npd), wt(npd), fscale(ncomp))
    fscale(1:ncomp) = 1.0  ! set scales to be 1.0
    cvars(1:6) = 0.0       ! set contol variables to default settings

 	! load initial samples from data file
 	open(11, file='data/'//trim(adjustl(casename))//'_samples.txt')
        do i=1,np
            read(11,*) f(i,1:ncomp), wt(i)
        enddo
    close(11)
     
    ! load variance states from data file
    open(31, file='data/'//trim(adjustl(casename))//'_variances.txt')
        read(31,*) tmp_input
        nvar = floor(tmp_input)
        allocate( vars(nvar) )
        do i=1,nvar
            read(31,*) vars(i)
        enddo
    close(31)

    ! compute the initial variances
	call compute_var(vf0)	
    
    ! ==============================
    ! start mixing
    omt = 0.0
    ivar = 1
    timecost = 0.0
	do while( omt < omt_end )
        call compute_var(vf)	! calculate and output variance function
        
        if (sqrt(vf/vf0) < vars(ivar)/vars(1)) then
            write(*,'(3e15.7)') omt, vf, vf/vf0
            
            ! save data
            write(tmp_str,'(f12.6)') vars(ivar)
            filename = 'data/'//trim(adjustl(casename))//'_fortran_'// &
                        trim(adjustl(model_name))//'_'//trim(adjustl(tmp_str))//'.txt'
            open(201, file=filename)
                do i = 1, np
                    write(201,'(5e15.7)')  f(i,1:ncomp), wt(i)
                enddo
            close(201)

            ivar = ivar + 1
            if (ivar > nvar) then
                goto 888
            endif
        endif
        
        call cpu_time(t0)

        if (model_name == "IEM") then
            call iem(np, npd, ncomp, wt, f, omdt)
        elseif (model_name == "MC") then
            call mcurl(np, npd, ncomp, wt, f, omdt)
        elseif (model_name == "KerM") then
            call kermix(np, npd, ncomp, wt, f, omdt, f(:,1))
        elseif (model_name == "EMST") then
            ! initialize if necessary
            if (omt == 0.0) then
                ! emst (mode=1) to initialize the state (age) variable state(1:np)
                call emst(1,np,ncomp,f(1:np,1:ncomp),state,wt,omdt,fscale,cvars,info)
                if( info /= 0 ) then
                    write(0,*)'emst initialization failed'
                    stop
                endif
            endif
            ! perform mixing
            call emst(2,np,ncomp,f(1:np,1:ncomp),state,wt,omdt,fscale,cvars,info)
            if( info /= 0 ) then
                write(*,*) 'emst failed with info = ', info
                stop
            endif
        endif

        call cpu_time(t1)
        timecost = timecost + (t1-t0)

        omt = omt + omdt
	end do
    
    888 continue
    write(*,*) 'time cost = ', timecost
    write(*,*) ''

    ! save cost
    filename = 'data/'//trim(adjustl(casename))//'_fortran_'// &
                        trim(adjustl(model_name))//'_costs.txt'
    open(301, file=filename, access='APPEND')
        write(301,'(I10 e15.7)') np, timecost
    close(301)
    
	stop

contains
    !-------------------------------------------------------
    ! Compute the scaled variances, vf.
    ! vf is defined as the trace of the covariance of the scaled compositions.
	subroutine compute_var(vf)
        real(kind(1.e0)), intent(out):: vf
        real(kind(1.e0)):: sumwt
        real(kind(1.e0)), allocatable:: mean(:), var(:)
        integer:: j, kk
	
        allocate( mean(ncomp), var(ncomp) )

        ! check that scales are strictly positive
        if( minval(fscale(1:ncomp)) <= 0. ) then
            write(0,*)'mixing_test: non-positive fscale in compute_var'
            stop
        endif

        ! check that weights are non-negative
        if( minval(wt(1:np)) <= 0. ) then
            write(0,*)'mixing_test: negative weight in compute_var'
            stop
        endif

        ! form means and variances (of un-scaled compositions)
        sumwt = sum(wt(1:np))
        
        mean(1:ncomp) = 0
        var(1:ncomp) = 0
        do j=1,ncomp
            mean(j) = dot_product( wt(1:np), f(1:np,j) ) / sumwt
            var(j) = dot_product( wt(1:np), (f(1:np,j)-mean(j))**2 ) / sumwt
        enddo

        ! scaled variances
        vf = sum( var(1:ncomp)/fscale(1:ncomp)**2 )

	end subroutine compute_var

end program mixing_test