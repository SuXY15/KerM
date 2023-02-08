    ! ==================================================
    ! Kernel Mixing Model by Xingyu Su, 2022/03/30.
    ! ==================================================
    subroutine kermix(np, npd, ncomp, wt, f, omdt, z)
    ! Mixing with single time scale (no differential diffusion)
    ! For more details, please refer to https://github.com/SuXY15/KerM 
    ! 
    ! Inputs:
    !   np     - number of particles
    !   npd    - leading dimension of f array (npd >= np, means some particles are not used)
    !   ncomp  - number of composition that equals to ns+1, where ns is the number of species
    !   wt(i)  - particle_i weight that is considered proportional to the mass of particle_i
    !   f(i,j) - mass fraction of species j of particle_i, j<=ncomp-1
    !            specific sensible enthalpy of particle_i, j=ncomp
    !   omdt   - normalized mixing time for all compositions (species and enthalpies)
    !   z(i)   - reference variable of particle_i (i.e., mixture fraction)
    !
    ! Outputs:
    !   f(i,j) - mass fraction of species j of particle_i, j<=ncomp-1
    !            specific sensible enthalpy of particle_i, j=ncomp
    !   wt(i)  - particle_i weight that is considered proportional to the mass of particle_i
    
    implicit none
    integer,          intent(in)    :: np, npd, ncomp
    real(kind(1.e0)), intent(in)    :: omdt, z(np)
    real(kind(1.e0)), intent(inout) :: f(npd,ncomp), wt(np)
    real(kind(1.e0)) :: sigma_k, var_k, coeffs, epsilon
    real(kind(1.e0)) :: mean_z, var_z, square_z, CDF(np)
    real(kind(1.e0)) :: xpq, dpq, fpq
    real(kind(1.e0)) :: wtsum, wtave, wtmax
    real(kind(1.e0)) :: ur(2), fbar(ncomp), a
    real(kind(1.e0)), allocatable :: centerB(:), countB(:), accumB(:)
    integer          :: i, p, q, nmix, imix, Nbin
        
    !  quick returns if no mixing required
    if( np <= 1 ) return;

    ! settings
    epsilon = 1e-24
    sigma_k = 0.25 ! the only model parameter, KerM(0.1)->EMST, KerM(1)->MC
    sigma_k = max(sigma_k, 1.0/np)
    var_k = sigma_k**2
    
    ! calculate z mean and z var and get weight statistics
    wtsum = 0.0
    wtmax = wt(1)
    mean_z = 0.0
    square_z = 0.0
    do i = 1,np ! single loop is the fastest
        wtsum = wtsum + wt(i)
        wtmax = max(wtmax, wt(i))
        mean_z = mean_z + wt(i) * z(i)
        square_z = square_z + wt(i) * z(i) * z(i)
    enddo
    mean_z = mean_z / wtsum
    square_z = square_z / wtsum
    var_z = square_z - mean_z**2
    wtave = wtsum / np
        
    ! get CDF of reference variable and emstimate coefficients
    if (np < 100) then
        call quickSort_CDF(z, wt, CDF, np)
        !* estimate with all possible pairs, O(N^2) complexity
        call get_coeffs(0, np, z, wt, CDF, var_z, var_k, coeffs) 
    else 
        Nbin = max(10, min(floor(sqrt(real(np))/2), 100))
        allocate( centerB(Nbin), countB(Nbin), accumB(Nbin) )
        call bucketSort_CDF(z, wt, CDF, np, Nbin, centerB, countB, accumB)
        accumB(:) = accumB(:) + 0.5*countB(:)
        
        !* estimate with np pairs, O(N) complexity, accurate enough
        call get_coeffs(1, np, z, wt, CDF, var_z, var_k, coeffs)  
        
        ! !* estimate with the bins PDF, O(Nbin^2) complexity, about 3x faster
        ! call get_coeffs(0, NBin, centerB, countB, accumB, var_z, var_k, coeffs)
    endif
            
    ! write(*,*) 'kermix: coeffs=', coeffs

    ! get number of mixing pairs
    call random_number( ur(1) )
    nmix = int( 1.5*omdt*np*coeffs + ur(1) )
    
    ! loop over particle-pair interactions
    do imix = 1, nmix
        ! select p:  marginal prob ~ wtave + w(p)
        do
            call random_number( ur )
            p = min( np , 1 + int( ur(1) * np ) )
            if( wt(p)+wtave > ur(2) * (wtmax+wtave) ) exit ! accept p
        end do
        ! select q:  conditional prob ~ w(p) + w(q)
        do
            call random_number( ur )
            q = min( np , 1 + int( ur(1) * np ) )
            if( wt(p)+wt(q) > ur(2) * (wt(p)+wtmax) ) exit ! accept q
        end do
        if( p == q .or. wt(p) <= 0. .or. wt(q) <= 0. ) cycle

        xpq = abs(CDF(p) - CDF(q))
        dpq = abs(z(p) - z(q))
        fpq = exp(-xpq**2/4/var_k)

        call random_number(ur)
        if ( ur(1) <= fpq ) then
            a = ur(2)
            fbar(1:ncomp) = (f(p,1:ncomp)*wt(p) + f(q,1:ncomp)*wt(q)) / (wt(p) + wt(q))
            
            f(p,1:ncomp) = f(p,1:ncomp) - a * (f(p,1:ncomp) - fbar(1:ncomp))
            f(q,1:ncomp) = f(q,1:ncomp) - a * (f(q,1:ncomp) - fbar(1:ncomp))
        endif
    enddo

    return
    end subroutine kermix



    subroutine get_coeffs(mode, np, z, wt, CDF, var_z, var_k, coeffs)
        implicit none
        integer,          intent(in)  :: mode     ! =0 for all possible pairs
                                                  ! =1 for np pairs
                                                  ! =2 for np*np pairs
        integer,          intent(in)  :: np       ! number of elements
        real(kind(1.e0)), intent(in)  :: z(np)    ! value of reference variable
        real(kind(1.e0)), intent(in)  :: wt(np)   ! weight of each element
        real(kind(1.e0)), intent(in)  :: CDF(np)  ! CDF of z
        real(kind(1.e0)), intent(in)  :: var_z    ! variance of z
        real(kind(1.e0)), intent(in)  :: var_k    ! sigma_k ** 2
        real(kind(1.e0)), intent(out) :: coeffs   ! coeffs to be calculated
        real(kind(1.e0)) :: xpq, dpq, fpq, prob, ur(2)
        real(kind(1.e0)) :: wtsum, wtave, wtmax, dvar_z, epsilon
        integer          :: i, p, q, nmix
        
        epsilon = 1e-24
        wtsum = sum( wt )
        wtmax = maxval( wt )
        wtave = wtsum / np
        if (mode == 0) then
            do p = 1,np
                do q = 1,np
                    if( p == q .or. wt(p) <= 0. .or. wt(q) <= 0. ) cycle
                    xpq = abs(CDF(p) - CDF(q))
                    dpq = abs(z(p) - z(q))
                    fpq = exp(-xpq**2/4/var_k)
                    prob = (wt(p)+wt(q))/(2.*wtsum) * fpq
                    dvar_z = dvar_z + prob * wt(p)*wt(q)/(wt(p)+wt(q)) * dpq**2 / wtsum
                enddo
            enddo
        else
            if (mode == 1) then
                nmix = np
            else
                nmix = np * np
            end if

            do i = 1,nmix
                ! select p:  marginal prob ~ wtave + w(p)
                do
                    call random_number( ur )
                    p = min( np , 1 + int( ur(1) * np ) )
                    if( wt(p)+wtave > ur(2) * (wtmax+wtave) ) exit ! accept p
                end do
                ! select q:  conditional prob ~ w(p) + w(q)
                do
                    call random_number( ur )
                    q = min( np , 1 + int( ur(1) * np ) )
                    if( wt(p)+wt(q) > ur(2) * (wt(p)+wtmax) ) exit ! accept q
                end do
                if( p == q .or. wt(p) <= 0. .or. wt(q) <= 0. ) cycle
                xpq = abs(CDF(p) - CDF(q))
                dpq = abs(z(p) - z(q))
                fpq = exp(-xpq**2/4/var_k)
                dvar_z = dvar_z + wt(p)*wt(q)/(wt(p)+wt(q)) * fpq * dpq**2 / wtsum
            enddo
            dvar_z = dvar_z * np / nmix
        endif
        coeffs = (var_z + epsilon) / (dvar_z + epsilon)
    end subroutine get_coeffs

    ! ==================================================
    ! bucket sort for CDF information
    ! ==================================================
    subroutine bucketSort_CDF(A, wtA, cdfA, N, Nbin, centerB, countB, accumB)
    
    implicit none
    integer,          intent(in)    :: N             ! size of array
    integer,          intent(in)    :: Nbin          ! number of bins / sections
    real(kind(1.e0)), intent(in)    :: A(N)          ! Array to be sorted (not changed)
    real(kind(1.e0)), intent(in)    :: wtA(N)        ! weight of each element (not changed)
    real(kind(1.e0)), intent(inout) :: cdfA(N)       ! Array CDF to be calculated
    real(kind(1.e0)), intent(inout) :: centerB(Nbin) ! centers of Bins
    real(kind(1.e0)), intent(inout) :: countB(Nbin)  ! counted weights of Bins
    real(kind(1.e0)), intent(inout) :: accumB(Nbin)  ! CDF of Bins' left boundary
    real(kind(1.e0)) :: wtsum
    real(kind(1.e0)) :: minA, maxA, stepA
    real(kind(1.e0)) :: minB(Nbin), maxB(Nbin)
    integer          :: i, j

    minA = minval(A) - 1e-10
    maxA = maxval(A) + 1e-10
    stepA = (maxA-minA) / Nbin

    do j=1,Nbin
        minB(j) = minA + stepA * (j-1)
        maxB(j) = minA + stepA * j
        centerB(j) = minA + stepA * real(j) - 0.5
    enddo

    wtsum = 0.0
    countB(:) = 0
    do i=1,N
        wtsum = wtsum + wtA(i)
        j = min(int((A(i)-minA)/stepA + 1), Nbin)
        countB(j) = countB(j) + wtA(i)
    enddo
    countB(:) = countB(:) / wtsum

    accumB(1) = 0
    do j=1,Nbin-1
        accumB(j+1) = accumB(j) + countB(j)
    enddo

    do i=1,N
        j = min(int((A(i)-minA)/stepA + 1), Nbin) 
        cdfA(i) = (A(i)-minB(j))/stepA * countB(j) + accumB(j)
    enddo
    
    end subroutine bucketSort_CDF



    ! ==================================================
    ! quick sort for CDF information
    ! ==================================================
    subroutine quickSort_CDF(A, wtA, cdfA, N)

    implicit none
    integer,          intent(in)    :: N       ! size of array
    real(kind(1.e0)), intent(in)    :: A(N)    ! Array to be sorted (not changed)
    real(kind(1.e0)), intent(in)    :: wtA(N)  ! weight of each element (not changed)
    real(kind(1.e0)), intent(inout) :: cdfA(N) ! Array CDF to be calculated
    integer          :: i
    integer          :: indx(N)
    real(kind(1.e0)) :: wtnow, wtsum, tmpA(N)

    tmpA(:) = A(:)
    
    do i=1,N
       indx(i) = i
    end do

    call quickSort(1, N, tmpA, indx)

    wtsum = sum(wtA)
    wtnow = 0.0
    do i=1,N
        wtnow = wtnow + wtA(indx(i)) / wtsum
        cdfA(indx(i)) = wtnow
    end do

        contains
        ! ==================================================
        ! quick sort with assistance array 
        recursive subroutine quickSort(l, r, A, nA)

        implicit none
        real(kind(1.e0)), intent(inout), dimension(:) :: A  ! Array to be sorted
        integer,          intent(inout), dimension(:) :: nA ! Assist array, i.e., index of A
        integer,          intent(in)                  :: l, r
        integer          :: i, j
        integer          :: ntmp1, ntmp2
        real(kind(1.e0)) :: tmp1, tmp2
        
        if (l>r) return;

        tmp2 = A(l); ntmp2 = nA(l);
        i = l; j = r;
        do while (i.ne.j)
            do while( (A(j)<=tmp2) .and. (i<j) )
                j = j-1;
            end do
            do while( (A(i)>=tmp2) .and. (i<j) )
                i = i+1;
            end do
            if (i<j) then
                 tmp1 =  A(i);  A(i) =  A(j);  A(j) =  tmp1;
                ntmp1 = nA(i); nA(i) = nA(j); nA(j) = ntmp1;
            endif
        end do
         A(l) = A(i);   A(i) =  tmp2;
        nA(l) = nA(i); nA(i) = ntmp2;
        call quickSort(l, i-1, A, nA);
        call quickSort(i+1, r, A, nA);
        
        end subroutine quickSort

    end subroutine quickSort_CDF



    ! ==================================================
    ! Kernel Mixing Model with differential diffusion by Xingyu Su, 2022/03/30.
    ! ==================================================
    subroutine kermix_dd(np, npd, ncomp, wt, f, omdt, z)
    ! Mixing with time scales of each composition
    ! For more details, please refer to https://github.com/SuXY15/KerM 
    ! 
    ! Inputs:
    !   np      - number of particles
    !   npd     - leading dimension of f array (npd >= np, means some particles are not used)
    !   ncomp   - number of composition that equals to ns+1, where ns is the number of species
    !   wt(i)   - particle_i weight that is considered proportional to the mass of particle_i
    !   f(i,j)  - mass fraction of species j of particle_i, j<=ncomp-1
    !             specific sensible enthalpy of particle_i, j=ncomp
    !   omdt(j) - normalized mixing times for each composition (species and enthalpies)
    !   z(i)    - reference variable of particle_i (i.e., mixture fraction)
    !
    ! Outputs:
    !   f(i,j)  - mass fraction of species j of particle_i, j<=ncomp-1
    !             specific sensible enthalpy of particle_i, j=ncomp
    !   wt(i)   - particle_i weight that is considered proportional to the mass of particle_i
    
    implicit none
    integer,          intent(in)    :: np, npd, ncomp
    real(kind(1.e0)), intent(in)    :: omdt(ncomp), z(np)
    real(kind(1.e0)), intent(inout) :: f(npd,ncomp), wt(np)

    real(kind(1.e0)) :: sigma_k, var_k, coeffs, epsilon
    real(kind(1.e0)) :: mean_z, var_z, square_z, CDF(np)
    real(kind(1.e0)) :: xpq, dpq, fpq
    real(kind(1.e0)) :: wtsum, wtave, wtmax
    real(kind(1.e0)) :: ur(2), fbar(ncomp), theta(ncomp), a
    real(kind(1.e0)), allocatable :: centerB(:), countB(:), accumB(:)
    integer          :: i, j
    integer          :: p, q, nmix, imix, Nbin
    real(kind(1.e0)) :: MassH( np, ncomp)
    !   MassH(i,j) j<=ncomp-1 - mass_j ( mass of species_j ) of particle_i, 
    !              j= ncomp   - sensible enthalpy of particle_i,            
    
    !  quick returns if no mixing required
    if( np <= 1 ) return;

    ! settings
    epsilon = 1e-24
    sigma_k = 0.25  ! the only model parameter, KerM(0.1)->EMST, KerM(1)->MC
    sigma_k = max(sigma_k, 1.0/np)
    var_k = sigma_k**2

    ! calculate z mean and z var and get weight statistics
    wtsum = 0.0
    wtmax = wt(1)
    mean_z = 0.0
    square_z = 0.0
    do i = 1,np ! single loop is the fastest
        wtsum = wtsum + wt(i)
        wtmax = max(wtmax, wt(i))
        mean_z = mean_z + wt(i) * z(i)
        square_z = square_z + wt(i) * z(i) * z(i)
    enddo
    mean_z = mean_z / wtsum
    square_z = square_z / wtsum
    var_z = square_z - mean_z**2
    wtave = wtsum / np

    ! get CDF of reference variable and emstimate coefficients
    if (np < 100) then
        call quickSort_CDF(z, wt, CDF, np)
        !* estimate with all possible pairs, O(N^2) complexity
        call get_coeffs(0, np, z, wt, CDF, var_z, var_k, coeffs) 
    else 
        Nbin = max(10, min(floor(sqrt(real(np))/2), 100))
        allocate( centerB(Nbin), countB(Nbin), accumB(Nbin) )
        call bucketSort_CDF(z, wt, CDF, np, Nbin, centerB, countB, accumB)
        
        !* estimate with np pairs, O(N) complexity, accurate enough
        call get_coeffs(1, np, z, wt, CDF, var_z, var_k, coeffs)  
        
        ! !* estimate with the bins PDF, O(Nbin^2) complexity, about 3x faster
        ! call get_coeffs(0, NBin, centerB, countB, accumB, var_z, var_k, coeffs)
    endif

    ! above contents are the same with the one without differential diffusion 
    ! ==================================================

    ! get number of mixing pairs
    call random_number( ur(1) )
    nmix = int( 1.5 * maxval(omdt) * np * coeffs + ur(1) )
    
    ! get mixing ratio coeffs theta
    call theta_calculate(ncomp, omdt, theta)

    ! re-normalize mass fractions
    do i = 1, np
        f(i,1:ncomp-1) = f(i,1:ncomp-1)/sum(f(i,1:ncomp-1))
    end do

    ! Get mass of compositions from compositions f(i,j) and mass weights wt(i)
    do j = 1, ncomp
        MassH(1:np, j ) = f(1:np, j) * wt(1:np)
    end do

    ! loop over particle-pair interactions
    do imix = 1, nmix
        ! select p:  marginal prob ~ wtave + w(p)
        do
            call random_number( ur )
            p = min( np , 1 + int( ur(1) * np ) )
            if( wt(p)+wtave > ur(2) * (wtmax+wtave) ) exit ! accept p
        end do
        ! select q:  conditional prob ~ w(p) + w(q)
        do
            call random_number( ur )
            q = min( np , 1 + int( ur(1) * np ) )
            if( wt(p)+wt(q) > ur(2) * (wt(p)+wtmax) ) exit ! accept q
        end do
        if( p == q .or. wt(p) <= 0. .or. wt(q) <= 0. ) cycle

        xpq = abs(CDF(p) - CDF(q))
        dpq = abs(z(p) - z(q))
        fpq = exp(-xpq**2/4/var_k)

        call random_number(ur)
        if ( ur(1) <= fpq ) then
            a = ur(2)
            fbar(1:ncomp) = ( MassH(p,1:ncomp)+MassH(q,1:ncomp) ) / (wt(p) + wt(q))
            
            MassH(p,:) = MassH(p, :) - a * theta(:) * (MassH(p, :) - wt(p)*fbar(:))
            MassH(q,:) = MassH(q, :) - a * theta(:) * (MassH(q, :) - wt(q)*fbar(:))

            wt(p) = sum(MassH(p,1:ncomp-1))
            wt(q) = sum(MassH(q,1:ncomp-1))
        endif
    enddo

    ! Get compositions f(i,j) from mass of compositions MassH(i,j) and mass weights wt(i)
    do j = 1, ncomp
        f(1:np, j) = MassH(1:np, j) / wt(1:np)
    end do

    return
    end subroutine kermix_dd



    ! ==================================================
    ! Calculate thetas for all scalars according to T.Yang_CNF_2020
    ! ==================================================
    subroutine theta_calculate(ncomp, omdt, theta)
    ! Input:  
    !   ncomp       -  number of composition
    !   omdt(j)     -  normalized mixing time for scalar_j
    !
    ! Output:     
    !   theta(j)    -  controlling parameter for uniform distribution of a_j
    implicit none
    integer                       :: ncomp
    real(kind(1.e0)), intent(in)  :: omdt(ncomp)
    real(kind(1.e0)), intent(out) :: theta(ncomp)
    real(kind(1.e0)) :: omdt_relative(ncomp), omdt_max
    integer          :: j
    
    omdt_max = maxval( omdt(1:ncomp) )
    omdt_relative(1:ncomp) = omdt(1:ncomp) / omdt_max
    
    do j = 1, ncomp
        theta(j) = ( 3.0 - sqrt(9.0-8.0*omdt_relative(j)) )/2.0
    enddo
    
    end subroutine theta_calculate
