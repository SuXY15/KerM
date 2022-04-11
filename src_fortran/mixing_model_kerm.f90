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
    real(kind(1.e0)) :: mean_z, var_z, dvar_z, CDF(np)
    real(kind(1.e0)) :: xpq, dpq, fpq
    real(kind(1.e0)) :: wtsum, wtave, wtmax
    real(kind(1.e0)) :: ur(2), fbar(ncomp), a
    integer          :: i, p, q, nmix, imix

    !  quick returns if no mixing required
    if( np <= 1 ) return;

    ! settings
    epsilon = 1e-24;
    sigma_k = 0.25; ! the only model parameter, KerM(0.1)->EMST, KerM(1)->MC
    sigma_k = max(sigma_k, 1.0/np);
    var_k = sigma_k**2;

    ! weights statistics
    wtsum = sum( wt )
    wtave = wtsum / np
    wtmax = maxval( wt )
    
    ! calculate z mean and z var
    do i = 1,np
        mean_z = mean_z + wt(i) * z(i);
    enddo
    mean_z = mean_z / wtsum;

    do i = 1,np
        var_z = var_z + wt(i) * (z(i)-mean_z)**2;
    enddo
    var_z = var_z / wtsum;

    ! get CDF of mixture fraction
    if (np<=100) then
        call quickSort_CDF(z, CDF, np)
    else if(np<=400) then
        call bucketSort_CDF(z, CDF, np, 20)
    else
        call bucketSort_CDF(z, CDF, np, 50)
    endif

    dvar_z = 0;
    ! calculate coeffs of kernel mixing
    if (np <= 50) then
        ! estimate with all possible pairs
        do p = 1,np
            do q = 1,np
                if( p == q .or. wt(p) <= 0. .or. wt(q) <= 0. ) cycle
                xpq = abs(CDF(p) - CDF(q))
                dpq = abs(z(p) - z(q))
                fpq = exp(-xpq**2/4/var_k)
                dvar_z = dvar_z + wt(p)*wt(q)/(wt(p)+wt(q)) * fpq * dpq**2 / wtsum / np
            enddo
        enddo
        coeffs = (var_z + epsilon) / (dvar_z + epsilon)
    else
        ! estimate with np pairs
        do i = 1,np
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
    endif
    coeffs = (var_z + epsilon) / (dvar_z + epsilon)
    
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
    real(kind(1.e0)) :: mean_z, var_z, dvar_z, CDF(np)
    real(kind(1.e0)) :: xpq, dpq, fpq
    real(kind(1.e0)) :: wtsum, wtave, wtmax
    real(kind(1.e0)) :: ur(2), fbar(ncomp), theta(ncomp), a
    integer          :: i, j
    integer          :: p, q, nmix, imix
    real(kind(1.e0)) :: MassH( np, ncomp)
    !   MassH(i,j) j<=ncomp-1 - mass_j ( mass of species_j ) of particle_i, 
    !              j= ncomp   - sensible enthalpy of particle_i,            
    
    !  quick returns if no mixing required
    if( np <= 1 ) return;

    ! settings
    epsilon = 1e-24;
    sigma_k = 0.25; ! the only model parameter, KerM(0.1)->EMST, KerM(1)->MC
    sigma_k = max(sigma_k, 1.0/np);
    var_k = sigma_k**2;

    ! weights statistics
    wtsum = sum( wt(1:np) )
    wtave = wtsum / np
    wtmax = maxval( wt )
    
    ! re-normalize mass fractions
    do i = 1, np
        f(i,1:ncomp-1) = f(i,1:ncomp-1)/sum(f(i,1:ncomp-1))
    end do

    ! Get mass of compositions from compositions f(i,j) and mass weights wt(i)
    do j = 1, ncomp
        MassH(1:np, j ) = f(1:np, j) * wt(1:np)
    end do

    ! calculate z mean and z var
    do i = 1,np
        mean_z = mean_z + wt(i) * z(i);
    enddo
    mean_z = mean_z / wtsum;

    do i = 1,np
        var_z = var_z + wt(i) * (z(i)-mean_z)**2;
    enddo
    var_z = var_z / wtsum;

    ! get CDF of mixture fraction
    if (np<=100) then
        call quickSort_CDF(z, CDF, np)
    else if(np<=400) then
        call bucketSort_CDF(z, CDF, np, 20)
    else
        call bucketSort_CDF(z, CDF, np, 50)
    endif

    ! calculate coeffs of kernel mixing
    dvar_z = 0;
    do i = 1,np
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
        dvar_z = dvar_z + wt(p)*wt(q)/(wt(p)+wt(q)) * fpq * dpq**2
    enddo
    coeffs = (var_z + epsilon) / (dvar_z + epsilon)

    ! get number of mixing pairs
    call random_number( ur(1) )
    nmix = int( 1.5 * maxval(omdt) * np * coeffs + ur(1) )
    
    ! get mixing ratio coeffs theta
    call theta_calculate(ncomp, omdt, theta)

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
            
            MassH(p, 1:ncomp) = MassH(p, 1:ncomp) - a * theta(1:ncomp) &
                                        * ( MassH(p, 1:ncomp) - wt(p) * fbar(1:ncomp) )
            MassH(q, 1:ncomp) = MassH(q, 1:ncomp) - a * theta(1:ncomp) &
                                        * ( MassH(q, 1:ncomp) - wt(q) * fbar(1:ncomp) )

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
    ! bucket sort for CDF information
    ! ==================================================
    subroutine bucketSort_CDF(A, cdfA, N, Nbin)
    
    implicit none
    integer,          intent(in)    :: N, Nbin
    real(kind(1.e0)), intent(in)    :: A(N)
    real(kind(1.e0)), intent(inout) :: cdfA(N)
    real(kind(1.e0)) :: minA, maxA
    real(kind(1.e0)) :: minB(Nbin), maxB(Nbin)
    integer          :: countB(Nbin), accumB(Nbin)
    integer          :: i, j

    minA = minval(A) - 1e-10
    maxA = maxval(A) + 1e-10

    do j=1,Nbin
      minB(j) = minA + (maxA-minA) / Nbin * (j-1)
      maxB(j) = minA + (maxA-minA) / Nbin * j
    end do

    countB(:) = 0
    do i=1,N
       j = min(int((A(i)-minA)/(maxA-minA)*Nbin + 1), Nbin)
       countB(j) = countB(j) + 1
    end do

    accumB(1) = 0
    do j=2,Nbin
       accumB(j) = accumB(j-1) + countB(j-1)
    end do

    do i=1,N
       j = min(int((A(i)-minA)/(maxA-minA)*Nbin + 1), Nbin)
       cdfA(i) = ((A(i)-minB(j))/(maxB(j)-minB(j))*countB(j) + accumB(j)) / N
    end do
    end subroutine bucketSort_CDF



    ! ==================================================
    ! quick sort for CDF information
    ! ==================================================
    subroutine quickSort_CDF(A, cdfA, N)

    implicit none
    integer,          intent(in)    :: N       ! size of array
    real(kind(1.e0)), intent(in)    :: A(N)    ! Array to be sorted (not changed)
    real(kind(1.e0)), intent(inout) :: cdfA(N) ! Array CDF to be calculated
    integer          :: i
    integer          :: indx(N)
    real(kind(1.e0)) :: tmpA(N)

    tmpA(:) = A(:)
    
    do i=1,N
       indx(i) = i
    end do

    call quickSort(1, N, tmpA, indx)

    do i=1,N
       cdfA(indx(i)) = real(i) / N
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
    end do
    
    end subroutine theta_calculate