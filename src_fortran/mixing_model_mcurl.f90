    subroutine mcurl(np, npd, ncomp, wt, f, omdt)
    ! Inputs:
    !   np     - number of particles
    !   npd    - leading dimension of f array (npd >= np, means some particles are not used)
    !   ncomp  - number of composition that equals to ns+1, where ns is the number of species
    !   wt(i)  - particle_i weight that is considered proportional to the mass of particle_i
    !   f(i,j) - mass fraction of species j of particle_i, j<=ncomp-1
    !            specific sensible enthalpy of particle_i, j=ncomp
    !   omdt   - normalized mixing time for all compositions (species and enthalpies)
    !
    ! Outputs:
    !   f(i,j) - mass fraction of species j of particle_i, j<=ncomp-1
    !            specific sensible enthalpy of particle_i, j=ncomp
    !   wt(i)  - particle_i weight that is considered proportional to the mass of particle_i
    
    implicit none
    integer,          intent(in)    :: np, npd, ncomp
    real(kind(1.e0)), intent(in)    :: omdt
    real(kind(1.e0)), intent(inout) :: f(npd,ncomp), wt(np)

    real(kind(1.e0)) :: wtsum, wtave, wtmax
    real(kind(1.e0)) :: ur(2), fbar(ncomp), a
    integer          :: p, q, nmix, imix

    !  quick returns if no mixing required
    if( np <= 1 ) return;

    ! weights statistics
    wtsum = sum( wt )
    wtave = wtsum / np
    wtmax = maxval( wt )
    
    ! get number of mixing pairs
    call random_number( ur(1) )
    nmix = int( 1.5*omdt*np + ur(1) )
    
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

        call random_number(a)
        fbar(1:ncomp) = (f(p,1:ncomp)*wt(p) + f(q,1:ncomp)*wt(q)) / (wt(p) + wt(q))
        
        f(p,1:ncomp) = f(p,1:ncomp) - a * (f(p,1:ncomp) - fbar(1:ncomp))
        f(q,1:ncomp) = f(q,1:ncomp) - a * (f(q,1:ncomp) - fbar(1:ncomp))
    enddo

    return
    end subroutine mcurl
