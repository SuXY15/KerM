    subroutine iem(np, npd, ncomp, wt, f, omdt)
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
    real(kind(1.e0)) :: wtsum
    real(kind(1.e0)) :: f_mean(ncomp)
    integer:: i,j

    !  quick returns if no mixing required
    if( np <= 1 ) return;
    
    wtsum = sum( wt(1:np) )
    do j=1,ncomp
        f_mean(j) = dot_product( f(1:np,j), wt(1:np) ) / wtsum
    enddo

    do i=1,np
        f(i,:) = f(i,:) - 0.5 * omdt * (f(i,:) - f_mean(:))
    enddo

    return
    end subroutine iem
