    subroutine mcmg(np, npd, ncomp, wt, f, omdt)
    ! The Mapping Closure (with Mapping of Gausian) mixing model
    ! Adopted from the MM-INTAS code from TNF workshop website:
    !   https://tnfworkshop.org/wp-content/uploads/2019/03/MM-INTAS.zip
    ! 
    ! Function :
    !    routine to advance scalars for one time step according to the
    !    mapping closure models for which the evolution of the scalar pdf
    !    is calculated by a time dependent mapping of a standard Gaussian
    !    Pope (1991), Valino, Ros & Dopazo (1991).
    !       - Pope (1991). Mapping closures for turbulent mixing and
    !         reaction.  Theoret. Comput. Fluid Dynamics, 2, 255-270.
    !       - Valino, L., Ros, J. and Dopazo, C. (1991). Monte Carlo
    !         implementation and analytic solution of an inert-scalar
    !         turbulent-mixing test problem using a mapping closure.
    !         Physics of Fluids A, 3, 2191-2198.
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
    
    USE blas_lapack_module, ONLY: dpbsv

    implicit none

!-----  input variables
    integer,          intent(in)    :: np, npd, ncomp
    real(kind(1.e0)), intent(in)    :: omdt

!-----  input/output variables
    real(kind(1.e0)), intent(inout) :: f(npd,ncomp), wt(np)

!-----  local variables
    INTEGER          :: i, ig, igadd, k, l
    real    (kind=8) :: hcphomdt, inc_wt, omega, wtsum, &
                        fold, fnew, ffold, ffnew, &
                        etai, etaph, etap, gph, fac

    integer         , ALLOCATABLE :: isp(:)                 ! (np)
    real    (kind=8), ALLOCATABLE :: cdf(:)                 ! (np+1)
    real    (kind=8), ALLOCATABLE :: eta(:)                 ! (np+1)
    real    (kind=8), ALLOCATABLE :: phi(:)                 ! (np)
    real    (kind=8), ALLOCATABLE :: bph(:)                 ! (np-1)

!-----  variables needed for call of LAPACK routine DPBSV
    integer         , PARAMETER :: kd = 1, nrhs = 1, ldab = kd+1

    integer          :: info
    real(kind=8), ALLOCATABLE :: ab(:,:)                ! (ldab,np)
    character(len=1) :: uplo

    real    (kind=8), PARAMETER :: rtpi = 0.39894228d0
    real    (kind=8), PARAMETER :: phimin = -1.E30, phimax = 1.E30
    real    (kind=8), PARAMETER :: tiny = 1.E-30, small = 1.E-10

    !  quick returns if no mixing required
    if( np <= 1 ) return;

    !--- set variables needed for DSBEV
    uplo = "U"

    hcphomdt = 0.5 * omdt

    ! weights statistics
    wtsum = sum( wt )

    ALLOCATE(isp(np))
    ALLOCATE(cdf(np+1))
    ALLOCATE(eta(np+1))
    ALLOCATE(phi(np))
    ALLOCATE(bph(np-1))
    ALLOCATE(ab(ldab,np))
    
    !--- loop over independent scalars
    SCALAR_LOOP: &
    DO k = 1, ncomp

        !--- calculate the mean and variance of f
        fold = sum( wt(:) * f(1:np,k) ) / wtsum
        ffold = sum( wt(:) * f(1:np,k) * f(1:np,k) ) / wtsum - fold**2

        !--- skip mixing if variance is very small
        IF (ffold < small) THEN
            CYCLE SCALAR_LOOP
        END IF

        !--- order particles on scalar value
        DO i = 1, np
            isp(i) = i
            phi(i) = f(i,k)
        END DO
        CALL sortip(np, phi, isp)

        !--- probabilities
        inc_wt = 0.
        DO i = 1, np-1
            ig     = isp(i)
            inc_wt = inc_wt + wt(ig)
            cdf(i) = inc_wt / wtsum
        END DO

        !--- to construct eta(1)
        cdf(np) = 0.5 * cdf(1)

        !--- to construct eta(np)
        cdf(np+1) = 0.5 * (cdf(np-1) + 1.)

        !--- mapped particles eta(1.5), eta(2.5),.... eta(1), eta(np)
        CALL geterf1(np+1, cdf , eta)

        !--- construct matrix coefficients

        !--- B_1+1/2
        etai   = eta(np)
        etaph  = eta(1)
        etap   = 0.5 * (eta(1) + eta(2))
        gph    = rtpi * exp(-0.5 * etaph**2)
        bph(1) = dble(np) * gph / (etap - etai)

        !--- B_i+1/2
        DO i = 2, np-2
            etai   = 0.5 * (eta(i-1) + eta(i  ))
            etaph  = eta(i)
            etap   = 0.5 * (eta(i  ) + eta(i+1))
            gph    = rtpi * exp(-0.5 * etaph**2)
            bph(i) = dble(np) * gph / (etap - etai)
        END DO

        !--- B_npt-1/2
        etai       = 0.5 * (eta(np-2) + eta(np-1))
        etaph      = eta(np-1)
        etap       = eta(np+1)
        gph        = rtpi * exp(-0.5 * etaph**2)
        bph(np-1) = dble(np) * gph / (etap - etai)

        !--- construct matrix (I-dt*A)
        ab(2,1) = 1. + bph(1) * hcphomdt
        DO i = 2, np-1
            ab(2,i) = 1. + (bph(i-1) + bph(i))*hcphomdt
        END DO
        ab(2,np) = 1. + bph(np-1)*hcphomdt

        DO i = 1, np-1
            ab(1,i+1) = -bph(i)*hcphomdt
        END DO

        !--- solve system
        CALL dpbsv(uplo, np, kd, nrhs, ab, ldab, phi, np, info)
        IF (info /= 0) THEN
            WRITE (*,*) "Error in subroutine MCMG"
            WRITE (*,*) "On return from DPBSV: INFO = ", info
            STOP 'MCMG: see output file.'
        END IF

        !--- new scalar values
        DO i = 1, np
            f(isp(i),k) = phi(i)
        END DO

        !--- calculate the mean and variance of f
        fnew = sum( wt(:) * f(1:np,k) ) / wtsum
        ffnew = sum( wt(:) * f(1:np,k) * f(1:np,k) ) / wtsum - fold**2

        !--- perform correction
        fac = sqrt(ffold/(ffnew + tiny)) * exp(-hcphomdt)
        
        DO i = 1, np
            f(i,k) = fold + (f(i,k) - fnew) * fac
            f(i,k) = MAX(MIN(f(i,k),phimax),phimin)
        END DO

    END DO SCALAR_LOOP

    DEALLOCATE(isp)
    DEALLOCATE(cdf)
    DEALLOCATE(eta)
    DEALLOCATE(phi)
    DEALLOCATE(bph)
    DEALLOCATE(ab)

    return
    end subroutine mcmg


!=======================================================================
!  subroutine geterf1( n, x, y )
!=======================================================================
!  method:
!    approximation to the inverse normal cumulative density function
!    according to the approximation by Abramowitz & Stegun formula
!    26.2.23 approximation error < 4.5E-04
!
!  limits :
!    x           [zero,one]
!    y           [-big,big] , y = -big for x < small
!                             y =  big for x > one - small
!
!    small       [2.8D-17]    chosen smallest number possible
!    big         [1000   ]    chosen arbitrary large number one order
!                             of magnitude larger than
!                             ERF1[ONE-SMALL+E]
!
!  input :
!    n           : number of function values
!    x           : probability
!
!  output :
!    f           : inverse error function of y
!----------------------------------------------------------------------
      SUBROUTINE geterf1(n, x, y)

      IMPLICIT NONE

!----- input variables
      integer         , INTENT(IN) :: n
      real    (kind=8), INTENT(IN) :: x(n)                    !(n)

!----- output variables
      real    (kind=8), INTENT(OUT) :: y(n)                   !(n)

!----- local variables
      integer          :: i
      real    (kind=8) :: xi, x1, x2, x3, yi, y1, t, t1, t2, f, f1, f2, s

      real    (kind=8), PARAMETER :: zero = 0.d0, &
                                     half = 0.5d0, &
                                     one = 1.d0, &
                                     two = 2.d0, &
                                     big = 1.D+03, &
                                     small = 2.8D-17

      real    (kind=8), PARAMETER :: c0 = 2.515517d0, &
                                     c1 = 0.802853d0, &
                                     c2 = 0.010328d0, &
                                     d1 = 1.432788d0, &
                                     d2 = 0.189269d0, &
                                     d3 = 0.001308d0

!-----  calculate inverse error function
      DO i = 1,n

!--------  endpoints clamped to -BIG, BIG
         IF     (x(i).LT.(zero + small)) THEN
            y(i) = -big
         ELSE IF (x(i).GT.(one  - small)) THEN
            y(i) =  big

!--------  normal calulation
         ELSE
            x1   = x(i)
            x2   = half - x1
            x3   = ABS(x2)
            xi   = half - x3
            s    = -SIGN(one, x2)

            t1   = LOG(xi)
            t2   = -two * t1
            t    = SQRT(t2)

            f1   = c0  + c1 * t + c2 * t2
            f2   = one + d1 * t + d2 * t2 + d3 * t*t2
            f    = f1 / f2

            y1   = t - f
            yi   = s * y1
            y(i) = yi
         END IF
      END DO

      END SUBROUTINE geterf1


!=======================================================================
!  subroutine inter1d
!=======================================================================
!  function:
!    interpolates in an array of function values and arguments
!    it is assumed that the arguments are sorted
!    extra endpoints can be passed using XIN, XMAX, FMIN, FMAX
!    values outsside [XMIN,XMAX] are clamped to FMIN or FMAX
!
!  input:
!    nin           : number of function values
!    xin           : list of arguments
!    xmin          : minimum argument value
!    xmax          : maximum argument value
!    fin           : list of function values
!    fmin          : minimum function value
!    fmax          : maximum function value
!    nout          : number of output values
!    xout          : list of arguments for interpolation
!
!  output:
!    fout          : new interpolated values
!-----------------------------------------------------------------------
      SUBROUTINE inter1d (nin, xin, xmin, xmax, fin, fmin, fmax, &
                          nout, xout, fout)

      IMPLICIT NONE

!-----  input variables
      integer         , INTENT(IN) :: nin, nout
      real    (kind=8), INTENT(IN) :: xmin, xmax, fmin, fmax
      real    (kind=8), INTENT(IN) :: xin(nin)                ! (nin)
      real    (kind=8), INTENT(IN) :: fin(nin)                ! (nin)
      real    (kind=8), INTENT(IN) :: xout(nout)              ! (nout)
                                     

!-----  output variables
      real    (kind=8), INTENT(INOUT) :: fout(nout)           ! (nout)

!-----  local variables
      integer          :: i, j
      real    (kind=8) :: f1, f2, fac1, fac2, x1, x2, dx, xmm, xxx, xouti

!=====  MAIN ACTION

!-----  limits
      xmm = xin(1  )
      xxx = xin(nin)

!-----  loop over points
      j = 2
      DO i = 1,nout
         xouti = xout(i)
         IF (xouti.LT.xmin) THEN
            fout(i) = fmin
         ELSE IF (xouti.GT.xmax) THEN
            fout(i) = fmax
         ELSE
            IF (xouti.LT.xmm) THEN
               x1 = xmin
               x2 = xmm
               f1 = fmin
               f2 = fin(1)
            ELSE IF (xouti.GT.xxx) THEN
               x1 = xxx
               x2 = xmax
               f1 = fin(nin)
               f2 = fmax
            ELSE
               DO
                  IF (xin(j).GT.xouti) THEN
                     x1 = xin(j-1)
                     x2 = xin(j  )
                     f1 = fin(j-1)
                     f2 = fin(j  )
                     EXIT
                  ELSE
                     j = j + 1
                  END IF
               END DO
               dx      =  x2    - x1
               fac1    = (x2    - xouti) / dx
               fac2    = (xouti - x1   ) / dx
               fout(i) = fac1 * f1 + fac2 * f2
            END IF
         END IF
      END DO
!-----  end of loop over points

      END SUBROUTINE inter1d


!=======================================================================
!  subroutine sortip
!=======================================================================
!  function:
!    routine to sort array RA into ascending order, and to rearrange
!    the array JB accordingly.  Modified version of numerical recipes
!    routine sort2.
!
!  input:
!    n      : number of array elements
!
!  input/output:
!    ra     : array to be sorted
!    jb     : integer array to be sorted accordingly
!-----------------------------------------------------------------------
      SUBROUTINE sortip(n, ra, jb)

      IMPLICIT NONE

!=====  VARIABLES

!-----  input variables
      integer         , INTENT(IN) :: n

!-----  input/output variables
      real    (kind=8), INTENT(INOUT) :: ra(n)                ! (n)
      integer         , INTENT(INOUT) :: jb(n)                ! (n)

!-----  local variables
      integer          :: i, j, ir, l, jjb
      real    (kind=8) :: rra

!=====  MAIN ACTION

      l  = n/2 + 1
      ir = n

      DO

         IF (l > 1) THEN
            l   = l - 1
            rra = ra(l)
            jjb = jb(l)
         ELSE
            rra    = ra(ir)
            jjb    = jb(ir)
            ra(ir) = ra(1 )
            jb(ir) = jb(1 )
            ir     = ir - 1

            IF (ir == 1) THEN
               ra(1) = rra
               jb(1) = jjb
               EXIT
            END IF
         END IF

         i = l
         j = l + l

         DO
            IF (j <= ir) THEN
               IF (j < ir) THEN
                  IF (ra(j) < ra(j+1)) j = j + 1
               END IF

               IF (rra < ra(j)) THEN
                  ra(i) = ra(j)
                  jb(i) = jb(j)
                  i     = j
                  j     = j  + j
               ELSE
                  j     = ir + 1
               END IF

            ELSE
               EXIT
            END IF
         END DO

         ra(i) = rra
         jb(i) = jjb

      END DO

      END SUBROUTINE sortip
