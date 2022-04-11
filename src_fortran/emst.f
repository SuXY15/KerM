      subroutine emst(mode,np,nc,f,state,wt,omdt,fscale,cvars,info)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  Fortran 90 subroutine to implement the EMST mixing model.

c  by Z.Ren and S.B. Pope, Cornell University, 2002.
c  Most lower level routines originally written by S. Subramaniam.
c  Ref: S.Subramaniam and S.B. Pope, Combust. Flame, 115: 487-514 (1998).

c  The EMST mixing model is applied to an ensemble of particles for a non-dimensional
c  time omdt. The ensemble is composed of np particles. Each particle has nc composition
c  variables (e.g., species mass fractions and enthalpy).  The array element f(i,j) contains
c  the j-th composition of the i-th particle. In addtion to the composition f(i,:),
c  the i-th particle has state (or age) variable state(i), and the numerical weight wt(i).

c Input:
c   mode  =1 -the state variables of all particles are initialized. No mixing is performed.
c         =2 -mixing is performed and the state variables are incremented.
c   np    -The number of particles in the ensemble. >=1 (integer)
c   nc    -The number of composition variables.     >=1 (integer)
c   f     -particle composition: f(i,j) is the j-th composition of the i-th particle
c   state -State (or age) variables.
c   wt    -particle weights. wt(i)>0 is the numerical weight of the i-th particle.
c   omdt  -normalized mixing time - see explanation below
c   fscale -scale factors for the compositions (fscale(j)>0) - see explanation below
c   cvar  -control variables for expert users.  Set cvar(1:6)=0. to obtain default values.

c Output:
c   f     -particle composition: unchanged for mode=1; value after mixing for mode=2.
c   state -age properties for particles: initialized for mode=1; incremented for mode=2.
c   info  = 0 for successful execution
c         < 0 for error in input

      use emst_subs, only: emst_init, emst_scl

      implicit none
      integer, intent(in) :: mode, np, nc
      real(kind(1.e0)), intent(in) :: wt(np), fscale(nc), omdt, cvars(6)
      real(kind(1.e0)), intent(inout) :: f(np,nc), state(np)
      integer, intent(out) :: info
    
c  Further explanations:
c
c  fscale:
c        The scale factors are used to define the scaled composition variables:
c        g(i,j)=f(i,j)/fscale(j),  and hence fscale(j) must be strictly positive.
c        The EMST is constructed based on the values of g.   Hence a relatively large value
c        of fscale(j) makes the EMST independent of f(:,j). 
c        Example of usage: set fscale=0.1 for all species, and set fscale=1.e16 for enthalpy, so
c        that the EMST is constructed solely from the species.  (Once the EMST has been formed,
c        the mixing does not depend on fscale(:).)  For numerical stability, fscale(:) should be
c        chosen so that the maximum range of the scaled variables is of order unity.
c
c  omdt:
c        In conventional usage, omdt = dt * C_phi * eps/k, where dt is the time step, C_phi is
c        the mixing model constant (usually C_phi=2) and k and eps are the turbulent kinetic
c        energy and its dissiapation rate.  In the case of a single composition (nc=1), the
c        composition variance decreases by the factor exp(-omdt) due to mixing. In general, 
c        the "variance function" (VF) is defined as the trace of the covariance matrix of the scaled
c        compositions, g.  Due to the mixing performed by emst, in expectation, the variance 
c	   function decreases by the factor exp(-omdt).
c
c  cvars:
c        All but expert users should set cvars(1:6)=0. which results in default values of various
c        control variables. The control variables are associated with the following local variable
c        names, and have the following meanings:
c
c  cvars(1)=icheck -  for cvars(1)=1. the input particle properties (wt, fscale,and state) are checked,
c                     otherwise (by default) they are not.  This can be used for debugging.
c
c  cvars(2)=Max_Ave - in some extreme circumstances (e.g., very small variance) mixing is performed
c                     by IEM instead of by EMST.  By default (cvars(2)=0.), the criterion for using
c                     IEM is based on the maximum range of the compositons.  Set cvars(2)=1. to use
c                     the average range instead.
c
c  cvars(3)=sdev_small - The threshold for performing IEM instead of EMST. (Setting cvars(3)=0. results
c                     in the default value sdev_small=1.e-4.)  IEM is used if: range < 4*sdev_samll,
c                     where  range=max_j( max_i( g(i,j))-min_i( g(i,j)))  for Max_Ave=0
c                     and    range=ave_j(max_i( g(i,j))-min_i( g(i,j)))   for Max_Ave=1.
c
c  cvars(4)=omdtmin - if omdt < omdtmin, then a rejection method is used. This means that more or
c                     less than the indicated amount of mixing is performed, but, in expectation,
c                     the correct amount is performed. Setting cvars(4)=0. results in the default
c                     value omdtmin=1.e-6.  Values in the range 1.e-7<omdtmin<1.e-2 can be set
c                     in cvars(4)
c  cvars(5)= fphthr - used in mixemst
c  cvars(6)= rphthr - used in mixemst
c
c  Notes:
c 1/ All reals in arguments are real(kind(1.e0)).
c 2/ In the calling routine, if f is dimensioned f(npd,ncd) with npd>np and ncd>=nc,
c    then the array section f(1:np,1:nc) should be used in the calling sequence.
c 3/ All routines are fixed format fortran 90.
c 4/ For mode=1, the only arguments references are: mode, np, state and info.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       local variables
        integer i,icheck
        real g(np,nc),cons(4)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   info=0
c check mode
         if(mode /= 1 .and. mode /=2) then
            write(0,*) 'emst: error input for mode= ', mode
            info=-1
	      return
         endif
         
c  check np
	   if( np < 0 ) then
            write(0,*) 'emst: error input for np= ', np
            info=-2	! must be non-negative
	      return
	   elseif( np == 0 ) then
	      return	! null ensemble
	   elseif( np==1 .and. mode==2 ) then
	      return	! singleton - no mixing
	   endif

c  get all the coeffient constants and model parameters for EMST-----------------

             cons(1) = 0.1666
             cons(2) = 0.1667
             cons(3) = 0.0178
             cons(4) = 0.3157

c-----perform mode 1, initialize the age properties for the particles in the ensemble.

        if (mode == 1) then
            call emst_init( np, cons(1), cons(2), cons(3), cons(4),
     1                    state )
            return
        endif


c-----perform mode 2. ---------------------------------------------------

      if (omdt == 0.0) then
	      return	! no mixing to perform
	    elseif( omdt < 0.0 ) then
	      write(0,*) 'emst: error input for omdt= ', omdt
            info=-7
	      return
	    endif

	   if( nc < 1 ) then
	      write(0,*) 'emst: error input for nc= ', nc
            info=-3
	      return
	   endif

c-------if icheck =1, check the inputs of wt,state,and fscale.
         icheck=nint(cvars(1))
         if(icheck == 1) then
	     if( minval(wt(1:np)) < 0. ) then
	         write(0,*) 'emst: negative particle weight'
	         info=-4
               return
           endif

	     if( maxval(state(1:np)) > cons(4) .or.
	1         minval(state(1:np)) < -cons(2) ) then
                
                write(0,*)'emst: particle age property out of range'
	          info=-5
                return
           endif

           if( minval(fscale(1:nc)) <= 0. ) then
              write(0,*) 'emst: non-positive value of fscale'
	        info=-6
              return
           endif
         endif

c------------------- scale the particle compositions -----------------------
c 
         do i=1,nc
           g(1:np,i)=f(1:np,i)/fscale(i)
         enddo

c  mix scaled compositions

	   call emst_scl(np,nc,g,state,wt,omdt,cvars,cons)

c-----------unscale the compositions before returning ----------------------
         do i=1,nc
            f(1:np,i)=g(1:np,i)*fscale(i)
         enddo

         return
         end subroutine emst