      module emst_subs

c  Fortran 90 subroutines which, together with emst, implement 
c  the EMST mixing model.

c  by Z.Ren and S.B. Pope, Cornell University, 2002.
c  Most lower level routines originally written by S. Subramaniam.
c  Ref: S.Subramaniam and S.B. Pope, Combust. Flame, 115: 487-514 (1998).

	contains
c
c////////////////////////////////////////////////////////////
c
      subroutine bmod(npd,n,nkmax,msubn,msub,imod,sumwtk,b)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    By Shankar
c    Subroutine bmod returns the model EMST B coefficients
c
c     given D=number of dimensions ;
c     N = number of equal weight nodes, wmin=1/N
c     2 parameters per dimension which are
c     f(inf)/fmax, y*, ( and fmax).
c     Normalise edwt by the ensemble weight sumwtk.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer,intent(in):: npd,n,nkmax,msubn,imod, msub(msubn,2)
      real(kind(1.e0)),intent(in)::sumwtk(msubn)
      real(kind(1.e0)),intent(out):: b(nkmax-1,msubn)

      integer:: iens,nk,ie
      real(kind(1.e0)):: wsub
c

c_________________________________________________________
c     imod determines the model type.
c     do for all edges, find edwt relative to sub-ens wt.
c     imod = 1 is B(w) = const.
c_________________________________________________________
      if(imod.eq.1)then
        do iens=1,msubn
          nk =  msub(iens,2)
          do ie=1,nk-1
            b(ie,iens) = 1.
          enddo
        enddo
c_________________________________________________________
c     imod = 2 is B(w) = 2*w
c_________________________________________________________
      elseif(imod.eq.2)then
        do iens = 1, msubn
          nk =  msub(iens,2)
          do ie=1,nk-1
            wsub = b(ie,iens)/sumwtk(iens)
            b(ie,iens) = 2.*wsub
          enddo
        enddo
      endif
      return
      end subroutine bmod
c
c////////////////////////////////////////////////////////////
c
             subroutine control_set(cvars,icheck,Max_Ave,
     1                      sdev_small,omdtmin,fphthr,rphthr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c return the control parameters:
c by Z.Ren, S.B.Pope 2002
c
c Input:
c      cvars
c Output:
c     icheck   -The control parameter which decide whether the inputs of particle properties(wt,fscale,and state)
c               should be checked or not.
c     Max_Ave  -The control parameter which decide the criteria for performing IEM based on the maximum range
c               or based on the average range.
c    sdev_small-The control parameter which decides whether EMST or IEM will
c               being executed.   Default sdev_small=1.e-4.
c     omdtmin  -The control parameter which decides whether the rejection method
c               will be used or not. Default omdtmin=1.e-6.
c    fphthr
c    rphthr :  -fphthr and rphthr are two parameters used to decide the dphthr.
c               default value for fphthr is 0.3. default value for rphthr is 1.e-5.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             implicit none
	       real(kind(1.e0)), intent(in) :: cvars(6)
	       integer, intent(out) :: icheck,Max_Ave
	       real(kind(1.e0)), intent(out) :: sdev_small, omdtmin,
	1			 fphthr,rphthr

              real(kind(1.e0)) sdev_def,omdt_def, fph_def,rph_def

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              sdev_def=1.e-4
              omdt_def=1.e-6
              fph_def =0.3
              rph_def =1.e-5

            if( nint(cvars(1)) /= 1)  then
               icheck=0
             else
               icheck=1
            endif

            if(nint(cvars(2)) /= 1)  then
               Max_Ave=0
             else
               Max_ave=1
             endif

            if(cvars(3) <= 0.) then
               sdev_small=sdev_def
            else
               sdev_small=cvars(3)
            endif

            if(cvars(4) <= 1.e-7 .or. cvars(4) > 1.e-2) then
              omdtmin  =omdt_def
            else
              omdtmin  =cvars(4)
            endif

            if(cvars(5) <= 0.)then
              fphthr=fph_def
            else
              fphthr=cvars(5)
            endif

            if(cvars(6) <= 0.)then
              rphthr=rph_def
            else
              rphthr=cvars(6)
            endif

            end subroutine control_set
c
c////////////////////////////////////////////////////////////
c
      subroutine covcal(mode,nd,n,ndimd, ndim, w, y, ymean, cov )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Find the means/covariance matrix of y
c    mode=1, return the means of y.
c        =2, return both the means and covariance of y.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none
       integer,intent(in):: mode, nd, n, ndimd,ndim
       real(kind(1.e0)),intent(in):: w(nd),y(nd,ndimd)
       real(kind(1.e0)),intent(out):: cov(ndimd,ndimd),ymean(ndimd)
       integer:: j,idim, jdim
       real(kind(1.e0))::  ave, s
       real(kind(1.d0)):: dsum
c
      dsum = 0.d0
c____________________________________________________________
c
c     First find the means of the y's
c____________________________________________________________
      do idim=1,ndim
        dsum=0.d0
        do j=1,n
          dsum = dsum+dble(w(j)*y(j,idim))
        enddo
        ave = dsum
c____________________________________________________________
c
c     w(i) are normalized such that they sum to unity.
c____________________________________________________________
        ymean(idim)=ave
      enddo
      if (mode .eq. 1) then
        return
      elseif (mode .eq. 2) then
        do idim=1,ndim
          do jdim = idim,ndim
            dsum = 0.d0
            do j=1,n
              dsum=dsum+dble( w(j)*(y(j,idim)-ymean(idim))*
     &           (y(j,jdim)-ymean(jdim)))
            enddo
            s=dsum
            cov(idim,jdim) = s
          enddo
        enddo
c____________________________________________________________
c     Fill in the lower half of the covariance matrix by symmetry.
c____________________________________________________________
        do idim=1,ndim
          do jdim = 1,idim-1
            cov(idim,jdim) = cov(jdim,idim)
          enddo
        enddo
      endif
c     
      return
      end subroutine covcal
c
c////////////////////////////////////////////////////////////
c
	subroutine edgewt( npd, np, ndim, nkmax, msubn, msub, mst, lc,
     $   rs, parent, isp, w, edwt, sumwt, wtnode, stack )
	integer npd, np, ndim, istart, nk, mst(nkmax-1,2,msubn),
     $     lc(npd), rs(npd), parent(npd,2), node, top, child,
     $     stack(npd), isp(npd), msub(msubn,2)
	real edwt(nkmax-1,msubn), w(npd), wtnode(npd), sumwt(msubn)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       By Shankar
c       Compute edge weights from the LC-RS representation of
c       the EMST in an efficient manner.
c	Latest version : assumes LC and RS store relative
c	LC and RS addresses in absolute locations.
c	Also stores all the edge weights of all edges of
c	all sub-ensembles simultaneously in the b arrays,
c	which are later overwritten by the B coefficients
c	themselves.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	Initialize edwt to zero.

	do iens = 1, msubn
          istart = msub(iens,1)
          nk= msub(iens,2)
          sumwt(iens) = 0.
          do i = 1, nk
            stack(i) = 0
            wtnode(i) = 0.
            sumwt(iens) = sumwt(iens) + w(isp(istart + i - 1))
          enddo

c	Start at the root node.

          node = mst(1,1,iens)
          top = 1
          stack(top) = node

c	Descend down the tree.
c	If there are children : pick the leftmost first.

  2       if( lc(node+istart-1).ne.0 )then

c	Check if this path has been traversed.

            if( wtnode(lc(node+istart-1)).eq.0. )then
              node = lc(node+istart-1)
              top = top + 1
              stack(top) = node
              goto 2
            else

c	If traversed, look for other children not accounted for

              ihroot = lc(node+istart-1)
  4           if( rs(ihroot+istart-1).ne.0 )then

c	Check if this path has been traversed.

                if( wtnode(rs(ihroot+istart-1)).eq.0. )then
                  node = rs(ihroot+istart-1)
                  top = top + 1
                  stack(top) = node
                  goto 2
                else
                  ihroot = rs(ihroot+istart-1)
                  goto 4
                endif
              else

c	Finished all paths at this level : pop current node off the
c       stack.
c	If nk > 2 need check for quitting only here.
	
                wtnode(node) = wtnode(node) +
     $             w( isp(istart + node - 1))
                child = node
                top = top - 1
                if(top.gt.0)then
                  node = stack(top)
                  wtnode(node) = wtnode(node) + wtnode(child)
                  goto 2
                else
                  goto 3
                endif
              endif
            endif
          else

c	Reached a leaf : add its weight and pop it off the stack.
	
            wtnode(node) = wtnode(node) + w(isp(istart + node - 1))
            child = node
            top = top - 1
            node = stack(top)
            wtnode(node) = wtnode(node) + wtnode(child)
            goto 2
          endif
  3       do ie = 1, nk-1
            i = mst(ie,1,iens)
            j = mst(ie,2,iens)
            if(j.eq.parent(i+istart-1,1)) i = j
            edwt(ie,iens) = wtnode(j)
            edwt(ie,iens) = min(edwt(ie,iens),sumwt(iens) -
     $         edwt(ie,iens))
          enddo
        enddo
	return
	end subroutine edgewt
c
c////////////////////////////////////////////////////////////
c
      subroutine emst_init(n,z0l,z0u,z1l,z1u,state)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inistat initializes the state of the state() array for mixemst
c input:
c     n - number of particles
c     z0l, z0u - lower and upper bounds of the Z0 r.v.
c     z1l, z1u - lower and upper bounds of the Z1 r.v.

c output:
c     state() - age property array of particles

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none
	integer, intent(in) :: n
	real(kind(1.e0)), intent(in) :: z0l,z0u,z1l,z1u
	real(kind(1.e0)), intent(inout) :: state(n)
        integer:: i, nt, int, notint
        real(kind(1.e0)):: p0,work(n)

      p0 = (z1l + z1u)/((z0l+z0u) + (z1l+z1u))

      call emst_rnu(state,n)

      nt = 0
      do i = 1, n
        if ( state(i) <= p0 ) then
          state(i) = 1.
          nt = nt + 1
        else
          state(i) = -1.
        endif
      enddo

      if( nt > 0 )call raninit(nt, z1l, z1u, work)
      if( nt < n )call raninit((n-nt), z0l, z0u, work(nt+1))

      int = 1
      notint = nt+1
      do i = 1, n
        if ( state(i) > 0. ) then
          state(i) = work(int)
          int = int + 1
        else
          state(i) = -work(notint)
          notint = notint + 1
        endif
      enddo

      return
      end subroutine emst_init
c
c////////////////////////////////////////////////////////////
c
      subroutine emst_scl(np,nc,g,state,wt,omdt,cvars,cons)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Perform mixing on the scale compositions g.

c  by Z.Ren and S.B. Pope, Cornell University, 2002

c Input:
c   np    -The number of particels in the ensemble. >=1 (integer)
c   nc    -The number of composition variables.     >=1 (integer)
c   g     -scaled particle composition
c   state -age properties for particles. The state variable state(i) is the age property
c          of the i-th particle.
c   wt    -particle weights. wt(i)>0 is the numerical weight of the i-th particle.
c   omdt  -normalized mixing time
c   cvars -control variables
c   cons  -constants in state equation
c Output:
c   g     -scaled particle composition after mixing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        integer, intent(in) ::  np, nc
        real(kind(1.e0)), intent(in) :: wt(np),  omdt, cvars(6),cons(4)
        real(kind(1.e0)), intent(inout) :: g(np,nc), state(np)

c       local variables

        integer:: i,j,imod,msubn,icheck,Max_Ave

        real(kind(1.e0))::  gtemp(np,nc),dphthr(nc)
        real(kind(1.e0))::  omtleft, MaxAvevar,range ,range_final,omdtm
        real(kind(1.e0))::  a(1),gtmax,sdev_small,omdtmin,fphthr,rphthr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c-----get the control parameters

         call control_set(cvars, icheck,Max_Ave,
     1           sdev_small,omdtmin,fphthr,rphthr)

         msubn   =1
         imod    =2
         gtmax   =2.5

         call Max_Avevar(Max_Ave,np,nc,g,wt,MaxAvevar,range)

c  -----------------for large omdt using IEM and return.
         range_final= range*exp(-omdt)

         if (range_final < sdev_small/10.) then
           call iem(np, nc, g, wt, omdt)
           call statinc(np,omdt,cons(1),cons(2),cons(3),cons(4),
     1                    state)
           return        ! all done
         endif

         omtleft =omdt
         do while ( omtleft > 0. ) !===keep mixing until omtleft=0 =====

            call Max_Avevar(Max_Ave,np,nc,g,wt,MaxAvevar,range)

c--------------perform IEM for small variance ---------------

            if(MaxAvevar < 1.e-12 .or. range < 4*sdev_small) then
                call iem(np, nc, g, wt, omtleft)
                call statinc(np,omtleft,cons(1),cons(2),cons(3),cons(4),
     1                     state)
	          return     ! all done
	      endif

c--------------perform EMST ------------------------

c   mix up to omdtm
c   deal with small omdt. if omdt<omdtmin: rejection method is used for EMST mixng model.
c
            if (omtleft <= omdtmin) then

                     gtemp(1:np,1:nc)=g(1:np,1:nc)

                     call mixemst (np,nc,state,gtemp,wt,2*omdtmin,cons,
     1                omdtm,gtemp,dphthr,fphthr,rphthr,msubn,imod,gtmax)

                     if(omdtm >= omtleft) then ! too much mixing
                         ! use rejection method to accept or reject
                         call emst_rnu(a,1)
                         if (a(1) <= omtleft/omdtm)
     1                      g(1:np,1:nc)=gtemp(1:np,1:nc)  ! accept

                         call statinc(np,omtleft,cons(1),cons(2),
     1                          cons(3),cons(4), state)
                         return	! all done
                     else
c					  accept mixing step; further mixing required
                        g(1:np,1:nc)=gtemp(1:np,1:nc)
                        call statinc(np,omdtm,cons(1),cons(2),
     1                          cons(3),cons(4), state)
                       omtleft=omtleft-omdtm
                     endif

            else	! regular case
                     call mixemst(np,nc,state,g,wt,omtleft,cons,omdtm,
     1                        g,dphthr,fphthr,rphthr,msubn,imod,gtmax)
c omdtm is the actual mixing time performed when we performed EMST Mixing on the Ensemble
c with desired mixing time omtleft.

                    call statinc(np,omdtm,cons(1),cons(2),
     1               cons(3),cons(4), state)

                    omtleft = omtleft - omdtm
            end if
         end do  !===============================

         return
	   end subroutine emst_scl

c
c////////////////////////////////////////////////////////////
c
	   subroutine emsti( ned, np, ndim, efdim, nefdim, nkmax, msubn,
     $   isp, msub, x, dismin, mst, itree, inear )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     by Shankar
c
c     Jan 4 1994 : edited to include difference betn ned and np.
c     Oct 18 : including effective dimension ideas.
c     Sep 3 1993 : changed to account for the difference in
c     the mst array - now stores all the trees.
c     and emsti computes all the trees at one shot.
c     Euclidean minimum-cost spanning tree
c     Given a set of vertices V finds the emst using Prim's
c     algorithm - evaluates cost function as it goes along
c     i.e. n(n-1)/2 evaluations of distances.
c     x and y store the locations of the n vertices.
c     array mst stores the n-1 edges of the minimum spanning tree
c     which is basically a list of n vertices.
c
c     v 1.2 : implemented the foll improvements
c
c     1. Passing parameters instead of using common blocks.
c     2. Indirect addressing.
c     3. Modified to form an EMST of a sub-ensemble of size
c     nk from an ensemble of size n.
c     4. Array  x contains the full ensemble : n particles, ndim.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none
       integer,intent(in)::ned,np,ndim,nefdim,nkmax,msubn
       integer,intent(in):: msub(msubn,2),isp(ned),efdim(ndim)
       integer,intent(out)::mst(nkmax-1,2,msubn),inear(ned),itree(ned)
       integer::iens,istart,i,iptr,iefdim,intree,nk,k,kloc,jnode,itmp,j
       integer::jptr,kptr,idim
       real(kind(1.e0))::  x(ned,ndim)
       real(kind(1.d0)):: dismin(ned),dmin,dist

c     Initialize the set of vertices U that form the tree at any
c     time to be the first vertex. also initialize the lowcost(dismin)
c     and inear arrays.

      do iens=1,msubn
        istart = msub(iens,1)
        nk= msub(iens,2)
        itree(1) = 1
        do i = 2, nk
          iptr = isp(i + (istart-1))
          dismin(i) = 0.d0
          do iefdim = 1, nefdim
            idim = efdim(iefdim)
            dismin(i) =  (dble(x(iptr,idim)-x(isp(istart),idim)))**2
     $         + dismin(i)
          enddo
          inear(i) = 1
          itree(i) = i
        enddo
c____________________________________________________________
c
c     Now pick the minimum lowcost vertex in V-U and connect it to its
c     inear node in the tree U.
c     k stores the current vertex under consideration.
c     Using indirect addressing the array itree is a pointer
c     from the from the set of vertices in the tree to the
c     initial ordering of the vertices, so that after k iterations
c     itree stores the k nodes in the tree and itree(k+1) to
c     itree(n) are not in the tree.
c
c____________________________________________________________
        intree=1
        do i = 2, nk
          dmin = dismin(itree(intree+1))
          k=itree(intree+1)
          kloc = intree+1
          do j = intree+2, nk
            jnode = itree(j)
            if(dismin(jnode).lt.dmin)then
              dmin = dismin(jnode)
              k=jnode
              kloc=j
            endif
          enddo
          mst(intree,1,iens) = inear(k)
          mst(intree,2,iens) = k
          intree = intree + 1
          itmp = itree(intree)
          itree(intree) = k
          itree(kloc) = itmp
          do j = intree+1, nk
c
c     For all vertices not in the tree update inear and dismin
c
            jnode = itree(j)
            jptr = isp(jnode + (istart-1))
            kptr = isp(k + (istart-1))
            dist = 0.d0
            do iefdim=1,nefdim
              idim = efdim(iefdim)
              dist = dist + (dble(x(jptr,idim)-x(kptr,idim)))**2
            enddo
            if(dist.lt.dismin(jnode))then
              dismin(jnode)=dist
              inear(jnode) = k
            endif
          enddo
        enddo
      enddo
      return
      end subroutine emsti
c
c////////////////////////////////////////////////////////////
c
 	subroutine iem(  np, nc, g, wt,omdt )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Perform IEM mixing on the particles in the ensemble.
c   by Z.Ren, S.B.Pope 2002
c
c  Input:
c	np  - The number of particles in the ensemble.
c     nc  - The number of composition variables
c	g	- particle compostions : g(i,j) is the j-th composition of the i-th
c		- particle in the ensemble
c     wt  - particle weights.
c     omdt-normalized mixing time.
c
c  Output:
c       g  -particle compositions after mixing.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
	integer, intent(in) :: np,nc
	real(kind(1.e0)), intent(in) :: wt(np), omdt
	real(kind(1.e0)), intent(inout) :: g(np,nc)

        integer :: i,j,k
	real(kind(1.d0)) :: arg, gmean,wtsum,decay
c
c  form sum of weights: return if sumwt <= 0 or np=0
c
	if( np == 0 ) return

      wtsum=0.0d0
      do i=1,np
	  wtsum=wtsum+wt(i)
      end do
c
	if( wtsum <= 0.d0 )  return

c
c  decay rate
c
	arg       = 0.5d0 * omdt
	decay     = 1.d0 - exp( -arg )
c
c  loop over compositions
c
	do k = 1,nc
        gmean=0.0
        do j=1,np
	     gmean    = gmean+g(j,k)*wt(j)
        end do
        gmean =gmean/wtsum
c
	  do i = 1, np
	    g(i,k)   = g(i,k) - decay * ( g(i,k) - gmean )
        end do
      end do
c
	return
	end subroutine iem
c
c////////////////////////////////////////////////////////////
c
	subroutine implct( npd, n, ndim, nkmax, msubn,
     $   msub, phit, phith, w, isp, mst, lc, rs,
     $   parent, b, eta, coefi, coefd, coefe,
     $   dphi, stack, itrav)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     By Shankar
c     Integration of the ODEs by an implicit, first-order method.
c     This method is unconditionally stable in the limit w -> 0.
c     Matrix solution by an direct method.
c     phi(t + h) = phi(t) - (h)*dphidt(t+h,phi(t+h)) + order(h*2)
c
c prodedure:
c     For all dimensions,
c     Start from root and descend to a leaf and work back up.
c     Back substitution :
c     Get root node soln.
c     For each child calculate the soln. recursively.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer mst(nkmax-1,2,msubn), parent(npd,2), lc(npd), rs(npd),
     $   child, stack(npd), top, itrav(npd), isp(npd), childa,
     $   msub(msubn,2)
      real phit(npd,ndim), phith(npd,ndim), w(npd), b(nkmax-1,msubn),
     $   coefd(npd), coefe(npd), coefi(npd), dphi(npd)
      
c____________________________________________________________
      do iens=1,msubn
        istart = msub(iens,1)
        nk = msub(iens,2)
        do i = 1,nk
          itrav(i) = 0
          coefi(i) = 0.
          coefd(i) = 0.
          coefe(i) = 0.
        enddo
        do idim =1,ndim

c     First get the Iks.

          iroot = mst(1,1,iens)
          do node=1,nk
            nodea = node+istart-1
            if(node.ne.iroot)then
              coefi(node) = coefi(node) + b(parent(nodea,2),iens)
              coefd(node) = coefd(node) - b(parent(nodea,2),iens)
     $           *( phit(isp(nodea),idim) - 
     $           phit(isp(istart+parent(nodea,1)-1),idim))
            endif
            child = lc(nodea)
            if(child.ne.0)then
              childa = child + istart - 1
              coefi(node) = coefi(node) + b(parent(childa,2),iens)
              coefd(node) = coefd(node) - b(parent(childa,2),iens)
     $           *( phit(isp(nodea),idim) - phit(isp(childa),idim))
              ihroot = child
              ihrota = childa
  40          if(rs(ihrota).ne.0)then
                child = rs(ihrota)
                childa = child + istart - 1
                coefi(node) = coefi(node) + b(parent(childa,2),iens)
                coefd(node) = coefd(node) - b(parent(childa,2),iens)
     $             *( phit(isp(nodea),idim) - phit(isp(childa),idim))
                ihroot = child
                ihrota = childa
                goto 40
              endif
            endif
c     
            coefi(node) = 1./(coefi(node)*eta + w(isp(nodea)))
          enddo
c     
c     Start at the root node.
c
          node = mst(1,1,iens)
          nodea = node + istart - 1
          top = 1
          stack(top) = node
c     
c     Descend down the tree.
c     If there are children : pick the leftmost first.
c     
  2       if( lc(nodea).ne.0 )then
c     
c     Check if this path has been traversed.
c     
            if( itrav(lc(nodea)).eq.0 )then
              node = lc(nodea)
              nodea = node + istart - 1
              itrav(node) = 1
              top = top + 1
              stack(top) = node
              goto 2
            else
c     
c     If traversed, look for other children not accounted for
c     
              ihroot = lc(nodea)
              ihrota = ihroot + istart - 1
  4           if( rs(ihrota).ne.0 )then
c     
c     Check if this path has been traversed.
c     
                if( itrav(rs(ihrota)).eq.0 )then
                  node = rs(ihrota)
                  nodea = node + istart - 1
                  itrav(node) = 1
                  top = top + 1
                  stack(top) = node
                  goto 2
                else
                  ihroot = rs(ihrota)
                  ihrota = ihroot + istart - 1
                  goto 4
                endif
              else
c     
c     Finished all paths at this level : pop current node off the stack.
c     All children' contributions are now added to the dead-end node.
c     If nk.gt.2 , then lc(iroot).ne.0 , so check for quitting only
c     here.
c     
                coefd(node) = coefi(node)*eta*coefd(node)/
     $             ( 1. - coefi(node)*eta*coefe(node) )
                if(node.ne.iroot)then
                  coefe(node) = coefi(node)*eta*b(parent(nodea,2),iens)
     $               / ( 1. - coefi(node)*eta*coefe(node) )
                endif
                child = node
                childa = child + istart - 1
                top = top - 1
                if(top.gt.0)then
                  node = stack(top)
                  nodea = node + istart - 1
                  coefd(node) = coefd(node) + 
     $               coefd(child)*b(parent(childa,2),iens)
                  coefe(node) = coefe(node) + 
     $               coefe(child)*b(parent(childa,2),iens)
                  goto 2
                else
                  goto 3
                endif
              endif
            endif
          else
c     
c     Reached a leaf : find its coefs and pop it off the stack.
c     
            coefd(node) = coefi(node)*coefd(node)*eta
            coefe(node) = coefi(node)*eta*b(parent(nodea,2),iens)
            child = node
            childa= child + istart - 1
            top = top - 1
            node = stack(top)
            nodea = node + istart -1
c     
c     To the leaf's parent add the child's contributions.
c     
            coefd(node) = coefd(node) + 
     $         coefd(child)*b(parent(childa,2),iens)
            coefe(node) = coefe(node) + 
     $         coefe(child)*b(parent(childa,2),iens)
            goto 2
          endif
  3       continue
c     
c     Start at the root node.
c     
          do i=1,nk
            itrav(i) = 0
          enddo
          node = mst(1,1,iens)
          nodea = node + istart - 1
          top = 1
          stack(top) = node
          dphi(node) = coefd(node)
          coefi(node) = 0.
          coefd(node) = 0.
c     
c     Descend down the tree.
c     If there are children : pick the leftmost first.
c     
  5       if( lc(nodea).ne.0 )then
c     
c     Check if this path has been traversed.
c     
            if( itrav(lc(nodea)).eq.0 )then
              iparent = node
              node = lc(nodea)
              nodea = node + istart - 1
              itrav(node) = 1
              dphi(node) = coefd(node) + coefe(node)*dphi(iparent)
              coefi(node) = 0.
              coefd(node) = 0.
              top = top + 1
              stack(top) = node
              goto 5
            else
c     
c     If traversed, look for other children not accounted for
c     
              ihroot = lc(nodea)
              ihrota = ihroot + istart - 1
  6           if( rs(ihrota).ne.0 )then
c     
c     Check if this path has been traversed.
c     
                if( itrav(rs(ihrota)).eq.0 )then
                  iparent = node
                  node = rs(ihrota)
                  nodea= node + istart - 1
                  itrav(node) = 1
                  dphi(node) = coefd(node) + coefe(node)*dphi(iparent)
                  coefi(node) = 0.
                  coefd(node) = 0.
                  top = top + 1
                  stack(top) = node
                  goto 5
                else
                  ihroot = rs(ihrota)
                  ihrota = ihroot + istart - 1
                  goto 6
                endif
              else
c     
c     Finished all paths at this level : pop current node off the stack.
c     If nk.gt.2 , then lc(iroot).ne.0 , so check for quitting only
c     here.
c     
                child = node
                childa = child + istart - 1
                top = top - 1
                if(top.gt.0)then
                  node = stack(top)
                  nodea = node + istart - 1
                  goto 5
                else
                  goto 7
                endif
              endif
            endif
          else
c     
c     Reached a leaf 
c     
            child = node
            childa = child + istart - 1
            top = top - 1
            node = stack(top)
            nodea = node + istart - 1
            goto 5
          endif
  7       continue
c     
          do i=1,nk
            itrav(i) = 0
            iabs= istart + i -1
            phith(isp(iabs),idim) = phit(isp(iabs),idim)+dphi(i)
            dphi(i) = 0.
            coefe(i) = 0.
            coefi(i) = 0.
            coefd(i) = 0.
            stack(i) = 0.
          enddo
         enddo
        enddo
      return
      end subroutine implct
c
c////////////////////////////////////////////////////////////
c
        subroutine Max_Avevar(Max_Ave,np,nc,g,wt,MaxAvevar,range)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calulate the averge variance or maximum variance of the ensemble
c  and the average difference or maximum difference of the ensemble in the compostion space 
c  with nc compositions according to the value of Max_Ave.

c  by Z.Ren, S.B.Pope 2002

c  Input:
c   Max_Ave =0 -calculate the maximum composition variance of g. Maximum variance=max(var(g(:,j))).
c              -calculate the maximum difference of g in composition space range.
c              -range=max_j(  max_i( f(i,j) ) - min_i( f(i,j) ) )
c           =1 -calculate the average composition variance of g. Average variance=sum(var(g(:,j)))/nc.
c              -calculate the average difference of g in composition space range.
c              -range=sum(  max_i( f(i,j) ) - min_i( f(i,j) ))/nc.
c        np -The number of particles in the ensemble.
c        nc -the number of composition variables.
c         g -particle compostions : g(i,j) is the j-th composition of the i-th
c	    -particle in the ensemble
c        wt -particle weights.
c  output:
c  MaxAvevar-maximum composition variance of the ensemble when Max_Ave=1.
c  MaxAvevar-average composition variance of the ensemble when Max_Ave=2.
c     range-maximum difference in composition space when Max_Ave=1.
c     range-average difference in composition space  when Max_Ave=2.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
	  integer, intent(in) :: Max_Ave,np,nc
	  real(kind(1.e0)), intent(in) ::  wt(np), g(np,nc)
	  real(kind(1.e0)), intent(out)::  MaxAvevar,range

        integer :: i,j
        real(kind(1.e0)) :: gmin(nc),gmax(nc)
        real(kind(1.d0)) ::  ave(nc),var(nc),sumwt,sumvar
c
c  Note: sums are accumulated in double precision.  Hence intrinsics (sum, matmul) are not used.

        do i=1,nc
            gmin(i)=minval(g(1:np,i))
            gmax(i)=maxval(g(1:np,i))
        end do

        sumwt=0.d0
        do i=1,np
          sumwt=sumwt+wt(i)
        end do

	  ave(1:nc)=0.d0
        var(1:nc)=0.d0

        do i=1,nc
          do j=1,np
            ave(i)=ave(i)+g(j,i)*wt(j)
          end do
          ave(i)=ave(i)/sumwt	! mean

          do j=1,np
            var(i)=var(i)+ wt(j)*(g(j,i)-ave(i))**2
          end do
          var(i)=var(i)/sumwt ! variance
	  end do

        if(Max_Ave == 0) then  ! take maxima
            MaxAvevar=maxval(var(1:nc))
            range=maxval(gmax(1:nc)-gmin(1:nc))
        else
	      MaxAvevar=sum(var(1:nc))/float(nc)
            range=sum( gmax(1:nc)-gmin(1:nc) ) /float(nc)
        endif

        end subroutine Max_Avevar 
c
c////////////////////////////////////////////////////////////
c
	  subroutine mixemst(np, ncomp,
     $   state, phi0, wt, dt,  cons,
     $    dtm, phim,dphthr,fphthr,rphthr,msubn,imod,gtmax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     mixemst - the EMST mixing module for intermittent mixing.
c     Z.Ren, S.B. Pope 2002
c     Shankar Subramaniam, January 1997
c
c algorithm :
c     for a flow chart describing the algorithm see notes on
c     Description of the Intermittent EMST Mixing Model.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c inputs:
c     np : The numble of particles in the ensemble.
c     ncomp : number of composition variables.
c     state:age properties for particles.
c     phi0 : array of composition values before mixing.
c     wt : array of particle weights, dimensioned np.
c     dt : time step over which mixing is desired
c     cons : model constants.
c Outputs:
c     dtm : time step over which mixing was actually performed subject
c           to the constraint on dphi.
c     phim : array of compositions after mixing over dtm.
c
cLocal parameter:
c     msubn : number of sub-ensembles (must be 1)
c     dphthr : threshold value for change in phi due to EMST mixing:
c     rphthr : parameter determining dphthr
c     fphthr : parameter determining dphthr
c     z0l, z0u, z1l, z1u : sample space intervals for Z_0 and Z_1
c     efdim : array of effective dimensions : ncomp
c     mst : array of edge indices for each sub-ensemble : np-1, 2,
c     itree : storage array for emsti
c     inear : -"-
c     msub : array storing starting particle index (istart) and
c     number of mixing particles (nk) for each sub-ensemble : msubn, 2
c     lc : leftmost child array : np
c     rs : right sibling array : np
c     parent : stores parent node and connecting edge : np, 2
c     stack : storage array for emsti
c     edptr : - " -
c     itrav : - " -
c     isp : cross reference array ( for randomly forming subensembles)
c     retained for historical reasons: np
c     (isp(i) = i for msubn = 1)
c     phi : array of compositions at initial time : np, ncomp
c     phihat : array of stdized compositions : np, ncomp
c     phitr : array of compositions to form the tree: np, ncomp
c     phivfn : array of compositions to calculate variance function:
c     np, ncomp
c     phitmp : - " -
c     phit : array of compositions after mixing : np, ncomp
c     b : array of B coefficients for each edge in each sub-ensemble
c     np - 1, msubn
c     sumwtk : sum of particle weights for each sub-ensemble : msubn
c     temp1 : temp storage array : np
c     wnm : particle weights normalized : np
c     coefi : I coefficient in matrix soln : np
c     coefd : D coefficient in matrix soln : np
c     coefe : E coefficient in matrix soln : np
c     q : orthogonal matrix of singular vectors: ncomp, ncomp
c     sigma : vector of singular values of cov matrix: ncomp
c     work : work array : ncomp
c     dphi : composition changes : np
c     expphi : mean scalar values : ncomp
c     covphi : scalar covariance matrix : ncomp, ncomp
c     avg : mean scalar values : ncomp
c     wtmin : minimum partice weight
c     state : particle age property : np
c     ds : incremental scaled time
c     gtmax : upper limit on normalized tree variance decay rate (relative
c	       to variance decay rate)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     type declaration of parameters passed

      integer   msubn, imod, np, ncomp
      real phi0(np, ncomp), wt(np), dt,  dtm, gtmax,
     $   phim(np, ncomp), state(np), cons(4),fphthr, rphthr



c     local type declarations for parameters

      integer dxiefl, jacfl

      real dphthr(ncomp),  z0l, z0u, z1l, z1u, 
     $       dxie(ncomp), jacob(ncomp,ncomp),
     $   tinit, tfin, qmin, qmax, qevolgl(5), p, fracin

      integer efdim(ncomp), mst(np-1,2,msubn), itree(np),
     $   inear(np), msub(msubn,2), lc(np), rs(np),
     $   parent(np,2), stack(np),
     $   edptr(np), itrav(np), isp(np),
     $   i, j, k, idim, nk, iefdim,    freac,
     $    mode, iens, istart, ntsub, it, iabs, nq

      real phi(np,ncomp), b(np-1,msubn), sumwtk(msubn),
     $   temp1(np), wnm(np), r(ncomp), coefi(np), coefd(np),
     $   coefe(np), q(ncomp,ncomp), sigma(ncomp), work(ncomp),
     $   dphi(np), sd(ncomp), mxsdev,
     $   phihat(np,ncomp), phitr(np,ncomp), phivfn(np,ncomp),
     $   phitmp(np,ncomp),
     $   expphi(ncomp), covphi(ncomp,ncomp), avg(ncomp),
     $   wtmin, covtmp(ncomp,ncomp),
     $   totwt, dtlim, dtage, deltat, varfn, varfnT,
     $   alpest, dphimx(ncomp),  ds, dphi0(np, ncomp),
     $   phrxmn(ncomp), phrxmx(ncomp), dphrx, ave, adev, var, sdev,
     $   svar, skew, curt,  t0, tf

      double precision dismin(np)
c----------------------------------------------------------------------
	z0l  = cons(1)
	z0u  = cons(2)
	z1l  = cons(3)
	z1u  = cons(4)

c     Find dphthr
c
        totwt = 0.0
        do i = 1, np
          do idim = 1, ncomp
            phi(i, idim) = phi0(i,idim)
          enddo
          wnm(i) = wt(i)
          totwt = totwt + wnm(i)
        enddo
        do  i = 1, np
          wnm(i) = wnm(i) / totwt
        enddo

c
        mxsdev = 0.0
        do idim = 1, ncomp
          call moment( 2,  np, wnm, phi(1,idim),  ave,
     $       adev, sdev, var, skew, curt )
        sd(idim) = sdev
        mxsdev=max(mxsdev,sdev)
        enddo
c
       do i = 1,ncomp
          dphthr(i) = fphthr*max(sd(i), rphthr*mxsdev)
       enddo

c
c     assign particle compositions from phi0 to phi
c     assign particle weights to wnm and normalize such that sum=1.0

      totwt = 0.0
      do i = 1, np
        do idim = 1, ncomp
          phi(i, idim) = phi0(i,idim)
        enddo
        wnm(i) = wt(i)
        totwt = totwt + wnm(i)
      enddo
      do  i = 1, np
        wnm(i) = wnm(i) / totwt
      enddo
      wtmin=wnm(1)
      do i = 1, np
        if(wnm(i).lt.wtmin)then
          wtmin=wnm(i)
        endif
      enddo
c

c

c     form fluctuating component by subtracting the means
c     (stored in avg)
c
        call stripm( 1, np, ncomp, np, avg, wnm, phi)

        do idim = 1,ncomp
          efdim(idim) = idim
        enddo

c


c     define the tree composition array
c

        do idim = 1, ncomp
          do i = 1, np
            phitr(i,idim) = phi(i,idim)
          enddo
        enddo
c
c     define the composition array for variance function calculation
c
        do idim = 1, ncomp
          do i = 1, np
            phivfn(i,idim) = phi(i,idim)
          enddo
        enddo
c____________________________________________________________
c     Form  subensembles randomly. temp1 is a work array.
c____________________________________________________________
c
      call subensi(np,ncomp,msubn,wnm,wtmin,isp,msub,temp1)
c
      call subint(np,state,msubn,isp,msub)
c
      call pcheck(np, ncomp, np, ncomp, phivfn, phitmp, wnm, isp,
     $   state, msubn, msub, varfn, varfnT, temp1, expphi,
     $   covtmp, z1l, z1u, gtmax)
c

c     form tree and generate edges and edge coefficients
c     form the mixing matrix C_nm. first get the model (B)
c     coefficients.
c     using a linear profile (imod=2).
c
      call emsti( np, np, ncomp, efdim, ncomp, np, msubn,
     $   isp, msub, phitr, dismin, mst, itree, inear )
c
      call tree( np, np, np, msubn, msub, mst, lc, rs,
     $   parent, stack, edptr)
c
      call edgewt( np, np, ncomp, np, msubn, msub, mst,
     $   lc, rs, parent, isp, wnm, b, sumwtk, temp1, stack )
c
      call bmod( np, np, np, msubn, msub, imod, sumwtk, b )
c
c     check the time step
c
       dtage =z1l
       deltat = amin1( dt,  dtage)
c     estimate the root
c
      call rtest(np, ncomp, np, phivfn, isp, msubn, msub, mst,
     $   b, varfn,  alpest)
c
      call covcal(2, np,np,ncomp,ncomp,wnm,phi,expphi,covphi)

c    find the root using rejection method
      call rtfini(np,  ncomp, phi, phitmp, phivfn, msubn,
     $   msub, isp, wnm, varfn, dphi0,  deltat,  alpest,
     $    efdim, dphi, mst, lc, rs, parent, b, coefi, coefd,
     $   coefe, stack, itrav,  expphi, covtmp, phihat, q,
     $   sigma, dphthr, dphimx,  work, dtm)


      call covcal(2, np,np,ncomp,ncomp,wnm,phi,expphi,covtmp)
c
      call stripm( 2, np, ncomp, np, avg, wnm, phi)

      do idim = 1, ncomp
        do i = 1, np
          phim(i, idim) = phi(i,idim)
        enddo
      enddo
c
      return
      end subroutine mixemst
c
c////////////////////////////////////////////////////////////
c
      subroutine moment( mode,  n, w, data,  ave, adev, sdev,
     $   var, skew, curt )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Subroutine to calculate the moments of a sample : cf Numerical Recipes.
c     modified for UNEQUAL weights.
c     modified so it calculates only the moments needed.
c     mode = 1 mean only
c     mode = 2 mean and variance
c     mode = 3 mean, variance and skewness
c     mode = 4 mean, variance, skewness and kurtosis.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer n, mode ,j
      real  ave, adev, sdev, var, skew, curt
      real data(n),w(n)
      double precision dsum,dwt,dfluc,p,dvar,dskew,dcurt,dflucw

      dsum=0.d0
      adev=0.
      var=0.
      skew=0.
      curt=0.
c
      dwt = 0.d0
      dfluc = 0.d0
      dflucw = 0.d0
      dsum= 0.d0
      dvar=0.d0
      dskew=0.d0
      dcurt=0.d0
c
      do j = 1,n
        dsum = dsum + dble(w(j)*data(j))
      enddo
      ave=sngl(dsum)
c
      if ( mode .eq. 1 ) then
        return
      elseif( mode .eq. 2 )then
        do j = 1,n
          dwt = dble(w(j))
          dfluc=dble(data(j)-ave)
          dflucw=dfluc*dwt
          dsum = dsum + abs(dflucw)
          dflucw = dflucw*dfluc
          dvar=dvar+ dflucw
        enddo
        adev=sngl(dsum)
        var = sngl(dvar)
        sdev=sqrt(var)
      elseif( mode .eq. 3 )then
        do j = 1,n
          dwt = dble(w(j))
          dfluc=dble(data(j)-ave)
          dflucw=dfluc*dwt
          dsum = dsum + abs(dflucw)
          dflucw = dflucw*dfluc
          dvar=dvar+ dflucw
          dflucw = dflucw*dfluc
          dskew=dskew+dflucw
        enddo
        adev=sngl(dsum)
        var = sngl(dvar)
        skew= sngl(dskew)
        sdev=sqrt(var)
        skew=skew/(sdev**3)
      elseif( mode .eq. 4 )then
        do j = 1,n
          dwt = dble(w(j))
          dfluc=dble(data(j)-ave)
          dflucw=dfluc*dwt
          dsum = dsum + abs(dflucw)
          dflucw = dflucw*dfluc
          dvar=dvar+ dflucw
          dflucw = dflucw*dfluc
          dskew=dskew+dflucw
          dflucw = dflucw*dfluc
          dcurt=dcurt+dflucw
        enddo
        adev=sngl(dsum)
        var = sngl(dvar)
        skew= sngl(dskew)
        curt = sngl(dcurt)
        sdev=sqrt(var)
          skew=skew/(sdev**3)
          curt=curt/(var**2) - 3.
      endif
      return
      end subroutine moment
c
c////////////////////////////////////////////////////////////
c
      subroutine pcheck(np, ncomp, npcd, ndimd, phi, phitmp, w, isp,
     $   state, msubn, msub, varfn, varfnT, work, expphi, covtmp, t1l,
     $   t1u, gtmax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     By Shankar, Z.Ren, S.B.Pope 2002
c
c     pcheck checks if the value of p is large enough to achieve
c     the required decay rate of the variance function without
c     exceeding the threshold decay rate for any of the scalar
c     variances.
c input:
c     np : number of particles
c     ncomp : number of compositions
c     npcd : dimensioned number of particles in the ensemble.
c     ndimd : dimensioned number of compositions
c     phi : array of compositions : npcd, ndimd
c     phitmp : array of compositions : npcd, ndimd
c     wnm : array of normalized particle weights : npcd
c     isp : cross reference array
c     state : particle age property array : np
c     msubn : number of sub-ensembles
c     msub : particle indices of sub-ensembles : msubn,2
c     varfn : variance function for the whole ensemble
c     varfnT : variance function for particles in the tree
c     expphi : means of phi : ndimd
c     covtmp : temp covariance matrix array : ndimd, ndimd
c     t1l, t1u : limits for the T1 r.v.
c     gtmax : limiting normalized tree variance decay rate
c output:
c     msub : with modified number of particles in tree
c     state: modified particle age property array:
c     workspace
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       implicit none

      integer np, ncomp, ndimd, npcd, isp(npcd), msubn, msub(msubn,2)

      real phi(npcd, ndimd), w(npcd), state(np), work(np),
     $   expphi(ndimd), covtmp(ndimd, ndimd), varfn, varfnT,
     $   phitmp(npcd, ndimd), t1l, t1u, gtmax

c     local declarations

      integer idim, iens, nt, i, nk, istart, notint, ntsub, it, iabs,
     $   ntadd, npsub, nwarn


      real pmin, p, pnew, vadd


c     given nt calculate var fn (t) and var fn_T (t)

      varfn = 0.
      call covcal(2, npcd, np, ndimd,ncomp, w, phi, expphi, covtmp)
      do idim = 1, ncomp
        varfn = varfn + covtmp(idim,idim)
      enddo
      p = 0.
      nt = 0
      do iens = 1, msubn
        istart = msub(iens,1)
        ntsub = msub(iens, 2)

        do i = istart, istart + ntsub - 1
          p = p + w(isp(i))
        enddo
        nt = nt + ntsub
      enddo
       it = 1
      do iens = 1, msubn
        istart = msub(iens,1)
        ntsub = msub(iens, 2)
        do i = 1, ntsub
          iabs = isp(istart + i - 1)
          work(it) = w(iabs)/p
          do idim = 1, ncomp
            phitmp(it,idim) = phi(iabs,idim)
          enddo
          it = it + 1
        enddo
      enddo

      call covcal(2, npcd, nt, ndimd, ncomp, work, phitmp, expphi,
     $   covtmp)
      
      varfnT = 0.
      do idim = 1, ncomp
        varfnT = varfnT + covtmp(idim,idim)
      enddo
      if ( varfnT .gt. 0. ) then
        pmin = varfn/(varfnT*gtmax)
      else
        pmin = 1.
      endif
      pmin = min(pmin, 1.)

          
      if ( p .lt. pmin) then

        do iens = 1, msubn
          istart = msub(iens,1)
          ntsub = msub(iens,2)
          if (iens .le. msubn-1) then
            nk = msub(iens+1,1) - msub(iens,1)
          else
            nk = np - msub(iens,1) + 1
          endif
          notint = nk - ntsub
          do i = 1, notint
            work(i) = abs(state(isp(istart + ntsub - 1 + i)))
          enddo
          if ( notint .gt. 1 ) then
            call sort2(notint,work,isp(istart+ntsub))
          endif
        enddo

        do while (p .lt. pmin)

c     add a particle

          if ( nt .le. (np-1) ) then

            istart = msub(1,1)

            call raninit(1, t1l, t1u, state(isp(istart+nt)))
            pnew = p + w(isp(istart+nt))
            varfnT = varfnT*p/pnew
            vadd = 0.
            do idim = 1, ncomp
              vadd = vadd + (phi(isp(istart+nt),idim)-
     $           expphi(idim))**2
              expphi(idim) = (expphi(idim)*p +
     $           w(isp(istart+nt))*
     $           phi(isp(istart+nt),idim))/pnew
            enddo
            varfnT = varfnT +
     $         vadd*p*w(isp(istart+nt))/(pnew*pnew)
            nt = nt + 1
            msub(1,2) = nt

            if ( varfnT .gt. 0. ) then
              pmin = varfn/(varfnT*gtmax)
            else
              pmin = 1.
            endif
            pmin = min(pmin, 1.)
            p = pnew
          else
            return
          endif
        enddo

      else
        return
      endif
      return
      end subroutine pcheck
c
c////////////////////////////////////////////////////////////
c
      subroutine raninit(n,a,b,x)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine to generate random samples from the following
c     distribution.

c     if x in [0,a) F(x) = 2*x/(a+b)

c     if x in [a,b] F(x) = 2*a/(a+b) + ( 2*b*(x-a) - (x^2-a^2))/(b^2-a^2)

c     solved by inverse function method

c     1. generate F uniform on [0,1]

c     2. for each value of F find x:

c     F in [0, 2*a/(a+b) ), x = F*(a+b)/2

c     F in [2*a/(a+b) , 1], x = b - sqrt( b^2 - (F*(b^2-a^2)+a^2))

c input:

c     n - number of samples
c     a,b - parameters of the distn

c output:
c     x - array of n samples
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
	integer, intent(in) :: n
	real(kind(1.e0)), intent(in) :: a, b
	real(kind(1.e0)), intent(out) :: x(n)

      integer i

      real f, pa, b2, a2

      pa = 2.*a/(a+b)
      b2 = b*b
      a2 = a*a

      call emst_rnu(x,n)

      do i = 1, n
        if ( x(i) .lt. pa) then
          x(i) = x(i)*(a+b)*0.5
        else
          x(i) = b - sqrt( (b2 - a2)*(1.-x(i)) )
        endif
      enddo

      return
      end subroutine raninit
c
c////////////////////////////////////////////////////////////
c
	subroutine emst_rnu( x, n )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  routine to generate pseudo-random numbers uniformly in the
c  exclusive interval [ 0 , 1 ].
c  By S.B. Pope ,november 1990

c   x(i) , i = 1, n  array in which random numbers are returned
c
c  method: see f. james, computer physics communications, vol.60
c	p.329, 1990.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	data i1, i2 / 12345, 67890 /
c
	dimension x(n)
c
	do 100 i = 1, n
c
	ik = i1 / 53668
	i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
	if( i1 .lt. 0 ) i1 = i1 + 2147483563
c
	ik = i2 / 52774
	i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
	if( i2 .lt. 0 ) i2 = i2 + 2147483399
c
	ix = i1 - i2
	if( ix .lt. 1 ) ix = ix + 2147483562
c
	x(i) = float( ix ) * 4.656612e-10
100     enddo
c
	return
	end subroutine emst_rnu
c
c////////////////////////////////////////////////////////////
c
      subroutine rtest(np, ncomp, npcd, phi, isp, msubn, msub, mst,
     $   b, varfn,  alpest)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     by Shankar Z.Ren. S.B.Pope
c     rtest estimates the value of the root to be used in the
c     root finder.
c  input:
c     np : number of particles
c     ncomp : number of compositions
c     npcd : dimensioned number of particles per cell
c     phi : particle composition array : npcd, ncomp
c     isp : particle locator array : npcd
c     msubn :  number of sub-ensembles
c     msub :  particle indices of sub-ensembles : msubn,2
c     mst : MST edge array : npcd-1, 2, msubn
c     b : edge coefficients : (npcd-1, msubn)
c     varfn : variance function for the ensemble
c  output:
c     alpest : estimate of the root
c workspace:
c     work() - real array size n
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer np, ncomp, npcd, msubn, msub(msubn,2),
     $   mst(npcd-1,2,msubn), isp(npcd)
      real phi(npcd,ncomp), varfn,  alpest, b(npcd-1,msubn)

c     Local declarations

      integer idim, iens, iedge, jabs, jnuabs, nk, istart
      real dvfndt, rbb, dphinu

c     1. Find change in variance function

      dvfndt = -varfn

c     2. Find r_beta_beta

      do iens = 1, msubn
        istart = msub(iens,1)
        nk = msub(iens,2)
        rbb = 0.
        do iedge = 1, nk-1
          jabs = istart + mst(iedge,1,iens) - 1
          jnuabs = istart + mst(iedge,2,iens) - 1
            do idim = 1, ncomp
            dphinu = phi(isp(jabs),idim) - phi(isp(jnuabs),idim)
            rbb = rbb + b(iedge,iens)*dphinu*dphinu
            enddo
        enddo
      enddo
      alpest = -0.5*dvfndt/rbb
      return
      end subroutine rtest
c
c////////////////////////////////////////////////////////////
c
      subroutine rtfini(np,  ncomp, phit, phitmp, phivfn,
     $   msubn, msub, isp, w, varfn, dphi0,  dt,
     $   alpest,  efdim, dphi, mst, lc, rs, parent, b, coefi,
     $   coefd, coefe, stack, itrav,  expphi, covtmp, phihat,
     $   q, sigma, dphthr, dphimx,
     $    work, dtm)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    rtfini finds the root to the variance function equation
c    by Z.Ren, S.B.Pope 2002
c    In this version, rejection method is used.

c input:
c     np : number of particles
c     ncomp : number of composition varirbales
c     phit (which is phi in mixemst; compositions in scalar space)
c       : np, ncomp
c     phitmp : particle composition array : np, ncomp
c     phivfn : particle composition array for variance fn calculation
c            : np, ncomp
c     msubn : number of sub-ensembles
c     msub :  particle indices of sub-ensembles : msubn,2
c     isp : cross reference array : np
c     w : array of normalized particle weights : np
c     varfn : variance function for the whole ensemble
c     dphi0 : particle composition change array : np, ncomp
c     dt : normalized time step over which mixing is reqd
c     alpest : estimate of the root
c     efdim : array of effective dimensions
c     dphi : particle composition change array : np
c     mst : MST edge list : np-1, 2, msubn
c     lc : leftmost child array : np
c     rs : right-sibling array : np
c     parent : parent node and connecting edge information : np,2
c     b : edge coefficient array : np-1, msubn
c     coefi : coefficient array : np
c     coefd : coefficient array : np
c     coefe : coefficient array : np
c     stack : temp array for node stacking : np
c     itrav : array for listing traversed nodes : np
c     expphi : means : ncomp
c     covtmp : temp covariance matrix array : ncomp, ncomp
c     phihat : array of stdized compositions : np, ncomp
c     q : orthogonal transformation matrix : ncomp, ncomp
c     sigma : vector of singular values : ncomp
c     dphthr(ncomp) : threshold values for composition changes
c     dphimx(ncomp) : max values of composition change


c output:
c      phit: particle composition array after mixing.
c      dtm : the time step over which mixing was actually performed subject
c            to the constraint on dphi

c workspace:
c     work() - real array size ncomp

c local variables:
c     aa: random number uniformly distributed in[0,1]
c     np : number of particles in the ensemble.
c     npsub : number of particles in the subensemble of interest
c     ntsub : number of particles in that subensemble in the tree
c     notint : number of particles in that subensemble NOT in the tree
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      implicit none

      integer np, ncomp, msubn, msub(msubn,2), isp(np),
     $    efdim(ncomp), mst(np-1, 2, msubn), lc(np),
     $   rs(np), parent(np,2), stack(np), itrav(np)

      real phit(np, ncomp), phitmp(np, ncomp),
     $   phivfn(np, ncomp),
     $   w(np), varfn, dphi0(np, ncomp),  dt,
     $   alpest, dphi(np), b(np-1, msubn), coefi(np),
     $   coefd(np), coefe(np), expphi(ncomp), covtmp(ncomp,ncomp),
     $   phihat(np, ncomp), q(ncomp, ncomp), sigma(ncomp),
     $     dphthr(ncomp), dphimx(ncomp), work(ncomp), dtm


c     local declarations
      integer idim, iens, isrch, npsub, ntsub, notint, i,
     $    accur, sucstp,  istart, iabs, iefdim, it,
     $	  ibrack, mxsrch, mbrack

      real ave, adev, sdev, svar, skew, curt, fh,
     $    eta2,
     $   etal, fl,  etflsp, etainc, varfnd, dvarfn, f,
     $   vfnprv, arg, del, fh0, fhl,  ainc
      real aa(1), dtmtemp

       mxsrch=20
       mbrack=2       ! mxsrch: number of times searching
                      !         in each cycle of bracketing
                      ! mbrack: total number of bracketing
                      !         cycles

c_____________________________________________________________________

c     1. Calculate variance fn and the reqd variance fn

      fh = 0.
      do idim = 1, ncomp
        call moment( 2, np, w, phivfn(1,idim),  ave,
     $     adev, sdev, svar, skew, curt )
        fh = fh + svar
      enddo
      varfn = fh
      vfnprv = fh
      arg = -dt
      varfnd = varfn*exp(arg)
      dvarfn = varfnd - varfn

      do idim = 1, ncomp
        do i = 1, np
          dphi0(i,idim) = 0.
        enddo
        dphimx(idim) = 0.
      enddo

c     3. Initialize root value etc

      ibrack= 0
      fh0   = fh
      sucstp = 0
      accur = 1
      isrch = 0

       eta2=1.08*alpest*dt
       ainc = 1.05
5001  continue

      ibrack = ibrack + 1
      fh = fh0 - varfnd
      etal=eta2

c     4. Try to Bracket the root

  101 call implct(np,np,ncomp,np,msubn,msub,phit,phitmp,
     $   w,isp,mst,lc,rs,parent,b,etal,coefi,coefd,coefe,dphi,
     $   stack,itrav)
c
      do iens = 1, msubn
        istart = msub(iens,1)
        if (iens .le. msubn-1) then
          npsub = msub(iens+1,1) - msub(iens,1)
        else
          npsub = np - msub(iens,1) + 1
        endif
        ntsub = msub(iens,2)
        notint = npsub - ntsub

        do idim = 1, ncomp
          do i = ntsub+1,npsub
            iabs = isp(istart + i - 1)
            phitmp(iabs,idim) = phit(iabs,idim)
          enddo
        enddo
      enddo


        do idim = 1, ncomp
          do i = 1, np
            phivfn(i,idim) = phitmp(i,idim)
          enddo
        enddo

      fl = 0.
c
      do idim = 1, ncomp

        call moment( 2, np, w, phivfn(1,idim),  ave,
     $     adev, sdev, svar, skew, curt )
        fl = fl + svar
      enddo
c
      fl = fl - varfnd
c  mixing is not enough

      if(fl.gt.0.) then
        if(isrch.lt.mxsrch) then
c  avoid fl~fh situation
          if(fl.gt.0.4*fh) then
            etal=1.10*etal
           else
            fhl  = fl / fh
            etal = ainc * etal / (1.0-fhl)
           endif
          isrch = isrch + 1
          goto 101
        elseif(ibrack .le. mbrack) then
         eta2 = 0.99 * eta2
	 ainc = ainc * 1.05
	  goto 5001
        else
         continue
        endif
      endif
c  enough mixing is reached!

c
c

777      continue
           
c     Accumulate dphi0

          it = 0
          do iens = 1, msubn
            istart = msub(iens,1)
            ntsub = msub(iens,2)

            do idim = 1, ncomp
              do i = 1,ntsub
                iabs = isp(istart + i - 1)
                dphi0(it + i,idim) = dphi0(it + i,idim) +
     $             abs( phitmp(iabs,idim) - phit(iabs,idim))
                dphimx(idim) = max(dphimx(idim), dphi0(it+i,idim))
              enddo
            enddo
            it = it + ntsub
          enddo

c    Check dphimx < dphthr
           do idim = 1, ncomp
            if ( dphimx(idim) .gt. dphthr(idim) ) then
              accur = 0
            endif
           enddo


         if ( accur .eq. 1  ) then
                 dtmtemp = log(varfn/(fl+ varfnd))

            if(dtmtemp.gt.dt)then
c   using rejection method.
                 call emst_rnu(aa,1)
                 if(aa(1) .lt. (dt/dtmtemp)) then
                   do idim = 1, ncomp
                   do i = 1, np
                     iabs = isp(i)
                     phit(iabs,idim) = phitmp(iabs,idim)
                   enddo
                   enddo
                 endif
                dtm=dt
             else
                 do idim = 1, ncomp
                 do i = 1, np
                    iabs = isp(i)
                    phit(iabs,idim) = phitmp(iabs,idim)
                 enddo
                 enddo
               dtm=dtmtemp
             endif
            return
         else
c if inaccurate, halve etal.
                 etal=etal/2.0
                 do idim = 1, ncomp
                  do i = 1, np
                  dphi0(i,idim) = 0.
                  enddo
                  dphimx(idim) = 0.
                enddo
                accur = 1

         call implct(np,np,ncomp,np,msubn,msub,phit,phitmp,
     $   w,isp,mst,lc,rs,parent,b,etal,coefi,coefd,coefe,dphi,
     $   stack,itrav)
c
      do iens = 1, msubn
        istart = msub(iens,1)
        if (iens .le. msubn-1) then
          npsub = msub(iens+1,1) - msub(iens,1)
        else
          npsub = np - msub(iens,1) + 1
        endif
        ntsub = msub(iens,2)
        notint = npsub - ntsub

        do idim = 1, ncomp
          do i = ntsub+1,npsub
            iabs = isp(istart + i - 1)
            phitmp(iabs,idim) = phit(iabs,idim)
          enddo
        enddo
      enddo


        do idim = 1, ncomp
          do i = 1, np
            phivfn(i,idim) = phitmp(i,idim)
          enddo
        enddo

      fl = 0.
c
      do idim = 1, ncomp

        call moment( 2,  np, w, phivfn(1,idim),  ave,
     $     adev, sdev, svar, skew, curt )
        fl = fl + svar
      enddo
c
      fl = fl - varfnd
         go to 777

      endif
       return
      end subroutine rtfini
c
c////////////////////////////////////////////////////////////
c

      subroutine sort2(n,ra,ira)
      integer ira(n)
      dimension ra(n)
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
          iira=ira(l)
        else
          rra=ra(ir)
          iira=ira(ir)
          ra(ir)=ra(1)
          ira(ir)=ira(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            ira(1)=iira
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            ira(i)=ira(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
        ira(i)=iira
      go to 10
      end subroutine sort2
c
c////////////////////////////////////////////////////////////
c

      subroutine statinc(n,ds,z0l,z0u,z1l,z1u,state)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  statinc updates the age property
c  input:
c     n - number of particles
c     ds - normalized time increment
c     z0l, z0u - lower and upper bounds of the Z0 r.v.
c     z1l, z1u - lower and upper bounds of the Z1 r.v.

c output:
c     state() - array of states for n particles

c workspace:
c     work() - real array size n
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
	integer, intent(in) :: n
	real(kind(1.e0)), intent(in) :: ds,z0l,z0u,z1l,z1u
	real(kind(1.e0)), intent(inout) :: state(n)

      integer :: i
      real(kind(1.e0)) :: work(n), p0, temp

      p0 = (z1l + z1u)/((z0l+z0u) + (z1l+z1u))

      call emst_rnu(work,n)

      do i = 1, n
        if ( abs(state(i)) <= ds ) then
          if ( state(i) < 0. ) then
            state(i) = z1l + (z1u-z1l)*work(i)
          else
            state(i) = - (z0l + (z0u-z0l)*work(i))
          endif
        else
          if ( state(i) < 0. ) then
            state(i) = state(i) + ds
          else
            state(i) = state(i) - ds
          endif
        endif
      enddo

      return 
      end subroutine statinc  
c
c////////////////////////////////////////////////////////////
c
    
	subroutine stripm( mode, npd, ndim, np, avg, w, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     by Shankar
c     modified to account for variable weights.
c     mode=1: calculate means over np particles and for ndim scalars
c     and strip calculated means from f array for ndim scalars.
c     mode=2: input f array and mean fields for ndim diemnsions;
c     add mean fields to f arrayfor ndim scalars.
c
c input:
c     mode	= 1 - calculate, subtract means from f array for ndim
c     scalars.
c     = 2 - add mean values to f array for ndim scalars.
c     npd	- leading dimension of particle array
c     ndim	- second dimension of particle array
c     np	- number of entries
c     avg     - mean value for ndim scalars (mode 2)
c     f	- particle properties
c output:
c     f	- value of field interpolated onto particle
c     avg     - mean value calculated for ndim scalars (mode 1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real f(npd, ndim), avg(ndim), w(npd)
      double precision    dsum
c

c
      if( mode .eq. 1 ) then
        do j = 1, ndim
          dsum     = 0.d0
          do i   = 1, np
            dsum   = dsum + dble( w(i)*f(i,j) )
          enddo
          avg(j)   = sngl(dsum)
          do i     = 1, np
            f(i,j) = f(i,j) -  avg(j)
          enddo
        enddo
c
      else
        do j = 1, ndim
          do i = 1, np
            f(i,j) = f(i,j) + avg(j)
          enddo
        enddo
      endif
c
      return
      end subroutine stripm
c
c////////////////////////////////////////////////////////////
c
      subroutine subensi(n,ndim,msubn,w,wtmin,isp,msub,randn)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     By Shankar
c     edited Sep 6, 1993 to swap only a pointer array called
c     isp - for subensemble pointer.
c     Form m(n) sub-ensembles of the ensemble of N particles.
c     Given N particles and their weights w(i).
c     Define wmax = max(w(i))
c     Define n = int(1/wmax) : effective number of particles.
c     Let m(n) =  int(log2(n))
c     Randomly order the N particles.
c     Define the kth sub-ensemble as follows:
c     add particles till the following inequality is satisfied
c
c     Wk >= ( 1 - sum(Wj, j=1,k))/( m - k)
c
c     where Wk is the sum of w(i) for all particles in the kth
c     sub-ensemble and the equality means,
c
c     1 = (m-k) Wk + sum( Wj, j=1,k)
c
c
c     N is denoted by the variable n.
c     n is denoted by the variable neff.
c     msubn is the number of sub-ensembles, m.
c     w is the array of weights.
c     msub defines the array of subensembles, msubn in number,
c     by the first entry pointing to the beginning of the sub-ensemble
c     and the second the length of the sub-ensemble in the array of
c     randomly  ordered particles x.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      integer msub(msubn,2),isp(n)
      real randn(n),w(n),wtmin
      double precision sumtot, rhswt

      do i=1,n
        isp(i) = i
      enddo
      if(msubn.eq.1)then
        msub(1,1)=1
        msub(1,2)=n
        return
      else
        call emst_rnu(randn,n)
c____________________________________________________________
c     At the ith chioce of a random number between 1 and N
c     given 0 < r(i) < 1., compute
c     iran = int(r(i)*(N-i-1)) + 1 s.t.  1 < iran < N - i
c     and swap the i + iran th particle with the ith one.
c     s.t. i + 1 < iswap < N
c
c     Note that you can only swap n-1 particles, nth particle must
c     swap with itself.
c____________________________________________________________
	do i=1,n-1
	  iran = int(randn(i)*(n-i-1)) + 1
	  iswap = iran + i
	  itemp  = isp(iswap)
	  isp(iswap) = isp(i)
	  isp(i) = itemp
        enddo
c____________________________________________________________
c
c     Now the array isp(n) contains the indices of the randomly
c     reordered particles.
c     Forming sub-ensembles ..
c     sumtot = total weight of selected particles in the k-1 previous
c     sub-ensembles.
c     sumwk = total weight of the j particles selected in the
c     k-th sub-ensemble.
c____________________________________________________________
	sumtot=0.d0
	jstart = 1
	do k = 1,msubn-1
	  sumwk = 0.d0
	  numpk = 0
	  rhswt = ( 1. - sumtot )/dble(float(msubn-k+1))
	  msub(k,1) = jstart
	  do j=jstart,n
	    if( (sumwk+dble(w(isp(j)))). le. (rhswt+0.5*wtmin) ) then
	      sumwk = sumwk + dble(w(isp(j)))
	      numpk = numpk + 1
	    else
	      goto 40
            endif
          enddo
  40	  msub(k,2) = numpk
	  sumtot = sumtot + sumwk
	  jstart = jstart + numpk
        enddo
	msub(msubn,1) = jstart
	msub(msubn,2) = n - jstart + 1
      endif
      return
      end subroutine subensi
c
c////////////////////////////////////////////////////////////
c
      subroutine subint(n,s,msubn,isp,msub)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     By Shankar
c     subint forms the subensemble containing the particles in the EMST.
c input:
c     n - number of particles
c     s() - state array for the n particles
c     msubn - number of sub-ensembles
c     isp() - particle locator array
c     msub(istart,nk)

c output:
c     msub(istart,nk1) for the one subensemble, nk0 = nk - nk1 + 1
c     isp() - particle locator array
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer n, msubn, isp(n), msub(msubn,2), nk0, nk1, i, itemp, iens,
     $   istart, nk, icur, nlast

      real s(n)

      do iens = 1, msubn
        istart = msub(iens,1)
        if (iens .le. msubn-1) then
          nk = msub(iens+1,1) - msub(iens,1)
        else
          nk = n - msub(iens,1) + 1
        endif
        nk0 = 0
        nk1 = 0
        nlast = nk - nk0
        icur = 1
        do i = 1, nk
          if ( s(isp(istart+icur-1)) .gt. 0. ) then
            nk1 = nk1 + 1
            icur = icur + 1
c     leave isp unchanged
          else
            itemp = isp(istart+nlast-1)
            isp(istart+nlast-1) = isp(istart+icur-1)
            isp(istart+icur-1) = itemp
            nk0 = nk0 + 1
            nlast = nlast - 1
          endif
        enddo
        msub(iens,2) = nk1
      enddo
      return
      end subroutine subint
c
c////////////////////////////////////////////////////////////
c
      subroutine tree( ned, np, nkmax, msubn, msub, mst, lc, rs,
     $   parent, stack, edptr )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     by Shankar
c     Jan 4 1994 : accounts for ned different from np.
c     Sep 2 1993, incorporates changes in mst array :
c     does all the sub-ensembles at one shot.
c     **Now the LC and RS arrays store the leftmost child
c     and rightsibling of the absolute value of node ( i.e.
c     1 < node < n ) in its absolute address location,
c     but the actual value of LC or RS is still a relative
c     address ( and needs istart of that sub-ensemble to
c     get the absolute address).
c     Construct the leftmost-child (LC), right-sibling (RS)
c     representation of the EMST from the edges. Also construct
c     the parent array. : second entry has edge information.
c     Stack is a workspace stack of length n.
c     edptr is a work array.
c     Extra info for model :
c     Distribution of the number of edges.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer ned, np, istart, nk, mst(nkmax-1,2,msubn), lc(ned),
     $   rs(ned), parent(ned,2), stack(ned), edptr(ned),
     $   top, child, ihroot, msub(msubn,2)

c     Initialize all arrays to zero.

      do i = 1, np
        lc(i) = 0
        rs(i) = 0
        parent(i,1) = 0
        parent(i,2) = 0
      enddo
      do iens = 1, msubn
        istart = msub(iens,1)
        nk = msub(iens,2)
        do i = 1,nk
          stack(i) = 0
          edptr(i) = i
        enddo
c     
c     *1	Initialize the node to be the root. Choose root to be
c               the first
c               node in the mst array.
c     
	node = mst(1,1,iens)
	top = 1
	stack(top) = node
	neleft = nk -1
	ieptr = 1
c     
c     *2	Given node search mst for existing edges and find a
c               child.
c     
  2	if(ieptr.le.neleft)then
  	  ie = edptr(ieptr)
	else
	  goto 8
	endif
  	if( node.eq.mst(ie,1,iens) ) then
  	  child = mst(ie,2,iens)
	  itmp = edptr(ieptr)
	  edptr(ieptr) = edptr(neleft)
	  edptr(neleft) = itmp
	  neleft = neleft - 1
	  goto 3
	endif
 	if( node.eq.mst(ie,2,iens) ) then
	  child = mst(ie,1,iens)
	  itmp = edptr(ieptr)
	  edptr(ieptr) = edptr(neleft)
	  edptr(neleft) = itmp
	  neleft = neleft - 1
	  goto 3
	endif
	ieptr = ieptr + 1
	goto 2
c     
c     *3	Found a child, is it the leftmost ?
c     
  3	if( lc(node+istart-1).eq.0 )then
c     
c     *4	If it is then fill it.
c     
	  lc(node+istart-1) = child
	  parent(child+istart-1,1) = node
	  parent(child+istart-1,2) = ie
	  node = child
	  top = top + 1
	  stack(top) = node
	  ieptr = 1
	  goto 2
	else
c     
c     *5	Leftmost child entry exists, look for empty right
c               sibling slots.
c               Make the LC entry the hor. root for the search.
c     
	  ihroot = lc(node+istart-1)	
  6	  if( rs(ihroot+istart-1).eq.0 ) then
c     *6
	    rs(ihroot+istart-1) = child
	    parent(child+istart-1,1) = node
	    parent(child+istart-1,2) = ie
	    node = child
	    top = top + 1
	    stack(top) = node
	    ieptr = 1
	    goto 2
	  else
c     *7	Look down the hor. line
	    ihroot = rs(ihroot+istart-1)
	    goto 6
	  endif
	endif
c     
c     *8	No children left : either hit a leaf or accounted for
c               all children.
c               If any edges are left then
c               Pop the stack and reset node.
c               else quit.
c     
  8	if( neleft.gt.0 )then
	  top = top - 1
	  node = stack(top)
	  ieptr = 1
	  goto 2
	endif
      enddo
      return
      end subroutine tree
c
c////////////////////////////////////////////////////////////
c
	
	end module emst_subs
