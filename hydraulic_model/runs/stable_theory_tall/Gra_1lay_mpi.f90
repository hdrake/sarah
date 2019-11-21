! Karl Helfrich, WHOI, khelfrich@whoi.edu
! solves the single-layer shallow water equations in flux form 
! this version is a direct descendant of the model introduced in Helfrich et al (1999, JFM)
!
!
! compilation commands
! Note the -Wl,... is necessary on Mac OSX due to stack overflow issues. Could also use -heap_arrays 
! option except that this produces a memory leak in v13.0 of the intel compilers
!
! mpif90 -fpp -fpe0 -O3 Gra_1lay_mpi.f90 -o Gra_1lay_mpi.x -Wl,-stack_size,0x10000000,-stack_addr,0xc0000000
!
! ifort -DSERIAL -fpp -fpe0 -O3 Gra_1lay_mpi.f90 -o Gra_1lay_ser.x -Wl,-stack_size,0x10000000,-stack_addr,0xc0000000

! on gfdlab15 (ubuntu under gfortran):
! gfortran -cpp -DSERIAL -ffpe-trap=zero,overflow,underflow -O3 Gra_1lay_mpi.f90 -o Gra_1lay_ser.x
!
! or for MPI
!
! mpif90 -cpp -ffpe-trap=zero,overflow,underflow -O3 Gra_1lay_mpi.f90 -o Gra_1lay_mpi.x



#ifndef SERIAL
#define USEMPI
#endif

!*************************************************************************** 
 module m_parameters
!***************************************************************************
 
 implicit none
 integer, parameter :: fpp = kind(1.d0)   ! set floating point precision here
 
 real(fpp), parameter :: eps=0.00000001_fpp,eps8=0.00000001_fpp
 real(fpp), parameter :: half=0.5_fpp, quarter=0.25_fpp
 real(fpp), parameter :: zero = 0._fpp, one=1._fpp
 real(fpp), parameter :: two=2._fpp, three=3._fpp, four=4._fpp 
 real(fpp), parameter :: five=5._fpp, six=6._fpp, seven=7._fpp
 real(fpp), parameter :: nine=9._fpp, ten=10._fpp, eleven=11._fpp
 real(fpp), parameter :: pi=4._fpp*atan(1._fpp), twopi=2._fpp*pi
 
 integer(kind=4), parameter :: mbc=2, meqn=3
 integer(kind=4) :: nx,ny
 integer(kind=4) :: dims(2) 
 integer(kind=4) :: ixbc, iybc
 integer(kind=4) :: nout, irestrt
 integer(kind=4) :: sx, ex, sy, ey 
 integer(kind=4) :: isxs, isxe, isys, isye
 real(fpp)       :: dx,dy
 real(fpp)       :: xlength,ylength
 real(fpp)       :: dt,t0,tend
 real(fpp)       :: gamma,cf,cb
 real(fpp)       :: u0, xb0, yb0, ay1, ax1, slpy
 real(fpp)       :: alpha, xc, sill, ramp1
 real(fpp)       :: initalpha, initbeta, initgamma
 character*80    :: runname, outfile, restart

 contains
  
   subroutine m_parameters_init  !============================
   
    character*80   :: fname
   
    call get_run_name(runname)  ! reads the runname from the command line
    fname = trim(runname) // '.ic'
    open(10,file=fname,status='unknown',form='formatted') 
      read(10,*)  nx, ny, xlength, ylength
      read(10,*)  dims(1), dims(2)
      read(10,*)  dt, tend, nout, irestrt
      read(10,*)  ixbc, iybc
      read(10,*)  gamma, cf, cb
      read(10,*)  alpha, xc, sill, ramp1
      read(10,*)  initalpha, initbeta, initgamma ! pchan stab
      read(10,'(a)') outfile
      read(10,'(a)') restart
    close(10)

	dx = xlength/real(nx,fpp)
	dy = ylength/real(ny,fpp)
		
  end subroutine m_parameters_init

  subroutine get_run_name(runname)  !=========================
! read the runname from the command line

    character*80, intent(out) :: runname
    character*80 :: tmpstr
    integer      :: iargc, ierr

    if(iargc().eq.0) then
       call getarg(0,tmpstr)
       print *, 'Stopping: correct usage: ',trim(tmpstr),' <run name>'
#ifdef USEMPI
       call MPI_FINALIZE(ierr)
#endif
       stop
    end if
    call getarg(1,runname)    
  end subroutine get_run_name
 
 end module m_parameters
!***************************************************************************
 
!***************************************************************************
 module m_variables
!***************************************************************************
  
  use m_parameters
  
  implicit none
  real(fpp), allocatable :: q(:,:,:)
  real(fpp), allocatable :: gradp(:,:,:),diss(:,:,:)
  real(fpp), allocatable :: uh(:,:), vh(:,:)
  real(fpp), allocatable :: b(:,:)
  real(fpp), allocatable :: cbxy(:,:)
  real(fpp), allocatable :: x(:), y(:)
  real(fpp), allocatable :: qlm1(:,:,:),qlm2(:,:,:),ql(:,:)
  real(fpp), allocatable :: qrm1(:,:,:),qrm2(:,:,:),qr(:,:)
  real(fpp), allocatable :: wrk1(:,:), wrk2(:,:), wrk3(:,:) 
  real(fpp), allocatable :: qwrk1(:,:,:),qwrk2(:,:,:)  

 contains
 
  subroutine m_variables_init
      allocate(q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,meqn))
      allocate(gradp(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2))
      allocate(diss(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2))
      allocate(uh(sx-mbc:ex+mbc,sy-mbc:ey+mbc))
      allocate(vh(sx-mbc:ex+mbc,sy-mbc:ey+mbc))
      allocate(b(sx-mbc:ex+mbc,sy-mbc:ey+mbc))
      allocate(cbxy(sx-mbc:ex+mbc,sy-mbc:ey+mbc))
      allocate(x(1-mbc:nx+mbc), y(1-mbc:ny+mbc))
      allocate(qwrk1(sx-mbc:ex+mbc,sy-mbc:ey+mbc,meqn))
      allocate(qwrk2(sx-mbc:ex+mbc,sy-mbc:ey+mbc,meqn))
      allocate(wrk1(nx,ny), wrk2(nx,ny), wrk3(nx,ny))      
      allocate(qlm1(3,sy-mbc:ey+mbc,meqn),qlm2(3,sy-mbc:ey+mbc,meqn))
      allocate(qrm1(3,sy-mbc:ey+mbc,meqn),qrm2(3,sy-mbc:ey+mbc,meqn))
      allocate(ql(sy-mbc:ey+mbc,meqn),qr(sy-mbc:ey+mbc,meqn))  
      
      q(:,:,:) = zero
      diss(:,:,:) = zero
      b(:,:) = zero
                        
  end subroutine m_variables_init
 
  subroutine m_variables_final 
      deallocate(x)
      deallocate(y)
      deallocate(q)
      deallocate(gradp,diss)
      deallocate(uh,vh,b,cbxy)
      deallocate(qlm1,qlm2,qrm1,qrm2)
      deallocate(ql,qr)
      deallocate(wrk1,wrk2,wrk3)
      deallocate(qwrk1,qwrk2)
  end subroutine m_variables_final
 
 end module m_variables
!***************************************************************************


!***************************************************************************
   program Gra_1lay_mpi
!***************************************************************************

   use m_parameters
   use m_variables

#ifdef USEMPI
   use mpi
#endif
  
   implicit none
   
! MPI-related variables. note some are still used in serial mode. 
      integer(kind=4) :: msgtype,source
      integer(kind=4) :: root,done,myid
      integer(kind=4) :: numprocs,comm2d,ierr,stride 
      integer(kind=4) :: nbrleft,nbrright,nbrbottom,nbrtop
      logical :: periods(2) 
#ifdef USEMPI
      integer(kind=4) :: status(MPI_STATUS_SIZE)
#endif

! local variables
      integer(kind=4) :: iosx, iosy, numx, numy      
      integer(kind=4) :: i, j, k, m, ii
      integer(kind=4) :: ni, no, nimax
      real(fpp) :: t, dtout, tstarto, tendo
      real(fpp) :: tmpi0, tmpi1, tmpi2, tmpi3
      real(fpp) :: tex,tbc,tflx,tvel,tgp,tdis,tsrc
      
!  assign some variables (might be changed later)
      data periods /.false., .false./ 
      data root,done /0,4/
      data tex,tbc,tflx,tvel,tgp,tdis,tsrc /7*zero/      

#ifdef USEMPI
! start the MPI
      call MPI_INIT( ierr ) 
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) 
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr ) 
#else
! or define some things here if in serial mode
      myid=0
      numprocs=1
      dims(:) = 1
#endif

! read the .ic file and initialize the parameters
      call m_parameters_init

           
#ifdef USEMPI
! check on BC specification 
      if (ixbc == 0) periods(1) = .true.
      if (iybc == 0) periods(2) = .true. 
! decomposition: only dims entries = 0 are modified by mpi_dims_create           
      if (dims(1).eq.0 .and. dims(2).ne.0) then ! 1d y-decomposition
          dims(1) = 1
          dims(2) = numprocs   
        elseif (dims(2).eq.0 .and. dims(1).ne.0) then ! 1d x-decomposition
          dims(1) = numprocs  
          dims(2) = 1
        elseif (dims(1)*dims(2) .ne. numprocs) then ! 2d decomposition
          dims(:) = 0
      endif
! Get a new communicator for the decomposition of the domain.  
      call MPI_DIMS_CREATE( numprocs, 2, dims, ierr ) 
      call MPI_CART_CREATE( MPI_COMM_WORLD, 2, dims, periods, &
                            .true., comm2d, ierr ) 
      if (myid == 0) print *, ' decomposition: dims(1:2) = ', dims(1), dims(2)
! Get my position in this communicator   
      call MPI_COMM_RANK( comm2d, myid, ierr ) 
! My neighbors are now +/- 1 with my rank.  
      call get_nbrs2d(comm2d,nbrleft,nbrright,nbrbottom,nbrtop) 
! Compute the decomposition       
      call decomp2d( comm2d, nx, ny, sx, ex, sy, ey )            
! Create a new, "strided" datatype for the exchange in the "non-contiguous" 
! y-direction, including mbc corner halo cells 
! MPI_TYPE_VECTOR(COUNT, BLOCKLENGTH, STRIDE, OLDTYPE, NEWTYPE,IERROR)
!      call mpi_Type_vector( ey-sy+1+2*mbc, 1, ex-sx+1+2*mbc, &
      call mpi_Type_vector( ey-sy+1+2*mbc, mbc, ex-sx+1+2*mbc, &
                 MPI_DOUBLE_PRECISION,stride,ierr)      		                                                                                     
      call mpi_Type_commit( stride, ierr ) 
#else
      sx=1
      ex=nx
      sy=1
      ey=ny
#endif

! allocate data arrays and initialize things used on each processes   
      call m_variables_init 
      
      do i=1-mbc,nx+mbc
	     x(i) = (real(i,fpp)-half) * dx
      enddo
      do j=1-mbc,ny+mbc
         y(j) = (real(j,fpp)-half) * dy
      enddo
      
! identify if this process includes a boundary  
      isxs=0; isxe=0; isys=0; isye=0
      if(sx==1)  isxs=1
      if(ex==nx) isxe=1
      if(sy==1)  isys=1
      if(ey==ny) isye=1

! read the .ic file and initialize q, b, ...
      call ic2dswe_pchan
! if a restart then replace q with q read from the restart file
      if (irestrt == 1) then  
         if (myid==0) call rddump(t0)
#ifdef USEMPI
         call MPI_BCAST(t0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(wrk1,nx*ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(wrk1,nx*ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(wrk3,nx*ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
         q(sx:ex,sy:ey,1) = wrk1(sx:ex,sy:ey)                     
         q(sx:ex,sy:ey,2) = wrk2(sx:ex,sy:ey)                     
         q(sx:ex,sy:ey,3) = wrk3(sx:ex,sy:ey)                     
      endif
      if (ixbc==2) call radinit
             
#ifdef USEMPI
      call MPI_BARRIER( MPI_COMM_WORLD, ierr ) 
      call MPI_BARRIER( comm2d, ierr ) 
#endif
 
! open the output file and write the initial state, etc
! collect b (cbxy) into wrk1 (wrk2) on root for writing
#ifdef USEMPI
      call collect(   b(sx:ex,sy:ey),wrk1,myid,numprocs,comm2d,root,done) 
      call collect(cbxy(sx:ex,sy:ey),wrk2,myid,numprocs,comm2d,root,done) 
#else
      wrk1(sx:ex,sy:ey) =    b(sx:ex,sy:ey)
      wrk2(sx:ex,sy:ey) = cbxy(sx:ex,sy:ey)
#endif
      if(myid == 0) then
         open(11,file=outfile,status='new',form='unformatted')
          write(11) nout, tend
          write(11) nx, ny, meqn, mbc
          write(11) dt, dx, dy, gamma, cf, cb
          write(11) ixbc, iybc
          write(11) x(1:nx)
          write(11) y(1:ny)
          write(11) wrk1
          write(11) wrk2
         close(unit=11)
      endif 
#ifdef USEMPI
      call collect(q(sx:ex,sy:ey,1),wrk1,myid,numprocs,comm2d,root,done)
      call collect(q(sx:ex,sy:ey,2),wrk2,myid,numprocs,comm2d,root,done)
      call collect(q(sx:ex,sy:ey,3),wrk3,myid,numprocs,comm2d,root,done)
#else
      wrk1(sx:ex,sy:ey) = q(sx:ex,sy:ey,1)
      wrk2(sx:ex,sy:ey) = q(sx:ex,sy:ey,2)
      wrk3(sx:ex,sy:ey) = q(sx:ex,sy:ey,3)
#endif
      if(myid == 0) then
        call out(t0)
        print *, ' done writing at T = ', t0
      endif

        
! =========================================================
! MAIN LOOP
! =========================================================     
!     # time increments between writing solution:
      dtout = (tend - t0)/real(nout,fpp)
      nimax = int(dtout/dt + eps)

 
      do no = 1, nout
	    tstarto = (no-1) * dtout + t0
	    tendo   = no*dtout + t0
        do ni = 1, nimax         
           t = tstarto + (ni-1)*dt
           m = (no-1)*nimax + ni    

#ifdef USEMPI             
           do k = 1,meqn
            call exchng2( q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,k), &
                    comm2d,stride,nbrleft,nbrright,nbrbottom,nbrtop ) 
           enddo
#endif           

           call applybcs(1, q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:) )
           
           call pgrad
           
           if(cb.gt.zero) call dissip_bot
           if(cf.gt.zero) call dissip_lat

 ! calculate u,v at t+dt/2
	       call vel(   q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:), &
                   gradp(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:), &
                    diss(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:), &
                      uh(sx-mbc:ex+mbc,sy-mbc:ey+mbc),   &
                      vh(sx-mbc:ex+mbc,sy-mbc:ey+mbc) )

           call qsource(1)

#ifdef USEMPI             
           do k = 1,meqn
            call exchng2( q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,k), &
                    comm2d,stride,nbrleft,nbrright,nbrbottom,nbrtop ) 
           enddo
#endif
            
           call applybcs(2, q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:))

	       call flux_old(q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:), &
                        uh(sx-mbc:ex+mbc,sy-mbc:ey+mbc),   &
                        vh(sx-mbc:ex+mbc,sy-mbc:ey+mbc) )

#ifdef USEMPI             
           do k = 1,meqn
            call exchng2( q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,k), &
                   comm2d,stride,nbrleft,nbrright,nbrbottom,nbrtop ) 
           enddo
#endif        

           call applybcs(2, q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:))

           call pgrad
           
           if(cb.gt.zero) call dissip_bot
           if(cf.gt.zero) call dissip_lat

           call qsource(2)

            if (ixbc==2 .and. (isxs==1 .or. isxe==1)) &
                call radup(q(sx-mbc:ex+mbc,sy-mbc:ey+mbc,:))

!        if(myid == 0) print *, no, ni


        enddo ! ni loop

!        if(myid == 0) print *, ' A '

!  collect q on root processor and write to file        
#ifdef USEMPI             
        call collect(q(sx:ex,sy:ey,1),wrk1,myid,numprocs,comm2d,root,done)
        call collect(q(sx:ex,sy:ey,2),wrk2,myid,numprocs,comm2d,root,done)
        call collect(q(sx:ex,sy:ey,3),wrk3,myid,numprocs,comm2d,root,done)
#else
        wrk1(sx:ex,sy:ey) = q(sx:ex,sy:ey,1)
        wrk2(sx:ex,sy:ey) = q(sx:ex,sy:ey,2)
        wrk3(sx:ex,sy:ey) = q(sx:ex,sy:ey,3)
#endif
        if(myid == 0) then
          call out(tendo)
          print *, ' done writing at t = ', tendo
        endif
                  
      enddo ! no loop
    
      call m_variables_final

#ifdef USEMPI      
      call MPI_Type_free( stride, ierr ) 
      call MPI_Comm_free( comm2d, ierr ) 
      call MPI_FINALIZE(ierr) 
#endif
      
      stop

end program Gra_1lay_mpi

! **************************************************************************** 
       subroutine ic2dswe_pchan
! **************************************************************************** 
!     # Set initial conditions for q.
       use m_parameters
       use m_variables 
       
       implicit none
    
       integer(kind=4) :: i,j,k
       real(fpp) :: yp, yc, hp, rn(2)

       y(:) = y(:) - half*ylength
           
	   CALL RANDOM_SEED     ! compiler uses date and time to generate seed
       
       b(:,:) = zero
       cbxy(:,:) = zero
       q(:,:,1) = eps
       q(:,:,2) = zero
       q(:,:,3) = zero
       
       do j=sy-mbc,ey+mbc
         do i=sx-mbc,ex+mbc
!           # set topography
            b(i,j) = sill*(dexp(-((x(i)-half*xlength)/(0.1*xlength))**2)) +alpha*y(j)**2
            cbxy(i,j) = cb
            CALL RANDOM_NUMBER(HARVEST = rn)	       
	        hp = ramp1*dsqrt(-two*dlog(rn(1)))*dcos(twopi*rn(2))
            if ((initalpha*y(j)**2 - initbeta*y(j) + initgamma + hp) .gt. eps) then
!             # set layer thickness d
              q(i,j,1) = initalpha*y(j)**2 - initbeta*y(j) + initgamma + hp
!             # print *, ' ', q(i,j,1)
!             # set down-channel layer transport v*d (where v is geostrophically balanced)
              q(i,j,2) = -(two*(initalpha + alpha)*y(j) - initbeta)
            endif              
          enddo
        enddo
             
       return
       end subroutine ic2dswe_pchan


! **********************************************************************
	subroutine vel(a,gp,d,ue,ve)
! **********************************************************************
       use m_parameters
       use m_variables, only: qwrk1,qwrk2
       
       implicit none
    
       integer(kind=4) :: i,j
       real(fpp) :: hdtdx,hdtdy,hdt,up,um,vp,vm
       real(fpp), intent(out) :: ue(sx-mbc:ex+mbc,sy-mbc:ey+mbc), &
                                 ve(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
       real(fpp), intent(in) :: a(sx-mbc:ex+mbc,sy-mbc:ey+mbc,meqn), &
                               gp(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2), &
                                d(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2)
     	
	  hdtdx = half*dt/dx
	  hdtdy = half*dt/dy
	  hdt   = half*dt
	  
      qwrk1(:,:,1) = a(:,:,2)/(a(:,:,1)+eps)  ! u
      qwrk2(:,:,1) = a(:,:,3)/(a(:,:,1)+eps)  ! v

	  do j = sy-1,ey+1
	     do i = sx-1,ex+1
	      up = dmax1(qwrk1(i,j,1), zero)*hdtdx
	      um = dmin1(qwrk1(i,j,1), zero)*hdtdx
	      vp = dmax1(qwrk2(i,j,1), zero)*hdtdy
	      vm = dmin1(qwrk2(i,j,1), zero)*hdtdy
	      qwrk1(i,j,2) = qwrk1(i,j,1) &
	                       -up*(qwrk1(i,j,1)  -qwrk1(i-1,j,1)) &
                           -um*(qwrk1(i+1,j,1)-qwrk1(i,j,1))   &
                           -vp*(qwrk1(i,j,1)  -qwrk1(i,j-1,1)) &
                           -vm*(qwrk1(i,j+1,1)-qwrk1(i,j,1))   &
                           +hdt*( gamma*qwrk2(i,j,1)-gp(i,j,1) &
                                +d(i,j,1)/(a(i,j,1)+eps) ) 
	      qwrk2(i,j,2) = qwrk2(i,j,1) &
	                       -up*(qwrk2(i,j,1)  -qwrk2(i-1,j,1)) &
                           -um*(qwrk2(i+1,j,1)-qwrk2(i,j,1))   &
                           -vp*(qwrk2(i,j,1)  -qwrk2(i,j-1,1)) &
                           -vm*(qwrk2(i,j+1,1)-qwrk2(i,j,1))   &
                           +hdt*(-gamma*qwrk1(i,j,1)-gp(i,j,2) &
                                +d(i,j,2)/(a(i,j,1)+eps) )
         enddo
      enddo
            
	  do j = sy-1,ey+1
	     do i = sx-1,ex+1
	       ue(i,j) = half*(qwrk1(i,j,2) + qwrk1(i-1,j,2))
	       ve(i,j) = half*(qwrk2(i,j,2) + qwrk2(i,j-1,2))
         enddo
      enddo

!  fix u.n = 0 for no flux boundary conditions if necessary
       if (ixbc==1 .and. isxs==1) ue(sx,:)   = zero
       if (ixbc==1 .and. isxe==1) ue(ex+1,:) = zero
       if (iybc==1 .and. isys==1) ve(:,sy)   = zero
       if (iybc==1 .and. isye==1) ve(:,ey+1) = zero

  	return
  	end subroutine vel


! **********************************************************************
	subroutine vel_hold(a,gp,d,ue,ve)
! **********************************************************************
       use m_parameters
       use m_variables, only: qwrk1
       
       implicit none
    
       integer(kind=4) :: i,j
       real(fpp) :: hdtdx,hdtdy,hdt,up,um,vp,vm
       real(fpp), intent(out) :: ue(sx-mbc:ex+mbc,sy-mbc:ey+mbc), &
                                 ve(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
       real(fpp), intent(in) :: a(sx-mbc:ex+mbc,sy-mbc:ey+mbc,meqn), &
                               gp(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2), &
                                d(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2)
       real(fpp) :: u1(sx-mbc:ex+mbc,sy-mbc:ey+mbc), &
                    v1(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
       real(fpp) :: uu(sx-mbc:ex+mbc,sy-mbc:ey+mbc), &
                    vv(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
     	
	  hdtdx = half*dt/dx
	  hdtdy = half*dt/dy
	  hdt   = half*dt
	  
      u1(:,:) = a(:,:,2)/(a(:,:,1)+eps)
      v1(:,:) = a(:,:,3)/(a(:,:,1)+eps)
      uu(:,:)=zero
      vv(:,:)=zero

	  do j = sy-1,ey+1
	     do i = sx-1,ex+1
	      up = dmax1(u1(i,j), zero)*hdtdx
	      um = dmin1(u1(i,j), zero)*hdtdx
	      vp = dmax1(v1(i,j), zero)*hdtdy
	      vm = dmin1(v1(i,j), zero)*hdtdy
	      uu(i,j) = u1(i,j)  &
                         -up*(u1(i,j)  -u1(i-1,j)) &
                         -um*(u1(i+1,j)-u1(i,j))   &
                         -vp*(u1(i,j)  -u1(i,j-1)) &
                         -vm*(u1(i,j+1)-u1(i,j))   &
                        + hdt*( gamma*v1(i,j)-gp(i,j,1) &
                                +d(i,j,1)/(a(i,j,1)+eps)) 
	      vv(i,j) = v1(i,j) &
                         -up*(v1(i,j)  -v1(i-1,j))  &
                         -um*(v1(i+1,j)-v1(i,j))    &
                         -vp*(v1(i,j)  -v1(i,j-1))  &
                         -vm*(v1(i,j+1)-v1(i,j))    &
                         + hdt*(-gamma*u1(i,j)-gp(i,j,2) &
                                +d(i,j,2)/(a(i,j,1)+eps))
         enddo
      enddo
            
	  do j = sy-1,ey+1
	     do i = sx-1,ex+1
	       ue(i,j) = half*(uu(i,j) + uu(i-1,j))
	       ve(i,j) = half*(vv(i,j) + vv(i,j-1))
         enddo
      enddo

!  fix u.n = 0 for no flux boundary conditions if necessary
       if (ixbc==1 .and. isxs==1) ue(sx,:) = zero
       if (ixbc==1 .and. isxe==1) ue(ex+1,:) = zero
       if (iybc==1 .and. isys==1) ve(:,sy) = zero
       if (iybc==1 .and. isye==1) ve(:,ey+1) = zero

  	return
  	end subroutine vel_hold


! ***********************************************************************
	subroutine flux_old(a,u,v)
! ***********************************************************************

!  Compute the nonlinear flux differences F_x = (uq)_x 
!                                     and G_y = (vq)_y
!  and update q(n+1) = q(n) - dt*[F_x + G_y]

!  This version follows R. J. LeVeque (1996, 
!  SIAM J. Numer. Anal. 33(2), 627-665). 
       use m_parameters
       use m_variables, only: qwrk1,qwrk2
       
       implicit none
    
 	   integer(kind=4), parameter :: ilim = 1
	   real(fpp), parameter :: eps2 = 0.0001_fpp

       real(fpp), intent(in) :: u(sx-mbc:ex+mbc,sy-mbc:ey+mbc), &
                                v(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
       real(fpp), intent(inout) :: a(sx-mbc:ex+mbc,sy-mbc:ey+mbc,meqn)
 
       integer(kind=4) :: i,j,m,isu,I2,I3,I4,J2,J3,J4
       real(fpp) :: philim
       real(fpp) :: hdt,dtdx,dtdy,hdtdx,hdtdy
       real(fpp) :: c,S,S1,waveud
       real(fpp) :: vpi,vmi,vpim,vmim,vp3,vm3
       real(fpp) :: upj,umj,upjm,umjm,up3,um3
       real(fpp) :: wave(meqn)
	   
	   dtdy  = dt/dy
	   dtdx  = dt/dx
	   hdtdy = half*dt/dy
	   hdtdx = half*dt/dx

       qwrk1(:,:,:) = zero  ! F x-directed flux
       qwrk2(:,:,:) = zero  ! G y-directed flux

! x-direction second-order with transverse corrections
          do j = sy, ey
            do i = sx, ex+1
              vpi  = dmax1(v(i,j+1),   zero)
              vmi  = dmin1(v(i,j),     zero)
              vpim = dmax1(v(i-1,j+1), zero)
              vmim = dmin1(v(i-1,j),   zero)              
!              isu = idnint(sign(one,u(i,j)))
!              I2  = i-1+(1-isu)/2
!              I3  = i-(1-isu)/2
!              I4  = i-isu
!              vp3 = dmax1(v(I3,j+1),   zero)
!              vm3 = dmin1(v(I3,j),     zero)              
              if(u(i,j).ge.zero) then
                 I2 = i-1
                 I3 = i
                 I4 = i-1
                 vp3 = vpi
                 vm3 = vmi
                else
                 I2 = i
                 I3 = i-1
                 I4 = i+1
                 vp3 = vpim
                 vm3 = vmim
              endif

! first-order with transverse corrections
              do m = 1, meqn
                 wave(m) = a(i,j,m) - a(i-1,j,m)
                 qwrk1(i,j,m) = qwrk1(i,j,m) + u(i,j)*a(I2,j,m)
                 c = hdtdx*u(i,j)*wave(m)
                 qwrk2(I3,j+1,m) = qwrk2(I3,j+1,m)-c*vp3
                 qwrk2(I3,j  ,m) = qwrk2(I3,j  ,m)-c*vm3
	          enddo
! check if layer depth near zero.  If not, use 2st-order scheme
              if(a(i,j,1).gt.eps2.and.a(i-1,j,1).gt.eps2) then
		        S1 = half*dabs(u(i,j))*(one-dtdx*dabs(u(i,j)))
                do m = 1, meqn
                  waveud = a(I4,j,m) - a(I4-1,j,m)
                  S = S1*philim(wave(m),waveud,ilim)*wave(m)
                  qwrk1(i,j,m) = qwrk1(i,j,m) + S
                  qwrk2(i  ,j+1,m) = qwrk2(i  ,j+1,m) + vpi  * dtdx*S
                  qwrk2(i-1,j+1,m) = qwrk2(i-1,j+1,m) - vpim * dtdx*S
                  qwrk2(i  ,j  ,m) = qwrk2(i  ,j  ,m) + vmi  * dtdx*S
                  qwrk2(i-1,j  ,m) = qwrk2(i-1,j  ,m) - vmim * dtdx*S
                enddo
	          endif
            enddo
          enddo

! y-direction second-order with transverse corrections
          do j = sy,ey+1
            do i = sx,ex
              upj  = dmax1(u(i+1,j),   zero)
              umj  = dmin1(u(i,j),     zero)
              upjm = dmax1(u(i+1,j-1), zero)
              umjm = dmin1(u(i,j-1),   zero)              
    !          isu = idnint(sign(one,v(i,j)))
    !          J2  = j-1+(1-isu)/2
    !          J3  = j-(1-isu)/2
    !          J4  = j-isu
    !          up3 = dmax1(u(i+1,J3), zero)
    !          um3 = dmin1(u(i,  J3), zero)
              if(v(i,j).ge.zero) then
                 J2 = j-1
                 J3 = j
                 J4 = j-1
                 up3 = upj
                 um3 = umj
                else
                 J2 = j
                 J3 = j-1
                 J4 = j+1
                 up3 = upjm
                 um3 = umjm
              endif

              do m = 1, meqn
                 wave(m) = a(i,j,m) - a(i,j-1,m)
                 qwrk2(i,j,m) = qwrk2(i,j,m) + v(i,j)*a(i,J2,m)
                 c = hdtdy*v(i,j)*wave(m)
                 qwrk1(i+1,J3,m) = qwrk1(i+1,J3,m)-c*up3
                 qwrk1(i  ,J3,m) = qwrk1(i  ,J3,m)-c*um3
              enddo
              if(a(i,j,1).gt.eps2.and.a(i,j-1,1).gt.eps2) then
                S1 = half*dabs(v(i,j))*(one-dtdy*dabs(v(i,j))) 
                do m = 1, meqn             
                   waveud = a(i,J4,m) - a(i,J4-1,m)
                   S = S1*philim(wave(m),waveud,ilim)*wave(m)
                   qwrk2(i,j,m) = qwrk2(i,j,m) + S
                   qwrk1(i+1,j  ,m) = qwrk1(i+1,j  ,m) + upj  * dtdy*S
                   qwrk1(i+1,j-1,m) = qwrk1(i+1,j-1,m) - upjm * dtdy*S
                   qwrk1(i  ,j  ,m) = qwrk1(i  ,j  ,m) + umj  * dtdy*S
                   qwrk1(i  ,j-1,m) = qwrk1(i  ,j-1,m) - umjm * dtdy*S
		        enddo
              endif
            enddo
          enddo
        
        do m = 1,meqn
          do j = sy,ey
            do i = sx,ex
              a(i,j,m) = a(i,j,m) &
                         - dtdx*(qwrk1(i+1,j,m)-qwrk1(i,j,m)) &
                         - dtdy*(qwrk2(i,j+1,m)-qwrk2(i,j,m))
        enddo; enddo; enddo
!        a(sx:ex,sy:ey,1:meqn) = a(sx:ex,sy:ey,1:meqn) &
!            - dtdx*(qwrk1(sx+1:ex+1,sy:ey,1:meqn)-qwrk1(sx:ex,sy:ey,1:meqn)) &
!            - dtdy*(qwrk2(sx:ex,sy+1:ey+1,1:meqn)-qwrk2(sx:ex,sy:ey,1:meqn))
  
        return
        end subroutine flux_old


! **********************************************************************
      subroutine pgrad
! **********************************************************************
! scaled with L = sqrt(gH)/f
       use m_parameters
       use m_variables, only: b,q,gradp
       
       implicit none
    
       integer(kind=4) :: i,j,k
       real(fpp) :: cx,cy,c1,gp,gm
       real(fpp) :: eta(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
      
      cx = half/dx
      cy = half/dy
      
	  eta(:,:) = q(:,:,1)+b(:,:)
      
      do j = sy-1, ey+1
         do  i = sx-1, ex+1
!  x-derivative
            c1 = two/(q(i-1,j,1)+q(i+1,j,1)+two*eps)
            gp = q(i+1,j,1)*c1
            gm = q(i-1,j,1)*c1
       	    gradp(i,j,1) = cx*(gp*(eta(i+1,j)-eta(i,j))+gm*(eta(i,j)-eta(i-1,j)))
!  y-derivative
	        c1 = two/(q(i,j-1,1)+q(i,j+1,1)+two*eps)
            gp = q(i,j+1,1)*c1
            gm = q(i,j-1,1)*c1
       	    gradp(i,j,2) = cy*(gp*(eta(i,j+1)-eta(i,j))+gm*(eta(i,j)-eta(i,j-1)))
         enddo
      enddo

!  fix for no-flux along y = ymin, ymax
      if (iybc==1 .and. isys == 1) then
       do i = sx-1,ex+1
	     gradp(i,1,2) = cy*(-three*eta(i,1)+four*eta(i,2)-eta(i,3))
       enddo
      endif
      if (iybc==1 .and. isye== 1) then
       do i = sx-1,ex+1
	     gradp(i,ny,2) = cy*( three*eta(i,ny)-four*eta(i,ny-1)+eta(i,ny-2))
       enddo
      endif
      
!  fix for no-flux along x = xmin, xmax
      if (ixbc == 1 .and. isxs==1) then
        do j = sy-1,ey+1
	      gradp(1,j,1) = cx*(-three*eta(1,j)+four*eta(2,j)-eta(3,j))
        enddo
      endif
      if (ixbc == 1 .and. isxe==1) then
        do j = sy-1,ey+1
	      gradp(nx,j,1) = cx*( three*eta(nx,j)-four*eta(nx-1,j)+eta(nx-2,j))
        enddo
      endif
      
      return
      end subroutine pgrad

! **********************************************************************
      subroutine pgrad_hold(h1,gradp)
! **********************************************************************
! scaled with L = sqrt(gH)/f
       use m_parameters
       use m_variables, only: b
       
       implicit none
    
       integer(kind=4) :: i,j,k
       real(fpp) :: cx,cy,c1,gp,gm
       real(fpp) :: eta(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
       real(fpp), intent(in)  :: h1(sx-mbc:ex+mbc,sy-mbc:ey+mbc)
       real(fpp), intent(out) :: gradp(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2)
      
      cx = half/dx
      cy = half/dy
      
	  eta(:,:) = h1(:,:)+b(:,:)
      
      do j = sy-1, ey+1
         do  i = sx-1, ex+1
!  x-derivative
            c1 = two/(h1(i-1,j)+h1(i+1,j)+two*eps)
            gp = h1(i+1,j)*c1
            gm = h1(i-1,j)*c1
       	    gradp(i,j,1) = cx*(gp*(eta(i+1,j)-eta(i,j))+gm*(eta(i,j)-eta(i-1,j)))
!  y-derivative
	        c1 = two/(h1(i,j-1)+h1(i,j+1)+two*eps)
            gp = h1(i,j+1)*c1
            gm = h1(i,j-1)*c1
       	    gradp(i,j,2) = cy*(gp*(eta(i,j+1)-eta(i,j))+gm*(eta(i,j)-eta(i,j-1)))
         enddo
      enddo

!  fix for no-flux along y = ymin, ymax
      if (iybc==1 .and. isys == 1) then
       do i = sx-1,ex+1
	     gradp(i,1,2) = cy*(-three*eta(i,1)+four*eta(i,2)-eta(i,3))
       enddo
      endif
      if (iybc==1 .and. isye== 1) then
       do i = sx-1,ex+1
	     gradp(i,ny,2) = cy*( three*eta(i,ny)-four*eta(i,ny-1)+eta(i,ny-2))
       enddo
      endif
      
!  fix for no-flux along x = xmin, xmax
      if (ixbc == 1 .and. isxs==1) then
        do j = sy-1,ey+1
	      gradp(1,j,1) = cx*(-three*eta(1,j)+four*eta(2,j)-eta(3,j))
        enddo
      endif
      if (ixbc == 1 .and. isxe==1) then
        do j = sy-1,ey+1
	      gradp(nx,j,1) = cx*( three*eta(nx,j)-four*eta(nx-1,j)+eta(nx-2,j))
        enddo
      endif
      
      return
      end subroutine pgrad_hold

! ***********************************************************************
      subroutine dissip_bot
! ***********************************************************************
       use m_parameters
       use m_variables, only: cbxy,q,diss 
       
       implicit none
    
       integer(kind=4) :: i,j
       real(fpp) :: u1,v1,s1
      
        do j = sy-1, ey+1
	      do i = sx-1, ex+1
           u1 = q(i,j,2)/(q(i,j,1)+eps)
           v1 = q(i,j,3)/(q(i,j,1)+eps)
	       s1 = sqrt(u1**2+v1**2)
	       if (q(i,j,1) .ge. 0.001_fpp) then	 
	         diss(i,j,1) = -cbxy(i,j)*u1*s1
	         diss(i,j,2) = -cbxy(i,j)*v1*s1
	       endif	   	 
          enddo
        enddo
      
      return
      end subroutine dissip_bot
      
! ***********************************************************************
      subroutine dissip_lat
! ***********************************************************************
       use m_parameters
       use m_variables, only: q,diss 
       
       implicit none
    
       integer(kind=4) :: i,j
       real(fpp) :: cx,cy,h
       real(fpp) :: cip,cim,cjp,cjm
       real(fpp) :: fpi,fmi,fpj,fmj
       
       if(cb==zero) diss(:,:,:)=zero
      
        cx = cf/dx**2
        cy = cf/dy**2     
        do j = sy-1, ey+1
	      do i = sx-1, ex+1
	 	 
	      h   = q(i,j,1)/(q(i,j,1)+eps)
	      fpi = q(i+1,j,1)/(q(i+1,j,1)+eps)*h
	      fmi = q(i-1,j,1)/(q(i-1,j,1)+eps)*h
	      fpj = q(i,j+1,1)/(q(i,j+1,1)+eps)*h
	      fmj = q(i,j-1,1)/(q(i,j-1,1)+eps)*h
	      cip = (q(i+1,j,1) - q(i,j,1)  ) / (q(i+1,j,1)+q(i,j,1)+eps)
	      cim = (q(i,j,1)   - q(i-1,j,1)) / (q(i-1,j,1)+q(i,j,1)+eps)
	      cjp = (q(i,j+1,1) - q(i,j,1)  ) / (q(i,j+1,1)+q(i,j,1)+eps)
	      cjm = (q(i,j,1)   - q(i,j-1,1)) / (q(i,j-1,1)+q(i,j,1)+eps)
	 
	      diss(i,j,1) = diss(i,j,1) + &
           cx*(fpi*(q(i+1,j,2)-q(i,j,2) - (q(i+1,j,2)+q(i,j,2))*cip)  &
              -fmi*(q(i,j,2)-q(i-1,j,2) - (q(i-1,j,2)+q(i,j,2))*cim)) &
         + cy*(fpj*(q(i,j+1,2)-q(i,j,2) - (q(i,j+1,2)+q(i,j,2))*cjp)  &
              -fmj*(q(i,j,2)-q(i,j-1,2) - (q(i,j-1,2)+q(i,j,2))*cjm))
               
	      diss(i,j,2) = diss(i,j,2) + &
           cx*(fpi*(q(i+1,j,3)-q(i,j,3) - (q(i+1,j,3)+q(i,j,3))*cip)  &
              -fmi*(q(i,j,3)-q(i-1,j,3) - (q(i-1,j,3)+q(i,j,3))*cim)) &
         + cy*(fpj*(q(i,j+1,3)-q(i,j,3) - (q(i,j+1,3)+q(i,j,3))*cjp)  &
              -fmj*(q(i,j,3)-q(i,j-1,3) - (q(i,j-1,3)+q(i,j,3))*cjm))   
	 
          enddo
        enddo
      
      return
      end subroutine dissip_lat

! ***********************************************************************
      subroutine dissip_hold(a,d)
! ***********************************************************************
       use m_parameters
       use m_variables, only: cbxy 
       
       implicit none
    
       integer(kind=4) :: i,j
       real(fpp) :: u1,v1,s1
       real(fpp) :: cx,cy,h,cip,cim,cjp,cjm,fpi,fmi,fpj,fmj
       real(fpp), intent(in)  :: a(sx-mbc:ex+mbc,sy-mbc:ey+mbc,meqn)
       real(fpp), intent(out) :: d(sx-mbc:ex+mbc,sy-mbc:ey+mbc,2)
      
      d(:,:,:) = zero
      
      if (cb .gt. zero) then                 
        do j = sy-1, ey+1
	      do i = sx-1, ex+1
           u1 = a(i,j,2)/(a(i,j,1)+eps)
           v1 = a(i,j,3)/(a(i,j,1)+eps)
	       s1 = sqrt(u1**2+v1**2)
	       if (a(i,j,1) .ge. 0.01_fpp) then	 
	         d(i,j,1) = -cbxy(i,j)*u1*s1
	         d(i,j,2) = -cbxy(i,j)*v1*s1
	       endif	   	 
          enddo
        enddo
      endif

      if (cf .gt. zero) then
        cx = cf/dx**2
        cy = cf/dy**2     
        do j = sy-1, ey+1
	      do i = sx-1, ex+1
	 	 
	      h   = a(i,j,1)/(a(i,j,1)+eps)
	      fpi = a(i+1,j,1)/(a(i+1,j,1)+eps)*h
	      fmi = a(i-1,j,1)/(a(i-1,j,1)+eps)*h
	      fpj = a(i,j+1,1)/(a(i,j+1,1)+eps)*h
	      fmj = a(i,j-1,1)/(a(i,j-1,1)+eps)*h
	      cip = (a(i+1,j,1) - a(i,j,1)  ) / (a(i+1,j,1)+a(i,j,1)+eps)
	      cim = (a(i,j,1)   - a(i-1,j,1)) / (a(i-1,j,1)+a(i,j,1)+eps)
	      cjp = (a(i,j+1,1) - a(i,j,1)  ) / (a(i,j+1,1)+a(i,j,1)+eps)
	      cjm = (a(i,j,1)   - a(i,j-1,1)) / (a(i,j-1,1)+a(i,j,1)+eps)
	 
	      d(i,j,1) = d(i,j,1) + &
           cx*(fpi*(a(i+1,j,2)-a(i,j,2) - (a(i+1,j,2)+a(i,j,2))*cip)  &
              -fmi*(a(i,j,2)-a(i-1,j,2) - (a(i-1,j,2)+a(i,j,2))*cim)) &
         + cy*(fpj*(a(i,j+1,2)-a(i,j,2) - (a(i,j+1,2)+a(i,j,2))*cjp)  &
              -fmj*(a(i,j,2)-a(i,j-1,2) - (a(i,j-1,2)+a(i,j,2))*cjm))
               
	      d(i,j,2) = d(i,j,2) + &
           cx*(fpi*(a(i+1,j,3)-a(i,j,3) - (a(i+1,j,3)+a(i,j,3))*cip)  &
              -fmi*(a(i,j,3)-a(i-1,j,3) - (a(i-1,j,3)+a(i,j,3))*cim)) &
         + cy*(fpj*(a(i,j+1,3)-a(i,j,3) - (a(i,j+1,3)+a(i,j,3))*cjp)  &
              -fmj*(a(i,j,3)-a(i,j-1,3) - (a(i,j-1,3)+a(i,j,3))*cjm))   
	 
          enddo
        enddo
      endif
      
      return
      end subroutine dissip_hold


! ************************************************************************
      subroutine qsource(itime)
! ************************************************************************
       use m_parameters
       use m_variables, only: qwrk1,q,gradp,diss 
       
       implicit none
    
       integer(kind=4), intent(in) :: itime
       real(fpp) :: hdt,c
      
       hdt = half*dt

       if (itime==1) then
         qwrk1(sx:ex,sy:ey,1) = q(sx:ex,sy:ey,2)
         qwrk1(sx:ex,sy:ey,2) = q(sx:ex,sy:ey,3)
         q(sx:ex,sy:ey,2) = q(sx:ex,sy:ey,2)+ &
            hdt*( gamma*qwrk1(sx:ex,sy:ey,2)- &
            q(sx:ex,sy:ey,1)*gradp(sx:ex,sy:ey,1)+diss(sx:ex,sy:ey,1)) 
       	 q(sx:ex,sy:ey,3) = q(sx:ex,sy:ey,3)+ &
            hdt*(-gamma*qwrk1(sx:ex,sy:ey,1)- &
            q(sx:ex,sy:ey,1)*gradp(sx:ex,sy:ey,2)+diss(sx:ex,sy:ey,2))
        else
         c = one/(one + (hdt*gamma)**2)
         qwrk1(sx:ex,sy:ey,1) = q(sx:ex,sy:ey,2)+ &
            hdt*(-q(sx:ex,sy:ey,1)*gradp(sx:ex,sy:ey,1)+diss(sx:ex,sy:ey,1))
         qwrk1(sx:ex,sy:ey,2) = q(sx:ex,sy:ey,3)+ &
            hdt*(-q(sx:ex,sy:ey,1)*gradp(sx:ex,sy:ey,2)+diss(sx:ex,sy:ey,2))
         q(sx:ex,sy:ey,2) = c*(qwrk1(sx:ex,sy:ey,1)+hdt*gamma*qwrk1(sx:ex,sy:ey,2))
         q(sx:ex,sy:ey,3) = c*(qwrk1(sx:ex,sy:ey,2)-hdt*gamma*qwrk1(sx:ex,sy:ey,1))
       endif
       
      return
      end subroutine qsource


! ********************************************************************************  
	subroutine applybcs( itime, a )
! ******************************************************************************** 
! periodic condtions (ixbc or iybc = 0) are handled by the halo exchange 
! when using MPI and here when in serial mode
! ixbc, iybc = 0  -- periodic
! ixbc, iybc = 1   --  no flux (wall)
! ixbc = 2  --  radiation (Orlanski)

    use m_parameters
    use m_variables, only : qlm1,qlm2,qrm1,qrm2,ql,qr
    
    implicit none
    integer(kind=4), intent(in) :: itime
    integer(kind=4) :: i,j,k,m
    real(fpp) :: anum, den, c
    real(fpp), intent(inout) :: a(sx-mbc:ex+mbc,sy-mbc:ey+mbc,1:meqn)
 
#ifdef SERIAL
 ! x and y periodic boundary conditions (serial version only) 
    if (ixbc==0) then    
         do j = sy-mbc, ey+mbc
            do i = 1,mbc
	          a(sx-i,j,1) = a(ex+1-i,j,1)
	          a(ex+i,j,1) = a(i,j,1)
	          a(sx-i,j,2) = a(ex+1-i,j,2)
	          a(ex+i,j,2) = a(i,j,2)
	          a(sx-i,j,3) = a(ex+1-i,j,3)
	          a(ex+i,j,3) = a(i,j,3)
            enddo
         enddo
    endif
    if (iybc==0) then  ! ystart	  do k = 1,meqn
	    do j = 1,mbc
          do i = sx-mbc, ex+mbc
	         a(i,sy-j,1) = a(i,ey+1-j,1)
	         a(i,ey+j,1) = a(i,j,1)
	         a(i,sy-j,2) = a(i,ey+1-j,2)
	         a(i,ey+j,2) = a(i,j,2)
	         a(i,sy-j,3) = a(i,ey+1-j,3)
	         a(i,ey+j,3) = a(i,j,3)
          enddo
        enddo
    endif
#endif
 

 ! y no-flux boundary conditions 
    if (isys==1 .and. iybc==1) then  ! ystart
      do i = sx,ex
	     do k = 1,mbc
	       a(i,1-k,1)  =  a(i,k,1)
	       a(i,1-k,2)  =  a(i,k,2)
	       a(i,1-k,3)  = -a(i,k,3)
         enddo
      enddo
    endif
    if (isye==1 .and. iybc==1) then  ! yend
      do i = sx,ex
	     do k = 1,mbc
	       a(i,ey+k,1) =  a(i,ey+1-k,1)
	       a(i,ey+k,2) =  a(i,ey+1-k,2)
	       a(i,ey+k,3) = -a(i,ey+1-k,3)
         enddo
      enddo
    endif
   
! x boundary conditions 
    if (isxs==1) then    ! xstart 
    
      if (ixbc==1) then  ! no-flux
        do j = sy-mbc,ey+mbc
	      do k = 1,mbc
	       a(sx-k,j,1)  =  a(k,j,1)
	       a(sx-k,j,2)  = -a(k,j,2)
	       a(sx-k,j,3)  =  a(k,j,3)
          enddo
        enddo
      endif     
      
      if (ixbc==2 .and. itime==1) then ! radiation
        do j = sy-mbc,ey+mbc

!         do 70 k=1,meqn
!           anum=qlm2(2,j,k)-a(1,j,k)
!           den=a(1,j,k)+qlm2(2,j,k)-two*qlm1(3,j,k)
!           if(den.eq.zero)then
!              ql(j,k) = qlm1(1,j,k)
!              go to 70
!           endif
!           c=anum/den
!           if(c.le.zero)then
!              c = zero
!             elseif(c.gt.one)then
!              c = one
!           endif
!           ql(j,k)=(qlm1(1,j,k)*(one-c)+two*c*a(1,j,k))/(one+c) 
!   70    continue

          do k=1,meqn
            c = (qlm2(2,j,k)-a(sx,j,k)) / &
               ( a(sx,j,k)+qlm2(2,j,k)-two*qlm1(3,j,k) + eps8)
            c = dmax1(zero,dmin1(one,c))
            ql(j,k)=(qlm1(1,j,k)*(one-c)+two*c*a(sx,j,k))/(one+c) 
          enddo 
       
        enddo ! j
        qlm2(:,:,:) = qlm1(:,:,:)
        qlm1(:,:,:) = a(sx-1:sx+1,:,:)
       endif  ! ixbc==2
          
    endif  ! isxs
    
    if (isxe==1) then    ! xend
    
      if (ixbc==1) then  ! no-flux
        do j = sy-mbc,ey+mbc
	     do k = 1,mbc
	       a(ex+k,j,1) =  a(ex+1-k,j,1)
	       a(ex+k,j,2) = -a(ex+1-k,j,2)
	       a(ex+k,j,3) =  a(ex+1-k,j,3)
         enddo
        enddo
      endif
          
      if (ixbc==2 .and. itime==1) then  ! radiation
        do j = sy-mbc,ey+mbc

!         do 50 k=1,meqn    
!           anum = qrm2(2,j,k)-a(nx,j,k)
!           den = a(nx,j,k)+qrm2(2,j,k)-2.d0*qrm1(1,j,k) 
!           if(den.eq.one)then
!              qr(j,k) = qrm1(3,j,k)
!              go to 50
!           endif
!           c=anum/den
!           if(c.le.zero)then
!              c = 0.d0
!             elseif(c.gt.one)then
!              c = 1.d0
!           endif
!           qr(j,k)=(qrm1(3,j,k)*(one-c)+two*c*a(nx,j,k))/(one+c) 
!   50    continue 
 
          do k=1,meqn    
            c = (qrm2(2,j,k)-a(ex,j,k)) / &
               ( a(ex,j,k)+qrm2(2,j,k)-two*qrm1(1,j,k)+eps8 ) 
            c = dmax1(zero,dmin1(one,c))
            qr(j,k)=(qrm1(3,j,k)*(one-c)+two*c*a(ex,j,k))/(one+c) 
          enddo 
        enddo ! j          
        qrm2(:,:,:)   = qrm1(:,:,:)
        qrm1(1:3,:,:) = a(ex-1:ex+1,:,:)  

      endif  ! isxe
        
    endif    
    
    return            
    end subroutine applybcs    

! ********************************************************************************  
	subroutine radinit
! ******************************************************************************** 
    use m_parameters
    use m_variables, only: q,qlm1,qlm2,qrm1,qrm2
    
    implicit none
    integer(kind=4) :: i,j,k
    
 !  set up radiation bc info
    if (isxs==1) then
       do k=1,meqn
          do j=sy-mbc,ey+mbc
             do i=1,3
               qlm1(i,j,k) = q(i-1,j,k)
               qlm2(i,j,k) = qlm1(i,j,k)
       enddo; enddo; enddo
     endif
     if (isxe==1) then
        do k=1,meqn
           do j=sy-mbc,ey+mbc
              do i=1,3
               qrm1(i,j,k) = q(nx-2+i,j,k)
               qrm2(i,j,k) = qrm1(i,j,k)         
        enddo; enddo; enddo
     endif  
     return            
     end subroutine radinit     

!     ==============================================================
      subroutine radup( a )
!     ==============================================================

    use m_parameters
    use m_variables,  only : ql,qr
    
    implicit none
    integer(kind=4) :: j, l, m
    real(fpp), intent(inout) :: a(sx-mbc:ex+mbc,sy-mbc:ey+mbc,1:meqn)
    
    if( isxs == 1 ) then
      do l = 1,mbc
         do m = 1,meqn
            do j = sy, ey
               a(1-l,j,m)  = ql(j,m)
            enddo
         enddo
      enddo
    endif
    
    if( isxe == 1 ) then
      do l = 1,mbc
         do m = 1,meqn
            do j = sy, ey
               a(nx+l,j,m) = qr(j,m)
            enddo
         enddo
      enddo
    endif
        
      return
      end subroutine radup


! ==================================================================
      subroutine out(time)
! ==================================================================
      use m_parameters, only: fpp
      use m_variables,  only: outfile,wrk1,wrk2,wrk3
      
      implicit none
      real(fpp), intent(in) :: time
      
      open(unit=11,file=outfile,access='append',form='unformatted')
        write(11) time
        write(11) wrk1
        write(11) wrk2
        write(11) wrk3
      close(unit=11)

      return
      end subroutine out
 
!   ********************************************************************
      subroutine rddump(time)
!   ********************************************************************
      use m_parameters, only: fpp,nx,ny
      use m_variables,  only: restart,wrk1,wrk2,wrk3
      
      implicit none
      real(fpp), intent(out) :: time
      integer(kind=4) :: n,k,nout,mx,my,meqn1,mbc1,ixbc,iybc
      real(fpp) :: tend,dt,dx,dy,gamma,cf,cb
      real(fpp) :: xx(1:nx),yy(1:ny)
      
      open(11,file=restart,status='old',form='unformatted')
       read(11) nout, tend
       read(11) mx, my, meqn1, mbc1
       read(11) dt, dx, dy, gamma, cf, cb
       read(11) ixbc, iybc
       read(11) xx
       read(11) yy
       read(11) wrk1
       read(11) wrk1
!   cycle through and read the last time from previous run
       do n = 1,nout+1 
         read(11,end=800) time
         read(11) wrk1(:,:)
         read(11) wrk2(:,:)
         read(11) wrk3(:,:)
       enddo
  800 close(unit=11)
  
      return
      end subroutine rddump
      
!     ==========================================**********************
      real(fpp) function philim(a,b,method)
!     ==========================================**********************
!     # Compute a limiter based on wave strengths a and b.
!     # method determines what limiter is used.
!     # a is assumed to be nonzero.
!       method = 1: minmod
!              = 2: superbee
!              = 3: van Leer
!              = 4: monotinized centered 

      use m_parameters
      
      implicit none
      
      integer(kind=4) :: method
      real(fpp), parameter :: eps2 = 0.0000000001_fpp
      real(fpp), intent(in) :: a,b
      real(fpp) :: r

      r = b/(a + sign(one, a)*eps2)
      
      select case (method)
        case (1)  
          philim = dmax1(zero, dmin1(one, r))
        case (2)  
          philim = dmax1(zero, dmin1(one, two*r), dmin1(two, r))
        case (3)  
          philim = (r + dabs(r)) / (one + dabs(r))
        case (4)  
          philim = dmax1(zero, dmin1((one + r)/two, two, two*r))
        case default
          philim = one
      end select
      
      return
      end function philim

     
#ifdef USEMPI
! ********************************************************************************  
      subroutine collect(a1,a2,myid,numprocs,comm2d,root,done) 
! ********************************************************************************  
! send/collect all parts of a distributed matrix to/on the root process
      use m_parameters
      use mpi
      
      implicit none
      
      integer(kind=4), intent(in) :: myid, numprocs
      integer(kind=4), intent(in) :: comm2d,root,done
      integer(kind=4) :: mx, my, i, j, k, ii
      integer(kind=4) :: iosx, iosy, numx, numy
      integer(kind=4) :: msgtype,source,ierr
      integer(kind=4) :: status(MPI_STATUS_SIZE)
      real(fpp), intent(in)  :: a1(sx:ex,sy:ey)
      real(fpp), intent(out) :: a2(nx,ny)
      
        if (myid .ne. 0) then    
!          call ret2root(a1(sx:ex,sy:ey),sx,ex,sy,ey, &
!                        comm2d,root,done) 
          mx = ex-sx+1
          my = ey-sy+1
          call MPI_SEND( sx, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
          call MPI_SEND( sy, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
          call MPI_SEND( mx, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
          call MPI_SEND( my, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
          do j = sy,ey
            call MPI_SEND( a1(sx,j), mx, MPI_DOUBLE_PRECISION, &
                        root, done, comm2d, ierr )
          enddo
        endif
        
        if (myid .eq. 0) then 
          a2(sx:ex,sy:ey) = a1(sx:ex,sy:ey)
          do ii=1, numprocs-1
!            source = ii
!            msgtype = done
            call MPI_RECV( iosx, 1, MPI_INTEGER, ii, &
                      done, comm2d, status, ierr )
            call MPI_RECV( iosy, 1, MPI_INTEGER, ii, &
                      done, comm2d, status, ierr )
            call MPI_RECV( numx, 1, MPI_INTEGER, ii, &
                      done, comm2d, status, ierr )
            call MPI_RECV( numy, 1, MPI_INTEGER, ii, &
                      done, comm2d, status, ierr )
            do j = 1,numy
              call MPI_RECV( a2(iosx,iosy+j-1),numx, &
                          MPI_DOUBLE_PRECISION,ii,done, & 
                          comm2d,status,ierr)
            enddo
          enddo
        endif	 

      end subroutine collect
#endif

#ifdef USEMPI      
! ********************************************************************************  
	subroutine ret2root( a,sx,ex,sy,ey,comm2d,root,done )
! ********************************************************************************  
!   send data from workers back to root processor
!	 
    use m_parameters, only: fpp
	use mpi	
    
	implicit none

	integer(kind=4), intent(in) :: sx,ex,sy,ey
	real(fpp),       intent(in) :: a(sx:ex,sy:ey) 
	integer(kind=4), intent(in) :: comm2d, root, done
	integer(kind=4) :: j, nx, ny
	integer(kind=4) :: status(MPI_STATUS_SIZE), ierr 
     
        nx = ex-sx+1
        ny = ey-sy+1
        call MPI_SEND( sx, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
        call MPI_SEND( sy, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
        call MPI_SEND( nx, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
        call MPI_SEND( ny, 1, MPI_INTEGER, root, done, &
                       comm2d, ierr )
        do j = sy,ey
          call MPI_SEND( a(sx,j), nx, MPI_DOUBLE_PRECISION, &
                        root, done, comm2d, ierr )
        enddo
 
 	return 
	end subroutine ret2root
#endif

#ifdef USEMPI      
! ********************************************************************************  
       subroutine exchng2( a, comm2d, stride, &
                           nbrleft, nbrright, nbrbottom, nbrtop  ) 
! ********************************************************************************  
       use mpi        
       use m_parameters, only: mbc, fpp, sx, ex, sy, ey
       
       implicit none
       real(fpp),    intent(inout) :: a(sx-mbc:ex+mbc, sy-mbc:ey+mbc) 
       integer(kind=4), intent(in) :: nbrleft, nbrright, nbrtop, nbrbottom
       integer(kind=4), intent(in) :: comm2d, stride
       integer(kind=4) :: status(MPI_STATUS_SIZE), ierr 
       integer(kind=4) :: i, j, numx

!  send/receive top/bottom first
       numx = mbc*(ex - sx + 1 + 2*mbc)
	   call MPI_SENDRECV( a(sx-mbc,ey-(mbc-1)), numx, MPI_DOUBLE_PRECISION, &
                          nbrtop, 0,  &
                          a(sx-mbc,sy-mbc), numx, MPI_DOUBLE_PRECISION,  &
                          nbrbottom, 0, comm2d, status, ierr ) 
	   call MPI_SENDRECV( a(sx-mbc,sy),  numx, MPI_DOUBLE_PRECISION, &
                          nbrbottom, 1,  &
                          a(sx-mbc,ey+1), numx, MPI_DOUBLE_PRECISION,  &
                          nbrtop, 1, comm2d, status, ierr )  
!  send/receive left/right and include the corner halo data
!  uses the vector datatype stridetype 
	   call MPI_SENDRECV( a(ex-(mbc-1),sy-mbc),   1, stride, nbrright, 0, &
                       a(sx-mbc,sy-mbc), 1, stride, nbrleft,  0, &
                       comm2d, status, ierr ) 
	   call MPI_SENDRECV( a(sx,sy-mbc),   1, stride, nbrleft,  1, &
                       a(ex+1,sy-mbc), 1, stride, nbrright, 1, &
                       comm2d, status, ierr ) 
       
	return 
	end subroutine exchng2
#endif

#ifdef USEMPI      
! ********************************************************************************  
       subroutine exchng2_hold( a, comm2d, stride, &
                           nbrleft, nbrright, nbrbottom, nbrtop  ) 
! ********************************************************************************  
       use mpi        
       use m_parameters, only: mbc, fpp, sx, ex, sy, ey
       
       implicit none
       real(fpp),    intent(inout) :: a(sx-mbc:ex+mbc, sy-mbc:ey+mbc) 
       integer(kind=4), intent(in) :: nbrleft, nbrright, nbrtop, nbrbottom
       integer(kind=4), intent(in) :: comm2d, stride
       integer(kind=4) :: status(MPI_STATUS_SIZE), ierr 
       integer(kind=4) :: i, j, numx

!  send/receive top/bottom first
       numx = ex - sx + 1  
       do j = 1, mbc
	   call MPI_SENDRECV( a(sx,ey-(j-1)),  numx, MPI_DOUBLE_PRECISION, &
                          nbrtop, 0,  &
                          a(sx,sy-j), numx, MPI_DOUBLE_PRECISION,  &
                          nbrbottom, 0, comm2d, status, ierr ) 
	   call MPI_SENDRECV( a(sx,sy+(j-1)),  numx, MPI_DOUBLE_PRECISION, &
                          nbrbottom, 1,  &
                          a(sx,ey+j), numx, MPI_DOUBLE_PRECISION,  &
                          nbrtop, 1, comm2d, status, ierr ) 
       enddo
 
!  send/receive left/right and include the corner halo data
!  uses the vector datatype stridetype 
       do i = 1,mbc
	call MPI_SENDRECV( a(ex,sy-mbc),   1, stride, nbrright, 0, &
                       a(sx-i,sy-mbc), 1, stride, nbrleft,  0, &
                       comm2d, status, ierr ) 
	call MPI_SENDRECV( a(sx,sy-mbc),   1, stride, nbrleft,  1, &
                       a(ex+i,sy-mbc), 1, stride, nbrright, 1, &
                       comm2d, status, ierr ) 
       enddo
       
	return 
	end subroutine exchng2_hold
#endif

#ifdef USEMPI      
! ********************************************************************************  
      subroutine get_nbrs2d( comm2d, nbrleft, nbrright, nbrbottom, nbrtop ) 
! ********************************************************************************  
! 
! This routine show how to determine the neighbors in a 2-d decomposition of 
! the domain. This assumes that MPI_Cart_create has already been called  
! 
      use mpi
      implicit none
      integer(kind=4), intent(in)  :: comm2d  
      integer(kind=4), intent(out) :: nbrleft, nbrright, nbrbottom, nbrtop  
      integer(kind=4) :: ierr 
 
      call MPI_Cart_shift( comm2d, 0,  1, nbrleft,   nbrright, ierr ) 
      call MPI_Cart_shift( comm2d, 1,  1, nbrbottom, nbrtop,   ierr ) 
 
      return 
      end subroutine get_nbrs2d
#endif

#ifdef USEMPI       
! ********************************************************************************  
      subroutine decomp2d( comm2d, nx, ny, sx, ex, sy, ey ) 
! ********************************************************************************  
      use mpi
      implicit none
      integer(kind=4), intent(in) :: comm2d, nx, ny  
      integer(kind=4), intent(out) :: sx, ex, sy, ey  

      integer(kind=4) :: dims(2), coords(2), ierr 
      logical :: periods(2) 
 
      call MPI_Cart_get( comm2d, 2, dims, periods, coords, ierr )  
      call decomp1d( nx, dims(1), coords(1), sx, ex ) 
      call decomp1d( ny, dims(2), coords(2), sy, ey ) 
 
      return 
      end subroutine decomp2d
#endif

#ifdef USEMPI             
! ********************************************************************************  
      subroutine decomp1d( n, numprocs, myid, s, e ) 
! ************************************************************************************  
! 
!  This file contains a routine for producing a decomposition of a 1-d array 
!  when given a number of processors.  It may be used in "direct" product 
!  decomposition.  The values returned assume a "global" domain in [1:n] 
!
      implicit none 
      integer(kind=4), intent(in)  :: n, numprocs, myid
      integer(kind=4), intent(out) :: s, e 
      integer(kind=4) :: nlocal, deficit 
 
      nlocal  = n / numprocs 
      s	      = myid * nlocal + 1 
      deficit = mod(n,numprocs) 
      s	      = s + min(myid,deficit) 
      if (myid .lt. deficit) then 
          nlocal = nlocal + 1 
      endif 
      e = s + nlocal - 1 
      if (e.gt.n .or. myid.eq.numprocs-1) e = n 
      return 
      end subroutine decomp1d
#endif





