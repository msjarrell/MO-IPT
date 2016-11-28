! ***************************************************************************
! ***************************************************************************
       module Global
       implicit none
       save
       integer Ndeg,LDA_DMFT,phsym,Norbs,Nm
       integer NQPTS,Nkp,N_PTS,Ndim
       integer, parameter :: Nmdos=1000
       !integer, parameter :: NQPTS=15625
       !integer, parameter :: Nkp=160
       real*8, allocatable :: nf(:),n0(:)                       
       real*8::sum1,sum2,delta1,dfac
       real*8, allocatable :: U_ab(:,:)
       complex*16, allocatable:: Ham0(:,:,:)
       complex*16, allocatable:: Ham_model(:,:,:)      
       real*8:: LATT_A,mu0,correc,ntot_input,correc2
       integer :: N,Ndos,idos,init,nloop,Jflag
       real*8:: U,t,J_H,eta,temp,beta,pi,frac                     
       real*8, allocatable :: nf_tmp(:),Zfac(:),slope(:)
       real*8, allocatable ::w(:),dw(:),ferm(:)            
       real*8, allocatable :: chi1(:,:),chi2(:,:)
       real*8, allocatable::ep_f(:),OrbE(:)          
       real*8:: mu_c,nf_new,ntot,mu_c_guess,alpha,shift,LDA_mu
       complex*16::ii,zero
       complex*16, allocatable::Gfscript(:,:),sigma(:,:),Gf(:,:)
       complex*16,allocatable::sigma_dyn(:,:)
       real*8:: ep1,ep2,ep3,ep4,ep5,b1,b2,db2,b3,db3,b4,db4,uj
       real*8, allocatable:: rr1(:),rr2(:),asym_p(:),epf_ref(:)       
       complex*16:: z1,z2,z3,z4,z5,z6
       integer  :: num_tasks,ierr,task_id
       integer :: jstart,jend,uniform_grid 


      end  module Global

! ***************************************************************************
! ***************************************************************************

      Program MULTI_ORBITAL_HUBBARD
!     Frequency grid from makegrid.f
!     Auxiliary routines are in funct.f
!     ep_f is DFT input

      use Global
      implicit none
      include 'mpif.h'
      real*8 r1,r2,r3,r4,r5,r6,xx,yy,mu0_tmp
      integer flag,io,jo,ie,acount,infot
      character*72 linemc
      integer i,j,k,info
      complex*16 gauss,z
      real*8, allocatable :: wl(:),muin(:)
      real*8 :: epf0,fermic,lutt_ei,lval,Lutt_val
      real t1,t2
      logical check
      external gauss,fermic,funcv2,Lutt_val,funcv3
      include'Glob_cons'
       
 1000	format(a72)


        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,num_tasks,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,task_id,ierr)

        if(task_id.eq.0) then
          write(6,*)
          write(6,*)'task_id=',task_id
!          write(6,*)'****************************************************************'
!          call system("banner MO-IPT")
!          write(6,*)'****************************************************************'
          write(6,*)
          write(6,*)
        end if


!	READ PAR.DAT INPUT FILE
 	open(unit=12,file='par.dat',status='unknown')
	read(12,*) b1,ep1,ep2,b2,db2,ep3! Grid parameters : see makegrid.f
	read(12,*) b3,db3,ep4,ep5	! Grid parameters : see makegrid.f
	read(12,*) U,J_H,dfac,correc,correc2,t,alpha,Jflag
        read(12,*) phsym,LDA_DMFT,init,idos,uniform_grid
	read(12,*) temp,frac,eta
	read(12,*) Norbs,Nkp,NQPTS,Ndim
	close(12)
        if(Norbs==0) then
           write(6,*) 'Number of orbitals must be greater than 0'
           stop
        end if
        if(Nkp==0) Nkp=1
        if(NQPTS==0) NQPTS=1

        if(Ndim==2) then
           N_PTS=((2*Nkp)+1)**2
        else
           !Default is 3-D
           N_PTS=((2*Nkp)+1)**3
        end if

        allocate(nf(Norbs),stat=info)
        allocate(n0(Norbs),stat=info)
        allocate(nf_tmp(Norbs),stat=info)
        allocate(Zfac(Norbs),stat=info)
        allocate(slope(Norbs),stat=info)
        allocate(ep_f(Norbs),stat=info)
        allocate(OrbE(Norbs),stat=info)
        allocate(rr1(Norbs),stat=info)
        allocate(rr2(Norbs),stat=info)
        allocate(asym_p(Norbs),stat=info)
        allocate(epf_ref(Norbs),stat=info)
        allocate(U_ab(Norbs,Norbs),stat=info)
        allocate(Ham0(NQPTS,Norbs,Norbs),stat=info)
        allocate(Ham_model(N_PTS,Norbs,Norbs),stat=info)
        allocate(muin(Norbs+1),stat=info)

  	open(unit=13,file='initials.dat',status='unknown')
        do io=1,Norbs 
          read(13,*) nf(io)
        end do
       	read(13,*) mu0,mu_c_guess,shift,ntot_input    		
 	close(13)
        open(unit=13,file='orb_energy.dat',status='unknown')
         do io=1,Norbs
        read(13,*) ep_f(io)
        end do
        close(13)


        call fill_U_matrix()


        open(unit=99,file='interaction_matrix.dat')      
        do io=1,Norbs
           do jo=1,Norbs
              write(99,*)U_ab(io,jo)
           end do
        end do
        close(99)
            
!***********************************************
! Solve for alpha - Mott pole weight not implemented in
!                   this release
           
           
!*******************************************************

!	THERE IS A CONSTRAINT IN THE ASYMMETRIC CASE: THE LUTTINGER
!	INTEGRAL. THE VARIABLE EI IS USED TO SATISFY THIS CONSTRAINT

!	DESCRIPTION OF VARIOUS VARIABLES
!	nf  : total correlated orbital occupation number.
!	ei  :Used to satisfy Luttinger's' theorem=-mu0 of KK paper(PRL V.77pp131)
!	ep_f =-U/2 for the particle-hole symmetric case
!       init = 0   Start a fresh calculation; no previous data available
!	     = 1   Supply previously computed data as input
!	     = 2   Supply previously computed data and stop after computing Transport
!	idos = 1   Hypercubic lattice Gaussian Density of states: see gauss(z) in funct.f
!	     = 2   Bethe lattice or Cayley tree: Semicircular density of states: see 
!		   sem(z) in funct.f	
!	     = 3   General density of states: see below and the subroutine gendos(z)
!	     = 4   LDA Tight - Binding input 
!	temp  : Temperature
!	


	
		
	if(temp.eq.0.d0) then
	  beta=0.d0
	else
	  beta=1.d0/temp
	end if

!	CREATE THE FREQUENCY GRID
	call makegrid		    ! Defined in the makegrid.f file
        
        acount=0
        infot=0
        allocate(wl(-Nm:Nm),stat=info)
        infot=infot+info
        allocate(chi1(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(chi2(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(Gfscript(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(sigma(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(Gf(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(sigma_dyn(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        if(infot.ne.0) write(6,*) 'Allocation failed in',task_id


!********************************************        
!**************************************************
         
	if(idos.gt.3) then
	   !READ INPUT FROM LDA-LMTO OUTPUT (HOPPING PARAMETERS AND
	   !ORBITAL ENERGIES
           call lattice_input()
        end if
          
           mu_c=mu_c_guess   
                 
            ntot=0.d0
            do io=1,Norbs
              ntot=ntot+nf(io)
            end do
          if(ntot.ne.ntot_input) ntot=ntot_input
             
              

 	  if(init.eq.0) then	! BEGIN WITH JUST THE Self-consistent
                                ! HARTREE Calculation

          !OR Get the occupancies from the LDA calculation and fix them
             sigma=0.d0
            
!       Guess self-energy
          do j=1,Norbs
           sum1=0.d0
           do k=1,Norbs
            if(k.ne.j) then
             sum1=sum1+U_ab(j,k)*nf(k)               
            end if   
           end do
           do i=-N,N
           sigma(j,i)=sum1 !+U**2*alpha**3/(4.d0*(w(i)+ii*eta))
           end do
         end do
          if((idos==1).or.(idos==2))then
           call findgf()
          else
           call ksummation()
          end if

          
       else				!BEGIN WITH OUTPUT OF A PREVIOUS RUN
	if(task_id==0) write(*,*)'start from old data'
        	 
        call readinput()
        
               
 
         if((idos==1).or.(idos==2))then
          call findgf()
         else
          call ksummation() 
         end if

       end if
        
           
!	INITIALIZATION OVER. NOW THE WORKS - GREEN FUNCTIONS
!	AND SELF ENERGIES AND DMFT.
	
!        we are not using outer loop(mu0). we change mu0 manually.

        call dmft_selfconsistency()
        lval=Lutt_val()

        if(task_id==0)write(6,*) 'mu0=',mu0,' LVAL=',lval,'mu_c=',mu_c 
        !if(task_id==0) then 
        !do io=1,Norbs
        ! write(6,*)io, ep_f(io)
        !end do    
        !endif  

!******************Quasi-particle-weights******************                  
        slope=0.0d0
        zfac=0.0d0
        do j=1,Norbs
          xx=0.0d0
          yy=0.0d0
          r1=0.0d0
          r2=0.0d0
          r3=0.0d0
          r4=0.0d0
          do i=-100,100
            xx=w(i)
            yy=dreal(sigma(j,i))
            r1=r1+xx
            r2=r2+xx**2
            r3=r3+yy
            r4=r4+xx*yy
          end do
          slope(j)=(r1*r2-201.d0*r4)/((r1**2)-201.d0*r2)
          zfac(j)=1.d0/(1.d0-slope(j))
        end do
        open(unit=73,file='quasi_wei.dat',status='unknown') 
        do j=1,Norbs
          write(73,*) j,zfac(j)
        end do
        close(73)
 
!******************************************************        

!       WRITE INTO FILES
        if(task_id==0) call writeoutput
        if(Norbs.le.6.and.task_id==0) call realtoImg()
        
!******************************************************        
        
        deallocate(nf)
        deallocate(n0)
        deallocate(nf_tmp)
        deallocate(Zfac)
        deallocate(slope)
        deallocate(ep_f)
        deallocate(OrbE)
        deallocate(rr1)
        deallocate(rr2)
        deallocate(asym_p)
        deallocate(epf_ref)
        deallocate(U_ab)
        deallocate(Ham0)
        deallocate(Ham_model)
        deallocate(chi1)
        deallocate(chi2)
        deallocate(Gfscript)
        deallocate(sigma)
        deallocate(Gf)
        deallocate(sigma_dyn)
        deallocate(wl)
        deallocate(muin)

        call MPI_FINALIZE(ierr)
	end
!	************************************************************
       
	  subroutine hartree()

	  use Global
	  implicit none
	  logical check
	  integer io,i,iter
	  real*8 muin(1),tmp_muc,broyden
	  external funcv1,broyden
	  include 'Glob_cons'


          if(task_id==0) write(6,*) 'Going into Hartree'

	  tmp_muc=mu_c
          !call broydn(muin,1,check,funcv1)
	  !mu_c=muin(1)
          mu_c=broyden(tmp_muc)
          call occupancies()

          if(task_id==0) write(6,*) 'Hartree done'
	  stop
          return
	  end       
!**************************************************************
	subroutine funcv1(nn,xx,fvec)
        use Global
	integer nn,io
        real*8 xx(nn),fvec(nn),nf_tot_tmp
    !    external Lutt_val
        include 'Glob_cons'
         mu_c=xx(1)
! 	call findgf()
        if((idos==1).or.(idos==2))then
         call findgf()
        else	  ! New Gf
         call ksummation() 
        end if

	 call occupancies()
	 nf_tot_tmp=0.d0
!         write(6,*) ntot
	 do io=1,Norbs
	  nf_tot_tmp=nf_tot_tmp+nf(io)
	 end do
	 fvec(1)=ntot-nf_tot_tmp
        if(task_id==0) write(*,*)'funcv1 mu_c=',mu_c,' fvec',fvec(1)
        

	return
	end

!********************************************
        subroutine funcv2(nn,xx,fvec)
        use Global
	integer nn
	real*8 xx(nn),fvec(nn),x1
        real*8 Lutt_val
	external Lutt_val
        include 'Glob_cons'

        mu0=xx(1)
!	x1=xx(1)
	 call dmft_selfconsistency()	
        fvec(1)=Lutt_val()

        if(task_id==0) then

        write(*,*)'--------------------------------------'
        write(*,*)'OUTER LOOP'
        write(*,*)'funcv2 mu0=',mu0,' fvec',fvec(1)
        write(*,*)'--------------------------------------'
        end if
	return
	end
 !***********************************************        
        subroutine funcv3(nn,xx,fvec)
        use Global
        integer nn
        real*8 xx(nn),fvec(nn),x1
        include 'Glob_cons'
        alpha=xx(1)
        fvec(1)=alpha**3-alpha**2+t**2/U**2
        if(task_id==0) write(6,*)'alpha=',alpha,'fvec=',fvec(1)
        return
        end

!*************************************************
	real*8 function Lutt_val()
        use Global
        implicit none
	real*8 r1,r2,r3,r4,r5,r6
        integer ie,jo
        real*8 Lutt_n(50)
	complex*16 dsig
!	IMPLEMENTATION OF EQ.(3.15) OF THESIS
	include 'Glob_cons'
!	Now evaluate the Luttinger integral eq.(3.15)

        if(task_id==0) then
        if(Norbs.gt.50) then
          write(6,*) 'Modify the array size of Lutt_n to Norbs+1'
          write(6,*) 'And recompile'
          stop
        end if
        end if


	do jo=1,Norbs
          r1=0.d0
          do ie=-N+1,0
	    dsig=0.5d0*((sigma(jo,ie+1)-sigma(jo,ie))/(w(ie+1)-w(ie))+&
     		(sigma(jo,ie+2)-sigma(jo,ie+1))/(w(ie+2)-w(ie+1)))
            r1=r1+dimag(dsig*Gf(jo,ie))*dw(ie)
	  end do
          Lutt_n(jo)=r1
        enddo 
	if(task_id==0) write(6,*) 'Luttinger values- Orb-1', Lutt_n(1),'Orb-3',Lutt_n(3)

        r1=0.d0 
        do jo=1,Norbs
          r1=r1+Lutt_n(jo)
        enddo   
        Lutt_val=r1/2.d0        ! Division by 2 in the paramagnetic case

	return
	end
!	************************************************************
	
        subroutine dmft_selfconsistency()
        use Global
        implicit none
	real*8 r1,r2,r3,r4,r5,r6
	integer nlmin,nlmax,info,io
        integer i,j
	real*8 rho,conv,tolf,lval,Lutt_val,mu_0tmp,nf_tot_tmp
        real*8 ::  muin(1),t5,t4,tmp_muc,broyden
        real*8,allocatable ::prho(:,:)
	external funcv1,Lutt_val,broyden
	logical check
	include 'Glob_cons'

        allocate(prho(Norbs,-Nm:Nm),stat=info)

 	nloop=1
	nlmin=2
	nlmax=3000
	tolf=2.D-4	! Tolerance for convergence; may be decreased if
			! greater accuracy is required
        conv=1.0 
! 	BEGIN DMFT SELF CONSISTENCY LOOPS
	do while(nloop.le.nlmin.or.((conv.gt.tolf).and. &
     					(nloop.le.nlmax))) 


	  do j=1,Norbs                                    !Orbitals
	     do i=-N,N
	        prho(j,i)=-dimag(Gf(j,i))/pi
	     enddo
	  enddo	                                      !Orbitals
 
          call evalchi		! Polarization propagator


	  if(task_id==0) write(6,*) 'Found Chi'
          call cpu_time(t4)
          call selfenergy	    ! Self energy
          call cpu_time(t5)
          if(task_id==0) write(6,*)'time2',t5-t4				
	  if(task_id==0)  write(6,*) 'Found Sigma'
            
          if((idos==1).or.(idos==2))then
            call findgf()
          !do i=-N,N
          !write(301,*) w(i),dreal(sigma(1,i)),-dimag(Gf(1,i))
          !end do
          !close(301)
          !stop
           else
            call ksummation()   
           end if     
         
        call occupancies()          
            


        if(phsym==0)then
                
           
	    
          nf_tot_tmp=0.d0
 	  do io=1,Norbs
	    nf_tot_tmp=nf_tot_tmp+nf(io)
	  end do
	  if(task_id==0)  write(6,*) 'Occupancy deviation=',ntot-nf_tot_tmp
	  if(dabs(ntot-nf_tot_tmp).gt.1.D-4) then
	    !muin(1)=mu_c  				
            !call broydn1(muin,1,check,funcv1)	! non-linear equation solver 
            !mu_c=muin(1)  					! 
	    tmp_muc=mu_c
            mu_c= broyden(tmp_muc)	! non-linear equation solver 
	  end if
          if(task_id==0) write(6,*)'mu_c',mu_c           
           if((idos==1).or.(idos==2))then
            call findgf()
           else
            call ksummation() 
           end if

          end if

          r1=Lutt_val()
          lval=r1


	  conv=0.0
	  r1=0.d0
	  r2=0.d0
          do j=1,Norbs			                      !Orbitals
	     do i=-N,N
	       rho=-dimag(Gf(j,i))/pi
	       r1=r1+dabs(rho-prho(j,i))*dw(i)
	       r2=r2+rho*dw(i)
	     end do
           enddo                                    ! Orbitals
	  conv=r1/r2
         
     
          if(task_id==0) write(6,*)&
          '**************************************************'
          if(task_id==0) write(6,206) nloop,conv,lval
          if(task_id==0) write(6,*)&
          '**************************************************'
	
          nloop=nloop+1
                

 
	end do			! END NLOOP


 206	format('nloop=',i4,'  conv=',f15.6,'  Lval=',f15.6)

        deallocate(prho)

 	return
	end
!	************************************************************
   	subroutine occupancies()

        use Global
 	real*8 r1,r2,r3,r4,r5,r6
        integer i,j,io,jo,ie,flag_occ
	include 'Glob_cons'


         do jo=1,Norbs      ! Orbitals
          r2=0.d0
          r1=0.d0
	  do ie=-N,N
            r2=r2+dw(ie)*(-dimag(Gf(jo,ie))/pi)*ferm(ie)
            r1=r1+dw(ie)*(-dimag(Gf(jo,ie))/pi)
	  end do
           nf(jo)=(r2/r1)
           
	   
            if(dabs(r1-1.d0).gt.0.02d0) then
             if(task_id==0) write(6,*) 'WARNING: Gf not normalized-',jo,r1
           end if	     	
         enddo                            ! Orbitals

         do jo=1,Norbs
           r1=0.d0
           r2=0.d0
           do ie=-N,N
	     r1=r1+(-dimag(Gfscript(jo,ie))/pi)*ferm(ie)*dw(ie)   
	     r2=r2+(-dimag(Gfscript(jo,ie))/pi)*dw(ie)
           end do
          
	  n0(jo)=(r1/r2)
          
          
	  if(dabs(r2-1.d0).gt.0.02d0) then
          if(task_id==0) write(6,*) 'WARNING: Gfscript not normalized-',jo,r2
          end if	     	
         end do

          if(phsym==1)then
            do jo=1,Norbs
              nf(jo)=0.5d0
              n0(jo)=0.5d0
             end do
           end if
          
	 if(task_id==0)  then
            if(Norbs==2) then
             write(6,*) 'nf=',nf(1),nf(2)
	     write(6,*) 'n0=',n0(1),n0(2)
             else if(Norbs==4) then
             write(6,*) 'nf=',nf(1),nf(2),nf(3),nf(4)
	     write(6,*) 'n0=',n0(1),n0(2),n0(3),n0(4)
             else if(Norbs.gt.4) then
             write(6,*) 'nf=',nf(1),nf(2),nf(3),nf(4), &
                              nf(5),nf(6)
	     write(6,*) 'n0=',n0(1),n0(2),n0(3),n0(4), &
                              n0(5),n0(6)
             end if
        end if
         
	 return
	 end    ! subroutine occupancies


!	************************************************************
   	subroutine evalchi
    	use Global
	real*8 r1,r2,r3,r4,r5,r6
        integer i,j,info
	real*8,allocatable :: rho1(:),rho2(:),chi(:)
        real*8 :: sgn
	external sgn
	include 'Glob_cons'

        allocate(rho1(-Nm:Nm),stat=info)
        if(info.ne.0) write(6,*) 'Allocation of rho1 failed'
        allocate(rho2(-Nm:Nm),stat=info)
        if(info.ne.0) write(6,*) 'Allocation of rho2 failed'
        allocate(chi(-Nm:Nm),stat=info)
        if(info.ne.0) write(6,*) 'Allocation of chi failed'



        do j=1,Norbs
	  do i=-N,N
	    rho1(i)=-ferm(i)*dimag(Gfscript(j,i))/pi
	    rho2(i)=-ferm(-i)*dimag(Gfscript(j,i))/pi
	  end do
          call convolve(Nm,rho1,rho2,chi)
          do i=-N,N
            chi1(j,i)=chi(i)
          end do
        enddo
!         do j=1,3,2
!           do i=-N,N
!        chi1(j,i)=chi1(j,i)!+alpha*0.5d0*(rho2(i)+rho1(-i))
!           end do
!         end do  

        !if(task_id==0) write(6,*) 'Chi1 done'
        do j=1,Norbs
	  do i=-N,N
	    rho1(i)=-ferm(-i)*dimag(Gfscript(j,i))/pi
	    rho2(i)=-ferm(i)*dimag(Gfscript(j,i))/pi
	  end do
          call convolve(Nm,rho1,rho2,chi)
        !if(task_id==0) write(6,*) 'Chi2 done,orb=',j
          do i=-N,N
            chi2(j,i)=chi(i)
          end do
        enddo

 
        !if(task_id==0) write(6,*) 'Chi2 done'
        deallocate(rho1)
        deallocate(rho2)
        deallocate(chi)

  100	return
	end


!	************************************************************
  	subroutine selfenergy

  	use Global
        implicit none
	real*8 :: r1,r2,r3,r4,r5,r6
        integer :: i,j,ko,ie,io,jo,lo,info,infot
        real*8,allocatable :: rho1(:),rho2(:),resig(:),&
     	rhosigma1(:),rhosigma2(:),rhosigma(:)
        real*8, allocatable:: muin(:) 
        logical :: check
	real*8,allocatable:: Afac(:),Bfac(:),Num1(:),denm(:)
        real*8, allocatable:: naa(:),Num2(:),Num3(:)
        real*8,allocatable :: Num4(:),Two_cr(:)
        complex*16 :: nsigma
        real*8,allocatable :: D_n(:,:),B_num1(:),Hart(:)
	real*8, allocatable ::  B_num2(:),Afac2(:)
	complex*16,allocatable :: sigmanew(:,:),sigma2(:,:,:),&
              reg_sigma(:,:),num(:,:)
        external funcv2
        real*8 :: sum4,sum5,sum6,alphanf,sum11,sum22
        real*8, allocatable :: Num5(:),Num6(:)
        real*8, allocatable :: Num7(:),Num8(:),Num9(:),FinalB(:)
        real*8 :: sum33,sum44,sum55,sum66,sum666,sum3
   
        infot=0
        allocate(rho1(-Nm:Nm),stat=info)
        infot=infot+info
        allocate(rho2(-Nm:Nm),stat=info)
        infot=infot+info
        allocate(resig(-Nm:Nm),stat=info)
        infot=infot+info
        allocate(rhosigma1(-Nm:Nm),stat=info)
        infot=infot+info
        allocate(rhosigma2(-Nm:Nm),stat=info)
        infot=infot+info
        allocate(rhosigma(-Nm:Nm),stat=info)
        infot=infot+info
        allocate(sigma2(Norbs,Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(reg_sigma(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(num(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(sigmanew(Norbs,-Nm:Nm),stat=info)
        infot=infot+info
        allocate(muin(Norbs),stat=info)
        infot=infot+info
        allocate(Afac(Norbs),stat=info)
        infot=infot+info
        allocate(Bfac(Norbs),stat=info)
        infot=infot+info
        allocate(Num1(Norbs),stat=info)
        infot=infot+info
        allocate(Denm(Norbs),stat=info)
        infot=infot+info
        allocate(naa(Norbs),stat=info)
        infot=infot+info
        allocate(Num2(Norbs),stat=info)
        infot=infot+info
        allocate(Num3(Norbs),stat=info)
        infot=infot+info
        allocate(Num4(Norbs),stat=info)
        infot=infot+info
        allocate(Two_cr(Norbs),stat=info)
        infot=infot+info
        allocate(D_n(Norbs,Norbs),stat=info)
        infot=infot+info
        allocate(B_num1(Norbs),stat=info)
        infot=infot+info
        allocate(B_num2(Norbs),stat=info)
        infot=infot+info
        allocate(Afac2(Norbs),stat=info)
        infot=infot+info
        allocate(Hart(Norbs),stat=info)
        infot=infot+info
        allocate(Num5(Norbs),stat=info)
        infot=infot+info
        allocate(Num6(Norbs),stat=info)
        infot=infot+info
        allocate(Num7(Norbs),stat=info)
        infot=infot+info
        allocate(Num8(Norbs),stat=info)
        infot=infot+info
        allocate(Num9(Norbs),stat=info)
        infot=infot+info
        allocate(FinalB(Norbs),stat=info)
        infot=infot+info

        if(infot.gt.0) then 
           write(6,*) 'Allocation failed in selfenergy'
           stop
        !else
        !   write(6,*) 'Allocation successful in selfenergy'
        end if


!        calculate sigma2 matrix
!        dfac=0.5d0
!        alphanf=alpha*dfac
        sigma2=zero
!	do io=1,1
!          do jo=2,2
            do io=1,Norbs    !1,3,2
              do jo=1,Norbs  !2,5,3
              if(jo.ne.io) then
            do i=-N,N
	      rho1(i)=(-dimag(Gfscript(io,-i))/pi)*ferm(i)
	      rho2(i)=chi1(jo,i)
	    end do

  	    call convolve(Nm,rho1,rho2,rhosigma1)
             
!            do i=-N,N
!           rhosigma1(i)=rhosigma1(i)!+alphanf*(alphanf*rho1(-i)+rho2(i))
!            end do
	    
             do i=-N,N
	     rho1(i)=(-dimag(Gfscript(io,-i))/pi)*ferm(-i)
	     rho2(i)=chi2(jo,i)
	    end do

	    call convolve(Nm,rho1,rho2,rhosigma2)
            
!            do i=-N,N
!           rhosigma2(i)=rhosigma2(i)!+alphanf*(alphanf*rho1(-i)+rho2(i))
!            end do
            	    

             do i=-N,N
	     rhosigma(i)=rhosigma1(i) + rhosigma2(i)
	    end do
            !write(6,*) 'CONVOLUTION DONE, NOW KKTRANSF'


            call kktransf(rhosigma,resig)
            !do i=-N,N
            !write(107,*) w(i),dw(i)
            !end do
            !close(107)
            !stop


            do i=-N,N
              sigma2(io,jo,i)=((U_ab(io,jo))**2)* &
                       (resig(i)-pi*ii*rhosigma(i))   
            end do
           end if
          end do ! End orbital index jo
        end do   ! End orbital index io
        

!       Now calculate the D_n, A and B factors
       
	call occupancies()

!*********************************************
!  calculation of A_{alpha} factor:
!**********************************************
!  calculate first term in the numerator

         Num1=0.d0
        do io=1,Norbs
            do jo=1,Norbs
              if(jo.ne.io) then
          Num1(io)=Num1(io)+((U_ab(io,jo)**2)*(nf(jo)*(1.d0-nf(jo))))
              end if
             end do
        end do

!**********************************************
! calculate second term in the numerator
          Two_cr=0.d0                   
           do io=1,Norbs
              do i=-N,N
       Two_cr(io)=Two_cr(io)+(ferm(i)*dw(i)*dimag(sigma(io,i)*Gf(io,i)))
              end do
         Two_cr(io)=-Two_cr(io)/(pi*(Norbs-1))
           end do
         do io=1,Norbs
!         write(120,*) Two_cr(io)
         end do   
!****************************************************
! Calculate 2nd term in the A factor(see eq 1 in the notes)
! \sum_{\beta .neq. \alpha}U_{\alpha \beta}\sum_{\gamma .neq. \beta .neq. \alpha!}U_{\alpha \gamma}<n_{\beta}n_{\gamma}>-<n_{beta}><n_{gamma}>!
!*************************************************************
!              Num2=0.d0
              Num3=0.d0
!              Num4=0.d0
              do io=1,Norbs
               do jo=1,Norbs
                  do ko=1,Norbs
          if((ko.ne.jo).and.(ko.ne.io).and.(io.ne.jo))then
           Num3(io)=Num3(io)+U_ab(io,jo)*U_ab(io,ko)*nf(ko)*nf(jo)
              end if
            end do
           end do
           end do
!***********************************************************************
             Num2=0.d0
             Num4=0.d0
              do io=1,Norbs
               do jo=1,Norbs
                  do ko=1,Norbs
          if((ko.ne.jo).and.(ko.ne.io).and.(io.ne.jo).and. &
        U_ab(jo,ko).ne.0.d0) then
           Num2(io)=Num2(io)+(U_ab(io,jo)*(U_ab(io,ko)/U_ab(jo,ko))* &
        Two_cr(jo))
              end if
            end do
           end do
           Num4(io)=Num2(io)-Num3(io)
           end do
!*****************************************************************************

!*****************************************************
!****************************************************
!  Calculate third term in the numerator
!  \sum_{\beta .neq. \alpha} U^3_{\alpha \beta}n_{\beta}(1-n_{\beta})
!*****************************************************
              B_num1=0.d0
             do io=1,Norbs
                do jo=1,Norbs
              if(jo.ne.io) then
          B_num1(io)=B_num1(io)+(U_ab(io,jo)**3)*(1.d0-nf(jo))*nf(jo)
              end if
              end do
             end do
!*******************************************************
!Calculate 3 and 4 th numerators terms in the Bterm.
! \sum_{\beta .neq. \alpha}U_{\alpha \beta} \sum_{\gamma .neq. \beta .neq.
! \alpha}U^2_{\alpha \gamma}(<n_{\beta}n_{\gamma}>-<n_{\beta}><n_{\gamma}>)+
! \sum_{\beta .neq. \alpha}U^2_{\alpha \beta} \sum_{\gamma .neq. \beta .neq.
! \alpha}U_{\alpha \gamma}(<n_{\beta}n_{\gamma}>(1-<n_{\beta}>)
!*********************************************
              Num5=0.d0
              Num7=0.d0
              do io=1,Norbs
                do jo=1,Norbs
                  do ko=1,Norbs
           if((ko.ne.jo).and.(ko.ne.io).and.(jo.ne.io).and. &
         (U_ab(jo,ko).ne.0.d0)) then
            Num5(io)=Num5(io)+(U_ab(io,jo)*(U_ab(io,ko)**2/ &
            U_ab(jo,ko))*Two_cr(jo))
            Num7(io)=Num7(io)+2*(((U_ab(io,jo)**2)* &
         (U_ab(io,ko)/U_ab(jo,ko))*(1.d0-nf(jo))*Two_cr(jo)))
                 end if
               end do
              end do
           end do
!************************************************************************
              Num6=0.d0
              do io=1,Norbs
                do jo=1,Norbs
                  do ko=1,Norbs
           if((ko.ne.jo).and.(ko.ne.io).and.(jo.ne.io))then
            Num6(io)=Num6(io)+U_ab(io,jo)*(U_ab(io,ko)**2)*nf(ko)*nf(jo)
             end if
               end do
              end do
              Num6(io)=-Num6(io)
            end do

!***********************************************************
! Calculate the 4 th numerator in the B term.
! \sum_{\alpha .neq. \alpha} U_{\alpha \beta} \sum_{\gamma .neq. \beta .neq.
! \alpha}U_{\alpha \gamma} \sum_{\eta .neq. \gamma .neq. \beta .neq. \alpha}
! U_{\alpha \eta}n_{\beta}<n_{\gamma}n_{\eta}>
! For the time being neglect the 3 particle correlation function.
!*************************************************************** 
             Num8=0.d0   
           do io=1,Norbs
               do jo=1,Norbs
                if(jo.ne.io) then
                  do ko=1,Norbs
            if((ko.ne.jo).and.(ko.ne.io)) then
                 do lo=1,Norbs
         if((lo.ne.ko).and.(lo.ne.jo).and.(lo.ne.io).and. &
         (U_ab(ko,lo).ne.0.d0)) then
         Num8(io)=Num8(io)+(U_ab(io,jo)*U_ab(io,ko)* &
          (U_ab(io,lo)/U_ab(ko,lo))*Two_cr(ko)*nf(jo))
            end if
            end do
            end if
            end do
            end if
            end do
            Num8(io)=-Num8(io)
            end do
!**********************************************
! Calculate 3 particle correlation function:
           Num9=0.d0
         do io=1,Norbs
           do jo=1,Norbs
            if(jo.ne.io) then
             do ko=1,Norbs
              if((ko.ne.jo).and.(ko.ne.io)) then
                do lo=1,Norbs
                if((lo.ne.ko).and.(lo.ne.jo).and.(lo.ne.io).and. &
        (U_ab(jo,ko).ne.0.d0).and.(U_ab(jo,ko).ne.0.d0).and. &
        (U_ab(jo,ko).ne.0.d0)) then 
         Num9(io)=Num9(io)+U_ab(io,jo)*U_ab(io,ko)*U_ab(io,lo)* &
    ((Two_cr(jo)*nf(lo)/U_ab(jo,ko))+(Two_cr(lo)*nf(ko)/U_ab(lo,jo))+ &
      (Two_cr(ko)*nf(jo)/U_ab(ko,lo))-(2.0d0*nf(jo)*nf(ko)*nf(lo)))
                           end if
                         end do
                        end if
                      end do
                    end if
                  end do
               end do



!*****************************************************
!  calculate the denominator in A factor  
         
               denm=0.d0
                           
               do io=1,Norbs
                 do jo=1,Norbs
                   if(jo.ne.io) then
          denm(io)=denm(io)+((U_ab(io,jo)**2)*(n0(jo)*(1.d0-n0(jo))))
                end if
               end do
            end do

!****************************************************
! calculate A factor
              Afac2= 0.d0
              Afac=0.d0
            do io=1,Norbs
                
             Afac2(io)= Num1(io)+correc*Num4(io)
             Afac(io)=Afac2(io)/denm(io) 
              if(Afac(io).lt.0.d0) then
                if(task_id==0) write(6,*) 'A is negative'
            end if
             end do
             
!**************************************************
!*************************************************
! calculation of B factor 
!************************************************
!************************************************
! calculate Hartree term 
              Hart=0.d0
              do io=1,Norbs
                do jo=1,Norbs
              if(jo.ne.io) then
           Hart(io)=Hart(io)+(U_ab(io,jo)*nf(jo))
              end if
             end do
           end do
!***************************************************
!  calculate first and third term in the numerator of B 
             B_num2=0.d0
             do io=1,Norbs
          B_num2(io)=(mu0-mu_c+ep_f(io)-Hart(io))/denm(io)             
            end do
!****************************************************
! Sum the total terms in the B term.
        FinalB=0.d0
        do io=1,Norbs
        FinalB(io)=B_num1(io)+Num5(io)+Num6(io)+Num7(io)+Num8(io)+&
       Num9(io)*correc2
        FinalB(io)=FinalB(io)/(Afac2(io)*denm(io))
        end do


!***************************************************
! calculate B factor
              Bfac=0.d0
              do io=1,Norbs

            Bfac(io)=B_num2(io)+FinalB(io)
!             Bfac(io)=0.d0
!          write(121,*)io,Bfac(io)
               
          end do

!****************************************************
         num=zero
         do io=1,Norbs
           do ie=-N,N
            do jo=1,Norbs
             if(jo.ne.io)then
            num(io,ie)=num(io,ie)+sigma2(io,jo,ie)
             end if
            end do
          end do
        end do

!          do ie=-N,N
!        write(107,*)ie,dimag(num(1,ie)),dreal(num(1,ie))
!          end do


!*****************************************************
! calculation of self-energy
           
!	  dfac=0.0d0
          sigma_dyn=cmplx(0.d0,0.d0)
         if(phsym==0)then
         do io=1,Norbs
!          if(nf(io).eq.0.50000) then
           do ie=-N,N
!          if(nf(io).ne.0.500) then
         sigma_dyn(io,ie)=(Afac(io)*num(io,ie))/ &
          (1.d0-(Bfac(io)*num(io,ie)))
           
         sigmanew(io,ie)=Hart(io)+sigma_dyn(io,ie)
          sigma(io,ie)=dfac*sigma(io,ie) &
         +(1.d0-dfac)*sigmanew(io,ie)
           if(dimag(sigma(io,ie)).gt.0.d0) then
             if(task_id==0) write(6,*) 'Acausal self energy',io,ie
	   end if
          end do
          
        end do
          else

           do io=1,Norbs
           do ie=-N,N
         sigma_dyn(io,ie)=(Afac(io)*num(io,ie))

         sigmanew(io,ie)=Hart(io)+sigma_dyn(io,ie)
          sigma(io,ie)=dfac*sigma(io,ie) &
         +(1.d0-dfac)*sigmanew(io,ie)
           if(dimag(sigma(io,ie)).gt.0.d0) then
             if(task_id==0)write(6,*) 'Acausal self energy',io,ie
           end if
          end do

        end do

         end if

        deallocate(rho1)
        deallocate(rho2)
        deallocate(resig)
        deallocate(rhosigma1)
        deallocate(rhosigma2)
        deallocate(rhosigma)
        deallocate(sigma2)
        deallocate(reg_sigma)
        deallocate(num)
        deallocate(sigmanew)
        deallocate(muin)
        deallocate(Afac)
        deallocate(Bfac)
        deallocate(Num1)
        deallocate(Denm)
        deallocate(naa)
        deallocate(Num2)
        deallocate(Num3)
        deallocate(Num4)
        deallocate(Two_cr)
        deallocate(D_n)
        deallocate(B_num1)
        deallocate(B_num2)
        deallocate(Afac2)
        deallocate(Hart)
        deallocate(Num5)
        deallocate(Num6)
        deallocate(Num7)
        deallocate(Num8)
        deallocate(Num9)
        deallocate(FinalB)

 299	return
	end


!	************************************************************

	subroutine convolve(Ntmp,rho1,rho2,rho3)
!			      /	
!	Calculate   rho3(w) = | dw' rho1(w') rho2(w+w')
!			      /

         use Global
	real*8 r1,r2,r3,r4,r5,r6
        integer i,j,info,infot,Ntmp
        integer llim,ulim,lj,locate,jmin,jmax,Nm2
        real*8 :: rho1(-Ntmp:Ntmp),rho2(-Ntmp:Ntmp), &
        rho3(-Ntmp:Ntmp)
        real*8, allocatable :: wl(:)
        real*8 rho2ji,wji,wjmin,wjmax
        real time1,time2,time3

        Nm2=2*N+1
        allocate(wl(Nm2),stat=info)
        infot=infot+info
        if(infot.ne.0) write(6,*) 'Allocation failed in',task_id


!        time1=secnds(0.0)

        do i=-N,N
            wl(i+N+1)=w(i)
        end do

        do i=-N,N
          wjmin=max(-w(i)-w(N),-w(N))
          wjmax=min(-w(i)+w(N),w(N))
          jmin=locate(2*N+1,wl,1,2*N+1,wjmin)-N-1+1
          jmax=locate(2*N+1,wl,1,2*N+1,wjmax)-N-1
          r1=0.d0
          llim=-N
          do j=jmin,jmax
            wji=w(j)+w(i)
            do while(wji.gt.w(llim+1))
              llim=min(llim+1,N)
            end do
            lj=llim
!           This implies that  ==> w(lj)<= w_j-w_i <=w(lj+1)
            rho2ji=rho2(lj)+(rho2(lj+1)-rho2(lj))*((wji-w(lj))/ &
                                  (w(lj+1)-w(lj)))
            r1=r1+rho1(j)*rho2ji*dw(j)
          end do
          rho3(i)=r1
        end do

!        time2=secnds(time1)
!        write(6,*) 'convol',time2

        deallocate(wl)

        return
        end

!       **************************************************************
	subroutine kktransf(rhosigma,rlsigma)
	use Global
        implicit none
	real*8 r1,r2,r3,r4,r5,r6
        integer i,j,info
	integer llim,ulim,Ntmp
	real*8 :: rhosigma(-Nm:Nm),rlsg(-Nm:Nm+1), &
                  rlsigma(-Nm:Nm),wo(-Nm:Nm)
        real*8 :: woj,dspj,resgj,spi,dspi,dwi,woji1,woji2

!		      oo		
!		       /			
!	rlsigma(w)= -P | dw' rhosigma(w')
!		       /     -----------	
!	             -oo       w' - w


!		KRAMERS-KRONIG

	llim=-N+1
	ulim=N		! calculate ReSig(w) for all w

	do j=llim,ulim

!	interlacing grid

	 wo(j)=.5D0*(w(j-1)+w(j))
 	 woj=wo(j)
	 dspj=rhosigma(j)-rhosigma(j-1)
	 resgj=0.0D0

 	
	  do i=-N,j-2

	  spi=rhosigma(i)
	  dspi=rhosigma(i+1)-rhosigma(i)
	  dwi=w(i+1)-w(i)
	  woji1=w(i)-woj
	  woji2=w(i+1)-woj
	  r1=dlog(woji2/woji1)


	  resgj=resgj-(spi*r1 + dspi )
	  resgj=resgj-(dspi/dwi)*(woj -w(i))*r1

	  end do

!	 skip the interval (j-1) to j
 	
	  do i=j,N-1

	  spi=rhosigma(i)
	  dspi=rhosigma(i+1)-rhosigma(i)
	  dwi=w(i+1)-w(i)
	  woji1=w(i)-woj
	  woji2=w(i+1)-woj
	  r1=dlog(woji2/woji1)

	  resgj=resgj-(spi*r1 + dspi )
	  resgj=resgj-(dspi/dwi)*(woj -w(i))*r1

	  end do


	  resgj=resgj - dspj
	  rlsg(j)=resgj
	 end do
	 rlsg(ulim+1)=rlsg(ulim)
	 rlsg(llim-1)=rlsg(llim)
	
	 do i=llim-1,ulim
	  rlsigma(i)=0.5d0*(rlsg(i)+rlsg(i+1))
	 end do


	 return
	 end


!	**************************************************************

	subroutine writeoutput
	use Global
        implicit none
        real*8 r1,r2,r3,r4,r5,r6
        integer i,j,io,jo
	include 'Glob_cons'
	
         open(unit=31,file='Gf1.dat',status='unknown')
         open(unit=40,file='Gfscript1.dat',status='unknown')
         open(unit=34,file='sig1up.dat',status='unknown')
         open(unit=35,file='sig1down.dat',status='unknown')
         if(Norbs.gt.2) then
         open(unit=32,file='Gf3.dat',status='unknown')
         open(unit=41,file='Gfscript3.dat',status='unknown')
         open(unit=36,file='sig2up.dat',status='unknown')
         open(unit=37,file='sig2down.dat',status='unknown')
         end if
         if(Norbs.gt.4) then
         open(unit=33,file='Gf5.dat',status='unknown')
         open(unit=42,file='Gfscript5.dat',status='unknown')
         open(unit=38,file='sig3up.dat',status='unknown')
         open(unit=39,file='sig3down.dat',status='unknown')
         end if


	do i=-N,N
           write(31,'(f12.6,e25.12,e25.12)')w(i), &
       -dimag(Gf(1,i))/pi,dreal(Gf(1,i))
           write(34,'(f20.12,e25.12,e25.12)') w(i) &
          ,dimag(sigma(1,i)),dreal(sigma(1,i))
   	   write(35,'(f20.12,e25.12,e25.12)') w(i) &
          ,dimag(sigma(2,i)),dreal(sigma(2,i))
   	   write(40,'(f12.6,e25.12,e25.12)')w(i), &
       -dimag(Gfscript(1,i))/pi,dreal(Gfscript(1,i))
	enddo
	close(31)
	close(34)
	close(35)
	close(40)

        if(Norbs.gt.2) then
	do i=-N,N
           write(32,'(f12.6,e25.12,e25.12)')w(i), &
       -dimag(Gf(3,i))/pi,dreal(Gf(3,i))
           write(36,'(f20.12,e25.12,e25.12)') w(i) &
          ,dimag(sigma(3,i)),dreal(sigma(3,i))
           write(37,'(f20.12,e25.12,e25.12)') w(i) &
          ,dimag(sigma(4,i)),dreal(sigma(4,i))
           write(41,'(f12.6,e25.12,e25.12)')w(i), &
       -dimag(Gfscript(3,i))/pi,dreal(Gfscript(3,i))
	enddo
        close(32)
        close(36)
        close(37)
        close(41)
        end if

        if(Norbs.gt.4) then
	do i=-N,N
           write(33,'(f12.6,e25.12,e25.12)')w(i), &
       -dimag(Gf(5,i))/pi,dreal(Gf(5,i))
           write(38,'(f20.12,e25.12,e25.12)') w(i) &
          ,dimag(sigma(5,i)),dreal(sigma(5,i))
           write(39,'(f20.12,e25.12,e25.12)') w(i) &
          ,dimag(sigma(6,i)),dreal(sigma(6,i))
           write(42,'(f12.6,e25.12,e25.12)')w(i),  &
       -dimag(Gfscript(5,i))/pi,dreal(Gfscript(5,i)) 
	enddo
	close(33)
        close(38)
        close(39)
        close(42)
        end if
  
       open(unit=43,file='initials.dat',status='unknown')
        do io=1,Norbs
   	write(43,*) nf(io)
    	end do
        write(43,*) mu0,mu_c,shift,ntot_input
 	close(43)
	return
	end

!	************************************************************

	complex*16 function gendos(z)

	use Global
	real*8 r1,r2,r3,r4,r5,r6
        integer i,j,info
	real*8 zr,zi,sgn
        real*8, allocatable::rho0(:),epg(:)
	complex*16 z

	external sgn
	include 'Glob_cons'

        allocate(rho0(Nmdos),stat=info)
        allocate(epg(Nmdos),stat=info)

	zr=dreal(z)
	zi=dimag(z)
	if(zi.eq.0.d0) zi=1.D-10


	r1=0.d0
	r2=0.d0
	do i=1,Ndos-1
	   r3=0.5d0*(rho0(i)+rho0(i+1))
	   r4=dlog(((zr-epg(i+1))**2+ zi**2)/((zr-epg(i))**2+ zi**2))
	   r5=datan((epg(i+1)-zr)/dabs(zi))-datan((epg(i)-zr)/dabs(zi))
	   r1=r1-0.5d0*r3*r4
	   r2=r2-sgn(zi)*r3*r5
	end do
	gendos=r1+ii*r2

        deallocate(rho0)
        deallocate(epg)

	return
	end

!	************************************************************

	subroutine lattice_input()
       	use Global
        implicit none
        integer ik,jo,ko,i,j,info
        real*8,allocatable ::  KX(:),KY(:),KZ(:),EK(:)
        real*8,allocatable :: H0_hyb(:,:),H0_RS(:,:)


        allocate(KX(-Nkp:Nkp),stat=info)
        allocate(KY(-Nkp:Nkp),stat=info)
        allocate(KZ(-Nkp:Nkp),stat=info)
        allocate(EK(N_PTS),stat=info)
        allocate(H0_hyb(Norbs,-Nm:Nm),stat=info)
        allocate(H0_RS(Norbs,-Nm:Nm),stat=info)

        if(LDA_DMFT==1)then ! Real material calculation
         open(unit=10,file='HK_data',status='unknown')             
	   do ik=1,NQPTS
	     do jo=1,Norbs
	       do ko=1,Norbs
                   read(10,*)Ham0(ik,jo,ko)
                  end do
                end do
              end do
            close(10)
           if(task_id==0) write(6,*) 'Hamiltonian constructed'

         else ! 2 or 3-D square lattice calculation

            if(task_id==0) then
            if(Ndim.ne.2) then
              write(6,*) 'Default lattice is 3-D cubic lattice.'
            else
              write(6,*) 'Lattice is 2-D square lattice.'
            end if
            end if

            open(unit=11,file='HK_constructed')      
            do i=1,Norbs
              do j=1,Norbs
               read(11,*) H0_hyb(i,j)
              end do
            end do
           do i=1,Norbs
                OrbE(i)=H0_hyb(i,i)
                ep_f(i)=OrbE(i)
            end do

            do i=1,Norbs
             do j=1,Norbs
               read(11,*) H0_RS(i,j)
             end do
           end do
           close(11)
            
          do i=-Nkp,Nkp
            KX(i)=(1.d0/dfloat(Nkp))*dfloat(i)
            KY(i)=KX(i)
            KZ(i)=KX(i)
          end do


          do ik=1,N_PTS
             if(Ndim==2) then
             Ek(ik)=dcos(pi*KX(ik))+dcos(pi*KY(ik)) 
             else
             Ek(ik)=dcos(pi*KX(ik))+dcos(pi*KY(ik))+dcos(pi*KZ(ik))
             end if
             do jo=1,Norbs
               do ko=1,Norbs
                  if(jo.eq.ko) then
                 Ham_model(ik,jo,ko)=-2.d0*H0_RS(jo,ko)*Ek(ik)+ep_f(jo)
                  else
                  Ham_model(ik,jo,ko)=-2.d0*H0_RS(jo,ko)*Ek(ik)+ &
                                       H0_hyb(jo,ko)
                  end if
               end do
             end do
          enddo !  ik for momentum k
          if(task_id==0) write(6,*) 'Hamiltonian constructed'

          end if

        deallocate(KX)
        deallocate(KY)
        deallocate(KZ)
        deallocate(EK)
        deallocate(H0_hyb)
        deallocate(H0_RS)
           
	 return
	 end


!	************************************************************
	subroutine readinput()
	use global
	implicit none
	integer i,j,N1,locate,sp,flag,k,jo,ko,info
        integer infot
	real*8 r1,r2,r3
        real*8 fermic,n0up,n0dn,wi,sgn
	real*8, allocatable :: wl(:),dwl(:)
	complex*16, allocatable:: sigl(:)
	character*100 str

	external fermic,locate

	include 'Glob_cons'


	  i=1
	  open(unit=24,file='sig1up.dat',status='unknown')
 23	  read(24,*,end=24) r1,r2,r3
	  i=i+1
	  goto 23
 24	  close(24)
 	  N1=i-1

        infot=0
        allocate(wl(N1),stat=info)
        infot=infot+info
        allocate(dwl(N1),stat=info)
        infot=infot+info
        allocate(sigl(N1),stat=info)
        infot=infot+info
        if(infot.ne.0) write(6,*) 'Allocation failed in',task_id

         j=0
         wl=0.d0
         sigl=0.d0        
         r1=0.d0
         r2=0.d0
         r3=0.d0
         N1=0
	  i=1
	  open(unit=24,file='sig1up.dat',status='unknown')
 11	  read(24,*,end=12) r1,r2,r3
 	  wl(i)=r1
	  sigl(i)=dcmplx(r3,r2)
	  i=i+1
	  goto 11
 12	  close(24)
 	  N1=i-1

	  do i=-N,N
	    wi=w(i)
	    j=locate(N1,wl,1,N1,wi)
	      sigma(1,i)=sigl(j)+(sigl(j+1)-sigl(j))* &
     		        ((wi-wl(j))/(wl(j+1)-wl(j)))
     	  end do
           
          j=0
          wl=0.d0
          sigl=0.d0
          r1=0.d0
          r2=0.d0
          r3=0.d0
          N1=0
 	  i=1
	  open(unit=25,file='sig1down.dat',status='unknown')
 13	  read(25,*,end=14) r1,r2,r3
 	  wl(i)=r1
	  sigl(i)=dcmplx(r3,r2)
	  i=i+1
	  goto 13
 14	  close(25)
 	  N1=i-1

	  do i=-N,N
	    wi=w(i)
	    j=locate(N1,wl,1,N1,wi)
	      sigma(2,i)=sigl(j)+(sigl(j+1)-sigl(j))* &
     		        ((wi-wl(j))/(wl(j+1)-wl(j)))
     	  end do
!          
        
        if(Norbs.gt.2) then
          j=0
          wl=0.d0
          sigl=0.d0
          r1=0.d0
          r2=0.d0
          r3=0.d0
          N1=0
          i=1
         open(unit=26,file='sig2up.dat',status='unknown')
 15       read(26,*,end=16) r1,r2,r3
          wl(i)=r1
          sigl(i)=dcmplx(r3,r2)
          i=i+1
          goto 15
 16       close(26)
          N1=i-1

          do i=-N,N
            wi=w(i)
            j=locate(N1,wl,1,N1,wi)
              sigma(3,i)=sigl(j)+(sigl(j+1)-sigl(j))* &
                       ((wi-wl(j))/(wl(j+1)-wl(j)))
          end do
          j=0
          wl=0.d0
          sigl=0.d0
          r1=0.d0
          r2=0.d0
          r3=0.d0
          N1=0  
          i=1
          open(unit=27,file='sig2down.dat',status='unknown')
 17       read(27,*,end=18) r1,r2,r3
          wl(i)=r1
          sigl(i)=dcmplx(r3,r2)
          i=i+1
          goto 17
 18       close(27)
          N1=i-1

          do i=-N,N
            wi=w(i)
            j=locate(N1,wl,1,N1,wi)
              sigma(4,i)=sigl(j)+(sigl(j+1)-sigl(j))* &
                       ((wi-wl(j))/(wl(j+1)-wl(j)))
          end do

          end if

          if(Norbs.gt.4)then
          j=0
          wl=0.d0
          sigl=0.d0
          r1=0.d0
          r2=0.d0
          r3=0.d0
          N1=0
          i=1
          open(unit=28,file='sig3up.dat',status='unknown')
 19       read(28,*,end=20) r1,r2,r3
          wl(i)=r1
          sigl(i)=dcmplx(r3,r2)
          i=i+1
          goto 19
 20       close(28)
          N1=i-1

          do i=-N,N
            wi=w(i)
            j=locate(N1,wl,1,N1,wi)
              sigma(5,i)=sigl(j)+(sigl(j+1)-sigl(j))* &
                       ((wi-wl(j))/(wl(j+1)-wl(j)))
          end do
           
          j=0
          wl=0.d0
          sigl=0.d0
          r1=0.d0
          r2=0.d0
          r3=0.d0
          N1=0
          i=1
          open(unit=29,file='sig3down.dat',status='unknown')
 21       read(29,*,end=22) r1,r2,r3
          wl(i)=r1
          sigl(i)=dcmplx(r3,r2)
          i=i+1
          goto 21
 22       close(29)
          N1=i-1

          do i=-N,N
            wi=w(i)
            j=locate(N1,wl,1,N1,wi)
              sigma(6,i)=sigl(j)+(sigl(j+1)-sigl(j))* &
                       ((wi-wl(j))/(wl(j+1)-wl(j)))
          end do

          end if
           
	  if(task_id==0) write(6,*) 'Reading done'

        deallocate(wl)
        deallocate(dwl)
        deallocate(sigl)

	return
	end


!	************************************************************
        subroutine findgf()
        use Global
        implicit none
        complex*16 sem,gauss,gamma,gendos
        integer io,i,jo,pw,ie
        real*8 m2,m4,r1,r2norm

        external gauss,sem,gendos

        include 'Glob_cons'


!       Calculate second and fourth moments of the Gaussian DOS
       
          if(idos.eq.1) then
           m2=t**2/2.d0
           m4=t**4*3.d0/4.d0
         else if(idos.eq.2) then
           m2=t**2/4.d0
           m4=t**4/8.d0
         end if

         if(idos==1)then
         do io=1,Norbs
          do i=-N,N
            z1=w(i)+ii*eta
            gamma=z1+mu_c-ep_f(io)-sigma(io,i)
          
          if(zabs(gamma).le.1.D6) then
              Gf(io,i)=sem(gamma)
             else
            Gf(io,i)=1.d0/gamma+m2/gamma**3   ! From moment expansion
            end if
        end do
      end do
        else if(idos==2)then

          do io=1,Norbs
          do i=-N,N
            z1=w(i)+ii*eta
            gamma=z1+mu_c-ep_f(io)-sigma(io,i)

!          if(zabs(gamma).le.1.D6) then
            Gf(io,i)=gauss(gamma)
            end do
         end do

          end if


        
       do io=1,Norbs
        do i=-N,N
        Gfscript(io,i)=(1.d0+(sigma(io,i)+ep_f(io)-mu_c+mu0)*Gf(io,i))
        Gfscript(io,i)=Gf(io,i)/Gfscript(io,i)
        if(dimag(Gfscript(io,i)).gt.0.d0)  &
            Gfscript(io,i)=dcmplx(dreal(Gfscript(io,i)),-1.D-9)
        end do
      end do

      return
      end

!******************************************************
!***************************************************************
        subroutine realtoImg()

        use Global
        implicit none
        integer info
        complex*16 :: msum1,msum2
        complex*16, allocatable :: Z(:),Gfwn(:,:),Sigmawn(:,:)
        integer i,j,k,ulim,nwn
        real*8 m2,m4,r1,r2
        real*8, allocatable :: Y(:),X(:)
        include 'Glob_cons'

        nwn=300

        allocate(z(-nwn:nwn),stat=info)
        allocate(Gfwn(Norbs,-nwn:nwn),stat=info)
        allocate(Sigmawn(Norbs,-nwn:nwn),stat=info)
        allocate(X(-nwn:nwn),stat=info)
        allocate(Y(-nwn:nwn),stat=info)

        if(Norbs.eq.2) then 
          ulim=1
        else if(Norbs.eq.4) then
          ulim=3
        end if


        open(unit=125,file='Matsubara_green.dat')
        open(unit=126,file='Matsubara_self.dat')

        do i=-nwn,nwn
          
         x(i)=0.0d0
         y(i)=pi*temp*(2*i+1)
         z(i)=dcmplx(x(i), y(i))
         do k=1,ulim,2
            msum1=dcmplx(0.d0,0.d0)
            msum2=dcmplx(0.d0,0.d0)
            do j=-N,N
         msum1=msum1+((-dimag(Gf(k,j))/pi)*dw(j))/(z(i)-w(j))
         msum2=msum2+((-dimag(sigma(k,j))/pi)*dw(j))/(z(i)-w(j))
             end do
            Gfwn(k,i)=msum1
            Sigmawn(k,i)=msum2
         
        end do
        if(Norbs.le.2) then
        write(125,*)dimag(z(i)),dimag(Gfwn(1,i))
        write(126,*)dimag(z(i)),dimag(Sigmawn(1,i))
        else if(Norbs.eq.4) then
        write(125,*)dimag(z(i)),dimag(Gfwn(1,i)),dimag(Gfwn(3,i))  
        write(126,*)dimag(z(i)),dimag(Sigmawn(1,i)),dimag(Sigmawn(3,i))        
        else if(Norbs.eq.6) then
        write(125,'(4f16.8,1x)')dimag(z(i)),dimag(Gfwn(1,i)), &
                               dimag(Gfwn(3,i)),dimag(Gfwn(5,i))  
        write(126,'(4f16.8,1x)')dimag(z(i)),dimag(Sigmawn(1,i)), &
                               dimag(Sigmawn(3,i)),dimag(Sigmawn(5,i))  
        end if
        end do         
        deallocate(z)
        deallocate(Gfwn)
        deallocate(Sigmawn)
        deallocate(X)
        deallocate(Y)
           
        return
        end
!****************************************************
      
        subroutine fill_U_matrix()

        use global
        implicit none
        integer io,jo
               
            
        if(Jflag==0) then
           do io=1,Norbs
              do jo=1,Norbs
                if(io.eq.jo) then
                   U_ab(io,jo)=0.d0
                else
                   U_ab(io,jo)=U
                end if
              end do
           end do
         else if(Jflag==1) then                       
           do io=1,Norbs
             do jo=1,Norbs
                if((io+jo.eq.3).or.(io+jo.eq.7)) then
                  U_ab(io,jo)=U
                else if(io+jo.eq.5) then
                  U_ab(io,jo)=U-2*(U/4.d0)
                else if(((io+jo.eq.4).or.(io+ &
                     jo.eq.6)).and.(io.ne.jo)) then
                  U_ab(io,jo)=U-3*(U/4.d0)
                else 
                  U_ab(io,jo)=0.d0
                end if
             end do
           end do
         else if(Jflag==2) then
           do io=1,Norbs
             do jo=1,Norbs
               if((io+jo.eq.3).or.(io+jo.eq.7)) then
                  U_ab(io,jo)=U
               else if(io+jo.eq.5) then
                  U_ab(io,jo)=U-2*(U/3.5d0)
               else if(((io+jo.eq.4).or.(io+ &
                    jo.eq.6)).and.(io.ne.jo)) then
                  U_ab(io,jo)=U-3*(U/3.5d0)
               else
                  U_ab(io,jo)=0.d0
               end if
             end do
           end do
         else if(Jflag==3) then
           do io=1,Norbs
             do jo=1,Norbs
               if((io+jo.eq.3).or.(io+jo.eq.7)) then
                  U_ab(io,jo)=U
               else if(io+jo.eq.5) then
                  U_ab(io,jo)=U-2*(U/3.2d0)
               else if(((io+jo.eq.4).or.(io+ &
                     jo.eq.6)).and.(io.ne.jo)) then
                   U_ab(io,jo)=U-3*(U/3.2d0)
               else
                   U_ab(io,jo)=0.d0
               end if
             end do
           end do
         else if(Norbs==6)then
               U_ab(1,1)=0.0d0
               U_ab(2,2)=0.0d0
               U_ab(3,3)=0.0d0
               U_ab(4,4)=0.0d0
               U_ab(5,5)=0.0d0
               U_ab(6,6)=0.0d0
!************************************
               U_ab(1,2)=U
               U_ab(1,3)=U-3.d0*J_H
               U_ab(1,4)=U-2.d0*J_H
               U_ab(1,5)=U_ab(1,3)
               U_ab(1,6)=U_ab(1,4)
!*****************************************
               U_ab(2,1)=U_ab(1,2) 
               U_ab(2,3)=U_ab(1,4)
               U_ab(2,4)=U_ab(1,3)
               U_ab(2,5)=U_ab(1,4)
               U_ab(2,6)=U_ab(1,3)
!*******************************************
               U_ab(3,1)=U_ab(1,3)
               U_ab(3,2)=U_ab(1,4)
               U_ab(3,4)=U_ab(1,2)
               U_ab(3,5)=U_ab(1,3)
               U_ab(3,6)=U_ab(1,4)
!*********************************************
               U_ab(4,1)=U_ab(1,4)
               U_ab(4,2)=U_ab(1,3)
               U_ab(4,3)=U_ab(1,2)
               U_ab(4,5)=U_ab(1,4)
               U_ab(4,6)=U_ab(1,3)
!*********************************************
               U_ab(5,1)=U_ab(1,3)
               U_ab(5,2)=U_ab(1,4)
               U_ab(5,3)=U_ab(1,3)
               U_ab(5,4)=U_ab(1,4)
               U_ab(5,6)=U_ab(1,2)
!********************************************
               U_ab(6,1)=U_ab(1,4)
               U_ab(6,2)=U_ab(1,3)
               U_ab(6,3)=U_ab(1,4)
               U_ab(6,4)=U_ab(1,3)
               U_ab(6,5)=U_ab(1,2)
!******************************************
             
         else if(Norbs==4)then
               U_ab(1,1)=0.0d0
               U_ab(2,2)=0.0d0
               U_ab(3,3)=0.0d0
               U_ab(4,4)=0.0d0
               !**********************
               U_ab(1,2)=U
               U_ab(1,3)=U-3.d0*J_H
               U_ab(1,4)=U-2.d0*J_H
               !***********************
               U_ab(2,1)=U_ab(1,2) 
               U_ab(2,3)=U_ab(1,4)
               U_ab(2,4)=U_ab(1,3)
               !**********************
               U_ab(3,1)=U_ab(1,3)
               U_ab(3,2)=U_ab(1,4)
               U_ab(3,4)=U_ab(1,2)
               !***********************
               U_ab(4,1)=U_ab(1,4)
               U_ab(4,2)=U_ab(1,3)
               U_ab(4,3)=U_ab(1,2)
               !***********************
             
         else if(Norbs==2)then
               U_ab(1,1)=0.0d0
               U_ab(2,2)=0.0d0
               !************************************
               U_ab(1,2)=U
               U_ab(2,1)=U_ab(1,2) 
             
           end if 

        return
        end
