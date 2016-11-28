!	************************************************************


	subroutine makegrid

	 use Global 
	implicit none
	real*8 r1,r2,r3,r4,r5,r6
        integer i,j,info
	integer grid,M1,M2,M3,M4,M5,M6,M7,lN,rN,i1,sp, &
     	i2,i3,j1,j2
	real*8 k,wmax,wmax1,tep,uw,fermic,dwi

	external grid,fermic

!	t=1.d0
	wmax=3.0d0*max(t,U/2.d0)
	wmax1=15.0d0*max(t,U/2.d0)
        

!       SET UP THE FREQUENCY GRID

        if(uniform_grid==1)then
         Nm=5000
         dwi=0.0025d0
         allocate(w(-Nm:Nm),stat=info)
         allocate(dw(-Nm:Nm),stat=info)
         allocate(ferm(-Nm:Nm),stat=info)
         w(0)=0.0d0
         N=Nm
         do i=-N,N
            dw(i)=dwi
            w(i)=dfloat(i)*dw(i)
   	   ferm(i)=fermic(w(i),beta)
          end do
 
        
        else
	  
!       First estimate N

        Nm=0
        k=ep2/b1
        M1=aint(dlog(b1/ep1)/k)
        if((M1/2*2-M1).ne.0) M1=M1+1
        Nm=Nm+M1

	M1=aint((b2-db2-b1)/ep2)
	if((M1/2*2-M1).ne.0) M1=M1+1
        Nm=Nm+M1

	M1=aint((2.d0*db2)/ep3)
	if((M1/2*2-M1).ne.0) M1=M1+1
        Nm=Nm+M1

	M1=aint((b3-db3-b2-db2)/ep2)
	if((M1/2*2-M1).ne.0) M1=M1+1
        Nm=Nm+M1

	M1=aint((2.d0*db3)/ep4)
	if((M1/2*2-M1).ne.0) M1=M1+1
        Nm=Nm+M1

	M1=aint((wmax-b3-db3)/ep2)
	if((M1/2*2-M1).ne.0) M1=M1+1
        Nm=Nm+M1

	M1=aint((wmax1-wmax)/ep5)
	if((M1/2*2-M1).ne.0) M1=M1+1
        Nm=Nm+M1

        Nm=Nm+200

        allocate(w(-Nm:Nm),stat=info)
        allocate(dw(-Nm:Nm),stat=info)
        allocate(ferm(-Nm:Nm),stat=info)
        

        w(0)=0.0d0
        k=ep2/b1
        M1=aint(dlog(b1/ep1)/k)
        if((M1/2*2-M1).ne.0) M1=M1+1
        do i = 1,M1
           w(i)=ep1*(dexp(i*k)-1.d0)
           w(-i)=-w(i)
           dw(i)=k*w(i)
           dw(-i)=dw(i)
        enddo
        dw(0)=dw(1)
        N=M1

	M2=grid(w(N),b2-db2,ep2,1.d0)
	N=N+M2

	M3=grid(w(N),w(N)+2.d0*db2,ep3,1.d0)
	N=N+M3

 	M4=grid(w(N),b3-db3,ep2,1.d0)
	N=N+M4

	M5=grid(w(N),w(N)+2.d0*db3,ep4,1.d0)
	N=N+M5
	
	M6=grid(w(N),wmax,ep2,1.d0)
	N=N+M6

	M7=grid(w(N),wmax1,ep5,1.d0)
	N=N+M7

	do i=-N,0
	  w(i)=-w(-i)
	  dw(i)=dw(-i)
	end do

	!if(task_id==0) write(6,*) M1,M2,M3,M4,M5,M6,M7
	if(task_id==0)  write(6,*) 'Number of frequencies=',N
        if(task_id==0)  then
        open(unit=12,file='grid.dat')
	do i=-N,N
	   write(12,*) w(i),dw(i)
	end do
	close(12)
        end if

!	Defining the fermi function as an array.
 400	do i=-N,N
	  ferm(i)=fermic(w(i),beta)
	end do

        end if

 401	return
	end

!	************************************************************
	
	integer function grid(aa,bb,ep,sign)
	 use Global

	integer M1
	real*8 aa,bb,ep,sign
	
	j=dint(sign)
	M1=aint(sign*(bb-aa)/ep)
	if((M1/2*2-M1).ne.0) M1=M1+1
	do i=N+1,N+M1
	    w(j*i)=aa+sign*ep*dfloat(i-N)
	    dw(j*i)=ep
	end do
	
	grid=M1

	return
	end

