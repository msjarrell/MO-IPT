c	**********************************************	  
	double precision function sgn(x)
	double precision x

	if(x.lt.0.D0) then
	  sgn=-1.D0
	else
	  sgn=1.D0
	end if

	return
	end

c	**********************************************	  
	double precision function theta(x)
	double precision x

	if(x.lt.0.D0) then
	  theta=0.D0
	else if(x.eq.0.d0) then
	  theta=0.5d0
	else
	  theta=1.d0
	end if

	return
	end

c	**********************************************	  



	double precision function fermic(x,beta)
	double precision x,beta
c	fermi function which avoids underflow/overflow.  If beta=0,
c	it is assumed that T=0 is meant!

	if(beta.eq.0.0D0) then
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0
	  else if(x.eq.0.0D0) then
	    fermic=0.5D0
	  else
	    fermic=0.0D0
	  end if
	else
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0/(dexp(beta*x)+1.0D0)
	  else
	    fermic=dexp(-beta*x)/(dexp(-beta*x)+1.0D0)
	  end if
	end if
	return
	end

c	**********************************************	  
	double precision function dfermic(x,beta)
	double precision x,beta
c
c	This function returns the minus derivative of the fermi function
c
c                               1
c	fermic(x,beta)= ----------------
c                       exp(-beta x) + 1
c
c       d fermic(x,beta)     -beta exp(-beta x)
c	----------------  = --------------------- = -dfermic(x,beta)
c              d x           (exp(-beta x) + 1)^2
c
c
	if(x.lt.0.0D0) then	  
	  dfermic=beta*dexp(beta*x)/(dexp(beta*x)+1.0D0)**2
        else
	  dfermic=beta*dexp(-beta*x)/(dexp(-beta*x)+1.0D0)**2
	end if

	return
	end

c	**********************************************	  
c
	INTEGER FUNCTION LOCATE(N,xx,ll,ul,x)
	integer N,ll,ul
	real*8 xx(1:N),x
c
c	given an array xx(1:n),and given a value x, returns a value j 
c 	such that x is between xx(j) and xx(j+1). xx(1:n) must be
c 	monotonic.

	integer jl,jm,ju


	if(x.lt.xx(ll)) then
	  locate=ll
c	  pause'input out of range in locate,left'
	  return
	end if

	if(x.gt.xx(ul)) then
	  locate=ul-1
c	  pause'input out of range in locate,right'
	  return
	end if

	jl=ll
	ju=ul
        do while((ju-jl) .gt. 1) 
	  jm=(ju+jl)/2
	  if (x .ge. xx(jm)) then
	    jl=jm
	  else
	    ju=jm
	 endif
	end do 
	locate=jl
c
	return
	end

c	****************************************************************
	double complex function gauss(z)
c	This block calculates -i*sqrt(pi)*w(z)
	double complex z,ii
	double precision pi2,x,y,t
	logical flag
	parameter(pi2=1.7724539D0)
	ii=dcmplx(0.D0,1.D0)

	t=1.d0
	if(dimag(z).lt.0.0D0) then
	  call wofz(dreal(z),-dimag(z),x,y,flag)
	  gauss=-pi2*ii*dcmplx(x,-y)
	else
	  call wofz(dreal(z),dimag(z),x,y,flag)
	  gauss=-pi2*ii*dcmplx(x,y)
	end if

	if(flag) write(6,*) 'error in cerfjd'
	return
	end


c	**************************************************
        complex*16 function sem(z)
c       This block calculates Hilbert transform for a semicircular 
c       DOS with t^*=t
c       In general sem=8*(z-sqrt(z**2-D**2/4))/D**2
        double complex z,z1,z2,ii,sem1,sem2,sz2
        double precision dband,pi,t,t2

        include 'Glob_cons'

	t=1.d0
!        t2=2.d0*t
	dband=t
        z1=dcmplx(dband,0.D0)
					             
        z2=z**2-z1**2
				                  
        sz2=zsqrt(z2)
												               
        sem1=2.D0/(z+sz2)
												                    
        sem2=2.D0/(z-sz2)
        
   
        if(dimag(sem1).le.0.D0) then
             sem=sem1
        else if(dimag(sem2).le.0.D0) then
             sem=sem2
        else 
             write(6,*) 'no causal root found'
             stop 
        end if
    
        return
        end  


c	*************************************************************
	
	  double complex function flat(z)
	  double complex z,ii
	  double precision rews,imws,rez,imz,pi
	  double precision dband,epsvh,zr,zi,acc,dnorm,rho0

	  common/dos/dband,epsvh,zr,zi,acc,dnorm,nnorm,idos

	
c	  Flat Band from -D/2 to D/2 where 1.D0
	  
	  include 'Glob_cons'

	  rez=dreal(z)
	  imz=dimag(z)
	  
	  if(dabs(rez).eq.dband/2.D0) rez=rez+1.D-6
	  if(dabs(imz).le.1.D-6) then
	    rews=dlog(dabs((rez+dband/2.D0)/(rez-dband/2.D0)))
	    if(dabs(rez).ge.dband/2.D0) then
	      imws=0.D0
	    else
	    imws=-pi
	    end if
	  else
	    rews=0.5D0*dlog(((rez+dband/2.D0)**2+imz**2)/
     &	         	   ((rez-dband/2.D0)**2+imz**2))

	    imws=-(datan((dband/2.D0-rez)/imz)+
     &		  datan((dband/2.D0+rez)/imz))
	  end if

 355	  flat=dcmplx(rews,imws)

 	  return
	  end

c	*************************************************************

C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
      SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
C
C  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
C  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
C  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
C  MEANS SQRT(-1).
C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
C  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
C  OF THE FUNCTION.
C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
C
C
C  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
C
C
C  PARAMETER LIST
C     XI     = REAL      PART OF Z
C     YI     = IMAGINARY PART OF Z
C     U      = REAL      PART OF W(Z)
C     V      = IMAGINARY PART OF W(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
C  PUT TO 0 UPON UNDERFLOW;
C
C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
C  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
C
*
*
*
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      LOGICAL A, B, FLAG
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     *           RMAXREAL = 0.5D+154,
     *           RMAXEXP  = 708.503061461606D0,
     *           RMAXGONI = 3.53711887601422D+15)
*
      FLAG = .FALSE.
*
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
*
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
      QRHO = X**2 + Y**2
*
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
*
      A     = QRHO.LT.0.085264D0
*
      IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
 10     CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
*
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
*
      ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
*
        B = (H.GT.0.0)
*
        IF (B) QLAMBDA = H2**KAPN
*
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
*
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
 11     CONTINUE
*
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
*
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
*
      END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
      IF (YI.LT.0.0) THEN
*
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
*
C
C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
*
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
*
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
*
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END

c	************************************************************

      SUBROUTINE LINT(XA,YA,Y2A,N,KLO,X,Y) 
      INTEGER N,KLO,KHI,K
      DOUBLE PRECISION XA(N),YA(N),Y2A(N),X,Y,H,B
      
      IF(KLO.EQ.0) THEN		! MEANS UNSET IN MAIN PROGRAM
      KLO=1 
      END IF
      KHI=N 
1     IF (KHI-KLO.GT.1) THEN 
        K=(KHI+KLO)/2 
        IF(XA(K).GT.X)THEN 
          KHI=K 
        ELSE 
          KLO=K 
        ENDIF 
      GOTO 1 
      ENDIF 
      H=XA(KHI)-XA(KLO) 
      IF (H.EQ.0.D0) then
         write(6,*) 'Bad XA input.' 
         stop
      END IF
      B=(X-XA(KLO))
      Y=YA(KLO)+Y2A(KLO)*B
      RETURN 
      END 

c	**********************************************	  
      SUBROUTINE LINE(X,Y,N,Y2) 
      PARAMETER (NMAX=5000) 
      INTEGER N,I
      DOUBLE PRECISION X(N),Y(N),Y2(N)

      DO I=1,N-1
        Y2(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      END DO
      Y2(N)=Y2(N-1)
        
      RETURN 
      END 

c	**********************************************	  

!  Inversion of a complex matrix a(m,m) 
       subroutine cmatinv(a,ai,m) 
       implicit none
       complex*16:: a(m,m),ai(m,m),zero 
       integer m,i,j,k 
       zero=dcmplx(0.d0,0.d0)
       ai=a
       do k=1,m 
        do j=1,m
        if(j.ne.k) then
            if (ai(k,k).eq.zero) ai(k,k)=dcmplx(1.d-9,1.D-10)
            ai(k,j)=ai(k,j)/ai(k,k) 
        endif
       end do 
       ai(k,k)=1.0d0/ai(k,k) 
       do i = 1,m 
        if (i.ne.k) then 
            do j=1,m 
                if (j.ne.k) ai(i,j)=ai(i,j)-ai(k,j)*ai(i,k)
             end do
         end if
        end do 
        do i=1,m 
        if (i.ne.k) ai(i,k)=-ai(i,k)*ai(k,k) 
       enddo
       end do 
       end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         real*8 function  broyden(tmp_muc)
         use Global
         implicit none
         include'mkl_rci.fi'
         integer M,L,K
         PARAMETER(L = 1)
         PARAMETER(M = 1)
         real*8 x(L),fjac(m,L),eps,f(m),D(L),U_(L)
         real*8 C,Y(L,L),Z(L,L),xold(L),tmp_muc
         integer ipiv(N),info,lwork
         PARAMETER(lwork=L*L)
         integer work(lwork)
         external fcn
         eps=1.0d-8
         x(1)=tmp_muc
         call fcn(m,l,x,f)

         if(DJACOBI (fcn,l,m,fjac,x,eps)
     .  /=TR_SUCCESS)then
        write(*,*) '|ERROR IN DJACOBI'
        endif

        call dgetrf(M,L,fjac,L,ipiv,info)

        if(info==0) then
        call dgetri(M,fjac,L,ipiv,work,lwork,info)
        else
        write(*,*) "problem with matrix inversion",info
        endif

        IF (info.NE.0) THEN
        stop 'Matrix inversion failed!'
        ELSE
        PRINT '(" Inverse Successful ")'
        ENDIF 


        D=-matmul(fjac,f) 

        x=x+D

!        write(*,"(2a,3F15.4)") "X values", "0", x


        k=0

        do
        call fcn(m,l,x,f)
        U_=matmul(fjac,f)
        C=dot_product(D,U_+D)

        Y=0.d0;Z=0.d0

        Y(1,:)=D
        Z(1,:)=U_

       Y=transpose(Y)
       fjac=fjac-(1/C)*matmul(matmul(Z,Y),fjac)
       k=k+1
       D=-matmul(fjac,f)
       xold=x
       x=x+D
       if(dot_product(x-xold,x-xold)<10.d-7) exit
       enddo

        broyden=x(1)

        

         return
         end 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         subroutine fcn(size_f,size_x,var,ff)
         use Global
         implicit none
         integer :: size_f,size_x,io
         real*8 var(size_x),ff(size_f)
         real*8 nf_tot_tmp
         include 'Glob_cons'
         mu_c=var(1)
         call findgf()
         call occupancies()
          nf_tot_tmp=0.d0
          do io=1,Norbs
          nf_tot_tmp=nf_tot_tmp+nf(io)
         end do

         ff(1)=ntot-nf_tot_tmp
         write(*,*)'funcv1 mu_c=',mu_c,' fvec',ff(1)              

        return
        end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
