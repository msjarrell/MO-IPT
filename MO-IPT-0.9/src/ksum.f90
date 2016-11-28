	subroutine ksummation()
	use Global
        implicit none
        include 'mpif.h'
	integer io,jo,ko,ik,ie,lo,mo,info
        integer i,j,IPT
	real*8 r1,r2,r3,r4,r5,r6
        complex*16 CZERO,x
        complex*16 INV1,INV2,INV3,det,dummy

        complex*16,allocatable:: GS(:,:,:),GI(:,:,:)
        complex*16,allocatable:: G_new(:,:),G1_new(:,:,:)
        complex*16,allocatable:: GS1(:,:),GI_new(:,:)
        complex*16,allocatable:: Gr(:,:,:)
        complex*16,allocatable:: Gr1(:,:,:)
        complex*16,allocatable:: c(:,:),cINV(:,:),XE(:,:)

        real t2,t1
	include 'Glob_cons'

        allocate(Gr(-N:N,Norbs,Norbs),stat=info)
        allocate(Gr1(-N:N,Norbs,Norbs),stat=info)
        allocate(c(Norbs,Norbs),stat=info)
        allocate(cINV(Norbs,Norbs),stat=info)
        allocate(GS(-Nm:Nm,Norbs,Norbs),stat=info)
        allocate(GI(-Nm:Nm,Norbs,Norbs),stat=info)
        allocate(G_new(Norbs,Norbs),stat=info)
        allocate(G1_new(-Nm:Nm,Norbs,Norbs),stat=info)
        allocate(GS1(Norbs,Norbs),stat=info)
        allocate(GI_new(Norbs,Norbs),stat=info)
        allocate(XE(Norbs,-Nm:Nm),stat=info)


        CZERO=dcmplx(0.d0,0.d0)
        if(task_id==0) write(6,*)N,Norbs
 
        if(LDA_DMFT==1)then
          do io=1,Norbs
	   do ie=-N,N
            z1=w(ie)+mu_c+shift+ii*eta      !-ep_f(io)+shift-LDA_mu+ii*eta     !+ep_f(io)+ii*eta 
            XE(io,ie)=z1-sigma(io,ie)
           end do
         end do
        else
           do io=1,Norbs
            do ie=-N,N
            z1=w(ie)+mu_c+ii*eta      !-ep_f(io)+shift-LDA_mu+ii*eta     !+ep_f(io)+ii*eta 
            XE(io,ie)=z1-sigma(io,ie)
           end do
          end do

         end if

         
          
           
          

          Gr=CZERO
          Gr1=CZERO
             c=CZERO
             cINV=CZERO
             x=CZERO
            call cpu_time(t1)

       jstart=-N+task_id*int((2*N)/num_tasks)+task_id
       jend=jstart+int((2*N)/num_tasks)
       if (task_id==num_tasks-1)then
         jend=N 
       endif
       if(LDA_DMFT==1)then
        if(Norbs==6)then
        do ie=jstart,jend  
            do ik=1,NQPTS ! ksum
 	         do jo=1,Norbs
 	           do ko=1,Norbs
		  c(jo,ko)=-Ham0(ik,jo,ko)
                  end do
                  c(jo,jo)=XE(jo,ie)+c(jo,jo)
                end do
	        
             x=c(1,1)*c(3,3)-c(1,3)*c(3,1)

            cINV(5,5)=c(5,5)-(c(5,1)*c(3,3)*c(1,5))/x &
          -(c(5,3)*c(1,1)*c(3,5))/x+(c(5,1)*c(1,3)*c(3,5))/x &
          +(c(5,3)*c(1,5)*c(3,1))/x

           cINV(5,5)=1.d0/cINV(5,5)

           cINV(5,1)=(c(3,1)*c(5,3)*cINV(5,5))/x &
         -(c(5,1)*c(3,3)*cINV(5,5))/x

           cINV(5,3)=(c(5,1)*c(1,3)*cINV(5,5))/x &
         -(c(5,3)*c(1,1)*cINV(5,5))/x

          cINV(1,5)=(c(1,3)*c(3,5)*cINV(5,5))/x &
        -(c(3,3)*c(1,5)*cINV(5,5))/x

           cINV(3,5)=(c(3,1)*c(1,5)*cINV(5,5))/x &
        -(c(1,1)*c(3,5)*cINV(5,5))/x

         cINV(1,1)=c(3,3)/x+(cINV(1,5)*cINV(5,1))/cINV(5,5)
         cINV(1,3)=-c(1,3)/x+(cINV(1,5)*cINV(5,3))/cINV(5,5)                        
         cINV(3,1)=-c(3,1)/x+(cINV(3,5)*cINV(5,1))/cINV(5,5)
         cINV(3,3)=c(1,1)/x+(cINV(3,5)*cINV(5,3))/cINV(5,5)
            
         cINV(2,2)=cINV(1,1)
         cINV(2,4)=cINV(1,3)
         cINV(2,6)=cINV(1,5)
         cINV(4,2)=cINV(3,1)
         cINV(4,4)=cINV(3,3)
         cINV(4,6)=cINV(3,5)
         cINV(6,2)=cINV(5,1)
         cINV(6,4)=cINV(5,3)
         cINV(6,6)=cINV(5,5)
  
       
        
                do lo=1,Norbs
                 do mo=1,Norbs
        Gr1(ie,lo,mo)=Gr1(ie,lo,mo)+(cINV(lo,mo)/dfloat(NQPTS))
                 end do
               end do
                  
	        end do
             end do


      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 
      call MPI_ALLREDUCE(Gr1,Gr,(2*N+1)*Norbs*Norbs, &
       MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)


            if(task_id==0)   write(6,*) '_____________'
	      call cpu_time(t2)
             if(task_id==0)  write(6,*)'time',t2-t1 
          
          do ie=-N,N
           Gf(1,ie)=Gr(ie,1,1)
           Gf(2,ie)=Gr(ie,2,2)
           Gf(3,ie)=Gr(ie,3,3)
           Gf(4,ie)=Gr(ie,4,4)
           Gf(5,ie)=Gr(ie,5,5)
           Gf(6,ie)=Gr(ie,6,6)
      
          end do    


          deallocate(c)
          deallocate(cINV)
               


          GI=CZERO         
          det=CZERO
          do ie=-N,N

        det=Gr(ie,1,1)*Gr(ie,3,3)-Gr(ie,1,3)*Gr(ie,3,1)
           
        GI(ie,5,5)=Gr(ie,5,5) &
       -(Gr(ie,5,1)*Gr(ie,3,3)*Gr(ie,1,5))/det &
       -(Gr(ie,5,3)*Gr(ie,1,1)*Gr(ie,3,5))/det &
       +(Gr(ie,5,1)*Gr(ie,1,3)*Gr(ie,3,5))/det &
       +(Gr(ie,5,3)*Gr(ie,1,5)*Gr(ie,3,1))/det
           
        GI(ie,5,5)=1.d0/GI(ie,5,5)
!        write(6,*)GI(ie,5,5),Gr(ie,5,5)
        
           
       GI(ie,5,1)=(Gr(ie,3,1)*Gr(ie,5,3)*GI(ie,5,5))/det &
        -(Gr(ie,5,1)*Gr(ie,3,3)*GI(ie,5,5))/det

       GI(ie,5,3)=(Gr(ie,5,1)*Gr(ie,1,3)*GI(ie,5,5))/det &
      -(Gr(ie,5,3)*Gr(ie,1,1)*GI(ie,5,5))/det

       GI(ie,1,5)=(Gr(ie,1,3)*Gr(ie,3,5)*GI(ie,5,5))/det &
      -(Gr(ie,3,3)*Gr(ie,1,5)*GI(ie,5,5))/det

       GI(ie,3,5)=(Gr(ie,3,1)*Gr(ie,1,5)*GI(ie,5,5))/det &
      -(Gr(ie,1,1)*Gr(ie,3,5)*GI(ie,5,5))/det

       GI(ie,1,1)=Gr(ie,3,3)/det &
      +(GI(ie,1,5)*GI(ie,5,1))/GI(ie,5,5)
       GI(ie,1,3)=-Gr(ie,1,3)/det &
      +(GI(ie,1,5)*GI(ie,5,3))/GI(ie,5,5)
       GI(ie,3,1)=-Gr(ie,3,1)/det &
      +(GI(ie,3,5)*GI(ie,5,1))/GI(ie,5,5)
       GI(ie,3,3)=Gr(ie,1,1)/det &
      +(GI(ie,3,5)*GI(ie,5,3))/GI(ie,5,5)          
        
         GI(ie,2,2)=GI(ie,1,1)
         GI(ie,2,4)=GI(ie,1,3)
         GI(ie,2,6)=GI(ie,1,5)
         GI(ie,4,2)=GI(ie,3,1)
         GI(ie,4,4)=GI(ie,3,3)
         GI(ie,4,6)=GI(ie,3,5)
         GI(ie,6,2)=GI(ie,5,1)
         GI(ie,6,4)=GI(ie,5,3)
         GI(ie,6,6)=GI(ie,5,5)

         end do
           
          deallocate(Gr)       
          
          do ie=-N,N
            do io=1,Norbs
              do jo=1,Norbs
                if(io.eq.jo) then
           GI(ie,io,jo)=GI(ie,io,jo)+sigma(io,ie)+mu0-mu_c
                end if 
               end do
              end do
           end do        
 
          Gs=CZERO
          det=CZERO
          do ie=-N,N
         
         det=GI(ie,1,1)*GI(ie,3,3)-GI(ie,1,3)*GI(ie,3,1)

         Gs(ie,5,5)=GI(ie,5,5) &
       -(GI(ie,5,1)*GI(ie,3,3)*GI(ie,1,5))/det &
       -(GI(ie,5,3)*GI(ie,1,1)*GI(ie,3,5))/det &
       +(GI(ie,5,1)*GI(ie,1,3)*GI(ie,3,5))/det &
       +(GI(ie,5,3)*GI(ie,1,5)*GI(ie,3,1))/det

         Gs(ie,5,5)=1.d0/Gs(ie,5,5)

       Gs(ie,5,1)=(GI(ie,3,1)*GI(ie,5,3)*Gs(ie,5,5))/det &
     -(GI(ie,5,1)*GI(ie,3,3)*Gs(ie,5,5))/det

       Gs(ie,5,3)=(GI(ie,5,1)*GI(ie,1,3)*Gs(ie,5,5))/det &
    -(GI(ie,5,3)*GI(ie,1,1)*Gs(ie,5,5))/det

       Gs(ie,1,5)=(GI(ie,1,3)*GI(ie,3,5)*Gs(ie,5,5))/det &
    -(GI(ie,3,3)*GI(ie,1,5)*Gs(ie,5,5))/det

       Gs(ie,3,5)=(GI(ie,3,1)*GI(ie,1,5)*Gs(ie,5,5))/det &
    -(GI(ie,1,1)*GI(ie,3,5)*Gs(ie,5,5))/det

       Gs(ie,1,1)=GI(ie,3,3)/det &
    +(Gs(ie,1,5)*Gs(ie,5,1))/Gs(ie,5,5)
    
       Gs(ie,1,3)=-GI(ie,1,3)/det &
    +(Gs(ie,1,5)*Gs(ie,5,3))/Gs(ie,5,5)

       Gs(ie,3,1)=-GI(ie,3,1)/det &
    +(Gs(ie,3,5)*Gs(ie,5,1))/Gs(ie,5,5)
       Gs(ie,3,3)=GI(ie,1,1)/det &
    +(Gs(ie,3,5)*Gs(ie,5,3))/Gs(ie,5,5)

         Gs(ie,2,2)=Gs(ie,1,1)
         Gs(ie,2,4)=Gs(ie,1,3)
         Gs(ie,2,6)=Gs(ie,1,5)
         Gs(ie,4,2)=Gs(ie,3,1)
         Gs(ie,4,4)=Gs(ie,3,3)
         Gs(ie,4,6)=Gs(ie,3,5)
         Gs(ie,6,2)=Gs(ie,5,1)
         Gs(ie,6,4)=Gs(ie,5,3)
         Gs(ie,6,6)=Gs(ie,5,5)

         end do

           do ie=-N,N
             Gfscript(1,ie)=Gs(ie,1,1)
             Gfscript(2,ie)=Gs(ie,2,2)
             Gfscript(3,ie)=Gs(ie,3,3)
             Gfscript(4,ie)=Gs(ie,4,4)
             Gfscript(5,ie)=Gs(ie,5,5)
             Gfscript(6,ie)=Gs(ie,6,6)
            end do
       
             
            else if(Norbs==4)then
 
            do ie=jstart,jend
             do ik=1,NQPTS ! ksum
                 do jo=1,Norbs
                   do ko=1,Norbs
                  c(jo,ko)=-Ham0(ik,jo,ko)
                  end do
                  c(jo,jo)=XE(jo,ie)+c(jo,jo)
                end do

                 x=c(1,1)*c(3,3)-c(1,3)*c(3,1)
                  cINV(1,1)=c(3,3)/x
                  cINV(2,2)=cINV(1,1)
                  cINV(3,3)=c(1,1)/x
                  cINV(4,4)=cINV(3,3)
                  cINV(1,3)=-c(1,3)/x
                  cINV(3,1)=cINV(1,3)
                  cINV(2,4)=cINV(1,3)
                  cINV(4,2)=cINV(1,3)
                do io=1,Norbs
                 do jo=1,Norbs
           Gr1(ie,io,jo)=Gr1(ie,io,jo)+cINV(io,jo)/NQPTS
                 end do
               end do

               
               end do

               end do
!             write(6,*) '_____________'

      call MPI_BARRIER(MPI_COMM_WORLD,IERR)

      call MPI_ALLREDUCE(Gr1,Gr,(2*N+1)*Norbs*Norbs, &
     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)


            if(task_id==0)   write(6,*) '_____________'
              call cpu_time(t2)
             if(task_id==0)  write(6,*)'time',t2-t1

             do ie=-N,N
           Gf(1,ie)=Gr(ie,1,1)
           Gf(2,ie)=Gr(ie,2,2)
           Gf(3,ie)=Gr(ie,3,3)
           Gf(4,ie)=Gr(ie,4,4)
           end do

          deallocate(c)
          deallocate(cINV)
       	
          GI_new=CZERO
          G_new=CZERO
          G1_new=CZERO         
!          det=CZERO
          do ie=-N,N
              
              G_new(1,1)=Gr(ie,1,1)
              G_new(2,2)=Gr(ie,2,2)
              G_new(3,3)=Gr(ie,3,3)
              G_new(4,4)=Gr(ie,4,4)
              G_new(1,2)=Gr(ie,1,2)
              G_new(1,3)=Gr(ie,1,3)
              G_new(1,4)=Gr(ie,1,4)
              G_new(2,1)=Gr(ie,2,1)
              G_new(2,3)=Gr(ie,2,3)
              G_new(2,4)=Gr(ie,2,4)
              G_new(3,1)=Gr(ie,3,1)
              G_new(3,2)=Gr(ie,3,2)
              G_new(3,4)=Gr(ie,3,4)
              G_new(4,1)=Gr(ie,4,1)
              G_new(4,2)=Gr(ie,4,2)
              G_new(4,3)=Gr(ie,4,3)
            
           call cmatinv(G_new,GI_new,Norbs)

              G1_new(ie,1,1)=GI_new(1,1)
              G1_new(ie,2,2)=GI_new(2,2)
              G1_new(ie,3,3)=GI_new(3,3)
              G1_new(ie,4,4)=GI_new(4,4)
              G1_new(ie,1,2)=GI_new(1,2)
              G1_new(ie,1,3)=GI_new(1,3)
              G1_new(ie,1,4)=GI_new(1,4)
              G1_new(ie,2,1)=GI_new(2,1)
              G1_new(ie,2,3)=GI_new(2,3)
              G1_new(ie,2,4)=GI_new(2,4)
              G1_new(ie,3,1)=GI_new(3,1)
              G1_new(ie,3,2)=GI_new(3,2)
              G1_new(ie,3,4)=GI_new(3,4)
              G1_new(ie,4,1)=GI_new(4,1)
              G1_new(ie,4,2)=GI_new(4,2)
              G1_new(ie,4,3)=GI_new(4,3)

               
            end do

            deallocate(Gr)

           do ie=-N,N
             do io=1,Norbs
              do jo=1,Norbs
                if(io.eq.jo) then
           G1_new(ie,io,jo)=G1_new(ie,io,jo)+sigma(io,ie)+mu0-mu_c
                end if
               end do
              end do
           end do
 
           GS1=CZERO
           GS=CZERO
!           det=CZERO
           GI_new=CZERO
 
            do ie=-N,N

              GS1(1,1)=G1_new(ie,1,1)
              GS1(2,2)=G1_new(ie,2,2)
              GS1(3,3)=G1_new(ie,3,3)
              GS1(4,4)=G1_new(ie,4,4)
              GS1(1,2)=G1_new(ie,1,2)
              GS1(1,3)=G1_new(ie,1,3)
              GS1(1,4)=G1_new(ie,1,4)
              GS1(2,1)=G1_new(ie,2,1)
              GS1(2,3)=G1_new(ie,2,3)
              GS1(2,4)=G1_new(ie,2,4)
              GS1(3,1)=G1_new(ie,3,1)
              GS1(3,2)=G1_new(ie,3,2)
              GS1(3,4)=G1_new(ie,3,4)
              GS1(4,1)=G1_new(ie,4,1)
              GS1(4,2)=G1_new(ie,4,2)
              GS1(4,3)=G1_new(ie,4,3)

           call cmatinv(GS1,GI_new,Norbs)


              GS(ie,1,1)=GI_new(1,1)
              GS(ie,2,2)=GI_new(2,2)
              GS(ie,3,3)=GI_new(3,3)
              GS(ie,4,4)=GI_new(4,4)
              GS(ie,1,2)=GI_new(1,2)
              GS(ie,1,3)=GI_new(1,3)
              GS(ie,1,4)=GI_new(1,4)
              GS(ie,2,1)=GI_new(2,1)
              GS(ie,2,3)=GI_new(2,3)
              GS(ie,2,4)=GI_new(2,4)
              GS(ie,3,1)=GI_new(3,1)
              GS(ie,3,2)=GI_new(3,2)
              GS(ie,3,4)=GI_new(3,4)
              GS(ie,4,1)=GI_new(4,1)
              GS(ie,4,2)=GI_new(4,2)
              GS(ie,4,3)=GI_new(4,3)


             end do


             do ie=-N,N
              Gfscript(1,ie)=Gs(ie,1,1)
              Gfscript(2,ie)=Gs(ie,2,2)
              Gfscript(3,ie)=Gs(ie,3,3)
              Gfscript(4,ie)=Gs(ie,4,4)
             end do

             else if(Norbs==2)then
        
            do ie=jstart,jend
             do ik=1,NQPTS ! ksum
                 do jo=1,Norbs
                   do ko=1,Norbs
                  c(jo,ko)=-Ham0(ik,jo,ko)
                  end do
                  c(jo,jo)=XE(jo,ie)+c(jo,jo)
                end do

                 cINV(1,1)=1.d0/c(1,1) 
                 cINV(2,2)=1.d0/c(2,2)               
                 cINV(1,2)=0.d0
                 cINV(2,1)=0.d0

                do io=1,Norbs
                 do jo=1,Norbs
            Gr1(ie,io,jo)=Gr1(ie,io,jo)+cINV(io,jo)/NQPTS
                 end do
               end do

               end do

               end do
!             write(6,*) '_____________'

      call MPI_BARRIER(MPI_COMM_WORLD,IERR)

      call MPI_ALLREDUCE(Gr1,Gr,(2*N+1)*Norbs*Norbs, &
     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)


            if(task_id==0)   write(6,*) '_____________'
              call cpu_time(t2)
             if(task_id==0)  write(6,*)'time',t2-t1

             do ie=-N,N
           Gf(1,ie)=Gr(ie,1,1)
           Gf(2,ie)=Gr(ie,2,2)
           end do

          deallocate(c)
          deallocate(cINV)
              
          GI_new=CZERO
          G_new=CZERO
          G1_new=CZERO
!          det=CZERO
          do ie=-N,N

              G_new(1,1)=Gr(ie,1,1)
              G_new(2,2)=Gr(ie,2,2)
              G_new(1,2)=Gr(ie,1,2)
              G_new(2,1)=Gr(ie,2,1)
            
               
            call cmatinv(G_new,GI_new,Norbs)

              G1_new(ie,1,1)=GI_new(1,1)
              G1_new(ie,2,2)=GI_new(2,2)
              G1_new(ie,1,2)=GI_new(1,2)
              G1_new(ie,2,1)=GI_new(2,1)
              
            end do

            deallocate(Gr)

           do ie=-N,N
             do io=1,Norbs
              do jo=1,Norbs
                if(io.eq.jo) then
           G1_new(ie,io,jo)=G1_new(ie,io,jo)+sigma(io,ie)+mu0-mu_c
                end if
               end do
              end do
           end do

           GS1=CZERO
           GS=CZERO
!           det=CZERO
           GI_new=CZERO

            do ie=-N,N

              GS1(1,1)=G1_new(ie,1,1)
              GS1(2,2)=G1_new(ie,2,2)
              GS1(1,2)=G1_new(ie,1,2)
              GS1(2,1)=G1_new(ie,2,1)
           
            call cmatinv(GS1,GI_new,Norbs)


              GS(ie,1,1)=GI_new(1,1)
              GS(ie,2,2)=GI_new(2,2)
              GS(ie,1,2)=GI_new(1,2)
              GS(ie,2,1)=GI_new(2,1)

              
             end do


             do ie=-N,N
              Gfscript(1,ie)=Gs(ie,1,1)
              Gfscript(2,ie)=Gs(ie,2,2)
             end do


           end if   

   else   ! if LDA_DMFT =0, Real lattice calculations
 
         if(Norbs==6)then
          do ie=jstart,jend
             do ik=1,N_PTS ! ksum
                 do jo=1,Norbs
                   do ko=1,Norbs
                  c(jo,ko)=-Ham_model(ik,jo,ko)
                  end do
                  c(jo,jo)=XE(jo,ie)+c(jo,jo)
                end do
        
             x=c(1,1)*c(3,3)-c(1,3)*c(3,1)

            cINV(5,5)=c(5,5)-(c(5,1)*c(3,3)*c(1,5))/x &
        -(c(5,3)*c(1,1)*c(3,5))/x+(c(5,1)*c(1,3)*c(3,5))/x &
        +(c(5,3)*c(1,5)*c(3,1))/x

           cINV(5,5)=1.d0/cINV(5,5)

           cINV(5,1)=(c(3,1)*c(5,3)*cINV(5,5))/x &
       -(c(5,1)*c(3,3)*cINV(5,5))/x

           cINV(5,3)=(c(5,1)*c(1,3)*cINV(5,5))/x &
       -(c(5,3)*c(1,1)*cINV(5,5))/x

          cINV(1,5)=(c(1,3)*c(3,5)*cINV(5,5))/x &
      -(c(3,3)*c(1,5)*cINV(5,5))/x

           cINV(3,5)=(c(3,1)*c(1,5)*cINV(5,5))/x &
      -(c(1,1)*c(3,5)*cINV(5,5))/x

         cINV(1,1)=c(3,3)/x+(cINV(1,5)*cINV(5,1))/cINV(5,5)
         cINV(1,3)=-c(1,3)/x+(cINV(1,5)*cINV(5,3))/cINV(5,5)                        
         cINV(3,1)=-c(3,1)/x+(cINV(3,5)*cINV(5,1))/cINV(5,5)
         cINV(3,3)=c(1,1)/x+(cINV(3,5)*cINV(5,3))/cINV(5,5)

         cINV(2,2)=cINV(1,1)
         cINV(2,4)=cINV(1,3)
         cINV(2,6)=cINV(1,5)
         cINV(4,2)=cINV(3,1)
         cINV(4,4)=cINV(3,3)
         cINV(4,6)=cINV(3,5)
         cINV(6,2)=cINV(5,1)
         cINV(6,4)=cINV(5,3)
         cINV(6,6)=cINV(5,5)


               do lo=1,Norbs
                do mo=1,Norbs
        Gr1(ie,lo,mo)=Gr1(ie,lo,mo)+(cINV(lo,mo)/dfloat(N_PTS))
                 end do
                  end do

                end do
             end do


      call MPI_BARRIER(MPI_COMM_WORLD,IERR)

      call MPI_ALLREDUCE(Gr1,Gr,(2*N+1)*Norbs*Norbs, &
     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)


            if(task_id==0)   write(6,*) '_____________'
              call cpu_time(t2)
             if(task_id==0)  write(6,*)'time',t2-t1

          do ie=-N,N
           Gf(1,ie)=Gr(ie,1,1)
           Gf(2,ie)=Gr(ie,2,2)
           Gf(3,ie)=Gr(ie,3,3)
           Gf(4,ie)=Gr(ie,4,4)
           Gf(5,ie)=Gr(ie,5,5)
           Gf(6,ie)=Gr(ie,6,6)
          end do


          deallocate(c)
          deallocate(cINV)



          GI=CZERO
          det=CZERO
          do ie=-N,N

        det=Gr(ie,1,1)*Gr(ie,3,3)-Gr(ie,1,3)*Gr(ie,3,1)

        GI(ie,5,5)=Gr(ie,5,5) &
     -(Gr(ie,5,1)*Gr(ie,3,3)*Gr(ie,1,5))/det &
     -(Gr(ie,5,3)*Gr(ie,1,1)*Gr(ie,3,5))/det &
     +(Gr(ie,5,1)*Gr(ie,1,3)*Gr(ie,3,5))/det &
     +(Gr(ie,5,3)*Gr(ie,1,5)*Gr(ie,3,1))/det

        GI(ie,5,5)=1.d0/GI(ie,5,5)
!        write(6,*)GI(ie,5,5),Gr(ie,5,5)


       GI(ie,5,1)=(Gr(ie,3,1)*Gr(ie,5,3)*GI(ie,5,5))/det &
      -(Gr(ie,5,1)*Gr(ie,3,3)*GI(ie,5,5))/det

       GI(ie,5,3)=(Gr(ie,5,1)*Gr(ie,1,3)*GI(ie,5,5))/det &
    -(Gr(ie,5,3)*Gr(ie,1,1)*GI(ie,5,5))/det

       GI(ie,1,5)=(Gr(ie,1,3)*Gr(ie,3,5)*GI(ie,5,5))/det &
    -(Gr(ie,3,3)*Gr(ie,1,5)*GI(ie,5,5))/det

       GI(ie,3,5)=(Gr(ie,3,1)*Gr(ie,1,5)*GI(ie,5,5))/det &
    -(Gr(ie,1,1)*Gr(ie,3,5)*GI(ie,5,5))/det

       GI(ie,1,1)=Gr(ie,3,3)/det &
    +(GI(ie,1,5)*GI(ie,5,1))/GI(ie,5,5)
       GI(ie,1,3)=-Gr(ie,1,3)/det &
    +(GI(ie,1,5)*GI(ie,5,3))/GI(ie,5,5)
       GI(ie,3,1)=-Gr(ie,3,1)/det &
    +(GI(ie,3,5)*GI(ie,5,1))/GI(ie,5,5)
       GI(ie,3,3)=Gr(ie,1,1)/det &
    +(GI(ie,3,5)*GI(ie,5,3))/GI(ie,5,5)

         GI(ie,2,2)=GI(ie,1,1)
         GI(ie,2,4)=GI(ie,1,3)
         GI(ie,2,6)=GI(ie,1,5)
         GI(ie,4,2)=GI(ie,3,1)
         GI(ie,4,4)=GI(ie,3,3)
         GI(ie,4,6)=GI(ie,3,5)
         GI(ie,6,2)=GI(ie,5,1)
         GI(ie,6,4)=GI(ie,5,3)
         GI(ie,6,6)=GI(ie,5,5)

         end do

          deallocate(Gr)

          do ie=-N,N
            do io=1,Norbs
              do jo=1,Norbs
                if(io.eq.jo) then
                GI(ie,io,jo)=GI(ie,io,jo)+sigma(io,ie)+mu0-mu_c
                end if
               end do
              end do
           end do

          Gs=CZERO
          det=CZERO
          do ie=-N,N

         det=GI(ie,1,1)*GI(ie,3,3)-GI(ie,1,3)*GI(ie,3,1)

         Gs(ie,5,5)=GI(ie,5,5) &
     -(GI(ie,5,1)*GI(ie,3,3)*GI(ie,1,5))/det &
     -(GI(ie,5,3)*GI(ie,1,1)*GI(ie,3,5))/det &
     +(GI(ie,5,1)*GI(ie,1,3)*GI(ie,3,5))/det &
    +(GI(ie,5,3)*GI(ie,1,5)*GI(ie,3,1))/det

         Gs(ie,5,5)=1.d0/Gs(ie,5,5)

       Gs(ie,5,1)=(GI(ie,3,1)*GI(ie,5,3)*Gs(ie,5,5))/det &
     -(GI(ie,5,1)*GI(ie,3,3)*Gs(ie,5,5))/det

       Gs(ie,5,3)=(GI(ie,5,1)*GI(ie,1,3)*Gs(ie,5,5))/det &
    -(GI(ie,5,3)*GI(ie,1,1)*Gs(ie,5,5))/det

       Gs(ie,1,5)=(GI(ie,1,3)*GI(ie,3,5)*Gs(ie,5,5))/det &
    -(GI(ie,3,3)*GI(ie,1,5)*Gs(ie,5,5))/det

       Gs(ie,3,5)=(GI(ie,3,1)*GI(ie,1,5)*Gs(ie,5,5))/det &
    -(GI(ie,1,1)*GI(ie,3,5)*Gs(ie,5,5))/det

       Gs(ie,1,1)=GI(ie,3,3)/det &
    +(Gs(ie,1,5)*Gs(ie,5,1))/Gs(ie,5,5)

       Gs(ie,1,3)=-GI(ie,1,3)/det &
    +(Gs(ie,1,5)*Gs(ie,5,3))/Gs(ie,5,5)

       Gs(ie,3,1)=-GI(ie,3,1)/det &
    +(Gs(ie,3,5)*Gs(ie,5,1))/Gs(ie,5,5)
       Gs(ie,3,3)=GI(ie,1,1)/det &
    +(Gs(ie,3,5)*Gs(ie,5,3))/Gs(ie,5,5)

         Gs(ie,2,2)=Gs(ie,1,1)
         Gs(ie,2,4)=Gs(ie,1,3)
         Gs(ie,2,6)=Gs(ie,1,5)
         Gs(ie,4,2)=Gs(ie,3,1)
         Gs(ie,4,4)=Gs(ie,3,3)
         Gs(ie,4,6)=Gs(ie,3,5)
         Gs(ie,6,2)=Gs(ie,5,1)
         Gs(ie,6,4)=Gs(ie,5,3)
         Gs(ie,6,6)=Gs(ie,5,5)

         end do

           do ie=-N,N
             Gfscript(1,ie)=Gs(ie,1,1)
             Gfscript(2,ie)=Gs(ie,2,2)
             Gfscript(3,ie)=Gs(ie,3,3)
             Gfscript(4,ie)=Gs(ie,4,4)
             Gfscript(5,ie)=Gs(ie,5,5)
             Gfscript(6,ie)=Gs(ie,6,6)
            end do

    

            else if(Norbs==4)then

            do ie=jstart,jend
             do ik=1,N_PTS ! ksum
                 do jo=1,Norbs
                   do ko=1,Norbs
                  c(jo,ko)=-Ham_model(ik,jo,ko)
                  end do
                  c(jo,jo)=XE(jo,ie)+c(jo,jo)
                end do

                 x=c(1,1)*c(3,3)-c(1,3)*c(3,1)
                  cINV(1,1)=c(3,3)/x
                  cINV(2,2)=cINV(1,1)
                  cINV(3,3)=c(1,1)/x
                  cINV(4,4)=cINV(3,3)
                  cINV(1,3)=-c(1,3)/x
                  cINV(3,1)=cINV(1,3)
                  cINV(2,4)=cINV(1,3)
                  cINV(4,2)=cINV(1,3)
             
               do io=1,Norbs
                 do jo=1,Norbs
           Gr1(ie,io,jo)=Gr1(ie,io,jo)+cINV(io,jo)/N_PTS
                 end do
               end do
                
               end do

               end do
!             write(6,*) '_____________'

      call MPI_BARRIER(MPI_COMM_WORLD,IERR)

      call MPI_ALLREDUCE(Gr1,Gr,(2*N+1)*Norbs*Norbs, &
     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)


            if(task_id==0)   write(6,*) '_____________'
              call cpu_time(t2)
             if(task_id==0)  write(6,*)'time',t2-t1

             do ie=-N,N
           Gf(1,ie)=Gr(ie,1,1)
           Gf(2,ie)=Gr(ie,2,2)
           Gf(3,ie)=Gr(ie,3,3)
           Gf(4,ie)=Gr(ie,4,4)
           end do

          deallocate(c)
          deallocate(cINV)
        
          GI_new=CZERO
          G_new=CZERO
          G1_new=CZERO
!          det=CZERO
          do ie=-N,N

              G_new(1,1)=Gr(ie,1,1)
              G_new(2,2)=Gr(ie,2,2)
              G_new(3,3)=Gr(ie,3,3)
              G_new(4,4)=Gr(ie,4,4)
              G_new(1,2)=Gr(ie,1,2)
              G_new(1,3)=Gr(ie,1,3)
              G_new(1,4)=Gr(ie,1,4)
              G_new(2,1)=Gr(ie,2,1)
              G_new(2,3)=Gr(ie,2,3)
              G_new(2,4)=Gr(ie,2,4)
              G_new(3,1)=Gr(ie,3,1)
              G_new(3,2)=Gr(ie,3,2)
              G_new(3,4)=Gr(ie,3,4)
              G_new(4,1)=Gr(ie,4,1)
              G_new(4,2)=Gr(ie,4,2)
              G_new(4,3)=Gr(ie,4,3)
             
                call cmatinv(G_new,GI_new,Norbs)

              G1_new(ie,1,1)=GI_new(1,1)
              G1_new(ie,2,2)=GI_new(2,2)
              G1_new(ie,3,3)=GI_new(3,3)
              G1_new(ie,4,4)=GI_new(4,4)
              G1_new(ie,1,2)=GI_new(1,2)
              G1_new(ie,1,3)=GI_new(1,3)
              G1_new(ie,1,4)=GI_new(1,4)
              G1_new(ie,2,1)=GI_new(2,1)
              G1_new(ie,2,3)=GI_new(2,3)
              G1_new(ie,2,4)=GI_new(2,4)
              G1_new(ie,3,1)=GI_new(3,1)
              G1_new(ie,3,2)=GI_new(3,2)
              G1_new(ie,3,4)=GI_new(3,4)
              G1_new(ie,4,1)=GI_new(4,1)
              G1_new(ie,4,2)=GI_new(4,2)
              G1_new(ie,4,3)=GI_new(4,3)


            end do

            deallocate(Gr)

           do ie=-N,N
             do io=1,Norbs
              do jo=1,Norbs
                if(io.eq.jo) then
           G1_new(ie,io,jo)=G1_new(ie,io,jo)+sigma(io,ie)+mu0-mu_c
                end if
               end do
              end do
           end do

           GS1=CZERO
           GS=CZERO
!           det=CZERO
           GI_new=CZERO

            do ie=-N,N

              GS1(1,1)=G1_new(ie,1,1)
              GS1(2,2)=G1_new(ie,2,2)
              GS1(3,3)=G1_new(ie,3,3)
              GS1(4,4)=G1_new(ie,4,4)
              GS1(1,2)=G1_new(ie,1,2)
              GS1(1,3)=G1_new(ie,1,3)
              GS1(1,4)=G1_new(ie,1,4)
              GS1(2,1)=G1_new(ie,2,1)
              GS1(2,3)=G1_new(ie,2,3)
              GS1(2,3)=G1_new(ie,2,3)
              GS1(2,4)=G1_new(ie,2,4)
              GS1(3,1)=G1_new(ie,3,1)
              GS1(3,2)=G1_new(ie,3,2)
              GS1(3,4)=G1_new(ie,3,4)
              GS1(4,1)=G1_new(ie,4,1)
              GS1(4,2)=G1_new(ie,4,2)
              GS1(4,3)=G1_new(ie,4,3)

           call cmatinv(GS1,GI_new,Norbs)


              GS(ie,1,1)=GI_new(1,1)
              GS(ie,2,2)=GI_new(2,2)
              GS(ie,3,3)=GI_new(3,3)
              GS(ie,4,4)=GI_new(4,4)
              GS(ie,1,2)=GI_new(1,2)
              GS(ie,1,3)=GI_new(1,3)
              GS(ie,1,4)=GI_new(1,4)
              GS(ie,2,1)=GI_new(2,1)
              GS(ie,2,3)=GI_new(2,3)
              GS(ie,2,4)=GI_new(2,4)
              GS(ie,3,1)=GI_new(3,1)
              GS(ie,3,2)=GI_new(3,2)
              GS(ie,3,4)=GI_new(3,4)
              GS(ie,4,1)=GI_new(4,1)
              GS(ie,4,2)=GI_new(4,2)
              GS(ie,4,3)=GI_new(4,3)


             end do


             do ie=-N,N
              Gfscript(1,ie)=Gs(ie,1,1)
              Gfscript(2,ie)=Gs(ie,2,2)
              Gfscript(3,ie)=Gs(ie,3,3)
              Gfscript(4,ie)=Gs(ie,4,4)
             end do


           else if(Norbs==2)then

            do ie=jstart,jend
             do ik=1,N_PTS ! ksum
                 do jo=1,Norbs
                   do ko=1,Norbs
                  c(jo,ko)=-Ham_model(ik,jo,ko)
                  end do
                  c(jo,jo)=XE(jo,ie)+c(jo,jo)
                end do

                 cINV(1,1)=1.d0/c(1,1)
                 cINV(2,2)=1.d0/c(2,2)
                 cINV(1,2)=0.d0
                 cINV(2,1)=0.d0

                do io=1,Norbs
                 do jo=1,Norbs
            Gr1(ie,io,jo)=Gr1(ie,io,jo)+cINV(io,jo)/N_PTS
                 end do
               end do

               end do

               end do
!             write(6,*) '_____________'
                              
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)

       call MPI_ALLREDUCE(Gr1,Gr,(2*N+1)*Norbs*Norbs, &
     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)


            if(task_id==0)   write(6,*) '_____________'
              call cpu_time(t2)
             if(task_id==0)  write(6,*)'time',t2-t1

             do ie=-N,N
           Gf(1,ie)=Gr(ie,1,1)
           Gf(2,ie)=Gr(ie,2,2)
           end do

          deallocate(c)
          deallocate(cINV)

          GI_new=CZERO
          G_new=CZERO
          G1_new=CZERO
!          det=CZERO
          do ie=-N,N

              G_new(1,1)=Gr(ie,1,1)
              G_new(2,2)=Gr(ie,2,2)
              G_new(1,2)=Gr(ie,1,2)
              G_new(2,1)=Gr(ie,2,1)


            call cmatinv(G_new,GI_new,Norbs)

              G1_new(ie,1,1)=GI_new(1,1)
              G1_new(ie,2,2)=GI_new(2,2)
              G1_new(ie,1,2)=GI_new(1,2)
              G1_new(ie,2,1)=GI_new(2,1)

            end do

            deallocate(Gr)

           do ie=-N,N
             do io=1,Norbs
              do jo=1,Norbs
                if(io.eq.jo) then
           G1_new(ie,io,jo)=G1_new(ie,io,jo)+sigma(io,ie)+mu0-mu_c
                end if
               end do
            end do
          end do 
              
           GS1=CZERO
           GS=CZERO
           GI_new=CZERO

            do ie=-N,N

              GS1(1,1)=G1_new(ie,1,1)
              GS1(2,2)=G1_new(ie,2,2)
              GS1(1,2)=G1_new(ie,1,2)
              GS1(2,1)=G1_new(ie,2,1)

            call cmatinv(GS1,GI_new,Norbs)


              GS(ie,1,1)=GI_new(1,1)
              GS(ie,2,2)=GI_new(2,2)
              GS(ie,1,2)=GI_new(1,2)
              GS(ie,2,1)=GI_new(2,1)


             end do


             do ie=-N,N
              Gfscript(1,ie)=Gs(ie,1,1)
              Gfscript(2,ie)=Gs(ie,2,2)
             end do

          end if

           end if


                    
        deallocate(GS)
        deallocate(GI)
        deallocate(G_new)
        deallocate(G1_new)
        deallocate(GS1)
        deallocate(GI_new)
        deallocate(XE)


	return
	end
!#####################################################################
