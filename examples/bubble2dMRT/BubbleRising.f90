      program MRT_LBM_MPF2D
!----------------------------------------------------------------------c
!     A  MRT LBM for multiphase flow problems
!----------------------------------------------------------------------c
!-------vars for LBM
      parameter(nq=9)
      integer, dimension (nq)::ex,ey,ieopp
      double precision, dimension (nq)::wt
!-------vars for Euler mesh
      integer nx,ny
      parameter(nx=121, ny=241)
      double precision xl ,yl 
      double precision, dimension (0:nx+1)::x
      double precision, dimension (0:ny+1)::y
      double precision dh
!-------vars for hydro on Euler mesh------------------------------------c
      double precision, dimension (0:nx+1,0:ny+1):: ux,uy,ux_pre,uy_pre,rho,pressure,visfluid
      double precision, dimension (0:nx+1,0:ny+1)::phi,phi_pre,phi_0,chpoten,nablaphix,nablaphiy,mchpoten
	  double precision, dimension (0:nx+1,0:ny+1)::psi,concentration,nablapsix,nablapsiy,nablaconx,nablacony,curv
	  double precision,dimension (0:nx+1,0:ny+1)::bodyforcex,bodyforcey,surtenforcex,surtenforcey,interforcex,interforcey,preforcex,preforcey,visforcex,visforcey
      double precision, dimension (0:nx+1,0:ny+1)::tau_ddf
!-------distribution function
      double precision, dimension (nq,0:nx+1,0:ny+1):: ddf,eddf
      double precision, dimension (nq,0:nx+1,0:ny+1):: ddg,eddg
	  double precision, dimension (nq,0:nx+1,0:ny+1):: ddh,eddh
!-------others parms
      double precision dt,tlim,uc,rho_low,rho_high,phi_low,phi_high,vis_low,vis_high,nu,mobility,surdiff,intthick,surtencoefficent,beta,kappa
	  double precision tau_ddg,tau_ddh,tau_ddf_low,tau_ddf_high
      integer istep
      double precision ctime, ptime, err, L2err,Loo,maxphi,minphi,maxmaxphi,minminphi,L2phi,mass1,mass2,mass3,HeavisideValue,phasevolume1,phasevolume2,phasevolume3,Euv2,Euvmax
!-------work array
      double precision, dimension (nq,0:nx+1,0:ny+1):: work3D
!----------------------------------------------------------------------c

!-------set LBM parms
      call setLBMParms2D(nq,ex,ey,ieopp,wt)
!-------generate Euler mesh
      call generateGrid2D(xl, yl, nx, ny, x, y, dh)
!-------set other parms
      call setOtherParms2D(xl,dh,dt,tlim,uc,rho_low,rho_high,phi_low,phi_high,vis_low,vis_high,nu,mobility,intthick,surdiff,surtencoefficent,beta,kappa,tau_ddg,tau_ddh,tau_ddf_low,tau_ddf_high)
!-------initialize for hydro macro vars
      call initHydroMacroVars2D(nx,ny,nq,x,y,tau_ddf,ux,uy,ux_pre,uy_pre,rho,pressure,phi,phi_pre,chpoten,mchpoten,&
      psi,concentration,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,bodyforcex,bodyforcey,surtenforcex,surtenforcey,interforcex,interforcey,preforcex,preforcey,&
      beta,kappa,phi_low,phi_high,intthick,dh,tau_ddf_low,tau_ddf_high,rho_low,rho_high,vis_low,vis_high,visfluid)
 
    open(2,file='x.plt',access='append') 
         WRITE(2,'(9F10.3)') ( x(i), i=1,nx )    
    close(2)        
      
!-------initialize for density
!------------------------------------------------------------------c
      
!-------initialize for Order parameter 
     call computeOrderEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,phi,chpoten,mchpoten,nablaphix,nablaphiy,mobility,intthick,eddg) 
     call computeOrderEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,phi,chpoten,mchpoten,nablaphix,nablaphiy,mobility,intthick,ddg)
     
     call computePsiEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,psi,eddh) 
     call computePsiEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,psi,ddh) 
     
     call computeDenEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,rho,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,eddf,dt) 
     call computeDenEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,rho,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,ddf,dt) 
     
!-------save results
      istep = 0
      ctime = istep*dt
      ptime = ctime*uc  
      
      call saveFlowField2D(ptime,ctime,uc,nx,ny,x,y,ux,uy,rho,pressure,phi,psi,concentration,phi_high,phi_low,chpoten,dh,surtenforcex,surtenforcey,preforcex,preforcey,bodyforcex,bodyforcey,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,curv,intthick,visforcex,visforcey)
!-----------------------------------------------------

       !do while(ptime<tlim)
		do while(istep.le.1000000)
         istep = istep + 1
         ctime = istep*dt  
         ptime = ctime*uc
       
      write(*,*)'-----------------------------------------------------istep=',istep
         
       call computeForce(nx,ny,nq,ex,ey,wt,chpoten,phi,rho,nablaphix,nablaphiy,nablaconx,nablacony,curv,concentration,bodyforcex,bodyforcey,surtenforcex,surtenforcey,interforcex,interforcey,preforcex,preforcey,dh,rho_high,rho_low,surtencoefficent,intthick)
       call collisionOrderDF2D(nq,ex,ey,wt,nx,ny,ux,ux_pre,uy,uy_pre,phi,phi_pre,chpoten,tau_ddg,dt,eddg,ddg,nablaphix,nablaphiy,mobility,intthick)
       !call collisionPsiDF2D(nq,ex,ey,wt,nx,ny,phi,psi,dt,eddh,ddh,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,intthick,tau_ddh)
       call collisionDenDF2D(nq,nx,ny,ex,ey,wt,ux,uy,rho,visfluid,pressure,tau_ddf,dt,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,visforcex,visforcey,nablaphix,nablaphiy,eddf,ddf,rho_low,rho_high)
       
       call streamOrderDF2D(nq,ex,ey,wt,nx,ny,ddg,work3D)
       !call streamPsiDF2D(nq,ex,ey,wt,nx,ny,ddh,work3D)
       call streamDenDF2D(nq,ex,ey,wt,nx,ny,ddf,work3D)
       
       !!call setOrderfBC2D(nq,nx,ny,ddg,eddg)
       !!call setPsifBC2D(nq,nx,ny,ddh,eddh) 

       
       call computeOrderMacro2D(nq,nx,ny,ex,ey,wt,ddg,eddg,ddh,phi,phi_pre,psi,concentration,chpoten,mchpoten,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,beta,kappa,dh,phi_low,phi_high,rho,ux,ux_pre,uy,uy_pre,istep,dt,x,y,tau_ddf,tau_ddf_low,tau_ddf_high,mobility,intthick,surtencoefficent,rho_low,rho_high,vis_low,vis_high,visfluid)
         
       !call explotation(nx,ny,nq,x,y,phi,psi,concentration,nablaphix,nablaphiy,dh,intthick,dt,ex,ey,wt,nablaconx,nablacony)
       
       call computeMacro2D(nq,ex,ey,wt,nx,ny,ddf,ux,uy,rho,pressure,nablaphix,nablaphiy,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,visforcex,visforcey,dt,rho_low,rho_high)
       
       !!call setMacroBC2D(uc,nx,ny,ux,uy,rho,pressure)
       call setMacroOrderBC2D(nx,ny,phi,phi_pre,psi,concentration,chpoten,phi_high,phi_low,dh)
       
        !if (mod(istep,100).eq.0) then
        if (mod(istep,100).eq.0) then
            call saveFlowField2D(ptime,ctime,uc,nx,ny,x,y,ux,uy,rho,pressure,phi,psi,concentration,phi_high,phi_low,chpoten,dh,surtenforcex,surtenforcey,preforcex,preforcey,bodyforcex,bodyforcey,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,curv,intthick,visforcex,visforcey)
        end if

 !        if(err.le.1e-9) stop

       call computeOrderEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,phi,chpoten,mchpoten,nablaphix,nablaphiy,mobility,intthick,eddg) 
       !call computePsiEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,psi,eddh) 
       call computeDenEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,rho,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,eddf,dt) 
       
      if (mod(istep,2000).eq.0) then
!      L2phi=0.d0
!      L2err=0.d0
!      maxphi=0.d0
!      do i=1,nx
!          do j=1,ny
!             phi_0(i,j,1)=0.50*(phi_low+phi_high)+0.50*(phi_high-phi_low)*tanh(2.d0*(15.0*1.0-sqrt((x(i)-50.0*1.0)**2+(y(j)-50.0*1.0)**2))/intthick)
!             phi_0(i,j,2)=0.50*(phi_low+phi_high)-0.50*(phi_high-phi_low)*tanh(2.d0*(30.0*1.0-sqrt((x(i)-50.0*1.0)**2+(y(j)-50.0*1.0)**2))/intthick)
!             phi_0(i,j,3)=1.0-phi(i,j,1)-phi(i,j,2)
!             phi_0(i,j,4)=0.0
!            
!             L2phi=L2phi+phi(i,j,3)**2
!             L2err=L2err+(phi(i,j,3)-phi_0(i,j,3))**2 
!             if(maxphi.lt.abs(phi(i,j,3))) maxphi=abs(phi(i,j,3))
!          end do          
!      end do
!    
222 format(f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14)    
    
!    open(2,file='err.plt',access='append')
!    write(2,*) ptime,sqrt(L2phi/real(nx)**2),maxphi
!    close(2)
!    
!    mass1=0.0
!    mass2=0.0
!    mass3=0.0
!      do i=1,nx
!          do j=1,ny
!            mass1=mass1+phi(i,j,4)*dh*dh
!            mass2=mass2+phi(i,j,2)*dh*dh
!            mass3=mass3+phi(i,j,3)*dh*dh
!          end do
!      end do
!    open(2,file='mass.plt',access='append')
!    write(2,222) ptime,mass1!,mass2,mass3
!    close(2)
!    
        phasevolume1=0.d0
        phasevolume2=0.d0
        phasevolume3=0.d0
        do i=1,nx
            do j=1,ny
                phasevolume1=phasevolume1+HeavisideValue(5.0*(phi(i,j)-0.5))*dh*dh
                !phasevolume2=phasevolume2+HeavisideValue(5.0*(phi(i,j,2)-0.5))*dh*dh
                !phasevolume3=phasevolume3+HeavisideValue(5.0*(phi(i,j,3)-0.5))*dh*dh
            end do                 
        end do  
            open(2,file='phasevolume.plt',access='append')
            write(2,222) ptime,phasevolume1!,phasevolume2,phasevolume3
            close(2)           
!        Euv2=0.0
!        Euvmax=0.0
!        do i=1,nx
!            do j=1,ny
!                Euv2=Euv2+ux(i,j)**2+uy(i,j)**2
!                if(Euvmax.le.sqrt(ux(i,j)**2+uy(i,j)**2)) Euvmax=sqrt(ux(i,j)**2+uy(i,j)**2)
!            end do
!        end do
!            open(2,file='Euv2.plt',access='append')
!            write(2,222) ptime,Euvmax,sqrt(Euv2/real(nx)/real(ny))
!            close(2)              
!            
      end if       


      end do

      end
!======================================================================c      
!-------set LBM parms
!----------------------------------------------------------------------c
      subroutine setLBMParms2D(nq,ex,ey,ieopp,wt)
!-------input parms
      integer nq
!-------output parms
      integer ex(nq), ey(nq), ieopp(nq)
      double precision wt(nq)
!----------------------------------------------------------------------c
      integer iq
      ex(1)= 1
      ey(1)= 0
      ex(2)= 0
      ey(2)= 1
      ex(3)=-1
      ey(3)= 0
      ex(4)= 0
      ey(4)=-1
      ex(5)= 1
      ey(5)= 1
      ex(6)=-1
      ey(6)= 1
      ex(7)=-1
      ey(7)=-1
      ex(8)= 1
      ey(8)=-1
      ex(9)= 0
      ey(9)= 0
!---------------------------------
      ieopp(1)=3
      ieopp(2)=4
      ieopp(3)=1
      ieopp(4)=2
      ieopp(5)=7
      ieopp(6)=8
      ieopp(7)=5
      ieopp(8)=6
      ieopp(9)=9
!---------------------------------
      do iq=1,nq
         if(iq.le.4) then
            wt(iq)=1.0/9.d0
         elseif(iq.le.8) then
            wt(iq)=1.0/36.d0
         else
            wt(iq)=4.0/9.d0
         endif
      end do

      end subroutine      
!======================================================================c
!-------generate Euler mesh
!----------------------------------------------------------------------c
      subroutine generateGrid2D(xl, yl, nx, ny, x, y, dh)
!-------input parms
      double precision xl, yl
      integer nx, ny
!-------output parms
      double precision x(0:nx+1), y(0:ny+1)
      double precision dh
!----------------------------------------------------------------------c
      integer i, j
      double precision dx, dy
      
      xl = real(nx-1)
      yl = real(ny-1)
      dx = xl/real(nx-1)
      dy = dx
      dh = dx
      do i=0,nx+1
         x(i)=(i-1)*dx
      end do
      do j=0,ny+1
         y(j)=(j-1)*dy-0.5
      end do
      
      end subroutine
!======================================================================c
!-------set other parms
!----------------------------------------------------------------------c
      subroutine setOtherParms2D(xl,dh,dt,tlim,uc,rho_low,rho_high,phi_low,phi_high,vis_low,vis_high,nu,mobility,intthick,surdiff,surtencoefficent,beta,kappa,tau_ddg,tau_ddh,tau_ddf_low,tau_ddf_high)
!-------input parms
      double precision xl,dh,dt,tlim,uc
!-------output parms
      double precision rho_low,rho_high,phi_low,phi_high,vis_low,vis_high,nu,mobility,intthick,surtencoefficent,beta,kappa,surdiff
	  double precision tau_ddg,tau_ddh,tau_ddf_low,tau_ddf_high
!----------------------------------------------------------------------c
    
      tlim = 2000.0

      dt = dh
      uc = 0.01d0  
      rho_low = 0.001d0
      rho_high= 1.d0
      phi_low = 0.d0
      phi_high= 1.d0
      mobility=0.1/10.0*1.0
	  surdiff=0.08
      intthick=4.d0*dh
      surtencoefficent=1.e-4*60.0/125.0
      beta=12.d0*surtencoefficent/intthick/(phi_high-phi_low)**4 
      kappa=3.d0*intthick*surtencoefficent/2.d0/(phi_high-phi_low)**2
      vis_low=0.6/3500.0
      vis_high=0.6/35.0  
       
      tau_ddg=3.d0*mobility/dt+0.5
	  tau_ddh=3.d0*surdiff/dt+0.5
      tau_ddf_low=3.d0*vis_low/dt+0.5d0
      tau_ddf_high=3.d0*vis_high/dt+0.5d0
      
      write(*,*) 'beta=',beta,'kappa=',kappa,'tau_ddg=',tau_ddg
      end subroutine
!======================================================================c
!-------initialize hydro vars
!----------------------------------------------------------------------c
      subroutine initHydroMacroVars2D(nx,ny,nq,x,y,tau_ddf,ux,uy,ux_pre,uy_pre,rho,pressure,phi,phi_pre,chpoten,mchpoten,psi,concentration,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,bodyforcex,bodyforcey,surtenforcex,surtenforcey,interforcex,interforcey,preforcex,preforcey,beta,kappa,phi_low,phi_high,intthick,dh,tau_ddf_low,tau_ddf_high,rho_low,rho_high,vis_low,vis_high,visfluid)
!-------input parms
      integer nx,ny,nq
      double precision beta,kappa,phi_low,phi_high,intthick,pi,dh,tau_ddf_low,tau_ddf_high,rho_low,rho_high,vis_low,vis_high
	  double precision x(0:nx+1),y(0:ny+1)
!-------output parms
      double precision, dimension(0:nx+1,0:ny+1)::tau_ddf,ux,uy,ux_pre,uy_pre,rho,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,interforcex,interforcey,preforcex,preforcey,visfluid,psi,concentration,tempphi,nablapsix,nablapsiy,nablaconx,nablacony
      double precision, dimension(0:nx+1,0:ny+1)::phi,phi_pre,chpoten,nablaphix,nablaphiy,mchpoten
!----------------------------------------------------------------------c
      integer i, j, k, iq

      pi=4.d0*atan(1.d0)
      
          do i=1,nx
              do j=1,ny
                  phi(i,j)=0.50-0.50*tanh(2.d0*(30.0-sqrt((x(i)-60.0)**2+(y(j)-60.0)**2))/intthick)
                  !phi(i,j)=0.50-0.50*tanh(2.d0*(60.0-y(j))/intthick)
                  ux(i,j)=0.d0
                  uy(i,j)=0.d0
              end do
          end do      
        
          do i=1,nx
              phi(i,0)=phi(i,ny)
              phi(i,ny+1)=phi(i,1)
              ux(i,0)=ux(i,ny)
              ux(i,ny+1)=ux(i,1)    
              uy(i,0)=uy(i,ny)
              uy(i,ny+1)=uy(i,1)		  
          end do      
      
          do j=1,ny
              phi(0,j)=phi(nx,j)
              phi(0,j)=phi(nx,j)
              phi(0,j)=phi(nx,j)
              phi(0,j)=phi(nx,j)
              phi(0,j)=phi(nx,j)
              phi(nx+1,j)=phi(1,j)
              ux(0,j)=ux(nx,j)
              ux(nx+1,j)=ux(1,j)
              uy(0,j)=uy(nx,j)
              uy(nx+1,j)=uy(1,j)          
          end do
      
          phi(0,0)=phi(nx,ny)
          phi(nx+1,0)=phi(1,ny)
          phi(0,ny+1)=phi(nx,1)
          phi(nx+1,ny+1)=phi(1,1)
          ux(0,0)=ux(nx,ny)
          ux(nx+1,0)=ux(1,ny)
          ux(0,ny+1)=ux(nx,1)
          ux(nx+1,ny+1)=ux(1,1)
          uy(0,0)=uy(nx,ny)
          uy(nx+1,0)=uy(1,ny)
          uy(0,ny+1)=uy(nx,1)
          uy(nx+1,ny+1)=uy(1,1)  
      
          do j=1,ny
             do i=1,nx
                nablaphix(i,j)=((phi(i+1,j+1)-phi(i-1,j+1))/6.0+2.0*(phi(i+1,j)-phi(i-1,j))/3.0+(phi(i+1,j-1)-phi(i-1,j-1))/6.0)/2.0/dh
                nablaphiy(i,j)=((phi(i+1,j+1)-phi(i+1,j-1))/6.0+2.0*(phi(i,j+1)-phi(i,j-1))/3.0+(phi(i-1,j+1)-phi(i-1,j-1))/6.0)/2.0/dh              
                chpoten(i,j) = 4.d0*beta*(phi(i,j)-phi_low)*(phi(i,j)-phi_high)*(phi(i,j)-0.5*(phi_low+phi_high))-kappa*(4.0*(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1))+(phi(i+1,j+1)+phi(i+1,j-1)+phi(i-1,j+1)+phi(i-1,j-1))-20.0*phi(i,j))/6.0/dh/dh
             end do
          end do  

          do i=1,nx
              chpoten(i,0)=chpoten(i,ny)
              chpoten(i,ny+1)=chpoten(i,1)
              nablaphix(i,0)=nablaphix(i,ny)
              nablaphix(i,ny+1)=nablaphix(i,1)
              nablaphiy(i,0)=nablaphiy(i,ny)
              nablaphiy(i,ny+1)=nablaphiy(i,1)          
          end do      
      
          do j=1,ny
              chpoten(0,j)=chpoten(nx,j)
              chpoten(nx+1,j)=chpoten(1,j)
              nablaphix(0,j)=nablaphix(nx,j)
              nablaphix(nx+1,j)=nablaphix(1,j)
              nablaphiy(0,j)=nablaphiy(nx,j)
              nablaphiy(nx+1,j)=nablaphiy(1,j)            
          end do
      
          chpoten(0,0)=chpoten(nx,ny)
          chpoten(nx+1,0)=chpoten(1,ny)
          chpoten(0,ny+1)=chpoten(nx,1)
          chpoten(nx+1,ny+1)=chpoten(1,1)
          nablaphix(0,0)=nablaphix(nx,ny)
          nablaphix(nx+1,0)=nablaphix(1,ny)
          nablaphix(0,ny+1)=nablaphix(nx,1)
          nablaphix(nx+1,ny+1)=nablaphix(1,1) 
          nablaphiy(0,0)=nablaphiy(nx,ny)
          nablaphiy(nx+1,0)=nablaphiy(1,ny)
          nablaphiy(0,ny+1)=nablaphiy(nx,1)
          nablaphiy(nx+1,ny+1)=nablaphiy(1,1)       
      
      do i=0,nx+1
          do j=0,ny+1
			  concentration(i,j)=4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
			  psi(i,j)=0.d0!concentration(i,j)
              concentration(i,j)=psi(i,j)*concentration(i,j)/(concentration(i,j)**2+1e-4)
              bodyforcex(i,j)=0.d0
              bodyforcey(i,j)=0.d0
              surtenforcex(i,j)=0.d0
              surtenforcey(i,j)=0.d0
              interforcex(i,j)=0.d0
              interforcey(i,j)=0.d0
              preforcex(i,j)=0.d0 
              preforcey(i,j)=0.d0
              pressure(i,j)=0.d0
              nablapsix(i,j)=0.d0
              nablapsiy(i,j)=0.d0
              nablapsix(i,j)=0.d0
              nablapsiy(i,j)=0.d0
              mchpoten(i,j)=chpoten(i,j)
          end do          
      end do
      
      do i=0,nx+1
          do j=0,ny+1
              rho(i,j)=rho_low*(1.d0-phi(i,j))+rho_high*phi(i,j)
              visfluid(i,j)=vis_low*(1.d0-phi(i,j))+vis_high*phi(i,j)
          end do          
      end do     
      
      do i=0,nx+1
          do j=0,ny+1
              tau_ddf(i,j)=3.d0*visfluid(i,j)/rho(i,j)/dh+0.5d0
          end do          
      end do      
      
      end subroutine      
!======================================================================c
!-------compute order parameter equilibrium
!----------------------------------------------------------------------c
     subroutine computeOrderEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,phi,chpoten,mchpoten,nablaphix,nablaphiy,mobility,intthick,eddg) 
!-------input parms
      integer nq, nx, ny
      integer ex(nq), ey(nq)
      double precision wt(nq),mobility,intthick,intnx,intny,ratio
      double precision, dimension(0:nx+1,0:ny+1)::ux,uy
      double precision, dimension(0:nx+1,0:ny+1)::phi,mchpoten,nablaphix,nablaphiy,chpoten
!-------output parms
      double precision,dimension(nq,0:nx+1,0:ny+1)::eddg
!----------------------------------------------------------------------c
      double precision eu,uv2
      integer i,j,k,iq
      
          do j=1,ny
             do i=0,nx+1
                uv2=ux(i,j)**2+uy(i,j)**2
                do iq=1,nq
                   eu=real(ex(iq))*ux(i,j)+real(ey(iq))*uy(i,j)         
                   eddg(iq,i,j)=phi(i,j)*wt(iq)+wt(iq)*3.d0*phi(i,j)*eu
                end do
             end do
          end do
     
    end subroutine   
!======================================================================c
!-------compute psi parameter equilibrium
!----------------------------------------------------------------------c
     subroutine computePsiEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,psi,eddh) 
!-------input parms
      integer nq, nx, ny
      integer ex(nq), ey(nq)
      double precision wt(nq)
      double precision, dimension(0:nx+1,0:ny+1)::ux,uy
      double precision, dimension(0:nx+1,0:ny+1)::psi
!-------output parms
      double precision,dimension(nq,0:nx+1,0:ny+1)::eddh
!----------------------------------------------------------------------c
      double precision eu,uv2
      integer i,j,k,iq
      
           do j=1,ny
             do i=0,nx+1
                uv2=ux(i,j)**2+uy(i,j)**2
                do iq=1,nq
                   eu=real(ex(iq))*ux(i,j)+real(ey(iq))*uy(i,j)         
                   eddh(iq,i,j)=psi(i,j)*wt(iq)+wt(iq)*3.d0*psi(i,j)*eu
                end do
             end do
          end do
     
    end subroutine  	
!======================================================================c
!-------compute density equilibrium
!----------------------------------------------------------------------c
      subroutine computeDenEDF2D(nq,ex,ey,wt,nx,ny,ux,uy,rho,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,eddf,dt) 
!-------input parms
      integer nq, nx, ny
      integer ex(nq), ey(nq)
      double precision wt(nq),dt
      double precision,dimension(0:nx+1,0:ny+1):: ux,uy,rho,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey
!-------output parms
      double precision eddf(nq,0:nx+1,0:ny+1)
!----------------------------------------------------------------------c
      double precision eu, uv2,gamma0,gamma1
      integer i,j,iq
      
      do j=1,ny
         do i=0,nx+1
            uv2=ux(i,j)**2+uy(i,j)**2
            do iq=1,nq
               eu=real(ex(iq))*ux(i,j)+real(ey(iq))*uy(i,j)
               gamma0=wt(iq)*(3.0d0*eu+4.50d0*eu**2-1.50d0*uv2)
               gamma1=wt(iq)*(1.d0+3.0d0*eu+4.50d0*eu**2-1.50d0*uv2)
               eddf(iq,i,j)=wt(iq)*pressure(i,j)+gamma0
               !+rho(i,j)/3.d0*wt(iq)*(3.0d0*eu+4.50d0*eu**2-1.50d0*uv2)&
               !-0.5d0*dt*(real(ex(iq))-ux(i,j))*(gamma0*preforcex(i,j)+gamma1*(surtenforcex(i,j)+bodyforcex(i,j)))&
               !-0.5d0*dt*(real(ey(iq))-uy(i,j))*(gamma0*preforcey(i,j)+gamma1*(surtenforcey(i,j)+bodyforcey(i,j)))
            end do
         end do
      end do
      end subroutine
!======================================================================c
!-------collision step for order parameter
!----------------------------------------------------------------------c
      subroutine collisionOrderDF2D(nq,ex,ey,wt,nx,ny,ux,ux_pre,uy,uy_pre,phi,phi_pre,chpoten,tau_ddg,dt,eddg,ddg,nablaphix,nablaphiy,mobility,intthick)
!-------input parms
      integer nq, nx, ny
      integer ex(nq), ey(nq)
      double precision wt(nq)
      double precision,dimension(0:nx+1,0:ny+1)::ux,uy,ux_pre,uy_pre
      double precision,dimension(0:nx+1,0:ny+1)::phi,phi_pre,chpoten,nablaphix,nablaphiy
      double precision tau_ddg,dt,mobility,intthick,intnx,intny,sumphi,sumnablaphix,sumnablaphiy,sumintnx,sumintny,lagx,lagy,sumphi2
      
      double precision eddg(nq,0:nx+1,0:ny+1),ddg(nq,0:nx+1,0:ny+1)
      double precision M(9,9),IM(9,9),MSIM(9,9),MSIM2(9,9),Sf(9),Sf2(9),DF(9)
      double precision intM(9,9),intM2(9,9),Fneq(9),MDF(9),Fneq2(9),MDF2(9)
!----------------------------------------------------------------------c
      integer i, j, k, iq, i1, j1, k1      
      
      M(1,:)=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
      M(2,:)=(/-1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0, -4.0/)
      M(3,:)=(/-2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, 4.0/)
      M(4,:)=(/1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0/)
      M(5,:)=(/-2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0/)
      M(6,:)=(/0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0, 0.0/)
      M(7,:)=(/0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0, 0.0/)
      M(8,:)=(/1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
      M(9,:)=(/0.0, 0.0, 0.0, 0, 1.0, -1.0, 1.0, -1.0, 0.0/)
      
      do iq=1,nq
           IM(iq,1)=M(1,iq)/9.d0
           IM(iq,2)=M(2,iq)/36.d0
           IM(iq,3)=M(3,iq)/36.d0
           IM(iq,4)=M(4,iq)/6.d0
           IM(iq,5)=M(5,iq)/12.d0
           IM(iq,6)=M(6,iq)/6.d0
           IM(iq,7)=M(7,iq)/12.d0
           IM(iq,8)=M(8,iq)/4.d0
           IM(iq,9)=M(9,iq)/4.d0
      end do
      
          do j=1,ny
             do i=0,nx+1
                tau_ddg=3.0*0.01+0.5 
                Sf(:)=(/1.0d0, 1.1d0, 1.1d0, 1.d0/tau_ddg, 1.d0/tau_ddg, 1.d0/tau_ddg, 1.0d0/tau_ddg, 1.2d0, 1.2d0/)
                do i1=1,9
                    do j1=1,9
                        intM(i1,j1)=Sf(i1)*M(i1,j1)
                    end do                
                end do
            
                do i1=1,9
                    do j1=1,9
                        MSIM(i1,j1)=0.d0
                        do k1=1,9
                            MSIM(i1,j1)=MSIM(i1,j1)+IM(i1,k1)*intM(k1,j1)
                        end do                    
                    end do                
                end do
            
                intnx=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                intny=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                                
                do i1=1,9
                    DF(i1)=wt(i1)*4.0/intthick*phi(i,j)*(1.d0-1.0*phi(i,j))*(real(ex(i1))*intnx+real(ey(i1))*intny)
                end do                  
            
                do i1=1,9
                    Fneq(i1)=0.d0
                    MDF(i1)=0.d0               
                    do j1=1,9
                        Fneq(i1)=Fneq(i1)+MSIM(i1,j1)*(ddg(j1,i,j)-eddg(j1,i,j))
                        MDF(i1)=MDF(i1)+MSIM(i1,j1)*DF(j1)*dt                  
                    end do                
                end do
            
                do iq=1,nq
                    ddg(iq,i,j)=ddg(iq,i,j)-Fneq(iq)+DF(iq)-0.5*MDF(iq)
                end do            
            
             end do
          end do
          !
          do iq=1,nq
              do j=1,ny
                  ddg(iq,0,j)=ddg(iq,nx,j)
                  ddg(iq,nx+1,j)=ddg(iq,1,j)
              end do
          end do
          
          do i=1,nx
                ddg(2,i,0)=ddg(4,i,1)
                ddg(5,i-1,0)=ddg(7,i,1)
                ddg(6,i+1,0)=ddg(8,i,1)
                
                ddg(4,i,ny+1)=ddg(2,i,ny)
                ddg(7,i+1,ny+1)=ddg(5,i,ny)
                ddg(8,i-1,ny+1)=ddg(6,i,ny)                
          end do          
      
      end subroutine
!======================================================================c
!-------collision step for psi parameter
!----------------------------------------------------------------------c
      subroutine collisionPsiDF2D(nq,ex,ey,wt,nx,ny,phi,psi,dt,eddh,ddh,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,intthick,tau_ddh)
!-------input parms
      integer nq, nx, ny
      integer ex(nq), ey(nq)
      double precision wt(nq)
      double precision,dimension(0:nx+1,0:ny+1)::phi,nablaphix,nablaphiy
      double precision,dimension(0:nx+1,0:ny+1)::psi,nablapsix,nablapsiy,nablaconx,nablacony
      double precision tau_ddh,dt,intthick,intnx,intny,lagx,lagy,nlagx,nlagy,ratio,source,gradsource
      
      double precision eddh(nq,0:nx+1,0:ny+1),ddh(nq,0:nx+1,0:ny+1)
      double precision M(9,9),IM(9,9),MSIM(9,9),MSIM2(9,9),Sf(9),Sf2(9),DF(9)
      double precision intM(9,9),intM2(9,9),Fneq(9),MDF(9),Fneq2(9),MDF2(9)
!----------------------------------------------------------------------c
      integer i, j, k, iq, i1, j1, k1      
      
      M(1,:)=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
      M(2,:)=(/-1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0, -4.0/)
      M(3,:)=(/-2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, 4.0/)
      M(4,:)=(/1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0/)
      M(5,:)=(/-2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0/)
      M(6,:)=(/0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0, 0.0/)
      M(7,:)=(/0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0, 0.0/)
      M(8,:)=(/1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
      M(9,:)=(/0.0, 0.0, 0.0, 0, 1.0, -1.0, 1.0, -1.0, 0.0/)
      
      do iq=1,nq
           IM(iq,1)=M(1,iq)/9.d0
           IM(iq,2)=M(2,iq)/36.d0
           IM(iq,3)=M(3,iq)/36.d0
           IM(iq,4)=M(4,iq)/6.d0
           IM(iq,5)=M(5,iq)/12.d0
           IM(iq,6)=M(6,iq)/6.d0
           IM(iq,7)=M(7,iq)/12.d0
           IM(iq,8)=M(8,iq)/4.d0
           IM(iq,9)=M(9,iq)/4.d0
      end do
      
          do j=1,ny
             do i=0,nx+1
                Sf(:)=(/1.0d0, 1.1d0, 1.1d0, 1.d0/tau_ddh, 1.d0/tau_ddh, 1.d0/tau_ddh, 1.0d0/tau_ddh, 1.2d0, 1.2d0/)
                do i1=1,9
                    do j1=1,9
                        intM(i1,j1)=Sf(i1)*M(i1,j1)
                    end do                
                end do
            
                do i1=1,9
                    do j1=1,9
                        MSIM(i1,j1)=0.d0
                        do k1=1,9
                            MSIM(i1,j1)=MSIM(i1,j1)+IM(i1,k1)*intM(k1,j1)
                        end do                    
                    end do                
                end do
            
                intnx=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                intny=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                
                source=(nablapsix(i,j)*intnx+nablapsiy(i,j)*intny)
                
                ratio=0.5
                gradsource=(1.0+ratio)*4.0*(1.0-2.0*phi(i,j))/intthick*psi(i,j)-ratio*source
                                    
                lagx=gradsource*intnx
                lagy=gradsource*intny
        
                               
                do i1=1,9
                    DF(i1)=wt(i1)*(real(ex(i1))*lagx+real(ey(i1))*lagy)
                end do                 
            
                do i1=1,9
                    Fneq(i1)=0.d0
                    MDF(i1)=0.d0               
                    do j1=1,9
                        Fneq(i1)=Fneq(i1)+MSIM(i1,j1)*(ddh(j1,i,j)-eddh(j1,i,j))
                        MDF(i1)=MDF(i1)+MSIM(i1,j1)*DF(j1)*dt                  
                    end do                
                end do
            
                do iq=1,nq
                    ddh(iq,i,j)=ddh(iq,i,j)-Fneq(iq)+DF(iq)-0.5*MDF(iq)
                end do            
            
             end do
          end do
          !
          do i=1,nx
                ddh(2,i,0)=ddh(4,i,1)
                ddh(5,i-1,0)=ddh(7,i,1)
                ddh(6,i+1,0)=ddh(8,i,1)
                
                ddh(4,i,ny+1)=ddh(2,i,ny)
                ddh(7,i+1,ny+1)=ddh(5,i,ny)
                ddh(8,i-1,ny+1)=ddh(6,i,ny)                
          end do
          
          do iq=1,nq
              !do i=1,nx
              !    ddh(iq,i,0)=ddh(iq,i,ny)
              !    ddh(iq,i,ny+1)=ddh(iq,i,1)
              !end do      
       
              do j=1,ny
                  ddh(iq,0,j)=ddh(iq,nx,j)
                  ddh(iq,nx+1,j)=ddh(iq,1,j)
              end do
      
              !ddh(iq,0,0)=ddh(iq,nx,ny)
              !ddh(iq,nx+1,0)=ddh(iq,1,ny)
              !ddh(iq,0,ny+1)=ddh(iq,nx,1)
              !ddh(iq,nx+1,ny+1)=ddh(iq,1,1)  
          end do
      
      end subroutine  	  
!======================================================================c
!-------collision step
!----------------------------------------------------------------------c
      subroutine collisionDenDF2D(nq,nx,ny,ex,ey,wt,ux,uy,rho,visfluid,pressure,tau_ddf,dt,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,visforcex,visforcey,nablaphix,nablaphiy,eddf,ddf,rho_low,rho_high)
!-------input parms
      integer nq, nx, ny
      integer ex(nq), ey(nq)
      double precision wt(nq),intviscoeffx,intviscoeffy,rho_low,rho_high
      double precision, dimension(0:nx+1,0:ny+1)::ux,uy,rho,visfluid,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,visforcex,visforcey,totalforcex,totalforcey
      double precision tau_ddf(0:nx+1,0:ny+1),dt,ref_rho
      double precision eddf(nq,0:nx+1,0:ny+1),ddf(nq,0:nx+1,0:ny+1)
      double precision,dimension (0:nx+1,0:ny+1)::nablaphix,nablaphiy      
      
      double precision M(9,9),IM(9,9),MSIM(9,9),Sf(9),DF(9),mgneq(9)
      double precision intM(9,9),Fneq(9),MDF(9)     
!----------------------------------------------------------------------c
	  double precision eu,uv2,gamma0,gamma1,force1,force2
      integer i, j, iq, jq, i1, j1, k

      M(1,:)=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
      M(2,:)=(/-1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0, -4.0/)
      M(3,:)=(/-2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, 4.0/)
      M(4,:)=(/1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0/)
      M(5,:)=(/-2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0/)
      M(6,:)=(/0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0, 0.0/)
      M(7,:)=(/0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0, 0.0/)
      M(8,:)=(/1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
      M(9,:)=(/0.0, 0.0, 0.0, 0, 1.0, -1.0, 1.0, -1.0, 0.0/)
      
      do iq=1,nq
           IM(iq,1)=M(1,iq)/9.d0
           IM(iq,2)=M(2,iq)/36.d0
           IM(iq,3)=M(3,iq)/36.d0
           IM(iq,4)=M(4,iq)/6.d0
           IM(iq,5)=M(5,iq)/12.d0
           IM(iq,6)=M(6,iq)/6.d0
           IM(iq,7)=M(7,iq)/12.d0
           IM(iq,8)=M(8,iq)/4.d0
           IM(iq,9)=M(9,iq)/4.d0
      end do
      
      do i=0,nx+1
          do j=1,ny
            Sf(:)=(/0.d0, 1.d0/tau_ddf(i,j), 1.d0/tau_ddf(i,j), 0.d0, 8.0*(2.d0-1.d0/tau_ddf(i,j))/(8.d0-1.d0/tau_ddf(i,j)), 0.d0, 8.0*(2.d0-1.d0/tau_ddf(i,j))/(8.d0-1.d0/tau_ddf(i,j)), 1.d0/tau_ddf(i,j), 1.d0/tau_ddf(i,j)/)
            do i1=1,9
                do j1=1,9
                    intM(i1,j1)=Sf(i1)*M(i1,j1)
                end do                
            end do
            
            do i1=1,9
                do j1=1,9
                    MSIM(i1,j1)=0.d0
                    do k=1,9
                        MSIM(i1,j1)=MSIM(i1,j1)+IM(i1,k)*intM(k,j1)
                    end do                    
                end do                
            end do
            
              do iq=1,nq
                  mgneq(iq)=0.d0
                  do jq=1,nq
                        mgneq(iq)=mgneq(iq)+MSIM(iq,jq)*(ddf(jq,i,j)-eddf(jq,i,j))    
                  end do                  
              end do
              intviscoeffx=0.d0
              intviscoeffy=0.d0              
              do iq=1,nq
                  intviscoeffx=intviscoeffx+real(ex(iq))*real(ex(iq))*mgneq(iq)*(rho_high*nablaphix(i,j)-rho_low*nablaphix(i,j))
                  intviscoeffx=intviscoeffx+real(ex(iq))*real(ey(iq))*mgneq(iq)*(rho_high*nablaphiy(i,j)-rho_low*nablaphiy(i,j))
                  intviscoeffy=intviscoeffy+real(ey(iq))*real(ex(iq))*mgneq(iq)*(rho_high*nablaphix(i,j)-rho_low*nablaphix(i,j))
                  intviscoeffy=intviscoeffy+real(ey(iq))*real(ey(iq))*mgneq(iq)*(rho_high*nablaphiy(i,j)-rho_low*nablaphiy(i,j))               
              end do
              visforcex(i,j)=-intviscoeffx*3.0/dt*visfluid(i,j)/rho(i,j)
              visforcey(i,j)=-intviscoeffy*3.0/dt*visfluid(i,j)/rho(i,j)
          end do          
      end do
      
      do j=1,ny
          visforcex(0,j)=visforcex(nx,j)
          visforcex(nx+1,j)=visforcex(1,j)
          visforcey(0,j)=visforcey(nx,j)
          visforcey(nx+1,j)=visforcey(1,j)          
      end do            
      
      do j=1,ny
         do i=1,nx
            Sf(:)=(/0.d0, 1.d0/tau_ddf(i,j), 1.d0/tau_ddf(i,j), 0.d0, 8.0*(2.d0-1.d0/tau_ddf(i,j))/(8.d0-1.d0/tau_ddf(i,j)), 0.d0, 8.0*(2.d0-1.d0/tau_ddf(i,j))/(8.d0-1.d0/tau_ddf(i,j)), 1.d0/tau_ddf(i,j), 1.d0/tau_ddf(i,j)/)

            totalforcex(i,j)=preforcex(i,j)+surtenforcex(i,j)+bodyforcex(i,j)+visforcex(i,j)
            totalforcey(i,j)=preforcey(i,j)+surtenforcey(i,j)+bodyforcey(i,j)+visforcey(i,j)            
            do iq=1,nq
               eu=real(ex(iq))*ux(i,j)+real(ey(iq))*uy(i,j)
               eforce=(3.d0*(real(ex(iq))-ux(i,j))+9.d0*eu*real(ex(iq)))*totalforcex(i,j)+(3.d0*(real(ey(iq))-uy(i,j))+9.d0*eu*real(ey(iq)))*totalforcey(i,j)       
               DF(iq)=wt(iq)*eforce*dt/rho(i,j)               
            end do
            
            do i1=1,9
                do j1=1,9
                    intM(i1,j1)=Sf(i1)*M(i1,j1)
                end do                
            end do
            
            do i1=1,9
                do j1=1,9
                    MSIM(i1,j1)=0.d0
                    do k=1,9
                        MSIM(i1,j1)=MSIM(i1,j1)+IM(i1,k)*intM(k,j1)
                    end do                    
                end do                
            end do
            
            do i1=1,9
                Fneq(i1)=0.d0
                MDF(i1)=0.d0
                do j1=1,9
                    Fneq(i1)=Fneq(i1)+MSIM(i1,j1)*(ddf(j1,i,j)-eddf(j1,i,j))
                    MDF(i1)=MDF(i1)+MSIM(i1,j1)*DF(j1)
                end do                
            end do
            
            do iq=1,nq
                ddf(iq,i,j)=ddf(iq,i,j)-Fneq(iq)+DF(iq)-0.5*MDF(iq)
            end do            
            
         end do
      end do        

      do iq=1,nq
          do j=1,ny
              ddf(iq,0,j)=ddf(iq,nx,j)
              ddf(iq,nx+1,j)=ddf(iq,1,j)
          end do          
      end do
      
      do i=1,nx
          ddf(2,i,0)=ddf(4,i,1)
          ddf(5,i-1,0)=ddf(7,i,1)-6.0*wt(7)*real(ex(7))*(-0.0)
          ddf(6,i+1,0)=ddf(8,i,1)-6.0*wt(8)*real(ex(8))*(-0.0)
          
          ddf(4,i,ny+1)=ddf(2,i,ny)
          ddf(7,i+1,ny+1)=ddf(5,i,ny)-6.0*wt(5)*real(ex(5))*0.00
          ddf(8,i-1,ny+1)=ddf(6,i,ny)-6.0*wt(6)*real(ex(6))*0.00
      end do        
                    
      end subroutine
!======================================================================c
!-------stream step
!----------------------------------------------------------------------c
      subroutine streamOrderDF2D(nq,ex,ey,wt,nx,ny,ddg,work3D)
!-------input parms
      integer nq,nx,ny
      integer ex(nq),ey(nq)
      double precision wt(nq)
      double precision ddg(nq,0:nx+1,0:ny+1),work3D(nq,0:nx+1,0:ny+1)
!-------output parms
!----------------------------------------------------------------------c
      integer i, j, i1, j1, iq, k
          do j=1,ny
             do i=1,nx
                do iq=1,nq
                   i1=i-ex(iq)
                   j1=j-ey(iq)
               
                   work3D(iq,i,j)=ddg(iq,i1,j1) 
                end do
             end do
          end do
      
          do j=1,ny
             do i=1,nx
                do iq=1,nq
                   ddg(iq,i,j)=work3D(iq,i,j)
                end do
             end do
          end do
      !
      end subroutine
!======================================================================c
!-------psi stream step
!----------------------------------------------------------------------c
      subroutine streamPsiDF2D(nq,ex,ey,wt,nx,ny,ddh,work3D)
!-------input parms
      integer nq,nx,ny
      integer ex(nq),ey(nq)
      double precision wt(nq)
      double precision ddh(nq,0:nx+1,0:ny+1),work3D(nq,0:nx+1,0:ny+1)
!-------output parms
!----------------------------------------------------------------------c
      integer i, j, i1, j1, iq, k
	  
          do j=1,ny
             do i=1,nx
                do iq=1,nq
                   i1=i-ex(iq)
                   j1=j-ey(iq)

               
                   work3D(iq,i,j)=ddh(iq,i1,j1) 
                end do
             end do
          end do
      
          do j=1,ny
             do i=1,nx
                do iq=1,nq
                   ddh(iq,i,j)=work3D(iq,i,j)
                end do
             end do
          end do

      end subroutine	  
!======================================================================c
!-------stream step
!----------------------------------------------------------------------c
      subroutine streamDenDF2D(nq,ex,ey,wt,nx,ny,ddf,work3D)
!-------input parms
      integer nq,nx,ny
      integer ex(nq),ey(nq)
      double precision wt(nq)
      double precision ddf(nq,0:nx+1,0:ny+1),work3D(nq,0:nx+1,0:ny+1)
!-------output parms
!----------------------------------------------------------------------c
      integer i, j, i1, j1, iq
      
      do j=1,ny
         do i=1,nx
            do iq=1,nq
               i1=i-ex(iq)
               j1=j-ey(iq)
               
               work3D(iq,i,j)=ddf(iq,i1,j1) 
            end do
         end do
      end do
      do j=1,ny
         do i=1,nx
            do iq=1,nq
               ddf(iq,i,j)=work3D(iq,i,j)
            end do
         end do
      end do    
      
    end subroutine
!======================================================================c
!-------set Order distribution function boundary condition
!----------------------------------------------------------------------c
      subroutine setOrderfBC2D(nq,nx,ny,ddg,eddg) 
!-------input parms
      integer nq, nx, ny
      double precision ddg(nq,0:nx+1,0:ny+1),eddg(nq,0:nx+1,0:ny+1)
      integer i,j,iq
                 
        do i=0,nx+1
            do iq=1,nq
                ddg(iq,i,1)=eddg(iq,i,1)+ddg(iq,i,2)-eddg(iq,i,2)
                ddg(iq,i,ny)=eddg(iq,i,ny)+ddg(iq,i,ny-1)-eddg(iq,i,ny-1)
            end do            
        end do           

      end subroutine   
!======================================================================c
!-------set Psi distribution function boundary condition
!----------------------------------------------------------------------c
      subroutine setPsifBC2D(nq,nx,ny,ddh,eddh) 
!-------input parms
      integer nq, nx, ny
      double precision ddh(nq,0:nx+1,0:ny+1),eddh(nq,0:nx+1,0:ny+1)
      integer i,j,iq
                 
        do i=0,nx+1
            do iq=1,nq
                ddh(iq,i,1)=eddh(iq,i,1)!+ddh(iq,i,2)-eddh(iq,i,2)
                ddh(iq,i,ny)=eddh(iq,i,ny)!+ddh(iq,i,ny-1)-eddh(iq,i,ny-1)
            end do            
        end do           

      end subroutine   	  
!======================================================================c
!-------set density distribution function boundary condition
!----------------------------------------------------------------------c
      subroutine setDenfBC2D(nq,nx,ny,ddf,eddf) 
!-------input parms
      integer nq, nx, ny
      double precision ddf(nq,0:nx+1,0:ny+1),eddf(nq,0:nx+1,0:ny+1)
      integer i,j,iq
         
        do j=1,ny
            do iq=1,nq
                ddf(iq,1,j)=eddf(iq,1,j)+ddf(iq,2,j)-eddf(iq,2,j)
                ddf(iq,nx,j)=eddf(iq,nx,j)+ddf(iq,nx-1,j)-eddf(iq,nx-1,j)
            end do            
        end do
        
        do i=1,nx
            do iq=1,nq
                ddf(iq,i,1)=eddf(iq,i,1)+ddf(iq,i,2)-eddf(iq,i,2)
                ddf(iq,i,ny)=eddf(iq,i,ny)+ddf(iq,i,ny-1)-eddf(iq,i,ny-1)
            end do            
        end do           
      end subroutine        
!======================================================================c
!-------compute macro vars for Oreder parameter
!----------------------------------------------------------------------c
      subroutine computeOrderMacro2D(nq,nx,ny,ex,ey,wt,ddg,eddg,ddh,phi,phi_pre,psi,concentration,chpoten,mchpoten,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,beta,kappa,dh,&
      phi_low,phi_high,rho,ux,ux_pre,uy,uy_pre,istep,dt,x,y,tau_ddf,tau_ddf_low,tau_ddf_high,mobility,intthick,surtencoefficent,rho_low,rho_high,vis_low,vis_high,visfluid) 
!-------input parms
      integer nq, nx, ny, istep
      integer ex(nq),ey(nq)
      double precision dh,beta,kappa,phi_low,phi_high, dt, pi,tau_ddf_low,tau_ddf_high,wt(nq),intthick,mobility,temp1,temp2,surtencoefficent,intnx,intny,lagx,lagy,aa,bb,rho_low,rho_high,vis_low,vis_high
      double precision ddg(nq,0:nx+1,0:ny+1),eddg(nq,0:nx+1,0:ny+1),ddh(nq,0:nx+1,0:ny+1),x(0:nx+1),y(0:ny+1)
!-------output parms
      double precision, dimension(0:nx+1,0:ny+1)::rho,ux,ux_pre,uy,uy_pre,tau_ddf,visfluid,psi,concentration,nablapsix,nablapsiy,nablaconx,nablacony
      double precision, dimension(0:nx+1,0:ny+1)::phi,phi_pre,chpoten,nablaphix,nablaphiy,mchpoten,nablaphi_x_y,nablaphi_x_x,nablaphi_y_y
!----------------------------------------------------------------------c
      integer i, j, k, iq
      
      pi=4.d0*atan(1.d0)
    
        do j=1,ny
            do i=1,nx           
                phi(i,j) = 0.d0           
                do iq=1,nq
                    phi(i,j) = phi(i,j) +ddg(iq,i,j)     
                end do
            end do
        end do       
      
          do i=1,nx
              phi(i,0)=phi(i,ny)
              phi(i,ny+1)=phi(i,1)
          end do      
      
          do j=1,ny
              phi(0,j)=phi(nx,j)
              phi(nx+1,j)=phi(1,j)
          end do
      
          phi(0,0)=phi(nx,ny)
          phi(nx+1,0)=phi(1,ny)
          phi(0,ny+1)=phi(nx,1)
          phi(nx+1,ny+1)=phi(1,1)
      

          do j=2,ny-1
             do i=1,nx
                !nablaphix(i,j)=((phi(i+1,j+1)-phi(i-1,j+1))/6.0+2.0*(phi(i+1,j)-phi(i-1,j))/3.0+(phi(i+1,j-1)-phi(i-1,j-1))/6.0)/2.0/dh
                !nablaphiy(i,j)=((phi(i+1,j+1)-phi(i+1,j-1))/6.0+2.0*(phi(i,j+1)-phi(i,j-1))/3.0+(phi(i-1,j+1)-phi(i-1,j-1))/6.0)/2.0/dh 
                nablaphi_x_y(i,j)=(phi(i+1,j+1)-phi(i+1,j-1)-phi(i-1,j+1)+phi(i-1,j-1))/4.0/dh/dh
                nablaphi_x_x(i,j)=((phi(i+1,j+1)-2.0*phi(i,j+1)+phi(i-1,j+1))+10.0*(phi(i+1,j)-2.0*phi(i,j)+phi(i-1,j))+(phi(i+1,j-1)-2.0*phi(i,j-1)+phi(i-1,j-1)))/12.0/dh/dh
                nablaphi_y_y(i,j)=((phi(i+1,j+1)-2.0*phi(i+1,j)+phi(i+1,j-1))+10.0*(phi(i,j+1)-2.0*phi(i,j)+phi(i,j-1))+(phi(i-1,j+1)-2.0*phi(i-1,j)+phi(i-1,j-1)))/12.0/dh/dh
                 
                nablaphix(i,j)=0.d0
                nablaphiy(i,j)=0.d0
                do iq=1,nq
                    nablaphix(i,j)=nablaphix(i,j)+wt(iq)*3.d0/dh*real(ex(iq))*phi(i+ex(iq),j+ey(iq))
                    nablaphiy(i,j)=nablaphiy(i,j)+wt(iq)*3.d0/dh*real(ey(iq))*phi(i+ex(iq),j+ey(iq))
                end do
				beta=12.d0*surtencoefficent*(1.d0-0.2d0*concentration(i,j)*0.6)/intthick/(phi_high-phi_low)**4
				kappa=3.d0*intthick*surtencoefficent*(1.d0-0.2d0*concentration(i,j)*0.6)/2.d0/(phi_high-phi_low)**2				
                chpoten(i,j) = 4.d0*beta*(phi(i,j)-phi_low)*(phi(i,j)-phi_high)*(phi(i,j)-0.5*(phi_low+phi_high))-kappa*(4.0*(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1))+(phi(i+1,j+1)+phi(i+1,j-1)+phi(i-1,j+1)+phi(i-1,j-1))-20.0*phi(i,j))/6.0/dh/dh
                aa=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                bb=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                
                mchpoten(i,j) = 4.d0*beta*(phi(i,j)-phi_low)*(phi(i,j)-phi_high)*(phi(i,j)-0.5*(phi_low+phi_high))-kappa*(aa**2*nablaphi_x_x(i,j)+2.0*aa*bb*nablaphi_x_y(i,j)+bb**2*nablaphi_y_y(i,j))+log((phi(i,j)+1e-8)/(1.0+1e-8-phi(i,j)))/(phi(i,j)+1e-8)/(1.0+1e-8-phi(i,j))*1e-3
             end do
          end do
          
          j=1
          do i=1,nx
                nablaphix(i,j)=((phi(i+1,j)-phi(i-1,j)))/2.0/dh
                nablaphiy(i,j)=nablaphiy(i,j+1)
                nablaphi_x_y(i,j)=nablaphi_x_y(i,j+1)
				beta=12.d0*surtencoefficent*(1.d0-0.2d0*concentration(i,j)*0.6)/intthick/(phi_high-phi_low)**4
				kappa=3.d0*intthick*surtencoefficent*(1.d0-0.2d0*concentration(i,j)*0.6)/2.d0/(phi_high-phi_low)**2					
                chpoten(i,j) = 4.d0*beta*(phi(i,j)-phi_low)*(phi(i,j)-phi_high)*(phi(i,j)-0.5*(phi_low+phi_high))-kappa*((phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/dh/dh+(phi(i,j+2)-2.d0*phi(i,j+1)+phi(i,j))/dh/dh)
                mchpoten(i,j) =mchpoten(i,j+1)
          end do
      
          j=ny
          do i=1,nx
                nablaphix(i,j)=((phi(i+1,j)-phi(i-1,j)))/2.0/dh
                nablaphiy(i,j)=nablaphiy(i,j-1)
                nablaphi_x_y(i,j)=nablaphi_x_y(i,j-1)
				beta=12.d0*surtencoefficent*(1.d0-0.2d0*concentration(i,j)*0.6)/intthick/(phi_high-phi_low)**4
				kappa=3.d0*intthick*surtencoefficent*(1.d0-0.2d0*concentration(i,j)*0.6)/2.d0/(phi_high-phi_low)**2					
                chpoten(i,j) = 4.d0*beta*(phi(i,j)-phi_low)*(phi(i,j)-phi_high)*(phi(i,j)-0.5*(phi_low+phi_high))-kappa*((phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/dh/dh+(phi(i,j)-2.d0*phi(i,j-1)+phi(i,j-2))/dh/dh)
                mchpoten(i,j) =mchpoten(i,j-1)
          end do          
      
      
          !do i=1,nx
          !    chpoten(i,0)=chpoten(i,ny)
          !    chpoten(i,ny+1)=chpoten(i,1)
          !    nablaphix(i,0)=nablaphix(i,ny)
          !    nablaphix(i,ny+1)=nablaphix(i,1)
          !    nablaphiy(i,0)=nablaphiy(i,ny)
          !    nablaphiy(i,ny+1)=nablaphiy(i,1)          
          !end do      
      
          do j=1,ny
              chpoten(0,j)=chpoten(nx,j)
              chpoten(nx+1,j)=chpoten(1,j)
              nablaphix(0,j)=nablaphix(nx,j)
              nablaphix(nx+1,j)=nablaphix(1,j)
              nablaphiy(0,j)=nablaphiy(nx,j)
              nablaphiy(nx+1,j)=nablaphiy(1,j)
              mchpoten(0,j)=mchpoten(nx,j)
              mchpoten(nx+1,j)=mchpoten(1,j)              
          end do
      
          !chpoten(0,0)=chpoten(nx,ny)
          !chpoten(nx+1,0)=chpoten(1,ny)
          !chpoten(0,ny+1)=chpoten(nx,1)
          !chpoten(nx+1,ny+1)=chpoten(1,1)
          !nablaphix(0,0)=nablaphix(nx,ny)
          !nablaphix(nx+1,0)=nablaphix(1,ny)
          !nablaphix(0,ny+1)=nablaphix(nx,1)
          !nablaphix(nx+1,ny+1)=nablaphix(1,1) 
          !nablaphiy(0,0)=nablaphiy(nx,ny)
          !nablaphiy(nx+1,0)=nablaphiy(1,ny)
          !nablaphiy(0,ny+1)=nablaphiy(nx,1)
          !nablaphiy(nx+1,ny+1)=nablaphiy(1,1)   

          do j=1,ny
             do i=1,nx           
                psi(i,j) = 0.d0           
                do iq=1,nq
                   psi(i,j) = psi(i,j) +ddh(iq,i,j)     
                end do
                concentration(i,j)=4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
                concentration(i,j)=psi(i,j)*concentration(i,j)/(concentration(i,j)**2+1e-4)
             end do
          end do       
          
          do j=2,ny-1
             do i=1,nx 
                nablapsix(i,j)=0.d0
                nablapsiy(i,j)=0.d0
                do iq=1,nq
                    nablapsix(i,j)=nablapsix(i,j)+wt(iq)*3.d0/dh*real(ex(iq))*psi(i+ex(iq),j+ey(iq))
                    nablapsiy(i,j)=nablapsiy(i,j)+wt(iq)*3.d0/dh*real(ey(iq))*psi(i+ex(iq),j+ey(iq))
                end do             
             end do
          end do
          
          j=1
          do i=1,nx
                nablapsix(i,j)=((psi(i+1,j)-psi(i-1,j)))/2.0/dh
                nablapsiy(i,j)=nablapsiy(i,j+1)
          end do
      
          j=ny
          do i=1,nx
                nablapsix(i,j)=((psi(i+1,j)-psi(i-1,j)))/2.0/dh
                nablapsiy(i,j)=nablapsiy(i,j-1)          
          end do            
          
          do j=1,ny
             do i=1,nx           
                temp1=4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
                
                intnx=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                intny=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                                    
                lagx=4.0/intthick*(1.d0-2.0*phi(i,j))*intnx*psi(i,j)
                lagy=4.0/intthick*(1.d0-2.0*phi(i,j))*intny*psi(i,j)
                
                nablaconx(i,j) = (nablapsix(i,j)-lagx)!*temp1/(temp1**2+1e-5)
                nablacony(i,j) = (nablapsiy(i,j)-lagy)!*temp1/(temp1**2+1e-5)
             end do
          end do
          
          !do j=2,ny-1
          !   do i=2,nx-1                           
          !      do iq=1,nq
          !          nablaconx(i,j)=nablaconx(i,j)+wt(iq)*3.d0/dh*real(ex(iq))*concentration(i+ex(iq),j+ey(iq))
          !          nablacony(i,j)=nablacony(i,j)+wt(iq)*3.d0/dh*real(ey(iq))*concentration(i+ex(iq),j+ey(iq))
          !      end do   
          !   end do
          !end do            
      
      do i=0,nx+1
          do j=0,ny+1
              rho(i,j)=rho_low*(1.d0-phi(i,j))+rho_high*phi(i,j)
              visfluid(i,j)=vis_low*(1.d0-phi(i,j))+vis_high*phi(i,j)
          end do          
      end do     
      
      do i=0,nx+1
          do j=0,ny+1
              tau_ddf(i,j)=3.d0*visfluid(i,j)/rho(i,j)/dh+0.5d0
          end do          
      end do        
      
      end subroutine
!====================================================================c
!-------compute force
!--------------------------------------------------------------------c  
     subroutine computeForce(nx,ny,nq,ex,ey,wt,chpoten,phi,rho,nablaphix,nablaphiy,nablaconx,nablacony,curv,concentration,bodyforcex,bodyforcey,surtenforcex,surtenforcey,interforcex,interforcey,preforcex,preforcey,dh,rho_high,rho_low,surtencoefficent,intthick)
     integer nx,ny,nq
     integer ex(nq),ey(nq)
     double precision wt(nq),dh,rho_high,rho_low,surtencoefficent,sigma_T,temp,intthick
     double precision,dimension (0:nx+1,0:ny+1)::rho,bodyforcex,bodyforcey,surtenforcex,surtenforcey,interforcex,interforcey,preforcex,preforcey
     double precision,dimension (0:nx+1,0:ny+1)::phi,chpoten,nablaphix,nablaphiy
     double precision,dimension (0:nx+1,0:ny+1)::sigma12,sigma13,sigma14,sigma23,sigma24,sigma34,nablaconx,nablacony,concentration,intnx,intny,curv
     
     
    do j=1,ny
        do i=1,nx
            intnx(i,j)=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
            intny(i,j)=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
            curv(i,j)=0.d0
        end do
    end do
    
    do j=2,ny-1
        do i=2,nx-1
            curv(i,j)=0.d0
            do iq=1,nq
                curv(i,j)=curv(i,j)+wt(iq)*3.d0/dh*real(ex(iq))*intnx(i+ex(iq),j+ey(iq))
                curv(i,j)=curv(i,j)+wt(iq)*3.d0/dh*real(ey(iq))*intny(i+ex(iq),j+ey(iq))
            end do
            if(abs(phi(i,j)).le.1e-5.or.abs(phi(i,j)-1.d0).le.1e-5) curv(i,j)=0.d0
        end do
    end do
     
     do i=1,nx
         do j=1,ny
             bodyforcex(i,j)=0.d0
             bodyforcey(i,j)=-(rho(i,j))*1e-4/60.0 

             
             surtenforcex(i,j)=chpoten(i,j)*nablaphix(i,j)
             surtenforcey(i,j)=chpoten(i,j)*nablaphiy(i,j)
			 !sigma_T=-0.2d0*surtencoefficent*0.6d0/(1.d0-0.6*concentration(i,j))
    !         temp=6.d0*phi(i,j)*(1.d0-phi(i,j))!**2/intthick
    !         surtenforcex(i,j)=temp*nablaconx(i,j)*sigma_T-temp*(intnx(i,j)*nablaconx(i,j)+intny(i,j)*nablacony(i,j))*intnx(i,j)*sigma_T-temp*curv(i,j)*intnx(i,j)*surtencoefficent*(1.d0+0.2d0*log(1.d0-concentration(i,j)*0.6))*4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
    !         surtenforcey(i,j)=temp*nablacony(i,j)*sigma_T-temp*(intnx(i,j)*nablaconx(i,j)+intny(i,j)*nablacony(i,j))*intny(i,j)*sigma_T-temp*curv(i,j)*intny(i,j)*surtencoefficent*(1.d0+0.2d0*log(1.d0-concentration(i,j)*0.6))*4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
			 
         end do         
     end do
     
      !do i=1,nx
      !    bodyforcex(i,0)=bodyforcex(i,ny)
      !    bodyforcex(i,ny+1)=bodyforcex(i,1)
      !    bodyforcey(i,0)=bodyforcey(i,ny)
      !    bodyforcey(i,ny+1)=bodyforcey(i,1)   
      !    surtenforcex(i,0)=surtenforcex(i,ny)
      !    surtenforcex(i,ny+1)=surtenforcex(i,1)  
      !    surtenforcey(i,0)=surtenforcey(i,ny)
      !    surtenforcey(i,ny+1)=surtenforcey(i,1)
      !    preforcex(i,0)=preforcex(i,ny)
      !    preforcex(i,ny+1)=preforcex(i,1)
      !    preforcey(i,0)=preforcey(i,ny)
      !    preforcey(i,ny+1)=preforcey(i,1)          
      !end do      
      
      do j=1,ny
          bodyforcex(0,j)=bodyforcex(nx,j)
          bodyforcex(nx+1,j)=bodyforcex(1,j)
          bodyforcey(0,j)=bodyforcey(nx,j)
          bodyforcey(nx+1,j)=bodyforcey(1,j)
          surtenforcex(0,j)=surtenforcex(nx,j)
          surtenforcex(nx+1,j)=surtenforcex(1,j)   
          surtenforcey(0,j)=surtenforcey(nx,j)
          surtenforcey(nx+1,j)=surtenforcey(1,j)
          preforcex(0,j)=preforcex(nx,j)
          preforcex(nx+1,j)=preforcex(1,j)
          preforcey(0,j)=preforcey(nx,j)
          preforcey(nx+1,j)=preforcey(1,j)          
      end do
      
      !bodyforcex(0,0)=bodyforcex(nx,ny)
      !bodyforcex(nx+1,0)=bodyforcex(1,ny)
      !bodyforcex(0,ny+1)=bodyforcex(nx,1)
      !bodyforcex(nx+1,ny+1)=bodyforcex(1,1)     
      !bodyforcey(0,0)=bodyforcey(nx,ny)
      !bodyforcey(nx+1,0)=bodyforcey(1,ny)
      !bodyforcey(0,ny+1)=bodyforcey(nx,1)
      !bodyforcey(nx+1,ny+1)=bodyforcey(1,1)
      !surtenforcex(0,0)=surtenforcex(nx,ny)
      !surtenforcex(nx+1,0)=surtenforcex(1,ny)
      !surtenforcex(0,ny+1)=surtenforcex(nx,1)
      !surtenforcex(nx+1,ny+1)=surtenforcex(1,1)
      !surtenforcey(0,0)=surtenforcey(nx,ny)
      !surtenforcey(nx+1,0)=surtenforcey(1,ny)
      !surtenforcey(0,ny+1)=surtenforcey(nx,1)
      !surtenforcey(nx+1,ny+1)=surtenforcey(1,1) 
      !preforcex(0,0)=preforcex(nx,ny)
      !preforcex(nx+1,0)=preforcex(1,ny)
      !preforcex(0,ny+1)=preforcex(nx,1)
      !preforcex(nx+1,ny+1)=preforcex(1,1)
      !preforcey(0,0)=preforcey(nx,ny)
      !preforcey(nx+1,0)=preforcey(1,ny)
      !preforcey(0,ny+1)=preforcey(nx,1)
      !preforcey(nx+1,ny+1)=preforcey(1,1)      
      
     end subroutine 
!====================================================================c
!-------compute macro vars
!--------------------------------------------------------------------c
      subroutine computeMacro2D(nq,ex,ey,wt,nx,ny,ddf,ux,uy,rho,pressure,nablaphix,nablaphiy,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,visforcex,visforcey,dt,rho_low,rho_high)
!-------input parms
      integer nq, nx, ny, istep
      integer ex(nq),ey(nq)
      double precision wt(nq),dt,rho_low,rho_high
      double precision x(0:nx+1),y(0:ny+1),ddf(nq,0:nx+1,0:ny+1)
      double precision,dimension (0:nx+1,0:ny+1)::nablaphix,nablaphiy
!-------output parms
      double precision, dimension(0:nx+1,0:ny+1)::ux,uy,rho,pressure,bodyforcex,bodyforcey,surtenforcex,surtenforcey,preforcex,preforcey,visforcex,visforcey
!----------------------------------------------------------------------c
      integer i, j
      pi=4.d0*atan(1.d0)
      do j=1,ny
         do i=1,nx
            pressure(i,j)=0.d0
            do iq=1,nq
                pressure(i,j)=pressure(i,j)+ddf(iq,i,j)
            end do              
            preforcex(i,j)=-1.d0/3.d0*pressure(i,j)*(rho_high*nablaphix(i,j)-rho_low*nablaphix(i,j))
            preforcey(i,j)=-1.d0/3.d0*pressure(i,j)*(rho_high*nablaphiy(i,j)-rho_low*nablaphiy(i,j))
            ux(i,j)  = 0.d0
            uy(i,j)  = 0.d0
            do iq=1,nq
               ux(i,j)  = ux(i,j) + real(ex(iq))*ddf(iq,i,j)
               uy(i,j)  = uy(i,j) + real(ey(iq))*ddf(iq,i,j)
            end do
            ux(i,j)  = ux(i,j)+0.5d0*(surtenforcex(i,j)+bodyforcex(i,j)+visforcex(i,j)+preforcex(i,j))*dt/rho(i,j)
            uy(i,j)  = uy(i,j)+0.5d0*(surtenforcey(i,j)+bodyforcey(i,j)+visforcey(i,j)+preforcey(i,j))*dt/rho(i,j)         
         end do
      end do
      !
      !do i=1,nx
      !    rho(i,0)=rho(i,ny)
      !    rho(i,ny+1)=rho(i,1)
      !    ux(i,0)=ux(i,ny)
      !    ux(i,ny+1)=ux(i,1)
      !    uy(i,0)=uy(i,ny)
      !    uy(i,ny+1)=uy(i,1)
      !    pressure(i,0)=pressure(i,ny)
      !    pressure(i,ny+1)=pressure(i,1)          
      !end do      
      
      do j=1,ny
          rho(0,j)=rho(nx,j)
          rho(nx+1,j)=rho(1,j)          
          ux(0,j)=ux(nx,j)
          ux(nx+1,j)=ux(1,j)
          uy(0,j)=uy(nx,j)
          uy(nx+1,j)=uy(1,j)
          pressure(0,j)=pressure(nx,j)
          pressure(nx+1,j)=pressure(1,j)          
      end do

      !rho(0,0)=rho(nx,ny)
      !rho(nx+1,0)=rho(1,ny)
      !rho(0,ny+1)=rho(nx,1)
      !rho(nx+1,ny+1)=rho(1,1)      
      !ux(0,0)=ux(nx,ny)
      !ux(nx+1,0)=ux(1,ny)
      !ux(0,ny+1)=ux(nx,1)
      !ux(nx+1,ny+1)=ux(1,1)   
      !uy(0,0)=uy(nx,ny)
      !uy(nx+1,0)=uy(1,ny)
      !uy(0,ny+1)=uy(nx,1)
      !uy(nx+1,ny+1)=uy(1,1)  
      !pressure(0,0)=pressure(nx,ny)
      !pressure(nx+1,0)=pressure(1,ny)
      !pressure(0,ny+1)=pressure(nx,1)
      !pressure(nx+1,ny+1)=pressure(1,1)        
      
    end subroutine
!======================================================================c
!-------set macro vars  boundary condition
!----------------------------------------------------------------------c
      subroutine setMacroOrderBC2D(nx,ny,phi,phi_pre,psi,concentration,chpoten,phi_high,phi_low,dh)
!-------input parms
      double precision phi_high,phi_low,dh,pi
      integer nx, ny
      double precision,dimension(0:nx+1,0:ny+1)::phi,phi_pre,chpoten
	  double precision,dimension(0:nx+1,0:ny+1)::psi,concentration
!----------------------------------------------------------------------c
      integer j, i

      pi=4.d0*atan(1.d0)
      
      do i=0,nx+1
         phi(i,ny) =1.0
         chpoten(i,ny) = (4.0*chpoten(i,ny-1)-chpoten(i,ny-2))/3.0
         phi(i,1) = 1.0
         chpoten(i,1) =(4.0*chpoten(i,2)-chpoten(i,3))/3.0

		 psi(i,ny)=0.d0!psi(i,ny-1)	
		 psi(i,1)=0.d0!psi(i,2)
      end do
      
      end subroutine      
!======================================================================c
!-------set macro vars  boundary condition
!----------------------------------------------------------------------c
      subroutine setMacroBC2D(uc,nx,ny,ux,uy,rho,pressure)
!-------input parms
      double precision uc
      integer nx, ny
      double precision,dimension(0:nx+1,0:ny+1)::ux,uy,rho,pressure
!----------------------------------------------------------------------c
      integer j, i

      !do j=1,ny
      !   ux(nx,j) = 0.d0
      !   uy(nx,j) = 0.d0
      !   ux(1,j) = 0.d0
      !   uy(1,j) = 0.d0
      !end do

      do i=0,nx+1
         ux(i,ny) = 0.00d0
         uy(i,ny) = 0.00d0
         ux(i,1) = -0.00d0
         uy(i,1) = -0.00d0
      end do
      end subroutine   
!======================================================================c
!-------save results
!----------------------------------------------------------------------c
      subroutine saveFlowField2D(ptime,ctime,uc,nx,ny,x,y,ux,uy,rho,pressure,phi,psi,concentration,phi_high,phi_low,chpoten,dh,surtenforcex,surtenforcey,preforcex,preforcey,bodyforcex,bodyforcey,nablaphix,nablaphiy,nablapsix,nablapsiy,nablaconx,nablacony,curv,intthick,visforcex,visforcey)
!-------input parms
      integer nx, ny
      double precision ptime,x(0:nx+1),y(0:ny+1),dh,ctime,uc,btux,btuy,btv,spcu,L2err,L2phi,L1err,L1phi,intthick,phi_high,phi_low,maxphi,minphi,uv2,intx,inty,sigma_T,temp
      double precision,dimension(0:nx+1,0:ny+1)::ux,uy,rho,pressure,surtenforcex,surtenforcey,preforcex,preforcey,bodyforcex,bodyforcey,psi,concentration,nablapsix,nablapsiy,nablaconx,nablacony,mfx,mfy,sfx,sfy,curv
      double precision,dimension(0:nx+1,0:ny+1)::phi,phi_0,nablaphix,nablaphiy,chpoten,visforcex,visforcey
!----------------------------------------------------------------------c
      integer i, j, ip, isp
      double precision pi
      character(len=64)::strtmp
      pi = 4.0d0*atan(1.d0)
      write(strtmp,888) ptime
888   format(f13.4)
      
      !spcu=0.d0
      !do i=1,nx
      !  do j=1,ny
      !      spcu=spcu+ux(i,j)**2+uy(i,j)**2
      !  end do        
      !end do
      !
      !write(*,*) 'spcu=',spcu
      
20    FORMAT(2(A,I6),A)
      
      open(2,file='log/flowfield_t'//trim(adjustl(strtmp))//'.plt')
      write(2,20) 'VARIABLES=X, Y, U, V, Phi, Psi, Psix, Psiy, Con, Conx, Cony, nablaCon, Nx, Ny, Pre, Vfx, Vfy'!
      write(2,*) 'ZONE T=FlowField I=',nx+2,', J=',ny,', F=POINT'

111   format(f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14f30.14)      
      
      do j=1,ny
         do i=0,nx+1
             intx=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
             inty=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
             temp=24.d0*phi(i,j)**2*(1.d0-phi(i,j))**2/intthick
             
             sigma_T=0.0!-0.2d0/(1.d0/0.6d0-concentration(i,j))
             mfx(i,j)=(nablaconx(i,j)-(intx*nablaconx(i,j)+inty*nablacony(i,j))*intx)*(-inty)*40.d0*sigma_T+(nablacony(i,j)-(intx*nablaconx(i,j)+inty*nablacony(i,j))*inty)*intx*40.d0*sigma_T
             mfy(i,j)=(nablaconx(i,j)-(intx*nablaconx(i,j)+inty*nablacony(i,j))*intx)*(-inty)*40.d0*sigma_T+(nablacony(i,j)-(intx*nablaconx(i,j)+inty*nablacony(i,j))*inty)*intx*40.d0*sigma_T
             sfx(i,j)=-curv(i,j)*(1.d0+0.20*log(1.d0-0.0d0*concentration(i,j)))*40.d0
             sfy(i,j)=-curv(i,j)*(1.d0+0.20*log(1.d0-0.0d0*concentration(i,j)))*40.d0

            write(2,111) x(i), y(j), ux(i,j), uy(i,j), phi(i,j), psi(i,j), nablapsix(i,j), nablapsiy(i,j), concentration(i,j), nablaconx(i,j), nablacony(i,j), sqrt(nablaconx(i,j)**2+nablacony(i,j)**2), intx, inty, pressure(i,j), visforcex(i,j), visforcey(i,j)
         end do
      end do

!      write(2,*) 'GEOMETRY'
!      write(2,*) 'F=POINT'
!      write(2,*) 'CS=GRID'
!      write(2,*) 'X=',100.0,'Y=',100.0
!      write(2,*) 'C=BLUE'
!      write(2,*) 'S=LOCAL'
!      write(2,*) 'L=SOLID'
!      write(2,*) 'PL=2'
!      write(2,*) 'LT=0.1'
!!      write(2,*) 'FC=CUST1,CLIPPING=CLIPTOVIEWPORT'
!      write(2,*) 'CLIPPING=CLIPTOVIEWPORT'
!      write(2,*) 'DRAWORDER=AFTERDATA'
!      write(2,*) 'MFC=""'
!      write(2,*) 'EP=72'
!      write(2,*) 'T=CIRCLE',40

      !L2err=0.d0
      !L2phi=0.d0
      !L1phi=0.d0
      !L1err=0.d0
      !maxphi=1.d0
      !minphi=0.d0
      !uv2=0.d0
      !do i=1,nx
      !    do j=1,ny
      !      phi_0(i,j)=0.50*(phi_low+phi_high)+0.50*(phi_high-phi_low)*tanh(2.d0*(real(nx-1)/4.0-sqrt((x(i)-real(nx-1)/2.0)**2+(y(j)-real(nx-1)/2.0)**2))/intthick)  
      !      !phi_0(i,j)=0.d0            
      !      !if(sqrt((x(i)-100.0)**2+(y(j)-100.0)**2).le.80.d0) then
      !      !    phi_0(i,j)=1.d0    
      !      !end if            
      !      !if(sqrt((x(i)-100.0)**2+(y(j)-100.0)**2).le.80.0.and.y(j).le.100.0.and.x(i).ge.92.0.and.x(i).le.108.0) then
      !      !    phi_0(i,j)=0.d0    
      !      !end if              
      !        !phi_0(i,j)=0.50*(0.0+1.0)+0.50*(1.0+0.0)*tanh(2.d0*(40.0-sqrt((x(i)-100.0)**2+(y(j)-150.0)**2))/3.0)
      !        L2phi=L2phi+phi_0(i,j)**2
      !        L2err=L2err+(phi(i,j)-phi_0(i,j))**2
      !        L1phi=L1phi+abs(phi_0(i,j))
      !        L1err=L1err+abs(phi(i,j)-phi_0(i,j))
      !        if(maxphi.lt.phi(i,j)) maxphi=phi(i,j)
      !        if(minphi.gt.phi(i,j)) minphi=phi(i,j)
      !        uv2=uv2+ux(i,j)**2+uy(i,j)**2
      !    end do          
      !end do
      !
      !L2err=sqrt(L2err/L2phi)
      !L1err=L1err/L1phi
      !uv2=sqrt(uv2/real(nx)/real(ny))
      !
      !write(*,*) 'L2err=',L2err,'L1err=',L1err
      !write(*,*) 'maxphi=',maxphi,'minphi=',minphi
      ! open(2,file='maxphi_Uerr.plt',access='append')
      ! write(2,111) ptime,uv2,maxphi,minphi,L2err
      ! close(2)        

    end subroutine
!=================================================================================================c
!	caculate  the  Heaviside function 
!-------------------------------------------------------------------------------------------------c
	function  HeavisideValue(Deltax)
	double precision HeavisideValue
	double precision  Deltax,pi
	
	pi=4.d0*atan(1.d0)
	
	if(Deltax.le.-2.d0) then
		HeavisideValue=0.d0
	else if(Deltax.ge.2.d0) then
		HeavisideValue=1.d0
    else
        HeavisideValue=0.5d0+Deltax/4.d0+2.d0/pi*sin(pi*Deltax/2.d0)
    end if
    
    end   
!=================================================================================================c
!	explotation 
!-------------------------------------------------------------------------------------------------c
subroutine explotation(nx,ny,nq,x,y,phi,psi,concentration,nablaphix,nablaphiy,dh,intthick,dt,ex,ey,wt,nablaconx,nablacony)
integer nx,ny,nq,ex(nq),ey(nq)
double precision x(0:nx+1),y(0:ny+1),dh,intthick,dt,signfun,cfl,wt(nq),intx,inty,xs,ys,dd,aa,bb,cc,ss1,ss2
double precision,dimension(0:nx+1,0:ny+1)::phi,psi,nablaphix,nablaphiy,nablaphi_x_y,nablaphi_x_x,nablaphi_y_y,nablaconx,nablacony
double precision,dimension(0:nx+1,0:ny+1)::eux,euy,concentration_new,concentration
integer ii,jj,kk,pp,qq,ll,ibeg,iend,jbeg,jend,ic,jc
double precision Con_int,tempArray(20),temp

!cfl=0.3
!do i=0,nx+1
!	do j=0,ny+1
!        intx=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
!        inty=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
!		eux(i,j)=signfun(phi(i,j)-0.5d0,intthick,dh)*intx*cfl
!		euy(i,j)=signfun(phi(i,j)-0.5d0,intthick,dh)*inty*cfl
!	end do
!end do
!
!
!do k=1,20
!    do i=2,nx-1
!	    do j=2,ny-1
!		    if(eux(i,j).ge.0.and.euy(i,j).ge.0) then
!			    concentration_new(i,j)=concentration(i,j)-eux(i,j)*dt/2.d0/dh*(3.d0*concentration(i,j)-4.d0*concentration(i-1,j)+concentration(i-2,j))-euy(i,j)*dt/2.d0/dh*(3.d0*concentration(i,j)-4.d0*concentration(i,j-1)+concentration(i,j-2))
!		    end if
!		
!		    if(eux(i,j).ge.0.and.euy(i,j).lt.0) then
!			    concentration_new(i,j)=concentration(i,j)-eux(i,j)*dt/2.d0/dh*(3.d0*concentration(i,j)-4.d0*concentration(i-1,j)+concentration(i-2,j))-euy(i,j)*dt/2.d0/dh*(-3.d0*concentration(i,j)+4.d0*concentration(i,j+1)-concentration(i,j+2))
!		    end if
!
!		    if(eux(i,j).lt.0.and.euy(i,j).ge.0) then
!			    concentration_new(i,j)=concentration(i,j)-eux(i,j)*dt/2.d0/dh*(-3.d0*concentration(i,j)+4.d0*concentration(i+1,j)-concentration(i+2,j))-euy(i,j)*dt/2.d0/dh*(3.d0*concentration(i,j)-4.d0*concentration(i,j-1)+concentration(i,j-2))
!		    end if
!
!		    if(eux(i,j).lt.0.and.euy(i,j).lt.0) then
!			    concentration_new(i,j)=concentration(i,j)-eux(i,j)*dt/2.d0/dh*(-3.d0*concentration(i,j)+4.d0*concentration(i+1,j)-concentration(i+2,j))-euy(i,j)*dt/2.d0/dh*(-3.d0*concentration(i,j)+4.d0*concentration(i,j+1)-concentration(i,j+2))
!		    end if		
!	    end do
!    end do

!do i=2,nx-1
!	    do j=2,ny-1
!            concentration(i,j)=concentration_new(i,j)
!	    end do
!    end do     
!end do


do i=2,nx-1
	    do j=2,ny-1
            concentration(i,j)=4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
            concentration(i,j)=psi(i,j)*concentration(i,j)/(concentration(i,j)**2+1e-4)
            if(phi(i,j).ge.5e-3.and.phi(i,j).le.0.995)  then
                intx=nablaphix(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                inty=nablaphiy(i,j)/sqrt(nablaphix(i,j)**2+nablaphiy(i,j)**2+1e-30)
                
                dd=intthick/4.d0*log(phi(i,j)/(1.d0-phi(i,j)))
                
                xs=x(i)-dd*intx
                ys=y(j)-dd*inty
                
                !write(*,*) 'xs=',xs,'ys=',ys, 'xy',sqrt((xs-80.d0)**2+(ys-80.d0)**2)

              ic=int(xs+1)
              jc=int(ys+1)
        
              ibeg=ic-2
              iend=ic+2
              jbeg=jc-2
              jend=jc+2
          
			    Con_int=0.d0
			    do ii=ibeg,iend
				    ll=ii-ibeg+1
				    tempArray(ll)=0.d0
				    do jj=jbeg,jend
					  
					    temp=concentration(ii,jj)
       
					    do qq=jbeg,jend
						    if(qq.ne.jj) then
							    temp=(ys-y(qq))/(y(jj)-y(qq))*temp
						    end if
					    end do
				    tempArray(ll)=tempArray(ll)+temp						
				    end do
			    end do
       
			    do ii=ibeg,iend
					    temp=tempArray(ii-ibeg+1)
					    do pp=ibeg,iend
						    if(pp.ne.ii) then
							    temp=(xs-x(pp))/(x(ii)-x(pp))*temp
						    end if
					    end do
				    Con_int=Con_int+temp
                end do
                concentration(i,j)=Con_int
            end if
	    end do
    end do     


    do j=1,ny
        do i=1,nx      
        nablaconx(i,j) =0.d0! (nablapsix(i,j)-lagx)*temp1/(temp1**2+1e-8)
        nablacony(i,j) =0.d0! (nablapsiy(i,j)-lagy)*temp1/(temp1**2+1e-8)
        end do
    end do
          
    do j=2,ny-1
        do i=2,nx-1                           
            do iq=1,nq
                nablaconx(i,j)=nablaconx(i,j)+wt(iq)*3.d0/dh*real(ex(iq))*concentration(i+ex(iq),j+ey(iq))
                nablacony(i,j)=nablacony(i,j)+wt(iq)*3.d0/dh*real(ey(iq))*concentration(i+ex(iq),j+ey(iq))
            end do
            nablaconx(i,j)=nablaconx(i,j)*4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
            nablacony(i,j)=nablacony(i,j)*4.d0*phi(i,j)*(1.d0-phi(i,j))/intthick
        end do
    end do   

    end subroutine
!=================================================================================================c
!	heavisidevalue 
!-------------------------------------------------------------------------------------------------c
function signfun(Deltax,intthick,dh)
double precision Deltax,intthick,signfun,dh,pi

pi=4.d0*atan(1.d0)

!heavisidevalue=tanh(2.d0*Deltax/intthick)
if(Deltax.lt.-dh) then
    signfun=-1.d0
else if(Deltax.gt.dh) then
    signfun=1.d0
else
    signfun=Deltax/dh+1.d0/pi*sin(pi*Deltax/dh)    
end if

end function

