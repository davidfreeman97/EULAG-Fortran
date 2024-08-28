      program moist_convec
c     Based on July 31st version of rce_crm_qmcl_shited_v2
c     Has single quantum state for each vertical column
c 8.7.24 - Changing syntax to use BLAS for speed


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Simple 2d periodic domain flow model
c
c  Author: W. Grabowski (grabow@ncar.ucar.edu) 
c          with help from P. Smolarkiewicz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc
cc Experiment: Try single mult operators with L96 
cc  nember of time steps:
cc    Was 20 in first run
      parameter(ntime0=500) ! 24/(4*60) hrs
      dimension precw(ntime0),dwtem(ntime0)
cc  major dimensions
cc HAVE TO BE CHANGED IN ALL ROUTINES:
      parameter(nx=4001,nz=51)
      parameter(nxz=nx*nz)
cc  grid
cc history tape: every 3hrs
      data dx,dz,dt,ntape /1.e3,500.,10.,720/ 
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      dimension xx(nx),zz(nz)

CC  MODEL VARIABLES
cc thermodynamics
      dimension theta(nx,nz),qv(nx,nz),qc(nx,nz),qr(nx,nz)
c      dimension qv_total(nx*nz, ntime0)
cc    Added in by DF for capturing thermo vars 
      dimension thermo_vars(8*nx*nz, ntime0)
cccccccccccccc
cc dynamics
      dimension ux(nx,nz),uz(nx,nz)     ! current time level
      dimension uxp(nx,nz),uzp(nx,nz)   ! previous time level
cc  forces for model variables
      dimension ft(nx,nz),fx(nx,nz),fz(nx,nz),
     *          fqv(nx,nz),fqc(nx,nz),fqr(nx,nz)
cc  advective velocities:
      dimension uxa(nx+1,nz),uza(nx,nz+1)

cc profiles  
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common /prof_a/ tau(nz)

cc required variables:
      dimension den(nx,nz),p(nx,nz),scr1(nx,nz),scr2(nx,nz)
      dimension clttemp(nx),atemp(nz),aqv(nz)


cc DF: For saving data 
      character(len=20) :: baseFileName
      character(len=30) :: fileName
cc constants
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

      fun23(a,b,c)=amax1(0.,amin1(1.,(b-a)/(b-c)))

c     DF: QMCl constants 
c      dimension phi_vec(500)
c      dimension B(500) 
c      CHARACTER training_text_file*13
c      CHARACTER eigenfunctions_text_file*18

c     Above are dummy vectors for storing data ported from matlab 
c     parameter(L = 50, ntime_training=25, D=nx*nz, N=500) 
c      parameter(L = 50, ntime_training=25, N_qm=500) 
      parameter(L = 50, ntime_training=5000, N_qm=5000) 
c      parameter(N_real = 500.0)
      integer D
      integer update_index
      real N_real
c      integer res_var

c     double precision training_vec(nxz*N)
      real training_vec(4001*51*5000*7)
      real T_bar_vec(nz*7)
      real T_bar(nz,7)
      real sigma_T(7)
c     ^^ 7 = number of distinct variables in training set
      real phi_vec(5000*50)
      integer i,ii,j,jj
c     ^For storing training data before reformatting 
      dimension S_matrix(L,L,4,nz)
c     DF: should have "res_var" in place of 4 to allow variation of parameter number 
c     DF: but currently throwing an error 
c     DF: Might no longer need U for now 
c     DF :Fixed for multi-var
      dimension U(L,L), trans_U(L,L)
c      dimension phi(number_training_pts,L)
      dimension phi(5000,L)
      dimension trans_phi(L,5000)

      dimension training_data_qm(nxz,N_qm, 7)
      dimension training_data_formatted(nx,nz,N_qm,7)
      dimension training_data_formatted_b(nx,7*nz,N_qm)

c      real rho(L,L) 
      dimension rho(L, nx)

cc    PCA stuff: 
      integer pca_res
      parameter(pca_res = 10)
      dimension training_pca(nx,pca_res, N_qm)
      real pca_coeff_mat_trans(pca_res, 4*nz)
      real pca_coeff_vec(4*nz*pca_res)

      real test_mat(3,2) 
      real test_vec(2,1)
      real test_vec_2(2,1) 

      test_mat = 0.0 
      test_vec = 0.0 
      test_mat(1,1) = 1.1 
      test_vec(1,1) = 1.1

      call dgemv('N',3,2,1.0,test_mat, 3
     &,test_vec,1,0.0,test_vec_2,1)

      training_pca = 0.0
cc    End PCA stuff

c     DF:print 
      print*, "Test Print"
c      call opngks
c      call gsclip(0)
c     D = nx*nz
      D = 7
      N_real = 5000.0
      update_index = 1
c     DF: Initialize quantum state rho 
c     DF: Fixed for multi-var
      do ii=1,nx 
      do i=1,L 
                 rho(i,ii) = 0 
      enddo
      rho(1,ii) = 1 

      enddo

      print*, "Created Rho"

c      rho = 0.0
c      rho(1,1) = 1.0

c     End of rho initialization 

c     DF: Import QMCl objects 
c      training_text_file = 'training2.txt'
c      eigenfunctions_text_file = 'eigenfunctions.txt'



c          open(unit=10, file = eigenfunctions_text_file, status = 'old')
c          do i =1,(L*N_qm)           
c          read(10,*) phi_vec(i)            
c         end do      
c          close(10)

      open(59,FILE="eigenfunctions.bin",FORM='UNFORMATTED',
     . ACTION='READ', STATUS='OLD', ACCESS='STREAM')
      read(59) phi_vec
      close(59)

          do i=1,L
              do j=1,N_qm
                  phi(j,i) = phi_vec(j + (i-1)*N_qm)
                  trans_phi(i,j) = phi_vec(j + (i-1)*N_qm)
              end do
          end do

      print*, "Loaded Eigenfunctions"

      open(77,FILE="T_bar.bin",FORM='UNFORMATTED',
     . ACTION='READ', STATUS='OLD', ACCESS='STREAM')
      read(77) T_bar_vec
      close(77)

          do i=1,nz
              do j=1,7
                  T_bar(i,j) = T_bar_vec(i + (j-1)*nz)
              end do
          end do

      print*, "Loaded T_bar"



      open(79,FILE="sigma_T.bin",FORM='UNFORMATTED',
     . ACTION='READ', STATUS='OLD', ACCESS='STREAM')
      read(79) sigma_T
      close(79)

      print*, "Loaded sigma_T"



      open(49,FILE="training_vec.bin",FORM='UNFORMATTED',
     . ACTION='READ', STATUS='OLD', ACCESS='STREAM')
c     training_vec_2 seems to be the working, updated course data
      read(49) training_vec
      close(49)

cc    PCA coeff initialization
      open(57,FILE="pca_coeff_vec.bin",FORM='UNFORMATTED',
     . ACTION='READ', STATUS='OLD', ACCESS='STREAM')
c     training_vec_2 seems to be the working, updated course data
      read(57) pca_coeff_vec
      close(57)
c           open(unit=12, file = training_text_file, status = 'old')
      do ii=1,pca_res
      do i=1,nz*4
            pca_coeff_mat_trans(ii,i) = pca_coeff_vec(
     .(ii-1)*nz*4 + i)
      end do
      end do

c     DF: Reformat training data to become a matrix 
c     DF: training data contains resolved and unresolved vars in 1st dim 
      do ii=1,7
      do i=1,N_qm
            do j=1,nx*nz
            training_data_qm(j,i,ii) = training_vec(
     .(ii-1)*nx*nz*N_qm + j +(i-1)*(nx*nz))
            end do
      end do
      end do

      do ii=1,7
      do i=1,N_qm
            do j=1,nx
            do jj=1,nz
            training_data_formatted(j,jj,i,ii) = training_vec(
     .(ii-1)*nx*nz*N_qm + j + (jj-1)*nx +(i-1)*(nx*nz))
            training_data_formatted_b(j, (7*jj)+ii, i) = 
     .training_data_formatted(j,jj,i,ii)
            enddo
            end do
      end do
      end do

      print*, "Set up training array"

c     DF:print
      print*, "Training data outputs"
      print*, training_data_qm(5,5,1)
      print*, training_data_qm(1,1,3)


c     End Import QMCl 
c     DF: Generate QMCl objects 
      CALL generate_qmcl_objects_qv(training_data_qm
     &, L, N_qm, nx, nz, 
     & phi, ntime_training, nx*nz, N_real, N_qm, S_matrix, U,
     & trans_U, pca_coeff_mat_trans, pca_res, training_pca)

c     DF:print 
      print*, "Generated QMCl Objects"
      print*, S_matrix(10,10,1,10)

c     End generate qmcl objects 
cc grid:
      time=0.
      dxi=1./dx
      dzi=1./dz
      dti=1./dt
      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo

      irstrt=1     ! 0-start from initial data, 1-start from tape

cc open tape to write to:
c      open(21,file='walker_bnch.d1-10',
c     .   form='unformatted',status='new')





cc open tape to read from:
      if(irstrt.eq.1)
     . open(31,file='walker_bnch.d1-10',
     .   form='unformatted',status='old')

cc    DF: Read and format training data  for QMCl 
                  
            
c      open(49,FILE="theta_qv_training.bin",FORM='UNFORMATTED',
c     . ACTION='READ', STATUS='OLD', ACCESS='STREAM')
c     read(49) training_vec
c      close(49)
c     Changed and moved this step up in the code
    

cc initialize moisture parameters:
      call moist_init

cc initialize model profiles:
      call prof_init

cc absorber (tau=1/time scale)
      zab = 16000.     ! height at which absorber begins
      t_abs=1000.       ! time scale
      towi=1./t_abs
      do k=1,nz
      tau(k)=towi*amax1(0.,zz(k)-zab)/(zz(nz)-zab)
      enddo

       do k=1,nz
       do i=1,nx
       den(i,k)=rho0(k)
       enddo
       enddo

            if(irstrt.eq.0) then



cc initial fields
       do k=1,nz
       do i=1,nx

        x1=xx(i)
        z1=zz(k)
        xc=xx(nx)/2.
        zc=zz(4)

        del=1.
        rad=sqrt((x1-xc)**2 + (z1-zc)**2)
        if(rad.gt.1.1e3) del=0.

       theta(i,k)=th_e(k) + 4.*del
       qv(i,k)=qv_e(k) + del*4.e-3
       qc(i,k)=0.
       qr(i,k)=0.

       ux(i,k)=ux_e(k)
       uz(i,k)=0.
       uxp(i,k)=ux_e(k)
       uzp(i,k)=0.

c        theta(i,k)=th_e(k) 
c        qv(i,k)=qv_e(k)
c        qc(i,k)=0.
c        qr(i,k)=0.

c        ux(i,k)=ux_e(k)
c        uz(i,k)=0.
c        uxp(i,k)=ux_e(k)
c        uzp(i,k)=0.

        ft(i,k)=0.
        fx(i,k)=0.
        fz(i,k)=0.

        fqv(i,k)=0.
        fqc(i,k)=0.
        fqr(i,k)=0.

        p(i,k)=0.
       enddo
       enddo
      
cc plot initial fields:
       call diagno_1(ux,uz,theta,nx,nz,scr1,scr2,rho0)
       call diagno_2(qv,qc,qr,nx,nz,scr1,scr2,rho0)
       call diagno_3(theta,qv,nx,nz,rho0,th_e,tm_e,a1,a2)
c       call plot_1(ux,uz,theta,nx,nz,xx,zz)
c       call plot_2(theta,qv,qc,qr,nx,nz,xx,zz,scr1,th_e,tm_e)

                   else

      ifl=2   ! number of files to read from tape
      call tape_rd(ifl,ux,uz,uxp,uzp,theta,qv,qc,qr,
     .             fx,fz,ft,fqv,fqc,fqr,p,nx,nz)
       do k=1,nz
       do i=1,nx
       den(i,k)=rho0(k)
       enddo
       enddo
                 endif


CCCC MARCH FORWARD IN TIME:
              ntime=ntime0
              do itime=1,ntime   ! TIME LOOP
               print*,'*** itime, time: ',itime,time
      print*, itime
cc extrapolate in time to get advective momentums:
       call velprd_1(ux,uxp,uxa,uz,uzp,uza,nx,nz,rho0)
          
cc save previous velocities:
       do i=1,nxz
        uxp(i,1)=ux(i,1)
        uzp(i,1)=uz(i,1)
       enddo

ccc surface flux
       call surfflux(theta,qv,ux,ft,fqv,th_e,nx,nz)

ccc radiative cooling:
       day=24.*3600.
       do k=1,nz
       do i=1,nx
       ft(i,k)=ft(i,k)-2.*1.5/day
       enddo
       enddo

ccc maintain wind:
       taunud=5.*24.*3600.   ! 5 day time scale
       do k=1,nz
       uxmean=0.
       do i=1,nx-1
       uxmean=uxmean+ux(i,k)/float(nx-1)
       enddo
       do i=1,nx
       fx(i,k)=fx(i,k) - 2.*(uxmean-ux_e(k))/taunud
       enddo
       enddo

cc add half of the force:
       do i=1,nxz
        theta(i,1)=theta(i,1)+.5*dt*ft(i,1)
        ux(i,1)   =   ux(i,1)+.5*dt*fx(i,1)
        uz(i,1)   =   uz(i,1)+.5*dt*fz(i,1)
        qv(i,1)   =   qv(i,1)+.5*dt*fqv(i,1)
        qc(i,1)   =   qc(i,1)+.5*dt*fqc(i,1)
        qr(i,1)   =   qr(i,1)+.5*dt*fqr(i,1)
       enddo

CC ADVECTION:
c liner: 1-iord=1, 0-iord prescribed inside mpdata
        liner=0
        if(itime/6*6.eq.itime) liner=1
c        if(itime/3*3.eq.itime) liner=1
c        if(itime/1*1.eq.itime) liner=1

cc advect velocities:
        call mpdat_2d(uxa,uza,   ux,den,1,liner)
        call mpdat_2d(uxa,uza,   uz,den,2,liner)
cc advect thermodynamic variables:
        call mpdat_2d(uxa,uza,theta,den,3,liner)
        call mpdat_2d(uxa,uza,   qv,den,4,liner)
        call mpdat_2d(uxa,uza,   qc,den,5,liner)
          call rain_fall(qr,tm_e,rho0,uza,nx,nz)
        call mpdat_2d(uxa,uza,   qr,den,6,liner)

cc save velocities after advection into advective velocities:
cc (used as scratch)
       do k=1,nz
       do i=1,nx
       uxa(i,k)=ux(i,k)
       uza(i,k)=uz(i,k)
       enddo
       enddo

cc finish thermodynamics
cc DF: calculation of forcings begins here for next cycle
c      open(unit = 92, file='therm_data_output')
c      write(92,*) theta, qv, qc, qr, ft, fqv, fqc, fqr
c      close(92)

cc    DF: writes qv output to big matrix, indexed by timestep
c      do k=1,nz
c      do i=1,nx 
c            qv_total(i+((k-1)*nx), itime) = qv(i,k)
c      enddo
c      enddo
      

       if (itime==1) then
      call thermo(theta,qv,qc,qr,ft,fqv,fqc,fqr)
       else
c     For now just using QC as the covariate input; change later
c     Still need to add code to read in training data and phi matrices

      call thermo_qmcl_qv(qv, qc, qr, theta, L
     & ,ux, uz, p, S_matrix, training_data_qm, 
     & phi, U, N_qm, 1.0, N_real, D, nx, nz, rho,
     & T_bar, sigma_T, update_index, trans_U, 
     & trans_phi, training_data_formatted_b, pca_res,
     & training_pca, pca_coeff_mat_trans)
       endif

       if(update_index.EQ.10) then
       update_index = 1
       else
       update_index = update_index + 1
       endif

c      do k=1,nz
c      do i=1,nx 
c            qv_total(i+((k-1)*nx), itime) = qv(i,k)
c      enddo
c      enddo
      
    
      
c      close(92)

cc DF: Replaced thermo with qmcl model
       
cc add absorber:
       if(zab.lt.zz(nz)) 
     1        call absor(ux,uz,theta,qv,qc,qr,ft,fqv,fqc,fqr)

cc add buoyancy
       epsb=rv/rg-1.
       do k=1,nz
       do i=1,nx
       scr1(i,k)=gg*( (theta(i,k)-th_e(k))/th0(k)
     *   + epsb*(qv(i,k)-qv_e(k))-qc(i,k)-qr(i,k) )
       enddo
       enddo
cc filter in vertical
       call integz(scr1,scr2,nx,nz)
c       call integxz(scr1,scr2,nx,nz)

cc apply
       do k=1,nz
       do i=1,nx
       uz(i,k) = uz(i,k)+.5*dt*scr1(i,k)
       enddo
       enddo

cc calculate pressure gradient force:
      epp=1.e-6
c      epp=1.e-7
      itp=100
      call gcrk_1(p,scr1,scr2,ux,uz,nx,nz,itp,epp)
      call prforc_1(p,scr1,scr2,ux,uz,nx,nz)
      do k=1,nz
      do i=1,nx
      ux(i,k)=scr1(i,k)
      uz(i,k)=scr2(i,k)
      enddo
      enddo

cc calculate velocity forces (using saved velocities after advection):
       do k=1,nz
       do i=1,nx
       fx(i,k)=(ux(i,k)-uxa(i,k))  *2./dt
       fz(i,k)=(uz(i,k)-uza(i,k))  *2./dt
       enddo
       enddo

cc update clock (in minutes...)
       time = time+dt/60. 

cc output and plot:
cc output and plot:
c       dtout=180.
       dtout=1.
       if(amod(time+.1*dt/60.,dtout).lt.0.5*dt/60.) then
cc plot selected fields:
c       call plot_1(ux,uz,theta,nx,nz,xx,zz)
c       call plot_2(theta,qv,qc,qr,nx,nz,xx,zz,scr1,th_e,tm_e)
cc analysis of output:
       print*,'   '
       call diagno_1(ux,uz,theta,nx,nz,scr1,scr2,rho0)
       call diagno_2(qv,qc,qr,nx,nz,scr1,scr2,rho0)
c       print*,' --- dwtemp, precw: ',a1,a2
       endif

cc  write data to tape:
        if(itime/ntape*ntape.eq.itime) then
      call tape_wr(ux,uz,uxp,uzp,theta,qv,qc,qr,
     .             fx,fz,ft,fqv,fqc,fqr,p,nx,nz)
       endif
      
             baseFileName = "qv_fromqm"

      write(fileName, '(A,I4,A)') baseFileName, itime, ".d1-10"

       open(93,file=
     .'/scratch/dfreeman/EULAG Fortran'//
     .'/VSCode Workspace'//
     .'/EULAG/Training Data 2/'//fileName,
     . form='unformatted',status='new')
         write(93) qv
       
        close(93)

             baseFileName = "qr_fromqm"

      write(fileName, '(A,I4,A)') baseFileName, itime, ".d1-10"

       open(101,file=
     .'/scratch/dfreeman/EULAG Fortran'//
     .'/VSCode Workspace'//
     .'/EULAG/Training Data 2/'//fileName,
     . form='unformatted',status='new')
         write(101) qr
       
        close(101)

                     baseFileName = "th_fromqm"

      write(fileName, '(A,I4,A)') baseFileName, itime, ".d1-10"

       open(103,file=
     .'/scratch/dfreeman/EULAG Fortran'//
     .'/VSCode Workspace'//
     .'/EULAG/Training Data 2/'//fileName,
     . form='unformatted',status='new')
         write(103) theta
       
        close(103)

                     baseFileName = "qc_fromqm"

      write(fileName, '(A,I4,A)') baseFileName, itime, ".d1-10"

       open(107,file=
     .'/scratch/dfreeman/EULAG Fortran'//
     .'/VSCode Workspace'//
     .'/EULAG/Training Data 2/'//fileName,
     . form='unformatted',status='new')
         write(107) qc
       
        close(107)




c       baseFileName = "qv_fromqm"

c       write(fileName, '(A,I4,A)') baseFileName, itime, ".d1-10"
 
c        open(93,file=
c      .'/scratch/dfreeman/EULAG Fortran'//
c      .'/VSCode Workspace'//
c      .'/EULAG/QM Output Data/'//fileName,
c      . form='unformatted',status='new')
c          write(93) qv
        
c         close(93)
 


      enddo


        stop
      end 


      subroutine velprd_1(ux,uxp,uxa,uz,uzp,uza,nx,nz,rho)
      dimension ux(nx,nz),uz(nx,nz)     
      dimension uxp(nx,nz),uzp(nx,nz)   
      dimension uxa(nx+1,nz),uza(nx,nz+1),rho(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

       do k=1,nz
       do i=2,nx
       uxa(i,k) =(0.75*(ux(i-1,k)+ux(i,k)) 
     .          - 0.25*(uxp(i-1,k)+uxp(i,k)) )*dt/dx * rho(k)
       enddo
cc cyclic in horizontal
       uxa(1,k) = uxa(nx,k)
       uxa(nx+1,k) = uxa(2,k)
       enddo
          
       do i=1,nx
       do k=2,nz
       uza(i,k) =(0.75*(uz(i,k-1)*rho(k-1)+uz(i,k)*rho(k)) 
     .          - 0.25*(uzp(i,k-1)*rho(k-1)+uzp(i,k)*rho(k)) )*dt/dz 
       enddo
cc zero flux in vertical
       uza(i,1) = - uza(i,2)
       uza(i,nz+1) = - uza(i,nz)
       enddo

       return
       end

      subroutine mpdat_2d(u1,u2,x,h,iflg,liner)
      parameter(nx=4001,nz=51)
      parameter(n1=nx+1,n2=nz+1)
      parameter(n1m=n1-1,n2m=n2-1)
      dimension u1(n1,n2m),u2(n1m,n2),x(n1m,n2m),h(n1m,n2m)
      common// v1(n1,n2m),v2(n1m,n2),f1(n1,n2m),f2(n1m,n2),
     *         cp(n1m,n2m),cn(n1m,n2m),
     *         mx(n1m,n2m),mn(n1m,n2m)
      real mx,mn
      parameter(iord0=2,isor=1,nonos=1,idiv=0)
      data ep/1.e-12/
c
      donor(y1,y2,a)=cvmgm(y2,y1,a)*a
      vdyf(x1,x2,a,r)=(abs(a)-a**2/r)*(abs(x2)-abs(x1))
     1                               /(abs(x2)+abs(x1)+ep)
      vcorr(a,b,y1,y2,r)=-0.125*a*b*y1/(y2*r)
      vcor31(a,x0,x1,x2,x3,r)= -(a -3.*abs(a)*a/r+2.*a**3/r**2)/3.
     1                         *(abs(x0)+abs(x3)-abs(x1)-abs(x2))
     2                         /(abs(x0)+abs(x3)+abs(x1)+abs(x2)+ep)
      vcor32(a,b,y1,y2,r)=0.25*b/r*(abs(a)-2.*a**2/r)*y1/y2
      vdiv1(a1,a2,a3,r)=0.25*a2*(a3-a1)/r
      vdiv2(a,b1,b2,b3,b4,r)=0.25*a*(b1+b2-b3-b4)/r
      pp(y)= amax1(0.,y)
      pn(y)=-amin1(0.,y)
 
      iord=iord0
      if(isor.eq.3) iord=max0(iord,3)
      if(liner.eq.1) iord=1

      do j=1,n2-1
        do i=1,n1
          v1(i,j) = u1(i,j)
        end do 
      end do 
      do i=1,n1-1
        do j=1,n2
          v2(i,j) = u2(i,j)
        end do 
      enddo

      if(nonos.eq.1) then
      do j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do i=1,n1m
      im=(i-1+(n1-i)/n1m*(n1-2))
      ip=(i+1    -i /n1m*(n1-2))
      mx(i,j)=amax1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
      mn(i,j)=amin1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
      end do
      end do
      endif
 
                         do 3 k=1,iord
 
      do 331 j=1,n2-1
      do 331 i=2,n1-1
  331 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      do j=1,n2-1
      f1(1 ,j)=f1(n1-1,j)
      f1(n1,j)=f1(2,j)
      enddo
      do 332 j=2,n2-1
      do 332 i=1,n1-1
  332 f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
      if (iflg.eq.6) then
        do i=1,n1m
          f2(i, 1)=donor(x(i,  1),x(i,  1),v2(i, 1))
          f2(i,n2)=donor(x(i,n2m),x(i,n2m),v2(i,n2))
        end do
      else
        do i=1,n1m
          f2(i, 1)=-f2(i,  2)
          f2(i,n2)=-f2(i,n2m)
        end do
      end if
  
      do 333 j=1,n2-1
      do 333 i=1,n1-1
  333 x(i,j)=x(i,j)-(f1(i+1,j)-f1(i,j)+f2(i,j+1)-f2(i,j))/h(i,j)
 
      if(k.eq.iord) go to 6

      do 49 j=1,n2-1
      do 49 i=1,n1
      f1(i,j)=v1(i,j)
   49 v1(i,j)=0.
      do 50 j=1,n2
      do 50 i=1,n1-1
      f2(i,j)=v2(i,j)
   50 v2(i,j)=0.
      do 51 j=2,n2-2
      do 51 i=2,n1-1
   51 v1(i,j)=vdyf(x(i-1,j),x(i,j),f1(i,j),.5*(h(i-1,j)+h(i,j)))
     *       +vcorr(f1(i,j), f2(i-1,j)+f2(i-1,j+1)+f2(i,j+1)+f2(i,j),
     *   abs(x(i-1,j+1))+abs(x(i,j+1))-abs(x(i-1,j-1))-abs(x(i,j-1)),
     *   abs(x(i-1,j+1))+abs(x(i,j+1))+abs(x(i-1,j-1))+abs(x(i,j-1))+ep,
     *                 .5*(h(i-1,j)+h(i,j)))
      if(idiv.eq.1) then
      do 511 j=2,n2-2
      do 511 i=2,n1-1
  511 v1(i,j)=v1(i,j)
     *    -vdiv1(f1(i-1,j),f1(i,j),f1(i+1,j),.5*(h(i-1,j)+h(i,j)))
     *    -vdiv2(f1(i,j),f2(i-1,j+1),f2(i,j+1),f2(i-1,j),f2(i,j),
     *                 .5*(h(i-1,j)+h(i,j)))
      endif
      do 52 j=2,n2-1
      do 52 i=2,n1-2
   52 v2(i,j)=vdyf(x(i,j-1),x(i,j),f2(i,j),.5*(h(i,j-1)+h(i,j)))
     *       +vcorr(f2(i,j), f1(i,j-1)+f1(i,j)+f1(i+1,j)+f1(i+1,j-1),
     *   abs(x(i+1,j-1))+abs(x(i+1,j))-abs(x(i-1,j-1))-abs(x(i-1,j)),
     *   abs(x(i+1,j-1))+abs(x(i+1,j))+abs(x(i-1,j-1))+abs(x(i-1,j))+ep,
     *                 .5*(h(i,j-1)+h(i,j)))
      i0=n1-2
      do j=2,n2-1
      v2(1,j)=vdyf(x(1,j-1),x(1,j),f2(1,j),.5*(h(1,j-1)+h(1,j)))
     *       +vcorr(f2(1,j), f1(1,j-1)+f1(1,j)+f1(2,j)+f1(2,j-1),
     *   abs(x(2,j-1))+abs(x(2,j))-abs(x(i0,j-1))-abs(x(i0,j)),
     *   abs(x(2,j-1))+abs(x(2,j))+abs(x(i0,j-1))+abs(x(i0,j))+ep,
     *                 .5*(h(1,j-1)+h(1,j)))
      v2(n1-1,j)=v2(1,j)
      enddo

      if(idiv.eq.1) then
      do 521 j=2,n2-1
      do 521 i=1,n1-1
  521 v2(i,j)=v2(i,j)
     *    -vdiv1(f2(i,j-1),f2(i,j),f2(i,j+1),.5*(h(i,j-1)+h(i,j)))
     *    -vdiv2(f2(i,j),f1(i+1,j),f1(i+1,j-1),f1(i,j-1),f1(i,j),
     *                 .5*(h(i,j-1)+h(i,j)))
      endif
      if(isor.eq.3) then
      do 61 j=2,n2-2
      do 61 i=3,n1-2
   61 v1(i,j)=v1(i,j)     +vcor31(f1(i,j),
     1        x(i-2,j),x(i-1,j),x(i,j),x(i+1,j),.5*(h(i-1,j)+h(i,j)))
      do j=2,n2-2
      v1(2,j)=v1(2,j)     +vcor31(f1(2,j),
     1        x(n1-2,j),x(1,j),x(2,j),x(3,j),.5*(h(1,j)+h(2,j)))
      v1(n1-1,j)=v1(n1-1,j) +vcor31(f1(n1-1,j),x(n1-3,j),x(n1-2,j),
     1                  x(n1-1,j),x(2,j),.5*(h(n1-2,j)+h(n1-1,j)))
      enddo
      do 62 j=2,n2-2
      do 62 i=2,n1-1
   62 v1(i,j)=v1(i,j)
     1 +vcor32(f1(i,j),f2(i-1,j)+f2(i-1,j+1)+f2(i,j+1)+f2(i,j),
     *   abs(x(i,j+1))-abs(x(i,j-1))-abs(x(i-1,j+1))+abs(x(i-1,j-1)),
     *   abs(x(i,j+1))+abs(x(i,j-1))+abs(x(i-1,j+1))+abs(x(i-1,j-1))+ep,
     *                   .5*(h(i-1,j)+h(i,j)))
      do 63 j=3,n2-2
      do 63 i=1,n1-1
   63 v2(i,j)=v2(i,j)     +vcor31(f2(i,j),
     1        x(i,j-2),x(i,j-1),x(i,j),x(i,j+1),.5*(h(i,j-1)+h(i,j)))
      do 64 j=3,n2-2
      do 64 i=2,n1-2
   64 v2(i,j)=v2(i,j)
     1 +vcor32(f2(i,j),f1(i,j-1)+f1(i+1,j-1)+f1(i+1,j)+f1(i,j),
     *   abs(x(i+1,j))-abs(x(i-1,j))-abs(x(i+1,j-1))+abs(x(i-1,j-1)),
     *   abs(x(i+1,j))+abs(x(i-1,j))+abs(x(i+1,j-1))+abs(x(i-1,j-1))+ep,
     *                   .5*(h(i,j-1)+h(i,j)))
      do 641 j=3,n2-2
      v2(1,j)=v2(1,j)
     1 +vcor32(f2(1,j),f1(1,j-1)+f1(2,j-1)+f1(2,j)+f1(1,j),
     *   abs(x(2,j))-abs(x(n1-2,j))-abs(x(2,j-1))+abs(x(n1-2,j-1)),
     *   abs(x(2,j))+abs(x(n1-2,j))+abs(x(2,j-1))+abs(x(n1-2,j-1))+ep,
     *                   .5*(h(1,j-1)+h(1,j)))
  641 v2(n1-1,j)=v2(1,j)
      endif
 
        do j=1,n2m
          v1( 1,j)=v1(n1m,j)
          v1(n1,j)=v1(  2,j)
        end do

      if (iflg.ne.6) then
        do i=1,n1m
          v2(i, 1)=-v2(i,  2)
          v2(i,n2)=-v2(i,n2m)
        end do
      end if

                  if(nonos.eq.1) then
c                 non-osscilatory option
      do 401 j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do 401 i=1,n1m
      im=(i-1+(n1-i)/n1m*(n1-2))
      ip=(i+1    -i /n1m*(n1-2))
      mx(i,j)=amax1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mx(i,j))
  401 mn(i,j)=amin1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mn(i,j))

      do 402 j=1,n2m 
      do 4021 i=2,n1-1
 4021 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      f1(1 ,j)=f1(n1m,j)
      f1(n1,j)=f1(2  ,j)
  402 continue
     
      do 403 i=1,n1m
      do 4031 j=2,n2m
 4031 f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
      if(iflg.ne.6) then
      f2(i, 1)=-f2(i,  2)
      f2(i,n2)=-f2(i,n2m)
      else
      f2(i, 1)=0.
      f2(i,n2)=0.
      endif
  403 continue

      do 404 j=1,n2m   
      do 404 i=1,n1m
      cp(i,j)=(mx(i,j)-x(i,j))*h(i,j)/
     1(pn(f1(i+1,j))+pp(f1(i,j))+pn(f2(i,j+1))+pp(f2(i,j))+ep)
      cn(i,j)=(x(i,j)-mn(i,j))*h(i,j)/
     1(pp(f1(i+1,j))+pn(f1(i,j))+pp(f2(i,j+1))+pn(f2(i,j))+ep)
  404 continue
      do 405 j=1,n2m 
      do 4051 i=2,n1m 
 4051 v1(i,j)=pp(v1(i,j))*
     1  ( amin1(1.,cp(i,j),cn(i-1,j))*pp(sign(1., x(i-1,j)))
     1   +amin1(1.,cp(i-1,j),cn(i,j))*pp(sign(1.,-x(i-1,j))) )
     2       -pn(v1(i,j))*
     2  ( amin1(1.,cp(i-1,j),cn(i,j))*pp(sign(1., x(i ,j )))
     2   +amin1(1.,cp(i,j),cn(i-1,j))*pp(sign(1.,-x(i ,j ))) )
      v1( 1,j)=v1(n1m,j)
      v1(n1,j)=v1( 2 ,j)
  405 continue

      do 406 i=1,n1m 
      do 406 j=2,n2m 
  406 v2(i,j)=pp(v2(i,j))*
     1  ( amin1(1.,cp(i,j),cn(i,j-1))*pp(sign(1., x(i,j-1)))
     1   +amin1(1.,cp(i,j-1),cn(i,j))*pp(sign(1.,-x(i,j-1))) )
     1       -pn(v2(i,j))*
     2  ( amin1(1.,cp(i,j-1),cn(i,j))*pp(sign(1., x(i ,j )))
     2   +amin1(1.,cp(i,j),cn(i,j-1))*pp(sign(1.,-x(i ,j ))) )
                  endif
    3                      continue
    6 continue
      return
      end   


      subroutine gcrk_1(p,pfx,pfz,u,w,n1,n3,itr,eps0)
      real p(*),pfx(*),pfz(*),u(*),w(*)
      parameter(nx=4001,nz=51)
      parameter(n=nx,l=nz)
      parameter(nn=n*l,nl=n*l)
      common// r(nn),qr(nn),ar(nn)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz)
      common/itero/ niter,nitsm,icount,eer,eem
      parameter (lord=3)
      dimension x(nn,lord),ax(nn,lord),ax2(lord),axar(lord),del(lord)
      dimension rho2d(nx,nz)
convergence test modes **************************************************
      logical ctest                                                     *
      data ctest/.false./                                               *
c     data ctest/.true./                                                *
      parameter (nplt=100)                                              *
      dimension err(0:nplt),xitr(0:nplt)                                *
      if(ctest) then                                                    *
      itr=6000/lord                                                     *
      ner=60                                                            *
      snorm=1./float(n*l)                                               *
      eps0=1.e-15                                                       *
      endif                                                             *
convergence test modes **************************************************
 
      do k=1,l
      do i=1,n
      rho2d(i,k)=rho0(k)
      enddo
      enddo

      eps=eps0*dti
      epa=1.e-30
      nlc=0
      itmn=5
cc iprc 0-no preconditioning, 1-with precon
      iprc=1

      call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,0)

      do k=1,nl
        r(k)=0.
       ar(k)=0.
       qr(k)=0.
      enddo
      do i=1,lord
       do k=1,nl
         x(k,i)=0.
        ax(k,i)=0.
       enddo
      enddo
      call prforc_1(p,pfx,pfz,u,w,n1,n3)
       call rhsdiv_1(pfx,pfz,rho2d,r,n1,n3,-1)
        call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1)
      eer0=0.
      eem0=-1.e15
      rl20=0.
      do k=1,nl
      eer0=eer0+qr(k)**2
      eem0=amax1(eem0,abs(qr(k)))
      rl20=rl20+r(k)**2
      enddo
      eer0=amax1(eer0,epa)
      eem0=amax1(eem0,epa)
      rl20=amax1(rl20,epa)
convergence test modes **************************************************
      if(ctest) then                                                    *
      do ier=0,nplt                                                     *
      err(ier)=eps                                                      *
      enddo                                                             *
      eer=-1.e15                                                        *
      do 3 k=1,nl                                                       *
    3 eer=amax1(eer,abs(r(k)))                                          *
      err(0)=eer                                                        *
      print 300,  err(0)                                                *
  300 format(4x,e11.4,' residual error at it=1')                        *
      endif                                                             *
convergence test modes **************************************************
       do k=1,nl
        x(k,1)=qr(k)
       enddo
      call laplc_1(x(1,1),ax(1,1),pfx,pfz,n1,n3)

      do 100 it=1,itr
       do i=1,lord
        ax2(i)=0.
        rax=0.
         do k=1,nl
          rax=rax+r(k)*ax(k,i)
          ax2(i)=ax2(i)+ax(k,i)*ax(k,i)
         enddo
        ax2(i)=amax1(epa,ax2(i))
        beta=-rax/ax2(i)
        dvmx=-1.e15
        rl2=0.
         do k=1,nl
          p(k)=p(k)+beta* x(k,i)
          r(k)=r(k)+beta*ax(k,i)
          dvmx=amax1(dvmx,abs(r(k)))
          rl2=rl2+r(k)*r(k)
         enddo
       if(dvmx.le.eps.and.it.ge.itmn) go to 200
       if(rl2.ge.rl20.and..not.ctest) go to 200
          rl20=amax1(rl2,epa)
       call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1)
       call laplc_1(qr,ar,pfx,pfz,n1,n3)
        nlc=nlc+1
         do ii=1,i
          axar(ii)=0.
           do k=1,nl
            axar(ii)=axar(ii)+ax(k,ii)*ar(k)
           enddo
          del(ii)=-axar(ii)/ax2(ii)
         enddo
        if(i.lt.lord) then
          do k=1,nl
            x(k,i+1)=qr(k)
           ax(k,i+1)=ar(k)
          enddo
           do ii=1,i
            do k=1,nl
              x(k,i+1)= x(k,i+1)+del(ii)* x(k,ii)
             ax(k,i+1)=ax(k,i+1)+del(ii)*ax(k,ii)
            enddo
           enddo
        else
          do k=1,nl
            x(k,1)=qr(k)+del(1)* x(k,1)
           ax(k,1)=ar(k)+del(1)*ax(k,1)
          enddo
           do ii=2,i
            do k=1,nl
              x(k,1 )= x(k,1)+del(ii)* x(k,ii)
              x(k,ii)=0.
             ax(k,1 )=ax(k,1)+del(ii)*ax(k,ii)
             ax(k,ii)=0.
            enddo
           enddo
        endif
convergence test modes **************************************************
      if(ctest) then                                                    *
      if(nlc/ner*ner.eq.nlc) then                                       *
      ier=nlc/ner                                                       *
      eer=-1.e15                                                        *
      do 50 k=1,nl                                                      *
   50 eer=amax1(eer,abs(r(k)))                                          *
      err(ier)=eer                                                      *
      endif                                                             *
      endif                                                             *
convergence test modes **************************************************
       enddo
  100 continue
  200 continue
      eer=0.
      eem=-1.e15
      do k=1,nl
      eer=eer+qr(k)**2
      eem=amax1(eem,abs(qr(k)))
      enddo
      eer=eer/eer0
      eem=eem/eem0
      niter=nlc
      nitsm=nitsm+niter
      icount=icount+1

convergence test modes **************************************************
      if(ctest) then                                                    *
      print 301, (err(ier),ier=1,nplt,1)                                *
  301 format(4x,5e11.4)                                                 *
      do 400 ier=0,nplt                                                 *
      xitr(ier)=ier*ner                                                 *
  400 err(ier)=alog10(err(ier)*dt )                                     *
      plmx=float(itr*lord)                                              *
c      call set(.1,.9,.1,.9,0.,plmx,-10.,0.,1)                           *
c      call labmod('(i4)','(f5.0)',4,4,2,2,20,20,0)                      *
c      call periml(4,10,5,2)                                             *
c      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)                         *
c      call curved(xitr,err,nplt+1)                                      *
c      i1=int(102.4+409.6)                                               *
c      call wtstr(cpux(i1),cpuy(50),'niter',3,0,0)                       *
c      call wtstr(cpux(17),cpuy(i1),'log e',3,90,0)                      *
c      call frame                                                        *
      endif                                                             *
convergence test modes **************************************************
      return
      end

      subroutine precon_1(rhs,p,r,c11,c33,d,iflg,jfl)
      parameter(nx=4001,nz=51)
      parameter(n=nx,l=nz,nl=nx*nz)
      dimension rhs(n,l),p(n,l),r(n,l),
     .          c11(n,l),c33(n,l),d(n,l)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
compute available storage 
      dimension e(n,0:l-1),f(n,0:l-1)
     .       ,px(n,l),dgh(n,l),po(n,l)

      data beta/-1.e15/
      data itr,line/2,1/

      if(iflg.eq.0) then
       do i=1,nl
        p(i,1)=rhs(i,1)
       enddo
      return
      endif
    
      omg=.7
      oms=1.-omg
      dxi2=0.25*dxi*dxi
      dzi2=0.25*dzi*dzi
      do i=1,nl
       c33(i,1)=d(i,1)*dzi2
       c11(i,1)=d(i,1)*dxi2
       dgh(i,1)=0.
        po(i,1)=0.
         p(i,1)=0.
         r(i,1)=0.
      enddo
      if(line.eq.1) then
       do k=1,l
         do i=2,n-1
          dgh(i,k)=c11(i+1,k)+c11(i-1,k)
         enddo
          dgh(1,k)=c11(2,k)+c11(n-1,k)
          dgh(n,k)=c11(2,k)+c11(n-1,k)
       enddo
      endif

      if(jfl.eq.0) then
      if(line.eq.0) then
      beta=-1.e15
      do i=1,nl
      beta=amax1(beta,abs(c11(i,1))/d(i,1))
      enddo
      beta=0.5/beta
      else
      beta=1.
      endif
      return
      endif
      beti=1./beta*(1-line)

      do 100 it=1,itr
      do i=1,nl
       r(i,1)=r(i,1)+d(i,1)*(beti*p(i,1)-rhs(i,1))
     .                  +dgh(i,1)*p(i,1)
      enddo
      do i=1,n
       e(i,0)=1.
       f(i,0)=0.
       dn=d(i,1)*beti+2.*c33(i,2)+dgh(i,1)
       e(i,1)=2.*c33(i,2)/dn
       f(i,1)=     r(i,1)/dn
      enddo
      do k=2,l-1
       do i=1,n
        dn=c33(i,k+1)+c33(i,k-1)*(1.-e(i,k-2))+d(i,k)*beti
     .                                + dgh(i,k)
        e(i,k)=                      c33(i,k+1)/dn
        f(i,k)=(c33(i,k-1)*f(i,k-2)+r(i,k))/dn
       enddo
      enddo
       do i=1,n
        dn=d(i,l)*beti+2.*(1.-e(i,l-2))*c33(i,l-1)
     .                                + dgh(i,l)
        p(i,l)=(r(i,l)+2.*f(i,l-2)*c33(i,l-1))/dn
        p(i,l-1)=f(i,l-1)/(1.-e(i,l-1))
       enddo
      do k=l-2,1,-1
       do i=1,n
        p(i,k)=e(i,k)*p(i,k+2)+f(i,k)
       enddo
      enddo


      if(line.eq.1) then
       do i=1,nl
        p(i,1)=oms*po(i,1)+omg*p(i,1)
       po(i,1)=     p(i,1)
       enddo
      endif

      if(it.eq.itr) go to 101
      do k=1,l
      do i=2,n-1
      px(i,k)=c11(i,k)*(p(i+1,k)-p(i-1,k))
      enddo
      px(1,k)=c11(1,k)*(p(2,k)-p(n-1,k))
      px(n,k)=c11(n,k)*(p(2,k)-p(n-1,k))
      enddo

      do k=1,l
      do 91 i=2,n-1
   91 r(i,k)=px(i+1,k)-px(i-1,k)
      r(1,k)=(px(2,k)-px(n-1,k))
      r(n,k)=(px(2,k)-px(n-1,k))
      enddo

  100 continue
  101 continue

      return
      end


      subroutine laplc_1(p,r,u,w,n1,l1)
      dimension p(n1,l1),r(n1,l1),u(n1,l1),w(n1,l1)
      parameter(nx=4001,nz=51)
      parameter(n=nx,l=nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz)
compute available storage in //
      parameter (nl=n*l)
      dimension px(n,l),pz(n,l)

      dxil=.5*dxi
      dzil=.5*dzi

compute pressure derivatives everywhere
      do 18 k=1,l
      do 1 i=2,n-1
    1 px(i,k)=     dxil*(p(i+1,k)-p(i-1,k))
      px(1,k)=dxil*(p(2,k)-p(n-1,k))
      px(n,k)=dxil*(p(2,k)-p(n-1,k))
   18 continue
      do 38 i=1,n
      do 3 k=2,l-1
    3 pz(i,k)=dzil*(p(i,k+1)-p(i,k-1))
      pz(i,1)= dzi*(p(i,2)-p(i,1))
   38 pz(i,l)= dzi*(p(i,l)-p(i,l-1))

compute interior pressure forces
      do 21 i=1,n
      do 10 k=2,l-1
      u(i,k)=-px(i,k)
   10 w(i,k)=-pz(i,k)
      w(i,1)=0.
      w(i,l)=0.
      u(i,1)=-px(i,1)
   21 u(i,l)=-px(i,l)

      do 99 k=1,l
      coef=rho0(k)
      do 99 i=1,n
      u(i,k)=coef*u(i,k)
   99 w(i,k)=coef*w(i,k)

compute laplacian
      do 911 k=1,l
      do 91  i=2,n-1
   91 r(i,k)=dxil*(u(i+1,k)-u(i-1,k))
      r(1,k)=dxil*(u(2,k)-u(n-1,k))
      r(n,k)=dxil*(u(2,k)-u(n-1,k))
  911 continue
      do 931 i=1,n
      do 93  k=2,l-1
   93 r(i,k)=r(i,k)+dzil*(w(i,k+1)-w(i,k-1))
      r(i,1)=r(i,1)+dzi *(w(i,2)-w(i,1)) 
      r(i,l)=r(i,l)+dzi *(w(i,l)-w(i,l-1))
  931 continue
      do 94 i=1,n
      do 94 k=1,l
   94 r(i,k)=-r(i,k)/rho0(k)

      return
      end

      subroutine prforc_1(p,pfx,pfz,u,w,n1,n3)
      dimension p(n1,n3),u(n1,n3),w(n1,n3),pfx(n1,n3),pfz(n1,n3)
      parameter(nx=4001,nz=51)
      parameter(n=nx,l=nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      parameter (nl=n*l)
      dimension px(n,l),pz(n,l)

      dxil=.5*dxi
      dzil=.5*dzi

compute pressure derivatives everywhere
      do 18 k=1,l
      do 1 i=2,n-1
    1 px(i,k)=     dxil*(p(i+1,k)-p(i-1,k))
      px(1,k)=dxil*(p(2,k)-p(n-1,k))
      px(n,k)=dxil*(p(2,k)-p(n-1,k))
   18 continue
      do 38 i=1,n
      do 3 k=2,l-1
    3 pz(i,k)=dzil*(p(i,k+1)-p(i,k-1))
      pz(i,1)= dzi*(p(i,2)-p(i,1))
   38 pz(i,l)= dzi*(p(i,l)-p(i,l-1))

compute interior pressure forces
      do 21 i=1,n
      do 10 k=2,l-1
      pfx(i,k)=u(i,k)-px(i,k)
   10 pfz(i,k)=w(i,k)-pz(i,k)
      pfz(i,1)=0.
      pfz(i,l)=0.
      pfx(i,1)=u(i,1)-px(i,1)
   21 pfx(i,l)=u(i,l)-px(i,l)

      return
      end

      subroutine rhsdiv_1(u,w,d,r,n1,l1,iflg)
      dimension u(n1,l1),w(n1,l1),d(n1,l1),r(n1,l1)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

      n=n1
      l=l1
      nl=n*l

      do 200 i=1,nl
  200 r(i,1)=0.

      dxil=.5*dxi
      dzil=.5*dzi

      do 1 i=2,n-1
      do 1 j=1,l
    1 r(i,j)=dxil*(u(i+1,j)*d(i+1,j)-u(i-1,j)*d(i-1,j))
      do 11 j=1,l
      r(1,j)=dxil*(u(2,j)*d(2,j)-u(n-1,j)*d(n-1,j))
      r(n,j)=dxil*(u(2,j)*d(2,j)-u(n-1,j)*d(n-1,j))
   11 continue
      do 3 k=2,l-1
      do 3 i=1,n
    3 r(i,k)=r(i,k)
     3        +dzil*(w(i,k+1)*d(i,k+1)-w(i,k-1)*d(i,k-1))
      do 13 i=1,n
      r(i,1)=r(i,1)+dzi*(w(i,2)*d(i,2)-w(i,1)*d(i,1))
   13 r(i,l)=r(i,l)+dzi*(w(i,l)*d(i,l)-w(i,l-1)*d(i,l-1))

      if(iflg.ne.0) then
      do 4 i=1,nl
    4 r(i,1)=iflg*r(i,1)/d(i,1)
      endif

      return
      end

       subroutine minmax(a,n,an,ax)
       dimension a(n)
       an= 1.e15
       ax=-1.e15
       do i=1,n
       an=amin1(a(i),an)
       ax=amax1(a(i),ax)
       enddo
       return
       end

      subroutine diagno_1(ux,uz,th,nx,nz,scr1,scr2,rho)
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),scr1(nx,nz),scr2(nx,nz)
      dimension rho(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
cc pressure solver diagnostics
      common/itero/ niter,nitsm,icount,eer,eem

      nxz=nx*nz
      do i=1,nx
      do k=1,nz
      scr2(i,k)=rho(k)
      enddo
      enddo

      print 200, time
 200  format(1x,' ****** analysis for time (min): ',f8.2)

      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max ux: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max uz: ',2e12.4)

      cour=0.
      do i=1,nxz
      cour=amax1(cour,abs(ux(i,1))*dt/dx+abs(uz(i,1))*dt/dz)
      enddo
      print 302,cour
 302  format(1x,' --> max courant: ',e12.4)

      call minmax(th,nxz,amn,amx)
      print 203,amn,amx
 203  format(1x,' --> min, max th: ',2e12.4)
      call rhsdiv_1(ux,uz,scr2,scr1,nx,nz,1)

      call minmax(scr1,nxz,amn,amx)
      print 204,amn,amx
 204  format(1x,' --> min, max div: ',2e12.4)

      nitav=nitsm/max0(icount,1)
      print 205, eer,eem,niter,nitav
  205 format(1x,'            eer,eem:',2e11.4/
     .       1x,'       niter, nitav:',2i4)

       if(cour.gt.1.) then
c       call clsgks
       stop 'courant'
       endif

       return
       end

      subroutine diagno_2(ux,uz,th,nx,nz,scr1,scr2,rho)
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),scr1(nx,nz),scr2(nx,nz)
      dimension rho(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
cc pressure solver diagnostics
      common/itero/ niter,nitsm,icount,eer,eem

      nxz=nx*nz
      do i=1,nx
      do k=1,nz
      scr2(i,k)=rho(k)
      enddo
      enddo

      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max qv: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max qc: ',2e12.4)

      call minmax(th,nxz,amn,amx)
      print 203,amn,amx
 203  format(1x,' --> min, max qr: ',2e12.4)

       return
       end

      subroutine diagno_3(th,qv,nx,nz,rho,th_e,tm_e,dwt,prew)
      dimension th(nx,nz),qv(nx,nz)
      dimension rho(nz),th_e(nz),tm_e(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

cc mean density weighted temp:
       sumg=0.
       do i=1,nx-1
       sum1=0.5*rho(1)   ! dz cancells out!
       sum2=0.5*rho(1)*th(i,1)*tm_e(1)/th_e(1)
       do k=2,nz
       sum1=sum1+rho(k)   
       sum2=sum2+rho(k)*th(i,k)*tm_e(k)/th_e(k)
       enddo
       sumg=sumg + sum2/sum1/float(nx-1)
       enddo
       dwt=sumg

cc mean precip water:
       sumg=0.
       do i=1,nx-1
       sum2=0.5*rho(1)*qv(i,1)*dz
       do k=2,nz
       sum2=sum2+rho(k)*qv(i,k)*dz
       enddo
       sumg=sumg + sum2/float(nx-1)
       enddo
       prew=sumg

       return
       end


      subroutine moist_init

create parametrs for the model:

cc mass, terminal velocity, diameter
      data ar,br,cr,dr /5.2e2,3.,130.,0.5/
      data as,bs,cs,ds /2.5e-2,2.,4.,0.25/
cc collection ef., alpha, beta
      data er,alphr,betr /0.8, 1., 2./
      data es,alphs,bets /0.2, .3, 3./
cc No
      data anor,anos /2*1.e7/
cc latent heats:
      data hlatv,hlats /2.53e6,2.84e6/
cc cloud droplet concentration (per cc)
      data dconc /200./ ! must be between 50 and 2000
cc limiting temperatures
      data tup,tdn /273.,263./   
cc gammas:
      data gamb1r,gambd1r /6.0,11.7/
      data gamb1s,gambd1s /2.0,2.56/
cc reference temperature and saturated vapor pressure:
      data tt0,ee0 /273.16,611./

      common/rain_p0/ ar,br,cr,dr,er,alphr,betr,gamb1r,gambd1r,anor
      common/rain_p1/ dconc,ddisp
      common/snow_p0/ as,bs,cs,ds,es,alphs,bets,gamb1s,gambd1s,anos
      common/temp_p/ tup,tdn
      common/latent/hlatv,hlats
      common/reference/ tt0,ee0

      common /const/ gg,cp,rg,rv
      data gg,cp,rg,rv   /9.72,1005.,287.,461./

cc check consistency
      if(dconc.lt.50..or.dconc.gt.2000.) then
      print*,' *** inconsistent droplet concentration. stop.'
      stop 'dconc'
      endif
cc calculate relative dispersion for Berry's autoconversion:
      ddisp=0.146-5.964e-2*alog(dconc/2000.)
      print 2075,anor,anos,dconc,ddisp
2075  format(1x,' N0 in raindrop distr.: ',e15.3/
     1 1x,' N0 in snowflake distr.: ',e15.3/
     1 1x,' Berry parameters of cloud droplet spectrum:'/
     1 1x,' droplet conc., relative disp.: ',2f12.4)

       return
       end

      subroutine rain_fall(qr,tm_e,rho,uza,nx,nz)
cc modify vertical advective velocity for rain fallout
      dimension qr(nx,nz),uza(nx,nz+1),tm_e(nz),rho(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common/rain_p0/ ar,br,cr,dr,er,alphr,betr,gamb1r,gambd1r,anor
      common/snow_p0/ as,bs,cs,ds,es,alphs,bets,gamb1s,gambd1s,anos
      common/temp_p/ tup,tdn

      parameter(nx1=4001)
      common /sprec/ precip(nx1)

      real lambdr,lambds,massr,masss

cc statement functions:
      alim01(fi)=amax1(0.,amin1(1.,fi))
      comb(tm,td,tu,ad,au)=
     1  alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad

      if(nx1.ne.nx) stop 'dim in rain_fall'
      gc3=dt*dzi

      do k=2,nz
        do i=1,nx
               dens= 0.5*(rho(k)+rho(k-1))
               qrv = 0.5*(qr(i,k)+ qr(i,k-1))
                coe_l=comb(tm_e(k),tdn,tup,0.,1.)   ! liquid part
                 qpr=qrv*coe_l         ! divide between rain and snow
                 qps=qrv-qpr           ! divide between rain and snow
         lambdr=(ar*anor*gamb1r/dens/(qpr+1.e-6))**(1./(1.+br)) 
         lambds=(as*anos*gamb1s/dens/(qps+1.e-6))**(1./(1.+bs)) 

        vtr=cr*gambd1r/gamb1r / lambdr**dr  ! terminal velocity
        vts=cs*gambd1s/gamb1s / lambds**ds  ! terminal velocity

               vtf=coe_l*vtr+(1.-coe_l)*vts   ! TERMINAL VELOCITY

               uza(i,k)=uza(i,k)-vtf*dens*gc3
        end do
      end do
ccc
CC LOWER AND UPPER BOUNDARIES:
cc lower:
        do i=1,nx
               dens=1.5*rho(1)-.5*rho(2)
               qrv =amax1(0.,1.5*qr(i,1)-.5*qr(i,2))
                coe_l=comb(tm_e(1),tdn,tup,0.,1.)   ! liquid part
                 qpr=qrv*coe_l         ! divide between rain and snow
                 qps=qrv-qpr           ! divide between rain and snow
         lambdr=(ar*anor*gamb1r/dens/(qpr+1.e-6))**(1./(1.+br))
         lambds=(as*anos*gamb1s/dens/(qps+1.e-6))**(1./(1.+bs))

        vtr=cr*gambd1r/gamb1r / lambdr**dr  ! terminal velocity
        vts=cs*gambd1s/gamb1s / lambds**ds  ! terminal velocity

               vtf=coe_l*vtr+(1.-coe_l)*vts   ! TERMINAL VELOCITY
cc
           precip(i)=dens*vtf*qrv*3600. ! in mm/hr
cc

               uza(i,1)=uza(i,1)-vtf*dens*gc3
        end do
cc upper:
        do i=1,nx
               dens=1.5*rho(nz)-.5*rho(nz-1)
               qrv =amax1(0.,1.5*qr(nz,1)-.5*qr(nz-1,2))
                coe_l=comb(tm_e(nz),tdn,tup,0.,1.)   ! liquid part
                 qpr=qrv*coe_l         ! divide between rain and snow
                 qps=qrv-qpr           ! divide between rain and snow
         lambdr=(ar*anor*gamb1r/dens/(qpr+1.e-6))**(1./(1.+br))
         lambds=(as*anos*gamb1s/dens/(qps+1.e-6))**(1./(1.+bs))

        vtr=cr*gambd1r/gamb1r / lambdr**dr  ! terminal velocity
        vts=cs*gambd1s/gamb1s / lambds**ds  ! terminal velocity

               vtf=coe_l*vtr+(1.-coe_l)*vts   ! TERMINAL VELOCITY

               uza(i,nz+1)=uza(i,nz+1)-vtf*dens*gc3
        end do

        return
        end
cc Potentially  use QMCl to replace this "thermo" subroutine 
cc This wouldn't involve worrying about coarsening 
cc
      subroutine thermo(th,qv,qc,qr,fth,fqv,fqc,fqr)
      parameter(nx=4001,nz=51)
      parameter(n=nx,l=nz)
      dimension th(n,l),qv(n,l),qc(n,l),qr(n,l),fth(n,l),
     1         fqv(n,l),fqc(n,l),fqr(n,l)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /const/ gg,cp,rg,rv

      common/rain_p0/ ar,br,cr,dr,er,alphr,betr,gamb1r,gambd1r,anor
      common/rain_p1/ dconc,ddisp
      common/snow_p0/ as,bs,cs,ds,es,alphs,bets,gamb1s,gambd1s,anos
      common/temp_p/ tup,tdn
      common/latent/hlatv,hlats
      common/reference/ tt0,ee0
      real lambdr,lambds,massr,masss

cc statement functions:
      alim01(fi)=amax1(0.,amin1(1.,fi))
      comb(tm,td,tu,ad,au)=
     1  alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad

condensation/evaporation

      pi=3.1416
      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

      do 100 k=1,l
      thetme=th_e(k)/tm_e(k)
      do 100 i=1,n
      coe_l=comb(tm_e(k),tdn,tup,0.,1.)   ! liquid contribution
      pre=1.e5*thetme**e
      tt=th(i,k)/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      esi=ee0*exp(b * delt)
      qvsw=a * esw /(pre-esw)
      qvsi=a * esi /(pre-esi)
      qvs=coe_l*qvsw + (1.-coe_l)*qvsi
ccc linearized condensation rate is next:
      cf1=thetme/th(i,k)
      cf1=cf1*cf1
      cf1=c*cf1*pre/(pre-esw)*d
      delta=(qv(i,k)-qvs)/(1.+qvs*cf1)
c--->
ccc one Newton-Raphson iteration is next:
      thn=th(i,k)+c*thetme*delta
      tt=thn/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      esi=ee0*exp(b * delt)
      qvsw=a * esw /(pre-esw)
      qvsi=a * esi /(pre-esi)
      qvs=coe_l*qvsw + (1.-coe_l)*qvsi
      fff=qv(i,k)-delta-qvs
      cf1=thetme/thn
      cf1=cf1*cf1
      cf1=c*cf1*pre/(pre-esw)*d
      fffp=-1.-qvs*cf1
      delta=delta-fff/fffp
ccc end of the iteration; if required, it can be repeated
c--->
      delta=amin1( qv(i,k), amax1(-qc(i,k),delta) )
      qv(i,k)=qv(i,k)-delta
      qc(i,k)=qc(i,k)+delta
      th(i,k)=th(i,k)+c*thetme*delta
      delta=amin1( qv(i,k), amax1(-qc(i,k),delta) )
      fqv(i,k)=-delta*2.*dti
      fth(i,k)=-c*thetme*fqv(i,k)
      fqc(i,k)=-fqv(i,k)
  100 continue

ccc remove trace of water variables:
      nl=n*l
      do i=1,nl
      qc(i,1)=cvmgm(0.,qc(i,1),qc(i,1)-1.e-9)
      qr(i,1)=cvmgm(0.,qr(i,1),qr(i,1)-1.e-10)
      enddo


compute moist forces update
      do 300 k=1,l
      thetme=th_e(k)/tm_e(k)
      do 300 i=1,n
       tt=th(i,k)/thetme
       pre=1.e5*thetme**e
       coe_l=comb(tm_e(k),tdn,tup,0.,1.)   ! liquid contribution
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      esi=ee0*exp(b * delt)
      qvsw=a * esw /(pre-esw)
      qvsi=a * esi /(pre-esi)

       ssw=qv(i,k) / qvsw      ! saturation ratio
       ssi=qv(i,k) / qvsi      ! saturation ratio

      qpr=qr(i,k)*coe_l                ! divide between rain and snow
      qps=qr(i,k)-qpr                  ! divide between rain and snow
      qcc=qc(i,k)*coe_l                ! divide between ice and water
      qci=qc(i,k)-qcc                  ! divide between ice and water

      lambdr=(ar*anor*gamb1r/rho0(k)/(qpr+1.e-6))**(1./(1.+br)) ! lambda
      lambds=(as*anos*gamb1s/rho0(k)/(qps+1.e-6))**(1./(1.+bs)) ! lambda

CC AUTOCONVERSION:
cc rain - Berry:
      del2=1.e3*rho0(k)*qcc
      autc=1./rho0(k)*1.67e-5*del2*del2 /
     1 (5. + .0366*dconc/(ddisp*(del2+1.E-6)))
cc snow:
       tc=tt-273.16
       times=amin1(1.e3,(3.56*tc+106.7)*tc+1.e3) ! time scale for
      auti=qci/times
      AUT = autc + auti

CC GROWTH:
      conr=anor/lambdr ! concentration
      cons=anos/lambds ! concentration

      massr=rho0(k)*(qpr+1.e-7) / conr  ! mass
      masss=rho0(k)*(qps+1.e-7) / cons  ! mass

      diamr=(massr/ar)**(1./br) ! diameter
      diams=(masss/as)**(1./bs) ! diameter

      rer=cr*diamr**(dr+1.)/2.e-5  ! Reynolds number
      res=cs*diams**(ds+1.)/2.e-5  ! Reynolds number

      ventr=amax1(1.,.78+.27*sqrt(rer))  ! ventilation factor
      vents=amax1(1.,.65+.39*sqrt(res))  ! ventilation factor

      thfun=1.e-7/(2.2*tm_e(k)/esw+2.2e-2/tm_e(k))  ! thermodynamic fun.

      g_acc_r=pi/4.*cr*diamr**(2.+dr)*er*alphr*rho0(k)*qc(i,k) ! growth
      g_acc_s=pi/4.*cs*diams**(2.+ds)*es*alphs*rho0(k)*qc(i,k) ! growth

      g_dep_r=4.*pi*diamr/betr*(ssw-1.)*ventr*thfun   ! growth/evap
      g_dep_s=4.*pi*diams/bets*(ssi-1.)*vents*thfun   ! growth/evap

      acc_r=conr * g_acc_r * qpr / (qpr + 1.e-9)
      acc_s=cons * g_acc_s * qps / (qps + 1.e-9)

       ACC= acc_r + acc_s  ! growth by accretion

      dep_r=conr * g_dep_r * qpr / (qpr + 1.e-9)
      dep_s=cons * g_dep_s * qps / (qps + 1.e-9)

       DEP= dep_r + dep_s  ! growth by deposition

      dcol=2.*(AUT + ACC)
      dcol=amin1(dcol,  2.*dti*qc(i,k)+fqc(i,k))
      devp=2.*DEP
      devp=amax1(devp, -2.*dti*qr(i,k)-dcol)
cc
      fqr(i,k)=devp+dcol
      fqc(i,k)=fqc(i,k)-dcol
      fqv(i,k)=fqv(i,k)-devp
      fth(i,k)=fth(i,k)+c*devp*thetme
  300 continue

      return
      end
cc    End of thermo subroutine, I believe  

      subroutine absor(ux,uz,th,qv,qc,qr,fth,fqv,fqc,fqr)
      parameter(nx=4001,nz=51)
      parameter(n=nx,l=nz)
      dimension ux(n,l),uz(n,l)
      dimension th(n,l),qv(n,l),qc(n,l),qr(n,l),fth(n,l),
     1         fqv(n,l),fqc(n,l),fqr(n,l)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common /prof_a/ tau(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

      do k=1,nz
      coe1=.5*dt*tau(k)
      coe2=1.+coe1
      do i=1,nx
      ux(i,k)=(ux(i,k)+coe1*ux_e(k))/coe2
      uz(i,k)=(uz(i,k)             )/coe2
      th(i,k)=(th(i,k)+coe1*th_e(k))/coe2
      qv(i,k)=(qv(i,k)+coe1*qv_e(k))/coe2
      qc(i,k)=(qc(i,k)             )/coe2
      qr(i,k)=(qr(i,k)             )/coe2

      fth(i,k)=fth(i,k) - tau(k)*(th(i,k)-th_e(k))
      fqv(i,k)=fqv(i,k) - tau(k)*(qv(i,k)-qv_e(k))
      fqc(i,k)=fqc(i,k) - tau(k)*(qc(i,k)-     0.)
      fqr(i,k)=fqr(i,k) - tau(k)*(qr(i,k)-     0.)
      enddo
      enddo
      return
      end

      subroutine integz(a,b,n1,n2)
      dimension a(n1,n2),b(n1,n2)

cc z smoothing:
      do k=2,n2-1
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do i=1,n1
      b(i,1 )=0.5*(a(i, 2)+a(i,   1))
      b(i,n2)=0.5*(a(i,n2)+a(i,n2-1))
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
      end

      subroutine integxz(a,b,n1,n2)
      dimension a(n1,n2),b(n1,n2)

cc z smoothing:
      do k=2,n2-1
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do i=1,n1
      b(i,1 )=0.5*(a(i, 2)+a(i,   1))
      b(i,n2)=0.5*(a(i,n2)+a(i,n2-1))
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

cc x smoothing:
      do i=1,n1
      im=i-1
      if(im.eq.0) im=n1-1
      ip=i+1
      if(ip.eq.n1+1) ip=2
      do k=1,n2
      b(i,k)=0.25*(a(im,k)+2.*a(i,k)+a(ip,k))
      enddo
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
      end

      subroutine prof_init
      parameter(nx=4001,nz=51)
      dimension xx(nx),zz(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)

      common /const/ gg,cp,rg,rv
      common/latent/hlatv,hlats
      common/reference/ tt0,ee0

      common /surf_ls/ thsrf(nx),qvsrf(nx)

      parameter(npin=23)
      dimension press(npin),temp(npin),zin(npin),vap(npin)
      common// tme1(nz),the1(nz),qve1(nz),ue1(nz)

       data press  / 
     1  1008.00, 991.25, 945.50, 893.79, 836.06, 772.82, 705.22,
     1   635.05, 564.48, 495.73, 430.71, 370.78, 316.72, 268.82,
     1   226.98, 190.82, 159.87, 133.55, 111.29,  92.56,  52.31,
     1    22.08,   9.32/
       data temp  / 
     1    25.26,  24.13,  21.04,  18.66,  16.50,  13.41,   9.06,
     1     3.73,  -1.51,  -6.97, -14.09, -22.44, -30.57, -39.60,
     1   -48.69, -57.40, -65.21, -72.58, -76.71, -74.98, -74.98,
     1   -74.98, -74.98/
       data vap  / 
     1   0.178E+02, 0.172E+02, 0.156E+02, 0.134E+02, 0.111E+02,
     1   0.888E+01, 0.631E+01, 0.487E+01, 0.396E+01, 0.200E+01,
     1   0.984E+00, 0.806E+00, 0.370E+00, 0.135E+00, 0.599E-01,
     1   0.258E-01, 0.123E-01, 0.582E-02, 0.367E-02, 0.589E-02,
     1   0.104E-02, 0.247E-02, 0.585E-02/

convert from temperature (deg C or K) into potential temperature
       do k=1,npin
c        temp(k)=temp(k)*(1.e3/press(k))**(rg/cp)
         temp(k)=(temp(k)+273.16)*(1.e3/press(k))**(rg/cp)
       enddo

        zin(1)=0.
        do k=2,npin
          km=k-1
          tempk =temp(k )*(1.e3/press(k ))**(-rg/cp)
     .                          * (1.+.6e-3*vap(k ))
          tempkm=temp(km)*(1.e3/press(km))**(-rg/cp)
     .                          * (1.+.6e-3*vap(km))
          delt=tempk-tempkm
          if (delt.gt.1.e-4) then
            tavi=alog(tempk/tempkm)/delt
          else
            tavi=1./tempk
          endif
          deltz=-rg/(tavi*gg) * alog(press(k)/press(km))
          zin(k)=zin(km)+deltz
        end do
cc define profiles:
      l=nz
cc surface data:
      iisn=1
      the1(1)=temp(iisn)
      tme1(1)=the1(1) * (1000./press(iisn))**(-rg/cp)
      qve1(1)=vap(iisn)*1.e-3
      ue1(1)=0.
c      print*,'DETERMINED SURFACE DATA'
cc higher levels - interpolate:
c      print*,'INTERPOLATION TO HIGHER LEVELS'
      do k=2,l
       zzz=(k-1)*dz
c        print*,'k=',k
        do kk=2,npin
          iisn=kk-1
          if(zin(kk).ge.zzz) go to 665
        enddo
        print*,'INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.'
        stop 'SOUNDING'
 665    continue
c        print*,'iisn=',iisn
        coe2=(zzz-zin(iisn))/(zin(iisn+1)-zin(iisn))
        the1(k)=coe2*temp(iisn+1) + (1.-coe2)*temp(iisn)
        qve1(k)=(coe2*vap(iisn+1) + (1.-coe2)*vap(iisn))*1.e-3
        presnl=coe2*press(iisn+1) + (1.-coe2)*press(iisn)
        tme1(k)=the1(k) * (1000./presnl)**(-rg/cp)
        ue1(k)=0.
      end do
c      print*,'ENVIRONMENTAL PROFILES'
c      do k=1,l
c      print 200,(k-1)*dz/1.e3,the1(k),tme1(k),qve1(k)*1.e3,ue1(k)
c 200    format(1x,'z,the,tme,qve,ue:',3f10.3,e12.3,f10.3)
c      enddo

cc no topography:
       do k=1,nz
        th_e(k)=the1(k)
        tm_e(k)=tme1(k)
        qv_e(k)=qve1(k)
        ux_e(k)= ue1(k)
       enddo
compute th00,tt00,pr00,rh00 and average stability for base state
      th00=the1(1)
      tt00=tme1(1)
      tvirt=tme1(1)*(1.+.6*qve1(1))
      rh00=press(1)*100./(rg*tvirt)
      pr00=press(1)*100.
      sum=0.
        do k=2,l-1
          sum = sum + (the1(k+1)-the1(k-1))/the1(k)
        enddo
      st=sum/(float(l-2)*2.*dz)

c      print*,'th00,tt00,pr00,rh00,st: ',th00,tt00,pr00,rh00,st
c      if(st.lt.1.e-10) stop 'st in init'

compute reference state vertical profiles
      cap=rg/cp
      capi=1./cap
      cs=gg/(cp*tt00*st)
        do k=1,nz
         zzz=(k-1)*dz
          th0(k)=th00*exp(st*zzz)
          exs=exp(-st*zzz)
          rho0(k)=rh00*exs*(1.-cs*(1.-exs))**(capi-1.)
        enddo
c      do k=1,l
c       print 207,(k-1)*dz/1.e3,th0(k),rho0(k)
c 207    format(1x,'z,th0,rho0',3f12.3)
c      enddo

cc surface data:
      a=rg/rv
      b=hlatv/(rv*tt0)
      c=hlatv/cp
      d=hlatv/rv
      e=-cp/rg
       pi=4.*atan(1.)
       do i=1,nx
c       thsrf(i)=th_e(1) + 1.  
       thsrf(i)=th_e(1) + 1. + 
     *     2.0*cos(2.*pi*float(i-1)/float(nx-1) - pi)
      thetme=th_e(1)/tm_e(1)
      pre=1.e5*thetme**e
      tt=thsrf(i)/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      qvsrf(i)=a*esw/(pre-esw)
      enddo

       do i=1,nx
       print*,'i,thsr,qvsr: ',thsrf(i),qvsrf(i)
       enddo

      return
      end

      subroutine surfflux(theta,qv,ux,fth,fqv,th_e,nx,nz)
      dimension theta(nx,nz),qv(nx,nz),ux(nx,nz),fth(nx,nz),fqv(nx,nz)      
      dimension th_e(nz),ss1(200),ss2(200)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /const/ gg,cp,rg,rv
cc ocean surface data:
      parameter(nx1=4001)
      common /surf_ls/ thsrf(nx1),qvsrf(nx1)

      dist(z,H)=(H-amin1(z,H))/H

      if(nz.gt.200) stop 'surfflux dim'

c calculate surface fluxes based on local conditions:
       drag=1.3e-3
       coeth=drag
       coeqv=drag
       hpbl=600.   ! prescribed pbl depth (m)

       do i=1,nx

         wind=ux(i,1)
           alf=drag*((thsrf(i)-theta(i,1))/th_e(1) 
     .          + .61*(qvsrf(i)-qv(i,1)))
            w_star=sqrt(abs(gg*hpbl*alf))
             wind=sqrt(wind*wind+w_star*w_star)
           wind=amax1(wind,2.)
        fl_sens=drag*wind*(thsrf(i)-theta(i,1))  ! sensible flux
        fl_lat =drag*wind*(qvsrf(i)-   qv(i,1))  ! latent flux

cc distribute in the vertical:
      do k=1,nz
cc note shift of flux positions:
      zz=(k-1)*dz+.5*dz
      fun=dist(zz,300.) 
      ss1(k)=fl_sens*fun
      ss2(k)=fl_lat *fun
      enddo

cc first level above the ground:
      fth(i,1)=fth(i,1)-(ss1(1)-fl_sens)/(0.5*dz) * 2.
      fqv(i,1)=fqv(i,1)-(ss2(1)- fl_lat)/(0.5*dz) * 2.
cc higher levels:
      do k=2,nz
      fth(i,k)=fth(i,k)-(ss1(k)-ss1(k-1))/dz *   2.
      fqv(i,k)=fqv(i,k)-(ss2(k)-ss2(k-1))/dz *   2.
      enddo

ccc print values for test:
c        print*,' @@@ box, s fluxes: ',i,fl_sens*1.e3,fl_lat*2.5e6
       enddo

      return
      end

      function cvmgm(a,ab,ac)
       if(ac.lt.0) then
        cvmgm=a
       else
        cvmgm=ab
       endif
       return
       end

      subroutine tape_wr(ux,uz,uxp,uzp,theta,qv,qc,qr,
     .                   fx,fz,ft,fqv,fqc,fqr,p,nx,nz)
      dimension theta(nx,nz),qv(nx,nz),qc(nx,nz),qr(nx,nz)
      dimension ux(nx,nz),uz(nx,nz)     ! current time level
      dimension uxp(nx,nz),uzp(nx,nz)   ! previous time level
      dimension ft(nx,nz),fx(nx,nz),fz(nx,nz),
     *          fqv(nx,nz),fqc(nx,nz),fqr(nx,nz),p(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

      write(21) time,ux,uz,uxp,uzp,theta,qv,qc,qr
      write(21) fx,fz,ft,fqv,fqc,fqr,p

      return
      end 

      subroutine tape_rd(ifl,ux,uz,uxp,uzp,theta,qv,qc,qr,
     .                   fx,fz,ft,fqv,fqc,fqr,p,nx,nz)
      dimension theta(nx,nz),qv(nx,nz),qc(nx,nz),qr(nx,nz)
      dimension ux(nx,nz),uz(nx,nz)     ! current time level
      dimension uxp(nx,nz),uzp(nx,nz)   ! previous time level
      dimension ft(nx,nz),fx(nx,nz),fz(nx,nz),
     *          fqv(nx,nz),fqc(nx,nz),fqr(nx,nz),p(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

       do if=1,ifl
      read(31,end=10) time,ux,uz,uxp,uzp,theta,qv,qc,qr
      read(31) fx,fz,ft,fqv,fqc,fqr,p
       print*,' ****** read data for time (min): ',time
       enddo
       return
 10    continue
      print*,'!! restarted with smaller number of files...'
      return
      end 
cc    DF: QMCL Subroutines 

      subroutine generate_qmcl_objects_qv(training_data
     &, L, number_training_pts, nx, nz, 
     & phi, ntime_training, D, N_real, N_qm, S_matrix, U,
     & trans_U, pca_coeff_matx_t, pca_dims, training_pca)
c     [Tentatively completed]
c     Note: still need 
c     Inputs: all these static variables and the dataset 
c     Output: U operator, mult operator array S_matrix

c     Import eigenfunctions and training dataset 

cQ: Can a subroutine generate a new variable, or can it only alter 
c currently existing variables
      real S_matrix(L,L,4, nz)
      integer pca_dims
c     DF: change "4" to "res_var" when possible
      real S(L,L), U(L,L), trans_U(L,L)
      integer nx, nz, L, D, number_training_pts
      integer k, e, j, i, m, N_qm
      real N_real, num_tr_pts_real
      real phi(number_training_pts, L)
      real phi_1(number_training_pts-1, L)
c      real tr_phi_1(L, number_training_pts-1)
      real phi_2(number_training_pts-1, L)
      real training_data(nx*nz, ntime_training, 7)
      real training_data_trunc(nx, nz*4, ntime_training)
      real pca_data_interm(nz*4, ntime_training)
      real pca_coeff_matx_t(pca_dims, nz*4)
c     Note: fourth argument in training_data dims corresponds to how many 
c     variables we are deciding to use in updating the state 
c     Each extra dim corresponds to another variable
      real Y(number_training_pts,1)
      real training_pca(nx,pca_dims, N_qm)

c     Creates Koopman operator 

      print*, "Began subroutine"

      num_tr_pts_real = REAL(number_training_pts)

cc 7.17.24 DF - Why was this even here?
      do m=1,ntime_training
      do i=1,nx
            do j=1,nz
                  do k=1,4
      training_data_trunc(i,j,k) = training_data(i+nz*(j-1),m,k)
                  enddo 
            enddo
      enddo
      enddo

c      training_data_trunc = reshape(training_data(:,:,1:4), 
c     &shape(training_data_trunc))
      
      do j=1,nx
      pca_data_interm = training_data_trunc(j,:,:)
      training_pca(j,:,:) = matmul(pca_coeff_matx_t, pca_data_interm)
      enddo

      do i=1,L 
            do j=1,L 
                  S(i,j) = 1
            enddo
      enddo

      print*, "phi vals" 
      print*, phi(1,1) 
      print*, phi(2,2)
      print*, phi(10,10)
      do i=1,(number_training_pts-1)
            do j=1,L
            phi_1(i,j) = phi(i,j)
            phi_2(i,j) = phi((i+1), j)
            end do
      end do
      print*,"phi_1 and phi_2 vals" 
      print*, phi_1(10,10) 
      print*, phi_2(1,1) 
      print*, phi_2(10,10)
c      tr_phi_1 = transpose(phi_1)
c      U = (1/num_tr_pts_real)*matmul(tr_phi_1, phi_2)
      U = (1/num_tr_pts_real)*matmul(transpose(phi_1), phi_2)
      do i=1,50
            do j=1,50
                  trans_U(i,j) = U(j,i)
            enddo
      enddo
      print*, "Generated U"
      print*, "U vals" 
      print*, U(1,1) 
      print*, U(10,10)
cc    Generates Y vectors 
      do k=1,4
c     k now iterates over variables, NOT over nx since 
c     we're using one box per height z
      do e=1,nz
            do i=1,number_training_pts
c     Need to figure out how this works with the training data
c     And change syntax 
c            Y(i,1) = training_data(nx*nz+nz*(k-1) + e,i) 
            Y(i,1) = training_data(1 + nx*(e-1),i,k) 

c     Note: make sure fourth argument corresponds to correct variable 
            end do
c      print*, "Got to generate S"
c      print*, "Y values"
c      print*, Y(1,1)
c      print*, Y(50,1)
      CALL generate_S(S, Y, phi, N_qm, L, N_real)
c      print*, "Finished generate S"
      do j=1,L
            do m = 1,L 
                  S_matrix(j,m,k,e) = S(j,m)
            enddo
      enddo
      enddo
      enddo


      end subroutine generate_qmcl_objects_qv
      
      subroutine thermo_qmcl_qv(qv, qc, qr,  
     & theta, L, ux, uz, p,
     & S_array_qv, training_data, phi, U, N_qm
     &, measurement_epsilon, N_real, D, nx, nz,
     & rho, T_bar, sigma_T, update_index, trans_U,
     & trans_phi, training_data_reshape, pca_terms,
     & pca_data, pca_coeff_mat_tpose)
c     phi is eigenfunctions matrix 
      real rho(L,nx)
      real U(L,L), rho_interm(L), trans_U(L,L)
      real vec_interm(L)
      real S_array_qv(L,L,4,nz)
      real S_qv(L,L),S_qr(L,L),S_qc(L,L),S_th(L,L)
      real training_data(nx*nz, ntime_training, 7)
      real phi(ntime_training, L)
      real trans_phi(L, ntime_training)
      real measurement_epsilon
      integer i,j,k,m
      integer kk,ll, mm, indexx
      integer nx, nz
      integer update_index
      real trace_val 
      real qv(nx,nz), qc(nx,nz), qr(nx,nz)
      real theta(nx,nz)
      real ux(nx,nz), uz(nx,nz), p(nx,nz)
      real covariate(nx*nz, 1)
      real N_real
      real T_bar(nz,7)
      real sigma_T(7)
      real sigma_unres(4,1)
      real T_bar_unres(4,1)
      real interm_prime(4,4)
      real closure_terms(4,1)
c      real covariate_array(3,3,3)
c     ^length and width of local nbhd, plus # of resolved vars (3) 
      integer D 
      integer pca_terms
      real covariate_array(D,nz)
      real qv_interm(L,L), qr_interm(L,L), qc_interm(L,L)
      real th_interm(L,L)
      real qv_interm_prime, qr_interm_prime
      real qc_interm_prime, th_interm_prime
      real diag1(L,1),diag2(L,1),diag3(L,1),diag4(L,1)
      real training_data_reshape(nx,D*nz,N_qm)
      real pca_data(nx, pca_terms, N_qm)
      real pca_coeff_mat_tpose(pca_terms, 4*nz)
      do i=1,4
      sigma_unres(i,1) = sigma_T(i)
      T_bar_unres(i,1) = T_bar(i,1)
      enddo 
      interm_prime = 0.0
      qv_interm = 0.0 
      qr_interm = 0.0
      qc_interm = 0.0
      th_interm = 0.0

      do j=1,nx
      print*, j
      rho_interm = rho(:, j)

      if(update_index.EQ.2) then

      do kk=1,nz
            covariate_array(1,kk)=(qv(j,kk) - T_bar(kk,1))/sigma_T(1)
            covariate_array(2,kk)=(qr(j,kk) - T_bar(kk,2))/sigma_T(2)
            covariate_array(3,kk)=(qc(j,kk) - T_bar(kk,3))/sigma_T(3)
            covariate_array(4,kk)=(theta(j,kk) - T_bar(kk,4))/sigma_T(4)
            covariate_array(5,kk)=(ux(j,kk) - T_bar(kk,5))/sigma_T(5)
            covariate_array(6,kk)=(uz(j,kk) - T_bar(kk,6))/sigma_T(6)
      enddo

      CALL update_rho_pca(rho_interm,phi,covariate_array,
     &  training_data, measurement_epsilon,L,
     &  N_qm, N_real, D, nx, nz, j, trans_phi, training_data_reshape,
     & pca_terms, pca_data(j,:,:), pca_coeff_mat_tpose)
      endif


      do i=1,nz

      S_qv = S_array_qv(1:L, 1:L, 1,i)
      S_qr = S_array_qv(1:L, 1:L, 2,i)
      S_qc = S_array_qv(1:L, 1:L, 3,i)
      S_th = S_array_qv(1:L, 1:L, 4,i)


      
      call dsymm('R','L',L,L,1.0d0, rho_interm,
     &L,S_qv,L,0.0d0, qv_interm,L) 
      call dsymm('R','L',L,L,1.0d0, rho_interm,
     &L,S_qr,L,0.0d0, qr_interm,L) 
      call dsymm('R','L',L,L,1.0d0, rho_interm,
     &L,S_qc,L,0.0d0, qc_interm,L) 
      call dsymm('R','L',L,L,1.0d0, rho_interm,
     &L,S_th,L,0.0d0, th_interm,L) 

      do k=1,L 
            diag1(k,1) = qv_interm(k,k)
            diag2(k,1) = qr_interm(k,k)
            diag3(k,1) = qc_interm(k,k)
            diag4(k,1) = th_interm(k,k)
      enddo

      interm_prime(1,1) = SUM(diag1) 
      interm_prime(2,2) = SUM(diag2) 
      interm_prime(3,3) = SUM(diag3) 
      interm_prime(4,4) = SUM(diag4) 
c     ^diagonal so vec prod can be done in BLAS

      call dsymv('L', 4, 1.0d0, interm_prime,
     &4, sigma_unres, 1,1.0d0, T_bar_unres,1)
c    T_bar_unres gets overwritten by the solutions

      qv(j,i) = T_bar_unres(1,1)
      qr(j,i) = T_bar_unres(2,1)
      qc(j,i) = T_bar_unres(3,1)
      theta(j,i) = T_bar_unres(4,1)
      enddo
      
c      CALL evolve_rho_fast(rho_interm, U, L, trans_U)
c     Update rho vectorized
      call dgemv('N', L,L,1.0,U, L,rho_interm,1,0.0,vec_interm,1)
      trace_value= 0.0
      do kk=1,L
            trace_value = trace_value + 
     &(vec_interm(kk))**2
      enddo 
      rho_interm = (1/sqrt(trace_value))*vec_interm

      do kk=1,L 
            rho(kk,j) = rho_interm(kk)
      enddo


      enddo

c      print*, "Covariate updated"

      end subroutine thermo_qmcl_qv


       
      subroutine generate_S(S, Y, phi, N_qm, L_qm, N_real)
       
      integer N_qm 
      real N_real
      integer L_qm
      real Y(N_qm,1)
      real phi(N_qm, L_qm)
      real S(L_qm,L_qm)
      real A_int(L_qm,L_qm)
      integer k, ii, jj, kk
      real intermediate_sum
                        
            
      intermediate_sum = 0.0

      do ii=1,L_qm
      do jj=1,L_qm
            do kk = 1,N_qm
            intermediate_sum = intermediate_sum 
     .+ (phi(kk,ii)*Y(kk,1))*phi(kk,jj)
            end do
            A_int(ii,jj) = (1.0/N_real)*intermediate_sum
            intermediate_sum = 0.0
      end do

      end do
      S = (A_int + transpose(A_int))/2.0
      return
      end subroutine generate_S
       
      subroutine trace(tr_val, A, size)
      
            implicit none
            integer size
            real A(size, size)
            real x, tr_val
            integer i 
      
            x=0.0
      do i = 1, size
      x = x + A(i,i)
      end do
      tr_val = x
      
      end subroutine trace
       
      subroutine evolve_rho_fast(rho_temp,U, spectral_res, trans_U) 
            integer spectral_res
            real U(spectral_res, spectral_res)
            real trans_U(spectral_res, spectral_res)
            real standin(spectral_res, spectral_res)
            real rho_temp(spectral_res,1)
            real intermediate(spectral_res, spectral_res)
            real tr_val 
            integer m,np,k
            integer lpa, lpb,lpc
            real alpha, beta

            m=50
            np=50
            k=50
            alpha = 1.0
            beta = 1.0
c           Should be 0^
            lpa = 50
            lpb=50
            lpc=50

           standin = 0.0
           intermediate = 0.0
c            call dsymm('R','L',spectral_res,spectral_res,
c     &1.0d0, trans_U,spectral_res,rho,spectral_res,0.0,
c     &standin,spectral_res) 
            call dgemm('N','N',m,
     &np,k, alpha,trans_U,
     &lpa,rho_temp,lpb,
     &beta, standin, lpc)

            call dgemm('N','N',m,
     &np,k, alpha,U,
     &lpa,standin,lpb,
     &beta, intermediate, lpc)
c            call dgemm('N','N',spectral_res,spectral_res,
c     &spectral_res, alpha,standin,spectral_res,U,spectral_res,
c     &beta, intermediate, spectral_res)
      tr_val = 0.0
 
            CALL trace(tr_val, intermediate, spectral_res)
            print*, "Trace Value" 
            print*, tr_val
            rho_temp = (1/tr_val)*intermediate

      end subroutine evolve_rho_fast
                     
      subroutine update_rho_pca(rho, phi, new_x_array,training_data,
     & epsilon, spectral_res, training_len, N_real, D, nx, nz, x_cord
     &, tr_phi, training_data_reshape, pca_terms, pca_data, C_T)
cc    x_cord and y_cord represent the spatial point at the center of the 
cc    neighborhood we are using to update with 
cc    For column updates, removed the "z_coord" argument
c     C_T is the transpose of the pca coefficients matrix
cc    D = number of resolved variable dimensions 
      integer spectral_res, training_len
      real N_real
      real rho(spectral_res)
      real new_rho(spectral_res)
      real G(spectral_res, spectral_res)
      real interm_G(spectral_res, spectral_res)
      real phi(training_len, spectral_res)
      real tr_phi(spectral_res, training_len)
      real phi_prime(training_len, spectral_res)
      real F(training_len)

ccccc          
      integer D
      integer nx, nz
      integer x_cord, z_cord
      real training_data(nx*nz, training_len, 7)
      real x_matrix(D*nz, training_len)
      real x_matrix_new(D*nz, training_len)
c      real x_matrix_newnew(D*nz, training_len)
      real new_x_array(D,nz)
      real new_x_spread(D*nz,training_len)
c     ^will be different for local nbhd kernel 
      real new_x(D*nz)
cc    Edit above dimensions when changing what variables/data used to train 
      real epsilon, tr_value
      integer i,j,k,l
      real dummy 
      integer m,np,kk,lpa,lpb,lpc
      real alpha, beta, alpha2
      real intermediate(spectral_res, spectral_res)
      real column_sums(training_len)

      real F_spread(training_len, spectral_res)
      real training_data_reshape(nx,D*nz,training_len)
      real training_interm(nz,training_len,D)
      
      real F_reshape(training_len, spectral_res)
      real diag1(spectral_res)

      integer pca_terms
      real pca_data(pca_terms, training_len)
      real C_T(pca_terms, D*nz)
      real new_x_pca(pca_terms)

      m=50
      np=50
      kk=training_len
      alpha = 1.0/N_real
      beta = 0.0
      lpa = 50
      lpb=training_len
      lpc=50

      alpha2 = 1.0


      
c      do i=1,D
c            do j=1,nz
c                  new_x((D*(j-1)) + i,1) = 
c     &new_x_array(i,j)
c            enddo
c      enddo

      new_x = reshape(new_x_array, shape(new_x))
c      new_x_pca = matmul(C_T, new_x(:,1))

      call dgemv('N',pca_terms,D*nz,alpha_2,C_T, pca_terms
     &,new_x,1,0.0,new_x_pca,1)
c      x_matrix = training_data_reshape(x_cord,:,:)
      new_x_spread = spread(new_x_pca(:), 2,training_len)
      x_matrix_new = (pca_data-new_x_spread)**2
c      x_matrix_newnew = x_matrix_new**2
c      column_sums = SUM(x_matrix_newnew, dim=1)
      column_sums = SUM(x_matrix_new, dim=1)

    
      F = EXP(column_sums*(-1.0/100.0))
      F_reshape = spread(F, 2, spectral_res)

      phi_prime = phi*F_reshape

      call dgemm('N','N',m,
     &np,kk, alpha,tr_phi,
     &lpa,phi_prime,lpb,
     &beta, G, lpc)

      call dgemv('N',spectral_res,spectral_res,1.0,G, spectral_res
     &,rho,1,0.0,new_rho,1)

      tr_value = 0.0

      do l=1,spectral_res
            diag1(k) = new_rho(l)**2
      enddo

      tr_value = SUM(diag1) 

c      do i=1,spectral_res
c      tr_value = tr_value + (new_rho(i))**2
c      enddo
      rho = new_rho/sqrt(tr_value)



           
c      G = (1/N_real)*matmul(tr_phi,phi_prime)
c     DF - try comment section
c      interm_G = matmul(rho, G)
c      new_rho = matmul(G, matmul(rho, G))
c      new_rho = matmul(G,matmul(rho,G))

      end subroutine update_rho_pca
       
       
             
ccccccccccc End of QMCL subroutines 