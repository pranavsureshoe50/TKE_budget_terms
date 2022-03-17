Program Main
  implicit none
  
  INCLUDE 'header'


! ------------------Input Paramters - Pranav ----------------------------------------------------------- 
  integer ( kind = 4 ), parameter :: nx = 402
  integer ( kind = 4 ), parameter :: ny = 1026
  integer ( kind = 4 ), parameter :: nz = 1538
  character(len=9), parameter :: meandir = 'meanField'
  character(len=6), parameter :: dir = '5pData'
 
  character(len=19), parameter :: edir = 'turbdis'
  integer ( kind = 4 ), parameter :: nx_start = 1
  integer ( kind = 4 ), parameter :: nx_end = 402

  integer ( kind = 4 ), parameter :: ny_start = 1
  integer ( kind = 4 ), parameter :: ny_end = 1026
  
  integer ( kind = 4 ), parameter :: nz_start = 1 
  integer ( kind = 4 ), parameter :: nz_end = 1538

  integer ( kind = 4 ), parameter :: nxx = 150
  integer ( kind = 4 ), parameter :: nyy = 426
  integer ( kind = 4 ), parameter :: nystart = 300

  integer ( kind = 4 ), parameter :: kmin5p = 260
  integer ( kind = 4 ), parameter :: kmax5p = 800
 

!-------------------Input parameters -------------------------------------------------------------------

  
  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx)
  real    ( kind = 8 ) ::  x5p(nxx), y5p(nyy) , dum1,dum2

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imax,jmax,kmax,tmp,kstart,kend,s1
  integer ( kind = 4 ) :: tag,nkrr,r
  real    ( kind = 8 ) :: time1,TOTALT,T_PRE,TSTEP1,DTM1,T_INI,nu,dx,dy,dz
  real    ( kind = 8 ) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  character(len=128)   :: buff,filename, ICFile
  integer (kind = 4)   :: nx_write, ny_write, nz_write
  
  logical file_exists, DEXIST
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: u5p,v5p,w5p,dens5p,p5p
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: u1,v1,w1,p1,dens1,baro_ver
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: uc,vc,wc,pc,dc
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: umean,vmean,wmean,pmean,densmean
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: diffu_st,diffu_sw,diffu_vr
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: stretch_ver, diss, int_diss
  REAL    ( kind = 8 ) , ALLOCATABLE, dimension(:) :: z5p

  real    ( kind = 8 ) :: theta(ny),z(nz),dtheta
  logical              :: save_vel_box,save_full_box, save_vel_full, save_full

  integer :: nk,kxx0,kxx1,kxx2,kxx3,KXX4

  REAL, DIMENSION(nx_start:nx_end)  :: x_w
  REAL, DIMENSION(ny_start:ny_end)  :: y_w
  REAL, DIMENSION(nz_start:nz_end)  :: z_w
!----------------------------------------------------------------------------------
  nx_write = (nx_end - nx_start) + 1
  ny_write = (ny_end - ny_start) + 1
  nz_write = (nz_end - nz_start) + 1
 
  nu = 0.001667
  nkrr = kmax5p - kmin5p +1 + 5*FLOOR((NZ-2-KMAX5P)/20.0)

  ALLOCATE( z5p(nkrr) ) 
  call read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nz,ru,rp,tag)
  k=1
  nk=1
  do while (k.le.nkrr) 
    
     if(k.le. (kmax5p-kmin5p+1)) then
      z5p(k) = zcg(kmin5p-1+k)
      k=k+1
     else
      
! 5p location             
      kxx0 = KMAX5P+NK*20-2
      KXX1 = KMAX5P+NK*20-1
      KXX2 = KMAX5P+NK*20
      KXX3 = KMAX5P+NK*20+1
      KXX4 = KMAX5P+NK*20+2

      z5p(k) = zcg(kxx0)
      if(k+1 .le. nkrr) z5p(k+1) = zcg(kxx1)
      if(k+2 .le. nkrr) z5p(k+2) = zcg(kxx2)
      if(k+3 .le. nkrr) z5p(k+3) = zcg(kxx3)
      if(k+4 .le. nkrr) z5p(k+4) = zcg(kxx4)

      k=k+5
      nk=nk+1
     endif      

  enddo
  do k= 1, nkrr
!     write(6,*) "z5p(",k,") = ", z5p(k)
  enddo

  do i=1,nxx
     x5p(i) = xc(i)
!     write(6,*) "x5p(",i,") = ", xc(i)
  enddo
  do j=1,nyy
     y5p(j) = yc(j+nystart)
!     write(6,*) 'y5p(',j,')=' , yc(j+nystart)
  enddo



! Allocate 5p data
  allocate(u5p(nxx,nyy,nkrr),v5p(nxx,nyy,nkrr),w5p(nxx,nyy,nkrr),dens5p(nxx,nyy,nkrr),p5p(nxx,nyy,nkrr))
  allocate(umean(nxx,nyy,nkrr),vmean(nxx,nyy,nkrr),wmean(nxx,nyy,nkrr),densmean(nxx,nyy,nkrr),pmean(nxx,nyy,nkrr))
  allocate(uc(nxx,nyy,nkrr),vc(nxx,nyy,nkrr),wc(nxx,nyy,nkrr),dc(nxx,nyy,nkrr),pc(nxx,nyy,nkrr))
  allocate(u1(nxx,nyy,nkrr),v1(nxx,nyy,nkrr),w1(nxx,nyy,nkrr),dens1(nxx,nyy,nkrr),p1(nxx,nyy,nkrr))
  allocate(diss(nxx,nyy,nkrr),int_diss(nxx,nyy,nkrr))


  umean=0.0
  vmean=0.0
  wmean=0.0
  pmean=0.0
  densmean=0.0
  uc = 0.0
  vc =0.0
  wc =0.0
  pc = 0.0
  dc =0.0 
  int_diss = 0.0
  T_PRE  = 0.0 
  TOTALT=0.0d0


  nx_write = nxx
  ny_write = nyy
  nz_write = nkrr 
 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        write(*,*) 'READING MEAN FIELD ..............'

        filename = trim(meandir) // "/UMEAN.5P"
        call read_post(filename,nx_write,ny_write,nz_write,umean,time1,dum1)

        filename = trim(meandir)// "/VMEAN.5P"
        call read_post(filename,nx_write,ny_write,nz_write,vmean,time1,dum1)


        filename = trim(meandir) // "/WMEAN.5P"
        call read_post(filename,nx_write,ny_write,nz_write,wmean,time1,dum1)

        filename = trim(meandir) // "/PMEAN.5P"
        call read_post(filename,nx_write,ny_write,nz_write,pmean,time1,dum1)


        filename = trim(meandir) // "/DENSMEAN.5P"
        call read_post(filename,nx_write,ny_write,nz_write,densmean,time1,dum1)

! ----------------------------------------------------------------------
        write(6,*) 'Done Read mean field ----------'

! ---------------------------------------------------------------------
! Looping over 5p files
        do r=res1,res2,pstride

! Check if files exist
! If not exist, variables of the next timestep will be assumed 
!to occupy 1000 steps instead of 500 (in case skipping 500)
!---------------------------------------------------------------------!           
           file_exists = .true.
           write(ICfile,'(a,i8.8,a)')trim(dir)//"/up_",r,".5p"
           INQUIRE(FILE=ICfile, EXIST=file_exists)
           IF(file_exists .eqv. .false.) THEN
              CYCLE
           ENDIF

           file_exists = .true.
           write(ICfile,'(a,i8.8,a)')trim(dir)//"/vp_",r,".5p"
           INQUIRE(FILE=ICfile, EXIST=file_exists)
           IF(file_exists .eqv. .false.) THEN
              CYCLE
           ENDIF

           file_exists = .true.
           write(ICfile,'(a,i8.8,a)')trim(dir)//"/wp_",r,".5p"
           INQUIRE(FILE=ICfile, EXIST=file_exists)
           IF(file_exists .eqv. .false.) THEN
              CYCLE
           ENDIF

           file_exists = .true.
           write(ICfile,'(a,i8.8,a)')trim(dir)//"/pp_",r,".5p"
           INQUIRE(FILE=ICfile, EXIST=file_exists)
           IF(file_exists .eqv. .false.) THEN
              CYCLE
           ENDIF

           file_exists = .true.
           write(ICfile,'(a,i8.8,a)')trim(dir)//"/densp_",r,".5p"
           INQUIRE(FILE=ICfile, EXIST=file_exists)
           IF(file_exists .eqv. .false.) THEN
              CYCLE
           ENDIF
!---------------------------------------------------------------------!           

  u5p=0.0
  v5p=0.0
  w5p=0.0
  p5p=0.0
  dens5p=0.0
  uc = 0.0
  vc = 0.0
  wc = 0.0
  diss = 0.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        write(*,*) 'READING 5P FILES ..............'

        write(ICfile,'(a,i8.8,a)')trim(dir)//"/up_",r,".5p"
        write(6,*) 'Reading u  file ...', ICfile
        call read_5p(ICfile,nxx,nyy,nkrr,u5p,time1,dtm1)

        write(ICfile,'(a,i8.8,a)')trim(dir)//"/vp_",r,".5p"
        write(6,*) 'Reading v  file ...',ICfile
        call read_5p(ICfile,nxx,nyy,nkrr,v5p,time1,dtm1)


        write(ICfile,'(a,i8.8,a)')trim(dir)//"/wp_",r,".5p"
        write(6,*) 'Reading w  file ...',ICfile
        call read_5p(ICFile,nxx,nyy,nkrr,w5p,time1,dtm1)

        write(ICfile,'(a,i8.8,a)')trim(dir)//"/pp_",r,".5p"
        write(6,*) 'Reading p  file ...',ICfile
        call read_5p(ICFile,nxx,nyy,nkrr,p5p,time1,dtm1)

        write(ICfile,'(a,i8.8,a)')trim(dir)//"/densp_",r,".5p"
        write(6,*) 'Reading dens file ...', ICfile
        call read_5p(ICFile,nxx,nyy,nkrr,dens5p,time1,dtm1)


        write(*,*) 'time,r :' ,  time1,r
 
        u1(:,:,:) = u5p(:,:,:)-umean(:,:,:)
        v1(:,:,:) = v5p(:,:,:)-vmean(:,:,:)
        w1(:,:,:) = w5p(:,:,:)-wmean(:,:,:)

        DO k=2,NKRR-1
        DO j=2,NYY-1
        DO i=2,NXX-1

           dz = z5p(k) - z5p(k-1)
           dy = y5p(j) - y5p(j-1)
           dx = x5p(i) - x5p(i-1)
   
          dudx = (u1(i,j,k) - u1(i-1,j,k))/dx
          dvdy = (v1(i,j,k) - v1(i,j-1,k))/dy
          dwdz = (w1(i,j,k) - w1(i,j,k-1))/dz

          dudy = (0.25*( u1(i,j,k)+u1(i,j+1,k)+u1(i-1,j+1,k)+u1(i-1,j,k) ) - &
             0.25*( u1(i,j,k)+u1(i,j-1,k)+u1(i-1,j-1,k)+u1(i-1,j,k) ) )/dy     
          dudz = (0.25*( u1(i,j,k)+u1(i,j,k+1)+u1(i-1,j,k+1)+u1(i-1,j,k) ) - &
             0.25*( u1(i,j,k)+u1(i,j,k-1)+u1(i-1,j,k-1)+u1(i-1,j,k) ) )/dz

          dvdx = (0.25*( v1(i,j,k)+v1(i,j-1,k)+v1(i+1,j-1,k)+v1(i+1,j,k) ) - &
              0.25*( v1(i,j,k)+v1(i,j-1,k)+v1(i-1,j-1,k)+v1(i-1,j,k) ) )/dx
          dvdz = (0.25*( v1(i,j,k)+v1(i,j-1,k)+v1(i,j-1,k+1)+v1(i,j,k+1) ) - &
              0.25*( v1(i,j,k)+v1(i,j-1,k)+v1(i,j-1,k-1)+v1(i,j,k-1) ) )/dz  

          dwdx = (0.25*( w1(i,j,k)+w1(i,j,k-1)+w1(i+1,j,k-1)+w1(i+1,j,k) ) - &
             0.25*( w1(i,j,k)+w1(i,j,k-1)+w1(i-1,j,k-1)+w1(i-1,j,k) ) )/dx
          dwdy = (0.25*( w1(i,j,k)+w1(i,j+1,k-1)+w1(i,j+1,k)+w1(i,j,k-1) ) - &
             0.25*( w1(i,j,k)+w1(i,j-1,k-1)+w1(i,j-1,k)+w1(i,j,k-1) ) )/dy   
          

          diss(i,j,k)   = nu * ( 2*dudx**2 + 2*dvdy**2 + 2*dwdz**2 &
                   + (dudy+dvdx)**2 + (dwdy+dvdz)**2 + (dudz+dwdx)**2 )

        ENDDO
        ENDDO
        ENDDO

        
        TSTEP1=TIME1-T_PRE
        if(r==res1) then
                TSTEP1=DTM1*pstride 
                T_INI = TIME1
        endif

                

        TOTALT=TOTALT+TSTEP1
        T_PRE=TIME1

        int_diss(:,:,:) =  int_diss(:,:,:) + diss(:,:,:)*TSTEP1
     
        WRITE(*,*)"TOTALT=",TOTALT,"TSTEP1=",TSTEP1

        if(TOTALT .GE. TTOTAL_LIM) exit

       enddo
 
! ---------------------------------------------------------------------
! End looping over 5p files 


    int_diss(:,:,:) =  int_diss(:,:,:)/TOTALT  

    CALL SYSTEM('mkdir ' // trim(edir) )
      
    WRITE(*,*)"MEAN DISS CALCULATED"

    filename  = trim(edir)  // '/TDISS.5P'
    call  write_5p(filename,nxx,nyy,nkrr,int_diss,T_INI,TIME1)

    filename  = trim(edir)  // '/Integ_TDISS.dat'
    call integxy(filename,int_diss,x5p,y5p,z5p,nxx,nyy,nkrr)

!    filename  = trim(edir)  // '/MKE.5P'    
!    call  write_5p(filename,nxx,nyy,nkrr,mke,T_INI,TIME1)
 
!   filename  =  trim(edir)  //  '/TKE.5P'    
!    call  write_5p(filename,nxx,nyy,nkrr,int_tke,T_INI,TIME1)
    
!    filename =    trim(edir)  //   '/tke.plt'
!    call write_plt_3d_1var(filename,int_tke,x5p,y5p,z5p,1,nxx,1,nyy,1,nkrr,nxx,nyy,nkrr)
    
!    filename =    trim(edir)  //  '/mke.plt'
!    call write_plt_3d_1var(filename,mke,x5p,y5p,z5p,1,nxx,1,nyy,1,nkrr,nxx,nyy,nkrr)




!   filename = 'meanField/vel_3D_full.plt'     
!   call  write_plt_3d_3vars(filename,umean,vmean,wmean,nxx,nyy,nkrr,x5p,y5p,z5p) 
!  do k=1,nkrr
!  do j=1,nyy
!      do i=1,nxx
! 
!       
!      enddo
!  enddo
!  enddo


  deallocate(u5p,v5p,w5p,p5p,dens5p)
  deallocate(umean,vmean,wmean,densmean,pmean)
  deallocate(uc,vc,wc,dc,pc)
  deallocate(u1,v1,w1,dens1,p1)
  deallocate(diss,int_diss)
  DEALLOCATE( z5p )
  stop
end Program Main



subroutine read_h(nz,ny,zch,ych,h)
  implicit none

  INTEGER nx,ny,nz,nzg,tag
  REAL ( kind = 8 ) :: zch(ny,nz),ych(ny,nz),h(ny,nz)


  INTEGER i,j,k

  OPEN(UNIT=1,FILE="./h_ridge.dat",STATUS='OLD',FORM='FORMATTED')

    do k=1,nz
        do j=1,ny
          read(1,*) zch(j,k), ych(j,k), h(j,k)
        enddo
    enddo
  close(1)

  write(6,*) "READ H DONE"

  return
end subroutine read_h

subroutine read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg,ru,rp,tag)
  implicit none

  INTEGER nx,ny,nz,nzg,tag
  REAL ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nzg)
  REAL ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nzg)
  REAL ( kind = 8 ) :: ru(nx),rp(nx)

  INTEGER i,j,k

  real rdelx,rdely,rdelz, dtheta
  real, allocatable, dimension(:) :: cug,cvg

  ALLOCATE(cug(nzg),cvg(nzg))

  ! ! READ GRID

  !OPEN(UNIT=1,FILE="x1_grid_for_masoud_conrad.in",STATUS='OLD',FORM='FORMATTED')
  OPEN(UNIT=1,FILE="./x1_grid.in",STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nx-1
     read(1,*) j, xu(i)
     !write(6,*) "xu(",i,") = ", xu(i)
  enddo
  close(1)
  xc(2:nx-1) = .5*(xu(1:nx-2)+xu(2:nx-1))
  xc(1 ) = 2.*xu(1  )-xc(2  )
  xc(nx) = 2.*xu(nx-1)-xc(nx-1)

!  do i= 2, nx-1
!     write(6,*) "xc(",i,") = ", xc(i)
!  enddo


  !OPEN(UNIT=1,FILE="x2_grid_for_masoud_conrad.in",STATUS='OLD',FORM='FORMATTED')
  OPEN(UNIT=1,FILE="./x2_grid.in",STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, ny-1
     read(1,*) j, yv(i)
     !write(6,*) "yv(",i,") = ", yv(i)
  enddo
  close(1)
  yc(2:ny-1) = .5*(yv(1:ny-2)+yv(2:ny-1))
  yc(1 ) = 2.*yv(1  )-yc(2  )
  yc(ny) = 2.*yv(ny-1)-yc(ny-1)

 ! do i= 2, ny-1
 !    write(6,*) "yc(",i,") = ", yc(i)
 ! enddo


  OPEN(UNIT=1,FILE="./x3_grid.in",STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nz-1
     read(1,*) j, zwg(i)
  enddo
  close(1)
  zcg(2:nz-1) = .5*(zwg(1:nz-2)+zwg(2:nz-1))
  zcg(1  )= 2.*zwg(1  )-zcg(2  )
  zcg(nz) = 2.*zwg(nz-1)-zcg(nz-1)

  !do i= 2, nz-1
  !   write(6,*) "zcg(",i,") = ", zcg(i)
  !enddo


  close(1)
  DEALLOCATE(cug,cvg)

  write(6,*) "READ GRID DONE"

  return
end subroutine read_grid

subroutine read_restart(filename,nx,ny,nz,var,time)
  implicit none
  
  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time,DTM1,grav

  write(6,*) "nx,ny,nz = ", nx,ny,nz

  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')  
  READ(19) I,J,K,JP 
  write(6,*) "I,J,K,JP = ", I,J,K,JP 
  DO K=1,NZ
!     write(6,*) " READ K = ", K
     READ(19) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  READ(19) nstep
  READ(19) TIME
  write(6,*) 'time=',time
  READ(19) DTM1,grav
  WRITE(*,*)"READING DONE"

  
  CLOSE(19)
  write(6,*) "READING RESTART FILE DONE"
  
  return
end subroutine read_restart


subroutine read_5p(filename,nx,ny,nz,var,time,dtm1)
  implicit none
  
  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nzg,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time,DTM1,grav
  integer (kind  = 4)  :: KMIN5p,KMAX5p

  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')  
  READ(19) I,J,K,JP 
  DO K=1,NZ
!     write(6,*) " READ K = ", K
     READ(19) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  READ(19) nstep
  READ(19) TIME
!  write(6,*) 'time=',time
  READ(19) DTM1,grav
  READ(19) KMIN5P,KMAX5P,NZG
  CLOSE(19)
  write(6,*) "READING 5p FILE DONE"
  
  return
end subroutine read_5p



subroutine write_5p(filename,nx,ny,nz,var,time1,time2)
  implicit none
  
  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nzg,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time1,time2
  integer (kind  = 4)  :: KMIN5p,KMAX5p

  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')  
  WRITE(19) nx,ny,nz 
  DO K=1,NZ
!     write(6,*) " READ K = ", K
     WRITE(19) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  WRITE(19) TIME1,TIME2
  CLOSE(19)
  
  return
end subroutine write_5p

subroutine read_post(filename,nx,ny,nz,var,time1,time2)
  implicit none
  
  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nzg,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time1,time2
  integer (kind  = 4)  :: nx1,ny1,nz1
    
  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')  
  READ(19) nx1,ny1,nz1 
  DO K=1,NZ
!     write(6,*) " READ K = ", K
     READ(19) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  write(6,*) " complete"
  READ(19) TIME1,TIME2
  CLOSE(19)
  
  return
end subroutine read_post




subroutine ascii_version

  integer (kind = 4) :: i,j,k,imax,jmax,kmax
  
  open(1,file="3d.vtk", form = 'formatted', status = 'unknown')

  imax = 10
  jmax = 10
  kmax = 10

  write(1,FMT='(a26)') "# vtk DataFile Version 2.0"
  write(1,FMT='(a9)') "FlowField"
  write(1,FMT='(a5)') "ASCII"
  write(1,FMT='(a23)') "DATASET STRUCTURED_GRID"
  write(1,FMT='(a10)', advance="no") "DIMENSIONS"
  write(1,"(a1, 3i7)") " ", imax, jmax, kmax
  write(1,FMT='(a6)', advance="no") "POINTS"
  write(1,"(a1, i15)", advance="no") " ", imax*jmax*kmax
  write(1,"(a6)") " float"

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1,*) real(i), real(j), real(k)
        end do
     end do
  end do

  close(1)

end subroutine ascii_version

subroutine binary_version

  integer (kind = 4) :: i,j,k,imax,jmax,kmax,s1
  character(len=128) :: buff
  real    (kind = 4), allocatable, dimension(:,:,:) :: vtkvar
  
  open(unit=1,file='3d.vtk',access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = 10
  jmax = 10
  kmax = 10

  allocate(vtkvar(imax,jmax,kmax))

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, " float"
  write(1) buff//char(10)

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1) real(i), real(j), real(k)
           vtkvar(i,j,k) = real(i)*real(j)*real(k)
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS vtkvar float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1) vtkvar(i,j,k)
        end do
     end do
  end do

  deallocate(vtkvar)
  close(1)
end subroutine binary_version

subroutine center_velocity(nx,ny,nz,Ucen,Uin,dir)
  !Passed Variables
  integer,intent(in)      :: dir,nx,ny,nz
  real (kind = 8),intent(in)         :: Uin(nx,ny,nz) 
  real (kind = 8),intent(out)        :: Ucen(nx,ny,nz) 


  !Local Variables
  integer              	:: i,j,k,err
  !Zero Output Array
  Ucen=0.d0
  !*************************************************
  !********************X1***************************
  !*************************************************
  if(dir.EQ.1) then
     !U
     do k=2,nz-1
        write(6,*) " CENTER_VELOCITY DIR1, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i-1,j,k))
           enddo
     	enddo
     enddo

     !*************************************************
     !********************X2***************************
     !*************************************************
  elseif (dir.EQ.2) then 
     !V
     !U
     do k=2,nz-1
        write(6,*) " CENTER_VELOCITY DIR2, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j-1,k))
           enddo
   	enddo
     enddo
     !*************************************************
     !********************X3***************************
     !*************************************************
  elseif (dir.EQ.3) then 
     !W
     !U
     do k=2,nz-1
        write(6,*) " CENTER_VELOCITY DIR3, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j,k-1))
           enddo
   	enddo
     enddo

  else
     !Invalid direction
     !write(*,'(a60,i2)') "INVALID DIRECTION IN center_velocities 
     !&                      dir must be 1,2,3.  dir= ", dir
  endif

  Ucen(:,1,:) = Ucen(:,ny-1,:)
  Ucen(:,ny,:) = Ucen(:,2,:)

  return
end subroutine center_velocity




subroutine write_plt_3d_3vars(filename,u,v,w,nx,ny,nz,xc,yc,zcg)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: u(2:nx-1,2:ny-1,2:nz-1),v(2:nx-1,2:ny-1,2:nz-1),w(2:nx-1,2:ny-1,2:nz-1)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz),ru1

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx-1
  jmax_plt = ny-1
  kmax_plt = nz-1

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
  ALLOCATE(yc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
  ALLOCATE(zc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
      
  kmod = 0
  DO k = 2, kmax_plt
     DO i = 2, imax_plt
        DO j = 2, jmax_plt
           xc_3d(i,j,k) = xc(i)!cos(yc(j))
           yc_3d(i,j,k) = yc(j)!sin(yc(j))
           zc_3d(i,j,k) = zcg(k)
        ENDDO
     ENDDO
  ENDDO
  
  kmod = 0
  write(6,*) "done1"
  ier = tecini(trim(title)//nulchar,'x,y,z,u,v,w'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt-1,jmax_plt-1,kmax_plt-1, &
              'BLOCK'//nulchar,nulchar)
  write(6,*) "done2"

  ! Write out the field data.
  itot = (imax_plt-1)*(jmax_plt-1)*(kmax_plt-1)
  ier = tecdat(itot,xc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,yc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,zc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,u(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,v(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,w(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_3d_3vars"
  deallocate(xc_3d,yc_3d,zc_3d)
  return 
end subroutine write_plt_3d_3vars



subroutine write_plt_3d_1var(filename,p,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end
  integer ( kind = 4 ) :: imax,jmax,kmax,s1,kstart,kend,nz,ny,nx

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: p(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,temp

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx_end-nx_start+1
  jmax_plt = ny_end-ny_start+1
  kmax_plt = nz_end-nz_start+1
 

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(yc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(zc_3d(imax_plt,jmax_plt,kmax_plt))

  ALLOCATE(temp(imax_plt,jmax_plt,kmax_plt))
      
  kmod = 0
  DO k = nz_start, nz_end
     DO i = nx_start, nx_end
        DO j = ny_start, ny_end

           xc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = xc(i)
           yc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = yc(j)
           zc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = zcg(k)
           
           temp(i+1-nx_start,j+1-ny_start,k+1-nz_start) = p(i,j,k)
!           write(*,*)  xc(i), yc(j) ,zcg(k) , p(i,j,k)
        ENDDO
     ENDDO
!     write(*,*) 'zc:',zcg(k)
  ENDDO


  ier = tecini(trim(title)//nulchar,'x,y,z,var'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt,jmax_plt,kmax_plt, &
              'BLOCK'//nulchar,nulchar)

  ! Write out the field data.
  itot = (imax_plt)*(jmax_plt)*(kmax_plt)
  ier = tecdat(itot,xc_3d,disdouble)
  ier = tecdat(itot,yc_3d,disdouble)
  ier = tecdat(itot,zc_3d,disdouble)
  ier = tecdat(itot,temp,disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_3d_1var"



  deallocate(xc_3d,yc_3d,zc_3d,temp)

  return 
end subroutine write_plt_3d_1var



subroutine write_plt_3d_5vars(filename,u,v,w,p,d,nx,ny,nz,xc,yc,zcg)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: u(2:nx-1,2:ny-1,2:nz-1),v(2:nx-1,2:ny-1,2:nz-1),w(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: p(2:nx-1,2:ny-1,2:nz-1),d(2:nx-1,2:ny-1,2:nz-1)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz),ru1

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx-1
  jmax_plt = ny-1
  kmax_plt = nz-1

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
  ALLOCATE(yc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
  ALLOCATE(zc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))

  kmod = 0
  DO k = 2, kmax_plt
     DO i = 2, imax_plt
        DO j = 2, jmax_plt
           xc_3d(i,j,k) = xc(i)!cos(yc(j))
           yc_3d(i,j,k) = yc(j)!sin(yc(j))
           zc_3d(i,j,k) = zcg(k)
        ENDDO
     ENDDO
  ENDDO
  kmod = 0

  write(6,*) "done1"
  ier = tecini(trim(title)//nulchar,'x,y,z,u,v,w,p,d'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt-1,jmax_plt-1,kmax_plt-1, &
              'BLOCK'//nulchar,nulchar)
  write(6,*) "done2"

  ! Write out the field data.
  itot = (imax_plt-1)*(jmax_plt-1)*(kmax_plt-1)
  ier = tecdat(itot,xc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,yc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,zc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,u(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,v(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,w(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,p(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ! Close the file
  ier = tecend()

  deallocate(xc_3d,yc_3d,zc_3d)
  write(6,*) "done write_plt_3d_5vars"
  return 
end subroutine write_plt_3d_5vars




subroutine write_3var_box(filename,u,v,w,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end
  integer ( kind = 4 ) :: imax,jmax,kmax,s1,kstart,kend,nz,ny,nx

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: u(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: v(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: w(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,tempu,tempv,tempw

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx_end-nx_start+1
  jmax_plt = ny_end-ny_start+1
  kmax_plt = nz_end-nz_start+1
 

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(yc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(zc_3d(imax_plt,jmax_plt,kmax_plt))

  ALLOCATE(tempu(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(tempv(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(tempw(imax_plt,jmax_plt,kmax_plt))
      
  kmod = 0
  DO k = nz_start, nz_end
     DO i = nx_start, nx_end
        DO j = ny_start, ny_end

           xc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = xc(i)
           yc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = yc(j)
           zc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = zcg(k)
           
           tempu(i+1-nx_start,j+1-ny_start,k+1-nz_start) = u(i,j,k)
           tempv(i+1-nx_start,j+1-ny_start,k+1-nz_start) = v(i,j,k)
           tempw(i+1-nx_start,j+1-ny_start,k+1-nz_start) = w(i,j,k)
!           write(*,*)  xc(i), yc(j) ,zcg(k) , p(i,j,k)
        ENDDO
     ENDDO
  ENDDO


  ier = tecini(trim(title)//nulchar,'x,y,z,u,v,w'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt,jmax_plt,kmax_plt, &
              'BLOCK'//nulchar,nulchar)

  ! Write out the field data.
  itot = (imax_plt)*(jmax_plt)*(kmax_plt)
  ier = tecdat(itot,xc_3d,disdouble)
  ier = tecdat(itot,yc_3d,disdouble)
  ier = tecdat(itot,zc_3d,disdouble)
  ier = tecdat(itot,tempu,disdouble)
  ier = tecdat(itot,tempv,disdouble)
  ier = tecdat(itot,tempw,disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_3d_box"



  deallocate(xc_3d,yc_3d,zc_3d,tempu,tempv,tempw)

  return 
end subroutine write_3var_box



subroutine write_5var_box(filename,u,v,w,p,d,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end
  integer ( kind = 4 ) :: imax,jmax,kmax,s1,kstart,kend,nz,ny,nx

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: u(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: v(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: w(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: p(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: d(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,tempu,tempv,tempw,tp,td

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx_end-nx_start+1
  jmax_plt = ny_end-ny_start+1
  kmax_plt = nz_end-nz_start+1
 

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(yc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(zc_3d(imax_plt,jmax_plt,kmax_plt))

  ALLOCATE(tempu(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(tempv(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(tempw(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(tp(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(td(imax_plt,jmax_plt,kmax_plt))
      
  kmod = 0
  DO k = nz_start, nz_end
     DO i = nx_start, nx_end
        DO j = ny_start, ny_end

           xc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = xc(i)
           yc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = yc(j)
           zc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = zcg(k)
           
           tempu(i+1-nx_start,j+1-ny_start,k+1-nz_start) = u(i,j,k)
           tempv(i+1-nx_start,j+1-ny_start,k+1-nz_start) = v(i,j,k)
           tempw(i+1-nx_start,j+1-ny_start,k+1-nz_start) = w(i,j,k)
           tp(i+1-nx_start,j+1-ny_start,k+1-nz_start) = p(i,j,k)
           td(i+1-nx_start,j+1-ny_start,k+1-nz_start) = d(i,j,k)
        ENDDO
     ENDDO
  ENDDO


  ier = tecini(trim(title)//nulchar,'x,y,z,u,v,w,p,d'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt,jmax_plt,kmax_plt, &
              'BLOCK'//nulchar,nulchar)

  ! Write out the field data.
  itot = (imax_plt)*(jmax_plt)*(kmax_plt)
  ier = tecdat(itot,xc_3d,disdouble)
  ier = tecdat(itot,yc_3d,disdouble)
  ier = tecdat(itot,zc_3d,disdouble)
  ier = tecdat(itot,tempu,disdouble)
  ier = tecdat(itot,tempv,disdouble)
  ier = tecdat(itot,tempw,disdouble)
  ier = tecdat(itot,tp,disdouble)
  ier = tecdat(itot,td,disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_3d_box"



  deallocate(xc_3d,yc_3d,zc_3d,tempu,tempv,tempw,tp,td)

  return 
end subroutine write_5var_box


subroutine integxy(filename,p,xu,yv,zw,nx,ny,nz)
  implicit none

  integer ( kind = 4 ) :: i,j,k
  integer ( kind = 4 ) :: imax,jmax,kmax,s1,kstart,kend,nz,ny,nx

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: p(nx,ny,nz)

  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz)
  real    ( kind = 8 ) :: integ_p(nz),dx,dy,lxx,lyy

  integ_p = 0.0
   do k=2,nz
   do i=2,nx
   do j=2,ny
         dx = xu(i)-xu(i-1)
         dy = yv(j)-yv(j-1)
         integ_p(k) = integ_p(k) + p(i,j,k)*dx*dy
   enddo
   enddo
   enddo
     lxx =  xu(nx)-xu(1)
     lyy = yv(ny)-yv(1)


     OPEN(UNIT=2,FILE=adjustl(trim(filename)),STATUS='NEW',FORM='FORMATTED')
     do k=2,nz
        write(2,'(2F20.14)')  0.5*(zw(k)+zw(k-1)) , integ_p(k)/(lxx*lyy)
     enddo
     CLOSE(2)
     return
  end  subroutine integxy
