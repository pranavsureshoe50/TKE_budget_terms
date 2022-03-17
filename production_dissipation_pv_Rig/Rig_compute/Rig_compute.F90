Program Main
  implicit none

! ------------------Input Paramters - Pranav ----------------------------------------------------------- 
  integer ( kind = 4 ), parameter :: nx = 402
  integer ( kind = 4 ), parameter :: ny = 1282
  integer ( kind = 4 ), parameter :: nz = 1282
  character(len=12), parameter :: time_res = '00000650.res'
 
  integer ( kind = 4 ), parameter :: nx_start = 2
  integer ( kind = 4 ), parameter :: nx_end = 300

  integer ( kind = 4 ), parameter :: ny_start = 550
  integer ( kind = 4 ), parameter :: ny_end = 750
  
  integer ( kind = 4 ), parameter :: nz_start = 2 
  integer ( kind = 4 ), parameter :: nz_end = 1281

!-------------------Input parameters -------------------------------------------------------------------

  logical :: full_u,full_w,full_v
  
  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx)

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imax,jmax,kmax,tmp,kstart,kend,s1
  integer ( kind = 4 ) :: tag
  real    ( kind = 8 ) :: time0,time1,time2,DTM1,grav,Re,dx,dy,dz,dvdx,dwdx
  character(len=128)   :: buff,filename
  integer (kind = 4)   :: nx_write, ny_write, nz_write

  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: u1,v1,w1,dens1
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: Rig, by_write,bz_write
  real    ( kind = 8 ) :: int2d_hl_1(nz),int2d_hr_1(nz)
  real    ( kind = 8 ) :: int2d_hl_2(nz),int2d_hr_2(nz)
  real    ( kind = 8 ) :: int2d_hl_3(nz),int2d_hr_3(nz)

  real    ( kind = 8 ) :: g, rho0, drdx, error
  logical              :: save_Rig


  REAL, DIMENSION(nx_start:nx_end)  :: x_w
  REAL, DIMENSION(ny_start:ny_end)  :: y_w
  REAL, DIMENSION(nz_start:nz_end)  :: z_w

! Input parameters - change -Pranav ---------------------------------------------- 
! --- Pranav ---- save parameters ----------------------------------------------------------
!----------------------------------------------------------------------------------
  nx_write = (nx_end - nx_start) + 1
  ny_write = (ny_end - ny_start) + 1
  nz_write = (nz_end - nz_start) + 1

  error = 1.0e-10
  g = 0.0098
  rho0 = 1.0
!
! Rig = N^2/(dvdx^2 + dwdx^2)
!

  call read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nz,ru,rp,tag)

  allocate(dens1(nx,ny,nz),v1(nx,ny,nz),w1(nx,ny,nz))

  allocate(Rig(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end))

  Rig = 0.0
  
  filename = './dens_' // time_res
  
  call read_restart(filename,nx,ny,nz,dens1,time1)

  filename = './v_' // time_res
  
  call read_restart(filename,nx,ny,nz,v1,time1)

  filename = './w_' // time_res
  
  call read_restart(filename,nx,ny,nz,w1,time1)



  do k=nz_start,nz_end
  do j=ny_start,ny_end
  do i=nx_start,nx_end

     dx = xc(i+1)-xc(i-1)     

     drdx = ( dens1(i+1,j,k)-dens1(i-1,j,k) )/dx
     dvdx = (0.5*( v1(i+1,j,k)+v1(i+1,j-1,k) ) - 0.5*( v1(i-1,j,k) + v1(i-1,j-1,k) ))/dx 
     dwdx = (0.5*( w1(i+1,j,k)+w1(i+1,j,k-1) ) - 0.5*( w1(i-1,j,k) + w1(i-1,j,k-1) ))/dx 

     Rig(i,j,k) = -(g/rho0)*drdx/(dvdx**2 + dwdx**2 + error)

     if (Rig(i,j,k).ge.10.0) Rig(i,j,k) = 10.0 
     if (Rig(i,j,k).le.-10.0) Rig(i,j,k) = -10.0 
  enddo
  enddo
  enddo
 
 

   write(6,*) "Over"
   filename = 'Rig.plt'       
   call write_plt_3d_1var(filename,Rig,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)


!--------------------------------------------------------------------!

  deallocate(dens1)

  deallocate(v1,w1,Rig)

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

  do i= 2, nx-1
!     write(6,*) "xc(",i,") = ", xc(i)
  enddo


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

  do i= 2, ny-1
!     write(6,*) "yc(",i,") = ", yc(i)
  enddo


  !OPEN(UNIT=1,FILE="x3_grid_for_masoud_conrad.in",STATUS='OLD',FORM='FORMATTED')
  OPEN(UNIT=1,FILE="./x3_grid.in",STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nz-1
     read(1,*) j, zwg(i)
     !write(6,*) "zwg(",i,") = ", zwg(i)
  enddo
  close(1)
  zcg(2:nz-1) = .5*(zwg(1:nz-2)+zwg(2:nz-1))
  zcg(1  )= 2.*zwg(1  )-zcg(2  )
  zcg(nz) = 2.*zwg(nz-1)-zcg(nz-1)

  do i= 2, nz-1
     write(6,*) "zcg(",i,") = ", zcg(i)
  enddo


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
  CLOSE(19)
  write(6,*) "READING RESTART FILE DONE"
  
  return
end subroutine read_restart

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


