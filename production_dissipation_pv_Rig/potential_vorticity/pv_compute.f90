Program Main
  implicit none
  include 'header'

!-------------------Input parameters -------------------------------------------------------------------
 logical :: full_box

  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx)

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imax,jmax,kmax,tmp,kstart,kend,s1
  integer ( kind = 4 ) :: tag,nx_write,ny_write,nz_write
  real    ( kind = 8 ) :: time0,time1,time2,DTM1,grav,Re
  character(len=128)   :: buff,filename

  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: u1,v1,w1,p1,dens1,pv_box
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x1(nx,ny,nz)
  real    ( kind = 8 ) :: pv(nx,ny,nz)
  logical :: save_pv 

 save_pv = .TRUE.
 full_box = .FALSE.

 omg_x = 0.0 ; omg_y = 0 ; omg_z = 0.0 

  write(*,*) 'reading grid ..'
  call read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nz,ru,rp,tag)

  allocate(w1(nx,ny,nz),u1(nx,ny,nz),v1(nx,ny,nz),p1(nx,ny,nz),dens1(nx,ny,nz))
  allocate(pv_box(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end))
  

    filename = './u_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,u1,time1)
    filename = './v_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,v1,time1)
    filename = './w_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,w1,time1)
    filename = './p_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,p1,time1)
    filename = './dens_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,dens1,time1)
 


    write(*,*) 'finished reading  restart files'

! Calculated the 3 comp of omega with first restart files time 
    call vorticity_3comps_cartesian(u1,v1,w1,omg_x,omg_y,omg_z,dq1x1,nx,ny,nz,xu,yv,zwg)
    write(6,*) "Calculation"
    deallocate(w1,u1,v1)
  
! Computed potential vorticity 
  call potential_vorticity(omg_x,omg_y,omg_z,dens1,pv,xc,yc,zcg,nx,ny,nz)
  write(*,*) 'finished computing diffusion of potential vorticity'

   do k=nz_start,nz_end
      do j=ny_start,ny_end
      do i=nx_start,nx_end
         pv_box(i,j,k) = pv(i,j,k)
      enddo
      enddo
   enddo


 

  if (save_pv)   then  
  filename = 'pv_box.plt'
  call write_plt_3d_1var(filename,pv_box,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  endif


  if (full_box)  then    
  filename = 'pv.plt'        
    call write_plt_3d_1var(filename,pv,xc,yc,zcg,2,nx-1,2,ny-1,2,nz-1,nx,ny,nz)
  endif

 deallocate(pv_box)
  !--------------------------------------------------------------------!



  stop
end Program Main

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
!     write(6,*) "zcg(",i,") = ", zcg(i)
  enddo


  ! ru(1:nx)=xu(1:nx)
  ! rp(1:nx)=xc(1:nx)

  close(1)

  write(6,*) "READ GRID DONE"

  return
end subroutine read_grid

subroutine read_restart(filename,nx,ny,nz,var,time)
  implicit none
  
  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time,DTM1,grav
  var = 0.0
  write(6,*) "nx,ny,nz = ", nx,ny,nz

  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')  
  READ(19) I,J,K,JP 
!  write(6,*) "I,J,K,JP = ", I,J,K,JP 
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





subroutine potential_vorticity(omg_x,omg_y,omg_z,dens,var,xc,yc,zcg,nx,ny,nz)
implicit none
  
!  INCLUDE 'header' 

  integer ( kind = 4 ) :: nx,ny,nz
  real ( kind = 8 ) :: ru1
  real    ( kind = 8 ) :: xu(nx),yv(ny),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: var(nx,ny,nz), dens(nx,ny,nz)
  real    (kind = 8)   :: drdx_p, drdx_m, drdy_p , drdy_m 
  real    (kind = 8)   :: drdz_p , drdz_m
  real    (kind = 8)   :: drdx,drdy,drdz,dx,dy,dz
  integer (kind = 4)   :: i,j,k  
  var =0.0

  do k=2,nz-1
!     write(6,*) " OMG: K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dx= xc(i) - xc(i-1)
           dy= yc(j) - yc(j-1)
           dz= zcg(k)- zcg(k-1)

! at u
    drdx_m = (dens(i,j,k) - dens(i-1,j,k))/dx
!at v
    drdy_m = (dens(i,j,k) - dens(i,j-1,k))/dy
!at w
    drdz_m = (dens(i,j,k) - dens(i,j,k-1))/dz      
 
    
           dx= xc(i+1) - xc(i)
           dy= yc(j+1) - yc(j)
           dz= zcg(k+1)- zcg(k)

! at u
    drdx_p = (dens(i+1,j,k) - dens(i,j,k))/dx
!at v
    drdy_p = (dens(i,j+1,k) - dens(i,j,k))/dy
!at w
    drdz_p = (dens(i,j,k+1) - dens(i,j,k))/dz      
          

    drdx = 0.5*(drdx_p + drdx_m)
    drdy = 0.5*(drdy_p + drdy_m)
    drdz = 0.5*(drdz_p + drdz_m)
 
         
    var(i,j,k) = omg_x(i,j,k)*drdx + omg_y(i,j,k)*drdy + &
                 omg_z(i,j,k)*drdz
 

   enddo
   enddo
   enddo


return
end subroutine potential_vorticity




subroutine vorticity_3comps_cartesian(u,v,w,omg_x,omg_y,omg_z,dq1x1,nx,ny,nz,xu,yv,zwg)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq2x2,dq3x3,dq1x1(nx,ny,nz)
  real    ( kind = 8 ) :: omgt,omgr,omgz

  real    ( kind = 8 ) :: xu(nx),yv(ny),zwg(nz)
  real    ( kind = 8 ) :: dx,dy,dz

  do k=2,nz-1
!     write(6,*) " OMG: K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dx= xu(i) - xu(i-1)
           dy= yv(j) - yv(j-1)
           dz= zwg(k)- zwg(k-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dy
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dx

           dq2x1=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i-1,j,k)+v(i-1,j-1,k)))/dx
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dy

           dq1x1(i,j,k)=(u(i,j,k)-u(i-1,j,k))/dx
                 
           omg_x(i,j,k) = dq3x2-dq2x3
           omg_y(i,j,k) = dq1x3-dq3x1
           omg_z(i,j,k) = dq2x1-dq1x2


        enddo
     enddo
  enddo

  return
end subroutine vorticity_3comps_cartesian


subroutine write_vtk_cartesian_3vars(filename,var1,var2,var3,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var1(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: var2(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: var3(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx-2
  jmax = ny-2
  kmax = nz-2

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, "float"
  write(1) buff//char(10)

  do k = 2,nz-1
     write(6,*) "GRID: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(xc(i)), real(yc(j)), real(zcg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS U float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var1(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS V float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var2(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS W float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var3(i,j,k))
        end do
     end do
  end do


  close(1)


end subroutine write_vtk_cartesian_3vars




subroutine write_plt_3d_3vars(filename,omg_x,omg_y,omg_z,dq1x1,dissipation,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x1(nx,ny,nz),dissipation(nx,ny,nz)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz),ru1

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,u_car,v_car,diss,w_car

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
  allocate(diss(2:nx-1,2:ny-1,2:nz-1))



  ru1 = 0.001278
         
  
  kmod = 0
  DO k = 2, kmax_plt

     write(6,*) "sudo diss: k = ", k
     DO i = 2, imax_plt
        DO j = 2, jmax_plt

           xc_3d(i,j,k) = xc(i)!cos(yc(j))
           yc_3d(i,j,k) = yc(j)!sin(yc(j))
           zc_3d(i,j,k) = zcg(k)
	   diss(i,j,k) = (ru1)*(omg_x(i,j,k)*omg_x(i,j,k) + omg_y(i,j,k)*omg_y(i,j,k) +omg_z(i,j,k)*omg_z(i,j,k))


        ENDDO
     ENDDO
  ENDDO



  write(6,*) "done1"
  ier = tecini(trim(title)//nulchar,'x,y,z,omg_x,omg_y,omg_z'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt-1,jmax_plt-1,kmax_plt-1, &
              'BLOCK'//nulchar,nulchar)
  write(6,*) "done2"

  ! Write out the field data.
  itot = (imax_plt-1)*(jmax_plt-1)*(kmax_plt-1)
  ier = tecdat(itot,xc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,yc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,zc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,omg_x(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,omg_y(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,omg_z(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  !ier = tecdat(itot,dq1x1(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  !ier = tecdat(itot,diss(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  !ier = tecdat(itot,dissipation(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)


  ! Close the file
  ier = tecend()
  deallocate(xc_3d,yc_3d,zc_3d,diss)
  write(6,*) "done write_plt_3d_3vars"
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


  ier = tecini(trim(title)//nulchar,'x,y,z,omgx,omgy,omgz'//nulchar, &
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

  
