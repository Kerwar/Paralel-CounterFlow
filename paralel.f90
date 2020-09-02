program main
  
  use type_m   , only: DP 
  use param_m  , only: param_t 
  use paralel_m, only: parll_t
  use coefun_m , only: ADI, write_tstep
  use point_m  , only: point_t, coefx_1, coefx_2, coefx_3, coefx_4, coefx_5, coefx_6, &
  coefx_7, coefx_8, coefx_9, coefy_1, coefy_2, coefy_3, coefy_4, coefy_5
  
  implicit none
  
  include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'
  
  type(param_t) :: param
  type(parll_t) :: PL
  
  real(DP) :: t_cpu_start, t_cpu_end
  
  type(point_t), allocatable :: field(:,:), field0(:,:), solution(:,:,:)
  ! Message Send/Recive Left/Riht/Bot/Top
  real(DP), allocatable :: msl(:), mrl(:), msr(:), mrr(:), &
  msb(:), mrb(:), mst(:), mrt(:), soltransf(:), everyerror(:)
  
  real(DP) :: tmperror, error
  ! Error for MPI
  integer :: ierr, ROW_COMM, COL_COMM
  
  integer :: status(MPI_STATUS_SIZE)
  
  ! Integer needed for indeces
  integer :: i, j, k, l, t, meshsize, rowsize, colsize, proc, index
  ! Indeces for data treatment in proc0
  integer, allocatable :: i_str_p(:), i_end_p(:), &
  j_str_p(:), j_end_p(:)
  
  ! Initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, PL%nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, PL%myproc, ierr)
  
  if (PL%myproc == 0) call cpu_time(t_cpu_start)
  
  call param%read_param()
  
  call PL%read_parll()
  
  if (PL%myproc == 0) allocate(solution(param%nx, param%ny, param%nt), &
  everyerror(PL%nprocs))
  
  ! Domain decomposition in processors
  PL%nrows = int(sqrt(float(param%ny)/float(param%nx)* &
  float(PL%nprocs)))
  if (PL%nrows == 0) PL%nrows = 1
  
  PL%ncols = int(float(PL%nprocs)/float(PL%nrows))
  
  do while (PL%nrows * PL%ncols /= PL%nprocs)
    PL%nrows = PL%nrows - 1
    PL%ncols = int(float(PL%nprocs)/float(PL%nrows))
  end do
  
  ! Domain decompositions
  
  PL%nxadd = mod(param%nx, PL%ncols)
  PL%nyadd = mod(param%ny, PL%nrows)
  
  PL%myrow = min(int(float(PL%myproc)/float(PL%ncols)), PL%ncols-1)
  PL%mycol = PL%myproc - PL%myrow * PL%ncols 
  
  !NEED TO CHECK IF THIS DECOMPOSITION WORKS FOR MORE THAN 2 PROC
  ! Set size of every mesh for each processor
  if (PL%mycol < PL%nxadd) then
    PL%istr = PL%mycol * &
    (int(float(param%nx)/float(PL%ncols)) + 1) + 1
    PL%iend = (PL%mycol + 1) * &
    (int(float(param%nx)/float(PL%ncols)) + 1) 
  else
    PL%istr = PL%nxadd + PL%mycol * &
    int(float(param%nx)/float(PL%ncols))  + 1
    PL%iend = PL%nxadd + (PL%mycol + 1) * &
    int(float(param%nx)/float(PL%ncols)) 
  end if
  
  if (PL%myrow < PL%nyadd) then
    PL%jstr = PL%myrow * &
    (int(float(param%ny)/float(PL%nrows)) + 1) + 1
    PL%jend = (PL%myrow + 1) * &
    (int(float(param%ny)/float(PL%nrows)) + 1) 
  else
    PL%jstr = PL%nyadd + PL%myrow * &
    int(float(param%ny)/float(PL%nrows)) + 1
    PL%jend = PL%nyadd + (PL%myrow + 1) * &
    int(float(param%ny)/float(PL%nrows)) 
  end if
  
  PL%nx = PL%iend - PL%istr + 1
  PL%ny = PL%jend - PL%jstr + 1
  
  allocate(soltransf(3 * PL%nx * PL%ny))
  
  if (PL%myproc == 0) then
    allocate(i_str_p(0:PL%nprocs-1), i_end_p(0:PL%nprocs-1), & 
    j_str_p(0:PL%nprocs-1), j_end_p(0:PL%nprocs-1))
    
    do j = 0, PL%nrows-1
      do i = 0, PL%nxadd-1
        i_str_p(i+j*PL%ncols) =  i    * (param%nx/PL%ncols + 1) + 1 
        i_end_p(i+j*PL%ncols) = (i+1) * (param%nx/PL%ncols + 1) 
      enddo
      do i = PL%nxadd, PL%ncols-1
        i_str_p(i+j*PL%ncols) = PL%nxadd + &
        i      * (param%nx/PL%ncols) + 1
        i_end_p(i+j*PL%ncols) = PL%nxadd + &
        (i + 1) * (param%nx/PL%ncols) 
      enddo
    enddo
    
    do i = 0, PL%ncols -1
      do j = 0, PL%nyadd-1
        j_str_p(j*PL%ncols+i) =  j    * (param%ny/PL%nrows + 1) + 1
        j_end_p(j*PL%ncols+i) = (j+1) * (param%ny/PL%nrows + 1) 
      enddo
      do j = PL%nyadd, PL%nrows-1
        j_str_p(j*PL%ncols+i) = PL%nyadd + &
        j     * (param%ny/PL%nrows) + 1
        j_end_p(j*PL%ncols+i) = PL%nyadd + &
        (j+1) * (param%ny/PL%nrows)
      enddo
    enddo
  end if
  !print *, PL%myproc, PL%myrow, PL%nrows, PL%mycol, &
  !PL%ncols, PL%istr, PL%iend, PL%jstr, PL%jend
  
  !Make Row and COl communicators in order to have easier comms
  
  call MPI_COMM_SPLIT(&
  MPI_COMM_WORLD, PL%myrow, PL%mycol, ROW_COMM, ierr)
  call MPI_COMM_RANK (ROW_COMM, PL%myid_row , ierr)
  call MPI_COMM_SIZE (ROW_COMM, PL%nodes_row, ierr)
  
  call MPI_COMM_SPLIT(&
  MPI_COMM_WORLD, PL%mycol, PL%myrow, COL_COMM, ierr)
  call MPI_COMM_RANK (COL_COMM, PL%myid_col , ierr)
  call MPI_COMM_SIZE (COL_COMM, PL%nodes_col, ierr)
  
  !print *, "Comss done!"
  !Set IDS of neighbours
  PL%myleft = PL%myid_col - 1
  if (PL%myleft < 0) PL%myleft = MPI_PROC_NULL
  
  PL%myright = PL%myid_col + 1
  if (PL%myright == PL%nodes_col) PL%myright = MPI_PROC_NULL
  
  ! We impose simmetry
  if (PL%nrows /= 1) then
    PL%mybot = PL%myid_row - 1
    if (PL%mybot < 0) PL%mybot = PL%nodes_row - 1
    
    PL%mytop = PL%myid_row + 1
    if (PL%mytop == PL%nodes_row) PL%mytop = 0
  else 
    PL%mybot = MPI_PROC_NULL
    
    PL%mytop = MPI_PROC_NULL
  end if
  
  if (PL%myproc == 0) call cpu_time(t_cpu_end)
  
  allocate(field(0:PL%nx+1, 0:PL%ny+1))
  allocate(field0(0:PL%nx+1, 0:PL%ny+1))
  
  !print *, "Memory allocated in the ", PL%myproc, " processor!"
  do i = 1, PL%nx
    field(i,0)%aT(:) = 0.0_DP
    field(i,0)%aF(:) = 0.0_DP
    field(i,0)%aZ(:) = 0.0_DP
    
    field(i,0)%aT(5) = 1.0_DP
    field(i,0)%aF(5) = 1.0_DP
    field(i,0)%aZ(5) = 1.0_DP
    do j = 1, PL%ny
      field(0,j)%aT(:) = 0.0_DP
      field(0,j)%aF(:) = 0.0_DP
      field(0,j)%aZ(:) = 0.0_DP
      
      field(0,j)%aT(5) = 1.0_DP
      field(0,j)%aF(5) = 1.0_DP
      field(0,j)%aZ(5) = 1.0_DP
      
      call field(i,j)%set_coef(param, i + PL%istr - 1, &
      j + PL%jstr - 1)
      if (PL%myproc == 1 ) then
        open(102, file = "testing.txt", position = "append") 
        if (associated(field(i,j)%setXc, coefx_1)) write(102,*) "1x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_2)) write(102,*) "2x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_3)) write(102,*) "3x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_4) )write(102,*) "4x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_5)) write(102,*) "5x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_6)) write(102,*) "6x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_7)) write(102,*) "7x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_8)) write(102,*) "8x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setXc, coefx_9)) write(102,*) "9x", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setYc, coefy_1)) write(102,*) "1y", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setYc, coefy_2)) write(102,*) "2y", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setYc, coefy_3)) write(102,*) "3y", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setYc, coefy_4)) write(102,*) "4y", i + PL%istr - 1, &
        j + PL%jstr - 1
        if (associated(field(i,j)%setYc, coefy_5)) write(102,*) "5y", i + PL%istr - 1, &
        j + PL%jstr - 1
        
        close(102)
      end if
      !if (PL%myproc == 0) print *, i,j, field(i,j)%T, field(i,j)%F, field(i,j)%Z
      
      field(PL%nx,j)%aT(:) = 0.0_DP
      field(PL%nx,j)%aF(:) = 0.0_DP
      field(PL%nx,j)%aZ(:) = 0.0_DP
      
      field(PL%nx,j)%aT(5) = 1.0_DP
      field(PL%nx,j)%aF(5) = 1.0_DP
      field(PL%nx,j)%aZ(5) = 1.0_DP
    end do
    
    field(i,PL%ny)%aT(:) = 0.0_DP
    field(i,PL%ny)%aF(:) = 0.0_DP
    field(i,PL%ny)%aZ(:) = 0.0_DP
    
    field(i,PL%ny)%aT(5) = 1.0_DP
    field(i,PL%ny)%aF(5) = 1.0_DP
    field(i,PL%ny)%aZ(5) = 1.0_DP
  end do
  
  allocate(msl(3 * PL%ny), mrl(3 * PL%ny))
  allocate(msr(3 * PL%ny), mrr(3 * PL%ny))
  allocate(msb(3 * PL%nx), mrb(3 * PL%nx))
  allocate(mst(3 * PL%nx), mrt(3 * PL%nx))
  
  print *, "My Proc:", PL%myproc, "My Row:", PL%myrow, "My Col:",PL%mycol, &
  "Istr:", PL%istr, "Iend:", PL%iend, "Jstr:",PL%jstr, "Jend:", PL%jend
  
  !print *, "Let the iter begin in the ", PL%myproc, " processor!"
  do t = 1, param%nt
    
    print *, "Processor:", PL%myproc, param%nt, PL%nx, PL%ny
    do i = 1, PL%nx
      do j = 1, PL%ny
        if (.not. associated(field(i,j)%setXc)) print *, i, j
        call field(i,j)%setXc(param)
        call field(i,j)%setYc(param)
      end do
    end do
    
    print *, "Ok let's try it", param%t0hs, param%r0hs
    !print *, "Coeficien functions done in the ", PL%myproc, " processor!"
    call ADI(field(:,:)%aT(1), field(:,:)%aT(2), field(:,:)%aT(3), &
    field(:,:)%aT(4), field(:,:)%aT(5), field(:,:)%sT, field(:,:)%T, &
    tmperror, 4, param%valpha)
    error = max(tmperror, 0.0_DP)
    
    call ADI(field(:,:)%aF(1), field(:,:)%aF(2), field(:,:)%aF(3), &
    field(:,:)%aF(4), field(:,:)%aF(5), field(:,:)%sF, field(:,:)%F, &
    tmperror, 4, param%valpha)
    error = max(tmperror, error)
    
    call ADI(field(:,:)%aZ(1), field(:,:)%aZ(2), field(:,:)%aZ(3), &
    field(:,:)%aZ(4), field(:,:)%aZ(5), field(:,:)%sZ, field(:,:)%Z, &
    tmperror, 4, param%valpha)
    error = max(tmperror, error)
    
    !print *, "System solved in the ", PL%myproc, " processor!, ", error
    
    do i = 0,1
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      if (mod(PL%myid_col, 2) == i) then
        do j = 1, PL%ny
          msl(j) = field(1,j)%T
          msl(j + PL%ny) = field(1,j)%F
          msl(j + PL%ny) = field(1,j)%Z
        end do
        
        call MPI_SEND(msl, size(msl), MPI_DOUBLE_PRECISION, &
        PL%myleft, 0, ROW_COMM, ierr)
      else
        call MPI_RECV(mrr, size(mrr), MPI_DOUBLE_PRECISION, &
        PL%myright, 0, ROW_COMM, status,  ierr)
        
        do j = 1, PL%ny
          field(PL%ny + 1,j)%T = mrr(j)
          field(PL%ny + 1,j)%F = mrr(j + PL%ny)
          field(PL%ny + 1,j)%Z = mrr(j + PL%ny)
        end do
      end if
      
      if (mod(PL%myid_col, 2) == i) then
        do j = 1, PL%ny
          msr(j) = field(PL%nx,j)%T
          msr(j + PL%ny) = field(PL%ny,j)%F
          msr(j + PL%ny) = field(PL%ny,j)%Z
        end do
        
        call MPI_SEND(msr, size(msr), MPI_DOUBLE_PRECISION, &
        PL%myright, 1, ROW_COMM, ierr)
      else
        call MPI_RECV(mrl, size(mrl), MPI_DOUBLE_PRECISION, &
        PL%myleft, 1, ROW_COMM, status, ierr)
        
        do j = 1, PL%ny
          field(0,j)%T = mrl(j)
          field(0,j)%F = mrl(j + PL%ny)
          field(0,j)%Z = mrl(j + PL%ny)
        end do
      end if
      
      if (mod(PL%myid_row, 2) == i) then
        do j = 1, PL%nx
          msb(j) = field(j,1)%T
          msb(j + PL%nx) = field(j,1)%F
          msb(j + PL%nx) = field(j,1)%Z
        end do
        
        call MPI_SEND(msb, size(msb), MPI_DOUBLE_PRECISION, &
        PL%mybot, 2, COL_COMM, ierr)
      else
        call MPI_RECV(mrt, size(mrt), MPI_DOUBLE_PRECISION, &
        PL%mytop, 2, COL_COMM, status, ierr)
        
        do j = 1, PL%nx
          field(j,PL%ny + 1)%T = mrt(j)
          field(j,PL%ny + 1)%F = mrt(j + PL%nx)
          field(j,PL%ny + 1)%Z = mrt(j + PL%nx)
        end do
      end if
      
      if (mod(PL%myid_row, 2) == i) then
        do j = 1, PL%nx
          mst(j) = field(j,PL%ny)%T
          mst(j + PL%nx) = field(j,PL%ny)%F
          mst(j + PL%nx) = field(j,PL%ny)%Z
        end do
        
        call MPI_SEND(mst, size(mst), MPI_DOUBLE_PRECISION, &
        PL%mytop, 3, COL_COMM, ierr)
      else
        call MPI_RECV(mrb, size(mrb), MPI_DOUBLE_PRECISION, &
        PL%mybot, 3, COL_COMM, status, ierr)
        do j = 1, PL%nx
          field(j,0)%T = mrb(j)
          field(j,0)%F = mrb(j + PL%nx)
          field(j,0)%Z = mrb(j + PL%nx)
        end do
      end if
      
    end do
    
    !print *, "Info sent and recieved in the ", PL%myproc, " processor!"
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    if (PL%myproc /= 0) then
      !print *, "Processor ", PL%myproc, " is up to duty"
      do j = 1, PL%ny
        do i = 1, PL%nx
          index = (j - 1) * PL%nx + i
          soltransf(                    index) = field(i,j)%T
          soltransf(    PL%nx * PL%ny + index) = field(i,j)%F
          soltransf(2 * PL%nx * PL%ny + index) = field(i,j)%Z
        end do
      end do
      print *, "Actual size of ", PL%myproc, ":", size(soltransf)
      !print *, "Sending info to main proc"
      !print *, soltransf
      call MPI_SsEND(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
      0, 1000 + PL%myproc, MPI_COMM_WORLD, ierr)
      !print *, "Proc" , PL%myproc, "done sending info to main prco"
    else 
      !print *, "Processor 0 is up to duty"
      
      do i = 0, PL%nrows-1
        do j = 0, PL%ncols-1
          
          proc = i * PL%nrows + j
          !print *, "We are in proc", proc, "if u ask proc 0", PL%nrows, PL%ncols
          
          rowsize = i_end_p(proc) - i_str_p(proc) + 1
          colsize = j_end_p(proc) - j_str_p(proc) + 1
          meshsize = rowsize * colsize
          
          if (i == 0 .and. j == 0) then 
            do k = 1,rowsize
              do l = 1,colsize
                
                solution(k ,l, t)%T = field(k,l)%T
                solution(k ,l, t)%F = field(k,l)%F
                solution(k ,l, t)%Z = field(k,l)%Z
              end do
            end do
            
          else
            
            deallocate(soltransf)
            allocate(soltransf( 3 * meshsize))
            print *, "Main proc idea of size of proc ", proc, ":", size(soltransf) 
            call MPI_RECV(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
            proc, 1000 + proc, MPI_COMM_WORLD, status, ierr)            
            !print *, "We are chilling in proc 0"
            do l = 1,colsize
              do k = 1,rowsize
                index = (l - 1) * rowsize + k
                
                solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%T = &
                soltransf(               index)
                solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%F = &
                soltransf(    meshsize + index)
                solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%Z = &
                soltransf(2 * meshsize + index)
              end do
            end do
          end if
          
          
        end do
      end do
    end if
    if (PL%myproc == 0) call write_tstep(solution(:,:,t), param, t)
    
    field0 = field
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    if (PL%myproc == 0) call write_tstep(solution(:,:,t), param, t)
  end do
  call MPI_FINALIZE(ierr)
end program main

! subroutine send(field, PL, msg, pos)
!   use type_m   , only: DP 
!   use param_m  , only: param_t 
!   use paralel_m, only: parll_t
!   use coefun_m
!   use point_m  , only: point_t
!   implicit none
!   class(point_t), intent(inout) :: field(:,:)
!   type(parll_t) , intent(in   ) :: PL
!   real(DP)      , intent(inout) :: msg(:)
!   character     , intent(in   ) :: pos
!   integer :: i, j, ierr

!   select case(pos)
!   case("l")
!     do j = 1, PL%ny
!       msg(j) = field(1,j)%T
!       msg(j + PL%nx) = field(1,j)%F
!       msg(j + PL%nx) = field(1,j)%Z
!     end do

!     call MPI_SEND(msg, size(msg), MPI_DOUBLE_PRECISION, PL%myleft, &
!     ROW_COMM, ierr)

!   case("r")
!     do j = 1, PL%ny
!       msg(j) = field(PL%nx,j)%T
!       msg(j + PL%nx) = field(PL%nx,j)%F
!       msg(j + PL%nx) = field(PL%nx,j)%Z
!     end do

!     call MPI_SEND(msg, size(msg), MPI_DOUBLE_PRECISION, PL%myright, &
!     ROW_COMM, ierr)

!   case("b")
!     do i = 1, PL%nx
!       msg(i) = field(i,1)%T
!       msg(i + PL%nx) = field(i,1)%F
!       msg(i + PL%nx) = field(i,1)%Z
!     end do

!     call MPI_SEND(msg, size(msg), MPI_DOUBLE_PRECISION, PL%mybot, &
!     COL_COMM, ierr)

!   case("t")
!     do i = 1, PL%nx
!       msg(i) = field(i,PL%ny)%T
!       msg(i + PL%nx) = field(i,PL%ny)%F
!       msg(i + PL%nx) = field(i,PL%ny)%Z
!     end do

!     call MPI_SEND(msg, size(msg), MPI_DOUBLE_PRECISION, PL%mytop, &
!     COL_COMM, ierr)
!   end select 

! end subroutine send

! subroutine recv(field, PL, msg, pos)
!   use type_m   , only: DP 
!   use param_m  , only: param_t 
!   use paralel_m, only: parll_t
!   use coefun_m
!   use point_m  , only: point_t
!   implicit none
!   class(point_t), intent(inout) :: field(:,:)
!   type(parll_t) , intent(in   ) :: PL
!   real(DP)      , intent(inout) :: msg(:)
!   character     , intent(in   ) :: pos
!   integer :: i, j

!   select case(pos)
!   case("l")
!     call MPI_RECV(msg, size(msg), MPI_DOUBLE_PRECISION, PL%myleft, &
!     ROW_COMM, ierr)

!     do j = 1, PL%ny
!       field(0,j)%T = msg(j)
!       field(0,j)%F = msg(j + PL%ny)
!       field(0,j)%Z = msg(j + PL%ny)
!     end do

!   case("r")
!     call MPI_RECV(msg, size(msg), MPI_DOUBLE_PRECISION, PL%myright, &
!     ROW_COMM, ierr)

!     do j = 1, PL%ny
!       field(PL%nx + 1,j)%T = msg(j)
!       field(PL%nx + 1,j)%F = msg(j + PL%ny)
!       field(PL%nx + 1,j)%Z = msg(j + PL%ny)
!     end do

!   case("b")
!     call MPI_RECV(msg, size(msg), MPI_DOUBLE_PRECISION, PL%mybot, &
!     COL_COMM, ierr)

!     do i = 1, PL%nx
!       field(i,0)%T = msg(i)
!       field(i,0)%F = msg(i + PL%nx)
!       field(i,0)%Z = msg(i + PL%nx)
!     end do
!   case("t")
!     call MPI_RECV(msg, size(msg), MPI_DOUBLE_PRECISION, PL%mytop, &
!     COL_COMM, ierr)

!     do i = 1, PL%nx
!       field(i,PL%ny + 1)%T = msg(i)
!       field(i,PL%ny + 1)%F = msg(i + PL%nx)
!       field(i,PL%ny + 1)%Z = msg(i + PL%nx)
!     end do
!   end select 

! end subroutine recv

