module data_m
  
  use type_m , only: DP 
  use point_m, only: point_t
  use param_m, only: param_t
  
  implicit none
  
  contains
  
  subroutine ploterror(u_error, v_error, p_error)
    real(DP), intent(in   ) :: u_error, v_error, p_error
    integer, save :: it
    real(DP), save :: u_max, v_max, p_max
    logical :: exist
    
    inquire(file = "err.dat", exist = exist)
    ! t x y u v p rho
    
    if (exist) then
      open (444, file = "err.dat", position = "append")
      u_max = max(u_max, u_error)
      v_max = max(v_max, v_error)
      p_max = max(p_max, p_error)
    else
      open (444, file = "err.dat")
      it = 0
      u_max = u_error
      v_max = v_error
      p_max = p_error
      
    end if
    
    it = it + 1
    write(444, "(I6, 3ES25.14)") it, u_error/u_max, v_error/v_max, p_error/p_max
    close(444)
  end subroutine ploterror
  
  subroutine plotinfo(points, param, t)
    real(DP)     , intent(in   ) :: t
    type(point_t), intent(in   ) :: points(:,:)
    type(param_t), intent(  out) :: param
    integer :: i,j
    
    logical :: exist
    
    inquire(file = "sol.dat", exist = exist)
    ! t x y u v p rho
    
    if (exist) then
      open (333, file = "sol.dat", position = "append")
      write(333,*)
      write(333,*)
    else
      open (333, file = "sol.dat")
      write(333, "(A1,A19,6A20)") "#", "Time", "X", "Y", "U", &
      "V", "Pressure", "Density"
    end if
    
    do i = 1, param%nx
      do j = 1,param%ny
        
        write(333,"(9ES20.10E4)") t, points(i,j)%x ,points(i,j)%y, &
        points(i,j)%u%c, points(i,j)%v%c, &
        points(i,j)%p%c, points(i,j)%rho, points(i,j)%uw, points(i,j)%vs
        
      end do
    end do
    
    close(333)
    
  end subroutine plotinfo
  
end module data_m