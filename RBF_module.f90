module RBF_mod
contains

  real(kind=8) function dist_sq(lon_1, lat_1, lon_2, lat_2)
    real(kind=8), intent(in) :: lon_1, lat_1, lon_2, lat_2
    real(kind=8)             :: earth_rad = 6371000 !in meters
    dist_sq = (cos(lon_1*pi/180)*cos(lat_1*pi/180) - cos(lon_2*pi/180)*cos(lat_2*pi/180))**2
    dist_sq = dist_sq + (sin(lon_1*pi/180)*cos(lat_1*pi/180) - sin(lon_2*pi/180)*cos(lat_2*pi/180))**2
    dist_sq = dist_sq + (sin(lat_1*pi/180) - sin(lat_2*pi/180))**2
    dist_sq = dist_sq*(earth_rad**2)
  end function dist_sq

  subroutine dist_sq_deriv(lon_1, lat_1, lon_2, lat_2, deriv)
    real(kind=8), intent(in) :: lon_1, lat_1, lon_2, lat_2
    real(kind=8),            :: lo_1, la_1, lo_2, la_2
    real(kind=8)             :: earth_rad = 6371000 !in meters
    real(kind=8), intent(out), dimension(2) :: deriv
    la_1 = lat_1*pi/180
    lo_1 = lon_1*pi/180
    la_2 = lat_2*pi/180
    lo_2 = lon_2*pi/180
    deriv(1) = -2*(cos(lo_1)*cos(la_1) - cos(lo_2)*cos(la_2))*sin(lo_1)
    deriv(1) = deriv(1) + 2*(sin(lo_1)*cos(la_1) - sin(lo_1)*cos(la_1))*cos(lo_1)
    deriv(2) = -2*(cos(lo_1)*cos(la_1) - cos(lo_2)*cos(la_2))*cos(lo_1)*sin(la_1)
    deriv(2) = deriv(2) - 2*(sin(lo_1)*cos(la_1) - sin(lo_1)*cos(la_1))*sin(lo_1)*sin(la_1)
    deriv(2) = deriv(2) + 2*(cos(la_1) - cos(la_1))*cos(la_1)
    deriv(1) = deriv(1)*(earth_rad**2)
    deriv(2) = deriv(2)*(earth_rad**2)
  end subroutine dist_sq_deriv

  real(kind=8) function gaussian(dist, epsilon) result(output)
    real(kind=8), intent(in)    :: dist
    real(kind=8), intent(in)    :: epsilon
    output = exp(-epsilon*(dist))
  end function gaussian

  subroutine point_set_distance(points_1, points_2, dist_mat)
    !fill dist_mat with distances between point set 1, and 2
    real(kind=8), intent(in), dimension(:,:) :: points_1, points_2
    real(kind=8), intent(out), dimension(:,:):: dist_mat
    integer                                  :: i, j, i_max, j_max
    i_max = size(points_1, 2)
    j_max = size(points_2, 2)
    do i = 1,i_max
      do j = i, j_max
        dist_mat(i,j) = dist_sq(points_1(1, i), points_1(2, i), points_2(1, j), points_2(2,j))
        dist_mat(j,i) = dist_mat(i,j)
      end do
    end do
  end subroutine

  subroutine matrix_eval(mat_in, mat_out, epsilon)
    !evaluate function for all elements of a matrix
    real(kind=8), intent(in), dimension(:,:)  :: mat_in
    real(kind=8), intent(out), dimension(:,:) :: mat_out
    real(kind=8), intent(in)                  :: epsilon
    integer                                   :: i_max, i, j, dummy

    i_max = size(mat_in,1)
    do i = 1,i_max
      do j = 1,i_max
        mat_out(i,j) = gaussian(mat_in(i,j), epsilon)
      end do
    end do
    call dpotrf('L', i_max, mat_out, i_max, dummy)
  end subroutine

  subroutine RBF_matrix(points, RBF_func, coefficients, epsilon, values)
    !Determines proper coefficients for an RBF decomposition
    !over a set of 2d points
    real(kind=8), intent(in), dimension(:,:)     :: points
    !real(kind=8), intent(in)                  :: RBF_func
    real(kind=8), intent(in)                     :: epsilon
    real(kind=8)                                 :: phase
    real(kind=8), dimension(:,:), allocatable    :: dist_mat !matrix for storing distances between nodes
    real(kind=8), dimension(:,:), allocatable    :: RF_mat   !matrix for storing values of RBF
    real(kind=8), dimension(:)                   :: values   !function to be interpolated evaluated at various points (overwritten with coefficients)
    integer                                      :: point_count
    integer                                      :: SUCCESS_FLAG

    point_count = size(points, 2)
    allocate(dist_mat(point_count, point_count))
    call point_set_distance(points, points, dist_mat)
    call matrix_eval(dist_mat, RF_mat, epsilon)
    call dpotrs('L', point_count, 1, RF_mat, point_count, values, point_count, SUCCESS_FLAG)
  end subroutine


end module RBF_mod
