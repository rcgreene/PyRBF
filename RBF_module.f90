module RBF_mod
contains

  real(kind=8) function distance(X_1, Y_1, X_2, Y_2)
    real(kind=8), intent(in) :: X_1, Y_1, X_2, Y_2
    distance = (X_1 - X_2)**2 + (Y_1 - Y_2)**2
    distance = SQRT(distance)
  end function distance

  real(kind=8) function gaussian(dist, epsilon) result(output)
    real(kind=8), intent(in)    :: dist
    real(kind=8), intent(in)    :: epsilon
    output = exp(-epsilon*(dist**2))
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
        dist_mat(i,j) = distance(points_1(1, i), points_1(2, i), points_2(1, j), points_2(2,j))
        dist_mat(j,i) = dist_mat(i,j)
      end do
    end do
  end subroutine

  subroutine matrix_eval(mat_in, mat_out, epsilon)
    !evaluate function for all elements of a matrix
    real(kind=8), intent(in), dimension(:,:)  :: mat_in
    real(kind=8), intent(out), dimension(:,:) :: mat_out
    real(kind=8), intent(in)                  :: epsilon
    integer                                   :: i_max, i, j

    i_max = size(mat_in,1)
    do i = 1,i_max
      do j = 1,i_max
        mat_out(i,j) = gaussian(mat_in(i,j), epsilon)
      end do
    end do
  end subroutine

  subroutine RBF_matrix(points, RBF_func, coefficients, epsilon, values)
    !Determines proper coefficients for an RBF decomposition
    !over a set of 2d points
    real(kind=8), intent(in), dimension(:,:)     :: points
    !real(kind=8), intent(in)                  :: RBF_func
    real(kind=8), intent(in)                     :: epsilon
    real(kind=8), dimension(:,:), allocatable    :: dist_mat !matrix for storing distances between nodes
    real(kind=8), dimension(:,:), allocatable    :: RF_mat   !matrix for storing values of RBF
    real(kind=8), dimension(:)                   :: values   !function to be interpolated evaluated at various points (overwritten with coefficients)
    integer                                      :: point_count
    integer                                      :: SUCCESS_FLAG

    point_count = size(points, 2)
    allocate(dist_mat(point_count, point_count))
    allocate(RF_mat(point_count, point_count))
    call point_set_distance(points, points, dist_mat)
    call matrix_eval(dist_mat, RF_mat, epsilon)
    call dposv('L', point_count, 1, RF_mat, point_count, values, point_count, SUCCESS_FLAG)
  end subroutine

!  subroutine vector_matrix(points, RBF_func, coefficients, epsilon, values)
    !RBF interpolation for convergence-free vector fields
!    real(kind=8), intent(in), dimension(:,:)    :: points
!    real(kind=8), intent(in)                    :: epsilon
!    real(kind=8), dimension(:,:), allocatable   :: dist_mat
!    real(kind=8), dimension(:,:), allocatable   :: RF_mat
!    real(kind=8), dimension(:)                  :: values
!    integer                                     :: point_count
!    integer                                     :: SUCCESS_FLAG

!    point_count = size(points, 2)
!    allocate(dist_mat(point_count, point_count))
!    allocate(dist_mat(2*point_count, 2*point_count))
!    call point_set_distance(points, points, dist_mat)
!    call vector_matrix_eval(dist_mat, RF_mat, epsilon)

  subroutine load_from_sample(filename, pressure_vals, lats, lons)
    character(:), intent(in)                                :: filename
    integer                                                 :: ncid
    real(kind = 8), dimension(:,:), allocatable, intent(in) :: pressure_vals
    real(kind = 8), dimension(:),   allocatable, intent(in) :: lons, lats
    real(kind = 8),

    use netcdf

    call nf90_open(filename, NF90_NOWRITE, ncid)
    allocate(pressure_vals(201,201))
    allocate(lons(201))
    allocate(lats(201))
    call nf90_get_var(ncid, 1, pressure_vals)
    call nf90_get_var(ncid, 2, lons)
    call nf90_get_var(ncid, 3, lats)
    call nf90_close(ncid)
  end subroutine

  subroutine construct_points(array, x, y, points, vals)
    real(kind=8), intent(in), dimension(:,:)               :: array
    real(kind=8), intent(in), dimension(:,)                :: x, y
    real(kind=8), intent(out), dimension(:,:), allocatable :: points
    real(kind=8), intent(out), dimension(:,), allocatable  :: vals
  end subroutine
end module RBF_mod
