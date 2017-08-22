module data_load
contains
use netcdf
  subroutine load_netcdf(path, file_id, ERR_FLAG)
    character, dimension(:,), intent(in) :: path
    integer, intent(out)                 :: file_id, ERR_FLAG

    ERR_FLAG = NF90_OPEN(path, NF_NOWRITE, file_id)

  end subroutine
end module
