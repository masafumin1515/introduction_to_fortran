module subprogs
  implicit none
contains
  subroutine gauss_jordan(a0, x, b, n)
    ! -- ガウス・ジョルダン法（部分pivot選択なし） --
    integer, intent(in) :: n  ! 配列の寸法
    real(8), intent(in) :: a0(n, n), b(n) ! ay = b
    real(8), intent(out) :: x(n) 
    integer i, k
    real(8) ar, a(n, n)
    a(:, :) = a0(:, :) ! a0をaにcopy
    x(:) = b(:) ! bをxにcopy
    do k = 1, n
      if (a(k, k) == 0.0d0) stop 'pivot = 0' ! pivotが0なら停止する
      ar = 1.0d0 / a(k, k)
      a(k, k) = 1.0d0
      a(k, k+1:n) = ar * a(k, k+1:n) ! k行のk+1列からn列にarをかける
      x(k) = ar * x(k)
      do i = 1, n
        if (i /= k) then
          a(i, k+1:n) = a(i, k+1:n) - a(i, k) * a(k, k+1:n)
          x(i) = x(i) - a(i, k) * x(k)
          a(i, k) = 0.0d0
        endif
      enddo
    enddo
  end subroutine gauss_jordan

  subroutine set_random_ad(a, b, x, n)
    ! nを取得 a,b,xを割付け aとbに乱数を設定
    real(8), allocatable, intent(out) :: a(:, :), b(:), x(:)
    integer  n
    write(*, '(a)', advance = 'no') ' input n : '
    read(*, *) n
    if (n < 1 .or. 100 < n) stop 'n must be 0 < n < 101'
    allocate (a(n, n), b(n), x(n))
    call random_number(a)
    call random_number(b)
    print *, a
    print *, b
  end subroutine set_random_ad
end module subprogs

program main
  use subprogs
  implicit none
  real(8), allocatable :: a(:, :), b(:), x(:), r(:)
  integer n
  call set_random_ad(a, b, x, n)
  call gauss_jordan(a, x, b, n)
  allocate (r(n)) ! 残差ベクトルの内積を出力
  r(:) = b(:) - matmul(a, x)
  write(*, *) 'Gauss-Jordan error = ', dot_product(r,r)
  deallocate(a, b, x)
end program main
