subroutine gauss_jordan_pv(a0, x, b, n)
  ! -- ガウス・ジョルダン法（部分ピボットあり）
  integer, intent(in) :: n
  real(8), intent(in) :: a0(n, n), b(n)
  real(8), intent(out) :: x(n)
  integer i, k, m
  real(8) ar, am, t, a(n, n), w(n) ! a, w は作業用の自動割付行列

  a(:, :) = a0(:, :)
  x(:) = b(:)
  do i = 1, n
    ! 部分pivot選択
    m = k
    am = abs(a(k , k))
    do i = k+1, n ! a(i, k) の絶対値が最大となるm行を探す
      if (abs(a(i, k)) > am) then
        am = abs(a(i, k))
        m = i
      endif
    enddo

    if (am == 0.0d0) stop 'A is singular' ! A が特異行列なら停止する
    if (k /= m) then ! k行とm行の入れ替え
      w(k:n) = a(k, k:n)
      a(k, k:n) = a(m, k:n)
      a(m, k:n) = w(k:n)
      t = x(k)
      x(k) = x(m)
      x(m) = t
    endif
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
end subroutine gauss_jordan_pv
