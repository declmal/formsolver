      subroutine print_mat(mat, nr, nc, trans)
      integer nr, nc, it, trans
      real*8 mat(nr,nc)
      if (trans.eq.0) then
        do it=1,nr
          write(*,*) mat(it,1:nc)
        enddo
      else
        do it=1,nc
          write(*,*) mat(1:nr,it)
        enddo
      endif
      write(*,*)
      return
      end
