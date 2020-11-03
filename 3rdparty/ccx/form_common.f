      subroutine print_mat(mat, nr, nc)
      integer nr, nc, it
      real*8 mat(nr,nc)
      do it=1,nr
        write(*,*) mat(it,:)
      enddo
      write(*,*)
      return
      end
