      subroutine gauss_roots_3d2(arr)

      real*8 arr(3,8)

      include "gauss.f"

      do i=1,3
        do j=1,8
          arr(i,j)=gauss3d2(i,j)
        enddo
      enddo

      return
      end
