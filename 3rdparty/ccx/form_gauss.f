      subroutine gauss_roots_3d2(ref)
      real*8 ref(3,8)
      include "gauss.f"
      do i=1,3
        do j=1,8
          ref(i,j)=gauss3d2(i,j)
        enddo
      enddo
      return
      end


      subroutine gauss_roots_3d3(ref)
      real*8 ref(3,27)
      include "gauss.f"
      do i=1,3
        do j=1,27
          ref(i,j)=gauss3d3(i,j)
        enddo
      enddo
      return
      end


      subroutine gauss_weights_3d2(ref)
      real*8 ref(8)
      include "gauss.f"
      do i=1,8
        ref(i)=weight3d2(i)
      enddo
      return
      end


      subroutine gauss_weights_3d3(ref)
      real*8 ref(27)
      include "gauss.f"
      do i=1,27
        ref(i)=weight3d3(i)
      enddo
      return
      end
