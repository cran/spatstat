      subroutine fexit(msg)
      integer nc
      character*(*) msg
      nc = len(msg)
      call fexitc(msg, nc)
      end
