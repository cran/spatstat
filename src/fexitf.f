C Output from Public domain Ratfor, version 1.0
      subroutine fexit(msg)
      integer nc
      character*(*) msg
      nc = len(msg)
      call fexitc(msg, nc)
      end
