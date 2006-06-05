# $Id: fexitf.r,v 1.2 2005/03/08 20:22:50 rolf Exp $
subroutine fexit(msg)
character*(*) msg
nc = len(msg)
call fexitc(msg, nc)
end
