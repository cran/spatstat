# $Id: aru.r,v 1.2 2005/03/08 20:18:22 rolf Exp $
subroutine aru(n,a,b,iseed,rrr)
implicit double precision(a-h,o-z)
dimension iseed(3), rrr(n)

w = b-a
do i = 1,n {
	call arand(iseed(1),iseed(2),iseed(3),rv)
	rrr(i) = a + w*rv
}

return
end
