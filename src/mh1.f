C Output from Public domain Ratfor, version 1.0
      subroutine mh1(par,period,xprop,yprop,mprop,ntypes, iseed,nrep,mre
     *p,p,q,npmax,nverb,x,y,marks,aux,npts,fixall)
      implicit double precision(a-h,o-z)
      dimension par(1), iseed(3)
      dimension x(1), y(1), marks(1), period(2)
      dimension xprop(1), yprop(1), mprop(1)
      integer aux(1)
      logical verb, fixall
      logical dummyl
      integer dummyi
      one = 1.d0
      zero = 0.d0
      m1 = -1
      verb = .not.(nverb.eq.0)
      qnodds = (one - q)/q
      dummyi = mprop(1) + ntypes + marks(1)
      dummyl = fixall
      dummyd = aux(1)
      dummyi = dummyi + 1
      dummyl = .not. dummyl
      dummyd = dummyd * dummyd
      irep = mrep
23000 if(irep .le. nrep)then
      if(verb .and. mod(irep,nverb).eq.0)then
      iprt = irep/nverb
      call intpr('irep/nverb=',-1,iprt,1)
      endif
      itype = 0
      call aru(1,zero,one,iseed,urp)
      if(urp.gt.p)then
      call aru(1,zero,one,iseed,urq)
      if(urq.gt.q)then
      u = xprop(irep)
      v = yprop(irep)
      call strauss(u,v,m1,x,y,npts,par,period,cifval)
      anumer = cifval
      adenom = qnodds*(npts+1)
      call aru(1,zero,one,iseed,bp)
      if(bp*adenom .lt. anumer)then
      npts = npts + 1
      itype = 1
      endif
      else
      if(npts .ne. 0)then
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      call strauss(x(ix),y(ix),ix,x,y,npts,par,period,cifval)
      anumer = qnodds * npts
      adenom = cifval
      call aru(1,zero,one,iseed,dp)
      if(dp*adenom .lt. anumer)then
      itype = 2
      endif
      endif
      endif
      else
      if(npts .ne. 0)then
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      u = xprop(irep)
      v = yprop(irep)
      call strauss(x(ix),y(ix),ix,x,y,npts,par,period,cvd)
      call strauss(u,v,ix,x,y,npts,par,period,cvn)
      call aru(1,zero,one,iseed,sp)
      if(sp*cvd .lt. cvn)then
      itype = 3
      endif
      endif
      endif
      if(itype .gt. 0)then
      if(itype.eq.1)then
      ix = npts
      x(ix) = u
      y(ix) = v
      if(npts .ge. npmax)then
      mrep = irep+1
      return
      endif
      else
      if(itype.eq.2)then
      call deathu(x,y,npts,ix)
      else
      x(ix) = u
      y(ix) = v
      endif
      endif
      endif
      irep = irep+1
      mrep = irep
      goto 23000
      endif
      return
      end
