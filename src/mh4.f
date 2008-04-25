C Output from Public domain Ratfor, version 1.0
      subroutine mh4(par,period,xprop,yprop,mprop,ntypes, iseed,nrep,mre
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
      call arand(iseed(1),iseed(2),iseed(3),urp)
      if(urp.gt.p)then
      call arand(iseed(1),iseed(2),iseed(3),urq)
      if(urq.gt.q)then
      u = xprop(irep)
      v = yprop(irep)
      mrk = mprop(irep)
      call straussm(u,v,mrk,m1,x,y,marks,npts,ntypes,par,period,cifval)
      anumer = cifval
      adenom = qnodds*(npts+1)
      call arand(iseed(1),iseed(2),iseed(3),bp)
      if(bp*adenom .lt. anumer)then
      npts = npts + 1
      itype = 1
      endif
      else
      if(npts .ne. 0)then
      call arand(iseed(1),iseed(2),iseed(3),xi)
      ix = 1 + int(npts*xi)
      mrki = marks(ix)
      call straussm(x(ix),y(ix),mrki,ix,x,y,marks,npts,ntypes,par,period
     *,cifval)
      anumer = qnodds * npts
      adenom = cifval
      call arand(iseed(1),iseed(2),iseed(3),dp)
      if(dp*adenom .lt. anumer)then
      itype = 2
      endif
      endif
      endif
      else
      if(npts .ne. 0)then
      call arand(iseed(1),iseed(2),iseed(3),xi)
      ix = 1 + int(npts*xi)
      u = xprop(irep)
      v = yprop(irep)
      mrki = marks(ix)
      mrk = mprop(irep)
      if(fixall .and. (mrk .ne. mrki))then
      itype = 0
      else
      call straussm(x(ix),y(ix),mrki,ix,x,y,marks,npts,ntypes,par,period
     *,cvd)
      call straussm(u,v,mrk,ix,x,y,marks,npts,ntypes,par,period,cvn)
      call arand(iseed(1),iseed(2),iseed(3),sp)
      if(sp*cvd .lt. cvn)then
      itype = 3
      endif
      endif
      endif
      endif
      if(itype .gt. 0)then
      if(itype.eq.1)then
      ix = npts
      x(ix) = u
      y(ix) = v
      marks(ix) = mrk
      if(npts .ge. npmax)then
      mrep = irep+1
      return
      endif
      else
      if(itype.eq.2)then
      call deathm(x,y,marks,npts,ix)
      else
      x(ix) = u
      y(ix) = v
      marks(ix) = mrk
      endif
      endif
      endif
      irep = irep+1
      mrep = irep
      goto 23000
      endif
23001 continue
      return
      end
