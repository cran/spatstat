C Output from Public domain Ratfor, version 1.0
      subroutine methas(nmbr,par,period,xprop,yprop,mprop,ntypes, iseed,
     *nrep,mrep,p,q,npmax,nverb,x,y,marks,aux,npts,fixall)
      implicit double precision(a-h,o-z)
      dimension par(1), iseed(3)
      dimension x(1), y(1), marks(1), period(2)
      dimension xprop(1), yprop(1), mprop(1)
      integer aux(1)
      logical verb, marked, fixall
      one = 1.d0
      zero = 0.d0
      m1 = -1
      verb = .not.(nverb.eq.0)
      marked = ntypes .gt. 1
      qnodds = (one - q)/q
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
      if(marked)then
      mrk = mprop(irep)
      else
      mrk = -1
      endif
      call cif(nmbr,u,v,mrk,m1,x,y,marks,npts,ntypes, par,period,cifval,
     *aux)
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
      if(marked)then
      mrki = marks(ix)
      else
      mrki = -1
      endif
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y, marks,npts,ntypes,par,perio
     *d,cifval,aux)
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
      if(marked)then
      mrki = marks(ix)
      mrk = mprop(irep)
      if(fixall .and. (mrk .ne. mrki))then
      itype = 0
      else
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y,marks,npts, ntypes,par,perio
     *d,cvd,aux)
      call cif(nmbr,u,v,mrk,ix,x,y,marks,npts, ntypes,par,period,cvn,aux
     *)
      call aru(1,zero,one,iseed,sp)
      if(sp*cvd .lt. cvn)then
      itype = 3
      endif
      endif
      else
      mrki = -1
      mrk= -1
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y,marks,npts, ntypes,par,perio
     *d,cvd,aux)
      call cif(nmbr,u,v,mrk,ix,x,y,marks,npts, ntypes,par,period,cvn,aux
     *)
      call aru(1,zero,one,iseed,sp)
      if(sp*cvd .lt. cvn)then
      itype = 3
      endif
      endif
      endif
      endif
      if(itype .gt. 0)then
      if(nmbr .eq. 8)then
      call updaux(itype,x,y,u,v,npts,ix,par,period,aux)
      endif
      if(itype.eq.1)then
      ix = npts
      x(ix) = u
      y(ix) = v
      if(marked)then
      marks(ix) = mrk
      endif
      if(npts .ge. npmax)then
      mrep = irep+1
      return
      endif
      else
      if(itype.eq.2)then
      call death(x,y,marks,marked,npts,ix)
      else
      x(ix) = u
      y(ix) = v
      if(marked)then
      marks(ix) = mrk
      endif
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
