C Output from Public domain Ratfor, version 1.0
      subroutine cif(nmbr,u,v,mark,ix,x,y,marks,npts,nmarks,par,period,t
     *par, cifval,aux)
      implicit double precision(a-h,o-z)
      dimension par(1), tpar(1), x(1), y(1), marks(1), period(2)
      integer aux(1)
      if(nmbr .eq. 1)then
      call strauss(u,v,ix,x,y,npts,par,period,cifval)
      if(tpar(1) .gt. 0)then
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
      endif
      else
      if(nmbr .eq. 2)then
      call straush(u,v,ix,x,y,npts,par,period,cifval)
      if(tpar(1) .gt. 0)then
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
      endif
      else
      if(nmbr.eq.3)then
      call sftcr(u,v,ix,x,y,npts,par,period,cifval)
      if(tpar(1) .gt. 0)then
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
      endif
      else
      if(nmbr.eq.4)then
      call straussm(u,v,mark,ix,x,y,marks,npts,nmarks,par,period,cifval)
      if(tpar(1) .gt. 0)then
      call lptm(u,v,mark,tpar,tval)
      cifval = cifval*tval
      endif
      else
      if(nmbr .eq. 5)then
      call straushm(u,v,mark,ix,x,y,marks,nmarks,npts,par,period,cifval)
      if(tpar(1) .gt. 0)then
      call lptm(u,v,mark,tpar,tval)
      cifval = cifval*tval
      endif
      else
      if(nmbr .eq. 6)then
      call dig1(u,v,ix,x,y,npts,par,period,cifval)
      if(tpar(1) .gt. 0)then
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
      endif
      else
      if(nmbr .eq. 7)then
      call dig2(u,v,ix,x,y,npts,par,period,cifval)
      if(tpar(1) .gt. 0)then
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
      endif
      else
      if(nmbr .eq. 8)then
      call geyer(u,v,ix,x,y,npts,par,period,cifval,aux)
      if(tpar(1) .gt. 0)then
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
      endif
      else
      call fexit("Cif number is greater than 8; bailing out.")
      endif
      endif
      endif
      endif
      endif
      endif
      endif
      endif
      return
      end
