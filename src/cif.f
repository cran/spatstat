C Output from Public domain Ratfor, version 1.0
      subroutine cif(nmbr,u,v,mark,ix,x,y,marks,npts,ntypes,par,period, 
     *cifval,aux)
      implicit double precision(a-h,o-z)
      dimension par(1), x(1), y(1), marks(1), period(2)
      integer aux(1)
      if(nmbr .eq. 1)then
      call strauss(u,v,ix,x,y,npts,par,period,cifval)
      else
      if(nmbr .eq. 2)then
      call straush(u,v,ix,x,y,npts,par,period,cifval)
      else
      if(nmbr.eq.3)then
      call sftcr(u,v,ix,x,y,npts,par,period,cifval)
      else
      if(nmbr.eq.4)then
      call straussm(u,v,mark,ix,x,y,marks,npts,ntypes,par,period,cifval)
      else
      if(nmbr .eq. 5)then
      call straushm(u,v,mark,ix,x,y,marks,ntypes,npts,par,period,cifval)
      else
      if(nmbr .eq. 6)then
      call dgs(u,v,ix,x,y,npts,par,period,cifval)
      else
      if(nmbr .eq. 7)then
      call diggra(u,v,ix,x,y,npts,par,period,cifval)
      else
      if(nmbr .eq. 8)then
      call geyer(u,v,ix,x,y,npts,par,period,cifval,aux)
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
