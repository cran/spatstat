      subroutine cif(nmbr,u,v,mark,ix,x,y,marks,npts,nmarks,par,period,
&     tpar, cifval,aux)
      implicit double precision(a-h,o-z)
      dimension par(1), tpar(1), x(1), y(1), marks(1), period(2)
      integer aux(1)
      if(.not.(nmbr .eq. 1))goto 23000
      call strauss(u,v,ix,x,y,npts,par,period,cifval)
      if(.not.(tpar(1) .gt. 0))goto 23002
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
23002 continue
      goto 23001
23000 continue
      if(.not.(nmbr .eq. 2))goto 23004
      call straush(u,v,ix,x,y,npts,par,period,cifval)
      if(.not.(tpar(1) .gt. 0))goto 23006
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
23006 continue
      goto 23005
23004 continue
      if(.not.(nmbr.eq.3))goto 23008
      call sftcr(u,v,ix,x,y,npts,par,period,cifval)
      if(.not.(tpar(1) .gt. 0))goto 23010
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
23010 continue
      goto 23009
23008 continue
      if(.not.(nmbr.eq.4))goto 23012
      call straussm(u,v,mark,ix,x,y,marks,npts,nmarks,par,period,cifval)
      if(.not.(tpar(1) .gt. 0))goto 23014
      call lptm(u,v,mark,tpar,tval)
      cifval = cifval*tval
23014 continue
      goto 23013
23012 continue
      if(.not.(nmbr .eq. 5))goto 23016
      call straushm(u,v,mark,ix,x,y,marks,nmarks,npts,par,period,cifval)
      if(.not.(tpar(1) .gt. 0))goto 23018
      call lptm(u,v,mark,tpar,tval)
      cifval = cifval*tval
23018 continue
      goto 23017
23016 continue
      if(.not.(nmbr .eq. 6))goto 23020
      call dgs(u,v,ix,x,y,npts,par,period,cifval)
      if(.not.(tpar(1) .gt. 0))goto 23022
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
23022 continue
      goto 23021
23020 continue
      if(.not.(nmbr .eq. 7))goto 23024
      call diggra(u,v,ix,x,y,npts,par,period,cifval)
      if(.not.(tpar(1) .gt. 0))goto 23026
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
23026 continue
      goto 23025
23024 continue
      if(.not.(nmbr .eq. 8))goto 23028
      call geyer(u,v,ix,x,y,npts,par,period,cifval,aux)
      if(.not.(tpar(1) .gt. 0))goto 23030
      call lpt(u,v,tpar,tval)
      cifval = cifval*tval
23030 continue
      goto 23029
23028 continue
      call fexit("Cif number is greater than 8; bailing out.")
23029 continue
23025 continue
23021 continue
23017 continue
23013 continue
23009 continue
23005 continue
23001 continue
      return
      end
