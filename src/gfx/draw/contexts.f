
      subroutine save_context(idev,iwin)
c
c     Save the context of the specified  device/window.
c
      integer context_id
c
      id = context_id(idev,iwin,1)
c
      if (id.gt.0) then
          call save_mcpak_context(id)
          call save_mcdraw_context(id)
      end if
c
      end


      subroutine load_context(idev,iwin)
c
c     Load the stored context for the specified device/window.
c
      common /prompt/ iprompt
      integer context_id
c
      id = context_id(idev,iwin,2)
c
      if (id.gt.0) then
c
c          if (iprompt.ne.0) write(6,'(a)')
c     $    'Restoring graphics context.  Note that x, y, and z '//
c     $    'are not reloaded'
c
          call restore_mcpak_context(id)
          call restore_mcdraw_context(id)
      end if
c
      end


      integer function context_id(idev,iwin,iopt)
c
c     Return the index corresponding to the specified idev and iwin.
c     If the pair are not found, create a new entry (iopt = 1), or return
c     with a value of -1 (opt = 2).  Also return -1 if NCMAX is exceeded.
c
c     For now, allow up to 100 contexts (static storage!)
c
      parameter (NCMAX = 100)
c
      integer nc,id(NCMAX),iw(NCMAX)
      data nc/0/
c
      do ic=1,nc
          if (id(ic).eq.idev.and.iw(ic).eq.iwin) then
              context_id = ic
              return
          end if
      end do
c
      if (iopt.eq.2.or.nc.ge.NCMAX) then
          context_id = -1
      else
          nc = nc + 1
          id(nc) = idev
          iw(nc) = iwin
          context_id = nc
      end if
c
      end
