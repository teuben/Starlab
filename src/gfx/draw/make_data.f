      program make_data
c
c     Make sample data for mcdraw.
c
      parameter (N = 250, ND = 5, PI = 3.14159)
c
      open(10,file='DATA.1')
      rewind 10
      do i = 1,N
          xx = 3.*i/float(N)
          write(10,*)i,xx,sin(2.*PI*xx),cos(2.*PI*xx),
     $                    exp(-xx)*sin(2.*PI*xx)
      end do
      close(10)
c
      open(10,file='DATA.2')
      rewind 10
      do j = 1,N/2
          yy = 3.*j/float(N/2)
          do i = 1,N/2
              xx = 2.*i/float(N/2)
              write(10,*)sin(2.*PI*xx)*cos(2.*PI*yy)
          end do
      end do
      close(10)
c
      end
