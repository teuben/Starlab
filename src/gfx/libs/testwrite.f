      character*80 foo
      open(42,file='foo')
      write(42,*)'Record 11111'
      write(42,*)'Record 22222'
      write(42,*)'Record 33333'
      write(42,*)'Record 44444'
      write(42,*)'Record 55555'
      rewind 42
      read(42,'(a)')foo
      write(42,*)'Record xxxxx'
      do i=1,100
          read(42,'(a)',err=999,end=999)foo
      end do
999   end
