program pnpreader

  use, intrinsic :: iso_c_binding

  IMPLICIT NONE

  INTEGER :: chp_handle = 100
  CHARACTER(len=200) :: diagdir1,diagdir2,chp_name1,chp_name2
  COMPLEX(C_DOUBLE_COMPLEX) :: temp_small1(1:8,1:16,1:16),temp_small2(1:8,1:16,1:16)
  COMPLEX(C_DOUBLE_COMPLEX) :: temp_big1(1:13,1:24,1:24),temp_big2(1:13,1:24,1:24)

  diagdir1 = "/pscratch/sd/e/echansen/DNA2MHDruns/mpifft1"
  diagdir2 = "/pscratch/sd/e/echansen/DNA2MHDruns/mpifft3"

  chp_name1 = "/padding.dat"
  chp_name2 = "/unpack.dat"

  OPEN(unit=chp_handle,file=trim(diagdir1)//trim(chp_name1),&
       form='unformatted', status='unknown',access='stream')
  READ(chp_handle) temp_big1
  CLOSE(chp_handle)  

  OPEN(unit=chp_handle,file=trim(diagdir2)//trim(chp_name1),&
       form='unformatted', status='unknown',access='stream')
  READ(chp_handle) temp_big2
  CLOSE(chp_handle)

  print *, "Read Temp bigs"

  OPEN(unit=chp_handle,file=trim(diagdir1)//trim(chp_name2),&
       form='unformatted', status='unknown',access='stream')
  READ(chp_handle) temp_small1
  CLOSE(chp_handle)

  OPEN(unit=chp_handle,file=trim(diagdir2)//trim(chp_name2),&
       form='unformatted', status='unknown',access='stream')
  READ(chp_handle) temp_small2
  CLOSE(chp_handle)

  print * , "Maxval Diff tempbigs",maxval(abs(temp_big1(:,:,:)-temp_big2(:,:,:))),maxloc(abs(temp_big1(:,:,:)-temp_big2(:,:,:)))

  print * , "Maxval Diff TempSmalls",maxval(abs(temp_small1-temp_small2)),maxloc(abs(temp_small1-temp_small2))

  print *, temp_small1(6,14,16)
  print *, temp_small2(6,14,16)

  print *, temp_big1(2,24,24)
  print *, temp_big2(2,24,24)

  print *, "Maxval tempbig1",maxval(abs(temp_big1)),maxloc(abs(temp_big1))
  print *, "Maxval tempbig2",maxval(abs(temp_big2)),maxloc(abs(temp_big2))

  
end program pnpreader
