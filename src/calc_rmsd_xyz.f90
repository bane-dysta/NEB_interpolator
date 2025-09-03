! written by jxzou at 20200612: RMSD calculator (atom-atom correspondence must be assured by users)
! updated by jxzou at 20200615: coordinates of the 1st gjf file will be fixed, coordinates of
!  the 2nd gjf file will be translated and rotated.
! modified to support xyz files

! Note: The result of this subroutine is almost exactly equal to with that of VMD.
! Why 'almost'? Because when VMD use .pdb files, coordinates in .pdb have only 3 digits, causing small error

! Note: You can specify the range of atoms to be compared. Other atoms will be
! translated and rotated correspondingly.

! --------------------------------------
! How to compile:
!
!  ifort calc_rmsd_xyz.f90 -mkl -o calc_rmsd_xyz
! --------------------------------------

! MAIN program of subroutine rmsd_wrapper
program main
 implicit none
 integer :: i, idx(4)
 integer, parameter :: iout = 6
 character(len=17) :: str
 character(len=240) :: fname1, fname2
 character(len=4) :: ext1, ext2
 logical :: alive

 i = iargc()
 if(.not. (i==2 .or. i==4)) then
  write(iout,'(/,A)') ' ERROR in subroutine rmsd: wrong command line arguments!'
  write(iout,'(A)') ' Format: ./calc_rmsd_xyz file1 file2 [range1] [range2]'
  write(iout,'(/,A)') ' Example 1: ./calc_rmsd_xyz a.gjf b.gjf'
  write(iout,'(A)') ' Example 2: ./calc_rmsd_xyz a.xyz b.xyz'
  write(iout,'(A,/)') ' Example 3: ./calc_rmsd_xyz a.gjf b.xyz 1-5 2-6'
  stop
 end if

 call getarg(1,fname1)
 call getarg(2,fname2)

 ! Check file extensions
 i = index(fname1, '.', back=.true.)
 if(i > 0) then
  ext1 = fname1(i:)
  call to_lower(ext1)
 else
  ext1 = ''
 end if
 
 i = index(fname2, '.', back=.true.)
 if(i > 0) then
  ext2 = fname2(i:)
  call to_lower(ext2)
 else
  ext2 = ''
 end if

 ! Check supported formats
 if(.not. (ext1 == '.gjf' .or. ext1 == '.xyz')) then
  write(iout,'(A)') 'ERROR: unsupported file format for file1: '//TRIM(fname1)
  write(iout,'(A)') 'Supported formats: .gjf, .xyz'
  stop
 end if

 if(.not. (ext2 == '.gjf' .or. ext2 == '.xyz')) then
  write(iout,'(A)') 'ERROR: unsupported file format for file2: '//TRIM(fname2)
  write(iout,'(A)') 'Supported formats: .gjf, .xyz'
  stop
 end if

 inquire(file=TRIM(fname1),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine rmsd: file '//TRIM(fname1)//' does not exist!'
  stop
 end if
 inquire(file=TRIM(fname2),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine rmsd: file '//TRIM(fname2)//' does not exist!'
  stop
 end if

 idx = 0
 str = ' '
 if(i == 4) then
  call getarg(3,str)
  i = index(str, '-')
  if(i == 0) then
   write(iout,'(A)') "ERROR in subroutine rmsd: error range1, no '-' found."
   write(iout,'(A)') 'range1 = '//TRIM(str)
   stop
  end if
  read(str(1:i-1),*) idx(1)
  read(str(i+1:),*) idx(2)

  call getarg(4,str)
  i = index(str, '-')
  if(i == 0) then
   write(iout,'(A)') "ERROR in subroutine rmsd: error range2, no '-' found."
   write(iout,'(A)') 'range2 = '//TRIM(str)
   stop
  end if
  read(str(1:i-1),*) idx(3)
  read(str(i+1:),*) idx(4)
 end if

 if(ANY(idx<0) .or. idx(1)>idx(2) .or. idx(3)>idx(4)) then
  write(iout,'(A)') 'ERROR in subroutine rmsd: input range not valid.'
  write(iout,'(A,2I8)') 'range1=', idx(1:2)
  write(iout,'(A,2I8)') 'range2=', idx(3:4)
  stop
 end if
 if(idx(2)-idx(1) /= idx(4)-idx(3)) then
  write(iout,'(A)') 'ERROR in subroutine rmsd: number of atoms in two ranges are not equal.'
  write(iout,'(A,2I8)') 'range1=', idx(1:2)
  write(iout,'(A,2I8)') 'range2=', idx(3:4)
  stop
 end if

 call rmsd_wrapper(fname1, fname2, ext1, ext2, idx)
 stop
end program main

! convert string to lowercase
subroutine to_lower(str)
 implicit none
 integer :: i, ic
 character(len=*), intent(inout) :: str
 
 do i = 1, len(str)
  ic = iachar(str(i:i))
  if(ic >= 65 .and. ic <= 90) str(i:i) = achar(ic + 32)
 end do
end subroutine to_lower

! calculate the RMSD value of two given molecules (gjf and xyz files are supported)
subroutine rmsd_wrapper(fname1, fname2, ext1, ext2, idx)
 implicit none
 integer :: i, natom1, natom2, natom
 integer, intent(inout) :: idx(4)
 integer, parameter :: iout = 6
 real(kind=8) :: rmsd_v, trans1(3), trans2(3), rotation(3,3)
 real(kind=8), allocatable :: coor1(:,:), coor2(:,:)
 character(len=240), intent(in) :: fname1, fname2
 character(len=4), intent(in) :: ext1, ext2

 ! We need to keep coordinates in file fname1 fixed, and coordinates in file
 ! fname2 to be changed. So coordinates of fname1 is hold in coor2, and
 ! coordinates of fname2 is hold in coor1

 ! Read first file
 if(ext1 == '.gjf') then
  call read_natom_from_gjf(fname1, natom2)
  allocate(coor2(3,natom2), source=0.0d0)
  call read_coor_from_gjf(fname1, natom2, coor2)
 else if(ext1 == '.xyz') then
  call read_natom_from_xyz(fname1, natom2)
  allocate(coor2(3,natom2), source=0.0d0)
  call read_coor_from_xyz(fname1, natom2, coor2)
 end if

 ! Read second file
 if(ext2 == '.gjf') then
  call read_natom_from_gjf(fname2, natom1)
  allocate(coor1(3,natom1), source=0.0d0)
  call read_coor_from_gjf(fname2, natom1, coor1)
 else if(ext2 == '.xyz') then
  call read_natom_from_xyz(fname2, natom1)
  allocate(coor1(3,natom1), source=0.0d0)
  call read_coor_from_xyz(fname2, natom1, coor1)
 end if

 if(ALL(idx==0)) idx = [1,natom2,1,natom1] ! compare all atoms

 if(natom2 < idx(2)) then
  write(iout,'(A)') 'ERROR in subroutine rmsd_wrapper: error range1.'
  write(iout,'(A,2I8,A1,I0)') 'natom2, range1=', natom2, idx(1), '-', idx(2)
  stop
 end if

 if(natom1 < idx(4)) then
  write(iout,'(A)') 'ERROR in subroutine rmsd_wrapper: error range2.'
  write(iout,'(A,2I8,A1,I0)') 'natom1, range2=', natom1, idx(3), '-', idx(4)
  stop
 end if

 natom = idx(2) - idx(1) + 1
 call rmsd(natom, coor1(:,idx(3):idx(4)), coor2(:,idx(1):idx(2)), rmsd_v,&
           trans1, trans2, rotation)
 write(iout,'(A,F12.6)') 'RMSD = ', rmsd_v

 ! translate and rotate other atoms(not compared in RMSD) in coor1
 ! also translate compared atoms
 do i = 1, natom1, 1
  if(i<idx(3) .or. i>idx(4)) then
   coor1(:,i) = coor1(:,i) + trans1
   coor1(:,i) = MATMUL(coor1(:,i), rotation)
  end if
  coor1(:,i) = coor1(:,i) - trans2
 end do ! for i

 ! Write output file based on original format
 if(ext2 == '.gjf') then
  call prt_coor_into_gjf(fname2, natom1, coor1, .false.)
 else if(ext2 == '.xyz') then
  call prt_coor_into_xyz(fname2, natom1, coor1, .false.)
 end if

 return
end subroutine rmsd_wrapper

! find the number of atoms from a given xyz file
subroutine read_natom_from_xyz(xyzname, natom)
 implicit none
 integer :: fid
 integer, parameter :: iout = 6
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: xyzname

 open(newunit=fid,file=trim(xyzname),status='old',position='rewind')
 read(fid,*) natom
 close(fid)

 if(natom <= 0) then
  write(iout,'(a)') 'ERROR in subroutine read_natom_from_xyz: invalid natom.'
  write(iout,'(a)') 'problematic file: '//trim(xyzname)
  write(iout,'(a,i0)') 'natom = ', natom
  stop
 end if

 return
end subroutine read_natom_from_xyz

! read Cartesian coordinates from a given xyz file
subroutine read_coor_from_xyz(xyzname, natom, coor)
 implicit none
 integer :: i, fid
 integer, parameter :: iout = 6
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: xyzname
 character(len=2) :: elem

 open(newunit=fid,file=trim(xyzname),status='old',position='rewind')
 read(fid,*) i ! read number of atoms
 if(i /= natom) then
  write(iout,'(a)') 'ERROR in subroutine read_coor_from_xyz: natom mismatch.'
  write(iout,'(a)') 'problematic file: '//trim(xyzname)
  write(iout,'(a,2i0)') 'expected, actual natom = ', natom, i
  stop
 end if

 read(fid,'(a)') buf ! skip comment line
 
 do i = 1, natom, 1
  read(fid,*) elem, coor(1:3,i)
 end do ! for i

 close(fid)
 return
end subroutine read_coor_from_xyz

! print a set of Cartesian coordinates into a .xyz file
! variable replace determines whether to replace the original file
subroutine prt_coor_into_xyz(xyzname, natom, coor, replace)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2) :: elem
 character(len=240) :: buf, new_xyz
 character(len=240), intent(in) :: xyzname
 logical, intent(in) :: replace

 i = index(xyzname, '.xyz', back=.true.)
 if(i > 0) then
  new_xyz = xyzname(1:i-1)//'_new.xyz'
 else
  new_xyz = trim(xyzname)//'_new.xyz'
 end if

 open(newunit=fid1,file=TRIM(xyzname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(new_xyz),status='replace')

 read(fid1,'(A)') buf ! copy natom
 write(fid2,'(A)') TRIM(buf)
 read(fid1,'(A)') buf ! copy comment line
 write(fid2,'(A)') TRIM(buf)

 do i = 1, natom, 1
  read(fid1,*) elem
  write(fid2,'(A2,1X,3F16.8)') elem, coor(1:3,i)
 end do ! for i

 if(replace) then
  close(fid1,status='delete')
  close(fid2)
  i = RENAME(TRIM(new_xyz), TRIM(xyzname))
 else
  close(fid1)
  close(fid2)
 end if

 return
end subroutine prt_coor_into_xyz

! find the number of atoms from a given gjf file
subroutine read_natom_from_gjf(gjfname, natom)
 implicit none
 integer :: i, fid, nblank
 integer, parameter :: iout = 6
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 nblank = 0
 open(newunit=fid,file=trim(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(a)',iostat=i) buf
  if(i /= 0) exit
  if(len_trim(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(a)') 'error in subroutine read_coor_from_gjf: wrong content in gjf.'
  write(iout,'(a)') 'problematic file: '//trim(gjfname)
  stop
 end if

 read(fid,'(a)') buf ! skip charge and mult

 natom = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  natom = natom + 1
 end do ! for while

 close(fid)
 return
end subroutine read_natom_from_gjf

! read Cartesian coordinates from a given gjf file
subroutine read_coor_from_gjf(gjfname, natom, coor)
 implicit none
 integer :: i, k, fid, nblank, ncol
 integer, parameter :: iout = 6
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 character(len=2), allocatable :: elem(:)

 nblank = 0
 open(newunit=fid,file=trim(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(a)',iostat=i) buf
  if(i /= 0) exit
  if(len_trim(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(a)') 'error in subroutine read_coor_from_gjf: wrong content in gjf.'
  write(iout,'(a)') 'problematic file: '//trim(gjfname)
  stop
 end if

 read(fid,'(a)') buf ! skip charge and mult
 read(fid,'(a)') buf ! read one line of coordinates
 call detect_ncol(buf, ncol) ! test the number of columns

 BACKSPACE(fid)
 allocate(elem(natom))

 if(ncol <= 4) then
  do i = 1, natom, 1
   read(fid,*) elem(i), coor(1:3,i)
  end do ! for i
 else ! more than 4 columns
  do i = 1, natom, 1
   read(fid,*) elem(i), k, coor(1:3,i)
  end do ! for i
 end if

 close(fid)
 return
end subroutine read_coor_from_gjf

! detect the number of columns in a line
subroutine detect_ncol(buf, ncol)
 implicit none
 integer :: i, k
 integer, intent(out) :: ncol
 character(len=30) :: str(10)
 character(len=240), intent(in) :: buf

 ncol = 0
 do i = 1, 10, 1
  read(buf,*,iostat=k) str(1:i)
  if(k /= 0) exit
 end do ! for i

 ncol = i - 1
 return
end subroutine detect_ncol

! print a set of Cartesian coordinates into a .gjf file
! variable replace determines whether to replace the original file
subroutine prt_coor_into_gjf(gjfname, natom, coor, replace)
 implicit none
 integer :: i, fid1, fid2, nblank, RENAME
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2) :: elem
 character(len=240) :: buf, new_gjf
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: replace

 i = index(gjfname, '.gjf', back=.true.)
 new_gjf = gjfname(1:i-1)//'_new.gjf'
 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(new_gjf),status='replace')

 nblank = 0
 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 read(fid1,'(A)') buf ! copy charge and mult
 write(fid2,'(A)') TRIM(buf)

 do i = 1, natom, 1
  read(fid1,*) elem
  write(fid2,'(A2,1X,3F16.8)') elem, coor(1:3,i)
 end do ! for i

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(replace) then
  close(fid1,status='delete')
  close(fid2)
  i = RENAME(TRIM(new_gjf), TRIM(gjfname))
 else
  close(fid1)
  close(fid2)
 end if

 return
end subroutine prt_coor_into_gjf

! calculate the geometry center (centeroid) of a set of coordinates,
! return the corresponding translational vector
! if moved = .True., also return translated coordinates
subroutine move_coor_to_center(natom, coor, trans)
 implicit none
 integer :: i
 integer, intent(in) :: natom
 real(kind=8), intent(inout) :: coor(3,natom)
 real(kind=8), intent(out) :: trans(3)

 do i = 1, 3
  trans(i) = -SUM(coor(i,:))
 end do ! for i

 trans = trans/DBLE(natom)

 do i = 1, natom, 1
  coor(:,i) = coor(:,i) + trans
 end do ! for i

 return
end subroutine move_coor_to_center

! calculated the RMSD value of two sets of coordinates
subroutine rmsd(natom, coor1, coor2, rmsd_v, trans1, trans2, rotation)
 implicit none
 integer :: i, j, lwork
 integer, intent(in) :: natom
 integer, parameter :: iout = 6
 real(kind=8), intent(inout) :: coor1(3,natom)
 real(kind=8), intent(in) :: coor2(3,natom)
 real(kind=8), intent(out) :: trans1(3), trans2(3), rotation(3,3)
 real(kind=8), intent(out) :: rmsd_v
 real(kind=8) :: d, H(3,3), u(3,3), vt(3,3), s(3), uvt(3,3), unity(3,3)
 real(kind=8), allocatable :: work(:), new_coor1(:,:), coor(:,:)

 ! translate coor1 to its origin
 call move_coor_to_center(natom, coor1, trans1)

 ! make a copy of coor2 as coor. translate coor to its origin, and keep coor2 fixed
 allocate(coor(3,natom), source=coor2)
 call move_coor_to_center(natom, coor, trans2)

 H = 0.0d0
 call dgemm('N', 'T', 3, 3, natom, 1.0d0, coor1, 3, coor, 3, 0.0d0, H, 3)

 ! SVD on H
 lwork = 30
 allocate(work(lwork), source=0.0d0)
 ! call dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
 call dgesvd('A', 'A', 3, 3, H, 3, s, u, 3, vt, 3, work, lwork, i)
 deallocate(work)

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine rmsd: MKL dgesvd error info /= 0.'
  write(iout,'(A,I0)') 'info = ', i
  stop
 end if

 uvt = MATMUL(u, vt) ! U(V^T)

 ! d = |U(V^T)|
 d = uvt(1,1)*(uvt(2,2)*uvt(3,3) - uvt(2,3)*uvt(3,2)) - &
     uvt(1,2)*(uvt(2,1)*uvt(3,3) - uvt(2,3)*uvt(3,1)) + &
     uvt(1,3)*(uvt(2,1)*uvt(3,2) - uvt(2,2)*uvt(3,1))
 unity = 0.0d0
 forall(i=1:3) unity(i,i) = 1.0d0
 if(d < 0.0d0) unity(3,3) = -1.0d0

 ! R = UD(V^T), using array uvt to hold UD
 uvt  = MATMUL(u, unity)
 rotation = MATMUL(uvt, vt)

 ! P' = PR
 allocate(new_coor1(natom,3), source=0.0d0)
 call dgemm('T', 'N', natom, 3, 3, 1.0d0, coor1, 3, rotation, 3, 0.0d0, new_coor1, natom)
 coor1 = TRANSPOSE(new_coor1)
 deallocate(new_coor1)

 rmsd_v = 0.0d0
 do i = 1, natom, 1
  do j = 1, 3, 1
   d = coor1(j,i) - coor(j,i)
   rmsd_v = rmsd_v + d*d
  end do ! for j
 end do ! for i

 deallocate(coor)
 rmsd_v = DSQRT(rmsd_v/DBLE(natom))
 return
end subroutine rmsd