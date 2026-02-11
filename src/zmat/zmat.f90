program zmatrix_generator
    implicit none
    
    integer, parameter :: maxatom = 1000
    integer, parameter :: maxneigh = 20
    real*8, parameter :: bohr2ang = 0.52917721067d0
    real*8, parameter :: ang2bohr = 1.0d0/bohr2ang
    real*8, parameter :: pi = 3.14159265358979d0
    real*8, parameter :: bondcrit = 1.15d0
    
    type :: atom_type
        character(len=2) :: name
        real*8 :: x, y, z
        integer :: atomic_num
    end type
    
    type(atom_type) :: atoms(maxatom)
    integer :: natom
    integer :: zmat(maxatom,3)
    integer :: nneigh(maxatom)
    integer :: neigh(maxneigh,maxatom)
    real*8 :: covalent_radii(118)
    character(len=200) :: input_file, output_file
    integer :: dot_pos
    
    call get_command_argument(1, input_file)
    if (len_trim(input_file) == 0) then
        write(*,*) 'Error: No input file specified'
        write(*,*) 'Usage: ./zmatrix_generator input.xyz'
        stop 1
    end if
    
    dot_pos = index(input_file, '.', back=.true.)
    if (dot_pos > 0) then
        output_file = input_file(1:dot_pos-1) // '.zmat'
    else
        output_file = trim(input_file) // '.zmat'
    end if
    
    call init_covalent_radii()
    call read_xyz_file()
    call generate_neighbor_list()
    call generate_zmatrix()
    call output_zmatrix()
    
contains

    subroutine init_covalent_radii()
        covalent_radii = 0.0d0
        covalent_radii(1) = 0.31d0   ! H
        covalent_radii(6) = 0.76d0   ! C
        covalent_radii(7) = 0.71d0   ! N
        covalent_radii(8) = 0.66d0   ! O
        covalent_radii(9) = 0.57d0   ! F
        covalent_radii(15) = 1.07d0  ! P
        covalent_radii(16) = 1.05d0  ! S
        covalent_radii(17) = 0.99d0  ! Cl
        covalent_radii(35) = 1.20d0  ! Br
        covalent_radii(53) = 1.39d0  ! I
    end subroutine

    integer function get_atomic_number(element)
        character(len=*), intent(in) :: element
        select case(trim(adjustl(element)))
            case('H'); get_atomic_number = 1
            case('C'); get_atomic_number = 6
            case('N'); get_atomic_number = 7
            case('O'); get_atomic_number = 8
            case('F'); get_atomic_number = 9
            case('P'); get_atomic_number = 15
            case('S'); get_atomic_number = 16
            case('Cl'); get_atomic_number = 17
            case('Br'); get_atomic_number = 35
            case('I'); get_atomic_number = 53
            case default; get_atomic_number = 6
        end select
    end function

    subroutine read_xyz_file()
        character(len=200) :: line
        character(len=2) :: element
        real*8 :: x, y, z
        integer :: i, ios
        
        open(10, file=trim(input_file), status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: Cannot open file ', trim(input_file)
            stop 1
        end if
        
        read(10,*) natom
        read(10,'(a)') line
        
        do i = 1, natom
            read(10,*) element, x, y, z
            atoms(i)%name = element
            atoms(i)%x = x * ang2bohr
            atoms(i)%y = y * ang2bohr
            atoms(i)%z = z * ang2bohr
            atoms(i)%atomic_num = get_atomic_number(element)
        end do
        
        close(10)
    end subroutine

    real*8 function atom_distance(i, j)
        integer, intent(in) :: i, j
        atom_distance = sqrt((atoms(i)%x - atoms(j)%x)**2 + &
                           (atoms(i)%y - atoms(j)%y)**2 + &
                           (atoms(i)%z - atoms(j)%z)**2)
    end function

    real*8 function atom_angle(i, j, k)
        integer, intent(in) :: i, j, k
        real*8 :: v1(3), v2(3), dot_prod, mag1, mag2
        
        v1(1) = atoms(i)%x - atoms(j)%x
        v1(2) = atoms(i)%y - atoms(j)%y
        v1(3) = atoms(i)%z - atoms(j)%z
        
        v2(1) = atoms(k)%x - atoms(j)%x
        v2(2) = atoms(k)%y - atoms(j)%y
        v2(3) = atoms(k)%z - atoms(j)%z
        
        dot_prod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        mag1 = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
        mag2 = sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2)
        
        if (mag1 < 1d-10 .or. mag2 < 1d-10) then
            atom_angle = 0.0d0
        else
            atom_angle = acos(max(-1.0d0, min(1.0d0, dot_prod/(mag1*mag2)))) * 180.0d0/pi
        end if
    end function

    real*8 function dihedral_angle(i, j, k, l)
        integer, intent(in) :: i, j, k, l
        real*8 :: v1(3), v2(3), v3(3), n1(3), n2(3)
        real*8 :: dot_prod, mag1, mag2, cross_prod
        
        v1(1) = atoms(i)%x - atoms(j)%x
        v1(2) = atoms(i)%y - atoms(j)%y
        v1(3) = atoms(i)%z - atoms(j)%z
        
        v2(1) = atoms(k)%x - atoms(j)%x
        v2(2) = atoms(k)%y - atoms(j)%y
        v2(3) = atoms(k)%z - atoms(j)%z
        
        v3(1) = atoms(l)%x - atoms(k)%x
        v3(2) = atoms(l)%y - atoms(k)%y
        v3(3) = atoms(l)%z - atoms(k)%z
        
        n1(1) = v1(2)*v2(3) - v1(3)*v2(2)
        n1(2) = v1(3)*v2(1) - v1(1)*v2(3)
        n1(3) = v1(1)*v2(2) - v1(2)*v2(1)
        
        n2(1) = v2(2)*v3(3) - v2(3)*v3(2)
        n2(2) = v2(3)*v3(1) - v2(1)*v3(3)
        n2(3) = v2(1)*v3(2) - v2(2)*v3(1)
        
        dot_prod = n1(1)*n2(1) + n1(2)*n2(2) + n1(3)*n2(3)
        mag1 = sqrt(n1(1)**2 + n1(2)**2 + n1(3)**2)
        mag2 = sqrt(n2(1)**2 + n2(2)**2 + n2(3)**2)
        
        if (mag1 < 1d-10 .or. mag2 < 1d-10) then
            dihedral_angle = 0.0d0
        else
            cross_prod = (n1(2)*n2(3) - n1(3)*n2(2))*v2(1) + &
                        (n1(3)*n2(1) - n1(1)*n2(3))*v2(2) + &
                        (n1(1)*n2(2) - n1(2)*n2(1))*v2(3)
            
            dihedral_angle = atan2(cross_prod, dot_prod) * 180.0d0/pi
        end if
    end function

    subroutine generate_neighbor_list()
        integer :: i, j
        real*8 :: dist, sum_cov
        
        nneigh = 0
        neigh = 0
        
        do i = 1, natom
            do j = 1, natom
                if (i == j) cycle
                
                dist = atom_distance(i, j)
                sum_cov = covalent_radii(atoms(i)%atomic_num) + &
                         covalent_radii(atoms(j)%atomic_num)
                
                if (dist < bondcrit * sum_cov) then
                    nneigh(i) = nneigh(i) + 1
                    if (nneigh(i) <= maxneigh) then
                        neigh(nneigh(i), i) = j
                    end if
                end if
            end do
        end do
    end subroutine

    subroutine generate_zmatrix()
        integer :: iatm, jatm, i1, i2, i3
        real*8 :: distmin, dist, tmpang
        logical :: found
        
        zmat = 0
        
        do iatm = 2, natom
            found = .false.
            do jatm = 1, iatm-1
                if (any(neigh(1:nneigh(iatm),iatm) == jatm)) then
                    i1 = jatm
                    found = .true.
                    exit
                end if
            end do
            
            if (.not. found) then
                distmin = 1d10
                do jatm = 1, iatm-1
                    dist = atom_distance(iatm, jatm)
                    if (dist < distmin) then
                        i1 = jatm
                        distmin = dist
                    end if
                end do
            end if
            zmat(iatm,1) = i1
            
            if (iatm >= 3) then
                found = .false.
                do jatm = 1, iatm-1
                    if (jatm == i1) cycle
                    
                    if (any(neigh(1:nneigh(i1),i1) == jatm)) then
                        if (natom > 3) then
                            tmpang = atom_angle(iatm, i1, jatm)
                            if (tmpang < 0.5d0 .or. tmpang > 179.5d0) cycle
                        end if
                        i2 = jatm
                        found = .true.
                        exit
                    end if
                end do
                
                if (.not. found) then
                    distmin = 1d10
                    do jatm = 1, iatm-1
                        if (jatm == i1) cycle
                        
                        dist = atom_distance(i1, jatm)
                        if (dist < distmin) then
                            if (natom > 3) then
                                tmpang = atom_angle(iatm, i1, jatm)
                                if (tmpang < 0.5d0 .or. tmpang > 179.5d0) cycle
                            end if
                            i2 = jatm
                            distmin = dist
                        end if
                    end do
                    
                    if (distmin == 1d10) then
                        write(*,*) 'Error: Cannot find suitable angle reference for atom', iatm
                        stop 1
                    end if
                end if
                zmat(iatm,2) = i2
                
                if (iatm >= 4) then
                    found = .false.
                    do jatm = 1, iatm-1
                        if (jatm == i1 .or. jatm == i2) cycle
                        
                        if (any(neigh(1:nneigh(i2),i2) == jatm)) then
                            tmpang = atom_angle(i1, i2, jatm)
                            if (tmpang < 0.5d0 .or. tmpang > 179.5d0) cycle
                            i3 = jatm
                            found = .true.
                            exit
                        end if
                    end do
                    
                    if (.not. found) then
                        distmin = 1d10
                        do jatm = 1, iatm-1
                            if (jatm == i1 .or. jatm == i2) cycle
                            
                            dist = atom_distance(i2, jatm)
                            if (dist < distmin) then
                                tmpang = atom_angle(i1, i2, jatm)
                                if (tmpang < 0.5d0 .or. tmpang > 179.5d0) cycle
                                i3 = jatm
                                distmin = dist
                            end if
                        end do
                        
                        if (distmin == 1d10) then
                            write(*,*) 'Error: Cannot find suitable dihedral reference for atom', iatm
                            stop 1
                        end if
                    end if
                    zmat(iatm,3) = i3
                end if
            end if
        end do
    end subroutine

    subroutine output_zmatrix()
        integer :: i, i1, i2, i3, ios
        real*8 :: bond_dist, bond_ang, dihed_ang
        
        open(20, file=trim(output_file), status='replace', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: Cannot create output file ', trim(output_file)
            stop 1
        end if
        
        write(20,'(a)') '# Z-Matrix format'
        write(20,'(a)') '# Atom  Bond_ref  Bond(A)  Angle_ref  Angle(deg)  Dihedral_ref  Dihedral(deg)'
        
        write(20,'(a)') trim(atoms(1)%name)
        
        if (natom >= 2) then
            i1 = zmat(2,1)
            bond_dist = atom_distance(2, i1) * bohr2ang
            write(20,'(a,i5,f12.6)') trim(atoms(2)%name), i1, bond_dist
        end if
        
        do i = 3, natom
            i1 = zmat(i,1)
            i2 = zmat(i,2)
            bond_dist = atom_distance(i, i1) * bohr2ang
            bond_ang = atom_angle(i, i1, i2)
            
            if (i >= 4) then
                i3 = zmat(i,3)
                dihed_ang = dihedral_angle(i, i1, i2, i3)
                write(20,'(a,i5,f12.6,i5,f10.4,i5,f10.4)') &
                    trim(atoms(i)%name), i1, bond_dist, i2, bond_ang, i3, dihed_ang
            else
                write(20,'(a,i5,f12.6,i5,f10.4)') &
                    trim(atoms(i)%name), i1, bond_dist, i2, bond_ang
            end if
        end do
        
        close(20)
    end subroutine

end program zmatrix_generator