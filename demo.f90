subroutine dft_pcg_atomic_relax(natom,coord,AL,atom_type,num_mov,force_tolerance,imov_at,iis_gupta)
    ! link this subroutine and libPCGopt.a
    ! atomMV_FFPCG and atomMV_VFF are provided in libPCGopt.a

    implicit none
    ! input output
    integer natom                ! input: number of atoms
    real*8 coord(3,natom)        ! input: atomic positions in cartesian coord
    real*8 AL(3,3), ALI(3,3)     ! input: cell lattices in Bohr unit
    integer atom_type(natom)     ! input: atomic number for each atom
    integer num_mov              ! input: the maximum number of movement
    real*8 force_tolerance       ! input: force criterion in Hartree/Bohr
    integer imov_at(3,natom)     ! input: 1,1,1, or 0,0,0, or 1 0 0 (etc), to control the movement of each atom 
    integer iis_gupta(natom)     ! input: if you don't have any metallic region, set this a 0-array simply,
                                 !        set "iis_gupta(ia) = 1" if the ia-th atom is a mettalic atom in a metallic region
    ! input output end

    real*8, parameter :: ang2bohr   =  1.8897259890       ! Angstrom to Bohr
    real*8, parameter :: hartree_ev = 27.2113961320       ! Hartree to eV
    real*8, parameter :: A_AU_1     =  0.5291772490       ! Bohr to Angstrom
    integer mov,i,j
    integer ia, ib, ic, inode, itodo_precond
    integer is_correction 
    real*8 xatom(3,natom)
    real*8 xatom_out(3,natom)
    real*8 fatom(3,natom), fatom_relax_p(3,natom), fatom_old_p(3,natom)
    real*8 tmp(3)
    real*8 reclat(3,3)
    real*8 E_tot
    real*8 force_max
    real*8 Etot_old2
    real*8 fatom_relax(3,natom)
    real*8 tolforce
    real*8 lat(3,3)
    real*8 xatom_old(3,natom)
    real*8 xatom_old2(3,natom),fatom_old2(3,natom)
    real*8 dx1,dx2,dx3,dx,dy,dz,dd,dd_ave
    real*8 fatom_old(3,natom),px(3,natom)
    
    call initial_pcg_relaxation()

    do mov = 1,num_mov
        !**********************************************************************
        ! This is a user plugin subroutine for DFT calculation
        ! input:
        !   atom_type(natom): atomic number for all atoms
        !   xatom: fractional coordinates
        !   AL: lattice 3x3, unit: Bohr
        call calculate_energy_force(atom_type, xatom, AL, E_tot, fatom)
        ! output:
        !   E_tot: total energy, unit: Hartree
        !   fatom: force, unit: Hartree/Bohr
        !**********************************************************************

        ! force criterion
        force_max = -100.d0
        do i = 1,natom
            if(dabs(fatom(1,i)*imov_at(1,i)).gt.force_max) force_max = dabs(fatom(1,i)*imov_at(1,i))
            if(dabs(fatom(2,i)*imov_at(2,i)).gt.force_max) force_max = dabs(fatom(2,i)*imov_at(2,i))
            if(dabs(fatom(3,i)*imov_at(3,i)).gt.force_max) force_max = dabs(fatom(3,i)*imov_at(3,i))
        enddo
        if(force_max .lt. force_tolerance) then
            write(*,*) "end, force tolforce"
            goto 2333  ! goto exit point
        endif
        fatom_relax(:,:) = fatom_relax(:,:)*imov_at(:,:)
        is_correction = 1
        ! PCG here, atomMV_FFPCG and atomMV_VFF are provided in the libPCGopt.a, directly link them
        call atomMV_FFPCG(is_correction,natom,E_tot,lat,reclat,xatom,xatom_old,fatom_relax,fatom_old,px,tolforce,inode,imov_at,itodo_precond,fatom_relax_p,fatom_old_p)
        if(itodo_precond.eq.1) then  ! calculate the preconditioning, recall it
            call atomMV_VFF(natom,atom_type,fatom,lat,reclat,xatom,xatom_out,mov,imov_at)
            call get_fatom_relax_p()

            call atomMV_FFPCG(is_correction,natom,E_tot,lat,reclat,1,xatom,xatom_old,fatom_relax,fatom_old,px,tolforce,inode, imov_at,itodo_precond,fatom_relax_p,fatom_old_p)
        endif    ! itodo_precond
        ! end PCG
        is_correction = 0
    enddo ! mov

2333    continue     ! exit point
    coord = matmul(lat,xatom)

    contains
            subroutine initial_pcg_relaxation()
                ! xatom_cartesian = matmul(lat,xatom_fraction)       
                ! xatom_fraction = matmul(transpose(reclat),xatom_cartesian)
                ! AL = lat * ang2bohr
                ! ALI = inverse(transpos(AL))
                ! reclat = inverse(transpos(lat))

                ! tolforce in hatree/bohr, AL in \AA
                tolforce = force_tolerance / hartree_ev / ang2bohr
                lat(1:3,1:3) = al(1:3,1:3) / ang2bohr
             
                ! AL->ALI, lat->reclat
                tmp(1:3) = 1
                ALI=transpose(AL)
                reclat=transpose(lat)
                call gaussj(ALI, 3, 3, tmp, 1, 1)
                call gaussj(reclat, 3, 3, tmp, 1, 1)
                xatom = matmul(transpose(reclat),coord)
            end subroutine initial_pcg_relaxation
                        
            subroutine get_fatom_relax_p()
                do i=1,natom
                    dx1=xatom_out(1,i)-xatom(1,i)
                    dx2=xatom_out(2,i)-xatom(2,i)
                    dx3=xatom_out(3,i)-xatom(3,i)
             
                    if(dx1.gt.0.5) dx1=dx1-1
                    if(dx1.lt.-0.5) dx1=dx1+1
                    if(dx2.gt.0.5) dx2=dx2-1
                    if(dx2.lt.-0.5) dx2=dx2+1
                    if(dx3.gt.0.5) dx3=dx3-1
                    if(dx3.lt.-0.5) dx3=dx3+1
             
                    dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
                    dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
                    dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
             
                    fatom_relax_p(1,i)=dx*imov_at(1,i)*A_AU_1
                    fatom_relax_p(2,i)=dy*imov_at(2,i)*A_AU_1
                    fatom_relax_p(3,i)=dz*imov_at(3,i)*A_AU_1
                enddo
                ! use symmetry operator here, if necessary
            end subroutine get_fatom_relax_p


            
end subroutine
