      module coordstorage
        implicit none
         real*8, allocatable, save :: coord_array(:,:)
         integer, allocatable, save :: idxandtype(:,:)
         character*3, allocatable, save :: name_array(:)
      end module coordstorage