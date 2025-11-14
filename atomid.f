c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module atomid  --  atomic properties for current atoms  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     tag       integer atom labels from input coordinates file
c     class     atom class number for each atom in the system
c     atomic    atomic number for each atom in the system
c     valnum    expected valence for each atom in the system
c     mass      atomic weight for each atom in the system
c     name      atom name for each atom in the system
c     tier      tier name (residue, motif, etc.) for each atom
c     story     descriptive type for each atom in the system
c     ntypes    number of unique atom types defined
c     atom_name array of element symbols for each atom type
c
      module atomid
      use sizes
      implicit none
      integer tag(maxatm)
      integer class(maxatm)
      integer atomic(maxatm)
      integer valnum(maxatm)
      real*8 mass(maxatm)
      character*3 name(maxatm)
      character*3 tier(maxatm)
      character*24 story(maxatm)

      
      integer, parameter :: ntypes = 2        
      character(len=2) :: atom_name(ntypes)

      ! Initialize element symbols
      data atom_name / 'H ', 'O ' /
      save
      end
