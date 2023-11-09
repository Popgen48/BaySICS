 module treeator


    implicit none

	contains



		

   subroutine tree_maker(n,dn,gr,en,gentime,samplinfo,ages,Ne,growth,events,taxa,nodes,coal_t,&
              &NodRanges,maxL,NM,MS,TM,Mig_labels,MigMat)

      use routines
      use irandom
      use luxury

      implicit none

      !External variables
      integer, intent(in):: n                             !Overall sample size
      integer, intent(in):: dn                            !# of demes in the sample
      integer, intent(in):: gr                            !Number of subsamples (by deme and age)
      integer, intent(in):: en                            !Number of events
      real(8), intent(in):: gentime                       !Generation time
      real(8), dimension(3,gr), intent(in):: samplinfo    !Array with the 1) sizes, 2) deme identities and 3) ages of temporal subsamples
      real(8), dimension(n), intent(in):: ages            !Ages of each taxa (individual) in the sample (ordered)
      real(8), dimension(dn), intent(in):: Ne             !Vector of population sizes of each population (deme)
      real(8), dimension(dn), intent(in):: growth         !Vector of growth rates for each block of coalescent simulations
      real(8), dimension(en,8), intent(in):: events       !Array with the info of ancient demographic events:
                                                        !  1:= time;
                                                        !  2:= Ne;
                                                        !  3:= new growth rate;
                                                        !  4:= 1st deme involved;
                                                        !  5:= 2nd deme involved;
                                                        !  6:= % of entring lineages from 1st pop (deme);
                                                        !  7:= % of entring lineages from 2nd pop (deme);
                                                        !  8:= label of the block
      integer, dimension(n),intent(out):: taxa                      !Taxa store the permutation order of the taxa
      integer, dimension(n-1), intent(out):: nodes                  !Nodes store the permutation order of the coalescent times
      real(8), dimension(n-1), intent(out):: coal_t                 !Array storing the coalescent times
      integer, dimension(2,n-1), intent(out):: NodRanges            !Array with the range of each of the final nodes of the tree
      real(8), intent(in):: maxL                                    !Maximum limit:= its the ceiling for coalescent times
      integer, intent(in):: NM                                      !Number of migration matrix
      integer, dimension(NM), intent(in):: MS                       !Array with the sizes of those matrix
      real, dimension(NM,2), intent(in):: TM                        !Timespans for the validity of the corresponding matrix
      integer, dimension(NM,dn+en), intent(in):: Mig_labels         !We only need the columns/rows with values (this are their labels)
      real(8), dimension(NM,en+dn,en+dn), intent(in):: MigMat       !Migration matrix
      !Internal variables
      !Integers and counters
      integer:: h, i, j, l                                          !Cycles counters
      integer:: k                                                   !Number of lineages
      integer:: x, a, b, c, d, e, f, z, zbis, nod                   !Counters
      !Arrays
      integer, dimension(n):: taxa2, taxa3, taxa4                   !Taxa store the permutation order of the taxa
      integer, dimension(n-1):: nodes2, nodes4, nodes5, nodes5B     !Nodes store the permutation order of the coalescent times
      real(8), dimension(n-1):: nodes3                              !Nodes store the permutation order of the coalescent times
      integer, dimension(n):: taxaC                                 !Taxa store the permutation order of the taxa
      integer, dimension(n-1):: nodesC                              !Nodes store the permutation order of the coalescent times
      integer, dimension(dn+en):: ni, gri                           !ni:= stores the sample size of all the blocks/populations; gri:= stores the grs of all the blocks/populations
      integer, dimension(dn+en,2):: demes_ranges                    !Array with the upper and lower boundaries of each block population (cell indicate positions in the large taxa array)
      real(8), dimension(dn+en,3,gr):: Msamplinfo                   !Master samplinfo stores the samplinfos of all blocks/populations
      integer, dimension(dn+en,n):: Mtaxa                           !Master taxa stores the taxas of all blocks/populations
      integer, dimension(dn+en,n):: MtaxaB                          !Master taxa stores the taxas of all blocks/populations with internal local numeration (serves for intrabock instructions)
      real(8), dimension(dn+en,n):: Mages                           !Master ages contains the ages of the blocks/populations
      integer, dimension(dn+en,n-1):: Mnodes                        !Master nodes contains the nodes of all the blocks/populations with continuous overall numeration (serves for interblock instructions)
      integer, dimension(dn+en,n-1):: MnodesB                       !Master nodes contains the nodes of all the blocks/populations with internal local numeration (serves for intrabock instructions)
      real(8), dimension(dn+en,n-1):: Mcoal_t                       !Master coal_t contains all the coal_t of all the blocks/populations
      real(8), allocatable, dimension(:,:):: samplinfoi             !A transient samplinfo for coalescing each block/population
      real(8), allocatable, dimension(:):: agesi                    !A transient ages for coalescing each block/population
      integer, allocatable, dimension(:):: taxai                    !A transient taxa for coalescing each block/population
      real(8), allocatable, dimension(:):: taxaii                   !Transient taxa with type real (for variate uses)
      integer, allocatable, dimension(:):: nodesi                   !A transient nodes for coalescing each block/population
      real(8), allocatable, dimension(:):: nodes2i                  !A transient nodes for coalescing each block/population
      real(8), dimension(3,(n-1)*(en+1)):: Nodages                  !Nodes information: 1) Age (a coal_t value); 2-3) the range it covers,...
                                                                    !...IMPORTANT! The position has 'n' substracted since their numeration continues that of taxa
      real(8), allocatable, dimension(:):: coal_ti                  !A transient coal_t for coalescing each block/population
      logical, dimension(2,en):: favor                              !An array that associatess one logical variable to each event: true if they already have the entring lineages chosen from another block instruction
      real, allocatable, dimension(:):: rvec                        !Random numbers vector
      real(8), dimension(n-1):: y                                   !Array used for sorting
      !Reals
      real(8):: bound
      real(8), allocatable, dimension(:):: x1, y1
      !Migration
      logical:: Access, shared                                      !Access: True if an entring lineage is suitable for accessing a block; shared: it lets a block to know if last block shared an entring block
      logical:: last                                                !Last:= Switch being .true. if the coalescing block is the last one (so truncation of times surpasing maxL is deactivated)
      logical:: switch, switch4
      integer:: receiver
      integer:: dummy
      real(8), allocatable, dimension(:,:):: emigrants            !1st column:= Numbers of the node/taxa that emigrate;
                                                                  !2nd column:= time;
                                                                  !3rd column:= donor;
                                                                  !4th column:= receiver population/block;
                                                                  !5th column:= 0 if it is a taxon and 1 if it is a node
      real(8), allocatable, dimension(:,:):: V
      real(8):: time, etime, Ptime
      logical:: mig, New_mig                                        !Switch indicating if migration is present; New_mig: if a new migration happened in a migrant lineage
      integer:: r1, r2, m1, m2, lineage                             !Storing the ranges of receiver and migrant lineages
      integer:: Mr1, Mr2, Mm1, Mm2                                  !Storing the ranges of receiver and migrant lineages in the 'M' arrays (Mtaxa,Mnodes,etc.)
      real(8), dimension(2):: node1, node2, node3                   !1st is the new name (id) of the node; 2nd is the age
      real(8):: nodeM1, nodeM3                                      !Store special nodes for migration
      integer:: pos
      integer:: A1, A2, TN1, TN2, Em, Lm, go_to, nbis
      logical:: older                                               !True if a migrant lineage coalesces outside the receiver block (and goes back to an ancestor)
      !*****************************************************************************************************************************
      !   S   T   A   R   T
      !*****************************************************************************************************************************
      !(0) Separating samplinfo, taxa, ages, nodes and coal_t in mini arrays (one for each pop/block)
      !-----------------------------------------------------------------------------------------------------------------------------
      if (ANY(MigMat/=0)) then
         mig=.true.
      else
         mig=.false.
      endif
      !0.1) Samplinfo, gr & n
      Msamplinfo=0
      ni=0
      gri=0
      x=1    ! x counts the deme number
      a=1
      do i = 1, gr, 1
         if (samplinfo(2,i)/=x) then
            x=x+1
            a=1
         endif
         ni(x) = ni(x) + NINT(samplinfo(1,i))
         gri(x) = gri(x) + 1
         Msamplinfo(x,1,a) = samplinfo(1,i)
         Msamplinfo(x,2,a) = samplinfo(2,i)
         Msamplinfo(x,3,a) = samplinfo(3,i)
         a=a+1
      enddo
      !0.2 Ages
      Mages=0
      do h = 1, dn, 1
         x=1
         do i = 1, gri(h), 1
            do l = 1, nint(Msamplinfo(h,1,i)), 1
               Mages(h,x) = Msamplinfo(h,3,i)
               x = x + 1
            enddo
         enddo
      enddo
      !0.3 Taxa, coal_t and nodes
      Mtaxa=0
      Mnodes=0
      Mcoal_t=0
      x=1
      do h = 1, dn, 1
         do i = 1, ni(h)-1, 1   !Only initializing mini arrays with ordered numbers
            Mtaxa(h,i)=x
            Mnodes(h,i)=x - (h-1) + n !Substraction of (h-1) adjust the count
            x=x+1
         enddo
         Mtaxa(h,ni(h))=x
         x=x+1
      enddo
      !---
      MnodesB=MAX(Mnodes-n,0)
      !---
      Nodages=0
      !---
      favor=.false.
      !-----------------------
      !X.1 Making demes_ranges
      x=0
      do i = 1, dn, 1
         demes_ranges(i,1)=x+1
         demes_ranges(i,2)=x+ni(i)
         x=x+ni(i)
      enddo
      !-----------------------------------------------------------------------------------------------------------------------------
      !(1) Coalescing each population (deme)
      !-----------------------------------------------------------------------------------------------------------------------------
      do h = 1, dn, 1
         !1.1) Making mini arrays needed for coalescing the block/population taking the info from their Masters
         allocate(samplinfoi(3,gri(h)), taxai(ni(h)), agesi(ni(h)), nodesi(ni(h)-1), coal_ti(ni(h)-1))
         !----------
         do i=1,3,1
            do j=1,gri(h),1
               samplinfoi(i,j)=Msamplinfo(h,i,j)
            enddo
         enddo
         do i=1,ni(h)-1,1
            agesi(i)=Mages(h,i)
            taxai(i)=i
            nodesi(i)=i
         enddo
         taxai(ni(h))=ni(h)
         agesi(ni(h))=Mages(h,ni(h))
         coal_ti=0
         !Now we turn on/off the switch for indicating if the coalescing block is the last one
         last=.false.
         if (dn==1) then
            last=.true.
         endif

         !1.2) Coalesce bastards! Coalesce! ha ha ha!
         !------------------------------------------------------------------------------------------------
         call AncTree(gentime,Ne(h),growth(h),ni(h),gri(h),samplinfoi,taxai,nodesi,coal_ti,maxL,last)
         !------------------------------------------------------------------------------------------------
         if (mig) then
            if (h==1) then
               allocate(emigrants(2*n,5))
               emigrants=0
            endif
            m1=COUNT(emigrants(:,1)/=0)
            call Migration(h,n,ni(h),gentime,agesi,taxai,nodesi,coal_ti,NM,dn+en,Mig_labels,MS,TM,MigMat,emigrants)
            m2=COUNT(emigrants(:,1)/=0)
         endif
         !1.3) Retrieveng the obtained information into the main arrays
         do i=1,ni(h)-1,1
            MtaxaB(h,i)=taxai(i)
            taxai(i)=Mtaxa(h,taxai(i))
            MnodesB(h,i)=nodesi(i)         !MnodesB should be filled here since info in nodesi is gonna be changed
            nodesi(i)=Mnodes(h,nodesi(i))  !Since the info in taxa and nodes is a permutation of their own numbers it is needed to use a third
         enddo                             !...array to move the info, however taxai and nodesi can be used without loss of information
         MtaxaB(h,ni(h))=taxai(ni(h))
         taxai(ni(h))=Mtaxa(h,taxai(ni(h)))
         !---
         if (mig) then
            if ((m2-m1)>0) then
               do i=m1+1,m2,1
                  if (emigrants(i,5)==0) then  !If it is a taxon
                     emigrants(i,1)=Mtaxa(h,NINT(emigrants(i,1)))
                  elseif (emigrants(i,5)==1) then !If it is a node
                     emigrants(i,1)=Mnodes(h,NINT(emigrants(i,1)))
                  endif
               enddo
            endif
         endif
         !---
         do i=1,ni(h)-1,1
            Mtaxa(h,i)=taxai(i)
            Mnodes(h,i)=nodesi(i)
            Mcoal_t(h,i)=coal_ti(i)
         enddo
         Mtaxa(h,ni(h))=taxai(ni(h))
         !...also filling the corresponding NodAges
         do i=1,ni(h)-1,1
            Nodages(1,Mnodes(h,i)-n)=coal_ti(MnodesB(h,i)) !IMPORTANT!:Always add 'n' to the position in Nodages to match the number of node
         enddo 

         deallocate(samplinfoi,taxai,agesi,nodesi,coal_ti)
      enddo
      !-----------------------------------------------------------------------------------------------------------------------------
      !(2) Coalescing blocks
      !-----------------------------------------------------------------------------------------------------------------------------
          nod=2*n-dn
          !----
          do i=dn+1,en+dn,1
             ni(i)=0
          enddo
          do h = 1, en, 1
             !ni(h+dn)=0
             do i = 4, 5, 1 !This for doing all next stuff for both populations/blocks that are sources of our block
                if (   (.not.favor(i-3,h))  .and.  ( (i==4) .or. (events(h,4)/=events(h,5)) )   ) then  !We enter only if favor was not received (other block already choose our lineages); and if it is the first one or if is the second one but is different
                    a = NINT(events(h,i)) !The number of entering (source) population/block
                    bound = events(h,1) !Time to the beginning of the block

                  !2.1) Making Mages and Mtaxa for the entring lineages of the entring population/block
                    allocate(rvec(1))
                    !------------------------------------
                    !Here it goes an instruction for locating the other block that shares the lineages of deme "a" so that we fill its arrays also
                    e=0
                    if (events(h,i+2)<1.0) then
                       if (h<en) then
                          do j=h+1,en,1
                             if (events(j,4)==a) then
                                e=j
                                favor(1,e)=.true.        ! e gets the position of the alternative block
                             endif
                             if (events(j,5)==a) then
                                e=j
                                favor(2,e)=.true.
                             endif
                          enddo
                       endif
                       e=e+dn !So that we can use e as the absolute position in Mtaxa and Mages
                    endif
                    !------------------------------------
                    d=COUNT(Mtaxa(h+dn,:)/=0.0)   !The counts the already filled cells (when favor has been done)
                    !------------------------------------
                    z=0 !z counts only the accepted lineages
                    zbis=0
                    c=0 !c indicates the position of the last node being older than boundary
                    do j = 1, ni(a), 1
                       Access=.false.                                       !This instruction is for deciding if the lineage is suitable for accessing the block, which happens if:
                       if (j==ni(a)) then                                               !1)The position reached is the last one, then the group after the last suitable lineage automatically is suitable
                          Access = .true.
                       elseif ( (j<ni(a)) .and. (Nodages(1,Mnodes(a,j)-n)>bound) ) then !2)If the age of the node is larger than boundary
                          Access = .true.
                       endif
                       if (Access) then
                          call ranlux(rvec,1)
                          if (rvec(1)<events(h,i+2)) then
                             !If the lineage succeded in entering the block---------------------------------------------------------------------------------------------------
                             z=z+1
                             if ( (j==1) .or. (c==j-1) ) then              !The lineage is a taxa (i.e. if the immediate previous node is also older than boundary, or is the first node)
                                   Mtaxa(h+dn,z+d)=Mtaxa(a,j)                          !So we record in Mtaxa of the block then umber of node of the entring block/population
                                if (Mtaxa(a,j)<=n) then        !In spite of this lineage is a taxa since the previous node it could be actually a node two blocks before so we have to be shure we are checking the proper age (Nodages or ages)
                                   Mages(h+dn,z+d)=MAX(bound,ages(Mtaxa(a,j)))         !We add "d" in this instructions for (re)starting to fill Mtaxa and Mages from the last position they wre filled (if a previous cycle was already run)
                                else
                                   Mages(h+dn,z+d)=MAX(bound,Nodages(1,Mtaxa(a,j)-n))
                                endif
                             else                                          !The lineage is a node
                                b=0
                                do l=c+1,j-1
                                   b=MAX(b,Mnodes(a,l))
                                enddo
                                Mtaxa(h+dn,z+d)=b       !We name the taxa cell of the block with the name (number) of the node dominating the entring group
                                Mages(h+dn,z+d)=bound
                                !---also recording the range of the node (it's gonna be used later)
                                if (a<=dn) then !We only get range of a node if belongs to a population (internal blocks ranges are found in section 3)
                                   Nodages(2,b-n)=c+1 + demes_ranges(a,1) - 1   !We want the position in taxa so we add the lowest boundary - 1
                                   Nodages(3,b-n)=j   + demes_ranges(a,1) - 1
                                endif
                             endif
                          else
                             !If it failed entering the block in the draw, necesarily there is an alternative block-----------------------------------------------------------
                             zbis=zbis+1
                             if ( (j==1) .or. (c==j-1) ) then              !The lineage is a taxa (i.e. if the immediate previous node is also older than boundary, or is the first node)
                                Mtaxa(e,zbis)=Mtaxa(a,j)                   !So we record in Mtaxa of the block the number of node of the entring block/population
                                if (Mtaxa(a,j)<=n) then              ! BUT it not necesarily is a taxa cuz could be a node inherited from a previous block (taken as taxa in the entring block)
                                   Mages(e,zbis)=MAX(bound,ages(Mtaxa(a,j)))   !We add "d" in this instructions for (re)starting to fill Mtaxa and Mages from the last position they were filled (if a previous cycle was run)
                                elseif (Mtaxa(a,j)>n) then
                                   Mages(e,zbis)=bound
                                endif
                             else                                          !The lineage is a node
                                b=0
                                do l=c+1,j-1
                                   b=MAX(b,Mnodes(a,l))
                                enddo
                                Mtaxa(e,zbis)=b       !We name the taxa cell of the block with the name (number) of the node dominating the entring group
                                Mages(e,zbis)=bound
                                !---also recording the range of the node (it's gonna be used later)
                                if (a<=dn) then !We only get range of a node if belongs to a population (internal blocks ranges are found in section 3)
                                   Nodages(2,b-n)=c+1 + demes_ranges(a,1) - 1   !We want the position in taxa so we add the lowest boundary - 1
                                   Nodages(3,b-n)=j   + demes_ranges(a,1) - 1
                                endif
                             endif
                          endif !---------------------------------------------------------------------------------------------------------------------------------------------
                          c=j !Recording the position of the last old node
                       endif
                    enddo
                    !
                    deallocate(rvec)
                    ni(h+dn)=ni(h+dn)+z !z is the number of lineages we filled in this lap
                    if (events(h,i+2)<1.0) then !The stuff regarding the alternate block ("e" and "zbis") only happens if there is alternate block (i.e. % of entring lineages is < 100%)
                       ni(e)=ni(e)+zbis
                    endif
                endif
             enddo

            !2.2) Now we're gonna make samplinfo
             if (ni(h+dn)>0) then
               !Instruction for sorting ages and taxa and
               allocate(taxai(ni(h+dn)),taxaii(ni(h+dn)),agesi(ni(h+dn)),nodesi(ni(h+dn)-1),coal_ti(ni(h+dn)-1))
               do i=1,ni(h+dn),1
                  agesi(i)=Mages(h+dn,i)
                  taxaii(i)=dble(Mtaxa(h+dn,i))
               enddo
               !--
               call SSORT(agesi, taxaii, ni(h+dn), 2)
               !--don't forget returning the sorted ages to Mages and Mtaxa
               do i=1,ni(h+dn),1
                  Mages(h+dn,i)=agesi(i)
                  Mtaxa(h+dn,i)=int(taxaii(i))
               enddo
               !Making samplinfo from ages
               e=0
               do i = 1, ni(h+dn), 1
                  access=.false.
                  if (i==1) then
                     access=.true.
                     e=e+1
                  elseif (  (i>1) .and. ( agesi(i) /= agesi(i-1) )  ) then
                     access=.true.
                     e=e+1
                  endif
                  Msamplinfo(h+dn,1,e)=Msamplinfo(h+dn,1,e)+1
                  if (access) then
                    Msamplinfo(h+dn,2,e)=h+dn
                    Msamplinfo(h+dn,3,e)=agesi(i)
                  endif
               enddo
               gri(h+dn)=e
               !--
               allocate(samplinfoi(3,e))
               do i=1,3,1
                  do j=1,e,1
                     samplinfoi(i,j)=Msamplinfo(h+dn,i,j)
                  enddo
               enddo
              !2.3) Preparing arrays for coalescing
               do i=1,ni(h+dn)-1
                  taxai(i)=i
                  nodesi(i)=i
               enddo
               taxai(ni(h+dn))=ni(h+dn)
               !--------------------
               last=.false.
               if (h==en) then
                  last=.true.
               endif
          !--------------------------------------------------------------------------------------------------------------------
          call AncTree(gentime,events(h,2),events(h,3),ni(h+dn),gri(h+dn),samplinfoi,taxai,nodesi,coal_ti,maxL,last)
          !--------------------------------------------------------------------------------------------------------------------
               !2.35 Migration!
               if (.NOT.last) then
                  if (mig) then
                     m1=COUNT(emigrants(:,1)/=0)
                     call Migration(h+dn,n,ni(h+dn),gentime,agesi,taxai,nodesi,coal_ti,NM,dn+en,Mig_labels,MS,TM,MigMat,emigrants)
                     m2=COUNT(emigrants(:,1)/=0)
                  endif
               endif
              !2.4) Then retrieving the obtained taxa, nodes, ages, and coal_t to their masters
               do i=1,ni(h+dn)-1,1
                  MtaxaB(h+dn,i)=taxai(i)
                  taxai(i)=Mtaxa(h+dn,taxai(i))
                  MnodesB(h+dn,i)=nodesi(i)         !Mnodes/MnodesB can be filled here
                  Mnodes(h+dn,i)=nodesi(i)+nod      !Since the info in taxa and nodes is a permutation of their own numbers it is needed to use a third
               enddo                                !...array to move the info, however taxai and nodesi can be used without loss of information
               MtaxaB(h+dn,ni(h+dn))=taxai(ni(h+dn))
               taxai(ni(h+dn))=Mtaxa(h+dn,taxai(ni(h+dn)))
               !---
               if (.NOT.last) then
                  if (mig) then
                     if ((m2-m1)>0) then
                        do i=m1+1,m2,1
                           if (emigrants(i,5)==0) then  !If it is a taxon
                              emigrants(i,1)=Mtaxa(h+dn,taxai(i))
                           elseif (emigrants(i,5)==1) then !If it is a node
                              emigrants(i,1)=Mnodes(h+dn,nodesi(i))
                           endif
                        enddo
                     endif
                  endif
               endif
               !---
               do i=1,ni(h+dn)-1,1
                  Mtaxa(h+dn,i)=taxai(i)
                  Mcoal_t(h+dn,i)=coal_ti(i)
               enddo
               Mtaxa(h+dn,ni(h+dn))=taxai(ni(h+dn))
               !...also filling the corresponding NodAges
               do i=1,ni(h+dn)-1,1
                  Nodages(1,Mnodes(h+dn,i)-n)=coal_ti(MnodesB(h+dn,i)) !IMPORTANT!:Always add 'n' to the position in Nodages to match the number of node
               enddo
               !------------------------------------------------
               deallocate(taxai,taxaii,agesi,nodesi,coal_ti,samplinfoi)
             endif
             nod=nod+MAX(0,ni(h+dn)-1)
          enddo
         !---------------------------------------------------------------------------------
         !2.5) Updating the original Nodes and Taxa arrays, only with the coalescents of the populations (not the internal blocks)
         x=0
         do i=1,dn
             do j=1,ni(i)-1
                 x=x+1
                 taxa(x)=Mtaxa(i,j)
                 nodes(x)=Mnodes(i,j)
             enddo
             x=x+1
             taxa(x)=Mtaxa(i,ni(i))
         enddo

        !-----------------------------------------------------------------------------------------------------------------------------
        !(3) THE MONSTER: The most terrible instruction ever: assembling all the trees from the populations/blocks into one single tree (aaaaahhhhhhh!!!)
        !-----------------------------------------------------------------------------------------------------------------------------
          shared=.false.
          do h = 1, en, 1 !We're gonna do this for all the blocks
           !3.1) Allocate arrays and get the info on them
              c=0
              do i=1,ni(h+dn),1 !This instruction is only for getting the size of the array (the range of the block)
                 if (Mtaxa(h+dn,i)<=n) then !its a taxa so only accounts for 1
                    c=c+1
                 else                        !its a node so it accounts for its range
                    c=c+NINT(Nodages(3,Mtaxa(h+dn,i)-n)-Nodages(2,Mtaxa(h+dn,i)-n))+1
                 endif
              enddo
              !
              if (c==0) then !Theres nothing to do
                  shared=.false.
              else           !Then everything has to be done
                 allocate(taxai(c),nodesi(c-1)) !The arrays are gonna work in the range of the block
                 !Mezzaninne) Get the starting point in taxa (it is the lowest boundary of the range of first entring population/block)
                 if (ni(NINT(events(h,4)))>0) then
                    a = demes_ranges(NINT(events(h,4)),1)
                 else
                    a = demes_ranges(NINT(events(h,5)),1)
                 endif
                 b = a+c !..& starting point of the next block
                 x=0
               !3.2) Get the corresponding taxa and nodes into taxai and nodesi
                 do i=1,ni(h+dn),1 !Have to walk all the taxa of the block & move them to taxai, also have to move nodes to nodesi
                    if (Mtaxa(h+dn,i)<=n) then         !It is an individual taxa
                       x=x+1
                       taxai(x)=Mtaxa(h+dn,i)
                    else                               !It is a node (group of taxa)
                       f=x+1 !Update first value of the range
                       do j = NINT(Nodages(2,Mtaxa(h+dn,i)-n)), NINT(Nodages(3,Mtaxa(h+dn,i)-n)-1), 1
                          x=x+1
                          if (shared) then
                             taxai(x)=taxaC(j)
                             nodesi(x)=nodesC(j)
                          else
                             taxai(x)=taxa(j)
                             nodesi(x)=nodes(j)
                          endif
                       enddo
                       x=x+1
                       if (shared) then
                          taxai(x)=taxaC(NINT(Nodages(3,Mtaxa(h+dn,i)-n)))
                       else
                          taxai(x)=taxa(NINT(Nodages(3,Mtaxa(h+dn,i)-n)))
                       endif
                       !Now we update the 2nd value of the range of this node
                       Nodages(2,Mtaxa(h+dn,i)-n) = f + a-1
                       Nodages(3,Mtaxa(h+dn,i)-n) = x + a-1
                    endif
                    if (i<ni(h+dn)) then
                       nodesi(x)=Mnodes(h+dn,i) !This instruction makes nodesi getting the block nodes, whose are acually filling the nodes not occupied by entring groups
                    endif
                 enddo

               !3.3) Updating the ranges of the new nodes
                 if (c>1) then
                    z=1   !z walks Mnodes
                    i=1   !i walks nodesi
                    do
                       if (nodesi(i)==Mnodes(h+dn,z)) then
                          z=z+1
                          x=0
                          do
                             if (   ( Nodages(1,nodesi(i+x)-n) > Nodages(1,nodesi(i)-n) )  .or.  (i+x==c-1)   ) exit !i.e. if we reach a node older than Mnodes(,z) or the last one exit
                             x=x+1
                          enddo
                          if ( Nodages(1,nodesi(i+x)-n) > Nodages(1,nodesi(i)-n) ) then !If we reached the end the final node is the next one (our range is ok IN NODES)
                             e=i+x-1
                          else
                             e=i+x       !If we reach the node that is also outside the range IN NODES (range is for taxa)
                          endif
                          x=0
                          do
                             if (   ( Nodages(1,nodesi(i-x)-n) > Nodages(1,nodesi(i)-n) )  .or.  (i-x==1)   ) exit !i.e. if we reach a node older than Mnodes(,z) or the last one exit
                             x=x+1
                          enddo
                          if ( Nodages(1,nodesi(i-x)-n) > Nodages(1,nodesi(i)-n) ) then !If we not reached the end
                             d=i-x+1
                          else
                             d=i-x       !If we reach the node that is also outside the range IN NODES (range is for taxa)
                          endif
                          Nodages(2,nodesi(i)-n)=d    + a-1
                          Nodages(3,nodesi(i)-n)=e+1  + a-1   !This is cuz the range is for taxa so if we reach node x, the taxa included is one more and count starts from 'a' (i.e. starting point of the block in taxa)
                       endif
                       if (i==c-1) exit
                       i=i+1
                    enddo
                 endif

                !3.3.2 Relocating the taxa/nodes that doesn't enter the block if second entring block is shared
                 if (h<en) then   !If the 2nd entring block is shared       !NOTE: Only 2nd entring block is checked since first one is adjusted by the previous block
                    if (events(h,5)==events(h+1,4)) then  !if theres a shared pop/block it corresponds to the 2nd entring pop/block of first block and 1st entring pop/block of next block
                       shared=.true.
                       taxaC=taxa
                       nodesC=nodes
                    else
                       shared=.false.
                    endif
                 endif

                !3.4) Retrieving taxai/nodesi to taxa and nodes
                 !So..
                 do i=1,c-1,1
                    taxa(a+i-1)=taxai(i)
                    nodes(a+i-1)=nodesi(i)
                 enddo
                 taxa(a+c-1)=taxai(c)

                !3.4.1 Update demes_ranges (this adjustment allows the next block starting at the right point (we temporarily shift the range of the entring shared block))
                 !3.4.1.1 updating demes_ranges of the 2nd entring block if shared
                 if (shared) then
                    demes_ranges(NINT(MIN(events(h,4),events(h,5))),1)=b  !Updating
                 endif
                 !3.3.3.2 Updating demes_ranges of our current block
                 demes_ranges(h+dn,1)=a
                 demes_ranges(h+dn,2)=b-1

                 deallocate(taxai,nodesi)
             endif
          enddo
          !Monster done!
      !-----------------------------------------------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------------------------------------------
      !(5) But wait, there's more: M I G R A T I O N
      !-----------------------------------------------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------------------------------------------
      !5.0 Sort the migrations from present to past
      if (mig) then
         E = COUNT(emigrants(:,1)/=0)
         if (E>0) then
            !---- This instruction is only to sort emigrants by age
            allocate(V(E,5))
            do i=1,E,1
               V(i,1)=emigrants(i,2)
               V(i,2)=i
            enddo
            CALL SSORT(V(:,1),V(:,2), E, 2)
            do i=1,E,1
               m1=NINT(V(i,2))  !m1 is only being recycled
               do j=1,5,1
                  V(i,j)=emigrants(m1,j)
               enddo
            enddo
            do i=1,E,1
               do j=1,5,1
                  emigrants(i,j)=V(i,j)
               enddo
            enddo
            !---
            !Cycle of migration events
            i=0
            do
                i=i+1
                !-------------------------------------------------------------
                do j=1,n-1,1  !nodes3 will be used several times
                   nodes3(j) = Nodages(1,nodes(j)-n)
                enddo
                nodes4=nodes
                taxa3=taxa
                !-------------------------------------------------------------
                Switch4=.true.
                if (emigrants(i,1)>n) then !If the migrant is a node
                   if ( Nodages(1,NINT(emigrants(i,1)-n)) == MAXVAL(nodes3) ) then !This means the migrant is already the MRCA
                      Switch4=.false.
                      i=E              !Since migration events are sorted by age, at the time we get one migrant being the MRCA, we are done with migration
                   endif
                endif
                !******************************************************************************************************************
                IF (Switch4) THEN
                    !******************************************************************************************************************
                    !5.1   G L O B A L   D O N O R   A D J U S T M E N T S
                    !------------------------------------------------------------------------------------------------------------------
                    b=NINT(emigrants(i,3))       !The donor population/block
                    lineage=NINT(emigrants(i,1)) !This is the lineage being erased cuz is the migrant or its parental node
                    go_to=lineage
                    !-------------------------Next instruction is only to find the parental node of the migrant
                    do j=1,n-1,1
                       taxa4(j)=Mtaxa(b,j)
                       nodes5B(j)=MnodesB(b,j)
                       nodes5(j)=Mnodes(b,j)
                    enddo
                    taxa4(n)=Mtaxa(b,n)
                    nbis=ni(b)
                    do j=1,en,1
                       ALLOCATE(nodes2i(nbis-1))   !nbis instead ni(k) because if more than 1 block (event) will be analysed nodes2i will store the block...
                       nodes2i=nodes5B             !...without being adjusted and ni(k) will store the size after the adjustment
                       if ((ANY(taxa4==lineage)).or.(ANY(nodes5==lineage)).or.(ANY(taxa4==go_to))) then  !This avoids entering blocks if the migrant & parental already coalesced
                           if (nbis>1) then   !If there is only one lineage in the block, theres no need to look for anything
                              if (  (ANY(taxa4==lineage)) .or.  (ANY(nodes5==lineage))  ) then
                                 if (ANY(taxa4==lineage)) then !i.e. if we got the migrant lineage
                                    pos = iMINMAXLOC1(1,n,taxa4,lineage,lineage,1,n)
                                    CALL GET_GOTO(.false.,nbis-1,nodes2i,pos,go_to,m1,m2)   !So go_to is the parental node
                                 elseif (ANY(nodes5==lineage)) then !i.e. if we got the parental lineage of the migrant lineage
                                    pos = iMINMAXLOC1(1,n-1,nodes5,lineage,lineage,1,n-1)
                                    CALL GET_GOTO(.true.,ni(b)-1,nodes2i,pos,go_to,m1,m2)    !So go_to is the parental node
                                 endif
                                 !----
                                 go_to=nodes5(go_to)
                                 !-------------------------Next instruction finds the sister lineage of the migrant
                                 pos = iMINMAXLOC1(1,n-1,nodes5,go_to,go_to,1,n-1)
                                 CALL Get_downstream_lin(nbis-1,pos,nodes2i,A1,TN1,A2,TN2)
                                 if (TN1==0) then
                                    A1=taxa4(A1)
                                 else
                                    A1=nodes5(A1)
                                 endif
                                 if (TN2==0) then
                                    A2=taxa4(A2)
                                 else
                                    A2=nodes5(A2)
                                 endif
                                 if (lineage==A1) then  !If the left incidence is actually the current emigrant
                                    Em=A2          !...then the other incidence is the one that inherit the name of the dissapeared node
                                    Lm=TN2
                                 elseif (lineage==A2) then !If the left incidence is actually the current emigrant
                                    Em=A1             !...then the other incidence is the one that inherit the name of the dissapeared node
                                    Lm=TN1
                                 endif
                              endif
                           else
                              go_to=lineage
                           endif
                           !------------------------------------------------------------------
                           if ((events(j,4)==b).or.(events(j,5)==b)) then !We found a block above the donor block
                              k=NINT(events(j,8)) !This is the actual id of the block
                              !----------------------------------------------------------------------------------------------
                              !-- --------------------------------------------------------------------------------------------
                              if ( ANY(Mtaxa(k,:)==lineage) ) then !i.e. The block has the migrant among its taxa
                                  do l=1,n-1,1      !First of all, we store the original taxa array (Mtaxa(k,:))
                                    taxa4(l)=Mtaxa(k,l)
                                    nodes5B(l)=MnodesB(k,l)
                                    nodes5(l)=Mnodes(k,l)
                                  enddo
                                  taxa4(n)=Mtaxa(k,n)
                                  !--
                                  do l=1,ni(k)-1,1 !Before everything else we adjust Ages: which is easy because we only have to erase the first cell..
                                     Mages(k,l)=Mages(k,l+1)  !..and shift every one to the left. This is because the erased lineage has the youngest possible..
                                  enddo                       !..ages since it's inherited from a downstream block
                                  Mages(k,ni(k))=0
                                  !----
                                  pos = iMINMAXLOC1(1,n,Mtaxa(k,:),lineage,lineage,1,n)
                                  !Adjusting taxa
                                  do l=pos,ni(k)-1,1
                                     Mtaxa(k,l)=Mtaxa(k,l+1)
                                     MtaxaB(k,l)=MtaxaB(k,l+1)
                                  enddo
                                  Mtaxa(k,ni(k))=0
                                  MtaxaB(k,ni(k))=0
                                  !Adjusting nodes
                                  if (ni(k)>1) then !If there's only one taxon there are no nodes to adjust
                                     if (pos==1) then
                                        do l=pos,ni(k)-2,1
                                          Mnodes(k,l)=Mnodes(k,l+1)
                                        enddo
                                        Mnodes(k,ni(k)-1)=0
                                     elseif (pos==ni(k)) then
                                        Mnodes(k,ni(k)-1)=0
                                     else
                                        if (Nodages(1,Mnodes(k,pos-1)-n)>Nodages(1,Mnodes(k,pos)-n)) then  !The parental node of the removed taxon is the one at right
                                           do l=pos,ni(k)-2,1
                                              Mnodes(k,l)=Mnodes(k,l+1)
                                           enddo
                                           Mnodes(k,ni(k)-1)=0
                                        else
                                           do l=pos-1,ni(k)-2,1
                                              Mnodes(k,l)=Mnodes(k,l+1)
                                           enddo
                                           Mnodes(k,ni(k)-1)=0
                                        endif
                                     endif
                                  endif
                                  nbis=ni(k)
                                  ni(k)=ni(k)-1
                                  !And its respective MtaxaB and MnodesB-------------------------------
                                  allocate(x1(ni(k)),y1(ni(k)))
                                  do l=1,ni(k),1
                                     x1(l)=MtaxaB(k,l)
                                  enddo
                                  y1=(/(l, l=1,ni(k))/)
                                  CALL SSORT(x1,y1,ni(k),2)
                                  do l=1,ni(k),1
                                     MtaxaB(k,NINT(y1(l)))=l
                                  enddo
                                  deallocate(x1,y1)
                                  !---
                                  if (ni(k)>1) then
                                     allocate(x1(ni(k)-1),y1(ni(k)-1))
                                     do l=1,ni(k)-1,1
                                        x1(l)=Nodages(1,Mnodes(k,l)-n)
                                     enddo
                                     y1=(/(l, l=1,ni(k)-1)/)
                                     CALL SSORT(x1,y1,ni(k)-1,2)
                                     do l=1,ni(k)-1,1
                                        MnodesB(k,NINT(y1(l)))=l
                                        Mcoal_t(k,l)=x1(l)
                                     enddo
                                     MnodesB(k,ni(k))=0
                                     Mcoal_t(k,ni(k))=0
                                     deallocate(x1,y1)
                                  else
                                     MnodesB(k,:)=0
                                     Mcoal_t(k,:)=0
                                     MnodesB(k,:)=0
                                     Mcoal_t(k,:)=0
                                  endif
                                  !----Done
                                  b=k   !Update the focal block to the last processed. From now on it will search for blocks having this as a source
                              !----------------------------------------------------------------------------------------------
                              elseif  ( ANY(Mtaxa(k,:)==go_to ) ) then         !i.e. The block has the parental of the migrant among its taxa
                              !----------------------------------------------------------------------------------------------
                                  do l=1,n-1,1      !First of all, we store the original taxa array (Mtaxa(k,:))
                                    taxa4(l)=Mtaxa(k,l)
                                    nodes5B(l)=MnodesB(k,l)
                                    nodes5(l)=Mnodes(k,l)
                                  enddo
                                  taxa4(n)=Mtaxa(k,n)
                                  !--
                                  pos = iMINMAXLOC1(1,n,Mtaxa(k,:),go_to,go_to,1,n)
                                  Mtaxa(k,pos)=Em
                                  !---Done
                                  b=k
                              endif
                              !----------------------------------------------------------------------------------------------
                              !----------------------------------------------------------------------------------------------
                           endif
                       endif
                       DEALLOCATE(nodes2i)
                    enddo
                    !******************************************************************************************************************
                    !5.2   S E T T I N G  T H E  R E C E I V E R  L I N E A G E  A N D  C O A L E S C E N T  T I M E
                    !------------------------------------------------------------------------------------------------------------------
                    !---------First we obtain the go_to (the parental node of the migrant) because its id is the one for the future coalescent
                    if (ANY(taxa==emigrants(i,1))) then !i.e. if the migrant lineage is a Taxon
                       pos = iMINMAXLOC1(1,n,taxa,INT(emigrants(i,1)),INT(emigrants(i,1)),1,n)
                       CALL GET_GOTO(.false.,n-1,nodes3,pos,go_to,m1,m2)   !So go_to is the parental node
                    elseif (ANY(nodes==emigrants(i,1))) then !i.e. if the migrant lineage is a Node
                       pos = iMINMAXLOC1(1,n-1,nodes,INT(emigrants(i,1)),INT(emigrants(i,1)),1,n-1)
                       CALL GET_GOTO(.true.,n-1,nodes3,pos,go_to,m1,m2)    !So go_to is the parental node
                    endif
                    go_to=nodes(go_to)
                    !+++
                    d=NINT(emigrants(i,4)) !d gets the number of receiver population/block
                    etime=emigrants(i,2)
                    CALL Get_mig_coal(INT(emigrants(i,1)),go_to,d,n,en,dn,ni,gentime,Ne,growth,events,Mages,&
                    &Mtaxa,MtaxaB,Mnodes,MnodesB,Mcoal_t,emigrants(i,2),lineage,Ptime,older)
                    time=emigrants(i,2)
                    !-------------------------------------------------------------
                    !5.3 Now we find out if an additional migration ocurred in the newly created lineage
                    CALL Single_Lin_Mig(etime+Ptime,etime ,INT(emigrants(i,4)),gentime,NM,&
                    &dn+en,Mig_labels,MS,TM,MigMat,New_mig,time,receiver)
                    if (New_mig) then !If the lineage emigrated again the we abort the present migration, update its migration event and start again
                      do j=E+1,i+1,-1
                         emigrants(j,:)=emigrants(j-1,:)
                      enddo
                      emigrants(i+1,1)=emigrants(i,1)
                      emigrants(i+1,2)=time
                      emigrants(i+1,3)=emigrants(i,4)
                      emigrants(i+1,4)=receiver
                      emigrants(i+1,5)=emigrants(i,5)
                      E=E+1
                    endif
                    !-------------------------------------------------------------
                    !-------------------------------------------------------------
                    !5.5 Defining the ranges of the migrant and receiver lineages both in
                    !    ..the main arrays (taxa, nodes) and in the local arrays (Mtaxa, Mnodes)
                    if (emigrants(i,1)>n) then  !It is a node
                       pos = iMINMAXLOC1(1,n-1,nodes,INT(emigrants(i,1)),INT(emigrants(i,1)),1,n-1)
                       CALL GET_GOTO(.true.,n-1,nodes3,pos,dummy,m1,m2)
                       Nodages(2,NINT(emigrants(i,1)-n))=m1
                       Nodages(3,NINT(emigrants(i,1)-n))=m2
                    else                        !It is a taxon
                       pos = iMINMAXLOC1(1,n,taxa,INT(emigrants(i,1)),INT(emigrants(i,1)),1,n)
                       m1=pos
                       m2=pos
                    endif
                    if (lineage>n) then   !If it is a node
                       pos = iMINMAXLOC1(1,n-1,nodes,lineage,lineage,1,n-1)
                       CALL GET_GOTO(.true.,n-1,nodes3,pos,dummy,r1,r2)
                       Nodages(2,lineage-n)=r1
                       Nodages(3,lineage-n)=r2
                    else
                       pos = iMINMAXLOC1(1,n,taxa,lineage,lineage,1,n)
                       r1=pos
                       r2=pos
                    endif
                    !------------------------------------
                    c=NINT(emigrants(i,4)) !'c' gets the number of receiver population/block
                    if ((m1>=r1).and.(m2<=r2)) then !i.e. in case the migrant is inside the receiver group, we consider the receiver the whole sample in the receiver pop/block
                       Mr1=MIN(1,ni(c))
                       Mr2=ni(c)
                    else
                       do j=1,en,1
                          if ((events(j,4)==c).or.(events(j,5)==c)) then !With this we get the upper bound of our present pop/block
                             bound=events(j,1)
                          endif
                       enddo
                       if (ni(c)>0) then
                          if (emigrants(i,2)<bound) then
                             pos = iMINMAXLOC1(1,n,Mtaxa(c,:),taxa(r1),taxa(r1),1,n)
                             Mr1=pos
                             Mr2=Mr1+(r2-r1)
                          else
                             Mr1=1
                             Mr2=ni(c)
                          endif
                       else
                          Mr1=0
                          Mr2=0
                       endif
                    endif
                    !--
                    c=NINT(emigrants(i,3)) !'c' gets the number of donor population/block
                    pos = iMINMAXLOC1(1,n,Mtaxa(c,:),taxa(m1),taxa(m1),1,n)
                    Mm1=pos
                    Mm2=Mm1+(m2-m1)
                    !------------------------------------------------------------------------------------------------------------------------------
                    !------------------------------------------------------------------------------------------------------------------------------
                    !5.6 Cutting and glueing the migrant lineages
                    if (r2<m1) then
                       !Preparing special nodes for main arrays
                       if (m2==n) then !The node is at the end of taxa, so the survivor of its parental nodes is automatically the one at left
                          node2(1)=nodes(m1-1)                 !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                          node2(2)=emigrants(i,2) !Time to the coalescent
                       else
                          if (Nodages(1,nodes(m1-1)-n)<Nodages(1,nodes(m2)-n)) then !The oldest one is the node that survive the merging of the two flanking nodes of the migrant
                             node1(1)=nodes(m2)                !Node1 is the resulting from colapsing flanking nodes of migrant node
                             node1(2)=Nodages(1,nodes(m2)-n)
                             node2(1)=nodes(m1-1)              !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                             node2(2)=emigrants(i,2) !Time to the coalescent
                          else
                             node1(1)=nodes(m1-1)              !Node1 is the resulting from colapsing flanking nodes of migrant node
                             node1(2)=Nodages(1,nodes(m1-1)-n)
                             node2(1)=nodes(m2)                 !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                             node2(2)=emigrants(i,2) !Time to the coalescent
                          endif
                       endif
                       node3(1)=nodes(r2) !Node3 is the former right-flanking node of the acceptor group; its age and range are the same and adjusted later
                       !Translating the portion before and including the receiver range
                       do j=1,r2-1,1
                          taxa2(j)=taxa(j)
                          nodes2(j)=nodes(j)
                       enddo
                       taxa2(r2)=taxa(r2)
                       !Now relocating the migrating range
                       do j=1,m2-m1,1
                          taxa2(r2+j)=taxa(m1+j-1)
                          nodes2(r2+j)=nodes(m1+j-1)
                       enddo
                       taxa2(r2+m2-m1+1)=taxa(m2)
                       !..and the in-betweeners now go after
                       if (m1>r2+1) then
                          do j=r2+1,m1-2,1
                             taxa2(j+m2-m1+1)=taxa(j)
                             nodes2(j+m2-m1+1)=nodes(j)
                          enddo
                          taxa2(m2)=taxa(m1-1)
                       endif
                       !..and finally the tail
                       if (m2<n) then
                          do j=m2+1,n-1,1
                             taxa2(j)=taxa(j)
                             nodes2(j)=nodes(j)
                          enddo
                          taxa2(n)=taxa(n)
                       endif
                       !-----
                       taxa=taxa2
                       nodes=nodes2
                       !Now inserting the nodes that required special treatment:
                       if (m2<n) then
                          nodes(m2)=NINT(node1(1))
                       endif
                       nodes(r2)=NINT(node2(1))
                       if (r2+(m2-m1+1)<n) then
                          if (node3(1)/=node2(1)) then
                             nodes(r2+(m2-m1+1))=NINT(node3(1))
                          endif
                       endif
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       if (i<E) then
                          !This is for Mtaxa and Mnodes
                          a=NINT(emigrants(i,4)) !receiver block
                          b=NINT(emigrants(i,3)) !donor block
                          !Preparing special nodes for Marrays (they can differ)
                          if ( (Mm2/=ni(b)) .and. (Mm1/=1) ) then !Only nodeM1 and nodeM3 need to get registered (node2 the node that desappears and is created gets the id from main arrays)
                             if (MnodesB(b,Mm1-1)<MnodesB(b,Mm2)) then !The oldest one is the node that survive the merging of the two flanking nodes of the migrant
                                nodeM1=Mnodes(b,Mm2)                 !Node1 is the resulting from colapsing flanking nodes of migrant node
                             else
                                nodeM1=Mnodes(b,Mm1-1)
                             endif
                          endif
                          if (Mr2<ni(a)) then
                             nodeM3=Mnodes(a,Mr2)
                          endif
                          !--Move the taxa at right of the receiver group
                          if (ni(a)>0) then
                             if (Mr2<ni(a)) then
                                c=ni(a)+(m2-m1+1)
                                Mtaxa(a,c)=Mtaxa(a,ni(a))
                                MtaxaB(a,c)=MtaxaB(a,ni(a))
                                do j=ni(a)-1,Mr2+1,-1  !From the original sample size in the receiver pop/block + migrant range, to the original sample size
                                   c=c-1
                                   Mtaxa(a,c)=Mtaxa(a,j)           !We send the taxa at right from the receiver group to the end
                                   MtaxaB(a,c)=MtaxaB(a,j)         !This is the only way to do it without using a third array
                                   Mnodes(a,c)=Mnodes(a,j)
                                enddo
                                Mnodes(a,c-1)=Mnodes(a,Mr2)
                             endif
                          endif
                          !--Get in the migrating group
                          c=Mr2
                          do j=Mm1,Mm2-1,1
                             c=c+1
                             Mtaxa(a,c)=Mtaxa(b,j)
                             Mnodes(a,c)=Mnodes(b,j)
                             MtaxaB(a,c)=ni(a)+j-Mm1+1  !It gets consecutive numbers continuing the sample size
                          enddo
                          Mtaxa(a,c+1)=Mtaxa(b,Mm2)
                          MtaxaB(a,c+1)=ni(a)+(Mm2-Mm1+1)
                          !--Now filling the hole left by the migrating group with cells at right
                          if (Mm2<ni(b)) then
                             c=Mm1
                             do j=Mm2+1,ni(b)-1,1
                                Mtaxa(b,c)=Mtaxa(b,j)
                                Mnodes(b,c)=Mnodes(b,j)
                                MtaxaB(b,c)=MtaxaB(b,j)
                                c=c+1
                             enddo
                             Mtaxa(b,c)=Mtaxa(b,ni(b))
                             MtaxaB(b,c)=MtaxaB(b,ni(b))
                          endif
                          !--
                          ni(a)=ni(a)+(Mm2-Mm1+1)
                          ni(b)=ni(b)-(Mm2-Mm1+1)
                          !--Mages
                          if (a>dn) then
                             pos = rMINMAXLOC1(1,en,events(:,8),DBLE(a),DBLE(a),1,en)
                             bound=events(pos,1) !Age of the block
                          else
                             bound=0
                          endif
                          do j=1,ni(a),1
                             Mages(a,j)=MIN(bound,ages(Mtaxa(a,j)))
                          enddo
                          if (b>dn) then
                             pos = rMINMAXLOC1(1,en,events(:,8),DBLE(b),DBLE(b),1,en)
                             bound=events(pos,1) !Age of the block
                          else
                             bound=0
                          endif
                          do j=1,ni(b),1
                             Mages(b,j)=MIN(bound,ages(Mtaxa(b,j)))
                          enddo
                          !...but wait. The special taxa and nodes had not been updated yet
                          !Now inserting the nodes that required special treatment:
                          if ((Mm2<ni(b)+(Mm2-Mm1+1)).and.(Mm1>1)) then
                             Mnodes(b,Mm1-1)=NINT(nodeM1)
                          endif
                          if ((Mr2<ni(a)).and.(Mr2>0)) then
                             if (older) then !This means: if the coalescent happened outside the block
                                zbis=MAXVAL(Mnodes)+1   !..then we don't name it as node2, we use a new name (id); this is because if we use node2 it could be different than the node2 at the main arrays
                                Mnodes(a,Mr2)=zbis     !..and the assignment of ages (which is given by name) could get wrong in posterior cycles
                                Nodages(1,zbis-n)=node2(2)+1
                             else
                                Mnodes(a,Mr2)=NINT(node2(1))
                             endif
                          endif
                          if (Mr2+(Mm2-Mm1+1)<ni(a)) then
                             Mnodes(a,Mr2+(Mm2-Mm1+1))=NINT(nodeM3)
                          endif
                       endif
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       !-------------------------------------------------------------
                       !5.8 Updating ranges and ages
                       if (m2<n) then
                          Nodages(1,NINT(node1(1)-n))=node1(2)  !For node1 (the result of collapsing flanking nodes of the migrant)
                       endif
                       Nodages(1,NINT(node2(1)-n))=node2(2)  !For node2 (the new coalescent of migrant and acceptor lineages)
                    elseif (m2<r1) then !--------------------------------------------------------------------------------------------------------------
                       IF (m1==1) THEN  !In case the migrant group is the first in taxa, then the surviving node is necessarily the one at right (in position m2)
                          node2(1)=nodes(m2)              !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migrant group
                          node2(2)=emigrants(i,2) !Time to the coalescent
                       ELSE
                          IF (Nodages(1,nodes(m1-1)-n)<Nodages(1,nodes(m2)-n)) THEN !The oldest one is the node that survive the merging of the two flanking nodes of the migrant
                             node1(1)=nodes(m2)                !Node1 is the resulting from colapsing flanking nodes of migrant node
                             node1(2)=Nodages(1,nodes(m2)-n)
                             node2(1)=nodes(m1-1)              !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                             node2(2)=emigrants(i,2) !Time to the coalescent
                          ELSE
                             node1(1)=nodes(m1-1)              !Node1 is the resulting from colapsing flanking nodes of migrant node
                             node1(2)=Nodages(1,nodes(m1-1)-n)
                             node2(1)=nodes(m2)                 !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                             node2(2)=emigrants(i,2) !Time to the coalescent
                          ENDIF
                       ENDIF
                       node3(1)=nodes(r1-1) !Node3 is the former left-flanking node of the acceptor group; its age and range are the same and adjusted later
                       if (m1>1) then
                          !Translating the portion before the receiver range
                          do j=1,m1-2,1
                             taxa2(j)=taxa(j)
                             nodes2(j)=nodes(j)
                          enddo
                          taxa2(m1-1)=taxa(m1-1)
                       endif
                       !..and the in-betweeners now go after
                       if (m2<r1-1) then
                          do j=m2+1,r1-2,1
                             taxa2(m1-m2+j-1)=taxa(j)
                             nodes2(m1-m2+j-1)=nodes(j)
                          enddo
                          taxa2(m1-m2+r1-2)=taxa(r1-1)
                       endif
                       !Now relocating the migrating range
                       do j=m1,m2-1,1
                          taxa2(r1-1-(m2-m1)+j-m1)=taxa(j)
                          nodes2(r1-1-(m2-m1)+j-m1)=nodes(j)
                       enddo
                       taxa2(r1-1)=taxa(m2)
                       !..and finally the tail
                       do j=r1,n-1,1
                          taxa2(j)=taxa(j)
                          nodes2(j)=nodes(j)
                       enddo
                       taxa2(n)=taxa(n)
                       !-----
                       taxa=taxa2
                       nodes=nodes2
                       !Now inserting the nodes that required special treatment:
                       if (m1>1) then
                          nodes(m1-1)=NINT(node1(1))  !(i)
                       endif
                       nodes(r1-1)=NINT(node2(1))     !(ii) Then we put the newly created node (the one joining receiver and migrant groups)
                       if (r1-(m2-m1+1)-1>0) then        !We restitute the node ar right of the receiver group; now it'll be at right of the receiver+migrant groups
                          if (node3(1)/=node2(1)) then   !...unless the node at right of the receiver group is also the node that dissapears (which happens)
                             nodes(r1-(m2-m1+1)-1)=NINT(node3(1))
                          endif
                       endif
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       if (i<E) then
                          !This is for Mtaxa and Mnodes
                          a=NINT(emigrants(i,4)) !receiver block
                          b=NINT(emigrants(i,3)) !donor block
                          !Preparing special nodes for Marrays (they can differ)
                          if ( (Mm2/=ni(b)) .and. (Mm1/=1) ) then !Only nodeM1 and nodeM3 need to get registered (node2 the node that desappears and is created gets the id from main arrays)
                             if (MnodesB(b,Mm1-1)<MnodesB(b,Mm2)) then !The oldest one is the node that survive the merging of the two flanking nodes of the migrant
                                nodeM1=Mnodes(b,Mm2)                 !Node1 is the resulting from colapsing flanking nodes of migrant node
                             else
                                nodeM1=Mnodes(b,Mm1-1)
                             endif
                          endif
                          if (Mr1>1) then
                             nodeM3=Mnodes(a,Mr1-1)   !Node3 is the former left-flanking node of the acceptor group
                          endif
                          !--Move the taxa at right of the receiver group
                          if(ni(a)>0) then
                             c=ni(a)+(m2-m1+1)
                             Mtaxa(a,ni(a)+(m2-m1+1))=Mtaxa(a,ni(a))
                             MtaxaB(a,ni(a)+(m2-m1+1))=MtaxaB(a,ni(a))
                             do j=ni(a)-1,Mr1,-1  !From the original sample size in the receiver pop/block + migrant range, to the original sample size
                                c=c-1
                                Mtaxa(a,c)=Mtaxa(a,j)           !We send the taxa at right from the receiver group to the end
                                MtaxaB(a,c)=MtaxaB(a,j)         !This is the only way to do it without using a third array
                                Mnodes(a,c)=Mnodes(a,j)
                             enddo
                          endif
                          !--Get in the migrating group
                          if (Mr1>0) then
                             c=Mr1-1
                          else
                             c=0
                          endif
                          do j=Mm1,Mm2-1,1
                             c=c+1
                             Mtaxa(a,c)=Mtaxa(b,j)
                             Mnodes(a,c)=Mnodes(b,j)
                             MtaxaB(a,c)=ni(a)+j-Mm1+1  !It gets consecutive numbers continuing the sample size
                          enddo
                          Mtaxa(a,c+1)=Mtaxa(b,Mm2)
                          MtaxaB(a,c+1)=ni(a)+(Mm2-Mm1+1)
                          !--Now filling the hole left by the migrating group with cells at right
                          if (Mm2<ni(b)) then
                             c=Mm1
                             do j=Mm2+1,ni(b)-1,1
                                Mtaxa(b,c)=Mtaxa(b,j)
                                Mnodes(b,c)=Mnodes(b,j)
                                MtaxaB(b,c)=MtaxaB(b,j)
                                c=c+1
                             enddo
                             Mtaxa(b,c)=Mtaxa(b,ni(b))
                             MtaxaB(b,c)=MtaxaB(b,ni(b))
                          endif
                          !--
                          ni(a)=ni(a)+(Mm2-Mm1+1)
                          ni(b)=ni(b)-(Mm2-Mm1+1)
                          !--Mages
                          if (a>dn) then
                             pos = rMINMAXLOC1(1,en,events(:,8),DBLE(a),DBLE(a),1,en)
                             bound=events(pos,1) !Age of the block
                          else
                             bound=0
                          endif
                          do j=1,ni(a),1
                             Mages(a,j)=MIN(bound,ages(Mtaxa(a,j)))
                          enddo
                          if (b>dn) then
                             pos = rMINMAXLOC1(1,en,events(:,8),DBLE(b),DBLE(b),1,en)
                             bound=events(pos,1) !Age of the block
                          else
                             bound=0
                          endif
                          do j=1,ni(b),1
                             Mages(b,j)=MIN(bound,ages(Mtaxa(b,j)))
                          enddo
                          !...but wait. The special taxa and nodes had not been updated yet
                          !Now inserting the nodes that required special treatment:
                          if ((Mm1>1).and.(Mm2<ni(b)+(Mm2-Mm1+1))) then !We update the left-flanking-migrating-group node only if:
                             Mnodes(b,Mm1-1)=NINT(nodeM1)                   !(1) There existed one (Mm1>1) and (2) The migrating group was not in the left tail
                          endif
                          if (Mr1>0) then
                             if (older) then !This means: if the coalescent happened outside the block
                                zbis=MAXVAL(Mnodes)+1   !..then we don't name it as node2, we use a new name (id); this is because if we use node2 it could be different than the node2 at the main arrays
                                Mnodes(a,Mr1+(Mm2-Mm1+1)-1)=zbis      !..and the assignment of ages (which is given by name) could get wrong in posterior cycles
                                Nodages(1,zbis-n)=node2(2)+1
                             else
                                Mnodes(a,Mr1+(Mm2-Mm1+1)-1)=NINT(node2(1))   !This is different than before because in Mtaxa and Mnodes the receiver group is not static when the migrant is at left
                             endif
                          endif
                          if (Mr1>1) then
                             Mnodes(a,Mr1-1)=NINT(nodeM3)
                          endif
                       endif
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       !-------------------------------------------------------------
                       !5.8 Updating ranges and ages
                       if (m1>1) then
                          Nodages(1,NINT(node1(1)-n))=node1(2)  !For node1 (the result of collapsing flanking nodes of the migrant)
                       endif
                       Nodages(1,NINT(node2(1)-n))=node2(2)  !For node2 (the new coalescent of migrant and acceptor lineages)
                    elseif ((m1>=r1).and.(m2<=r2)) then !i.e. if the migrating node is insidethe receiver node
                       if (m2==n) then !The node is at the end of taxa, so the survivor of its parental nodes is automatically the one at left
                          node2(1)=nodes(m1-1)        !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                          node2(2)=emigrants(i,2) !Time to the coalescent
                       elseif (m1==1) then  !In case the migrant group is the first in taxa, then the surviving node is necessarily the one at right (in position m2)
                          node2(1)=nodes(m2)          !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migrant group                           node2(2)=emigrants(i,2) !Time to the coalescent
                       else
                          if (Nodages(1,nodes(m1-1)-n)<Nodages(1,nodes(m2)-n)) THEN !The oldest one is the node that survive the merging of the two flanking nodes of the migrant
                             node1(1)=nodes(m2)                !Node1 is the resulting from colapsing flanking nodes of migrant node
                             node1(2)=Nodages(1,nodes(m2)-n)
                             node2(1)=nodes(m1-1)              !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                             node2(2)=emigrants(i,2) !Time to the coalescent
                          else
                             node1(1)=nodes(m1-1)              !Node1 is the resulting from colapsing flanking nodes of migrant node
                             node1(2)=Nodages(1,nodes(m1-1)-n)
                             node2(1)=nodes(m2)                 !Node2 is the new one: the coalescent of the migrant and acceptor lineages, it get the name of the node dissapeared from the flanking ones to the migran group
                             node2(2)=emigrants(i,2) !Time to the coalescent
                          endif
                       endif
                       !---
                       if (m1>1) then
                          !Translating the portion before the migrant range (including a part of the receiver)
                          do j=1,m1-1,1
                             taxa2(j)=taxa(j)
                             nodes2(j)=nodes(j)
                          enddo
                       endif
                       !..and the second half of the receiver group goes next
                       if (m2<=r2-1) then
                          if (m2<r2-1) then
                             do j=m2+1,r2-1,1
                                taxa2(m1-m2+j-1)=taxa(j)
                                nodes2(m1-m2+j-1)=nodes(j)
                             enddo
                          endif
                          taxa2(m1-m2+r2-1)=taxa(r2)
                       endif
                       !Now relocating the migrating range at right of the receiver group
                       do j=m1,m2-1,1
                          taxa2(r2-m2+j)=taxa(j)
                          nodes2(r2-m2+j)=nodes(j)
                       enddo
                       taxa2(r2)=taxa(m2)
                       !..and finally the tail
                       if (r2<n) then
                          nodes2(r2)=nodes(r2)
                          do j=r2+1,n-1,1
                             taxa2(j)=taxa(j)
                             nodes2(j)=nodes(j)
                          enddo
                          taxa2(n)=taxa(n)
                       endif
                       !-----
                       taxa=taxa2
                       nodes=nodes2
                       !Now inserting the nodes that required special treatment:
                       if (m1>1) then
                          nodes(m1-1)=NINT(node1(1))
                       endif
                       nodes(r2-(m2-m1+1))=NINT(node2(1))
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       if (i<E) then
                          !This is for Mtaxa and Mnodes
                          a=NINT(emigrants(i,4)) !receiver block
                          b=NINT(emigrants(i,3)) !donor block
                          !Preparing special nodes for Marrays (they can differ)
                          if ( (Mm2/=ni(b)) .and. (Mm1/=1) ) then !Only nodeM1 and nodeM3 need to get registered (node2 the node that desappears and is created gets the id from main arrays)
                             if (MnodesB(b,Mm1-1)<MnodesB(b,Mm2)) then !The oldest one is the node that survive the merging of the two flanking nodes of the migrant
                                nodeM1=Mnodes(b,Mm2)                 !Node1 is the resulting from colapsing flanking nodes of migrant node
                             else
                                nodeM1=Mnodes(b,Mm1-1)
                             endif
                          endif
                          !--Here the receiver group is not in the receiver population, so we only put the migrant group aside
                          c=ni(a)
                          do j=Mm1,Mm2-1,1
                             c=c+1
                             Mtaxa(a,c)=Mtaxa(b,j)
                             Mnodes(a,c)=Mnodes(b,j)
                             MtaxaB(a,c)=ni(a)+j-Mm1+1  !It gets consecutive numbers continuing the sample size
                          enddo
                          Mtaxa(a,c+1)=Mtaxa(b,Mm2)
                          MtaxaB(a,c+1)=ni(a)+(Mm2-Mm1+1)
                          !--Now filling the hole left by the migrating group with cells at right
                          if (Mm2<ni(b)) then
                             c=Mm1
                             do j=Mm2+1,ni(b)-1,1
                                Mtaxa(b,c)=Mtaxa(b,j)
                                Mnodes(b,c)=Mnodes(b,j)
                                MtaxaB(b,c)=MtaxaB(b,j)
                                c=c+1
                             enddo
                             Mtaxa(b,c)=Mtaxa(b,ni(b))
                             MtaxaB(b,c)=MtaxaB(b,ni(b))
                          endif
                          !--
                          ni(a)=ni(a)+(Mm2-Mm1+1)
                          ni(b)=ni(b)-(Mm2-Mm1+1)
                          !--Mages
                          if (a>dn) then
                             pos = rMINMAXLOC1(1,en,events(:,8),DBLE(a),DBLE(a),1,en)
                             bound=events(pos,1) !Age of the block
                          else
                             bound=0
                          endif
                          do j=1,ni(a),1
                             Mages(a,j)=MIN(bound,ages(Mtaxa(a,j)))
                          enddo
                          if (b>dn) then
                             pos = rMINMAXLOC1(1,en,events(:,8),DBLE(b),DBLE(b),1,en)
                             bound=events(pos,1) !Age of the block
                          else
                             bound=0
                          endif
                          do j=1,ni(b),1
                             Mages(b,j)=MIN(bound,ages(Mtaxa(b,j)))
                          enddo
                          !...but wait. The special taxa and nodes had not been updated yet
                          !Now inserting the nodes that required special treatment:
                          if ((Mm2<ni(b)+(Mm2-Mm1+1)).and.(Mm1>1)) then
                             Mnodes(b,Mm1-1)=NINT(nodeM1)
                          endif
                          if (ni(a)-(Mm2-Mm1+1)>=1) then
                             if (older) then !This means: if the coalescent happened outside the block
                                zbis=MAXVAL(Mnodes)+1   !..then we don't name it as node2, we use a new name (id); this is because if we use node2 it could be different than the node2 at the main arrays
                                Mnodes(a,ni(a)-(Mm2-Mm1+1))=zbis      !..and the assignment of ages (which is given by name) could get wrong in posterior cycles
                                Nodages(1,zbis-n)=node2(2)+1
                             else
                                Mnodes(a,ni(a)-(Mm2-Mm1+1))=NINT(node2(1))
                             endif
                          endif
                       endif
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       !5.8 Updating ranges and ages
                       if ((m1>1).and.(m2<n)) then
                          Nodages(1,NINT(node1(1)-n))=node1(2)  !For node1 (the result of collapsing flanking nodes of the migrant)
                       endif
                       Nodages(1,NINT(node2(1)-n))=node2(2)  !For node2 (the new coalescent of migrant and acceptor lineages)
                       !-------------------------------------------------------------
                    endif
                    !-----------------------------------------------------------------------------------
                    !-----------------------------------------------------------------------------------
                    if (i<E) then
                        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        !  A D J U S T M E N T S   I N   E M I G R A N T S
                        !The required adjustments are:
                        !(1) Update the name of the node that dissapears (the incidence of the migrant), in case it is a migrant too
                        if (ANY(emigrants(:,1)==node2(1))) then !Remember that node2 gets the node that dissapears/is_created

                           pos = iMINMAXLOC1(1,n-1,nodes4,INT(node2(1)),INT(node2(1)),1,n-1) !Since we are dealing with a node that changes name, it is always a node
                           call Get_downstream_lin(n-1,pos,nodes3,A1,TN1,A2,TN2)
                           if (TN1==0) then
                              A1=taxa3(A1)
                           else
                              A1=nodes4(A1)
                           endif
                           if (TN2==0) then
                              A2=taxa3(A2)
                           else
                              A2=nodes4(A2)
                           endif
                           if (emigrants(i,1)==A1) then  !If the left incidence is actually the current emigrant
                              Em=A2          !...then the other incidence is the one that heritage the name of the dissapeared node
                              Lm=TN2
                           elseif (emigrants(i,1)==A2) then !If the left incidence is actually the current emigrant
                              Em=A1             !...then the other incidence is the one that heritage the name of the dissapeared node
                              Lm=TN1
                           endif
                           pos = rMINMAXLOC1(1,2*n,emigrants(:,1),node2(1),node2(1),1,2*n)
                           emigrants(pos,1)=Em     !...then the other incidence is the one that inherit the name of the dissapeared node
                           emigrants(pos,5)=Lm
                           j=i+1
                           do    !Before we have to make a check to verify that: The lineage that is going to inherit the migration event was not already committed for another migration event
                              if ( (emigrants(j,1)==Em) .and. (j/=pos) ) then !This means that
                                 do l=MIN(pos,j),E,1 !This instruction is actually erasing the x-th row of emigrants (x=youngest migration event between position(1) and the other)
                                    emigrants(l,:)=emigrants(l+1,:)
                                 enddo
                                 E=E-1
                              endif
                              if (j>=E) exit
                              j=j+1
                           enddo
                        endif
                        !(2) Update the name of the lineage that receives the migrant, in case it is a migrant too
                        do j=i+1,E,1
                           if (emigrants(i,1)==lineage) then
                              pos = rMINMAXLOC1(1,2*n,emigrants(:,1),DBLE(lineage),DBLE(lineage),1,2*n)
                              if (emigrants(pos,2)>emigrants(i,2)) then !What this means is that the migrant lineage coalesced another migrant lineage at a time younger than the migration of it
                                 emigrants(pos,1)=node2(1)              !So the actual emigrant is not the previous one anymore but the newly created node
                              endif
                           endif
                        enddo
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       !5.9 Cleaning the tails of the arrays and creating B arrays (MtaxaB, MnodesB)
                       if (ni(a)>0) then
                          do j=ni(a)+1,n,1
                             Mtaxa(a,j)=0
                             MtaxaB(a,j)=0
                             Mages(a,j)=0
                          enddo
                          do j=ni(a),n-1,1
                             Mnodes(a,j)=0
                             MnodesB(a,j)=0
                             Mcoal_t(a,j)=0
                          enddo
                       else
                           Mtaxa(a,:)=0
                           MtaxaB(a,:)=0
                           Mages(a,:)=0
                           Mnodes(a,:)=0
                           MnodesB(a,:)=0
                           Mcoal_t(a,:)=0
                       endif
                       if (ni(b)>0) then
                          do j=ni(b)+1,n,1
                             Mtaxa(b,j)=0
                             MtaxaB(b,j)=0
                             Mages(b,j)=0
                          enddo
                          do j=ni(b),n-1,1
                             Mnodes(b,j)=0
                             MnodesB(b,j)=0
                             Mcoal_t(b,j)=0
                          enddo
                       else
                          Mtaxa(b,:)=0
                          MtaxaB(b,:)=0
                          Mages(b,:)=0
                          Mnodes(b,:)=0
                          MnodesB(b,:)=0
                          Mcoal_t(b,:)=0
                       endif
                       !-----------------
                       !---MtaxaB ; in MtaxaB from the donor pop/block there are missing taxa so we have to make numering continuous
                       if (ni(b)>0) then    !Not necessary for MtaxaB of the receiver because it only gets new numbers that are added to the count
                          allocate(x1(ni(b)),y1(ni(b)))
                          do j=1,ni(b),1
                             x1(j)=MtaxaB(b,j)
                          enddo
                          y1=(/(j, j=1,ni(b))/)
                          CALL SSORT(x1,y1,ni(b),2)
                          do j=1,ni(b),1
                             MtaxaB(b,NINT(y1(j)))=j
                          enddo
                          deallocate(x1,y1)
                       endif
                       !---MnodesB
                       if (ni(a)>1) then
                          allocate(x1(ni(a)-1),y1(ni(a)-1))
                          do j=1,ni(a)-1,1
                             x1(j)=Nodages(1,Mnodes(a,j)-n)
                          enddo
                          y1=(/(j, j=1,ni(a)-1)/)
                          CALL SSORT(x1,y1,ni(a)-1,2)
                          do j=1,ni(a)-1,1
                             MnodesB(a,NINT(y1(j)))=j
                             Mcoal_t(a,j)=x1(j)
                          enddo
                          deallocate(x1,y1)
                       endif
                       !---
                       if (ni(b)>1) then
                          allocate(x1(ni(b)-1),y1(ni(b)-1))
                          do j=1,ni(b)-1,1
                             x1(j)=Nodages(1,Mnodes(b,j)-n)
                          enddo
                          y1=(/(j, j=1,ni(b)-1)/)
                          CALL SSORT(x1,y1,ni(b)-1,2)
                          do j=1,ni(b)-1,1
                            MnodesB(b,NINT(y1(j)))=j
                            Mcoal_t(b,j)=x1(j)
                          enddo
                          deallocate(x1,y1)
                       endif
                       !++++++++++++++++++++++++++++++++++++++++++++++++++++
                    endif
                ENDIF
                !******************************************************************************************************************
                !******************************************************************************************************************
                IF (i==E) EXIT
            enddo
            deallocate(emigrants)
         endif
      endif   
      !-----------------------------------------------------------------------------------------------------------------------------
      !(4) Filling Coal_t (Have to get times, sort them and rename nodes with the proper positions)
      !-----------------------------------------------------------------------------------------------------------------------------
      do i=1,n-1,1
         coal_t(i)=Nodages(1,nodes(i)-n)
         NodRanges(1,i)=NINT(Nodages(2,nodes(i)-n))
         NodRanges(2,i)=NINT(Nodages(3,nodes(i)-n))
         y(i)=i
      enddo
      call SSORT(coal_t,y,n-1,2)
      do i=1,n-1,1
         nodes(NINT(y(i)))=i
         y(i)=i
      enddo
      !Finally, getting the ranges of the nodes that do not have them yet
      do i=1,n-1,1
         if ((NodRanges(1,i)==0).or.(NodRanges(2,i)==0)) then
            !Upper bound
            b=h
            switch=.false.
            do
               if (b==n-1) then
                  switch=.true.
               elseif ( coal_t(nodes(b+1)) > coal_t(nodes(i)) ) then
                  switch=.true.
               endif
               if (switch) exit
               b=b+1
            enddo
            !Lower bound
            a=h
            switch=.false.
            do
               if (a==1) then
                  switch=.true.
               elseif ( coal_t(nodes(a-1)) > coal_t(nodes(i)) ) then
                  switch=.true.
               endif
               if (switch) exit
               a=a-1
            enddo
            NodRanges(1,i)=a
            NodRanges(2,i)=b+1
         endif
      enddo



    return


   end subroutine tree_maker




   subroutine Migration(h,bigN,n,gentime,ages,taxa,nodes,coal_t,nm,dnen,Mig_labels,MS,TM,MigMat,emigrants)

      use routines
      use irandom
      use luxury

      implicit none

      !Externals
      integer, intent(in):: h                                   !Id of the focal population/block
      integer, intent(in):: bigN                                !Overall sample size (for sizing correctly emigrants array)
      integer, intent(in):: n                                   !Sample size (number of lineages coalesced in the focal population/block)
      real(8), intent(in):: gentime                             !Generation time
      real(8), dimension(n), intent(in):: ages                  !Same as always
      integer, dimension(n), intent(in):: taxa                  !Same as always
      integer, dimension(n-1), intent(in):: nodes               !Same as always
      real(8), dimension(n-1), intent(in):: coal_t              !Same as always
      integer, intent(in):: nm                                  !Number of migration matrix
      integer, intent(in):: dnen                                !Size of the migration matrix (it is actually not the effective size)
      integer, dimension(nm,dnen), intent(in):: mig_labels      !Since the whole matrix can be space-consuming, we only enter the (labels) columns/rows with values 
      integer, dimension(nm), intent(in):: MS                   !Array with the sizes of those matrix
      real, dimension(nm,2), intent(in):: TM                    !Timespans for the validity of the corresponding matrix
      real(8), dimension(nm,dnen,dnen), intent(in):: MigMat     !Migration matrix
      real(8), dimension(2*bigN,5):: emigrants                  !1st col:= Ids of the nodes/taxa that emigrate; 2nd col:= time; 3rd col:= receiver population/block; 4th col:= 0 if taxon and 1 if node
      !Internals
      integer:: i, j, k, m, o, z, z0
      real(8):: lo, up, mu, u, leng, length                     !Boundaries of the time window where a specific migration matrix applies; mu:= intensity of migration (lambda of a Poisson density)
      real, dimension(1):: RVEC                                 !Random number
      integer:: a, b, winner, go_to
      integer:: nl, na, lin, ad, s, me
      integer, allocatable, dimension(:,:):: list, tlist        !They are lists of the migration events only storing the donor and receiver populations
      real(8):: mbound, nbound
      logical:: switch, younger
      integer:: dummy1, dummy2
      !--------------------------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------------------------
      !START
      z=COUNT(emigrants(:,1)/=0)
      z0=z
      do i=1,NM,1
         if (ANY(mig_labels(i,:)==h)) then  !If our population appears in the labels, i.e. is considered for migration
            lo=TM(i,1)
            up=TM(i,2)
            nl = n - COUNT(coal_t<lo) - COUNT(ages>lo) !The number of lineages entring the time window is actually n minus: one for each coalescent
            length = nl*(up-lo)  !First we get the lenght of the tree in the interval
            do j=1,n-1,1
               if ((coal_t(j)>lo).and.(coal_t(j)<up)) then  !..which is the number of lineages times the lenght of the interval
                  length = length - (up-coal_t(j))         !..minus the difference of the upper bound and a node for each none
               endif
            enddo
            do j=1,n,1                                       !It is also needed to substract the age of a taxa
               if ((ages(j)>lo).and.(ages(j)<up)) then
                  length = length - (MIN(ages(j),up)-lo)      !The minimum is in case the taxa is actually older than the time window
               endif
            enddo
            !This is only to find the row/column in the matrix that corresponds to our population/block
            do j=1,MS(i),1
               if (mig_labels(i,j)==h) then
                  s=j
               endif
            enddo
            !Now finding out if emigration ocurrs
            do j=1,MS(i),1  !For each column..
               if (MigMat(i,s,j)/=0) then
                  mu=MigMat(i,s,j)*length/gentime
                  me=random_Poisson(real(mu),.true.)
                  if (.NOT.allocated(list)) then
                     allocate(list(me,2))
                     a=0
                  else
                     a=SIZE(list,1)
                  endif
                  allocate(tlist(me+a,2))
                  do k=1,a,1
                     tlist(k,1)=list(k,1)
                     tlist(k,2)=list(k,2)
                  enddo
                  if (me>0) then        !If succesful then
                     do k=1,me,1
                        tlist(a+k,1)=s                !1st column: id of the donor population/block
                        tlist(a+k,2)=mig_labels(i,j)  !2nd column: the id of the receiver population/block
                     enddo
                  endif
                  deallocate(list)
                  allocate(list(me+a,2))
                  list=tlist
                  deallocate(tlist)
               endif
            enddo
            !-------------------
            !Now producing those emigrations
            if (allocated(list)) then
               k=SIZE(list,1)
            else
               k=0
            endif
            if (k>0) then  !This is for not doing anything if there were no events!
               na=COUNT((ages<up).and.(ages>lo))   !Number of taxa whose ages fell in the time window
               do j=1,k,1
                  CALL RANLUX (RVEC, 1)  !This random is to decide in what time band (inter-arrival times) the migration happened
                  u=RVEC(1)
                  mbound=lo
                  a = rMINMAXLOC2(-1,n-1,coal_t,lo,dble(1000000000000000.0),1,n-1)
                  b = rMINMAXLOC2(-1,n,ages,lo,dble(1000000000000000.0),1,n)
                  lin=nl            !Next instructions works as follows:
                  switch=.false.    !Each inter-arrival time band (between nodes and also between ages of taxa if any) is screened
                  leng=0
                  do                !to count the lineages crossing, in order to get their lenghts' sum (len), and randomly choose the inter-arrival band with probability len/length
                     if ((b==0).or.(b>n)) then !If there are no ages in the time window
                        if ((a>n-1).or.(lo>coal_t(n-1))) then !This is if we are in the lineage being the MRCA so there are no upper nodes anymore
                           leng=leng+lin*(up-mbound)
                           nbound=mbound
                           mbound=up
                        else
                           leng=leng+lin*(MIN(coal_t(a),up)-mbound)
                           nbound=mbound
                           mbound=MIN(coal_t(a),up)
                        endif
                        a=a+1
                        ad=-1
                     elseif (coal_t(a)<ages(b)) then
                        leng=leng+lin*(MIN(coal_t(a),up)-mbound)
                        nbound=mbound
                        mbound=MIN(coal_t(a),up)
                        a=a+1
                        ad=-1
                     elseif (coal_t(a)>ages(b)) then
                        leng=leng+lin*(MIN(up,ages(b))-mbound)
                        nbound=mbound
                        mbound=MIN(ages(b),up)
                        b=b+1
                        ad=1
                     endif
                     if (u<leng/length) then !This means the migration point fell here!
                        !So we need to choose randomly the lineage
                        CALL RANLUX (RVEC, 1)
                        u=RVEC(1)
                        winner=NINT(lin*u+0.5) !Finally we now which lineage is
                        m=0
                        o=1    !This instruction works by scanning every taxa and node to see if they belong to the subset that cross the band
                        do     !..and stopping when we get the chosen one (so it is only to get the identity of the chosen one)
                           if (MOD(o,2)==0) then !This is for nodes
                              if (coal_t(nodes(o/2))<mbound) then !This 'if' is an assurance for not considering lineages that start at an older point than the time band
                                 younger=.true.
                                 CALL GET_GOTO(.true.,n-1,dble(nodes),o/2,go_to,dummy1,dummy2)
                              else
                                 younger=.false.
                              endif
                           else                  !This is for taxa
                              if (ages(taxa((o+1)/2))<mbound) then
                                 younger=.true.
                                 CALL GET_GOTO(.false.,n-1,dble(nodes),(o+1)/2,go_to,dummy1,dummy2)
                              else
                                 younger=.false.
                              endif
                           endif
                           if (younger) then
                              if ((MOD(o,2)==0).and.(o/2==go_to)) then  !This means that we got the MRCA lineage, so it goes up to the upper boundary of the block
                                 m=m+1
                              elseif (coal_t(nodes(go_to))>=mbound) then !mbound has the upper limit since it was updated before leaving the previous instruction
                                 m=m+1                      !The rationale here is that if our taxa/node goes to someone older than the time band (not the migration time window but the inter-nodes/tip time bands)
                              endif                         !then our node/taxa is candidate for the emigration event since it is a lineage crossing the chosen time band
                           endif
                           if (m==winner) exit
                           o=o+1
                        enddo
                        !Ladies and gentleman the leaving lineages is:
                        if (MOD(o,2)==0) then !Nodes
                           z=z+1
                           emigrants(z,1) = nodes(o/2)
                           CALL RANLUX (RVEC, 1)
                           u=RVEC(1)
                           emigrants(z,2) = u*(mbound-nbound) + nbound
                           emigrants(z,3) = h
                           emigrants(z,4) = list(j,2)
                           emigrants(z,5) = 1   !One means node
                        else                  !Taxa
                           z=z+1
                           emigrants(z,1) = taxa((o+1)/2)
                           CALL RANLUX (RVEC, 1)
                           u=RVEC(1)
                           emigrants(z,2) = u*(mbound-nbound) + nbound
                           emigrants(z,3) = h
                           emigrants(z,4) = list(j,2)
                           emigrants(z,5) = 0   !Zero means taxa
                        endif
                        !And modify everything
                        switch=.true.
                     endif
                     lin=lin+ad
                     if (switch) exit
                  enddo
               enddo
            endif
            if (allocated(list)) then
               deallocate(list)
            endif
         endif
      enddo
      !--------------------------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------------------------
      !Last instruction is for avoiding repeated lineages migrating (if repeated we keep the youngest migration event)
      z=COUNT(emigrants(:,1)/=0)
      if (z>z0) then
        do i=z0+1,z-1,1
          do j=i+1,z,1
            if (emigrants(j,1)/=0) then !this is because further instructions can erase rows so not necessarily they have data
              if ((emigrants(i,1)==emigrants(j,1)).and.(emigrants(i,5)==emigrants(j,5)).and.(emigrants(i,3)==emigrants(j,3))) then !If we find a repetition (both taxa, both nodes; and both share donor block)
                if (emigrants(i,2)>emigrants(j,2)) then !If i-th is older then we put j-th in i and erase j-th
                   emigrants(i,:)=emigrants(j,:)
                endif
                emigrants(j,1)=-1 !we only label the repeated ones with a (-1) because no taxon or node can have this number
              endif
            endif
          enddo
        enddo
        !--
        a=0
        do i=z,z0+1,-1
           if (emigrants(i,1)==-1) then
              do k=i,z-1,1
                 emigrants(k,:)=emigrants(k+1,:)
              enddo
              emigrants(z-a,:)=0
              a=a+1
           endif
        enddo
      endif


      return

   end subroutine Migration




   subroutine Get_mig_coal(mig_lin,go_to,h,n,en,dn,ni,gentime,Ne,growth,events,&
                           &Mages,Mtaxa,MtaxaB,Mnodes,MnodesB,Mcoal_t,time,winner,Ptime,older)

      use irandom
      use routines
      use luxury

      implicit none

      !Externals
      INTEGER, INTENT(IN):: mig_lin                              !Id of the inmigrant lineage
      INTEGER, INTENT(IN):: go_to                                !Id of the parental lineage of the inmigrant, locally (in the M arrays)
      INTEGER, INTENT(IN):: h                                    !Id of the current receiver population/block
      INTEGER, INTENT(IN):: n                                    !Overall sample size
      INTEGER, INTENT(IN):: en                                   !Number of events
      INTEGER, INTENT(IN):: dn                                   !# of demes in the sample
      INTEGER, DIMENSION(dn+en),INTENT(INOUT):: ni               !ni:= stores the sample size of all the blocks/populations
      REAL(8), INTENT(IN):: gentime                              !Generation time
      REAL(8), DIMENSION(dn), INTENT(IN):: Ne                    !Vector of population sizes of each population (deme)
      REAL(8), DIMENSION(dn), INTENT(IN):: growth                !Vector of growth rates for each block of coalescent simulations
      REAL(8), DIMENSION(en,8), INTENT(IN):: events              !Array with the info of ancient demographic events (same as always)
      REAL(8), DIMENSION(dn+en,n), INTENT(INOUT):: Mages         !Master ages contains the ages of all the blocks/populations
      INTEGER, DIMENSION(dn+en,n), INTENT(INOUT):: Mtaxa         !Master taxa:= the taxa of all blocks/populations
      INTEGER, DIMENSION(dn+en,n), INTENT(INOUT):: MtaxaB        !Master taxa:= the taxa of the blocks/populations with local numeration
      INTEGER, DIMENSION(dn+en,n), INTENT(INOUT):: Mnodes        !Master nodes:= the nodes of all blocks/populations with continuous overall numeration (for interblock instructions)
      INTEGER, DIMENSION(dn+en,n), INTENT(INOUT):: MnodesB       !Master nodes contains the nodes of the blocks/populations in local numeration
      REAL(8), DIMENSION(dn+en,n-1), INTENT(INOUT):: Mcoal_t     !Master coal_t contains all the coal_t of all the blocks/populations
      REAL(8), INTENT(INOUT):: time                              !Time to the coalescent of the migrant lineage
      INTEGER, INTENT(OUT):: winner                              !The chosen lineage
      REAL(8), INTENT(OUT):: Ptime                               !Time between the arrival and the coalescent of the migrant lineage
      LOGICAL, INTENT(OUT):: older                               !An indicator to be true when the lineage coalesced outside the receiver block
      !Internals
      INTEGER:: ih                                               !Same as h but for internal use (it can change when lineage migrate not only among blocks but to older blocks)
      INTEGER:: i, j, a, b, c, d, k, r, z, dummy, dummy2
      REAL(8):: e, f
      REAL(8), ALLOCATABLE, DIMENSION(:):: ages                  !Same as always
      INTEGER, ALLOCATABLE, DIMENSION(:):: taxa, taxaB           !Same as always
      INTEGER, ALLOCATABLE, DIMENSION(:):: nodes, nodesB         !Same as always
      REAL(8), ALLOCATABLE, DIMENSION(:):: coal_t                !Same as always
      REAL(8), ALLOCATABLE, DIMENSION(:):: nodes2                !Same as always
      REAL(8):: lambda
      REAL(8), ALLOCATABLE, DIMENSION(:,:):: list                !List of time bins where them igrant lineage can coalesce
      REAL, DIMENSION(1):: RVEC                                  !Random number
      REAL(8):: p1, p2, bound, itime, time1, time2, u
      INTEGER:: Nblock, nlin
      LOGICAL:: coalesced, done, last, leave, switch
      INTEGER:: loc_h, lineage                                   !loc_h:= Location of the block h; lineage:= temporal carrier of winner lineage
      INTEGER:: pos1, pos2                                       !For using MIN/MAXLOC
      LOGICAL:: id                                               !True if the lineage is a taxon, false if it is a node
      !--------------------------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------------------------
      !START
        older=.false.
        ih=h
        Ptime=0
        DO  !Cycle moving around blocks (in case the coalescent doesn't happen in the receiver block but upstream)
           itime=time
           !--------------------------------------
           !This is only to get the upper boundary of the block/population
           switch=.false.
           i=1  !We are gonna screen every event
           DO
              IF ((events(i,4)==ih).or.(events(i,5)==ih)) THEN
                 bound=events(i,1)
                 switch=.true.
              ENDIF
              IF ((i==en).and.(.NOT.switch)) THEN !This is in case no bound is found that means our block is the last one
                 bound=1000000000
                 switch=.true.
              ENDIF
              IF (switch) EXIT
              i=i+1
           ENDDO
           IF (ni(ih)>0) THEN
               !This is only to get the location in events of the block h
               IF (ih>dn) THEN
                  i=1
                  done=.false.
                  DO
                     IF (events(i,8)==ih) THEN
                        loc_h=i
                        done=.true.
                     ENDIF
                     IF (done) EXIT
                     i=i+1
                  ENDDO
               ENDIF
               !--------------------------------------
               !So we prepare everything for getting the coalescent time
               IF (ALLOCATED(taxa)) THEN
                   DEALLOCATE(ages,taxa,nodes,coal_t)
               ENDIF
               !----
               !+++++++
               ALLOCATE(ages(ni(ih)),taxa(ni(ih)),taxaB(ni(ih)),nodes(ni(ih)-1),nodesB(ni(ih)-1),coal_t(ni(ih)-1))
               DO i=1,ni(ih)-1,1
                  ages(i)=Mages(ih,i)
                  taxa(i)=Mtaxa(ih,i)
                  taxaB(i)=MtaxaB(ih,i)
                  nodes(i)=Mnodes(ih,i)
                  nodesB(i)=MnodesB(ih,i)
                  coal_t(i)=Mcoal_t(ih,i)
               ENDDO
               ages(ni(ih))=Mages(ih,ni(ih))
               taxa(ni(ih))=Mtaxa(ih,i)
               taxaB(ni(ih))=MtaxaB(ih,i)
               !------------------------------------------------------------------------
               if (ni(ih)>1) then
                  a = COUNT(coal_t<itime)  !The number of lineages crossing the timeline is actually n minus one for each coalescent
                  b = COUNT(ages>itime)    !Number of taxa older than time
                  nlin = ni(ih)-a-b
               else
                  nlin=1
               endif
               !Now we make the list of lineages changes:
               IF (ALLOCATED(list)) THEN
                  DEALLOCATE(list)
               ENDIF
               c = COUNT((coal_t>itime).and.(coal_t<bound))  !The number of coalescents in the time bin
               d = COUNT((ages>itime).and.(ages<bound))      !The number of ages in the time bin
               ALLOCATE(list(c+d+1,2))       !Number of timebands is the number of events happening (coalescent or new sample) in the timewindow + 1
               DO i=1,c,1                    !List gets in the first column the upper bound of the time band and in the second column the change in number of lineages..
                  list(i,1)=coal_t(a+i)      !...ocurring at that time (+1 is because there is a new ancient sample or -1 is a coalescent)
                  list(i,2)=-1
               ENDDO
               DO i=1,d,1
                  list(c+i,1)=ages(ni(ih)-b+i)
                  list(c+i,2)=+1
               ENDDO                                            !ni(h)-a+b-1 is the number of events (coalescent or ancient taxa appearing) before time
               list(c+d+1,1)=bound
               list(c+d+1,2)=0
               CALL SSORT(list(:,1),list(:,2),c+d+1,2)
               !--------------------------------------
               !We get the coalescent time
               lambda = nlin
               i=1
               time1=itime
               !We get the Ne and the growth rate to adjust times
               IF (ih>dn) THEN
                  e=events(loc_h,2)
                  f=events(loc_h,3)
               ELSE
                  e=Ne(ih)
                  f=growth(ih)
               ENDIF
               !-- The next instruction works as follows: the obtained time is going to be checked against the time bin, if it crosses the time bin upper boundary,
               !..a new time will be drawn (with adjusted number of lineages)
               DO
                  !-----------------So we get the time
                  time2 = random_exponential()*gentime/lambda  !Coalescent happens faster with intensity increased once per lineage present (instead of the pairwise combinations)
                  CALL Ne_adjuster(e,f,time1,time2)
                  !------------------
                  IF (time2<list(i,1)) THEN  !If the time doesnt cross the next bound, we'r done!
                     done=.true.
                     CALL RANLUX (RVEC, 1)
                     u=RVEC(1)
                     winner=NINT(real(lambda)*u+0.5)
                     time = time2
                  ELSE                         !If our time excedds the lenght in our time bin (so it crosses to the next one)
                     IF (i==c+d+1) THEN        !If there are no older events we keep the time
                        done=.true.
                        winner=0
                        time=time2
                     ELSE
                        done=.false.
                        time1=list(i,1)
                        lambda=lambda+list(i,2)
                     ENDIF
                  ENDIF
                  IF (done) EXIT
                  i=i+1
               ENDDO
               i=1
               j=1
               done=.false.
               if (winner>0) then !This is if we got a coalescent inside the time frame
                  if (ni(ih)>1) then
                     CALL Get_lineage(winner,ni(ih),ages,taxaB,taxa,nodesB,nodes,coal_t,time,lineage,id)
                  else
                     lineage=taxa(1)
                  endif
               endif
               !--------------------------------------
               IF (time<bound) THEN
                  coalesced=.true.
                  !--------------------------------------
                  !We get the info to find out if during the time between migration event and the coalescent there was a further migration event
                  Ptime=Ptime+MIN(bound,time)-itime
                  winner=lineage
               ELSE
                  coalesced=.false.
                  time=bound       !Now the time to the immigration is the age of the new block (which will be found in the next section)
               ENDIF
               DEALLOCATE(ages,taxa,taxaB,nodes,nodesB,coal_t)
           ELSE
               coalesced=.false.
               time=bound
           ENDIF
           !******************************************************************************************************************
           ! G L O B A L   R E C E I V E R   A D J U S T M E N T S
           !------------------------------------------------------------------------------------------------------------------
           IF (h/=ih) THEN  !This means that our current block is not the receiver block, it inherited the migrant lineage from a downstream block where it failed to coalesce
               !(A)First with ages-------------------
               pos1 = rMINMAXLOC2(-1,n,Mages(ih,:),itime,DBLE(100000000000000.0),1,n)
               IF (pos1==0) THEN !If no one was older than itime then we just put the age at the end of the array (all the taxa are same or younger)
                  pos1=ni(ih)+1
               ELSE
                  DO i=ni(ih)+1,pos1+1,-1  !If position is intermediate we have to make room for the new one (which is the last at right of the ones younger or = to itime)
                     Mages(ih,i)=Mages(ih,i-1)
                  ENDDO
               ENDIF
               Mages(ih,pos1)=itime
               !(B)Then with taxa and taxaB----------
               DO i=1,ni(ih),1       !First we update the names of the taxaB that are older than bound (they only get their own name +1)
                  IF (MtaxaB(ih,i)>=pos1) THEN
                     MtaxaB(ih,i)=MtaxaB(ih,i)+1
                  ENDIF
               ENDDO
               !Now we move the stuff at right of the receiver lineage to leave one single space for the migrant lineage
               IF ((.NOT.coalesced).OR.(ni(ih)==1))THEN    !First clause means that the lineage just cross the block without coalescing; Second clause means there is only one lineage
                  Mtaxa(ih,ni(ih)+1)=mig_lin
                  MtaxaB(ih,ni(ih)+1)=pos1
               ELSEIF (ni(ih)>1) THEN
                  IF  (((.NOT.id).AND.(Mtaxa(ih,ni(ih))==winner)).OR.((id).AND.(Mnodes(ih,ni(ih)-1)==winner)))  THEN !First clause means the receiver is a taxon & in the last position at right
                     Mtaxa(ih,ni(ih)+1)=mig_lin
                     MtaxaB(ih,ni(ih)+1)=pos1
                  ELSE         !If the lineage got its coalescent (and is not in the extreme) here then:
                     IF (.NOT.id) THEN                                   !If the receiver is a taxon then the right limit of its range is only its position
                        pos2 = iMINMAXLOC1(1,n,Mtaxa(ih,:),winner,winner,1,n)  !We get its position
                        r=pos2
                     ELSE                                                !If the receiver is a node we need to find its range:
                        ALLOCATE(nodes2(ni(ih)-1)) !This is because nodes arrays measure n-1 (we've not increased ni(ih) yet)
                        DO i=1,ni(ih)-1,1
                           nodes2(i)=Mnodes(ih,i)
                        ENDDO
                        pos2 = rMINMAXLOC1(1,ni(ih)-1,nodes2,DBLE(winner),DBLE(winner),1,ni(ih)-1)  !We get its position
                        CALL GET_GOTO(.true.,ni(ih)-1,nodes2,pos2,dummy,dummy2,r)
                        DEALLOCATE(nodes2)
                     ENDIF
                     !Now we move all the stuff at right exactly one place (the difference of the receiver being a taxon or node is only where we stop -r+2-)
                     DO i=ni(ih)+1,r+2,-1
                        Mtaxa(ih,i)=Mtaxa(ih,i-1)
                        MtaxaB(ih,i)=MtaxaB(ih,i-1)
                     ENDDO
                     Mtaxa(ih,r+1)=mig_lin
                     MtaxaB(ih,r+1)=pos1   !We place it at right of the receiver lineage
                  ENDIF
               ENDIF
               !(C)Afterwards with nodes-------------
               IF (coalesced) THEN
                  pos2 = rMINMAXLOC2(-1,n-1,Mcoal_t(ih,:),time,DBLE(100000000000000.0),1,n-1)

                  IF (ALL(Mcoal_t(ih,:)<time)) THEN
                     pos2=ni(ih)
                  ENDIF

               ELSE
                  pos2 = rMINMAXLOC2(-1,n-1,Mcoal_t(ih,:),bound+1,DBLE(100000000000000.0),1,n-1)
                  IF (pos2==0) THEN
                     pos2=ni(ih)
                  ENDIF
               ENDIF
               !IF (ALL(Mcoal_t(ih,:)<time)) THEN
               !   pos2=ni(ih)
               !ENDIF
               DO i=1,ni(ih)-1,1       !First we update the names of the nodesB that are older than bound (they only get their own name +1)
                  IF (MnodesB(ih,i)>=pos2) THEN
                     MnodesB(ih,i)=MnodesB(ih,i)+1
                  ENDIF
               ENDDO
               !First we find the node than changes
               IF (ni(ih)>0) THEN   !In an empty block, nodes will stay empty, so modifications only for non-empty blocks
                  IF ( (.NOT.coalesced) .OR. (ni(ih)==1) ) THEN
                     Mnodes(ih,ni(ih))=go_to
                     MnodesB(ih,ni(ih))=pos2
                  ELSEIF (ni(ih)>1) THEN
                     IF ( ((.NOT.id).AND.(Mtaxa(ih,ni(ih))==winner)) .OR. ((id).AND.(Mnodes(ih,ni(ih)-1)==winner)) )  THEN
                        Mnodes(ih,ni(ih))=go_to
                        MnodesB(ih,ni(ih))=pos2
                     ELSE
                        DO i=ni(ih),r+1,-1   !So we move the stuff at right exactly one place (the difference of the receiver being a taxon or node is only where we stop)
                           Mnodes(ih,i)=Mnodes(ih,i-1)
                           MnodesB(ih,i)=MnodesB(ih,i-1)
                        ENDDO
                        Mnodes(ih,r)=go_to
                        MnodesB(ih,r)=pos2
                     ENDIF
                  ENDIF
               ENDIF
               !..and finally coal_t
               DO i=ni(ih),pos2+1,-1   !So we move the stuff at right exactly one place (the difference of the receiver being a taxon or node is only where we stop)
                  Mcoal_t(ih,i)=Mcoal_t(ih,i-1)
               ENDDO
               IF (coalesced) THEN
                  Mcoal_t(ih,pos2)=time
               ELSE
                  IF (ni(ih)>0) THEN
                     Mcoal_t(ih,pos2)=bound+1
                  ENDIF
               ENDIF
               !---Done
               ni(ih)=ni(ih)+1   !Don't remember if we have now 1 sample more in the block, sample size is +1
           ENDIF
           !******************************************************************************************************************
           !-----------------------------------------------------------------------------------------
           !We have to find the new block accepting our immigrant lineage
           IF (.NOT.coalesced) THEN   !If the immigrant lineage failed to coalesce in the current block we have to find what block(s) take lineages from the current one...
              older=.true.
              !---
              done=.false.       !..and define the new current acceptor block
              last=.false.
              leave=.false.
              i=1
              DO     !This cycle walks all the events searching for the ones that takes our current one as a source of lineages
                 if ((events(i,4)==ih).or.(events(i,5)==ih)) then !If we find a block with our focal block (id=h) then:
                    if (.not.last) then   !Last is a switch used only to indicate if it's the first (=false) or second time (=true) we find a block receiving lineages from the focal one (h)
                       a=NINT(events(i,8))      !The only difference is that if is the first time (last=false) the name and proportion of lineages is registered in "a" and "p1" respectively
                       if (events(i,4)==ih) then
                          p1=events(i,6)
                       elseif (events(i,5)==ih) then
                          p1=events(i,7)
                       endif
                       last=.true.
                    else !-----------------> if is the second time (last=true) we get a block receiving lineages from the focal one (h)
                       b=NINT(events(i,8))         !...then the name and proportion of lineages is registered in "b" and "p2" respectively
                       if (events(i,4)==ih) then
                          p2=events(i,6)
                       elseif (events(i,5)==ih) then
                          p2=events(i,7)
                       endif
                       leave=.true.          !If we entered here, it is the second time and we are done finding blocks receving lineages (there are maximum 2)
                    endif
                    !----------------------------------------> if the block accepting lineages from the focal (h) get 100% of them then we enter the first if below without searching further
                    if ((events(i,4)==events(i,5)).or.((events(i,4)==ih).and.(events(i,6)==1.0))&
                    &.or.((events(i,5)==ih).and.(events(i,7)==1.0))) then
                       NBlock=a
                       leave=.true.   !...and since we are done we activate the exit pass
                    elseif (leave) then   !If we have the exit pass activated (leave) without visiting the previous if-option, then we have to acceptor blocks and have to sort the receiver one
                       CALL RANLUX (RVEC, 1)
                       u=RVEC(1)
                       if (p1/(p1+p2)<u) then !So here we get randomly who get our immigrant lineage
                          NBlock=a
                       else
                          NBlock=b
                       endif
                    endif
                 endif
                 IF (leave) EXIT
                 i=i+1
              ENDDO
           ELSE  !If coalesced we apply this adjustments to remove the new lineages in all upstream (ancestral) blocks
              z=ih
              pos1 = rMINMAXLOC1(1,en,events(:,8),DBLE(ih),DBLE(ih),1,en)      !1. First we get the number of event our block is
              DO i=MAX(pos1+1,1),en,1  ! MAX(ih-dn+1,1) scans every block starting in the next one to our current block (ih+1, but in the events account is ih-dn+1), or the 1st block if ih is a population (ih<=dn)
                 IF ((events(i,4)==z).or.(events(i,5)==z)) then
                    k=NINT(events(i,8))
                    IF (ANY(Mtaxa(k,:)==lineage)) THEN  !Second clause means: if the coalescent time is younger than the age of the block...
                       pos1 = iMINMAXLOC1(1,n,Mtaxa(k,:),lineage,lineage,1,n)          !...otherwise the lineage we want to update -the receiver- will survive until the start of the block and no update would be necessary
                       Mtaxa(k,pos1)=go_to
                    ENDIF
                    z=k
                 ENDIF
              ENDDO
              done=.true.
           ENDIF
           IF (done) EXIT
           ih=NBlock
        ENDDO
        !--------------------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------------------

        
      return


   end subroutine Get_mig_coal




   subroutine Get_lineage(winner,n,ages,taxaB,taxa,nodesB,nodes,coal_t,time,lineage,id)

      use irandom

      implicit none

      !Externals
      integer, intent(in):: winner
      integer, intent(in):: n                             !Sample size (number of lineages coalesced in the focal population/block)
      real(8), dimension(n), intent(in):: ages            !Same as always
      integer, dimension(n), intent(in):: taxaB           !Taxa array with numbers only in the range 1-n
      integer, dimension(n), intent(in):: taxa            !Taxa array with the actual names of the lineages taken as taxa (wich could be nodes names)
      integer, dimension(n-1), intent(in):: nodesB        !Nodes array with numbers in the range 1-(n-1) (internal local -intrablock- numeration)
      integer, dimension(n-1), intent(in):: nodes         !Nodes array with the actual names of the nodes (external -interblock- numeration)
      real(8), dimension(n-1), intent(in):: coal_t        !Master coal_t contains all the coal_t of all the blocks/populations
      real(8), intent(in):: time                          !The time at what the coalescent happens
      integer, intent(out):: lineage                      !Id of the taxa/node that coalesces with the migrant lineage
      logical, intent(out):: id                           !False if the lineage is a taxon; True if it is a node
      !Internals
      integer:: m, o
      integer:: go_to
      logical:: younger
      integer:: dummy1,dummy2
      !-------------------------------------------------------------------------------------------------------------------------
      !-------------------------------------------------------------------------------------------------------------------------
      !START
      m=0
      o=1    !This instruction works by scanning every taxa and node to see if they belong to the subset that cross the time
      do     !..and stopping when we get the chosen one (so it is only to get the identity of the chosen one)
         if (mod(o,2)==0) then !This is for nodes
            if (coal_t(nodesB(o/2))<time) then !This 'if' is an assurance for not considering lineages that start at an older point than the time band
               younger=.true.
               call GET_GOTO(.true.,n-1,dble(nodesB),o/2,go_to,dummy1,dummy2)
            else
               younger=.false.
            endif
         else                  !This is for taxa
            if (ages(taxaB((o+1)/2))<time) then
               younger=.true.
               call GET_GOTO(.false.,n-1,dble(nodesB),(o+1)/2,go_to,dummy1,dummy2)
            else
               younger=.false.
            endif
         endif
         if (younger) then
            if ((mod(o,2)==0).and.(o/2==go_to)) then  !This means that we got the MRCA lineage, so it goes up to the upper boundary of the block
               m=m+1
            elseif (coal_t(nodesB(go_to))>=time) then
               m=m+1                      !The rationale here is that if our taxa/node goes to someone older than the time
            endif
         endif   !then our node/taxa is candidate for the inmigration event since it is a lineage crossing the chosen time band
         if (m==winner) exit
         o=o+1
      enddo
      !Ladies and gentleman the leaving lineages is:
      if (mod(o,2)==0) then !Nodes
         lineage = nodes(o/2)
         id=.true.
      else                  !Taxa
         lineage = taxa((o+1)/2)
         id=.false.
      endif


      return

   end subroutine Get_lineage




   subroutine Get_downstream_lin(n,x,nodes,A1,TN1,A2,TN2)

      implicit none
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------------
      !THE PURPOSE OF THIS PROGRAM IS TO FIND THE TWO LINEAGES FOR WHICH OUR CURRENT NODE (X) IS THE PARENTAL NODE. IN OTHER WORDS THIS ROUTINE FINDS THE TwO
      !DOWNSTREAM INCIDENCES OF THE X NODE
      !-------------------------------------------------------------------------------------------------------------------------------------------------------
      integer, intent(in):: n                          !Size of the array (which is the sample size minus one)
      integer, intent(in):: x                          !The position of the focal node
      real(8), dimension(n), intent(in):: nodes        !Nodes array with their ages instead of names
      integer, intent(out):: A1                        !Position of the left incidence
      integer, intent(out):: TN1                       !The left incidence is a taxa if TRUE and a node if FALSE
      integer, intent(out):: A2                        !Position of the right incidence
      integer, intent(out):: TN2                       !The right incidence is a taxa if TRUE and a node if FALSE

      integer:: i
      logical:: switch
      real(8):: Maximum
      !---------------------------
      !For the left-hand indicence
      if (x==1) then
         A1=1
         TN1=0
      else
         if (nodes(x-1)>nodes(x)) then !If the immediate previous nodeis already older than x it means that the incident lineage comes from the taxon that isl ocated between both nodes (taxon x)
            A1=x
            TN1=0   !It is a taxon
         else
            TN1=1   !It is a node
            i=x
            switch=.false.
            Maximum=0
            do
               i=i-1
               if (nodes(i)>nodes(x)) then   !If the node we reach is older than x we're done
                  switch=.true.
               else
                  if (Maximum<nodes(i)) then  !Since all the nodes falling in this side of the clause are younger than x, if we get someone older than maximum it is the current candidate
                     A1=i                     !So currently A1 gets that position
                     Maximum=nodes(i)         !An Maximum is updated
                  endif
               endif
               if ((switch).or.(i==1)) exit
            enddo
         endif
      endif
      !---------------------------
      !For the right-hand indicence
      if (x==n) then
         A2=n+1
         TN2=0
      else
         if (nodes(x+1)>nodes(x)) then !If the immediate previous nodeis already older than x it means that the incident lineage comes from the taxon that isl ocated between both nodes (taxon x)
            A2=x+1
            TN2=0   !It is a taxon
         else
            TN2=1  !It is a node
            i=x
            switch=.false.
            Maximum=0
            do
               i=i+1
               if (nodes(i)>nodes(x)) then   !If the node we reach is older than x we're done
                  switch=.true.
               else
                  if (Maximum<nodes(i)) then !Since all the nodes falling in this side of the clause are younger than x, if we get someone older than maximum it is the current candidate
                     A2=i                  !So currently A1 gets that position
                     Maximum=nodes(i)      !An Maximum is updated
                  endif
               endif
               if ((switch).or.(i==n)) exit
            enddo
         endif
      endif


      return

   end subroutine Get_downstream_lin




   subroutine Single_Lin_Mig(Top,Bottom,h,gentime,nm,dnen,Mig_labels,MS,TM,MigMat,mig,time,receiver)

        use irandom
        use luxury

        implicit none

        !Externals
        real(8), intent(in):: Top                               !Upper bound of the new lineage (the coalescent time of the migrant lineage)
        real(8), intent(in):: Bottom                            !Lower bound of the new lineage (the time of arrival of the migrant lineage)
        integer, intent(in):: h                                 !Id of the focal population/block
        real(8), intent(in):: gentime                           !Generation time
        integer, intent(in):: nm                                !Number of migration matrix
        integer, intent(in):: dnen                              !Size of the migration matrix (it is actually not the effective size)
        integer, dimension(nm,dnen), intent(in):: mig_labels    !We only consider the columns/rows with values (this are their labels)
        integer, dimension(nm), intent(in):: MS                 !Array with the sizes of those matrix
        real, dimension(nm,2), intent(in):: TM                  !Timespans for the validity of the corresponding matrix
        real(8), dimension(nm,dnen,dnen), intent(in):: MigMat   !Migration matrix
        logical, intent(out):: mig                              !Indicates if migration ocurred
        real(8), intent(out):: time                             !Time to the migration event
        integer, intent(out):: receiver                         !Receiver population
        !Internals
        integer:: i, j, k, m
        real, dimension(1):: RVEC                     !Random number
        real(8):: lo, up, mu, u, length, one          !Boundaries of the migration time window; mu= intensity of migration (lambda of a Poisson density)
        integer:: counter, s, me
        integer, allocatable, dimension(:,:):: list
        real(8):: T
        logical:: done
        !--------------------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------------------
        !START
        allocate(list(SUM(MS),2))
        s=0
        T=Top*2
        mig=.false.
        do i=1,NM,1
           list=0
           lo=TM(i,1)
           up=TM(i,2)
           if ((Top>lo).and.(Bottom<up)) then
              length = MIN(Top,up)-MAX(Bottom,lo)  !First we get the lenght of the new lineage in the interval
              if ((ANY(mig_labels(i,:)==h)).and.(length>0)) then  !If our population appears in the labels, i.e. is considered for migration
                 !This is only to find the row/column in the matrix that corresponds to our population/block
                 do j=1,MS(i),1
                    if (mig_labels(i,j)==h) then
                       k=j
                    endif
                 enddo
                 !Now finding out if emigration ocurrs
                 counter=1
                 do j=1,MS(i),1  !For each column..
                    if (MigMat(i,k,j)/=0) then
                       mu=MigMat(i,k,j)*length/gentime
                       me=random_Poisson(real(mu),.true.)
                       s=s+me !Here s count the number of putative migration events
                       if (me>0) then    !If succesful then
                          list(counter,1)=s                !1st column: number of emigration events
                          list(counter,2)=mig_labels(i,j)  !2nd column: the id of the receiver population/block
                          counter=counter+1
                       endif
                    endif
                 enddo
              endif
              !Getting the time
              m=SUM(list(:,1))  !This is the overall number of putative migrations obtained (for this migration matrix)
              if (m>0) then
                 one=1.0
                 do j=1,m
                    call RANLUX (RVEC, 1)
                    u=RVEC(1)  !In this way we obtain the lowest of their times (of putative migration events)
                    one=MIN(u,one)
                 enddo
                 time = length*one + MAX(Bottom,lo)
                 if (time<T) then !This happens only if the current time is bigger than the previous one
                    T=time
                    !Getting the receiver population
                    call RANLUX (RVEC, 1)
                    u=RVEC(1)
                    j=1
                    one=0
                    done=.false.
                    do
                       one=one+list(j,1)
                       if (u*m<one) then
                          done=.true.
                       endif
                       if (done) exit
                       j=j+1
                    enddo
                    receiver=list(j,2)
                 endif
                 mig=.true.
              endif
           endif
        enddo


        return

   end subroutine Single_Lin_Mig




   subroutine Tree(gentime,Ne,growth,st_point,n,maxL,last,taxa,nodes,coal_times)

      use irandom

      !External variables
      real(8), intent(in):: Ne                                  !Current Ne (the one in the present population/block)
      real(8), intent(in):: growth                              !Growth rate when exponential growth is considered (0=constant Ne)
      real(8), intent(in):: gentime                             !Generation time
      real(8), intent(in):: st_point                            !The absolute time corresponding with the time zero of this sub-tree
      integer, intent(in):: n                                   !Sample size     
      real(8), intent(in):: maxL                                !Maximum limit:= its the ceiling for coalescent times
      logical, intent(in):: last                                !Switch for indicating if the coalescing block is the last one
      
      integer, intent(out), dimension(n):: taxa                 !Taxa store the permutation order of the taxa
      integer, intent(out), dimension(n-1):: nodes              !Nodes store the permutation order of the coalescent times
      real(8), intent(out), dimension(n-1):: coal_times         !Array storing the coalescent times
      !Internal variables
      integer:: i, j
      real:: lambda
      real(8):: t1
      !----------------------------------------------------------------------------------------------
      !Start
      !Getting the inter-node times
      coal_times = 0
      !---Instruction for i=1
      lambda = REAL((n)*(n-1))/2
      coal_times(1) = random_exponential()*gentime/lambda
      call Ne_adjuster(DBLE(Ne),growth,st_point,coal_times(1))
      !---Now for the rest
      if (n>1) then
         do i = 2, n-1, 1
            if ((.NOT.last).and.(coal_times(i-1)>maxL)) then !This instruction avoid getting new times if they are going to be shaved anyway
               do j=i, n-1                                                 !so the times are assigned in an increasing order
                  coal_times(j)=coal_times(i-1)+j
               enddo
               go to 100
            else
               lambda = REAL((n-i+1)*(n-i))/2
               coal_times(i) = random_exponential()*gentime/lambda
               t1=coal_times(i-1)
               call Ne_adjuster(DBLE(Ne),growth,t1,coal_times(i))
            endif
         enddo
      endif
      !Permutation of the taxa labels
100   call random_order(taxa, n)
      !Permutation of the coalescent times
      call random_order(nodes,n-1)
      !-----------------------------------
      
      
      return

   end subroutine Tree




   subroutine AncTree(gentime,Ne,growth,n,st,samplinfo,taxa,nodes,coal_t,maxL,last)

      use irandom

      !External variables
      real(8), intent(in):: gentime                            !Generation time
      real(8), intent(in):: Ne                              !Current Ne (the one in the present population/block)
      real(8), intent(in):: growth                          !growth rate when exponential growth is considered (0=constant Ne)
      integer, intent(in):: n                               !Overall sample size
      integer, intent(in):: st                              !Number of temporal subsamples
      real(8), dimension(3,st), intent(in):: samplinfo      !Array with the 1) sizes, 2) deme identities and 3) ages of temporal subsamples
      integer, intent(inout), dimension(n):: taxa           !Taxa store the permutation order of the taxa, order only have the integers in order for permutation purposes
      integer, intent(inout), dimension(n-1):: nodes        !Nodes store the permutation order of the coalescent times
      real(8), intent(out), dimension(n-1):: coal_t         !Array storing the coalescent times
      real(8), intent(in):: maxL                            !Maximum limit:= its the ceiling for coalescent times
      logical, intent(in):: last                            !Switch for indicating if the coalescing block is the last one
      !Internal variables
      integer:: k, i, j !dummy                              !K (master K:= indicator of usable lineages when a coalescent time is being obtained
      logical:: switch                                      !Switch for enabling/disabling the shaving of nodes overpasing the next subpopulation age
      real(8):: irange, bound                               !Range is the range of taxa available for coalescing; bound:= boundary for shaving nodes
      real(8), dimension(3,st):: isamplinfo                 !Samplinfo for internal use (can be modified)
      integer, allocatable, dimension(:):: taxita           !Same as taxa but for getting a partial tree
      integer, allocatable, dimension(:):: noditos          !Same as nodes but for a partial tree
      real(8), allocatable, dimension(:):: coal_tito        !Array storing the coalescent times of a subsample
      !---------------------------------------------------------------------------------------------------------------------------------------------------
      !Start
      !---------------------------------------------------------------------------------------------------------------------------------------------------
      isamplinfo = samplinfo
      !Initializing things
      coal_t = 0
      !------------------------Getting the coalescent tree of the first temporal subsample (t=0)
      k = NINT(isamplinfo(1,1))
      irange = isamplinfo(1,1)
      allocate(taxita(k),noditos(k-1),coal_tito(k-1))
      if (k>1) then
         call Tree(gentime,Ne,growth,DBLE(0.0),k,maxL,last,taxita,noditos,coal_tito)
         !IMPORTANT: Here the starting point is zero since we are starting a block
         coal_tito = coal_tito + samplinfo(3,1) !...and so we need to add samplinfo(3,1)
         !------------------------and putting the tree in the big arrays
         do i=1,NINT(irange-1),1
            taxa(i)=taxita(i)
            nodes(i)=noditos(i)
            coal_t(i)=coal_tito(i)
         enddo
         taxa(NINT(irange))=taxita(NINT(irange))
         !--------------------------
      endif
      !-----------------
      if (st>1) then
           k=1
           do j=1,NINT(irange-1),1
              if ( coal_t(nodes(j)) > isamplinfo(3,2) ) then !When the coalescent time is older than the age of the 2nd subsample
                 nodes(j) = 0
                 k = k + 1
              endif
           enddo
           deallocate(taxita,noditos,coal_tito)
           !----------------------------Here goes the main part-----------------------------------------------------------------
           do i=2,st,1
              !----This switch is for disabling/enabling the shaving of nodes overpassing next supopulation age
              switch=.false.
              if (i<st) then
                 switch=.true.    !False if we are using last (older) temporal subsample in this lap
                 bound=isamplinfo(3,i+1)
              endif
              !--------------------------
              k=k+INT(isamplinfo(1,i))
              !--------------------------------------------------------------------------------------------------
              call partial_tree(gentime,Ne,growth,n,k,irange,samplinfo(3,1),isamplinfo(1,i),&
              &isamplinfo(3,i),bound,taxa,nodes,coal_t,switch,maxL,last)
           enddo
           !--------------------------------------------------------------------------------------------------------------------
      endif

      
      return

   end subroutine AncTree




   subroutine partial_tree(gentime,Ne,growth,n,k,irange,sinfo0i,sinfo1i,sinfo3i,bound,taxa,nodes,coal_t,switch,maxL,last)

      USE irandom

      !External variables
      real(8), intent(in):: gentime                         !Generation time
      real(8), intent(in):: Ne                              !Current Ne (the one in the present population/block)
      real(8), intent(in):: growth                          !growth rate when exponential growth is considered (0=constant Ne)
      integer, intent(in):: n                               !Overall sample size
      integer, intent(inout):: k                            !Number of entring lineages (updated with the new subsample taxa)
      real(8), intent(inout):: irange                       !Range of taxa available for coalescing
      real(8), intent(in):: sinfo0i                         !Age of the block (sinfo0i - sinfo3i) is the pertinent starting point for adjusting  coal time when variable Ne
      real(8), intent(in):: sinfo1i                         !Number of lineages in the current temporal subsample
      real(8), intent(in):: sinfo3i                         !Age of the current temporal subsample
      real(8), intent(in):: bound                           !Boundary above which nodes should be cut off
      integer, intent(inout), dimension(n):: taxa           !Taxa store the permutation order of the taxa, order only have the integers in order for permutation purposes
      integer, intent(inout), dimension(n-1):: nodes        !Nodes store the permutation order of the coalescent times
      real(8), intent(inout), dimension(n-1):: coal_t       !Array storing the coalescent times
      logical, intent(in):: switch                          !Switch enabling/disabling the shaving of nodes older than next subsample age (or employed boundary)
      real(8), intent(in):: maxL                            !Maximum limit:= its the ceiling for coalescent times
      logical, intent(in):: last                            !Switch for indicating if the coalescing block is the last one
      !Internal variables
      integer:: i, j, g, x                                  !Counters
      integer, dimension(k):: taxita                        !Same as taxa but for getting a partial tree
      integer, dimension(k-1):: noditos                     !Same as nodes but for a partial tree
      real(8), dimension(k-1):: coal_tito                   !Array storing the coalescent times of a subsample
      integer, dimension(k):: gr_size, Up, Lo               !Arrays storing the size and limits of the groups (regresenting lineages)
      integer, dimension(n):: taxaB                         !Array for permutating groups in taxa
      integer, dimension(n-1):: nodesB
      !---------------------------------------------------------------------------------------------------------------------------
      !Start
      if (k-sinfo1i==irange) then      !If there are as many lineages to coalesce as taxa (i.e. all the times surpased the age of next subsample)
         irange=irange+sinfo1i !...so only a permutation is needed
         !---------------------------------------------------
         call Tree(gentime,Ne,growth,sinfo3i-sinfo0i,k,maxL,last,taxita,noditos,coal_tito)
         !Notice that sinfo3i is the starting point (but the point t=0 is at sinfo0i) for adjusting times when variable Ne
         coal_tito = coal_tito + sinfo0i
         !------------------------------------------
         do i=1,NINT(irange),1
            taxaB(i)=taxa(taxita(i))
         enddo
         do i=1,NINT(irange-1),1
            taxa(i)=taxaB(i)
            nodes(i)=noditos(i)
         enddo
         taxa(NINT(irange))=taxaB(NINT(irange))
      else
         !1) Oganizing the info about groupings without new samples
         !1a.)get the info of the groups
         i=0  ! Counter of cycles (cells)
         g=0  ! counter of groups
         do
            i=i+1          !Counter of cycles (of cells)
            if (i==irange) then
               g=g+1  !Counter of groups gotten
               Up(g)=i    !Upper limit of that group is i
            elseif (nodes(i)==0) then
               g=g+1  !Counter of groups gotten
               Up(g)=i    !Upper limit of that group is i
            endif
            if (i==irange) exit
         enddo
         gr_size(1)=Up(1) !The size of first group is Up(1)
         Lo(1)=1          !...and the lower limit of group 1 is 1
         do i=2,g,1
            gr_size(i)=Up(i)-Up(i-1)  !Asigning lower limits and sizes of the groups
            Lo(i)=Up(i-1)+1
         enddo
         !----------------Now defining the new taxa as groups of 1-element
         irange=irange+sinfo1i !the number of of cells in taxa already used (permutated)
         if (sinfo1i>=1) then
            do i=1,NINT(sinfo1i),1
               gr_size(NINT(k-sinfo1i+i))=1
               Up(k-NINT(sinfo1i)+i)=NINT(irange-sinfo1i+i)
               Lo(k-NINT(sinfo1i)+i)=NINT(irange-sinfo1i+i)
               g=g+1
            enddo
         endif
         !2) Getting the (sub)coalescent-tree of the k lineages
         !---------------------------------------------------
         call Tree(gentime,Ne,growth,sinfo3i-sinfo0i,k,maxL,last,taxita,noditos,coal_tito)
         !---------------------------------------------------
         !Notice that sinfo3i is the starting point (but the point t=0 is at sinfo0i) for adjusting times when variable Ne
         coal_tito=coal_tito+sinfo0i
         !3) Applying the obtained permutation of lineages (groups) to the actual groups on taxa
          !3a) Permutating groups (bringing to taxaB the groups in the order defined by taxita, and bringing the respective nodes with us)
         nodesB=0    ! This complex instruction works as follows:
         x=1   !x counts each cell being ocupied
         do i=1,g,1               !g is actually equal to k = # of permuted groups/lineages
            do j=Lo(taxita(i)),Up(taxita(i))-1,1              !(1) j moves along the groups that are going to fill taxaB (their position is given by taxita)
               if (gr_size(taxita(i))>1) then                  !(2) So if we are in the j-th group taken, its label is given by taxita (so the # of group being moved is taxita(j))
                  nodesB(x)=nodes(j)                               ! and Lo(taxita(j)) and Up(taxita(j)) provide the boundaries of the group being moved
               endif                                !We enter the if above if we have to move also the node (but the number in the corresponding node...
               taxaB(x)=taxa(j)                     !...is the consecutive order of coal times being used by being fixed)
               x=x+1
            enddo
            taxaB(x)=taxa(Up(taxita(i)))
            x=x+1
         enddo
         !3b) Setting the nodes--------------------------------------------------------------------------------------------------------
         noditos=noditos+(NINT(irange)-k)   ! range-k = # of fixed nodes, here we adjust the numbers in nodes for continuing the order
         j=0                                ! of the already fixed nodes
         do i=1,NINT(irange-1),1            !So new permuted nodes (and coalescent times) goes only into the empty spaces (0)
            if (nodesB(i)==0.0) then
               j=j+1
               nodesB(i)=noditos(j)
            endif
         enddo
         !---------------
         !returning to taxa from taxaB
         do i=1,NINT(irange-1),1
            taxa(i)=taxaB(i)
            nodes(i)=nodesB(i)
         enddo
         taxa(NINT(irange))=taxaB(NINT(irange))
         !--------------------------------------------
      endif
      !4a) Asigning coalescent times from coal_tito and anulating the ones that are older than the next temporal subsample
      do i=1,k-1,1
         coal_t(NINT(irange-k+i))=coal_tito(i)
      enddo
      k=1
      if (switch) then   ! Enabling/disabling the shaving of nodes older than next temporal subsample age (or bound)
         do i=1,NINT(irange-1),1
            if (coal_t(nodes(i))>bound) then
               nodes(i)=0
               k=k+1                     !k for the next cycle is the number of zeros in nodes + 1
            endif
         enddo
      endif
      !---------------------------------------------------------------------------------------------------------------------------

      
      return

   end subroutine partial_tree




   subroutine order_tree(n,taxa,ages,nodes,coal_t,MasterTree)

      use routines

      implicit none

      integer, intent(in):: n                                   !Overall sample size
      integer, intent(in), dimension(n):: taxa                  !Taxa store the permutation order of the taxa, order only have the integers in order for permutation purposes
      real(8), intent(in), dimension(n):: ages                  !Ages of each taxa (individual) in the sample (ordered)
      integer, intent(in), dimension(n-1):: nodes               !Nodes store the permutation order of the coalescent times
      real(8), intent(in), dimension(n-1):: coal_t              !Array storing the coalescent times
      real(8), intent(out), dimension(3,2*n-1):: MasterTree     !The tree info ready for being used by the graphicer

      integer:: i
      real(8), dimension(2*n-1):: Mtree, branch_len             !An array with ages in the odd cells and coal_times in the even cells, and branch_len := branch lenghts
      integer, dimension(2*n-1):: go_to                         !Array containing the information about of where the current node/taxa goes in the tree (master_tree array)
      integer, dimension(n-1):: nodeB, order
      integer:: A, B, UpB, LoB
      !-----------------------------------------------------------------------------------------------------------------------------------------------------
      !Start
      !Gettin the values in "go_to"
      !---------------------------------------------------------------------------------
      do i = 1, n-1, 1
         Mtree(2*i)=coal_t(nodes(i))
         Mtree(2*i-1)=ages(taxa(i))
      enddo
      Mtree((2*n)-1)=ages(taxa(n))
      ! 1st: for the taxa       We send the current taxa to the closest node(at most two nodes for choosing: the neighbors)
      go_to(1) = 2             !First one has no choice: goes to the first node
      do i = 3, (2*n)-3, 2
         if (Mtree(i+1)>Mtree(i-1)) then   !So we compare ages (coalescent times) of the nodes and go for the youngest
            go_to(i) = i - 1
         elseif (Mtree(i+1)<Mtree(i-1)) then
            go_to(i) = i + 1
         endif
      enddo
      go_to((2*n)-1) = (2*n)-2 !Last one has no choice: goes to the last node
      ! 2nd: for the nodes
      nodeB = nodes     ! For nodes it is necesary to use the "backwards window"
      do i = 1, n-1, 1
        order(i) = i    !The procedure consists in:
      enddo               !1. Choose one by one the nodes in increasing order of their age (coalescent time)
      do i=1,n-2,1        !2. Using an array (order) to track the usable nodes (by erasing from the array the already used -younger- ones)
                          !3. Use MAXLOC/MINLOC to get the closest usable node from the remaining in the array "order" but we use "nodeB" instead "nodes"
         B = iMINMAXLOC1(-1,n-1,nodeB,1,n,1,n-1)
         nodeB(B) = 0
         if (B>1) then
            UpB = iMINMAXLOC1(1,n-1,order,1,n,1,B-1)
         else
            UpB=0
         endif
         if (B<n-1) then
            LoB = iMINMAXLOC1(-1,n-1,order,1,n,B+1,n-1)
         else
            LoB=0
         endif
         order(B)=0
         if (UpB==0) then                 !All this is for preventing:   1) Trying to get a cell of a node in the location "0" (which happens if the node is in the edge)
            go_to(2*B) = 2*LoB                                        !  2) Trying to use a node whith a value is "0" (which happens if all the nodes at one side of our node were already used)
         elseif (LoB==0) then
            go_to(2*B) = 2*UpB
         else
            if (nodeB(UpB)==0) then                           !If UpB or LoB gives zero means that every node at one side are not usable, and necessarily the other side has one,
               go_to(2*B) = 2*LoB                             !..so the location given by the other bound (UpB or LoB) is the one.
            elseif (nodeB(LoB)==0) then
               go_to(2*B) = 2*UpB
            else
               if ( nodes(LoB)>nodes(UpB) ) then              !If there are usable nodes in both sides we keep the youngest one
                  go_to(2*B) = 2*UpB
               elseif ( nodes(LoB)<nodes(UpB) ) then
                  go_to(2*B) = 2*LoB
               endif
            endif
         endif
      enddo
      A = iMINMAXLOC1(1,n-1,nodes,1,n-1,1,n-1)   !This is for the oldest node (the root of the tree, or the MRCA of  everyone)
      go_to(2*A) = 2*A    !it send that node to itself
      !---------------------------------------------------------------------------------
      !Finally, we assemble MasterTree
      do i=1,2*n-1,1
         branch_len(i) = Mtree(go_to(i)) - Mtree(i)
      enddo
      if (MINVAL(branch_len)<0.0) then
         write(*,*) 'Negative branch detected'
         stop
      endif
      do i=1,2*n-1,1
         MasterTree(1,i)=Mtree(i)
         MasterTree(2,i)=go_to(i)
         MasterTree(3,i)=branch_len(i)
      enddo


      return

   end subroutine order_tree




   subroutine Ne_adjuster(Ne,growth,t1,t2)

      implicit none

      !External variables              !TIME IS ALWAYS BACKWARDS!!
      real(8), intent(in):: Ne         !Ne in the actual block (at the starting time of the block)
      real(8), intent(in):: growth     !Constant (k) of an exponential growth: Nt = N0 e**(kt)
      real(8), intent(in):: t1         !The time where the target time (t2) time departs from
      real(8), intent(inout):: t2      !The coalescent time we want to adjust
      !Internal variables
      real(8):: Nt
      !-------------------------------------------------------------------
      if (growth<0) then
         Nt=Ne*EXP(growth*t1)
         t2=LOG(t2*2*Nt*ABS(growth)+1)/ABS(growth)
         t2=t1+t2
      elseif (growth>0) then
         Nt=Ne*EXP(growth*t1)
         t2=ABS(LOG(ABS((-1*t2*2*Nt*growth)+1))/growth)
         t2=t1+t2
      elseif (growth==0) then
         t2=t1+t2*2*DBLE(Ne)
      endif


      return

   end subroutine Ne_adjuster




   subroutine GET_GOTO(mode,n,vector,x,go_to,lo,up)
   !This function is aimed to find the lowest of two neighbors the right and the left closest neighbors
   !of a defined point (x) in an array (vector) with the condition that they have to be the most immediate neighbor having
   !larger numbers than the one in position x. The output is a number (go_to) with the position of the most immediate
   !neighbor with a value higher than the reference one. It serves to find the node where the reference node (x) goes-to
   !in a tree with format taxa/nodes/coalescent_times

      implicit none

      LOGICAL, INTENT(IN):: mode                   !False if it is a taxa, true if it is a node
      INTEGER, INTENT(IN):: n                      !Size of the array
      REAL(8), DIMENSION(n), INTENT(IN):: vector   !Array
      INTEGER, INTENT(IN):: x                      !Reference position
      INTEGER, INTENT(OUT):: go_to                 !Output
      INTEGER, INTENT(OUT):: lo, up                !Range of the node

      INTEGER:: i,j
      REAL(8):: a,b
      !----
      if (mode) then
         if (vector(x)==MAXVAL(vector)) then
            go_to=x
            lo=1
            up=n+1
         else
            !Right-hand neighbor
            if (x<n) then
               i=x
               do
                  i=i+1
                  if ((vector(i)>vector(x)).or.(i==n)) exit
               enddo
               if (vector(i)>vector(x)) then  !If the exit of the do-cycles is because we found an older node then:
                  a=vector(i)                 !'a' gets the age of the older node (for later comparison to decide the go_to)
                  up=i                        !..and the upper (right) limit is 'i'
               else                       !Only way to not exiting the cycle by finding an older node is because we get the last node so:
                  a=MAXVAL(vector)*1.1    !'a' gets an age older than anyone so that the go_to be assigned to the other node (the left one)
                  up=n+1                  !..and the right limit is the end (one place at right of the node)
               endif
            else                      !If (x==n) our node is the last one,..
               a=MAXVAL(vector)*1.1   !..'a' gets an age older than anyone to assure that the go_to is the other (the left) incidence
               up=n+1                 !..and the right limit is the last taxon (n+1)
            endif
            !Left-hand neighbor
            if (x>1) then
               j=x
               do
                  j=j-1
                  if ((vector(j)>vector(x)).or.(j==1)) exit
               enddo
               if (vector(j)>vector(x)) then   !If the exit of the do-cycles is because we found an older node then:
                  b=vector(j)                  !'b' gets its age
                  lo=j+1                       !..and the lower (left) limit is thetaxon at right of that node
               else                         !If the exit was because we got the last node at left then:
                  b=MAXVAL(vector)*1.1      !'b' gets an age older than anyone so that the go_to be equal to the other node (the left)
                  lo=1                      !..and the left limit is the first taxon
               endif
            else                         !If x is the firts node:
               b=MAXVAL(vector)*1.1     !'b' gets an age older than anyone so that the go_to be necessarily the other (the right) node
               lo=1                     !..and the left limit is the first taxon
            endif
            !Now the goto
            if (a<b) then
               go_to=i
            elseif (a>b) then
               go_to=j
            endif
         endif
      else
         if (x==1) then
            go_to=1
         elseif (x==n+1) then
            go_to=n
         else
            if (vector(x-1)>vector(x)) then
               go_to=x
            elseif (vector(x-1)<vector(x)) then
               go_to=x-1
            endif
         ENDIF
      endif


      return

   end subroutine GET_GOTO




   subroutine mutate(n,mutrate,MasterTree,mut_Mtree)

      use irandom
    
      implicit none
      !Variables
      !Externals
      integer, intent(in):: n
      real(16), intent(in):: mutrate
      real(8), intent(in), dimension(3,(2*n)-1):: MasterTree
      real(8), intent(out), dimension(4,(2*n)-1):: mut_Mtree
      !Internals
      integer:: i
      real(16):: mu
      !----------------------------------------------------------------------
      !Start
      do i = 1, (2*n)-1, 1
         mu = MasterTree(3,i)*mutrate
         mut_Mtree(4,i) = random_Poisson(real(mu),.true.)
      enddo

      mut_Mtree(1,:)=MasterTree(1,:)
      mut_Mtree(2,:)=MasterTree(2,:)
      mut_Mtree(3,:)=MasterTree(3,:)
    
    
      return

   end subroutine mutate










end module treeator
