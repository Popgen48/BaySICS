!****************************************************************************
!  PROGRAM: Coalescent version 7.6.6 - Console - Intel
!
!  PURPOSE:  It performs coalescent simulations of single or multiple loci, 
!  following a custom design that allows an arbitrary demographic history
!  that can be both stepwise or with exponential growth/decrease. It also
!  allows arbitrary population structures in time with both split and 
!  merging of populations and migration among any teminal or ancestral 
!  populations. The program allows full integration of ancient DNA or 
!  other type of heterochronous sampling, and can perform an automatic
!  calibration of radiocarbon dates. 
!  For the simulations and the empirical data (a DNA alignment in 
!  FASTA format), this program can compute and save: 
!   1. Summary statistics (simulations and empirical alignment)
!   2. Coalescent times (only for simulations)    
!   3. FASTA files (only for simulations)
!   4. Arlequin files (both simulations and empirical alignment)
!   5. Neighbor-Joining trees (both simulations and empirical alignment)
! 
!  For multilocus simulation different DNA fragments can be simulated 
!  independently (including their own mutation rate and marker type 
!  specifications. In this case (multilocus), the summary statistics can 
!  be computed in two ways: making the average over loci, or reporting all  
!  summary statistics separatedly. For multilocus simulation, the 
!  Neighbor-Joining trees can be inferred and recorded for every individual	
!  locus or alltogether from distances pooled over all loci.
!   	
!  This program was originally designed for ABC analyses but can be used for 
!  other purposes.
!  	
!****************************************************************************

Program Coalescent766
       !!GUI stuff**********************************************
       !use ifport
       !use ifqwin
       !use graphs   
       !!GUI stuff**********************************************

      use treeator
      use routines
      use luxury
        
      implicit none

      !Variables-----------------------
      
       !!GUI stuff**********************************************
       !integer(4) ret
       !logical:: retl
       !!---Custom options
       !integer:: draftsize
       !integer:: graphsize
       !real:: treesize                                    !Size of the coalescent graph
       !integer:: gtype                                    !Type of coalescent graph: 1.dendrogram; 2.cladogram       
       !real, allocatable, dimension(:,:):: demecolor      !For color assignment to demes
       !!GUI stuff**********************************************
       
      !---General
      integer:: i, j, k, l, m, z                          !Counters
      integer:: a, b, ierror                              !Multipurpose
      integer:: n, dn, en, gr, g, maxg                    !Overall sample size, demes nr, events nr, sampling groups nr, st groups nr, max row size of a NJ tree 
      integer:: NM, r, lociNr                             !Number of migration matrixes, number of priors, number of loci        
      character(30):: nameA, nameB, nameC, nameD          !Input files names
      character(30), ALLOCATABLE, DIMENSION(:):: nameL    !Names of the loci       
      integer:: s_int, l_int                              !Interval of simulations/loci for displaying nr of simulation/locus 
      integer:: w_int                                     !Interval for writing results
      integer:: simul_nr                                  !Number of simulations  
      logical:: radio                                     !Radiocarbon dates calibration choice  
      logical:: appendSim                                 !To append/not the simulations to an existing table
      logical, dimension(5):: dofiles                     !Options for: (1)fasta, (2)arlequin, (3)coal-times, (4)summary statistics, & (5)Neighbor-Joining files      
      integer, dimension(2):: dtype                       !Type of Neighbot-Joining distance (cell 1: model (p-dist, Jukes-Cantor...), cell 2 (average, separate trees)
      integer, dimension(3):: datatype                    !Cell 1: Type of data (DNA single/DNA multi/SNPs), cell 2: average/independent (SuS calculation)
      logical:: empal                                     !True if empirical alignment (in fasta) is provided
      logical:: mimic                                     !True to get the simulated alignments mimic the '?'/'n' pattern of the empirical one
      logical, dimension(12):: SuS                        !Array with the info regarding whose summary statistics are being calculated and whose are not
      integer:: seed1      
      !---Intersubprograms 
      real(8), allocatable, dimension(:,:):: samplinfo, isamplinfo       !Sampling info: [1st col.] sample group sizes, [2nd col] deme number, [3rd col] group ages 
      real(8), allocatable, dimension(:):: Ne, iNe                       !Vectors of population sizes 
      real(8), allocatable, dimension(:):: growth, igrowth               !Vectors of growth rates 
      real(8), allocatable, dimension(:):: bid                           !Identifier of population blocks for Ne/growth
      real(8), allocatable, dimension(:,:):: events, ievents             !Info about the internal blocks
      real(8), allocatable, dimension(:,:,:):: MigMat, iMigMat           !Matrix of migration rates 
      real, allocatable, dimension(:,:):: TM                             !Intervals of time for application of the migration matrix
      integer, allocatable, dimension(:):: MS                            !Sizes of migration matrix 
      integer, allocatable, dimension(:,:):: mlbls                       !Id's of the relevant columns in the migration matrix 
      character(30), allocatable, dimension(:):: labels, ilabels         !Labels of organisms (sequences)
      integer, allocatable, dimension(:,:):: stgroups, istgroups         !Statistical groups (for SuS computation)
      integer, allocatable, dimension(:):: gsizes                        !Sizes of the statistical gorups
      real, allocatable, dimension(:,:):: demesranges                    !Range that each deme occupies in taxa       
      real(8):: gen, igen                                                !Generation time 
      real(8):: sexratio, isexratio, adjust                              !Sexratio, in male-female ratios
      integer:: ploidy                                                   !Ploidy:= (1) Haploid; (2) Diploid; (3) Haplo-Diploid 
      
      integer, allocatable, dimension(:):: nb                            !Number of nucleotides
      integer, allocatable, dimension(:):: marker                        !Marker type:= (1)Autosomal; (2)X-linked; (3)Y-linked; (4)Mitochondrial 
      real(16), allocatable, dimension(:):: mutrate, imutrate            !Mutation rate 
      real(16), allocatable, dimension(:,:,:):: sM, isM                  !Substitution matrix (rows/columns: A, C, G, T)
      real(16), allocatable, dimension(:,:):: ACGT, iACGT                !Parameters Pi(A), Pi(C), Pi(G), Pi(T)
      real(16), allocatable, dimension(:):: gamma1, igamma1              !Gamma parameter  
      real(16), allocatable, dimension(:):: invariant, iinvariant        !Proportion of invariant sites      
      integer, allocatable, dimension(:):: gcat                          !Number of gamma categories
      real(8), allocatable, dimension(:,:):: gammacat                    !The probabilities of the discrete gamma categories (all loci)
      real(8), allocatable, dimension(:):: igammacat                     !The probabilities of the discrete gamma categories (single locus)
      
      real(8), allocatable, dimension(:,:):: PrInfo, iPrInfo             !Array carrying the information about priors
      real(8), allocatable, dimension(:):: Priors                        !Array with the random values simulated from the priors          
      integer, allocatable, dimension(:):: taxa                          !Taxa store the permutation order of OTUs (order has the sorted integers for permutation)  
      real(8), allocatable, dimension(:):: ages                          !Ages of each taxa (individual) in the sample (ordered)
      integer, allocatable, dimension(:):: nodes                         !Nodes store the permutation order of the coalescent times 
      real(8), allocatable, dimension(:):: coal_t                        !Array storing the coalescent times 
      integer, allocatable, dimension(:,:):: NodRanges                   !Array with the range of each of the final nodes of the tree 
      real(8), allocatable, dimension(:,:):: MasterTree, mut_Mtree       !The tree info ready for mutation introduction, and the graphicer
      logical:: first                                                    !To control writing files de novo or appending 
      !Summary statistics and ABC stuff
      integer, allocatable, dimension(:,:):: eAlign, sAlign              !Empirical alignment and simulated alignment respectively   
      integer, allocatable, dimension(:,:,:):: MASTEReAlign              !To store all the empirical alignments 
      real(8), allocatable, dimension(:,:):: NJts                        !Neighbor-Joining tree(s), stored in g rows each with a  single tree 
      real(8), allocatable, dimension(:,:,:):: obsNJts                   !Neighbor-Joining tree(s), stored in g rows each with a  single tree 
      integer:: SuStN, SuStNr                                            !Single-locus number of summary statistics, overall number of summary statistics  
      character(12):: char1                                              !Its a recipient of numbers       
      character(17), allocatable, dimension(:):: headings, iheadings     !Headings of the results file (all)
      character(17), allocatable, dimension(:):: prheadings              !Headings of the results file (priors/parameters)
      real(8), allocatable, dimension(:):: irow, row, obsrow             !A single row of the reference table, subrows to carry the SuSt
      real(8), allocatable, dimension(:,:):: table                       !The reference table and the table with coalescent times 
      real(8), allocatable, dimension(:,:,:):: tableC                    !The reference table and the table with coalescent times 
      real(8), allocatable, dimension(:,:,:):: NJtable                   !The file of NJ trees 
      integer, allocatable, dimension(:,:):: dij, Sij, trij              !Distances matrixes
	  !Parallel
	  real:: t1, t2
      integer:: ithread
      character(5):: char 
      !----------------------------------------------------------------------------------------------------------------------
      !   S    T    A    R    T
      !----------------------------------------------------------------------------------------------------------------------
      call cpu_time(t1)
	  
       !!GUI stuff**********************************************
       !ret=DELETEMENUQQ(2,0)
       !ret=DELETEMENUQQ(2,0)
       !ret=DELETEMENUQQ(3,0) 
       !ret=DELETEMENUQQ(3,0)
       !retl=INSERTMENUQQ(3,0,$MENUENABLED, 'About'C, WINABOUT)
       !call inimainwindow(nameA)         
       !draftsize=300                ! Size of the boxes showing the sampling draft and the design of the simulations
       !graphsize=400                ! Size of the coalescent tree graph (not of the tree itself)
       !!GUI stuff**********************************************
      
	  write(*,*) 'Enter thread number'      
      read(*,*) ithread
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(UNIT=1, FILE='iparallel.txt', STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file
      openif: if (ierror==0) then                                              ! ierror=0 means Open was succesfull
         readloop: do i=1,ithread,1 
            read(1,*,IOSTAT=ierror) nameA, seed1 
            if (ierror/=0) exit   !EXIT if not valid 
         enddo readloop            
      else
         write(*,*) 'The seeds file could not be opened :('                
      endif openif
      close(UNIT=1)
      write(*,*) 'Seed in run ', ithread, ' is ', seed1
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  call read_info(radio,s_int,l_int,w_int,simul_nr,appendSim,dofiles,dtype,empal,mimic,SuS,datatype)
      !call read_info(nameA,radio,gtype,treesize,s_int,w_int,simul_nr,appendSim,dofiles,dtype,empal,mimic,SuS,datatype,seed1) 
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      nameB=TRIM(nameA)//'.dat'
      call preread_input(nameB,gr,dn,en,n,r,NM,lociNr)
      !-----------------------------------------------------------------------------------------------------------------------  
      allocate(nameL(lociNr),samplinfo(3,gr),isamplinfo(3,gr),Ne(dn),iNe(dn))
      allocate(growth(dn),igrowth(dn),bid(dn),events(en,8),ievents(en,8))
      allocate(MigMat(NM,dn+en,dn+en),iMigMat(NM,dn+en,dn+en),MS(NM),TM(NM,2))
      allocate(mlbls(NM,dn+en),labels(n),ilabels(n),stgroups(2,n),istgroups(2,n))
       !!GUI stuff***
       !allocate(demesranges(dn+en,2),demecolor(dn+en,3))  !demecolor only if graphic stuff
       !!GUI stuff***      
      allocate(PrInfo(r,9),iPrInfo(r,9),Priors(r),prheadings(r))       
      allocate(nb(lociNr),marker(lociNr),mutrate(lociNr),imutrate(lociNr))      
      allocate(sM(4,4,lociNr),isM(4,4,lociNr),ACGT(4,lociNr),iACGT(4,lociNr))
      allocate(gamma1(lociNr),igamma1(lociNr),gcat(lociNr),invariant(lociNr),iinvariant(lociNr))
      if (dofiles(3)) then
         allocate(tableC(w_int,n-1,lociNr))           
      endif 
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
      call read_input(nameB,nameL,dn,gr,en,n,radio,Ne,growth,bid,samplinfo,events,NM,MS,TM,mlbls,MigMat,gen,sexratio,&
                   &ploidy,lociNr,nb,marker,mutrate,sM,ACGT,gamma1,gcat,invariant,labels,stgroups,r,PrInfo,prheadings)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate( taxa(n),ages(n),nodes(n-1),coal_t(n-1),NodRanges(2,n-1) )  
      allocate( MasterTree(3,(2*n)-1),mut_Mtree(4,(2*n)-1),gammacat(lociNr,MAXVAL(gcat)) ) 
      !------------------------------------------------------------------------------------
      !------------------------------------------------------------------------------------
      !---1.Computing the number of SuSt-----
      g=MAXVAL(stgroups)
      if (dofiles(4)) then   !If summary statistics are required 
         !We get the number of SuSt          
         a=COUNT(SuS(1:8))
         b=COUNT(SuS(9:12))
         SuStN=(a*g + b*(g*(g-1)/2))              !Single-locus number of chosen summary statistics
         SuStNr=SuStN*(1+(lociNr-1)*datatype(3))  !Overall number of summary statistics
      else                                                
         SuStNr=0    
      endif 
      !---2.Now we prepare the headings of the table-----      
      allocate(headings(SuStNr+r+1),irow(SuStN),row(SuStNr+r+1),obsrow(SuStNr+r+1))      
      headings(SuStNr+2:SuStNr+r+1)=prheadings 
      if ((lociNr>1).and.(datatype(3)==1)) then  !Multilocus+independent
         allocate(iheadings(SuStN+1)) 
         headings(1)='Simulation'
         z=2         
         do i=1,lociNr,1
            write( char1, '(I12)' )  i           
            call SuSLabels(g,SuS,SuStN+r+1,iheadings)
            do j=1,SuStN,1
               headings(z)='LOC'//TRIM(ADJUSTL(char1))//'-'//TRIM(iheadings(j+1))
               z=z+1
            enddo
         enddo
         deallocate(iheadings)
      else
         call SuSLabels(g,SuS,SuStNr+r+1,headings)
      endif    
      !---3.Computing the summary statistics of the empirical alignment-----
      allocate(table(w_int,SuStNr+r+1))
      table=0 
	  write( char, '(I5)' )  ithread
      nameD=TRIM(nameA)//'-'//TRIM(ADJUSTL(char))//'.csv'
	  !nameD=TRIM(nameA)//'.csv'  
      if ((empal).and.(dofiles(4))) then !empal:= empirical FASTA available; dofiles(4):= SuSt required 
         obsrow=0                        !read FASTA->compute SuSt 
         z=2
         do l=1,lociNr,1 
            allocate(eAlign(n,nb(l))) 
            !--Each locus comes in a different FASTA file
            if (lociNr>1) then 
               nameC=TRIM(nameA)//'-'//TRIM(ADJUSTL(nameL(l)))//'.fas'    
            else
               nameC=TRIM(nameA)//'.fas'
            endif
            !--Reading the file
            call readFASTA(nameC,n,nb(l),eAlign)
            !--Now calculating the empirical SuSt
            call SummaryStatistics(n,nb(l),eAlign,SuStN,stgroups,SuS,irow)
            !--Pouring them in obsrow
            if (  (lociNr==1) .or. ((lociNr>1).and.(datatype(3)==0))  ) then  !Single locus or multilocus+average
               do i=1,SuStN,1 
                  obsrow(i+1)=obsrow(i+1)+irow(i)/REAL(lociNr,8)  !The division by lociNr turns the sum into the average
               enddo
            elseif ( (lociNr>1) .and. (datatype(3)==1) ) then  !Multilocus + independent SuSt    
               do j=1,SuStN,1
                  obsrow(z)=irow(j)
                  z=z+1
               enddo    
            endif 
            !--
            deallocate(eAlign) 
         enddo
         !Finally writing out the empirical SuSt          
         if (.not.appendSim) then
            table(1,:)=obsrow
            call write_results(.true.,nameD,w_int,1,1+SuStNr+r,headings,table)
         endif 
      else  
         if (.not.appendSim) then 
            call write_results(.true.,nameD,w_int,0,1+SuStNr+r,headings,table)
         endif    
      endif    
      !---4.Computing the Neighbor-Joining tree of the empirical alignment-----   
      if (dofiles(5)) then  !dofiles(5):=NJ-tree required (FASTA present is assumed)
         allocate(dij(n,n),Sij(n,n),trij(n,n),gsizes(g)) 
         do i=1,g,1
            m=COUNT(stgroups==i) 
            gsizes(i)=4*m-4
		 enddo
		 maxg=MAXVAL(gsizes)
		 !dtype(2)=0: For NJ-trees being calculated from distances pooled over loci--------
		 if (dtype(2)==0) then	
             allocate(NJts(g,maxg),obsNJts(g,maxg,1),NJtable(w_int*g,maxg,1))
             NJtable=0
             if (empal) then
                obsNJts=0  
                dij=0
                Sij=0
                trij=0
                do l=1,lociNr,1 
                   allocate(eAlign(n,nb(l))) 
                   !--Each locus comes in a different FASTA file
                   if (lociNr>1) then 
                      nameC=TRIM(nameA)//'-'//TRIM(ADJUSTL(nameL(l)))//'.fas'    
                   else
                      nameC=TRIM(nameA)//'.fas'
                   endif
                   !--Reading the file
                   call readFASTA(nameC,n,nb(l),eAlign)
                   !--Now calculating the empirical distance of the l-thieth locus
                   call NJDistanceMatrix(n,nb(l),eAlign,stgroups,dtype(1),dij,Sij,trij) 
                   !--
                   deallocate(eAlign) 
			    enddo 
                !With the composite distances we make the NJ-tree   
                call NJOnlyTree(n,stgroups,maxg,g,dtype(1),dij,Sij,trij,NJts) 
                !--Pouring them in obsNJts
                do i=1,g,1
                   obsNJts(i,:,1)=NJts(i,:)
                enddo 
                !Finally writing out the tree file  
                if (.not.appendSim) then
                   do i=1,g,1
                      NJtable(i,:,1)=NJts(i,:)
                   enddo    
                   call write_results2(.true.,0,nameA,1,maxg,w_int*g,g,gsizes,NJtable(:,:,1))
                endif  
             else
                if (.not.appendSim) then
                   call write_results2(.true.,0,nameA,0,maxg,w_int*g,g,gsizes,NJtable(:,:,1))
                endif
			 endif 
			 
		 !dtype(2)=1: For NJ-trees being calculated for each locus separatedly--------------	 
		 elseif (dtype(2)==1) then
		     allocate(NJts(g,maxg),obsNJts(g,maxg,lociNr),NJtable(w_int*g,maxg,lociNr))
             NJtable=0
             if (empal) then
                obsNJts=0    
                do l=1,lociNr,1 
                   allocate(eAlign(n,nb(l))) 
                   !--Each locus comes in a different FASTA file
                   if (lociNr>1) then 
                      nameC=TRIM(nameA)//'-'//TRIM(ADJUSTL(nameL(l)))//'.fas'    
                   else
                      nameC=TRIM(nameA)//'.fas'
                   endif
                   !--Reading the file
                   call readFASTA(nameC,n,nb(l),eAlign)
                   !--Now calculating the empirical NJ-tree
                   call NJtree(n,nb(l),eAlign,stgroups,maxg,g,dtype(1),NJts) 
                      !call NJDistanceMatrix(n,nb(l),eAlign,stgroups,dtype(1),dij,Sij,trij) 
                      !call NJOnlyTree(n,stgroups,maxg,g,dtype(1),dij,Sij,trij,NJts) 
                   !--Pouring them in obsNJts
                   do i=1,g,1
                      obsNJts(i,:,l)=NJts(i,:)
                   enddo 
                   !Finally writing out the tree file  
                   if (.not.appendSim) then
                      do i=1,g,1
                         NJtable(i,:,l)=NJts(i,:)
                      enddo    
                      call write_results2(.true.,0,nameL(l),1,maxg,w_int*g,g,gsizes,NJtable(:,:,l))
                   endif         
                   !--
                   deallocate(eAlign) 
			    enddo 
             else
                if (.not.appendSim) then
                   do l=1,lociNr,1  
                      call write_results2(.true.,0,nameL(l),0,maxg,w_int*g,g,gsizes,NJtable(:,:,l))
                   enddo
                endif
			 endif   
         !----------------------------------------------------------------------------------			 
		 endif
      endif      
      !---5.Arlequin file-----       
      if ((empal).and.(dofiles(2))) then  
         do l=1,lociNr,1
            allocate(eAlign(n,nb(l)))  
            !--Each locus comes in a different FASTA file
            if (lociNr>1) then 
               nameC=TRIM(nameA)//'-'//TRIM(ADJUSTL(nameL(l)))//'.fas'    
            else
               nameC=TRIM(nameA)//'.fas'
            endif
            !--Reading the file
            call readFASTA(nameC,n,nb(l),eAlign)
            !Then making the files 
            if (lociNr>1) then 
               nameC=TRIM(nameA)//'-'//TRIM(ADJUSTL(nameL(l)))    
            else
               nameC=TRIM(nameA)
            endif
            call do_arlequin(nameC,n,nb(l),stgroups,eAlign,labels,0) 
            !--
            deallocate(eAlign) 
         enddo
      endif
      !----------------------------------------------------------------------------------
      !----------------------------------------------------------------------------------
      if ( (empal) .and. ( (dofiles(2)).or.(dofiles(4)).or.(dofiles(5)) ) .and. (mimic) ) then
         z=MAXVAL(nb) 
         allocate(MASTEReAlign(lociNr,n,z)) 
         MASTEReAlign=0
         do l=1,lociNr,1
            allocate(eAlign(n,nb(l))) 
            !--Each locus comes in a different FASTA file
            if (lociNr>1) then 
               nameC=TRIM(nameA)//'-'//TRIM(ADJUSTL(nameL(l)))//'.fas'    
            else
               nameC=TRIM(nameA)//'.fas'
            endif
            !--Reading the file
            call readFASTA(nameC,n,nb(l),eAlign)
            !--..and pouring down to MASTER
            do j=1,n,1
               do i=1,nb(l),1    
                  MASTEReAlign(l,j,i)=eAlign(j,i)       
               enddo    
            enddo
            deallocate(eAlign)
         enddo   
      endif
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !S I M U L A T I O N --------------------------------------------------------------------------------------------------
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call RLUXGO(4,seed1,0,0) 
      gammacat=0
      first=.true.
      do i=1,simul_nr,1         
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          isamplinfo=samplinfo
          iNe=Ne
          igrowth=growth            
          ievents=events
          iMigMat=MigMat
          igen=gen
          isexratio=sexratio
          imutrate=mutrate
          isM=sM          
          iACGT=ACGT
          igamma1=gamma1
          iinvariant=invariant
          iPrInfo=PrInfo
          ilabels=labels
          istgroups=stgroups
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
          !Assigning priors 
          row=0
          if (r>0) then
             call Get_Priors(lociNr,r,gr,isamplinfo,dn,iNe,igrowth,en,ievents,NM,iMigMat,igen,isexratio,imutrate,isM,iACGT,igamma1,iinvariant,iPrInfo,Priors)             
             do j=1,r,1
                row(SuStNr+1+j)=Priors(j)
             enddo
          endif 
          !Resolving 'Match' commands 
          if ( (ANY(ievents(:,2)==-777777)).or.(ANY(ievents(:,3)==-777777)).or.(ANY(igrowth==-777777)) ) then
             call Get_Matches(igen,dn,iNe,igrowth,en,ievents)
          endif            
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !Making the stuff moronproof
          call order_samplinfo(gr,dn,n,isamplinfo,ilabels,istgroups)                        
          !Rescaling the age of the samples: everything is calendar time
          m=1
          do j=1,gr,1
             do k=1,isamplinfo(1,j),1    !Filling the array "ages"
                ages(m)=isamplinfo(3,j)
                m=m+1
             enddo
          enddo
          !Adjusting events
          if (en>1) then                 !First is for sorting events according to age
             call OrderEvents(dn,en,ievents)
          endif                          !Second adjustment is to name blocks after their order (they are ordered by age not by number..            
          call OrderEvents2(dn,en,ievents)   !...It is necessary since The Monster (in tree_building) assumes they are so. 
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    
           !!GUI stuff**********************************************
           if ((MOD(i,s_int)==0).or.(i==1))  then  
              write(*,*) 'Simulation', i
           !   call Simuldraft(i,nameA,n,dn,en,ievents,gr,isamplinfo,ilabels,iNe,igrowth,demesranges,demecolor,draftsize)              
		   endif 
           !if (i==1) then 
           !   call Messagetext(1,seed1)
           !endif
           !!GUI stuff**********************************************
          
		  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
          !  L  O  C  U  S    C  Y  C  L  E 
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
          a=MOD(i-1,w_int)+1  !a=simulations since the last recording
          z=2
          do l=1,lociNr,1
             allocate( eAlign(n,nb(l)), sAlign(n,nb(l)), igammacat(gcat(l)) ) 
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
             !Adjustments for marker and ploidy type 
             if ((ploidy==1).and.(marker(l)==1)) then              !Autosomal haploid
                adjust = isexratio*(1-isexratio)*4                  
             elseif ((ploidy==2).and.(marker(l)==1)) then          !Autosomal diploid
                adjust = isexratio*(1-isexratio)*8  
             elseif (ploidy==3) then                               !Haplo-diploid   
                adjust = (isexratio*(1-isexratio)*9)/(1+isexratio) 
             elseif ((ploidy==2).and.(marker(l)==2)) then          !X-linked
                adjust = (isexratio*(1-isexratio)*9)/(1+isexratio)      
             elseif ((ploidy==2).and.(marker(l)==3)) then          !Y-linked
                adjust = isexratio      
             elseif ((ploidy==2).and.(marker(l)==4)) then          !Mitochondrial
                adjust = 1-isexratio                                   
             endif  
             ievents(:,2)=ievents(:,2)*adjust
             iNe=iNe*adjust     
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             !Calling the routine that makes the coalescent tree 
             call tree_maker(n,dn,gr,en,igen,isamplinfo,ages,iNe,igrowth,ievents,taxa,nodes,coal_t,NodRanges,ievents(en,1),NM,MS,TM,mlbls,iMigMat)
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             call order_tree(n,taxa,ages,nodes,coal_t,MasterTree) 
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          
              !!GUI stuff**********************************************
              !!Making the graphic
              if ((MOD(i,s_int)==0).or.(i==1))  then 
                 write(*,*) 'Locus #', l 
              !   call graphic(nameL(l),i,n,MasterTree,gtype,graphsize,treesize)
              !   call graphlabels(n,dn,taxa,ilabels,en,ievents,demesranges,demecolor,graphsize,treesize)
              endif 
              !!GUI stuff**********************************************  
          
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             !Calling the subprogram that put mutations into the tree
             call mutate(n,imutrate(l)*nb(l),MasterTree,mut_Mtree)
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             !Making alignment, fasta, arlequin files, and computing summary statistics 
             if (i>1) then
                do j=1,gcat(l),1 
                    igammacat(j)=gammacat(l,j)  
                enddo
             endif
             !---Extracting eAlign from MASTER
             if ( (empal) .and. ( (dofiles(2)).or.(dofiles(4)).or.(dofiles(5)) ) .and. (mimic) ) then
                do j=1,n,1
                   do k=1,nb(l),1    
                      eAlign(j,k)=MASTEReAlign(l,j,k)       
                   enddo    
                enddo
             endif
             !-
             call do_alignment(nameL(l),n,nb(l),mut_Mtree,NodRanges,istgroups,mimic,eAlign,isM(:,:,l),iACGT(:,l),gamma1(l),gcat(l),igammacat,invariant(l),dofiles,i,taxa,ilabels,sAlign)
             !---
             if (i==1) then
                do j=1,gcat(l),1 
                    gammacat(l,j)=igammacat(j)  
                enddo
             endif
             !Pasing out the coalescent times
             if (dofiles(3)) then
                tableC(a,:,l)=coal_t 
             endif 
             !Computing the summary statistics
             if ((dofiles(4)).and.(ANY(SuS))) then
                call SummaryStatistics(n,nb(l),sAlign,SuStNr+r+1,stgroups,SuS,irow)
             endif
             !---Recording the simulated SuSt
             if (  (lociNr==1) .or. ((lociNr>1).and.(datatype(3)==0))  ) then !Single locus or multilocus+average
                do j=1,SuStN,1 
                   row(j+1)=row(j+1)+irow(j)/REAL(lociNr,8)  !The division by lociNr turns the sum into the average
                enddo
             elseif ( (lociNr>1) .and. (datatype(3)==1) ) then  !Multilocus + independent SuSt    
               do j=1,SuStN,1
                  row(z)=irow(j)
                  z=z+1
               enddo    
             endif 
             !Computing the Neighbor Joining tree
             if (dofiles(5)) then
                if (dtype(2)==0) then 
                   if (l==1) then
                      dij=0
                      Sij=0
                      trij=0
                   endif   
                   call NJDistanceMatrix(n,nb(l),sAlign,stgroups,dtype(1),dij,Sij,trij)
                elseif (dtype(2)==1) then    
                   call NJtree(n,nb(l),sAlign,stgroups,maxg,g,dtype(1),NJts)                 
                   do k=1,g,1
                      NJtable((a-1)*g+k,:,l)=NJts(k,:)
                   enddo 
                endif
             endif  
             !---  
             deallocate(eAlign,sAlign,igammacat)
          enddo
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
          !  E  N  D    O  F    L  O  C  U  S    C  Y  C  L  E
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !RECORDING AND WRITING OUT RESULTS
          !For SuSt
          row(1)=i
          do j=1,SuStNr+r+1,1
             table(a,j)=row(j)
          enddo
          !For NJ-trees (if multilocus composite)
          if ((dofiles(5)).and.(dtype(2)==0)) then
             call NJOnlyTree(n,stgroups,maxg,g,dtype(1),dij,Sij,trij,NJts) 
             do k=1,g,1
                NJtable((a-1)*g+k,:,1)=NJts(k,:)
             enddo
          endif           
          !Write out, if w_int sim complete or run finished
          if ((a==w_int).or.(i==simul_nr)) then 
             !Writing SuSt 
             call write_results(.false.,nameD,w_int,a,1+SuStNr+r,headings,table)
             !Writing the coalescent times
             if (dofiles(3)) then 
                do l=1,lociNr,1 
                   call write_results1(first,appendSim,nameL(l),a,n-1,tableC(:,:,l)) 
                enddo
             endif  
             !Writing the NJ-trees
             if (dofiles(5)) then
                if (dtype(2)==0) then 
                   call write_results2(.false.,i-a+1,nameA,a,maxg,w_int*g,g,gsizes,NJtable(:,:,1))       
                else    
                   do l=1,lociNr,1
                      call write_results2(.false.,i-a+1,nameL(l),a,maxg,w_int*g,g,gsizes,NJtable(:,:,l))
                   enddo
                endif
             endif
             !------
             if (first) then
                first=.false.  !This result in 'first' being true only the first time 
             endif    
          endif  
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo
      !----Appending the empirical SuS & NJt if appending results option was enabled
      deallocate(table)
      if ((empal).and.(dofiles(4)).and.(appendSim)) then !If simulations are appended empirical SuS go at the end         
         allocate(table(1,SuStNr+r+1))
         table(1,:)=0.0D00
         do i=1,SuStNr+1,1
            table(1,i)=obsrow(i)
         enddo         
         call write_results(.false.,nameD,w_int,1,1+SuStNr+r,headings,table)
         deallocate(table)
      endif
      !--
      if (dofiles(5)) then  !NJ-tree 
         deallocate(NJtable)
      endif   
      if ((empal).and.(dofiles(5)).and.(appendSim)) then !If simulations are appended empirical NJt goes at the end
         allocate(NJtable(1*g,maxg,lociNr))
         do l=1,lociNr,1 
            do i=1,g,1
               NJtable(i,:,l)=obsNJts(i,:,l)
            enddo    
            call write_results2(.false.,0,nameL(l),1,maxg,g,g,gsizes,NJtable(:,:,l))
         enddo  
         deallocate(NJtable)
      endif
      !**************************************************************************************************
      !**************************************************************************************************
      !**************************************************************************************************
      !---DEALLLOCATION FESTIVAL 
      if ( (empal) .and. ( (dofiles(2)).or.(dofiles(4)).or.(dofiles(5)) ) .and. (mimic) ) then          
         deallocate(MASTEReAlign) 
      endif
      if (dofiles(5)) then  !NJ-tree 
         deallocate(dij,Sij,trij,gsizes)
         deallocate(NJts,obsNJts)
      endif  
      if (dofiles(3)) then  !Coal times
         deallocate(tableC)    
      endif 
      deallocate(nameL,samplinfo,isamplinfo,Ne,iNe,growth,igrowth,bid,events,ievents)
      deallocate(MigMat,iMigMat,MS,TM,mlbls,labels,ilabels,stgroups,istgroups)
      deallocate(PrInfo,iPrInfo,Priors,prheadings)      
      deallocate(nb,marker,mutrate,imutrate,sM,isM,ACGT,iACGT)      
      deallocate(gamma1,igamma1,gcat,invariant,iinvariant)
      deallocate(taxa,ages,nodes,coal_t,NodRanges,MasterTree,mut_Mtree,gammacat)      
      deallocate(headings,irow,row,obsrow)
         
       !!GUI stuff**********************************************
       !call Messagetext(3,seed1) 
       !deallocate(demesranges,demecolor)
       !CALL ErrorMessage(9999) 
       !ret = CLICKMENUQQ(loc(WINEXIT))
       !!GUI stuff**********************************************
	  
	  call cpu_time(t2)
      write(*,*) NameA, "Time: ", t2-t1, "sec"
      !read(*,*) i
	        
       
      stop
      
EndPRogram Coalescent766
        
          