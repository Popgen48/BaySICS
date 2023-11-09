module routines
    
              
   implicit none

	contains


	  
	

 
  subroutine read_info(radio,s_int,l_int,w_int,simul_nr,appendSim,dofiles,distype,empal,mimic,SuS,datatype)
      !read_info(fname,radio,gtype,treesize,s_int,w_int,simul_nr,appendSim,dofiles,distype,empal,mimic,SuS,datatype,seed1)
        
    implicit none
        
    !Externals
    !character(30), intent(out):: fname              !Name of the input file
    logical, intent(out):: radio                    !Choice for radiocarbon dates calibration
    !integer, intent(out):: gtype                    !This is for choosing type of coalescent graph (dendrogram or cladogram)
    !real, intent(out):: treesize                    !Size of the coalescent graph
    integer, intent(out):: s_int                    !Simulation interval for number display
	integer, intent(out):: l_int                    !Loci interval for locus number display 
    integer, intent(out):: w_int                    !Interval for writing results
    integer, intent(out):: simul_nr                 !Number of simulations  
    logical, intent(out):: appendSim                !Yes/No for appending the simulations to an existing table    
    logical, dimension(5), intent(out):: dofiles    !Options for: (1)fasta, (2)arlequin, (3)coal-times, (4)summary statistics, & (5)Neighbor-Joining files
    integer, dimension(2), intent(out):: distype    !Type of distance for the NJ tree: 1.Nr. of diffs; 2.p-distance; 3.Jukes-Cantor; 4.Kimura2-P    
    logical, intent(out):: empal                    !True if empirical alignment (in fasta) is provided
    logical, intent(out):: mimic                    !True if simulated alignments mimic the '?'/'n' pattern of empirical files 
    logical, dimension(12), intent(out):: SuS       !Choices of summary statistics
    integer, dimension(3), intent(out):: datatype   !Cell 1: DNA/SNPs (1/2), Cell 2: multilocus?YES/NO(=2/1), Cell 3: average/independent (1/2) [if cell 2=2]
    !integer, intent(out):: seed1                    !Seed for random number generation 
    !Internals
    integer:: i, ierror 
    character(30):: text 
    character(12):: char12
    character(12):: char3
    !---------------------------------------------------------------------------------
    open(UNIT=1, FILE='info.dat', STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file         
    openif: if (ierror==0) then                                              
       !read(1,*,IOSTAT=ierror) text, fname
       !Radiocarbon dates calibration
       read(1,*,IOSTAT=ierror) text, char12
       radio=yesorno(char3)       
       !!Dendrogram or cladogram graphic
       !read(1,*,IOSTAT=ierror) text, char12
       !if ((TRIM(char12)=='Dendro').or.(TRIM(char12)=='DENDRO').or.(TRIM(char12)=='dendro')) then
       !   gtype=1
       !elseif ((TRIM(char12)=='Clado').or.(TRIM(char12)=='CLADO').or.(TRIM(char12)=='clado')) then
       !   gtype=2 
       !else
       !   call ErrorMessage(1100)    
       !endif          
       !!More stuff
       !read(1,*,IOSTAT=ierror) text, treesize   
       read(1,*,IOSTAT=ierror) text, s_int, l_int  
       read(1,*,IOSTAT=ierror) text, w_int
       read(1,*,IOSTAT=ierror) text, simul_nr
       !Appending the simulations to an existing table
       READ(1,*,IOSTAT=ierror) text, char3
       appendSim=yesorno(char3)       
       
       !-------------------------------------
       !Type of data (DNA/SNPs)
       read(1,*,IOSTAT=ierror) text, char12
       datatype=0
       if ((TRIM(char12)=='DNA').or.(TRIM(char12)=='dna')) then
          datatype(1)=1
       elseif ((TRIM(char12)=='SNP').or.(TRIM(char12)=='SNPS').or.(TRIM(char12)=='SNPs').or.(TRIM(char12)=='snp').or.(TRIM(char12)=='snps')) then
          datatype(1)=2
       else
          call ErrorMessage(1101)  
       endif 
         if (datatype(1)==1) then !If DNA: multilocus (yes/no)
            read(1,*,IOSTAT=ierror) text, char3 
            if ((TRIM(char3)=='YES').or.(TRIM(char3)=='Yes').or.(TRIM(char3)=='yes')) then   
               datatype(2)=2
            elseif ((TRIM(char3)=='NO').or.(TRIM(char3)=='No').or.(TRIM(char3)=='no')) then     
               datatype(2)=1
            else
               call ErrorMessage(2020)   
            endif
         endif
         
       !-------------------------------------
       !Dofiles (creation of files)
       !---YES|NO to create FASTA files
       read(1,*,IOSTAT=ierror) text, char3  
       dofiles(1)=yesorno(char3)  
       !---YES|NO to create arlequin files
       read(1,*,IOSTAT=ierror) text, char3  
       dofiles(2)=yesorno(char3)
       !---YES|NO to create coalescent times files
       read(1,*,IOSTAT=ierror) text, char3  
       dofiles(3)=yesorno(char3) 
       !---YES|NO to create summary statistics files
       read(1,*,IOSTAT=ierror) text, char3  
       dofiles(4)=yesorno(char3)
         if ((datatype(1)==1).and.(datatype(2)==2).and.(dofiles(4))) then  !If (DNA) & (multilocus) & (SuSt)
            read(1,*,IOSTAT=ierror) text, char12   !SuSt are averaged or taken separately for each locus
            if ((TRIM(char12)=='AVERAGE').or.(TRIM(char12)=='Average').or.(TRIM(char12)=='average')) then
               datatype(3)=0       !This is needed to calculate the number of columns in the reference table
            elseif ((TRIM(char12)=='SEPARATE').or.(TRIM(char12)=='Separate').or.(TRIM(char12)=='separate')) then
               datatype(3)=1 
            else
               call ErrorMessage(1102)    
            endif  
         endif
       !YES|NO to create Neighbor-Joining trees
       read(1,*,IOSTAT=ierror) text, char3  
       dofiles(5)=yesorno(char3)
         if (dofiles(5)) then
            read(1,*,IOSTAT=ierror) text, distype(1)
            if (datatype(1)==1) then !If DNA 
               if (datatype(2)==2) then !If multilocus
                  read(1,*,IOSTAT=ierror) text, char12
                  if ((TRIM(char12)=='COMPOSITE').or.(TRIM(char12)=='Composite').or.(TRIM(char12)=='composite')) then
                     distype(2)=0       !This allows to calculate the number of columns in the reference table
                  elseif ((TRIM(char12)=='SEPARATE').or.(TRIM(char12)=='Separate').or.(TRIM(char12)=='separate')) then
                     distype(2)=1 
                  else
                     call ErrorMessage(1102)    
                  endif 
               else
                  distype(2)=0  !The default is 0 (it is needed to process single-locus)       
               endif   
            endif    
         else
            distype=0 
         endif  
       
       !Now read the options of summary statistics
       if (dofiles(4)) then
          read(1,*,IOSTAT=ierror) text
          do i=1,12,1
             read(1,*,IOSTAT=ierror) text, char3
             SuS(i)=yesorno(char3) 
          enddo
       endif
       
       !Empirical alignment
       read(1,*,IOSTAT=ierror) text, char3
       empal=yesorno(char3)   
       !Mimic choice
       if (empal) then
          read(1,*,IOSTAT=ierror) text, char3
          mimic=yesorno(char3)          
       endif   
       
       !The random number seed
       !read(1,*,IOSTAT=ierror) text, seed1
       if (ierror/=0) then
          call ErrorMessage(1002)  
       endif
    else
       call ErrorMessage(1003)  
    endif openif         
    close(UNIT=1)
        
     
    return
     
  end subroutine read_info 

  
  
  logical function yesorno(chardata)
       
     implicit none
       
     character(12), intent(in):: chardata
     !---
     if ((TRIM(chardata)=='YES').or.(TRIM(chardata)=='Yes').or.(TRIM(chardata)=='yes')) then   
        yesorno=.true.
     elseif ((TRIM(chardata)=='NO').or.(TRIM(chardata)=='No').or.(TRIM(chardata)=='no')) then     
        yesorno=.false.
     else
        call ErrorMessage(2020)   
     endif
       
       
     return
       
   end function yesorno   
   
    
    
    
  SUBROUTINE preread_input(file_name,gr,dn,en,n,r,NM,lociNr)
        
    IMPLICIT NONE
        
    !Externals
    CHARACTER(30), INTENT(IN):: file_name       !Name of the input file
    INTEGER, INTENT(OUT):: dn                   !Nr of subpopulations (demes)
    INTEGER, INTENT(OUT):: gr                   !Nr of sampling groups
    INTEGER, INTENT(OUT):: en                   !Nr of demographic events
    INTEGER, INTENT(OUT):: n                    !Overall sample size (over all samples) 
    INTEGER, INTENT(OUT):: r                    !Nr of priors  
    INTEGER, INTENT(OUT):: NM                   !Nr of migration matrixes 
    INTEGER, INTENT(OUT):: lociNr               !Nr of loci
    !Internals
    INTEGER:: ierror 
    !---------------------------------------------------------------------------------
    OPEN(UNIT=1, FILE=file_name, STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file
       openif: IF (ierror==0) THEN                                              
          READ(1,*,IOSTAT=ierror) gr, dn, en, n, r, NM, lociNr
       ELSE
          CALL ErrorMessage(1005)
       ENDIF openif               
    CLOSE(UNIT=1)
    
     
    RETURN
     
  END SUBROUTINE preread_input 


    
    
  SUBROUTINE read_input(file_name,lociNames,dn,gr,en,n,radio,Ne,growth,bid,samplinfo,events,NM,MS,TM,Mtags,MigMat,gentime,&
    &sexratio,ploidy,lociNr,nb,marker,mutrate,sM,ACGT,gamma1,gammacat,invariant,labels,stgroups,r,PrInfo,headings)
        
    IMPLICIT NONE
        
    !Externals
    CHARACTER(30), INTENT(IN):: file_name                       !Name of the input file
    CHARACTER(30), DIMENSION(lociNr), INTENT(OUT):: lociNames   !Name of the different loci
    INTEGER, INTENT(IN):: dn                                    ! # of subpopulations (demes)
    INTEGER, INTENT(IN):: gr                                    ! # of sampling groups
    INTEGER, INTENT(IN):: en                                    ! # of demographic events
    INTEGER, INTENT(IN):: n                                     ! Total samlpe size (over all samples) 
    LOGICAL, INTENT(IN):: radio                                 ! True if automatic radiocarbon dating
    REAL(8), DIMENSION(dn), INTENT(OUT):: Ne             ! Vector of population sizes (one for each deme)  
    REAL(8), DIMENSION(dn), INTENT(OUT):: growth         ! Growth rate of the corresponding population (deme)
    REAL(8), DIMENSION(dn), INTENT(OUT):: bid            ! Identifier of population number        
    REAL(8), DIMENSION(3,gr), INTENT(OUT):: samplinfo    ! An array with: 1st col:= subsample sizes; 2nd col:= deme; and 3rd col:= sampling time (if ancient, 0 if not)    
    REAL(8), DIMENSION(en,8), INTENT(OUT):: events       ! Array with the info of ancient demographic events: 
                                                          !(1) time;
                                                          !(2) Ne; 
                                                          !(3) Growth rate; 
                                                          !(4) 1st deme involved; 
                                                          !(5) 2nd deme involved; 
                                                          !(6) % of migrating lineages from 1st deme; 
                                                          !(7) % of migrating lineages from 2nd deme; 
                                                          !(8) Block id
    
    INTEGER, INTENT(IN):: NM                                   !# of migration matrixes (external, for allocating the array)
    INTEGER, DIMENSION(NM), INTENT(OUT):: MS                   !Array with the sizes of those matrix
    REAL, DIMENSION(NM,2), INTENT(OUT):: TM                    !Timespans for the validity of the corresponding matrix   
    INTEGER, DIMENSION(NM,dn+en):: Mtags                       !To save space, we only record the columns/rows with values: this gets the labels of those columns/rows
    REAL(8), INTENT(OUT), DIMENSION(NM,dn+en,dn+en):: MigMat   !Array with the migration matrix (last two dimensions), 1st dimension index matrixes
    
    REAL(8), INTENT(OUT):: gentime                             !Generation time (number of years per generation) 
    REAL(8), INTENT(OUT):: sexratio                            !Proportion of males in the population [0.0-1.0]
    INTEGER, INTENT(OUT):: ploidy                              !Ploidy:= (1) Haploid; (2) Diploid; (3) Haplo-Diploid 
    
    INTEGER, INTENT(IN):: lociNr                               !# of independent loci
    INTEGER, DIMENSION(lociNr), INTENT(OUT):: nb               !Number of nucleotides in the fragment   
    INTEGER, DIMENSION(lociNr), INTENT(OUT):: marker           !Type of marker:= (1) Autosomal; (2) X-linked; (3) Y-linked; (4) Mitochondrial   
    REAL(16), DIMENSION(lociNr), INTENT(OUT):: mutrate         !Mutation rate in substitutions/site/year
    REAL(16), DIMENSION(4,4,lociNr), INTENT(OUT):: sM          !Substitution matrix (rows/columns: A, C, G, T)
    REAL(16), DIMENSION(4,lociNr), INTENT(OUT):: ACGT          !Proportion of A, C, G and T in the fragment
    REAL(16), DIMENSION(lociNr), INTENT(OUT):: gamma1          !Shape parameter of the gamma distribution
    INTEGER, DIMENSION(lociNr), INTENT(OUT):: gammacat         !Number of categories of gamma
    REAL(16), DIMENSION(lociNr), INTENT(OUT):: invariant       !Proportion of invariant sites        
    
    CHARACTER(30), DIMENSION(n), INTENT(OUT):: labels       !Labels of the samples, i.e. the names of the real sequences they represent in the simulations 
    INTEGER, DIMENSION(2,n), INTENT(OUT):: stgroups         !Statistical groups (for SuS computation)
    INTEGER, INTENT(IN):: r                                 !# of priors
    REAL(8), DIMENSION(r,9), INTENT(OUT):: PrInfo           !Array with the info of the priors  
    CHARACTER(17), DIMENSION(r):: headings                  !Headings corresponding to the (prior-sampled) parameters of the reference table 
    !Internals        
    INTEGER:: dn2, gr2, en2, n2, r2, NM2, lociNr2           !Same as dn, gr, en, n, r, and NM but obtained from input file for verification
    CHARACTER(12), DIMENSION(3):: charA                     !For capturing the info of the sampling and demes section
    CHARACTER(12), DIMENSION(8):: charB                     !For capturing the info of the events section
    CHARACTER(12), DIMENSION(2):: char2                     !For capturing pairs of values
    CHARACTER(12), DIMENSION(4):: char4                     !For capturing four values
    CHARACTER(15):: char15                                  !For capturing text
    CHARACTER(12):: char1                                   !For capturing a value
    CHARACTER(12), ALLOCATABLE, DIMENSION(:):: charALLO             !For capturing groups of variable size        
    REAL(8):: y                                                     !Multipurpose
    INTEGER:: i, j, k, l, x, z, ierror, dummy                       !Counters, multipurpose                                                
    INTEGER, ALLOCATABLE, DIMENSION(:):: matrixsize                 !Array with the sizes of the migration matrix
    CHARACTER(12), DIMENSION(8):: charPrInfo                        !For reading the priors info  
    !-----------------------------------------------------------------------------------------------------------------------
    !  S  T  A  R  T
    !-----------------------------------------------------------------------------------------------------------------------
    z=0
    !-------------------------------------------------------------------------------------------------------------
    OPEN(UNIT=1, FILE=file_name, STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file
    !-------------------------------------------------------------------------------------------------------------
    openif: IF (ierror==0) THEN                                              ! ierror=0 means Open was succesfull
       WRITE(*,*)'Succesfull!'                                              
       !---------------------------------------------------------------------------------------------------------
       !(1) Reading the number of demes, sampling groups and events
       READ(1,*,IOSTAT=ierror) gr2, dn2, en2, n2, r2, NM2, lociNr2 
       WRITE(*,1000) gr2
       1000 FORMAT(I5,' sampling groups')
       WRITE(*,1001) dn2
       1001 FORMAT(I5,' demes sampled')
       WRITE(*,1002) en2
       1002 FORMAT(I5,' events')  
       WRITE(*,1003) n2
       1003 FORMAT(I7,' individuals') 
       WRITE(*,1004) r2
       1004 FORMAT(I5,' priors') 
       WRITE(*,1005) NM2
       1005 FORMAT(I5,' migration matrixes')
       IF (lociNr2==1) THEN
          WRITE(*,1006) lociNr2
          1006 FORMAT(I9,' locus') 
       ELSE 
          WRITE(*,1007) lociNr2
          1007 FORMAT(I9,' loci')  
       ENDIF     
       !---------------------------------------------------------------------------------------------------------
       !(2) Reading sampling values (SAMPLINFO)   
       DO i=1,gr,1
          READ(1,*,IOSTAT=ierror)  (charA(l), l=1,3) 
          READ(charA(1),*) samplinfo(1,i)
          READ(charA(2),*) samplinfo(2,i)
          IF ((TRIM(charA(3))=='Prior').OR.(TRIM(charA(3))=='prior').OR.(TRIM(charA(3))=='PRIOR')) THEN
             samplinfo(3,i)=-666666
             !--
             z=z+1
             WRITE( char1, '(I12)' )  i
             headings(z)='Age_'//TRIM(ADJUSTL(char1))
          ELSE  
             READ(charA(3),*) y
             samplinfo(3,i)= calendar_date(y,radio)
          ENDIF 
          IF (ierror/=0) STOP
       ENDDO
       DO i=1,gr,1
          IF (samplinfo(3,i)==-666666) THEN 
             WRITE(*,1021) samplinfo(:,i)
          ELSE
             WRITE(*,1020) samplinfo(:,i) 
          ENDIF    
       ENDDO  
       1020 FORMAT(' ',F4.0,' ind from deme ',F4.0,' w/age ', F7.1)
       1021 FORMAT(' ',F4.0,' ind from deme ',F4.0,' w/age set by a prior')
       IF (n2 /= SUM(samplinfo(1,:))) THEN    
          CALL ErrorMessage(1007) 
       ENDIF    
       IF (MAXVAL(samplinfo(2,:))/=dn) THEN
          CALL ErrorMessage(1016) 
       ENDIF
            
       !---------------------------------------------------------------------------------------------------------           
       !(3) Reading population sizes and growth rates
       DO i=1,dn,1                
          READ(1,*,IOSTAT=ierror)  (charA(l), l=1,3) 
          IF ((TRIM(charA(1))=='Prior').OR.(TRIM(charA(1))=='prior').OR.(TRIM(charA(1))=='PRIOR')) THEN
             Ne(i)=-666666
             !--
             z=z+1
             WRITE( char1, '(I12)' )  i
             headings(z)='Ne_'//TRIM(ADJUSTL(char1))
          ELSEIF ((TRIM(charA(1))=='Match').OR.(TRIM(charA(1))=='match').OR.(TRIM(charA(1))=='MATCH')) THEN
             Ne(i)=-777777
          ELSE  
             READ(charA(1),*) Ne(i) 
          ENDIF 
          IF ((TRIM(charA(2))=='Prior').OR.(TRIM(charA(2))=='prior').OR.(TRIM(charA(2))=='PRIOR')) THEN
             growth(i)=-666666
             !--
             z=z+1
             WRITE( char1, '(I12)' )  i
             headings(z)='Growth_'//TRIM(ADJUSTL(char1))
          ELSEIF ((TRIM(charA(2))=='Match').OR.(TRIM(charA(2))=='match').OR.(TRIM(charA(2))=='MATCH')) THEN
             growth(i)=-777777
          ELSE   
             READ(charA(2),*) growth(i) 
          ENDIF 
          READ(charA(3),*) bid(i) 
          IF (ierror/=0) STOP
          !--
          IF ((Ne(i)/=-666666).AND.(growth(i)/=-666666).AND.(growth(i)/=-555555)) THEN  !value/value
             WRITE(*,1010) Ne(i), growth(i), bid(i)             
          ELSEIF ((Ne(i)/=-666666).AND.(growth(i)==-666666)) THEN                       !value/prior
             WRITE(*,1011) Ne(i), bid(i) 
          ELSEIF ((Ne(i)/=-666666).AND.(growth(i)/=-555555)) THEN                       !value/match
             WRITE(*,1012) Ne(i), bid(i)  
          ELSEIF ((Ne(i)==-666666).AND.(growth(i)/=-666666).AND.(growth(i)/=-555555)) THEN  !prior/value
             WRITE(*,1013) growth(i), bid(i)             
          ELSEIF ((Ne(i)==-666666).AND.(growth(i)==-666666)) THEN                           !prior/prior
             WRITE(*,1014) bid(i) 
          ELSEIF ((Ne(i)==-666666).AND.(growth(i)/=-555555)) THEN                           !prior/match
             WRITE(*,1015) bid(i)                   
          ELSEIF ((Ne(i)==-666666).AND.(growth(i)/=-666666).AND.(growth(i)/=-555555)) THEN  !match/value
             WRITE(*,1016) growth(i), bid(i)             
          ELSEIF ((Ne(i)==-666666).AND.(growth(i)==-666666)) THEN                           !match/prior
             WRITE(*,1017) bid(i) 
          ENDIF    
          1010 FORMAT(' Ne is',F11.1,' with growth rate of ',F11.6, ' in population # ', F3.0)
          1011 FORMAT(' Ne is',F11.1,' with growth rate sampled from a prior in population # ', F3.0)
          1012 FORMAT(' Ne is',F11.1,' with growth rate matched (computed) in population # ', F3.0)
          1013 FORMAT(' Ne is sampled from a prior, with growth rate of ',F11.6, ' in population # ', F3.0)
          1014 FORMAT(' Ne is sampled from a prior, and so is the growth rate in population # ', F3.0)
          1015 FORMAT(' Ne is sampled from a prior, and the growth rate matched (computed) in population # ', F3.0)
          1016 FORMAT(' Ne is matched (computed), with growth rate of ',F11.6, ' in population # ', F3.0)
          1017 FORMAT(' Ne is matched (computed), with growth rate sampled from a prior in population # ', F3.0)
       ENDDO
       IF (MAXVAL(bid)/=dn) THEN
          CALL ErrorMessage(1017)
       ENDIF 
       DO i=1,dn,1
          IF ((Ne(i)<1).AND.(Ne(i)/=-666666)) THEN
             CALL ErrorMessage(1018)
          ENDIF
       ENDDO 
            
       !---------------------------------------------------------------------------------------------------------            
       !(4) Reading events
       IF (en>0) THEN  
          WRITE(*,*) '-----------------------------------------------------------------------'
          WRITE(*,*) '-----------------------------------------------------------------------'
          WRITE(*,*)'Demographic events' 
          DO i=1,en,1 
             READ(1,*,IOSTAT=ierror)  (charB(l), l=1,8)  
             !1st column, time 
             IF ((TRIM(charB(1))=='Prior').OR.(TRIM(charB(1))=='prior').OR.(TRIM(charB(1))=='PRIOR')) THEN   
                events(i,1)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  i
                headings(z)='Ev_time_'//TRIM(ADJUSTL(char1))
             ELSE  
                READ(charB(1),*) events(i,1) 
             ENDIF  
             !2nd column, Ne  
             IF ((TRIM(charB(2))=='Prior').OR.(TRIM(charB(2))=='prior').OR.(TRIM(charB(2))=='PRIOR')) THEN   
                events(i,2)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  i
                headings(z)='Ev_Ne_'//TRIM(ADJUSTL(char1))
             ELSEIF ((TRIM(charB(2))=='Match').OR.(TRIM(charB(2))=='match').OR.(TRIM(charB(2))=='MATCH')) THEN   
                events(i,2)=-777777
             ELSE
                READ(charB(2),*) events(i,2) 
             ENDIF  
             !3rd column, growth  
             IF ((TRIM(charB(3))=='Prior').OR.(TRIM(charB(3))=='prior').OR.(TRIM(charB(3))=='PRIOR')) THEN   
                events(i,3)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  i
                headings(z)='Ev_growth_'//TRIM(ADJUSTL(char1))
             ELSEIF ((TRIM(charB(3))=='Match').OR.(TRIM(charB(3))=='match').OR.(TRIM(charB(3))=='MATCH')) THEN   
                events(i,3)=-777777
             ELSE
                READ(charB(3),*) events(i,3) 
             ENDIF  
             !4th column, 1st deme
             READ(charB(4),*) events(i,4)  
             !5th column, 2st deme
             READ(charB(5),*) events(i,5)                  
             !6th column, 1st deme contribution  
             IF ((TRIM(charB(6))=='Prior').OR.(TRIM(charB(6))=='prior').OR.(TRIM(charB(6))=='PRIOR')) THEN   
                events(i,6)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  i
                headings(z)='1st_prop_lin_'//TRIM(ADJUSTL(char1))
             ELSE
                READ(charB(6),*) events(i,6) 
             ENDIF  
             !7th column, 2nd deme contribution  
             IF ((TRIM(charB(7))=='Prior').OR.(TRIM(charB(7))=='prior').OR.(TRIM(charB(7))=='PRIOR')) THEN   
                events(i,7)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  i
                headings(z)='2nd_prop_lin_'//TRIM(ADJUSTL(char1))
             ELSE
                READ(charB(7),*) events(i,7) 
             ENDIF
             !8th column, Block id
             READ(charB(8),*) events(i,8) 
             !------------------------------- 
             WRITE(*,*) '___________________________________'
             WRITE(*,1049) i
             WRITE(*,*) '___________________________________'
             1049 FORMAT(' Event number ',I2) 

             IF (events(i,1)==-666666) THEN  !Time (prior)
                WRITE(*,*) 'Time to the event sampled from a prior' 
             ELSE                            !Time (value)
                WRITE(*,1050) events(i,1)  
                1050 FORMAT('Time to the event = ',F12.1)
             ENDIF
             
             IF (events(i,2)==-666666) THEN       !Ne (prior)
                WRITE(*,*) 'New Ne sampled from a prior'                
             ELSEIF (events(i,2)==-555555) THEN   !Ne (match)
                WRITE(*,*) 'New Ne matched (computed)'
             ELSE                                 !Ne (value)
                WRITE(*,1051) events(i,2) 
                1051 FORMAT(' New Ne =',F10.0)
             ENDIF
             
             IF (events(i,3)==-666666) THEN       !Growth rate (prior)
                WRITE(*,*) 'New growth rate sampled from a prior'                
             ELSEIF (events(i,3)==-555555) THEN   !Growth rate (match)
                WRITE(*,*) 'New growth rate matched (computed)'
             ELSE                                 !Growth rate (value)
                WRITE(*,1052) events(i,3) 
                1052 FORMAT(' New growth rate =',F12.9)
             ENDIF
                          
             IF (events(i,4)==events(i,5)) THEN
                WRITE(*,1053) events(i,4)
                WRITE(*,1056) events(i,6)*100
             ELSE
                WRITE(*,1054) events(i,4)  
                WRITE(*,1055) events(i,5)
                WRITE(*,1057) events(i,6)*100
                WRITE(*,1058) events(i,7)*100
             ENDIF  
             WRITE(*,1060) events(i,8)                  
             1053 FORMAT(' Deme involved =',F10.0)
             1054 FORMAT(' 1st deme involved =',F4.0)
             1055 FORMAT(' 2nd deme involved =',F4.0)
             1056 FORMAT(' Taking',F9.2,'% of its lineages')
             1057 FORMAT(' Taking',F9.2,'% of lineages from 1st deme')
             1058 FORMAT(' Taking',F9.2,'% of lineages from 2nd deme')                   
             1060 FORMAT(' Number of the coalescing block =',F4.0) 
          ENDDO 
          WRITE(*,*) '___________________________________'
       ENDIF 
       DO i=1,en,1
          IF ((events(i,1)<1).AND.(events(i,1)/=-777777).AND.(events(i,1)/=-666666)) THEN
             CALL ErrorMessage(1019)
          ENDIF
          IF ((events(i,2)<1).AND.(events(i,2)/=-777777).AND.(events(i,2)/=-666666)) THEN
             CALL ErrorMessage(1020)
          ENDIF
          IF ((events(i,4)<1).OR.(events(i,4)>dn+en)) THEN
             CALL ErrorMessage(1021)
          ENDIF
          IF ((events(i,5)<1).OR.(events(i,5)>dn+en)) THEN
             CALL ErrorMessage(1021)
          ENDIF
          IF (((events(i,6)<0.0).OR.(events(i,6)>1.0)).AND.(events(i,6)/=-777777).AND.(events(i,6)/=-666666)) THEN
             CALL ErrorMessage(1022)
          ENDIF
          IF (((events(i,7)<0.0).OR.(events(i,7)>1.0)).AND.(events(i,7)/=-777777).AND.(events(i,7)/=-666666)) THEN
             CALL ErrorMessage(1022)
          ENDIF
       ENDDO
       IF ( (ANY(events(:,8)>dn+en)).OR.(ANY(events(:,8)<dn)) ) THEN
          CALL ErrorMessage(1023)
       ENDIF
       DO i=1,en,1
          IF ( (events(i,4)==events(i,8)) .OR. (events(i,5)==events(i,8)) ) THEN
             CALL ErrorMessage(1025)
          ENDIF
       ENDDO
            
       !---------------------------------------------------------------------------------------------------------
       !(5) Reading migration matrix
       IF (NM>0) THEN 
          MigMat=0 
          DO i=1,NM,1                        !For each matrix
             READ(1,*,IOSTAT=ierror) MS(i)                            !The size (doesn't accept prior) 
             READ(1,*,IOSTAT=ierror) (char2(l), l=1,2)                !The time span (check for priors)
             IF ((TRIM(char2(1))=='Prior').OR.(TRIM(char2(1))=='prior').OR.(TRIM(char2(1))=='PRIOR')) THEN   
                TM(i,1)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  i
                headings(z)='MigMat_'//TRIM(ADJUSTL(char1))//'_time1'
             ELSE  
                READ(char2(1),*) TM(i,1) 
             ENDIF  
             IF ((TRIM(char2(2))=='Prior').OR.(TRIM(char2(2))=='prior').OR.(TRIM(char2(2))=='PRIOR')) THEN   
                TM(i,2)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  i
                headings(z)='MigMat_'//TRIM(ADJUSTL(char1))//'_time2'
             ELSE  
                READ(char2(2),*) TM(i,2) 
             ENDIF  
             !--
             READ(1,*,IOSTAT=ierror) (Mtags(i,j), j=1,MS(i))  !The labels (don't accept priors)
             !--
             ALLOCATE(charALLO(MS(i)))
             DO j=1,MS(i),1   !This cycles over rows
                READ(1,*,IOSTAT=ierror) dummy, (charALLO(k),k=1,MS(i)) 
                DO k=1,MS(i),1
                   IF ((TRIM(charALLO(k))=='Prior').OR.(TRIM(charALLO(k))=='prior').OR.(TRIM(charALLO(k))=='PRIOR')) THEN   
                      MigMat(i,Mtags(i,j),Mtags(i,k))=-666666
                      !--
                      z=z+1
                      WRITE( char1, '(I12)' )  i
                      headings(z)='MigMat_'//TRIM(ADJUSTL(char1))//'-'
                      WRITE( char1, '(I12)' )  j
                      headings(z)=TRIM(headings(z))//'-'//TRIM(ADJUSTL(char1))
                      WRITE( char1, '(I12)' )  k
                      headings(z)=TRIM(headings(z))//TRIM(ADJUSTL(char1))
                   ELSE  
                      READ(charALLO(k),*) MigMat(i,Mtags(i,j),Mtags(i,k)) 
                   ENDIF      
                ENDDO
             ENDDO 
             DEALLOCATE(charALLO)
             !-----------------------------------------------------------------
             WRITE(*,*)'Migration matrix', i 
             WRITE(*,*)'It applies to the timespan ', TM(i,1), TM(i,2)
             WRITE(*,*) (Mtags(i,j), j=1,MS(i))
             DO j=1,MS(i),1
                WRITE(*,1030) Mtags(i,j), (MigMat(i,j,k), k=1,MS(i))    
             ENDDO 
             1030 FORMAT(I3,100F9.5)
          ENDDO 
       ENDIF  
       !Checking
       IF (NM<0) THEN 
          CALL ErrorMessage(1013)
       ENDIF
       DO i=1,NM,1
          DO j=1,dn+en
             DO k=1,dn+en
                IF ((MigMat(i,j,k)<0).AND.(MigMat(i,j,k)/=-666666)) THEN
                   CALL ErrorMessage(1014)
                ENDIF
             ENDDO
          ENDDO
       ENDDO 
            
       !---------------------------------------------------------------------------------------------------------
       !(6) Generation time, marker type, sex ration, mutation rate, fragment size, substitution matrix
       !...and gamma+invariant parameters
       !Generation time
       WRITE(*,*) ' '
       WRITE(*,*) ' '       
       READ(1,*,IOSTAT=ierror) char1              
       IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
          gentime=-666666
          !--
          z=z+1
          headings(z)='gen_time'
       ELSE  
          READ(char1,*) gentime 
       ENDIF 
       IF (gentime==-666666) THEN
          WRITE(*,*) 'The generation time is sampled from a prior'
       ELSE 
          WRITE(*,1301) gentime    
          1301 FORMAT(' Generation time = ',F7.2)
       ENDIF
       
       !Sex ratio 
       READ(1,*,IOSTAT=ierror) char1              
       IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
          sexratio=-666666
          !--
          z=z+1
          headings(z)='sex_ratio'
       ELSE  
          READ(char1,*) sexratio 
       ENDIF
       IF (sexratio==-666666) THEN
          WRITE(*,*) 'The sex ratio is sampled from a prior'
       ELSE 
          WRITE(*,1302) sexratio    
          1302 FORMAT(' Sex ratio (M/Total) = ',F10.6)
       ENDIF
       
       !---Ploidy
       READ(1,*,IOSTAT=ierror) char1    
       IF ((TRIM(char1)=='Haploid').OR.(TRIM(char1)=='HAPLOID').OR.(TRIM(char1)=='haploid')) THEN
          ploidy=1 
          WRITE(*,*) 'The organism is HAPLOID'
       ELSEIF ((TRIM(char1)=='Diploid').OR.(TRIM(char1)=='DIPLOID').OR.(TRIM(char1)=='diploid')) THEN   
          ploidy=2
          WRITE(*,*) 'The organism is DIPLOID'
       ELSEIF ((TRIM(char1)=='Haplo-Diploid').OR.(TRIM(char1)=='HAPLO-DIPLOID').OR.(TRIM(char1)=='haplo-diploid')) THEN   
          ploidy=3  
          WRITE(*,*) 'The organism is HAPLO-DIPLOID'
       ENDIF
       
       !--Some checking
       IF ((gentime<0).AND.(gentime/=-666666)) THEN 
          CALL ErrorMessage(1008)
       ENDIF
       IF (((sexratio<0).OR.(sexratio>1)).AND.(sexratio/=-666666)) THEN 
          CALL ErrorMessage(10085)  
       ENDIF
       !-----------------------------------------------------------------------
       !---Cycle over loci
       DO l=1,lociNr,1    
          READ(1,*,IOSTAT=ierror) lociNames(l) !This is the name of the locus 
          !Fragment size
          READ(1,*,IOSTAT=ierror) nb(l) 
          !---Type of molecular marker
          READ(1,*,IOSTAT=ierror) char15    
          IF ((TRIM(char15)=='AUTOSOMAL').OR.(TRIM(char15)=='Autosomal').OR.(TRIM(char15)=='autosomal')) THEN
             marker(l)=1 
             WRITE(*,*) 'Locus #', l, ' is AUTOSOMAL' 
          ELSEIF ((TRIM(char15)=='X-Linked').OR.(TRIM(char15)=='X-linked').OR.(TRIM(char15)=='X-LINKED').OR.(TRIM(char15)=='x-linked')) THEN   
             marker(l)=2
             WRITE(*,*) 'Locus #', l, ' is X-LINKED'
          ELSEIF ((TRIM(char15)=='Y-Linked').OR.(TRIM(char15)=='Y-linked').OR.(TRIM(char15)=='Y-LINKED').OR.(TRIM(char15)=='y-linked')) THEN 
             marker(l)=3 
             WRITE(*,*) 'Locus #', l, ' is Y-LINKED'
          ELSEIF ((TRIM(char15)=='Mitochondrial').OR.(TRIM(char15)=='MITOCHONDRIAL').OR.(TRIM(char15)=='mitochondrial')) THEN   
             marker(l)=4 
             WRITE(*,*) 'Locus #', l, ' is MITOCHONDRIAL'
          ENDIF             
          !---Mutation rate
          READ(1,*,IOSTAT=ierror) char1              
          IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
             mutrate(l)=-666666
             !--
             z=z+1
             WRITE( char1, '(I12)' )  l
             headings(z)='mut_rate-LOCUS'//TRIM(ADJUSTL(char1))
             WRITE(*,*) 'Mutation rate is sampled from a prior'
          ELSE  
             READ(char1,*) mutrate(l) 
             WRITE(*,1404) mutrate(l) 
             1404 FORMAT('Mutation rate = ',F15.11)             
          ENDIF             
          !---Substitution matrix
          READ(1,*,IOSTAT=ierror) char1, char1, char1, char1   !This is meaningless text
          DO i=1,4,1
             READ(1,*,IOSTAT=ierror) char1, (char4(j), j=1,4)  
             DO j=1,4,1
                IF ((TRIM(char4(j))=='Prior').OR.(TRIM(char4(j))=='prior').OR.(TRIM(char4(j))=='PRIOR')) THEN   
                   sM(i,j,l)=-666666
                   !--
                   z=z+1
                   WRITE( char1, '(I12)' )  l
                   headings(z)='rate-LOCUS'//TRIM(ADJUSTL(char1))
                   WRITE( char1, '(I12)' )  i
                   headings(z)=TRIM(ADJUSTL(headings(z)))//'-'//TRIM(ADJUSTL(char1))
                   WRITE( char1, '(I12)' )  j
                   headings(z)=headings(z)//TRIM(ADJUSTL(char1))
                ELSE  
                   READ(char4(j),*) sM(i,j,l) 
                ENDIF
             ENDDO 
          ENDDO            
          !---Nucleotide equilibrium frequencies
          READ(1,*,IOSTAT=ierror) char1, char1, char1, char1   !This is meaningless text
          READ(1,*,IOSTAT=ierror) (char4(i), i=1,4)  
          DO i=1,4,1
             IF ((TRIM(char4(i))=='Prior').OR.(TRIM(char4(i))=='prior').OR.(TRIM(char4(i))=='PRIOR')) THEN   
                ACGT(i,l)=-666666
                !--
                z=z+1
                WRITE( char1, '(I12)' )  l
                headings(z)='pi-LOCUS'//TRIM(ADJUSTL(char1))
                WRITE( char1, '(I12)' )  i
                headings(z)=TRIM(ADJUSTL(headings(z)))//TRIM(ADJUSTL(char1))
             ELSE  
                READ(char4(i),*) ACGT(i,l) 
             ENDIF 
          ENDDO            
          !---Gamma parameter & number of categories
          READ(1,*,IOSTAT=ierror) char1              
          IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
             gamma1(l)=-666666
             !--
             z=z+1
             WRITE( char1, '(I12)' )  l
             headings(z)='gamma_par-LOCUS'//TRIM(ADJUSTL(char1))
          ELSE  
             READ(char1,*) gamma1(l) 
          ENDIF 
          IF (gamma1(l)/=0) THEN
             READ(1,*,IOSTAT=ierror) gammacat(l)
          ELSE
             gammacat(l)=0 
          ENDIF    
          !---Proportion of invariant sites
          READ(1,*,IOSTAT=ierror) char1  
          IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
            invariant(l)=-666666
            !--
            z=z+1
            WRITE( char1, '(I12)' )  l
            headings(z)='invariant-LOCUS'//TRIM(ADJUSTL(char1))
          ELSE  
             READ(char1,*) invariant(l) 
          ENDIF 
          !---Checking
          IF ((mutrate(l)<0).AND.(mutrate(l)/=-666666)) THEN 
             CALL ErrorMessage(1009)
          ENDIF
          IF ((nb(l)<=1).AND.(nb(l)/=-666666)) THEN 
             CALL ErrorMessage(1010)
          ENDIF
          IF ((gamma1(l)<0.0).AND.(gamma1(l)/=-666666)) THEN 
             CALL ErrorMessage(1011)
          ENDIF
          IF ((invariant(l)<0.0).AND.(invariant(l)/=-666666)) THEN 
             CALL ErrorMessage(1012)
          ENDIF
       ENDDO !End of cycle of loci
                            
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !Getting the heading name of extra-priors
       if (r>z) then
          do i=1,r-z,1
             WRITE( char1, '(I8)' )  i
             headings(z+i)='HYPER-PRIOR-'//TRIM(ADJUSTL(char1)) 
          enddo    
       endif    
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         
                        
       !---------------------------------------------------------------------------------------------------------
       !(7) Reading labels
       readloop1: DO i=1,n,1 
          READ(1,*,IOSTAT=ierror) labels(i), stgroups(1,i), stgroups(2,i)
          IF (ierror/=0) EXIT          ! EXIT if not valid
       ENDDO readloop1
            
       !---------------------------------------------------------------------------------------------------------
       !(8) Reading Priors
       readloop2: DO i=1,r,1
          READ(1,*,IOSTAT=ierror) (charPrInfo(j), j=1,8)  
          !(1) Retrieving the code of the PDF------------
          IF ((TRIM(charPrInfo(1))=='UNIFORM').OR.(TRIM(charPrInfo(1))=='Uniform').OR.(TRIM(charPrInfo(1))=='uniform')) THEN
             PrInfo(i,1)=1 
          ELSEIF ((TRIM(charPrInfo(1))=='LOG-UNIFORM').OR.(TRIM(charPrInfo(1))=='Log-Uniform').OR.&
          &(TRIM(charPrInfo(1))=='Log-uniform').OR.(TRIM(charPrInfo(1))=='log-Uniform').OR.(TRIM(charPrInfo(1))=='log-uniform')&
          &.OR.(TRIM(charPrInfo(1))=='LOGUNIFORM').OR.(TRIM(charPrInfo(1))=='LogUniform').OR.(TRIM(charPrInfo(1))=='Loguniform')&
          &.OR.(TRIM(charPrInfo(1))=='logUniform').OR.(TRIM(charPrInfo(1))=='loguniform')) THEN
             PrInfo(i,1)=2        
          ELSEIF ((TRIM(charPrInfo(1))=='NORMAL').or.(TRIM(charPrInfo(1))=='Normal').or.(TRIM(charPrInfo(1))=='normal')) THEN 
             PrInfo(i,1)=3
          ELSEIF ((TRIM(charPrInfo(1))=='LOG-NORMAL').OR.(TRIM(charPrInfo(1))=='Log-Normal').OR.&
          &(TRIM(charPrInfo(1))=='Log-normal').OR.(TRIM(charPrInfo(1))=='log-Normal').OR.(TRIM(charPrInfo(1))=='log-normal')&
          &.OR.(TRIM(charPrInfo(1))=='LOGNORMAL').OR.(TRIM(charPrInfo(1))=='LogNormal').OR.(TRIM(charPrInfo(1))=='Lognormal')&
          &.OR.(TRIM(charPrInfo(1))=='logNormal').OR.(TRIM(charPrInfo(1))=='lognormal')) THEN
             PrInfo(i,1)=4
          ELSEIF ((TRIM(charPrInfo(1))=='EXPONENTIAL').OR.(TRIM(charPrInfo(1))=='Exponential').OR.&
          &(TRIM(charPrInfo(1))=='exponential')) THEN  
             PrInfo(i,1)=5 
          ELSEIF ((TRIM(charPrInfo(1))=='GAMMA').OR.(TRIM(charPrInfo(1))=='Gamma').OR.(TRIM(charPrInfo(1))=='gamma')) THEN  
             PrInfo(i,1)=6
          ELSEIF ((TRIM(charPrInfo(1))=='GEOMETRIC').OR.(TRIM(charPrInfo(1))=='Geometric').OR.&
          &(TRIM(charPrInfo(1))=='geometric')) THEN  
             PrInfo(i,1)=7
          ELSEIF ((TRIM(charPrInfo(1))=='BINOMIAL').OR.(TRIM(charPrInfo(1))=='Binomial').OR.&
          &(TRIM(charPrInfo(1))=='binomial')) THEN  
             PrInfo(i,1)=8 
          ELSEIF ((TRIM(charPrInfo(1))=='BETA').OR.(TRIM(charPrInfo(1))=='Beta').OR.(TRIM(charPrInfo(1))=='beta')) THEN  
             PrInfo(i,1)=9
          ELSEIF ((TRIM(charPrInfo(1))=='POISSON').OR.(TRIM(charPrInfo(1))=='Poisson').OR.(TRIM(charPrInfo(1))=='poisson')) THEN  
             PrInfo(i,1)=10
          ELSE
             CALL ErrorMessage(1006) 
          ENDIF 
          !(2) Retrieving the sign-----------------------
          IF (TRIM(charPrinfo(2))=='-') THEN
             PrInfo(i,2)=-1.0
          ELSEIF (TRIM(charPrinfo(2))=='+') THEN
             PrInfo(i,2)=1.0 
          ELSE
             CALL ErrorMessage(1015) 
          ENDIF   
          !(3-7) Retrieving the parameters---------------
          DO j=3,7,1 
             IF ((charPrInfo(j)(1:5)=='Prior').or.(charPrInfo(j)(1:5)=='prior').or.(charPrInfo(j)(1:5)=='PRIOR')) THEN
                READ(charPrInfo(j)(6:12),*) x 
                PrInfo(i,j)=-555555000-x   
             ELSE
                READ(charPrInfo(j),*) PrInfo(i,j)  
             ENDIF
          ENDDO  
          !(8) Retrieving the type of prior--------------
          IF (TRIM(charPrInfo(8))=='NOVAR') THEN
             PrInfo(i,8)=0.0
          ELSEIF (TRIM(charPrInfo(8))=='PAR') THEN   
             PrInfo(i,8)=1.0
          ELSE
             CALL ErrorMessage(10061) 
          ENDIF  
       ENDDO readloop2
            
    ELSE
       CALL ErrorMessage(1005)           
    ENDIF openif
        
    CLOSE(UNIT=1)
    WRITE(*,*) ' '
        
    IF (ALLOCATED(matrixsize)) THEN
       DEALLOCATE(matrixsize)
    ENDIF
        
        
    RETURN  
           
  END SUBROUTINE read_input
    
    
    
    
  subroutine readFASTA(name,n,nb,ealign)
   
    implicit none
        
    !Externals
    character(30), intent(in):: name                                !Name of the input file
    integer, intent(in):: n                                         !Overall sample size
    integer, intent(in):: nb                                        !Lenght of the fragment (# of nucleotide positions)
    integer, dimension(n,nb), intent(out):: ealign                  !Alignment   
    !Internals
    character(1), allocatable, dimension(:,:):: ialign              !Alignment   
    integer:: i, j, ierror
    character(30):: label                                           !Names of the sequences (in fact not used)
    !-----------------------------------------------------------------------------------------------------------------------
    !  S  T  A  R  T
    !-----------------------------------------------------------------------------------------------------------------------
    allocate(ialign(n,nb))
    !----------------     
    OPEN(UNIT=1, FILE=name, STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file
    !-------------------------------------------------------------------------------------------------------------
    openif: if (ierror==0) then                                              ! ierror=0 means Open was succesfull
       write(*,*)'Opening file'                                              
       !---------------------------------------------------------------------------------------------------------
       !(1) Reading everything
       readloop: do i = 1, n, 1                   
          read(1,*,IOSTAT=ierror) label 
          read(1,100,IOSTAT=ierror) (ialign(i,j), j=1,nb)
          100 format(100000(A1))                    
          if (ierror/=0) exit          ! EXIT if not valid 
       enddo readloop
       !---------------------------------------------------------------------------------------------------------            
    else
       CALL ErrorMessage(1024)                
    endif openif 
    do i=1,n,1
       do j=1,nb,1
          if ((ialign(i,j)=='A').or.(ialign(i,j)=='a')) then
             ealign(i,j)=1
          elseif ((ialign(i,j)=='C').or.(ialign(i,j)=='c')) then
             ealign(i,j)=2
          elseif ((ialign(i,j)=='G').or.(ialign(i,j)=='g')) then
             ealign(i,j)=3
          elseif ((ialign(i,j)=='T').or.(ialign(i,j)=='t')) then
             ealign(i,j)=4
          elseif (ialign(i,j)=='-') then
             ealign(i,j)=-1     
          else    
             ealign(i,j)=0  
          endif    
       enddo
    enddo
    deallocate(ialign)
            

    return
   
  end subroutine readFASTA
    
        
        
    
  SUBROUTINE Get_Priors(lociNr,r,gr,samplinfo,dn,Ne,growth,en,events,NM,MigMatrix,gentime,&
  &sexratio,mutrate,sM,ACGT,gamma1,invariant,PrInfo,Priors)
   
    implicit none
   
    !Externals
    integer, intent(in):: lociNr                                   !# of loci
    integer, intent(in):: r                                        !# of priors
    integer, intent(in):: gr                                       !# of sampling groups 
    real(8), dimension(3,gr), intent(inout):: samplinfo            !Age, deme and # of sampling groups
    integer, intent(in):: dn                                       !# of demes
    real(8), dimension(dn), intent(inout):: Ne                     !Ne of all the demes
    real(8), dimension(dn), intent(inout):: growth                 !Growth rate of all the demes 
    integer, intent(in):: en                                       !# of events/blocks
    real(8), dimension(en,8), intent(inout):: events               !Events/blocks info
    integer, intent(in):: NM                                       !# of migration matrix
    real(8), dimension(NM,dn+en,dn+en), intent(inout):: MigMatrix  !All them igration matrix      
    real(8), intent(inout) :: gentime                              !Generation time
    real(8), intent(inout):: sexratio                              !Sex ratio (M/L e.g. 0.25 means 25% males)
    
    real(16), dimension(lociNr), intent(inout):: mutrate           !Mutation rate    
    real(16), dimension(4,4,lociNr), intent(inout):: sM            !Substitution matrix    
    real(16), dimension(4,lociNr), intent(inout):: ACGT            !Nucleotides equilibroum ratios (proportions)
    real(16), dimension(lociNr), intent(inout):: gamma1            !Gamma parameter (shape)
    real(16), dimension(lociNr), intent(inout):: invariant         !Ratio of invariant sites
    
    real(8), dimension(r,8), intent(inout):: PrInfo                !Array with the info of the priors
    real(8), dimension(r), intent(out):: Priors                    !Array with the random values simulated from the priors
    !Internals
    integer:: i, j, k, l, z, q, a
    real(8):: x                                                    !X gets the number of cross-referred prior  
    logical, dimension(r):: done                                   !Indicator of a prior already sampled (if needed by another one) 
    logical, dimension(5):: y                                      !'y' register which of the 5 parameters of each prior are ready
    !--------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------
    !We get the priors values first 
    !---------------------------------------------------------------------
    q=r
    done=.false.
    do
       do i=1,r,1
          if (.NOT.done(i)) then
             if ((ANY(INT(PrInfo(i,:)/10000)==-555555)).or.(ANY(INT(PrInfo(i,:)/1000)==-555555))&
             &.or.(ANY(INT(PrInfo(i,:)/100)==-555555)).or.(ANY(INT(PrInfo(i,:)/10)==-555555))) then                 
                !----This is only for getting the bloody reference of prior (number)
                do j=3,7,1
                   if  ((INT(PrInfo(i,j)/10)==(-555555)).or.(INT(PrInfo(i,j)/100)==(-555555))&
                   &.or.(INT(PrInfo(i,j)/1000)==(-555555)).or.(INT(PrInfo(i,j)/10000)==(-555555))) then 
                      a=int(log10(ABS(PrInfo(i,j))))
                      a=int(a-5) !a stores now the number of cifers in the prior reference
                      x=dble(10**a)*int(PrInfo(i,j)/dble(10**a))
                      x=ABS(PrInfo(i,j))+x
                      if (done(int(x))) then
                         PrInfo(i,j)=Priors(int(x)) 
                         y(j-2)=.true. 
                      else
                         y(j-2)=.false. 
                      endif
                   else
                      y(j-2)=.true. 
                   endif    
                enddo
                if (ALL(y)) then  !i.e. if the 5 possible parameters required by i-th prior are ready
                   Priors(i)=xrandom(PrInfo(i,:))
                   done(i)=.true.
                endif    
             else    
                Priors(i)=xrandom(PrInfo(i,:)) 
                done(i)=.true.
             endif
          endif    
       enddo 
       q=COUNT(done)
       if (q==r) exit
    enddo
    !-------------------------------------------------- 
    !..and then we assign the priors: 
    !--------------------------------------------------
    !Checking samplinfo
    z=1
    do i=1,gr,1
       do j=1,3,1    
          if (samplinfo(j,i)==-666666) then
             samplinfo(j,i)=Priors(z) 
             z=z+1
          endif    
       enddo          
    enddo    
    !--------------------------------------------------
    !Checking Ne and growth
    do i=1,dn,1
       if (Ne(i)==-666666) then
          Ne(i)=Priors(z) 
          z=z+1
       endif
       if (growth(i)==-666666) then
          growth(i)=Priors(z) 
          z=z+1
       endif
    enddo  
    !--------------------------------------------------
    !Checking events
    do i=1,en,1
       do j=1,8,1    
          if (events(i,j)==-666666) then                 
             events(i,j)=Priors(z)                 
             if ((j==6).or.(j==7)) then  !If the Prior is in the proportion of lineages we have to adjust the proportion in the block that taes the remaining
                do k=1,en,1              !So we check block by block...
                   do l=4,5,1            !Until we get the one that has the same entring block
                      if ((k==i).and.(events(k,l)==events(i,j-2))) then
                         events(k,l+2)=1-Priors(z)       
                      endif    
                   enddo
                enddo    
             endif  
             z=z+1
          endif    
       enddo          
    enddo 
    !Checking migration matrix
    do i=1,NM,1
       do j=1,dn+en,1 
          do k=1,dn+en,1
             if (MigMatrix(i,j,k)==-666666) then 
                MigMatrix(i,j,k)=Priors(z) 
                z=z+1
             endif
          enddo
       enddo
    enddo
    !Generation time and sexratio
    if (gentime==-666666) then 
       gentime=Priors(z) 
       z=z+1
    endif       
    if (sexratio==-666666) then 
       sexratio=Priors(z) 
       z=z+1
    endif 
    !Substitution-related and other parameters
    do l=1,lociNr,1 
       if (mutrate(l)==-666666) then 
          mutrate(l)=Priors(z) 
          z=z+1
       endif 
    enddo
    do l=1,lociNr,1 
       do i=1,4,1 
          do j=1,4,1    
             if (sM(i,j,l)==-666666) then 
                sM(i,j,l)=Priors(z) 
                z=z+1
             endif
          enddo   
       enddo 
    enddo
    do l=1,lociNr,1
       do i=1,4,1    
          if (ACGT(i,l)==-666666) then 
             ACGT(i,l)=Priors(z) 
             z=z+1
          endif
       enddo
       ACGT(:,l)=ACGT(:,l)/SUM(ACGT(:,l))
    enddo
    do l=1,lociNr,1
       if (gamma1(l)==-666666) then 
          gamma1(l)=Priors(z) 
          z=z+1
       endif  
    enddo
    do l=1,lociNr,1
       if (invariant(l)==-666666) then 
          invariant(l)=Priors(z) 
          z=z+1
       endif
    enddo
       
       
    return
          
  END SUBROUTINE Get_Priors
   
    
    
   
    SUBROUTINE Get_Matches(gen,dn,Ne,growth,en,events)
   
       implicit none
   
       !Variables
       real(8), intent(in):: gen
       integer, intent(in)::dn
       real(8), dimension(dn), intent(in)::Ne
       real(8), dimension(dn), intent(inout):: growth
       integer, intent(in):: en
       real(8), dimension(en,8):: events
       
       integer:: i, j
       integer, dimension(1):: pos
       integer:: a, b
       real(8):: Ne1, Ne2
       real(8):: g1, g2
       real(8):: t, t1, t2
       real(8):: factor
       !Start
       !First we check in the growth column
       do i=1,dn,1
          if (growth(i)==-777777) then
             a=0
             b=0
             do j=1,en !This instruction is only to get the id's of the blocks that have our current block(i-th) as an entring block
                if ((events(j,4)==i).or.(events(j,5)==i)) then
                   if (a==0) then    
                      a=j
                   else
                      b=j
                   endif
                endif   
             enddo 
             if (events(a,4)/=events(a,5)) then !In case the ancient a-th block takesl ineages also from another block/population, we adjust its proportion of Ne
                if (events(a,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of our entring block 
                   factor=events(a,6) !Which is in column 6 if our entring block id is in cloumn 4
                else
                   factor=events(a,7) !..or it is in column 7 if our entring block id was in col 5
                endif
                factor=factor/(events(a,6)+events(a,7))
             else
                factor=1
             endif
             Ne1=events(a,2)*factor !The Ne is in events(:,2)
             if (b/=0) then
                if (events(b,4)/=events(b,5)) then !In case the ancient a-th block takesl ineages also from another block/population, we adjust its proportion of Ne
                   if (events(b,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of our entring block 
                      factor=events(b,6) !Which is in column 6 if our entring block id is in cloumn 4
                   else
                      factor=events(b,7) !..or it is in column 7 if our entring block id was in col 5
                   endif
                   factor=factor/(events(b,6)+events(b,7))
                else
                   factor=1
                endif
                Ne2=events(b,2)*factor !The Ne is in events(:,2)
             else
                Ne2=0 
             endif
             Ne2=Ne1+Ne2
             Ne1=Ne(i)
             t=events(a,1)/gen !-events(i,1)
             !Finnally we calculate the bloody number:
             growth(i)=(log(Ne2)-log(Ne1))/t  
          endif
       enddo    
       !-----------------------------------------------------------------------------
       do i=1,en,1  !We check all the events looking for a Match code (-7777);
          if (events(i,2)==-777777) then  !If the Match is in column 2, it is a Ne so: 
             a=NINT(events(i,4))      !a & b gets the id of the blocks whose Ne have to match (the entring blocks of our block)
             b=NINT(events(i,5))
             if (a<=dn) then    !if the entring block is a population
                g1=growth(a)    !we get growth, Ne and age of the population
                Ne1=Ne(a)
                t1=events(i,1)
             else               !if the entring block is an internal block, the data is in events   
                pos=MAXLOC(events(:,8), Mask= events(:,8)==a) !The entring block is the 
                g1=events(pos(1),3)
                Ne1=events(pos(1),2)
                t1=events(i,1)-events(pos(1),1)
             endif
             if (b/=a) then    !The same made with 'a' is done for 'b' except if it is = a
                 if (b<=dn) then
                    g2=growth(b)
                    Ne2=Ne(b)
                    t2=events(i,1)
                 else
                    pos=MAXLOC(events(:,8), Mask= events(:,8)==b)
                    g2=events(pos(1),3)
                    Ne2=events(pos(1),2)
                    t2=events(i,1)-events(pos(1),1)
                 endif
             else              !if b=a we destroy the secont part of the sum (cuz is repeated)
                g2=1
                Ne2=1
                t2=1
                factor=0
             endif    
             events(i,2) = events(i,6)*(Ne1*(2.718281828459**(t1*g1))) 
             events(i,2) = events(i,2) + events(i,7)*(Ne2*(2.718281828459**(t2*g2))) 
          endif   
          if (events(i,3)==-777777) then !if Match is in column 3 it is a growth rate
             a=0
             b=0
             do j=i+1,en !This instruction is only to get the id's of the blocks that have our current block(i-th) as an entring block
                if ((events(j,4)==events(i,8)).or.(events(j,5)==events(i,8))) then
                   if (a==0) then    
                      a=j
                   else
                      b=j
                   endif
                endif   
             enddo 
             if (events(a,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of our entring block 
                factor=events(a,6) !Which is in column 6 if our entring block id is in cloumn 4
             else
                factor=events(a,7) !..or it is in column 7 if our entring block id was in col 5
             endif
             Ne1=events(a,2)*factor !The Ne is in events(:,2)
             if (b/=0) then
                if (events(b,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of ou entring block 
                   factor=events(b,6) !Which is in column 6 if our entring block id is in cloumn 4
                else
                   factor=events(b,7) !..or it is in column 7 if our entring block id was in col 5
                endif
                Ne2=events(b,2)*factor !The Ne is in events(:,2)
             else
                Ne2=0 
             endif
             Ne2=Ne1+Ne2
             Ne1=events(i,2)
             t=(events(a,1)-events(i,1))/gen
             !Finnally we calculate the bloody number:
             events(i,3)=(log(Ne2)-log(Ne1))/t
          endif    
       enddo
       
       
   
       return
   
    END SUBROUTINE Get_Matches
   
    
    
   
    REAL(8) FUNCTION xrandom(info)
    
    
       USE irandom
       USE luxury
   
       IMPLICIT NONE
     
       !Variables
       REAL(8), DIMENSION(8), INTENT(IN):: info       !Info of the probability distribution of priors
       
       INTEGER:: i
       REAL(8):: x
       REAL, DIMENSION(1):: RVEC
       
       !----------------------------------
       SELECT CASE(int(info(1)))
       CASE(1) !Uniform
          CALL RANLUX (RVEC,1)
          x=RVEC(1)*(info(4)-info(3))+info(3) + info(5) !info(5)=offset
       CASE(2) !Log-uniform
          CALL RANLUX (RVEC,1)
          x=EXP(RVEC(1)*(info(4)-info(3))+info(3)) + info(5) !info(5)=offset   
       CASE(3) !Normal
          x=random_normal()*info(4)+info(3) +info(5) !info(3)=mean, info(4)=standar deviation, info(5)=offset
       CASE(4) !Log-normal
          x=EXP(random_normal()*info(4)+info(3)) + info(5) !info(3)=ln-mean, info(4)=ln-standard deviation, info(5)=offset  
       CASE(5) !Exponential
          x=random_exponential()*info(3) + info(5)   !we use lambda parametrized as the rate (rate=1/lambda), info(5)=offset
       CASE(6) !Gamma
          x=random_gamma(real(info(3)),.true.)
          x=x*info(4)+info(5) !info(3)=shape,info(4)=scale,info(5)=offset         
       CASE(7) !Geometric
          i=0
          do
             i=i+1 
             CALL RANLUX (RVEC,1)
             if (RVEC(1)<info(3)) exit    
          enddo    
          x=i + info(5) !info(3)=p,info(5)=offset 
       CASE(8) !Binomial 
          x=random_binomial2(int(info(3)),real(info(4)),.true.) + info(5) !info(3)=n,info(4)=p,info(5)=offset 
          !x=random_binomial1(int(info(3)),real(info(4)),.true.) + info(5) !info(3)=n,info(4)=p,info(5)=offset 
       CASE(9) !Beta
          x=random_beta(real(info(3)),real(info(4)),.true.) + info(5) !info(3)=alpha,info(4)=beta,info(5)=offset
       CASE(10) !Poisson
          x=random_Poisson(real(info(3)),.true.) + info(5) !info(3)=lambda,info(5)=offset
       ENDSELECT
       x=x*info(2) !+ or -
       xrandom = MAX(info(6),MIN(info(7),x))   !Last checking to boundaries. Here always in natural scale
       !-----------------------------     
           
       
       RETURN          
                
    END FUNCTION xrandom        

    
    
    
    SUBROUTINE order_samplinfo(gr,dn,n,samplinfo,labels,groups)
    
    
        IMPLICIT NONE
    
        !Variables
        INTEGER, INTENT(IN):: gr,dn,n                         !Nr of subsamples, nr of populations, overall sample size
        REAL(8), DIMENSION(3,gr), INTENT(INOUT):: samplinfo   !Info of the sampling: 1st row:=sizes of subsamples; 2nd row:=population; 3rd:=ages of subsamples 
        CHARACTER(30), DIMENSION(n), INTENT(INOUT):: labels   !Names of the samples
        INTEGER, DIMENSION(2,n), INTENT(INOUT):: groups       !Assignation for statistical groups 
        !Internals
        INTEGER:: i,j,l,m,k,r,r2,ref
        LOGICAL::found
        REAL(8), DIMENSION(gr,3):: sinfo
        REAL(8), DIMENSION(gr,3):: isinfo
        CHARACTER(30), DIMENSION(n):: ilabels
        INTEGER, DIMENSION(2,n):: igroups    
        REAL(8), ALLOCATABLE, DIMENSION(:):: x, y
    
        !----------------------------
        DO i=1,gr,1
           sinfo(i,:)=samplinfo(:,i)
        ENDDO
        isinfo=sinfo
        ilabels=labels
        igroups=groups
        !--------
        IF (gr>1) THEN         
           r=0  != # of individuals up to present population
           l=0  !l walks the subsamples (rows of samplinfo)
           m=0
           DO j=1,dn        
              l=l+m                     !This is adding all the overall (before the current, j)
              m=COUNT(sinfo(:,2)==j)    !This obtains the number of sampling groups in population j
              IF (m>1) THEN
                 k=l+2
                 found=.false.
                 DO 
                    IF ( sinfo(k-1,3) > sinfo(k,3) ) THEN
                       found=.true.
                    ENDIF    
                    IF ((k==l+m).or.(found)) EXIT
                    k=k+1
                 ENDDO 
                 IF (found) THEN
                    ALLOCATE(x(m),y(m))
                    y=(/(i,i=1,m)/)
                    DO k=1,m,1
                       x(k)=sinfo(k+l,3) !X will contain the ages (the variable to be sorted)
                    ENDDO   
                    CALL SSORT(x,y,m,2)
                    ref=r
                    DO k=1,m,1
                       r2=ref 
                       isinfo(k+l,3)=x(k)            !X will contain the ages (the variable to be sorted)
                       isinfo(k+l,1)=sinfo(NINT(y(k))+l,1) !These is the subsample size                   
                       !Now the labels
                       DO i=1,NINT(y(k))-1,1
                          r2=r2+nint(sinfo(i+l,1))
                       ENDDO
                       DO i=1,NINT(sinfo(NINT(y(k))+l,1)),1
                          ilabels(r+i)=labels(r2+i)
                          igroups(:,r+i)=groups(:,r2+i)
                       ENDDO
                       r=r+nint(isinfo(k+l,1))
                    ENDDO   
                    DEALLOCATE(x,y)
                 ELSE
                    DO i=1,m
                       r=r+nint(sinfo(l+i,1))       !r saves the number of individuals up to (j-th) population
                    ENDDO   
                 ENDIF
              ELSE
                 r=r+nint(sinfo(l+1,1)) 
              ENDIF
           ENDDO
        ENDIF 
        labels=ilabels
        DO i=1,gr,1
           samplinfo(:,i)=isinfo(i,:)
        ENDDO
        groups=igroups
    
    
        RETURN
    
    ENDSUBROUTINE order_samplinfo
    
      
    
    
SUBROUTINE SSORT (X, Y, N, KFLAG) 
    !***BEGIN PROLOGUE SSORT 
    !***PURPOSE Sort an array and optionally make the same interchanges in 
    ! an auxiliary array. The array may be sorted in increasing 
    ! or decreasing order. A slightly modified QUICKSORT 
    ! algorithm is used. 
    !***LIBRARY SLATEC 
    !***CATEGORY N6A2B 
    !***TYPE SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I) 
    !***KEYWORDS SINGLETON QUICKSORT, SORT, SORTING 
    !***AUTHOR Jones, R. E., (SNLA) C Wisniewski, J. A., (SNLA) 
    !***DESCRIPTION 
    ! 
    ! SSORT sorts array X and optionally makes the same interchanges in 
    ! array Y. The array X may be sorted in increasing order or 
    ! decreasing order. A slightly modified quicksort algorithm is used. 
    ! 
    ! Description of Parameters 
    ! X - array of values to be sorted (usually abscissas) 
    ! Y - array to be (optionally) carried along 
    ! N - number of values in array X to be sorted 
    ! KFLAG - control parameter 
    ! = 2 means sort X in increasing order and carry Y along. 
    ! = 1 means sort X in increasing order (ignoring Y) 
    ! = -1 means sort X in decreasing order (ignoring Y) 
    ! = -2 means sort X in decreasing order and carry Y along. 
    ! 
    !***REFERENCES R. C. Singleton, Algorithm 347, An efficient algorithm 
    ! for sorting with minimal storage, Communications of 
    ! the ACM, 12, 3 (1969), pp. 185-187. 
    !***REVISION HISTORY (YYMMDD) 
    ! 761101 DATE WRITTEN 
    ! 761118 Modified to use the Singleton quicksort algorithm. (JAW) 
    ! 890531 Changed all specific intrinsics to generic. (WRB) 
    ! 890831 Modified array declarations. (WRB) 
    ! 891009 Removed unreferenced statement labels. (WRB) 
    ! 891024 Changed category. (WRB) 
    ! 891024 REVISION DATE from Version 3.2 
    ! 891214 Prologue converted to Version 4.0 format. (BAB) 
    ! 900315 CALLs to XERROR changed to CALLs to XERMSG. (THJ) 
    ! 901012 Declared all variables; changed X,Y to SX,SY. (M. McClain) 
    ! 920501 Reformatted the REFERENCES section. (DWL, WRB) 
    ! 920519 Clarified error messages. (DWL)
    ! 920801 Declarations section rebuilt and code restructured to use 
    ! IF-THEN-ELSE-ENDIF. (RWC, WRB) 
    !***END PROLOGUE SSORT 
    ! .. Scalar Arguments .. 
    INTEGER KFLAG, N 
    ! .. Array Arguments .. 
    REAL(8) X(*), Y(*) 
    ! .. Local Scalars .. 
    REAL(8) R, T, TT, TTY, TY 
    INTEGER I, IJ, J, K, KK, L, M, NN 
    ! .. Local Arrays .. 
    INTEGER IL(21), IU(21) 
    ! .. External Subroutines .. 
    ! None 
    ! .. Intrinsic Functions .. 
    INTRINSIC ABS, INT 
    !***FIRST EXECUTABLE STATEMENT SSORT 
    NN = N 
    IF (NN .LT. 1) THEN
       write(*,*) 'The number of values to be sorted is not positive.' 
       RETURN 
    ENDIF 
    ! 
    KK = ABS(KFLAG) 
    IF (KK.NE.1 .AND. KK.NE.2) THEN
       write(*,*) 'The sort control parameter, K, is not 2, 1, -1, or -2.'
       RETURN 
    ENDIF 
    ! 
    ! Alter array X to get decreasing order if needed 
    ! 
    IF (KFLAG .LE. -1) THEN
       DO 10 I=1,NN
          X(I) = -X(I) 
10     CONTINUE 
    ENDIF 
    ! 
    IF (KK .EQ. 2) GO TO 100  !C Sort X only 
    ! 
    M = 1 
    I = 1 
    J = NN 
    R = 0.375E0 
    ! 
20  IF (I .EQ. J) GO TO 60 
    IF (R .LE. 0.5898437E0) THEN
       R = R+3.90625E-2 
    ELSE 
       R = R-0.21875E0 
    ENDIF 
    ! 
30  K = I   ! C Select a central element of the array and save it in location T 
    ! 
    IJ = I + INT((J-I)*R) 
    T = X(IJ)  ! C If first element of array is greater than T, interchange with T  
    !
    IF (X(I) .GT. T) THEN
       X(IJ) = X(I)
       X(I) = T 
       T = X(IJ) 
    ENDIF 
    L = J    ! C If last element of array is less than than T, interchange with T 
    ! 
    IF (X(J) .LT. T) THEN
       X(IJ) = X(J) 
       X(J) = T 
       T = X(IJ)  ! C If first element of array is greater than T, interchange with T 
       ! 
       IF (X(I) .GT. T) THEN 
          X(IJ) = X(I) 
          X(I) = T 
          T = X(IJ) 
       ENDIF 
    ENDIF  ! C Find an element in the second half of the array which is smaller 
          ! than T 
          ! 
40  L = L-1 
    IF (X(L) .GT. T) GO TO 40    ! C Find an element in the first half of the array which is greater 
                                 ! than T 
                                 ! 
50  K = K+1
    IF (X(K) .LT. T) GO TO 50    ! C Interchange these elements 
                                   ! 
    IF (K .LE. L) THEN
       TT = X(L)
       X(L) = X(K)
       X(K) = TT
       GO TO 40 
    ENDIF     ! C Save upper and lower subscripts of the array yet to be sorted 
                ! 
    IF (L-I .GT. J-K) THEN
       IL(M) = I 
       IU(M) = L 
       I = K 
       M = M+1 
    ELSE 
       IL(M) = K 
       IU(M) = J 
       J = L 
       M = M+1 
    ENDIF 
    GO TO 70        ! C Begin again on another portion of the unsorted array 
                      ! 
60  M = M-1 
    IF (M .EQ. 0) GO TO 190 
    I = IL(M) 
    J = IU(M) 
    ! 
70  IF (J-I .GE. 1) GO TO 30 
    IF (I .EQ. 1) GO TO 20 
    I = I-1
    ! 
80  I = I+1
    IF (I .EQ. J) GO TO 60
    T = X(I+1)
    IF (X(I) .LE. T) GO TO 80 
    K = I 
    ! 
90  X(K+1) = X(K) 
    K = K-1 
    IF (T .LT. X(K)) GO TO 90 
    X(K+1) = T 
    GO TO 80    ! C Sort X and carry Y along 
                             ! 
100 M = 1
    I = 1 
    J = NN
    R = 0.375E0 
    ! 
110 IF (I .EQ. J) GO TO 150 
    IF (R .LE. 0.5898437E0) THEN
       R = R+3.90625E-2 
    ELSE 
       R = R-0.21875E0 
    ENDIF
    ! 
120 K = I    ! C Select a central element of the array and save it in location T 
               ! 
    IJ = I + INT((J-I)*R) 
    T = X(IJ) 
    TY = Y(IJ)   ! C If first element of array is greater than T, interchange with T 
                   ! 
    IF (X(I) .GT. T) THEN
       X(IJ) = X(I) 
       X(I) = T
       T = X(IJ)
       Y(IJ) = Y(I)
       Y(I) = TY
       TY = Y(IJ)
    ENDIF 
    L = J     ! C If last element of array is less than T, interchange with T
                      ! 
    IF (X(J) .LT. T) THEN
       X(IJ) = X(J)
       X(J) = T 
       T = X(IJ)
       Y(IJ) = Y(J)
       Y(J) = TY
       TY = Y(IJ)    ! C If first element of array is greater than T, interchange with T 
                       ! 
       IF (X(I) .GT. T) THEN
          X(IJ) = X(I) 
          X(I) = T 
          T = X(IJ) 
          Y(IJ) = Y(I) 
          Y(I) = TY 
          TY = Y(IJ)
       ENDIF 
    ENDIF      ! C Find an element in the second half of the array which is smaller 
                 ! than T 
                 ! 
130 L = L-1 
    IF (X(L) .GT. T) GO TO 130    ! C Find an element in the first half of the array which is greater 
                                    ! than T 
                                    ! 
140 K = K+1 
    IF (X(K) .LT. T) GO TO 140    ! C Interchange these elements 
                                    ! 
    IF (K .LE. L) THEN
       TT = X(L) 
       X(L) = X(K)
       X(K) = TT
       TTY = Y(L) 
       Y(L) = Y(K) 
       Y(K) = TTY 
       GO TO 130
    ENDIF               ! C Save upper and lower subscripts of the array yet to be sorted 
                          ! 
    IF (L-I .GT. J-K) THEN
       IL(M) = I 
       IU(M) = L 
       I = K 
       M = M+1
    ELSE
       IL(M) = K 
       IU(M) = J 
       J = L 
       M = M+1 
    ENDIF 
    GO TO 160      ! C Begin again on another portion of the unsorted array
                     ! 
150 M = M-1
    IF (M .EQ. 0) GO TO 190
    I = IL(M) 
    J = IU(M) 
    ! 
160 IF (J-I .GE. 1) GO TO 120
    IF (I .EQ. 1) GO TO 110 
    I = I-1 
    ! 
170 I = I+1
    IF (I .EQ. J) GO TO 150 
    T = X(I+1) 
    TY = Y(I+1)
    IF (X(I) .LE. T) GO TO 170
    K = I 
    ! 
180 X(K+1) = X(K) 
    Y(K+1) = Y(K) 
    K = K-1
    IF (T .LT. X(K)) GO TO 180 
    X(K+1) = T 
    Y(K+1) = TY 
    GO TO 170        ! C Clean up 
                       ! 
190 IF (KFLAG .LE. -1) THEN
       DO 200 I=1,NN
          X(I) = -X(I)
200    CONTINUE
    ENDIF
   
    RETURN 

END SUBROUTINE SSORT

    
      
   
    subroutine OrderEvents(dn,en,events)
   
             
          implicit none
          
          !External variables
          integer, intent(in):: dn
          integer, intent(in):: en
          real(8), dimension(en,8), intent(inout):: events
          !Internal variables
          integer:: i, j
          real(8), dimension(en,8):: tevents
          real, dimension(2,dn+en):: blockranges  
          real(8), dimension(en):: x, y                       !Array used for sorting arrays  
          integer:: a, b
          integer, dimension(1):: loc1, loc2
          real:: Mrange1, Mrange2, c
          logical:: switch
          !----------------------------------------------------------------------------
          !START
          !----------------------------------------------------------------------------
          tevents=events
          y=(/(j,j=1,en)/)
          x=events(:,1)
          switch=.false.
          do i=2,en,1
              if (x(i-1)>x(i)) then
                 switch=.true.    
              endif    
          enddo    
          if (switch) then          
              call SSORT(x,y,en,2)
              do j=1,en,1
                 events(j,1)=x(j)
                 events(j,2)=tevents(NINT(y(j)),2)
                 events(j,3)=tevents(NINT(y(j)),3)
                 events(j,4)=tevents(NINT(y(j)),4)
                 events(j,5)=tevents(NINT(y(j)),5)
                 events(j,6)=tevents(NINT(y(j)),6)
                 events(j,7)=tevents(NINT(y(j)),7)
                 events(j,8)=tevents(NINT(y(j)),8)
              enddo          
          endif    
          !Fourth events entring pop/blocks (so that at each sub-tree coalescing they always came from the left)
          blockranges=0
          do i=1,dn
             blockranges(1,i)=i-1
             blockranges(2,i)=i
          enddo    
          do i=1,en,1
             a=int(events(i,4))
             b=int(events(i,5)) 
             if (a/=b) then 
                if (a>dn) then !Getting the place of the block in events array
                   loc1=MINLOC(events(:,8),Mask=events(:,8)>=a)
                else
                   loc1=a 
                endif
                Mrange1=SUM(blockranges(:,a))
                if (b>dn) then
                   loc2=MINLOC(events(:,8),Mask=events(:,8)>=b)
                else
                   loc2=b 
                endif
                Mrange2=SUM(blockranges(:,b))
                if (Mrange1>Mrange2) then
                   c=real(events(i,4))
                   events(i,4)=events(i,5)
                   events(i,5)=c          
                endif 
             endif
             blockranges(1,NINT(events(i,8)))=blockranges(1,NINT(events(i,4)))
             blockranges(2,NINT(events(i,8)))=blockranges(2,NINT(events(i,5)))
          enddo 
           
        
        return            
           
    end subroutine OrderEvents
    
    
    

   subroutine OrderEvents2(dn,en,events)
    
      implicit none
    
      !Externals
      integer, intent(in)::   dn
      integer, intent(in)::   en
      real(8), dimension(en,8), intent(inout):: events
      !Internals  
      integer:: i
      real(8), allocatable, dimension(:):: list
      logical:: switch
      !--------------------------------
      if (en>1) then
         switch=.false.  
         i=1
         do 
            if (events(i,8)/=i+dn) then
               switch=.true.    
            endif 
            if ((switch).or.(i==en)) exit
            i=i+1
         enddo 
         if (switch) then
            allocate(list(en))
            do i=1,en,1
               list(NINT(events(i,8)-dn))=i+dn
            enddo   
            do i=1,en,1
               events(i,8)=i+dn
               if (events(i,4)>dn) then
                  events(i,4)=list(NINT(events(i,4)-dn))
               endif
               if (events(i,5)>dn) then
                  events(i,5)=list(NINT(events(i,5)-dn))
               endif
            enddo  
            deallocate(list)
         endif
      endif
    
    
      return
    
    end subroutine OrderEvents2
    
      
 
    
    
    subroutine do_alignment(locusName,n,nb,mtree,NodRanges,groups,mimic,obs_align,SM,&
    &ACGT,shape,gcat,gammacat,inv,dofiles,sim_n,taxa,labels,sim_align)
      
       use irandom
       use luxury
       
       implicit none
             
       !Externals
       character(30), intent(in):: locusName 
       integer, intent(in):: n                                   !Number of samples
       integer, intent(in):: nb                                  !number of nucleotides, size of the fragment
       real(8), dimension(4,(2*n)-1), intent(in):: mtree         !Tree information: 1.Ages of nodes/taxa, 2.Go-to, 3.Branch lenght 4.Number of mutations
       integer, dimension(2,n-1), intent(in):: NodRanges         !Array with the range of each of the final nodes of the tree 
       integer, dimension(2,n), intent(in):: groups              !A two column array with the information of the statistical groups 
       logical, intent(in):: mimic                               !True if simulated align. imitate the n/? of the obs. alignment
       integer, dimension(n,nb), intent(in):: obs_align          !Alignment of observed sequences 
       real(16), dimension(4,4), intent(in):: sM                 !Substitution matrix (4x4: A,C,G,T)
       real(16), dimension(4), intent(in):: ACGT                 !Parameters Pi(A), Pi(C), Pi(G), Pi(T)
       real(16), intent(in):: shape                              !Parameter of the gamma distribution  
       integer, intent(in):: gcat                                !Number of gamma categories
       real(8), dimension(gcat), intent(inout):: gammacat        !The probabilities of the discrete gamma categories       
       real(16), intent(in):: inv                                !Proportion of invariant sites         
       logical, dimension(5), intent(in):: dofiles               !En/Dis:abling the creation of: (1)fasta, (2)arlequin, (3)coal-times, (4)SuSt, & (5) NJ-trees (only 1 & 2 used here) 
       integer, intent(in):: sim_n                               !Number of simulation
       integer, dimension(n), intent(in):: taxa                  !Same as in the tree building routines
       character(30), dimension(n), intent(in) :: labels         !Individuals labels 
       integer, dimension(n,nb), intent(out):: sim_align         !Alignment of simulated sequences        
       !Internals
       integer:: i, j, h, l                                      !Counters
       integer:: e, k, f, a                                      !Multipurpose
       integer:: b, m, ss                                        !Specific purpose 
       real(8):: g, u, x, y, tot                                 !Multipurpose 
       real(8):: ref1, ref2, tarp, df   
       integer:: nes                                             !Number of effective sites (after discarding invariant)
       logical:: done
       real, dimension(1):: rv1
       real, dimension(2):: rv2
       real, allocatable, dimension(:):: rvec                    !Array of random numbers
       real(8), allocatable, dimension(:):: r_array, r_array2    !Array of times and rates for discrete gamma
       integer, allocatable, dimension(:):: ss_list, ss_list2    !List of positions of segregating sites 
       real(8), allocatable, dimension(:):: rs_list, rs_list2    !List of positions of segregating sites 
       integer, allocatable, dimension(:,:):: alignT, alignTr    !alignT:= alignment of taxa; 
       integer, allocatable, dimension(:,:):: alignN             !alignN:= alignment of nodes 
       real(8), dimension(4,4):: Q                               !Transition probabilities matrix
       real(8), dimension((2*n)-1):: alltimes, at2               !Array with times of both nodes and taxa
       integer, allocatable, dimension(:,:):: alignF             !Full alignment  
       integer, allocatable, dimension(:):: sites                !Used to distribute the segregating sites 
       character, allocatable, dimension(:):: iseq      
       character(30):: file_name                                 !Name of the input file      
       character(12):: file_num                                  !Name with the number of file made character
       integer:: ierror  
       !--------------------------------------------------------
       !   S   T   A   R   T
       !--------------------------------------------------------
       m=int(SUM(mtree(4,:)))  
       !--------------------------------------------------------
       !1. Compute the discrete gamma distribution
       !We use the k categories approach of Yang(1994) [J. Mol Evol.39:306-314] 
       if ((sim_n==1).and.(shape/=0.0)) then
          allocate(r_array(gcat),r_array2(gcat)) 
          r_array(1)=0.0
          df=2*real(shape,8) 
          !1.1. We compute the boundaries of the gamma pdf  
          !(a) We set the initial reference bounds
          tarp=dble(gcat-1)/dble(gcat)
          ref2=10  
          do
             x=ref2 
             call chi2NC(x,df,dble(0.0),y,ierror)    
             if (y>tarp) exit
             ref2=2*ref2        
          enddo 
          !(b) We search X* such that P(Gamma<=X*)=tarp (e.g. 0.25,0.5,0.75 for 4 categories)                       
          do i=gcat-1,1,-1
             ref1=0 
             tarp=dble(i)/dble(gcat)   !The target probability 
             !Now we iterate until getting within 0.00001 distance from target p (tarp)               
             do
                x=(ref1+ref2)/2 
                call chi2NC(x,df,dble(0.0),y,ierror)
                if (y>tarp) then
                   ref2=x     !ref1 stays the same
                elseif (y<tarp) then
                   ref1=x     !ref2 stays the same
                endif
                if (ABS(y-tarp)<0.0000001) exit 
             enddo
             r_array(i+1)=x
             ref2=x
          enddo            
          !(c) We assign the average values
          do i=1,gcat,1
             !First we get the width of the 10000 bins 
             if (i<gcat) then
                u=(r_array(i+1)-r_array(i))/10000 
             else    
                u=1.0
             endif
             tot=0             
             do j=1,10000,1
                !Evaluate the gamma function (alpha=shape, beta=2)
                x=r_array(i)+j*u 
                y=real((shape-1),8)*log(x)-x/2.0-lngamma(dble(shape))-real(shape*0.6931471805599453,8)  ! y = ln(gamma_pdf(x;shape,scale=2))
                tot=tot+exp(y)
             enddo   
             r_array2(i)=tot/10000  !The mean value
          enddo
          !(d) Finally we normalize them to sum 1.0
          r_array2=r_array2/sum(r_array2) 
          gammacat(1)=r_array2(1)
          do i=2,gcat,1
             gammacat(i)=r_array2(i)+gammacat(i-1)   !It's not an error: 2nd term has the sum of the last 2 terms
          enddo
          deallocate(r_array,r_array2)         
       endif
       !---------------------------------------------------------------------------------------------------------------------------
       if (m>0) then  !All the stuff below only makes sense if there are mutations in the tree: otherwise we skip them
          !------------------------------------------------------------------------------------------------------------------------    
          k=2*n-1
          alltimes=mtree(1,:)
          at2=(/(i,i=1,k)/)
          call SSORT(alltimes,at2,k,-2)
          !*********************************
          !We transfer the information to the substitution matrix Q
          Q=0
          Q=0
          do j=1,4,1
             do i=1,4,1
                Q(i,j)=REAL(sM(i,j),8)
             enddo
          enddo    
          !------------------------------------------------------------------------------------------------------------------------
          !1) Discarding the invariant sites---------------------------------------------------------------------------------------
          nes=NINT(nb*(1-inv))    
          !2) Choosing the segregating sites---------------------------------------------------------------------------------------
          allocate(ss_list(m),ss_list2(m),rs_list(m),rs_list2(m))   
          if (shape==0.0) then   !Shape parameter set to 0.0 means no heterogeneity accross sites 
             do i=1,m,1
                CALL RANLUX(rv1,1)  
                j=NINT(rv1(1)*real(nes)+0.5)  !j gets the position for the mutation             
                ss_list(i)=j
             enddo
          else 
             !2. Asign the mutations  
             do i=1,m,1
                !2.1 We get the category 
                CALL RANLUX(rv2,2)   
                j=1 
                done=.false.
                do     
                   if (rv2(1)<=gammacat(j)) then 
                      done=.true.    
                   endif
                   if (done) exit
                   j=j+1
                enddo
                !2.2 We get the position inside the category
                k=nint((real(nes)/real(gcat)))  !This is the size of the classes
                k=k*(j-1)+NINT(real(k)*rv2(2)+0.5)   !Now we get the absolute position
                ss_list(i)=k
             enddo 
          endif   
          rs_list=dble(ss_list)
          rs_list2=(/(i,i=1,m)/)
          CALL SSORT(rs_list,rs_list2,m,2)
          j=1           !This instruction substituties an array like this: 23,34,55,55,67,85,85,85,114,...
          u=rs_list(1)  !for one like this:  1,2,3,3,4,5,5,5,6,...
          ss_list(1)=1
          do i=2,m,1
             if (rs_list(i)==u) then
                ss_list(i)=j
             else
                j=j+1 
                u=rs_list(i)
                ss_list(i)=j 
             endif    
          enddo   
          ss=j   !This is the number of segregating sites
          !We return the ss_list to its random order to be used by the next instructions
          do i=1,m,1
            ss_list2(INT(rs_list2(i)))=ss_list(i)   
          enddo    
          !------------------------------------------------------------------------------------------------------------------------ 
          !2) Making the alignment only with segregating sites---------------------------------------------------------------------
          allocate(alignT(n,ss),alignN(n-1,ss),rvec(ss)) 
          !Asigning original sequences at random                      
          call ranlux(rvec,ss)         !1:=A, 2:=C, 3:=G, 4:=T
          do i=1,ss,1
             if ((rvec(i)>0).and.(rvec(i)<=ACGT(1))) then
                a=1
             elseif ((rvec(i)>ACGT(1)).and.(rvec(i)<=ACGT(1)+ACGT(2))) then
                a=2 
             elseif ((rvec(i)>ACGT(1)+ACGT(2)).and.(rvec(i)<=ACGT(1)+ACGT(2)+ACGT(3))) then
                a=3
             elseif (rvec(i)>=1-ACGT(4)) then
                a=4 
             endif 
             alignT(:,i)=a
             alignN(:,i)=a
          enddo  
          a=1
          g=MAXVAL(mtree(1,:))  
          do h=1,2*n-1,1   !We move along the entries of mTree
             e=NINT(at2(h))   !at2 has the times (both coal and ages sorted in decreasing order)
             g=mtree(1,e)
             if (mtree(4,e)>0) then 
                do j=1,int(mtree(4,e)),1   !For each mutation in the branch
                   !Finding the nucleotide in the ancestral node to our current branch
                   if (MOD(e,2)==0) then         !It is a node
                      f=alignN(int(real(e)/2),ss_list2(a))
                   else if (MOD(e,2)==1) then    !It is a taxa 
                      f=alignT(int(real(e+1)/2),ss_list2(a))
                   endif    
                   !Now we get the substitution
                   b=Rij(f,Q)
                   !We assign the mutated nucleotide to all taxa/nodes downstream the current lineage
                   if (MOD(e,2)==1) then
                      alignT(int(real(e+1)/2),ss_list2(a))=b
                   elseif (MOD(e,2)==0) then
                      do l=NodRanges(1,int(real(e)/2)),NodRanges(2,int(real(e)/2))-1,1
                         alignT(l,ss_list2(a))=b 
                         alignN(l,ss_list2(a))=b
                      enddo
                      alignT(NodRanges(2,int(real(e)/2)),ss_list2(a))=b
                   endif
                   !---...and move on
                   a=a+1
                enddo
             endif
          enddo
          !------------------------------------------------------------------------------------------------------------------------
          !3) Making the full alignment--------------------------------------------------------------------------------------------
          allocate(alignF(n,nb),alignTr(n,ss),sites(nb))
          alignTr=0 
          do j=1,ss,1
             do i=1,n,1     
                alignTr(taxa(i),j)=alignT(i,j)
             enddo
          enddo   
          call random_order(sites,nb) 
          do i=1,nb,1
             if (sites(i)<=ss) then                 
                alignF(:,i)=alignTr(:,sites(i)) 
             else
                call ranlux(rv1,1)    
                if ((rv1(1)>0).and.(rv1(1)<=ACGT(1))) then
                   b=1
                elseif ((rv1(1)>ACGT(1)).and.(rv1(1)<=ACGT(1)+ACGT(2))) then
                   b=2 
                elseif ((rv1(1)>ACGT(1)+ACGT(2)).and.(rv1(1)<=ACGT(1)+ACGT(2)+ACGT(3))) then
                   b=3
                elseif (rv1(1)>=1-ACGT(4)) then
                   b=4 
                endif
                alignF(:,i)=b
             endif
          enddo 
          deallocate(alignTr)   
       !-------------------------------------------------------------------------------------------------------------------------
       else !If not mutations occurred then we make the alignment without polymorphic sites
       !-------------------------------------------------------------------------------------------------------------------------
          allocate(alignF(n,nb)) 
          do i=1,nb,1
             call ranlux(rv1,1)    
             if ((rv1(1)>0).and.(rv1(1)<=ACGT(1))) then
                b=1
             elseif ((rv1(1)>ACGT(1)).and.(rv1(1)<=ACGT(1)+ACGT(2))) then
                b=2 
             elseif ((rv1(1)>ACGT(1)+ACGT(2)).and.(rv1(1)<=ACGT(1)+ACGT(2)+ACGT(3))) then
                b=3
             elseif (rv1(1)>=1-ACGT(4)) then
                b=4 
             endif
             alignF(:,i)=b
          enddo
       !-------------------------------------------------------------------------------------------------------------------------   
       endif
       !-------------------------------------------------------------------------------------------------------------------------
       !Now the instruction for the simulated alignment mirror the unknowns' pattern of the empirical alignment
       if (mimic) then
          do j=1,nb,1 
             do i=1,n,1
                if (obs_align(i,j)==0) then
                   alignF(i,j)=0   
                elseif (obs_align(i,j)==-1) then
                   alignF(i,j)=-1    
                endif       
             enddo        
          enddo
       endif   
       do j=1,nb,1
          do i=1,n,1  
             sim_align(i,j)=alignF(i,j)
          enddo
       enddo   
       !---------------------------------------------------------------------------- 
       !3.1) Writting the fasta alignment for taxa----------------------------------
       if (dofiles(1)) then   
          allocate(iseq(nb))
          write(file_num,'(I12)')  sim_n
          file_name = 'sim'//trim(adjustl(locusName))//'-'//trim(adjustl(file_num))//'.fas'
          OPEN(UNIT=1, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)       
          openif0: IF (ierror==0) THEN 
             do i=1,n,1
                write(1,112,iostat=ierror) '>', labels(i) 
                112  format(1x,T1,A1,T2,A30) 
                do j=1,nb,1
                   select case (alignF(i,j))
                    case(-1)   
                       iseq(j)='-'      
                    case(0)   
                       iseq(j)='?'    
                    case(1)   
                       iseq(j)='A'
                    case(2)   
                       iseq(j)='C'
                    case(3)   
                       iseq(j)='G'
                    case(4)   
                       iseq(j)='T'                      
                   end select 
                enddo  
                write(1,125,iostat=ierror) (iseq(l), l=1,nb)
                125  format(1x,T1,100000A1)
             enddo 
          else
              write(1,*)'Failed to write the simulated alignment'
          endif openif0
          CLOSE(UNIT=1)  
          deallocate(iseq)
       endif
       !---------------------------------------------------------------------------- 
       !3.2) Making the Arlequin files----------------------------------------------
       if (dofiles(2)) then
          call do_arlequin(locusName,n,nb,groups,alignF,labels,sim_n)
       endif
       !---------------------------------------------------------------------------------------------------------------------------
       !Goodbye
       if (m>0) then
          deallocate(ss_list,ss_list2,rs_list,rs_list2,alignT,alignN,rvec,sites) 
       endif
       deallocate(alignF) 
  
     
      
       return
      
    end subroutine do_alignment
     
     
     
     
    subroutine do_arlequin(locusName,n,nb,groups,alignF,labels,sim_n)
      
       implicit none
             
       !Externals
       character(30), intent(in):: locusName                     !It gets the name of the fragment/gene
       integer, intent(in):: n                                   !Number of samples
       integer, intent(in):: nb                                  !number of nucleotides, size of the fragment       
       integer, dimension(2,n), intent(in):: groups              !A two column array with the information of the statistical groups 
       integer, dimension(n,nb), intent(in):: alignF             !Full alignment  
       character(30), dimension(n), intent(in) :: labels         !Individuals labels 
       integer, intent(in):: sim_n                               !Number of simulation
       !Internals
       integer:: i, j, h, a, f                                   !Multipurpose
       character, allocatable, dimension(:,:):: arlign           !Alignment for creating arlequin files
       character(30), allocatable, dimension(:) :: ilabels       !Samples labels
       character(30):: file_name                                 !Name of the input file      
       character(12):: file_num                                  !Name with the number of file made character
       !---------------------------------------------------------------------------- 
       !   S   T   A   R   T
       !---------------------------------------------------------------------------- 
       !3.2) Making the Arlequin files----------------------------------------------
       f=COUNT(groups>0)
       allocate(arlign(f,nb),ilabels(f))
       a=1
       do h=1,MAXVAL(groups),1
          do i=1,n,1 
             if ( (groups(1,i)==h) .or. (groups(2,i)==h) ) THEN
                ilabels(a)=labels(i) 
                do j=1,nb,1 
                   select case (alignF(i,j))
                   case(-1)   
                      arlign(a,j)='-'    
                   case(0)   
                      arlign(a,j)='?'     
                   case(1)   
                      arlign(a,j)='A'
                   case(2)   
                      arlign(a,j)='C'
                   case(3)   
                      arlign(a,j)='G'
                   case(4)   
                      arlign(a,j)='T'
                   end select
                enddo 
                a=a+1
             endif 
          enddo    
       enddo 
       write(file_num,'(I12)')  sim_n
       if (sim_n==0) then
          file_name = TRIM(ADJUSTL(locusName))//'-'//TRIM(ADJUSTL(file_num))//'.arp'
       else
          file_name = 'sim'//TRIM(ADJUSTL(locusName))//'-'//TRIM(ADJUSTL(file_num))//'.arp'
       endif    
       !--------------------------------
       call writearpfile(file_name,n,f,nb,ilabels,arlign,groups)
       !--------------------------------
       deallocate(arlign,ilabels)
       !Goodbye
       
      
       return
      
    end subroutine do_arlequin
     
      
     
     
    integer function Rij(ni,Q)
    
       use luxury
       
       implicit none
       !Variables
       real(8), dimension(4,4), intent(in):: Q
       integer, intent(in):: ni
       
       integer:: i
       real(8):: a
       real, dimension(1):: rv
       real(8), dimension(4):: vec
       logical:: done
       !------------
       a=SUM(Q(ni,:))
       vec(1)=Q(ni,1)/a
       vec(2)=(Q(ni,1)+Q(ni,2))/a
       vec(3)=1-(Q(ni,4)/a)
       vec(4)=1.0
       CALL RANLUX(rv,1)   
       i=1 
       do     
          if (rv(1)<=vec(i)) then 
             done=.true.    
          endif
          if (done) exit
          i=i+1
       enddo
       
       Rij=i   
       
              
       return
    
    end function Rij
    
    
    
    
    SUBROUTINE writearpfile(file_name,n,SeqN,nb,labels,arlign,groups)
        
        implicit none
        
        !Variables
        !----------------------------------------------------------------------------------------
        character(30), intent(in):: file_name               !Name of the arp file to be created
        integer, intent(in):: n                             !Number of samples (original)
        integer, intent(in):: SeqN                          !Number of sequences in the alignment (not necesarily the same than n) 
        integer, intent(in):: nb
        character(30), dimension(SeqN), intent(in):: labels
        character, dimension(SeqN,nb), intent(in):: arlign
        integer, dimension(2,n), intent(in):: groups   
        integer:: h, i, j, a, b
        integer:: ierror 
        !character(30+SeqN), dimension(SeqN)::Lines
        !----------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------------------        
        !**************************************************************************************************************************
        !
        OPEN(UNIT=1, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror) ! Opening the file
        openif: IF (ierror==0) THEN                                              ! ierror=0 means Open was succesfull
            !
              !WRITE(*,*)'Succesfull!'           
            !---WRITING DATA------------------
              WRITE(1,*) '# An Arlequin file created by BaySICS'
              WRITE(1,*) '[Profile]'
              WRITE(1,*) '          Title= "Dataset"'
              WRITE(1,*) '          DataType=DNA'
              WRITE(1,91) '          NbSamples=', MAXVAL(groups)
91          FORMAT(1X,A20,I3)            
              WRITE(1,*) '          GenotypicData=0'
              WRITE(1,*) '          MissingData="?"'
              WRITE(1,*) '          LocusSeparator=NONE'
              WRITE(1,*) '[Data]'
              WRITE(1,*) '          [[Samples]]'
            !---here we go-
            a=1
            DO h=1,MAXVAL(groups), 1
               !---h-th sample-
               b=COUNT(groups==h)
                 WRITE(1,92) '            SampleName="Group',h,'"'
92             FORMAT(1X,A29,I3,A1)   
                 WRITE(1,93) '            SampleSize=',b
93             FORMAT(1X,A23,I3)
                 WRITE(1,*) '            SampleData={'   
               writeloop0: DO i=a,a+b-1,1                   
                     WRITE(1,110,iostat=ierror) labels(i),' 1 ',(arlign(i,j), j=1,nb)    
                     a=a+1
                     IF (ierror/=0) EXIT   ! EXIT if not valid
               ENDDO writeloop0 
               WRITE(1,*) '}'
            ENDDO   
                         
            !
        ENDIF openif
        !**************************************************************************************************************************
        !--------------------------------------------------------------------------------------------------------------------------
110     FORMAT(1X,T1,A30,T34,A3,T39,100000A1)
        
        CLOSE(UNIT=1)
        
        
        return
        
    ENDSUBROUTINE writearpfile
       
       
        
        
        
        
       
          
   SUBROUTINE ErrorMessage(code)
    
     !USE ifport
     !USE dflib
     !USE ifqwin

     IMPLICIT NONE
  
     INTEGER, INTENT(IN):: code
     INTEGER:: ret, ret2       
     !------------------------------
     SELECT CASE(code)
      CASE(1002)
        write(*,*) 'The simulations program could not read the options of the "info" file.&
        & Please, check the names and content of the info file.'
        STOP
      CASE(1100)
        write(*,*) 'The simulations program could not read the type of tree graph.&
        & Please, write "dendro" or "clado" in the info file.'
        STOP
      CASE(1101)
        write(*,*) 'The simulations program could not read the type of data.&
        & Please, write DNA_single, DNA_multi, or SNPs in the info file.'
        STOP  
      CASE(1102)
        write(*,*) 'After summary statistics choice, a choice of AVERAGE/INDEPENDENT is required &
        & when multilocus DNA is chosen.' 
        STOP   
      CASE(1003)
        write(*,*) 'The simulations program could not open the info file.&
        & Please check the file name of the infom file and its location in the same directory.'
        STOP      
      CASE(1005)
        write(*,*) 'The input file could not be open or was not in the folder.&
        & Please check that input file is a ".dat" file and is located in the folder.'
        STOP
      CASE(1006)
        write(*,*) 'The probability distribution of the prior was not recognized.&
        & Please, chech the input file and consult the handbook if necessary.'
        STOP
      CASE(10061)
        write(*,*) 'The status of the prior (8th column) should be either NOVAR (noise variable) or PAR (parameter).&
        &Please, add one of those two commands.'
        STOP   
      CASE(1007)
        write(*,*) 'The sample sizes provided at the head of the input file and calculated&
        & from the sampling section do not coincide.'
        STOP
      CASE(1008)
        write(*,*) 'The generation time cannot be negative. Please, check the input file.'
        STOP
      CASE(10085)
        write(*,*) 'The sex ratio has to be a proportion [0.0-1.0]. Please, check the input file.'
        STOP  
      CASE(1009)
        write(*,*) 'The mutation rate cannot be negative. Please, check the input file.'
        STOP
      CASE(1010)
        write(*,*) 'Please, check the number of nucleotides (lenght) of the DNA fragment.'
        STOP
      CASE(1011)
        write(*,*) 'The proportion of invariant sites should be between zero (inclusive) to one.'
        STOP
      CASE(1012)
        write(*,*) 'Shape parameter of (gamma parameter) should be positive.'
        STOP
      CASE(1013)
        write(*,*) 'The number of migration matrix should be a non-negative integer.'
        STOP
      CASE(1014)
        write(*,*) 'The migration rates should be given as probabilities&
        & (i.e. numbers between 0.0 and 1.0).'
        STOP
      CASE(1015)
        write(*,*) 'The name of the prior density should be followed by a "-" or "+" symbol indicating&
        & if the sampled parameter should be mirrored or not.'
        STOP
      CASE(1016)
        write(*,*) 'The number of populations provided at the head of the input file and gotten from the&
        & sampling section do not coincide.'C, 'ERROR!'
        STOP
      CASE(1017)
        write(*,*) 'The labels of the populations should be consecutive integers starting&
        & in 1 (1, 2, 3, ...).'
        STOP
      CASE(1018)
        write(*,*) 'The effective population sizes should be positive and larger than 1.0.'
        STOP
      CASE(1019)
        write(*,*) 'The times to the events should be positive.'C, 'ERROR!'
        STOP
      CASE(1020)
        write(*,*) 'The effective population sizes of the blocks should be larger than 1.'
        STOP
      CASE(1021)
        write(*,*) 'The id number of the events/blocks should continue the populations numbers.'
        STOP
      CASE(1022)
        write(*,*) 'The columns 6th-7th of events should contain a proportion in [0.0-1.0].'
        STOP
      CASE(1023)
        write(*,*) 'The labels of the events/blocks should continue the populations numeration.'
        STOP
      CASE(1024)
        write(*,*) 'The FASTA file was not found. Please, check the name and location of the alignment file.'C,&
        & 'ERROR'
        STOP
      CASE(1025)
        write(*,*) 'The blocks cannot be the source of lineages of themlseves.'
        STOP 
      CASE(1026)
        write(*,*) 'The Neighbor Joining tree cannot be done for less than N=3.&
        & Check the sizes of the statistical groups or deactivate the NJ tree option.'
        STOP   
      CASE(9999)
        write(*,*) 'THE ANALYSIS HAS FINISHED SUCCESFULLY! :)'
     END SELECT          
        
     RETURN

   END SUBROUTINE ErrorMessage

   
   
   
    
   subroutine write_results(first,fname,a,c,b,headings,table)
   
        implicit none
        
        !External Variables
        logical, intent(in):: first                           !If it's the first call to the program 
        character(30), intent(in):: fname                     !Name of the file
        integer, intent(in):: a                               !Table size: rows   
        integer, intent(in):: c                               !Effective number of rows to write   
        integer, intent(in):: b                               !Table size: columns   
        character(17), dimension(b), intent(in):: headings    !First row of the output file
        real(8), dimension(a,b), intent(in):: table           !The data 
        
        !Internal Variables
        integer:: i, j                                     !Counters
        integer:: ierror                                   !Indicator of a malfunction of opening/writing
        logical:: found                                    !Indicator of the exitence of the file
        character(30):: iname  
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        iname=TRIM(fname)
        inquire(FILE=iname, EXIST=found) 
        if (found) then
           if (first) then   !We only want to enter here the first time
              OPEN(UNIT=1, FILE=iname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror) 
           else
              OPEN(UNIT=1, FILE=iname, STATUS='OLD', ACTION='WRITE', ACCESS='APPEND', IOSTAT=ierror)  
           endif       
        else
           OPEN(UNIT=1, FILE=iname, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)                
        endif 
        openif0: if (ierror==0) then
           if (first) then  !This is because we only want to enter here the first time
              write(1,1024,iostat=ierror) (headings(i),i=1,b)
              1024  format(1x,10000(A17))
           endif
           do i=1,c,1     
              write(1,1040,iostat=ierror) INT(table(i,1)), (table(i,j),j=2,b) 
              1040 format(1x,I15,10000(E17.9)) 
           enddo
        else
           write(*,*)'Failed to write the table'
        endif openif0
        CLOSE(UNIT=1)
        
        
        return
        
   end subroutine write_results 
   
   
   
   
   
   subroutine write_results1(first,appendSim,fname,a,b,table)
   !This subroutine writes out a coalescent times file
   
        implicit none
        
        !External Variables
        logical, intent(in):: first                     !True only the first call
        logical, intent(in):: appendSim                 !True if appending the results to an existing table
        character(30), intent(in):: fname               !Name of the file
        integer, intent(in):: a                         !Number of rows in table
        integer, intent(in):: b                         !Number of columns in table
        real(8), dimension(a,b), intent(in):: table     !The content of the file
        !Internals
        integer:: i, j                                  !Counters
        integer:: ierror                                !Indicator of a malfunction of opening/writing
        logical:: found                                 !Indicator of the exitence of the file
        character(30):: iname                           !The actual name of the output file
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
        !Writing the output to the results file
        iname=TRIM(fname)//'_ct.csv'
        !Determining if this is the first call (so we create/reeplace the file)
        if ((.not.appendSim).and.(first)) then      !Enter here on first call and only if not appending  
           inquire(FILE=iname, EXIST=found)
           if (found) then
              open(UNIT=1, FILE=iname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
           else
              open(UNIT=1, FILE=iname, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)
           endif 
        else  
           OPEN(UNIT=1, FILE=iname, STATUS='OLD', ACTION='WRITE', ACCESS='APPEND', IOSTAT=ierror)    
        endif    
        openif0: IF (ierror==0) THEN                  
            writeloop0: DO i = 1, a, 1
                write(1,1040,iostat=ierror)   (table(i,j), j=1,b) 
                1040  format(1x,100000(E17.9))
                IF (ierror/=0) EXIT   ! EXIT if not valid
            enddo writeloop0 
        else
           write(*,*)'Failed to write the table'
        endif openif0
        CLOSE(UNIT=1)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        return
        
   end subroutine write_results1
   
   
   
   
   subroutine write_results2(first,sim,fname,a,b,c,g,gsizes,table)
   !This subroutine writes out the file with the Neighbor Joining trees
   
        implicit none
        
        !External Variables
        logical, intent(in):: first                         !True if this is the first call (to first create the file or append)         
        integer, intent(in):: sim                           !Number of simulation of the first one
        character(30), intent(in):: fname                   !Name of the file
        integer, intent(in):: a                             !Number of simulations in the file
        integer, intent(in):: b                             !Width of the table file below 
        integer, intent(in):: c                             !Size of the table array (may not coincide with the simulations)
        integer, intent(in):: g                             !Number of trees per simulation        
        integer, dimension(g), intent(in):: gsizes          !Sizes of trees rows  
        real(8), dimension(c,b), intent(in):: table         !The trees coded in special format
        !Internals
        character(12):: number                              !To write the simulation number
        integer:: i, j, k, l                                !Counters
        integer:: ierror                                    !Indicator of a malfunction of opening/writing
        logical:: found                                     !Indicator of the exitence of the file
        character(42):: iname                               !The actual name of the output file
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
        !Writing the output to the results file
        iname=TRIM(fname)//'-NJtrees.tree'
        !Determining if this is the first call (so we create/reeplace the file)
        if (first) then !This means this is the first call to the ruoutine
           inquire(FILE=iname, EXIST=found)
           if (found) then
              open(UNIT=1, FILE=iname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
           else
              open(UNIT=1, FILE=iname, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)   
           endif 
        else
           OPEN(UNIT=1, FILE=iname, STATUS='OLD', ACTION='WRITE', ACCESS='APPEND', IOSTAT=ierror)    
        endif    
        openif0: IF (ierror==0) THEN 
            if (first) then
               write(1,*,iostat=ierror) 'Neighbor_joining_trees'
            endif
            k=1
            writeloop0: DO i=1,a,1
               write(number,'(I12)') sim+i-1 
               write(1,*,iostat=ierror) 'Simulation_number_', adjustl(number)
               writeloop1: DO j=1,g,1
                  write(1,1040,iostat=ierror) (table(k,l), l=1,gsizes(j))
                  k=k+1
                  1040  format(1x,100000(E17.9))
                  IF (ierror/=0) EXIT   ! EXIT if not valid
               enddo writeloop1 
            enddo writeloop0
        else
           write(*,*)'Failed to write the table'
        endif openif0
        CLOSE(UNIT=1)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        return
        
   end subroutine write_results2
   
   
   
   
   SUBROUTINE SuSLabels(g,SuS,nr,headings)
   
      implicit none
   
      !Externals
      integer, intent(in):: g                                  !Number of statistical groups
      logical, dimension(12), intent(in):: SuS                 !Array with the choices of SuS
      integer, intent(in):: nr                                 !Lenght of the row  (# of columns in table)
      character(17), dimension(nr), intent(inout):: headings   !Labels of columns of the reference table    
      !Internals
      integer:: i,j, z
      character(3):: name, name2   
      !-----------------------------
      headings(1)='Simulation'
      z=2
      do i=1,g,1
         write(name,'(I3)')  i   
         if (Sus(1)) then        
            headings(z)='HapTypes'//TRIM(ADJUSTL(name))
            z=z+1
         endif
         if (Sus(2)) then        
            headings(z)='PrivHaps'//TRIM(ADJUSTL(name))
            z=z+1
         endif
         if (Sus(3)) then        
            headings(z)='SegSites'//TRIM(ADJUSTL(name))
            z=z+1
         endif
         if (Sus(4)) then        
            headings(z)='PairDiff'//TRIM(ADJUSTL(name))
            z=z+1
         endif
         if (Sus(5)) then        
            headings(z)='NucDiver'//TRIM(ADJUSTL(name))
            z=z+1
         endif
         if (Sus(6)) then        
            headings(z)='GenDiver'//TRIM(ADJUSTL(name))
            z=z+1
         endif
         if (Sus(7)) then        
            headings(z)='TajimasD'//TRIM(ADJUSTL(name))
            z=z+1
         endif
         if (Sus(8)) then        
            headings(z)='FusFs'//TRIM(ADJUSTL(name))
            z=z+1
         endif
      enddo
      do i=1,g-1,1
         do j=i+1,g,1
            write(name,'(I3)')  i  
            write(name2,'(I3)')  j  
            if (Sus(9)) then        
               headings(z)='PairDifG'//TRIM(ADJUSTL(name))//TRIM(ADJUSTL(name2))
               z=z+1
            endif
            if (Sus(10)) then        
               headings(z)='Fst'//TRIM(ADJUSTL(name))//TRIM(ADJUSTL(name2))
               z=z+1
            endif
            if (Sus(11)) then        
               headings(z)='ShareHap'//TRIM(ADJUSTL(name))//TRIM(ADJUSTL(name2))
               z=z+1
            endif
            if (Sus(12)) then        
               headings(z)='SharFreq'//TRIM(ADJUSTL(name))//TRIM(ADJUSTL(name2))
               z=z+1
            endif
         enddo
      enddo
   
      return
      
   END SUBROUTINE SuSLabels
   
   
   
   
   
   SUBROUTINE SummaryStatistics(n,nb,alignF,nr,groups,SuS,row)     !  alignS,taxa,ages,nodes,coal_t,sd) 
   
      implicit none
      
      !External variables
      integer, intent(in):: n                          !Sample size
      integer, intent(in):: nb                         !Number of nucleotides in the analyzed fragment
      integer, dimension(n,nb), intent(in):: alignF    !The alignment
      integer, intent(in):: nr                         !Lenght of the table row   
      integer, dimension(2,n), intent(in):: groups     !Two columna with the stat. groups info 
      logical, dimension(12), intent(in):: SuS         !Summary statistics choices   
      real(8), dimension(nr), intent(out):: row        !Vector of summary statistics 
      !Internal variables
      integer:: i, j, k, l, h, z, a, b, c                 !Counters and multipurpose 
      integer:: m, m1, m2                                 !Counters for statistical groups
      integer:: ss                                        !Number of segregating sites (width of align)
      integer, allocatable, dimension(:,:):: alignS       !Alignment of segregating sites 
      integer, dimension(n,n):: dij, Sij                  !Pairwise differences matrix, effective nr. of ocmparisons of dij 
      integer, allocatable, dimension(:,:):: Ndij         !A transient dij array containing 2 statistical groups (for pairwise comparisons)   
      integer, allocatable, dimension(:,:):: minidij      !Same as dij but only for the taxa of the group
      integer, allocatable, dimension(:,:):: miniSij      !Same as Sij but only for the taxa of the group
      integer, allocatable, dimension(:,:,:):: alldij     !Array containing the arrays dij of each statistical group
      integer, allocatable, dimension(:,:):: minialign    !Alignment including only the samples of the group (it is a subgroup of align)  
      real, allocatable, dimension(:):: freq              !Array with the haplotype frequencies of one group (actually the union of two groups for pairwise fst)  
      real, allocatable, dimension(:,:):: freq2           !Array with allele frequencies of 2 groups
      real(8), allocatable, dimension(:):: GDivA          !Array with the values of genetic diversity (1-sum(pi**2)) for each statistical group
      real(8), allocatable, dimension(:):: pi             !average number of pairwise differences
      integer:: g                                         !Number of statistical groups 
      integer, allocatable, dimension(:):: HapNr          !number of haplotypes in each group
      logical, dimension(n):: Selector                    !Indicates what samples are included for calculating certain summary statistics and whose are not  
      logical, allocatable, dimension(:,:):: Selector2    !Indicates what samples are included for calculating certain summary statistics and whose are not
      integer, allocatable, dimension(:):: s              !Number of segregating sites (in a statistical group)
      real(8):: GD                                        !Gene diversity 
      real(8):: x                                         !Multipurpose 
      !----------------------------------------------------------------------------------------------------------------
      !  S  T  A  R  T  
      !----------------------------------------------------------------------------------------------------------------
      call get_var_nr(n,nb,alignF,ss)
      allocate(alignS(n,ss))
      call get_var_algn(n,nb,alignF,ss,alignS) 
      !MAKING THE DISTANCES MATRIX
      dij=0
      Sij=0
      do j=2,n,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
         do i=1,j-1,1
            do k=1,nb,1  
               if ((alignF(i,k)/=0).and.(alignF(j,k)/=0)) then !No contribution of unknown (n/?) sites 
                  Sij(i,j)=Sij(i,j)+1
                  if (alignF(i,k)/=alignF(j,k)) then 
                     dij(i,j)=dij(i,j)+1
                  endif   
               endif
            enddo
         enddo    
      enddo
      do i=1,n-1,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
         do j=i+1,n,1
            Sij(j,i)=Sij(i,j)
            dij(j,i)=dij(i,j) 
         enddo
      enddo      
      !----------------------------------------------------------------------------------------------------------------
      !MAKING THE SUMMARY STATISTICS FOR EACH GROUP      
      g=MAXVAL(groups)
      if (SuS(7)) then
         allocate(s(g))    
      endif
      if ((SuS(7)).or.(SuS(8))) then
         allocate(pi(g))    
      endif
      if ((SuS(9)).or.(SuS(10))) then
         allocate(GDivA(g),alldij(g,2*n,2*n)) 
      endif 
      if ((SuS(8)).or.(SuS(11))) then
         allocate(HapNr(g))    
      endif
      !---------------------------------------------------------------------------------------------------------------
      !SINGLE GROUP SUMMARY STATISTICS
      z=1 
      do h=1,g,1          
         !First, making Selector (the array with the info for including/excluding samples in minidij)
         !------------------------------------------------------------------------------------------
         Selector=.false.
         do i=1,n,1
             if ( (groups(1,i)==h) .or. (groups(2,i)==h) ) then 
                Selector(i)=.true.
             endif
         enddo  
         !--------------------------------------
         m=COUNT(Selector)   
         !--------------------------------------    
         !Second, making mini distances matrix 
         !------------------------------------------------------------------------------------------ 
         allocate(minidij(m,m),miniSij(m,m),freq(m))
         freq=0
         minidij=0
         a=1
         b=2  
         do i=1,n-1,1
            do j=i+1,n,1
               if ((Selector(i)).and.(Selector(j))) then
                  minidij(a,b)=dij(i,j)
                  minidij(b,a)=dij(i,j)
                  miniSij(a,b)=Sij(i,j)
                  miniSij(b,a)=Sij(i,j)
                  if (b<m) then
                     b=b+1
                  elseif (b==m) then !If we reached the end of the row, then.. 
                     a=a+1           !..we start a new row  
                     b=a+1           !..and the column is the b+1 
                  endif    
               endif
            enddo 
         enddo  
         !In case we have pair diff. calculated among groups
         if (SuS(9)) then
            do j=1,m,1
               do i=1,m,1 
                  alldij(h,i,j)=minidij(i,j) 
               enddo
            enddo   
         endif 
         !Third, making mini alignment  
         !------------------------------------------------------------------------------------------
         allocate(minialign(m,ss))
         a=0
         i=0         
         do 
            i=i+1 
            if (Selector(i)) then
                a=a+1
                minialign(a,:)=alignS(i,:)                
            endif
            if ((a==m).or.(i==n)) exit
         enddo
         !------------------------------------------------------------------------------------------         
         !Finally, ladies and gentleman...the summary statistics (aplauses)         
         !1) HAPTYPES--------------------------------------         
         if ((SuS(1)).or.(SuS(8)).or.(SuS(11))) then            
            a=HapTypes(m,minidij)
            if (SuS(1)) then
               row(z)=a
               z=z+1
            endif
            if ((SuS(8)).or.(SuS(11))) then
               HapNr(h)=a      
            endif
         endif        
         !2) PRIVHAPS--------------------------------------
         if (SuS(2)) then
            row(z)=PrivHaps(n,m,dij,minidij,h,groups)            
            z=z+1
         endif
         !3) SEGSITES--------------------------------------
         if ((SuS(3)).or.(SuS(7))) then
            call get_var_nr(m,ss,minialign,a)
            if (SuS(3)) then
               row(z)=a
               z=z+1
            endif
            if (SuS(7)) then
               s(h)=a    
            endif    
         endif
         !4) PAIRDIFF--------------------------------------
         if ((SuS(4)).or.(SuS(7)).or.(SuS(8))) then
            x=SUM(dble(minidij))/dble(m*(m-1))
            if (SuS(4)) then
               row(z)=x*(real(m)/real(m-1))
               z=z+1
            endif
            if ((SuS(7)).or.(SuS(8))) then
               pi(h)=x     
            endif        
         endif
         !5) NUCDIVER--------------------------------------
         if (SuS(5)) then
            x=0 
            do j=2,m,1
               do i=1,j-1,1
                  x = x + dble(minidij(i,j))/dble(miniSij(i,j)) 
               enddo
            enddo    
            row(z)=x*2/dble(m*(m-1))
            z=z+1
         endif
         !6) GENDIVER--------------------------------------
         if ((SuS(6)).or.(SuS(10))) then
            call frequencies(m,minidij,freq)              
            GD=GenDiver(m,freq)*real(m)/real(m-1)
            if (SuS(6)) then
               row(z)=GD
               z=z+1
            endif
            if (SuS(10)) then
                GDivA(h)=GD    
            endif 
         endif 
         !7) TAJIMASD--------------------------------------
         if (SuS(7)) then
            row(z)=TajimasD(m,pi(h),s(h))
            z=z+1
         endif 
         !8) FU'S FS---------------------------------------
         if (SuS(8)) then
            row(z)=Fs(m,HapNr(h),pi(h))
            z=z+1 
         endif
         !-------------------------------------------------         
         deallocate(minidij,miniSij,freq,minialign) 
      enddo
      !INTER-GROUP SUMMARY STATISTICS 
      !------------------------------------------------------------------------------------------
      if ((SuS(9)).or.(SuS(10).or.(SuS(11)))) then
         do i=1,g-1,1
            do j=i+1,g,1       
               if ((SuS(9)).or.(SuS(10))) then 
                  m1=0
                  do k=1,n,1
                     if ((groups(1,k)==i).or.(groups(2,k)==i)) then 
                        m1=m1+1
                     endif
                  enddo
                  m2=0
                  do k=1,n,1
                     if ( (groups(1,k)==j) .or. (groups(2,k)==j) ) then 
                        m2=m2+1
                     endif
                  enddo
                  allocate(minidij(m1,m2)) 
                  minidij=0
                  a=1
                  b=1
                  do k=1,n,1
                     do l=1,n,1
                        if ( ((groups(1,k)==i).or.(groups(2,k)==i)).and.((groups(1,l)==j).or.(groups(2,l)==j)) ) then
                           minidij(a,b)=dij(k,l)  
                           if (b<m2) then
                               b=b+1
                           elseif (b==m2) then !If we reached the end of the row, then.. 
                               a=a+1           !..we start a new row  
                               b=1             !..and the column is the b+1 
                           endif    
                        endif
                     enddo 
                  enddo
                  !9) PAIRDIFF/POP----------------------------------
                  !The average number of pairwise differences between the populations
                  if (SuS(9)) then
                     row(z)=SUM(dble(minidij))/real(m1*m2)                       
                     z=z+1 
                  endif 
                  !10) FST-------------------------------------------
                  !Fst following...
                  if (SuS(10)) then
                     allocate(Ndij(m1+m2,m1+m2),freq(m1+m2))
                     do l=1,m1,1
                        do k=1,m1,1
                           Ndij(k,l)=alldij(i,k,l) 
                        enddo
                     enddo
                     Ndij=alldij(i,:,:)
                     do k=m1+1,m1+m2,1
                        do l=1,m1,1
                           Ndij(k,l)=minidij(l,k-m1) 
                        enddo
                     enddo   
                     do k=m1+1,m1+m2,1
                        do l=1,m1,1
                           Ndij(l,k)=minidij(l,k-m1) 
                        enddo
                     enddo
                     do k=m1+1,m1+m2,1
                        do l=m1+1,m1+m2,1
                           Ndij(k,l)=alldij(j,k-m1,l-m1) 
                        enddo
                     enddo
                     call frequencies(m1+m2,Ndij,freq)
                     row(z)=fst(i,j,g,m1+m2,GDivA,freq)
                     z=z+1 
                     deallocate(minidij,Ndij,freq)
                  else
                     deallocate(minidij) 
                  endif   
               endif
               !-------------------------------------------------
               if ((SuS(11)).or.(SuS(12))) then
                  Selector=.false.
                  do k=1,n,1
                     if ( (groups(1,k)==i).or.(groups(2,k)==i).or.(groups(1,k)==j).or.(groups(2,k)==j) ) then 
                        Selector(k)=.true.
                     endif
                  enddo  
                  m=COUNT(Selector)                    
                  allocate(minidij(m,m),miniSij(m,m))
                  minidij=0
                  miniSij=0
                  a=1
                  b=2
                  do k=1,n-1,1
                     do l=k+1,n,1  
                        if ((Selector(k)).and.(Selector(l))) then
                           minidij(a,b)=dij(k,l)
                           minidij(b,a)=dij(k,l)
                           miniSij(a,b)=Sij(k,l)
                           miniSij(b,a)=Sij(k,l)
                           if (b<m) then
                              b=b+1
                           elseif (b==m) then !If we reached the end of the row, then.. 
                              a=a+1           !..we start a new row  
                              b=a+1           !..and the column is the b+1 
                           endif    
                        endif
                     enddo 
                  enddo
                  a=HapNr(i)
                  b=HapNr(j)
                  c=HapTypes(m,minidij) 
                  !11) SHAREHAP-------------------------------------
                  !Number of shared haplotypes
                  if (SuS(11)) then
                     row(z)=a+b-c
                     z=z+1  
                  endif 
                  !12) SHAREFREQ-------------------------------------
                  !Sum allele frequencies differences
                  if (SuS(12)) then
                     allocate(freq2(c,2),Selector2(m,2))
                     Selector2=.false.
                     l=1
                     do k=1,n,1
                        if ( (groups(1,k)==i).or.(groups(2,k)==i).or.(groups(1,k)==j).or.(groups(2,k)==j) ) then 
                           if ( (groups(1,k)==i).or.(groups(2,k)==i) ) then 
                              Selector2(l,1)=.true.
                           endif
                           if ( (groups(1,k)==j).or.(groups(2,k)==j) ) then 
                              Selector2(l,2)=.true.
                           endif 
                           l=l+1
                        endif
                     enddo
                     !--
                     call shared_frequencies(m,c,Selector2,minidij,freq2) 
                     !--
                     x=0.0 
                     do l=1,c,1
                        x=x+MIN(freq2(l,1),freq2(l,2))       
                     enddo 
                     row(z)=x
                     z=z+1
                     deallocate(freq2,Selector2)
                  endif 
                  !--------------------------------------------------
                  deallocate(minidij,miniSij)
               endif 
            enddo
         enddo          
      endif
      
      deallocate(alignS) 
      if (allocated(s)) then
         deallocate(s)
      endif 
      if (allocated(pi)) then
         deallocate(pi)
      endif 
      if (allocated(HapNr)) then
         deallocate(HapNr)
      endif 
      if (allocated(GDivA)) then
         deallocate(GDivA)
      endif
      if (allocated(alldij)) then
         deallocate(alldij)
      endif
       
           
          
      return
      
   ENDSUBROUTINE SummaryStatistics
   
   
   
   SUBROUTINE get_var_nr(n,nb,alignment,segsites)
       
           implicit none
       
           !Variables
           integer, intent(in):: n
           integer, intent(in):: nb
           integer, dimension(n,nb), intent(in):: alignment          
           integer, intent(out):: segsites
           !Internals
           integer:: i, j, z, ref 
           logical:: done, done2
           !-------------------------------------------------------
           z=0
           do j=1,nb,1
              !This is only to get a reference nucleotide (in case the 1st is 0)
              i=1 
              done=.false.
              do
                 if (alignment(i,j)/=0) then
                    ref=alignment(i,j)
                    done=.true.
                 endif   
                 if ((done).or.(i==n)) exit !If we exit by i=n there is no defined nucleotide in the column
                 i=i+1
              enddo 
              !Now we compare the reference to each nucleotide
              if (done) then
                 done2=.false.
                 i=1
                 do
                    if ((alignment(i,j)/=0).and.(alignment(i,j)/=ref)) then
                       done2=.true. 
                       z=z+1
                    endif    
                    if ((done2).or.(i==n)) exit
                    i=i+1
                 enddo 
              endif
           enddo
           segsites=z            
           
           return       
       
   END SUBROUTINE get_var_nr
   
   
      
   SUBROUTINE get_var_algn(n,nb,alignment,ss,algn_red)
       
           implicit none
       
           !External variables
           integer, intent(in):: n
           integer, intent(in):: nb
           integer, dimension(n,nb), intent(in):: alignment          
           integer, intent(in):: ss
           integer, dimension(n,ss), intent(out):: algn_red
           !Internal variables
           integer:: i, j, z, ref 
           logical:: done, done2
           !Start
           z=1
           do j=1,nb,1
              !This is only to get a reference nucleotide (in case the 1st is 0)
              i=1 
              done=.false.
              do
                 if (alignment(i,j)/=0) then
                    ref=alignment(i,j)
                    done=.true.
                 endif   
                 if ((done).or.(i==n)) exit !If we exit by i=n there is no defined nucleotide in the column
                 i=i+1
              enddo 
              !Now we compare the reference to each nucleotide
              if (done) then
                 done2=.false.
                 i=1
                 do
                    if ((alignment(i,j)/=0).and.(alignment(i,j)/=ref)) then
                       done2=.true. 
                       algn_red(:,z)=alignment(:,j)
                       z=z+1
                    endif    
                    if ((done2).or.(i==n)) exit
                    i=i+1
                 enddo 
              endif
           enddo
       
       
           return
       
   END SUBROUTINE get_var_algn
   
   
   
   INTEGER FUNCTION HapTypes(n,d) 
   
      implicit none
      
      !External
      integer, intent(in):: n
      integer, dimension(n,n), intent(in):: d  !Matrix of distances 
      !Internal
      integer, dimension(n):: array            !1-dimension array with the taxa
      integer:: i                              !Counters
      integer:: m, r                           !m:= counter of places walked in array; r:= counter of haplotypes 
      !********************************************************************************************************************      
      !NUMBER OF HAPLOTYPES
      do i=1,n,1
         array(i)=i 
      enddo    
      m=0
      r=0
      do 
         m=m+1 
         if (array(m)/=0) then
            r=r+1 
            do i=m+1,n,1 
               if (array(i)/=0) then
                  if (d(array(m),array(i))==0) then
                     array(i)=0
                  endif   
               endif
            enddo
         endif
         if (m==n) exit    
      enddo    
      
      HapTypes=r
      
      
      return
      
   END FUNCTION HapTypes
   
   
   
   INTEGER FUNCTION PrivHaps(n,m,d,md,cur,groups)
   
      implicit none
      
      !External
      integer, intent(in):: n                       !Overall sample size
      integer, intent(in):: m                       !Sample size of the current group
      integer, dimension(n,n), intent(in):: d       !Matrix of distances (of the whole sampling)              
      integer, dimension(m,m), intent(in):: md      !Matrix of distances (of the current group)
      integer, intent(in):: cur                     !Current group being analyzed 
      integer, dimension(2,n), intent(in):: groups  !A two column array with the information of the groups of samples whose summary statistics are gonna be calculated upon
      !Internal
      integer:: i,j                                 !Counters
      integer:: a                                   !'a' counts the positions 
      logical:: new                                 !Switch being true if the haplotype is not repeated (the same) from a previous one 
      logical:: SameHap                             !Switch being true if the haplotypes being compared are equal
      !********************************************************************************************************************      
      !NUMBER OF PRIVATIVE HAPLOTYPES
      PrivHaps=0
      a=0
      do i=1,m,1
         new = .true.  !This is only for checking if the current haplotype to compare (the i-th) isn't the same that  
         if (i>1) then   !..another previous in the list. So comparison only takes place if 'new' is true. 
            do j=i-1,1,-1
               if (md(i,j)==0) then
                   new = .false.
               endif    
            enddo    
         endif
         !First we have to localize the position of our current haplotype (no matter if the comparison takes place or not)
         do
            a=a+1                                               !'a' is the position of the current haplotype in d
            if ((groups(1,a)==cur).or.(groups(2,a)==cur)) exit
         enddo
         !Now the comparison of our current haplotype (the one in position a)
         if (new) then
            SameHap = .false.  
            j=0
            do 
               j=j+1 
               if (  ((groups(1,j)/=0).and.(groups(1,j)/=cur)) .or. ((groups(2,j)/=0).and.(groups(2,j)/=cur))  ) then
                  if (d(a,j)==0) then
                     SameHap = .true. 
                  endif
               endif
               if ((SameHap).or.(j==n)) exit
            enddo
            if (.not.SameHap) then
               PrivHaps=PrivHaps+1 
            endif    
         endif    
      enddo   
         
      return 
      
   END FUNCTION PrivHaps
    
   
   
   SUBROUTINE frequencies(n,d,freq) 
   
      implicit none
      
      !External
      integer, intent(in):: n                  !Sample size
      integer, dimension(n,n), intent(in):: d  !Matrix of distances 
      real, dimension(n), intent(out)::freq    !Array with the haplotype frequencies  
      !Internal
      integer, dimension(n):: array            !1-dimension array with the taxa
      integer:: i                              !Counters
      integer:: m                              !m:= counter of places walked in array;  
      !********************************************************************************************************************      
      !START
      do i=1,n,1
         array(i)=i 
      enddo    
      freq=1
      do m=1,n,1 
         if (array(m)/=0) then
            do i=m+1,n,1 
               if (array(i)/=0) then
                  if (d(array(m),array(i))==0) then
                     array(i)=0
                     freq(i)=0
                     freq(m)=freq(m)+1
                  endif   
               endif
            enddo
         endif 
      enddo 
      freq=freq/SUM(freq)
      
          
      return
      
   END SUBROUTINE frequencies
   
   
   
   SUBROUTINE shared_frequencies(n,a,selector,d,freq) 
   
      implicit none
      
      !External
      integer, intent(in):: n                          !Sample size
      integer, intent(in):: a                          !Number of alleles in both groups   
      logical, dimension(n,2), intent(in):: Selector   !Cells are true if they belong to group 1 (col 1) or 2 (col 2) 
      integer, dimension(n,n), intent(in):: d          !Matrix of distances 
      real, dimension(a,2), intent(out)::freq          !Array with the haplotype frequencies  
      !Internal
      integer, dimension(n):: array            !1-dimension array with the taxa
      integer:: i                              !Counters
      integer:: m, z                           !m:counter for array; k:counter for freq 
      real:: n1, n2
      !********************************************************************************************************************        
      !START
      array=(/(i,i=1,n)/) 
      z=0
      do m=1,n,1
         if (array(m)/=0) then
            z=z+1  !z walks the rows of freq 
            if (Selector(m,1)) then 
               freq(z,1)=1
            else   
               freq(z,1)=0 
            endif
            if (Selector(m,2)) then 
               freq(z,2)=1
            else   
               freq(z,2)=0 
            endif 
            !--- 
            do i=m+1,n,1 
               if (array(i)/=0) then   
                  if (d(array(m),array(i))==0) then
                     array(i)=0   !deactivate this individual 
                     if (Selector(i,1)) then  !If i is in group 1
                        freq(z,1)=freq(z,1)+1
                     endif  
                     if (Selector(i,2)) then  !If i is in group 2
                        freq(z,2)=freq(z,2)+1
                     endif
                  endif   
               endif
            enddo
         endif   
      enddo
      n1=SUM(freq(:,1))
      n2=SUM(freq(:,2))
      freq(:,1)=freq(:,1)/n1
      freq(:,2)=freq(:,2)/n2
          
      return
      
   END SUBROUTINE shared_frequencies
       
       
   
   REAL(8) FUNCTION GenDiver(n,freq)
   
       implicit none
       
       !Externals
       integer, intent(in):: n                  !Sample size, the array is set to this lenght cuz is the max # of haplotypes present
       real, dimension(n), intent(in):: freq    !Array with the haplotype frequencies
       !Internals
       integer:: i
       real(8):: x
       !------------------------------
       x=0
       do i=1,n,1
          x=x+freq(i)**2 
       enddo    
       GenDiver=1-x
      
       return
   
   END FUNCTION GenDiver
   
   
   
   REAL(8) FUNCTION fst(a,b,g,n,GDivA,freq)
   
       implicit none
       
       !External variables
       integer, intent(in):: a
       integer, intent(in):: b
       integer, intent(in):: g
       integer, intent(in):: n
       real(8), dimension(g), intent(in):: GDivA
       real, dimension(n),intent(in):: freq
       !Internal variables
       real(8):: Hs, Ht
       !-------------------------
       !Start
       !------------
       Hs = 0.5*(GDivA(a)+GDivA(b))
       Ht = GenDiver(n,freq)*real(n,8)/real(n-1,8)
       if (Ht==0.0) then
          fst=-1000000000
       else    
          fst=1-(Hs/Ht)  !=Ht-Hs/Ht
       endif
     
       return   
   
   END FUNCTION fst
    
   
   
   FUNCTION TajimasD(n,pi,s)
   
       implicit none
       
       !External variables
       integer, intent(in):: n
       real(8), intent(in):: pi
       integer, intent(in):: s
       real(8):: TajimasD
       !Internal variables
       integer:: i
       real(8):: a1, a2, b1, b2, c1, c2, e1, e2
       !-------------------------
       !Start
       !------------
       a1=0
       do i=1,n-1
          a1=a1+1/dble(i)
       enddo
       !------------
       a2=0
       do i=1,n-1
          a2=a2+1/dble(i**2)
       enddo 
       !------------
       b1 = (n+1)/dble(3*(n-1)) 
       !------------
       b2 = 2*dble(n**2+n+3)/dble((9*n)*(n-1))
       !------------
       c1 = b1-a1**(-1)
       !------------
       c2 = b2 - dble(n+2)/(a1*dble(n)) + a2/a1**2
       !------------
       e1 = c1/a1
       !------------
       e2 = c2/(a1**2+a2)
       !Finnaly :
       if (S==0) then
          TajimasD = 0.0
       else   
          TajimasD = ( pi-dble(S)/a1 )/(e1*dble(S)+e2*dble(S)*(dble(S)-1))**0.5
       endif
       
       
       return
   
   END FUNCTION TajimasD
    
    
    
   FUNCTION Fs(n,k0,pi)
       !This function computes the statistic Fs of Fu:
       !Fu, Y-X (1997) Genetics 147:915-925.
   
       implicit none
       
       !External variables
       integer, intent(in):: n   !sample size       
       integer, intent(in):: k0  !Number of observed haplotypes
       real(8), intent(in):: pi  !Average number of pairwise differences
       real(8):: Fs
       !Internal variables
       integer:: i, j
       real(8):: Sp, Sn, Sk     !Correspond to Sum(|Sk|theta_k_pi), Sn(theta), S'
       real(8):: lSn            ! ln(Sn(theta))
       real(8):: x, y
       integer, allocatable, dimension(:,:):: S_table
       !-------------------------
       !Start
       !------------
       if (k0<=1) then
          Fs=1000000000       !Should be Infinity but we can't deal with that shit
       elseif (pi==0) then
          Fs=-1000000000      !Should be -Infinity, see comment above
       else             
          if (n<20) then 
             Sn=1  
             do i=0,n-1,1
                Sn=Sn*(pi+i)            
             enddo   ! Sn=theta(theta+1)...(theta+n-1): Original Fu's (1997) paper is wrong!       
             !First we create the table of Stirling numbers 
             allocate(S_table(n,n))
             call StirlingN(n,S_table)
             !Now we get Sum |Sk|theta**k
             Sk=0
             do j=k0,n,1   !For all k>k0 (alleles number) 
                x=ABS(S_table(n,j))  
                y=pi**j
                Sk=Sk+x*y       
             enddo
             Sp=Sk/Sn
             !---------------------
             !Alternate option if 1-Sp is too small
             if ((1-Sp)<0.1) then
                Sk=0
                do j=1,k0-1,1   !For all k>k0 (alleles number) 
                   x=ABS(S_table(n,j))  
                   y=pi**j
                   Sk=Sk+x*y       
                enddo
                Sp=Sk/Sn 
                Fs=log(1-Sp)-log(Sp) 
             else
                Fs=log(Sp)-log(1-Sp)  !Fs=logit(Sp)    
             endif
             deallocate(S_table)
          else  
             lSn=0  
             do i=0,n-1,1
                lSn=lSn+log(pi+i)            
             enddo
             !Now we get Sum |Sk|theta**k
             Sp=0
             do i=k0,n,1
                Sp = Sp + exp(lnStirling(n,i) + i*log(pi) - lSn)
             enddo
             !---------------------
             !Alternate option if 1-Sp is too small
             if ((1-Sp)<0.1) then
                Sp=0
                do i=1,k0-1,1
                   Sp = Sp + exp(lnStirling(n,i) + i*log(pi) - lSn)
                enddo
                Fs=log(1-Sp)-log(Sp) 
             else
                Fs=log(Sp)-log(1-Sp)  !Fs=logit(Sp)    
             endif
          endif 
       endif
       
       
              
       return
   
   END FUNCTION Fs
    
           
   
   subroutine StirlingN(n,table)
      !This function creates the table of Stirling numbers up to n
      !Rows correspond to values of n, and columns to values of k
      implicit none
      
      !External variables
      integer, intent(in):: n
      integer, dimension(n,n), intent(out):: table  !table represents rows=n, columns=k
      !Internal variables
      integer:: i, j
      integer, allocatable, dimension(:,:):: itable  !table represents rows=n, columns=k
      !---------
      allocate(itable(n+1,n+1))
      itable=0
      do i=1,n+1,1
         itable(i,i)=1    
      enddo    
      do i=2,n,1     !i walks on n (rows), corresponding values [2,n] in rows [3,n+1]
         do j=1,n-1,1     !j walks on k (columns), corresponding values [1,n-1] in coluns [2,n]
            itable(i+1,j+1) = (i-1)*itable(i,j+1) + itable(i,j)   !Formula: s(n+1,k)=n*s(n,k)+s(n,k-1)    
         enddo
      enddo   
      !Finally pouring out the values in the table
      table=0
      do i=1,n,1
         do j=1,n,1
            table(i,j)=itable(i+1,j+1)    
         enddo
      enddo   
      deallocate(itable)
      
      
      return
   
   end subroutine StirlingN
   
   
   
   function lnStirling(n,k)
      !Returns the natural; logarithm of a Stirling number of first kind   
      !Requires functions lngamma, digamma, and trigamma
      use irandom
   
      implicit none
      
      integer, intent(in):: n
      integer, intent(in):: k
      real(8):: lnStirling
      !
      real(8):: ni, ki
      real(8):: x, x0, t0, B, gt, lbin
      !------------------------------
      if ((n<0).or.(k<0).or.(n<k)) then
         write(*,*) 'n and k for Stirling numbers should be positive' 
         stop     
      elseif ((n==0).and.(k==0)) then
         lnStirling=0
      elseif ((n>0).and.(k==0)) then
         write(*,*) 'n and k should be positive to get log(Stirling numbers)' 
         stop 
      elseif ((n>=1).and.(k==1)) then
         x=n 
         lnStirling=lngamma(x) 
      elseif (n==k) then
         lnStirling=0 
      elseif ((n>1).and.(k>1).and.(n>k)) then
         ni = real(n-1,8) 
         ki = real(k-1,8)
         call phiprimeroot(ni,ki,0.1D00,ni*ki,0.0000001D00,x0)   !Here we get the root of phi' (digamma(x+n+1)-digamma(x+1)-k/x)
         t0 = ki/(ni-ki)
         B = lngamma(x0+ni+1) - lngamma(x0+1) - ki*log(x0) -ni*log(t0+1) + ki*log(t0)         
         gt = 1/x0 * sqrt(  ki*(ni-ki)/ni/( trigamma(x0+ni+1) - trigamma(x0+1) + ki/x0/x0 )  ) 
         lbin = lngamma(ni+1.D00) - lngamma(ki+1.D00) - lngamma(ni-ki+1.D00) !This is log of the binomial coefficient(ni,ki)
         lnStirling = B + log(gt) + lbin  
      endif
      
      
      return
   
   end function lnStirling


   
   
   subroutine phiprimeroot(n,m,lo,up,eps,x) 
   
      implicit none
      
      real(8), intent(in):: n
      real(8), intent(in):: m
      real(8), intent(in):: lo
      real(8), intent(in):: up
      real(8), intent(in):: eps
      real(8), intent(out):: x
      
      real(8):: a, b, c, fA, fB, fC, slope
      logical:: done
      integer:: k
      !-------------------------------
      a=lo
      b=up        
      fA=phi(a,n,m)
      fB=phi(b,n,m)      
      if (((fA>0).and.(fB>0)).or.((fA<0).and.(fB<0))) then !The zero is not in the interval
         if (fA>fB) then
            x=b
         else
            x=a
         endif
      else   
         done=.false. 
         k=1
         do
            if ((fA>-eps).and.(fA<eps)) then !if phi(a) is in [0.0-eps,0.0+eps]
               x=a  
               done=.true.
            elseif ((fB>-eps).and.(fB<eps)) then  
               x=b
               done=.true.
            else
               if (k<15) then 
                  c=(a+b)/2
               else 
                  slope=(fB-fA)/(b-a) !Approach by a stright line
                  c=a-fA/slope
               endif
               fC=phi(c,n,m)
               if ((fC>-eps).and.(fC<eps)) then
                  x=c 
                  done=.true.    
               elseif ((fA<0.0).and.(fB>0.0)) then
                  if (fC>0.0) then
                     b=c
                     fB=fC
                  else
                     a=c 
                     fA=fC
                  endif
               elseif ((fA>0.0).and.(fB<0.0)) then    
                  if (fC>0.0) then
                     a=c
                     fA=fC
                  else
                     b=c 
                     fB=fC
                  endif  
               endif   
            endif    
            if (done) exit
            k=k+1
         enddo  
      endif
      
      
      return
      
   end subroutine phiprimeroot
   
   
   
   
   real(8) function phi(x,n,m)
   
      implicit none
      
      real(8), intent(in):: x
      real(8), intent(in):: n
      real(8), intent(in):: m
      
      integer:: iflag
      !------------------------
      phi=digamma(x+n+1.D00,iflag)-digamma(x+1.D00,iflag)-m/x
   
      
      return
   
   end function
   
   
   
   
   function digamma ( x, ifault )

!*****************************************************************************80
!
!! DIGAMMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2016
!
!  Author:
!
!    Original FORTRAN77 version by Jose Bernardo.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jose Bernardo,
!    Algorithm AS 103:
!    Psi ( Digamma ) Function,
!    Applied Statistics,
!    Volume 25, Number 3, 1976, pages 315-317.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the digamma function.
!    0 < X.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = 8 ) DIGAMMA, the value of the digamma function at X.
!
  implicit none

  real ( kind = 8 ), parameter :: c = 8.5D+00
  real ( kind = 8 ), parameter :: euler_mascheroni = 0.57721566490153286060D+00
  real ( kind = 8 ) digamma
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
!
!  Check the input.
!
  if ( x <= 0.0D+00 ) then
    digamma = 0.0D+00
    ifault = 1
    return
  end if
!
!  Initialize.
!
  ifault = 0
!
!  Approximation for small argument.
!
  if ( x <= 0.000001D+00 ) then
    digamma = - euler_mascheroni - 1.0D+00 / x + 1.6449340668482264365D+00 * x
    return
  end if
!
!  Reduce to DIGAMA(X + N).
!
  digamma = 0.0D+00
  x2 = x

  do while ( x2 < c )
    digamma = digamma - 1.0D+00 / x2
    x2 = x2 + 1.0D+00
  end do
!
!  Use Stirling's (actually de Moivre's) expansion.
!
  r = 1.0D+00 / x2

  digamma = digamma + log ( x2 ) - 0.5D+00 * r

  r = r * r

  digamma = digamma &
    - r * ( 1.0D+00 / 12.0D+00 &
    - r * ( 1.0D+00 / 120.0D+00 &
    - r * ( 1.0D+00 / 252.0D+00 &
    - r * ( 1.0D+00 / 240.0D+00 &
    - r * ( 1.0D+00 / 132.0D+00 ) ) ) ) )

  return
end
subroutine psi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! PSI_VALUES returns some values of the Psi or Digamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[x]
!
!    or
!
!      PolyGamma[0,x]
!
!    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
!
!    PSI(1) = -Euler's constant.
!
!    PSI(X+1) = PSI(X) + 1 / X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.5772156649015329D+00, &
    -0.4237549404110768D+00, &
    -0.2890398965921883D+00, &
    -0.1691908888667997D+00, &
    -0.6138454458511615D-01, &
     0.3648997397857652D-01, &
     0.1260474527734763D+00, &
     0.2085478748734940D+00, &
     0.2849914332938615D+00, &
     0.3561841611640597D+00, &
     0.4227843350984671D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end






   function trigamma(x)
   
      implicit none
      
      real(8), intent(in):: x
      real(8):: trigamma
      
      real(8):: y 
      !-------------------------
      y = 5/(66*(x**11)) - 691/(2730*(x**13)) + 7/(6*(x**15))
      trigamma = 1/x + 1/(2*(x**2)) + 1/(6*(x**3)) - 1/(30*(x**5)) + 1/(42*(x**7)) - 1/(30*(x**9))
      
      return
      
   end function trigamma


   
   
   
   

  subroutine NJtree(n,nb,alignF,groups,nr,g,distype,NJts) 
   
    implicit none
      
    !External Variables
    integer, intent(in):: n                                !Sample size
    integer, intent(in):: nb                               !Number of nucleotides in the analyzed fragment
    integer, dimension(n,nb), intent(in):: alignF          !Full alignment 
    integer, dimension(2,n), intent(in):: groups           !A two rows array with the information of the groups of samples whose summary statistics are gonna be calculated upon 
    integer, intent(in):: nr                               !Maximum length of a tree-containing rows (trees are stored one-tree-per-row)      
    integer, intent(in):: g                                !Number of statistical groups 
    integer, intent(in):: distype                          !Type of distance:1.nr. of differences, 2.p-distance, 3.Jukes-Cantor, 4.Kimura 2-parameters
    real(8), intent(out), dimension(g,nr):: NJts           !Neighbor-Joining tree(s), stored in g rows each with a  single tree 
    !Internal Variables
    integer:: h, i, j, k, c, a, b                          !Counters and multipurpose
    integer:: z5                                           !Counters for cells walked in row arrays
    real(8), allocatable, dimension(:):: row5              !Vector of summary statistics (Neighbor-Joining tree)
    integer:: m                                            !Sample size of a specific statistical gorup
    integer, allocatable, dimension(:,:):: dij, Sij        !Array with the distances between sequences (distance is the number of different nucleotides)
    integer, allocatable, dimension(:,:):: trij            !Arrays with the porportions of transitions and transversions (respectively) between seuqences
    integer, allocatable, dimension(:,:):: minidij, miniSij   !Same as dij but only for the group
    integer, allocatable, dimension(:,:):: minitrij           !Same as dij but only for the transverison of the group
    real(8), allocatable, dimension(:,:):: Dm                 !Distances matrix (only for a statistical group))
    logical, dimension(n):: selector                          !Indicates what samples are included for calculating certain summary statistics and whose are not
    real(8), allocatable, dimension(:,:):: unrootree          !The NJ tree: 1st col are where it goes, 2nd branch lenght (the position in the array is the identifier)
    !----------------------------------------------------------------------------------------------------------------
    !  S  T  A  R  T  
    !---------------------------------------------------------------------------------------------------------------- 
    allocate(dij(n,n),Sij(n,n))
    !---
    NJts=0
    !---
    !MAKING THE DISTANCES MATRIX 
    dij=0
    Sij=0
    if (distype<4) then         
      do j=2,n,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do i=1,j-1,1
          if ((groups(1,i)==groups(1,j)).or.(groups(1,i)==groups(2,j))&
            &.or.(groups(2,i)==groups(1,j)).or.(groups(2,i)==groups(2,j))) then 
            do k=1,nb,1  
               if ((alignF(i,k)/=0).and.(alignF(j,k)/=0)) then !No contribution of unknown (n/?) sites 
                  Sij(i,j)=Sij(i,j)+1
                  if (alignF(i,k)/=alignF(j,k)) then 
                     dij(i,j)=dij(i,j)+1
                  endif   
               endif
            enddo
          endif
        enddo    
      enddo
      do i=1,n-1,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do j=i+1,n,1
          Sij(j,i)=Sij(i,j)
          dij(j,i)=dij(i,j) 
        enddo
      enddo  
    elseif (distype==4) then  
      allocate(trij(n,n)) 
      trij=0
      do j=2,n,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do i=1,j-1,1
          if ((groups(1,i)==groups(1,j)).or.(groups(1,i)==groups(2,j))&
          &.or.(groups(2,i)==groups(1,j)).or.(groups(2,i)==groups(2,j))) then 
            do k=1,nb,1  
               if ((alignF(i,k)/=0).and.(alignF(j,k)/=0)) then !No contribution of unknown (n/?) sites 
                  Sij(i,j)=Sij(i,j)+1
                  if ((alignF(i,k)*alignF(j,k)==3).or.(alignF(i,k)*alignF(j,k)==8)) then  !!Here we count A<->G + C<->T differences (transitional)
                     dij(i,j)=dij(i,j)+1
                  elseif (alignF(i,k)/=alignF(j,k)) then   !if the sites are different (and are not transitions) => transversional diff.
                     trij(i,j)=trij(i,j)+1   
                  endif   
               endif
            enddo
          endif
        enddo    
      enddo
      do i=1,n-1,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do j=i+1,n,1
          Sij(j,i)=Sij(i,j)
          dij(j,i)=dij(i,j) 
          trij(j,i)=trij(i,j) 
        enddo
      enddo
    endif
    !----------------------------------------------------------------------------------------------------------------
    !MAKING THE TREES OF EACH GROUP      
    do h=1,g,1          
       !First, making Selector (the array with the info for including/excluding samples in minidij)
       !------------------------------------------------------------------------------------------
       Selector=.false.
       do i=1,n,1
          if ( (groups(1,i)==h) .or. (groups(2,i)==h) ) then 
             Selector(i)=.true.
          endif
       enddo  
       !--------------------------------------
       m=COUNT(Selector)   
       !--------------------------------------   
       if (m<3) then
         call ErrorMessage(1026)  
       endif
       !Second, making mini distances matrix 
       !------------------------------------------------------------------------------------------ 
       allocate(minidij(m,m),miniSij(m,m),minitrij(m,m),Dm(m,m))
       minidij=0
       miniSij=0
       a=1
       b=2  
       do i=1,n-1,1
          do j=i+1,n,1
             if ((Selector(i)).and.(Selector(j))) then
                minidij(a,b)=dij(i,j)
                minidij(b,a)=dij(i,j)
                miniSij(a,b)=Sij(i,j)
                miniSij(b,a)=Sij(i,j)
                if (distype==4) then                     
                   minitrij(a,b)=trij(i,j)
                   minitrij(b,a)=trij(i,j)
                endif 
                if (b<m) then
                   b=b+1
                elseif (b==m) then !If we reached the end of the row, then.. 
                   a=a+1           !..we start a new row  
                   b=a+1           !..and the column is the b+1 
                endif    
             endif
          enddo 
       enddo 
       !---------------
       if (distype<4) then
          call Get_distance(distype,m,minidij,minitrij,miniSij,Dm)
       elseif (distype==4) then
          call Get_distance(distype,m,minidij,minitrij,miniSij,Dm)    
       endif    
       !---------------
       deallocate(minidij,minitrij,miniSij)    
       !----------------------------------------------------------
       !14) NEIGHBOR JOININ-------------------------------
       c=4*m-4 
       allocate(unrootree(2*m-2,2),row5(c)) 
       !+++++++++++++++++++++++++++
       call NJ2(m,Dm,unrootree)
       !call NJ(m,Dm,unrootree)
       !+++++++++++++++++++++++++++         
       z5=1 
       do i=1,2*m-2,1
          row5(z5)=unrootree(i,1)   
          z5=z5+1
       enddo
       do i=1,2*m-2,1
          row5(z5)=unrootree(i,2)   
          z5=z5+1
       enddo
       !---
       do i=1,c,1
          NJts(h,i)=row5(i)
       enddo
       !---
       deallocate(Dm,unrootree,row5)  
       !-------------------------------------------------
       !----------------------------------------------------------
    enddo
    deallocate(dij,Sij)
    if (distype==4) then  
       deallocate(trij)
    endif 
 
          
      
    return   
   
  end subroutine NJtree

  
  
  
  
  subroutine NJDistanceMatrix(n,nb,alignF,groups,distype,dij,Sij,trij) 
   
    implicit none
      
    !External Variables
    integer, intent(in):: n                           !Sample size
    integer, intent(in):: nb                          !Number of nucleotides in the analyzed fragment
    integer, dimension(n,nb), intent(in):: alignF     !Full alignment 
    integer, dimension(2,n), intent(in):: groups      !A two rows array with the information of the groups of samples whose summary statistics are gonna be calculated upon 
    integer, intent(in):: distype                     !Type of distance:1.nr. of differences, 2.p-distance, 3.Jukes-Cantor, 4.Kimura 2-parameters
    integer, dimension(n,n), intent(inout):: dij      !Array with the distances between sequences 
    integer, dimension(n,n), intent(inout):: Sij      !Array with the distances between sequences 
    integer, dimension(n,n), intent(inout):: trij     !Arrays with the porportions of transitions and transversions between sequences
    !Internal Variables
    integer:: i, j, k                                 !Counters and multipurpose
    integer, dimension(n,n):: idij, iSij, itrij       !dij, Sij and trij for internal use
    !----------------------------------------------------------------------------------------------------------------
    !  S  T  A  R  T  
    !---------------------------------------------------------------------------------------------------------------- 
    !MAKING THE DISTANCES MATRIX 
    idij=0
    iSij=0
    if (distype<4) then         
      do j=2,n,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do i=1,j-1,1
          if ((groups(1,i)==groups(1,j)).or.(groups(1,i)==groups(2,j))&
            &.or.(groups(2,i)==groups(1,j)).or.(groups(2,i)==groups(2,j))) then 
            do k=1,nb,1  
               if ((alignF(i,k)/=0).and.(alignF(j,k)/=0)) then !No contribution of unknown (n/?) sites 
                  iSij(i,j)=iSij(i,j)+1
                  if (alignF(i,k)/=alignF(j,k)) then 
                     idij(i,j)=idij(i,j)+1
                  endif   
               endif
            enddo
          endif
        enddo    
      enddo
      do i=1,n-1,1      !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do j=i+1,n,1
          iSij(j,i)=iSij(i,j)
          idij(j,i)=idij(i,j) 
        enddo
      enddo  
    elseif (distype==4) then  
      itrij=0
      do j=2,n,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do i=1,j-1,1
          if ((groups(1,i)==groups(1,j)).or.(groups(1,i)==groups(2,j))&
          &.or.(groups(2,i)==groups(1,j)).or.(groups(2,i)==groups(2,j))) then 
            do k=1,nb,1  
               if ((alignF(i,k)/=0).and.(alignF(j,k)/=0)) then !No contribution of unknown (n/?) sites 
                  iSij(i,j)=iSij(i,j)+1
                  if ((alignF(i,k)*alignF(j,k)==3).or.(alignF(i,k)*alignF(j,k)==8)) then  !!Here we count A<->G + C<->T differences (transitional)
                     idij(i,j)=idij(i,j)+1
                  elseif (alignF(i,k)/=alignF(j,k)) then   !if the sites are different (and are not transitions) => transversional diff.
                     itrij(i,j)=itrij(i,j)+1   
                  endif   
               endif
            enddo
          endif
        enddo    
      enddo
      do i=1,n-1,1       !This double "do" assures that all the combinations of i,j with i=/j are probed only once
        do j=i+1,n,1
          iSij(j,i)=iSij(i,j)
          idij(j,i)=idij(i,j) 
          itrij(j,i)=itrij(i,j) 
        enddo
      enddo
    endif
    !Done! ->output is dij, Sij, trij
    dij = dij + idij
    Sij = Sij + iSij
    trij = trij + itrij
    
    
    return   
   
  end subroutine NJDistanceMatrix
  
  
  
  
  
  subroutine NJOnlyTree(n,groups,nr,g,distype,dij,Sij,trij,NJts) 
   
    implicit none
      
    !External Variables
    integer, intent(in):: n                                !Sample size
    integer, dimension(2,n), intent(in):: groups           !A two rows array with the information of the groups of samples whose summary statistics are gonna be calculated upon 
    integer, intent(in):: nr                               !Maximum length of a tree-containing rows (trees are stored one-tree-per-row)      
    integer, intent(in):: g                                !Number of statistical groups 
    integer, intent(in):: distype                          !Type of distance:1.nr. of differences, 2.p-distance, 3.Jukes-Cantor, 4.Kimura 2-parameters
    integer, dimension(n,n), intent(in):: dij              !Array with the distances between sequences 
    integer, dimension(n,n), intent(in):: Sij              !Array with the distances between sequences 
    integer, dimension(n,n), intent(in):: trij             !Arrays with the porportions of transitions and transversions between sequences    
    real(8), intent(out), dimension(g,nr):: NJts           !Neighbor-Joining tree(s), stored in g rows each with a  single tree 
    !Internal Variables
    integer:: h, i, j, c, a, b                             !Counters and multipurpose
    integer:: z5                                           !Counters for cells walked in row arrays
    real(8), allocatable, dimension(:):: row5              !Vector of summary statistics (Neighbor-Joining tree)
    integer:: m                                            !Sample size of a specific statistical gorup
    integer, allocatable, dimension(:,:):: minidij, miniSij   !Same as dij but only for the group
    integer, allocatable, dimension(:,:):: minitrij           !Same as dij but only for the transverison of the group
    real(8), allocatable, dimension(:,:):: Dm                 !Distances matrix (only for a statistical group))
    logical, dimension(n):: selector                          !Indicates what samples are included for calculating certain summary statistics and whose are not
    real(8), allocatable, dimension(:,:):: unrootree          !The NJ tree: 1st col are where it goes, 2nd branch lenght (the position in the array is the identifier)
    !----------------------------------------------------------------------------------------------------------------
    !  S  T  A  R  T  
    !---------------------------------------------------------------------------------------------------------------- 
    !MAKING THE TREES OF EACH GROUP      
    do h=1,g,1          
       !First, making Selector (the array with the info for including/excluding samples in minidij)
       !------------------------------------------------------------------------------------------
       Selector=.false.
       do i=1,n,1
          if ( (groups(1,i)==h) .or. (groups(2,i)==h) ) then 
             Selector(i)=.true.
          endif
       enddo  
       !--------------------------------------
       m=COUNT(Selector)   
       !--------------------------------------   
       if (m<3) then
         call ErrorMessage(1026)  
       endif
       !Second, making mini distances matrix 
       !------------------------------------------------------------------------------------------ 
       allocate(minidij(m,m),miniSij(m,m),minitrij(m,m),Dm(m,m))
       minidij=0
       miniSij=0
       a=1
       b=2  
       do i=1,n-1,1
          do j=i+1,n,1
             if ((Selector(i)).and.(Selector(j))) then
                minidij(a,b)=dij(i,j)
                minidij(b,a)=dij(i,j)
                miniSij(a,b)=Sij(i,j)
                miniSij(b,a)=Sij(i,j)
                if (distype==4) then                     
                   minitrij(a,b)=trij(i,j)
                   minitrij(b,a)=trij(i,j)
                endif 
                if (b<m) then
                   b=b+1
                elseif (b==m) then !If we reached the end of the row, then.. 
                   a=a+1           !..we start a new row  
                   b=a+1           !..and the column is the b+1 
                endif    
             endif
          enddo 
       enddo 
       !---------------
       if (distype<4) then
          call Get_distance(distype,m,minidij,minitrij,miniSij,Dm)
       elseif (distype==4) then
          call Get_distance(distype,m,minidij,minitrij,miniSij,Dm)    
       endif    
       !---------------
       deallocate(minidij,minitrij,miniSij)    
       !----------------------------------------------------------
       !14) NEIGHBOR JOININ-------------------------------
       c=4*m-4 
       allocate(unrootree(2*m-2,2),row5(c)) 
       !+++++++++++++++++++++++++++
       call NJ2(m,Dm,unrootree)
       !call NJ(m,Dm,unrootree)
       !+++++++++++++++++++++++++++         
       z5=1 
       do i=1,2*m-2,1
          row5(z5)=unrootree(i,1)   
          z5=z5+1
       enddo
       do i=1,2*m-2,1
          row5(z5)=unrootree(i,2)   
          z5=z5+1
       enddo
       !---
       do i=1,c,1
          NJts(h,i)=row5(i)
       enddo
       !---
       deallocate(Dm,unrootree,row5)  
       !-------------------------------------------------
       !----------------------------------------------------------
    enddo
    
 
          
      
    return   
   
  end subroutine NJonlyTree
  
  
   
   
   
   subroutine Get_distance(distype,m,dij,trij,Sij,Dm)
   
      implicit none
      
      !Externals
      integer, intent(in):: distype               !Type of distance:1.# of differences, 2.p-distance, 3.Jukes-Cantor, 4.Kimura 2-p
      integer, intent(in):: m                     !Sample size (matrixes size)
      integer, dimension(m,m), intent(in):: dij   !Number of differences (transitions if distype=4)
      integer, dimension(m,m), intent(in):: trij  !Number of transversional differences
      integer, dimension(m,m), intent(in):: Sij   !Number of effective nucleotides in each comparison
      real(8), dimension(m,m), intent(out):: Dm   !Distances matrix
      !Internals
      integer:: i, j
      real(8):: p, q
      !-----------------------------------------------------------------
      !Start
      Dm=0
      select case(distype)
      case(1)
         do j=2,m,1       !This assures that all the combinations of i,j with i=/j are probed only once
            do i=1,j-1,1 
               Dm(i,j)=dble(dij(i,j)) 
            enddo
         enddo
      case(2)
         do j=2,m,1       !This assures that all the combinations of i,j with i=/j are probed only once
            do i=1,j-1,1 
               Dm(i,j)=dble(dij(i,j))/dble(Sij(i,j)) 
            enddo
         enddo
      case(3)
         do j=2,m,1       !This assures that all the combinations of i,j with i=/j are probed only once
            do i=1,j-1,1   
               Dm(i,j)=-0.75*log(1.0-1.33333333333333*(dble(dij(i,j))/dble(Sij(i,j)))) 
            enddo
         enddo   
      case(4)         
         do j=2,m,1       !This assures that all the combinations of i,j with i=/j are probed only once
            do i=1,j-1,1 
               p=dble(dij(i,j))/dble(Sij(i,j))
               q=dble(trij(i,j))/dble(Sij(i,j)) 
               Dm(i,j)=-0.5*log(1-2*p-q)-0.25*log(1-2*q)
            enddo
         enddo 
      end select
      do j=2,m,1      
         do i=1,j-1,1             
            Dm(j,i)=Dm(i,j) 
         enddo
      enddo
      
      
      return
   
   endsubroutine Get_distance
   
   
   
   
   subroutine NJ(n,dij,unrootree)
 
        implicit none

        !External
        integer, intent(in):: n                                 ! Number of taxa
        real(8), dimension(n,n), intent(in):: dij               ! Matrix of distances among taxa
        real(8), dimension(n+n-2,2), intent(out):: unrootree    ! The NJ tree: 1st col are where it goes, 2nd branch lenght (the position in the array is the identifier)
        !Internals 
        integer:: i, j, k, m, x
        integer:: a, b, c, d
        real(8):: S1, L1x,L2x, d1, d2, dist
        integer, dimension(n):: id
        real(8), allocatable, dimension(:,:):: idij, tdij    
        logical, allocatable, dimension(:):: Selector
        logical:: switch
        !--------------------------------------------------------
        !START
        unrootree=0
        !----------------------------------------------------------------
        !(1) Getting the lenght of the star tree
        S1=SUM(dij)/(2*(n-1)) 
        !----------------------------------------------------------------
        !(2) Joining the first neighbor
        allocate(idij(n,n),tdij(n,n))
        idij=real(dij)
        tdij=idij
        id=(/(j,j=1,n)/)
        m=n
        do i=1,n-2,1
           call Get_Neighbors(m,idij,a,b,L1x,L2x)
           !Retrieving the obtained info into the tree
           unrootree(id(a),1)=n+i
           unrootree(id(b),1)=n+i
           unrootree(id(a),2)=L1x
           unrootree(id(b),2)=L2x
           !--
           id(a)=n+i
           do j=b,n-1,1
              id(j)=id(j+1)    
           enddo
           !------------------
           deallocate(idij)
           if (allocated(Selector)) then
              deallocate(Selector)
           endif 
           allocate(idij(n-i,n-i),Selector(n-i+1))
           !Now calculating again idij with a-b (neighbors) as one taxa
           Selector=.true.
           Selector(b)=.false. !Position a is gonna carry the info of a-b
           !--
           tdij(a,:)=tdij(a,:)-L1x  
           tdij(:,a)=tdij(:,a)-L1x  
           !--
           idij=0
           c=1
           d=2
           do j=1,m-1,1
              do k=j+1,m,1
                 if ((Selector(j)).and.(Selector(k))) then
                    idij(c,d)=tdij(j,k)
                    idij(d,c)=idij(c,d)
                    if (d<m-1) then       !This instructions reset the row
                       d=d+1              !Remember the row has m-2 because matrix size will be m-1
                    elseif (d==m-1) then  !and we are leaving one reserved for the composite neighbors
                       c=c+1           
                       d=c+1           
                    endif    
                 endif
              enddo 
           enddo 
           deallocate(tdij)
           allocate(tdij(n-i,n-i))
           tdij=idij
           m=m-1
        enddo 
        !Finally we check if there is one remaining
        !At this point it is possible that one node is not conected with anyone else, if so it goes to the (2n-2)-th
        i=1       !So this do-cycle gets the node that has its go-to empty
        switch=.false. 
        do 
           if (unrootree(i,1)==0) then !We locate the neighbor without node assigned (so that it get the last node assigned)
              switch=.true.   
           endif        
           if (switch) exit
           i=i+1 
        enddo
        x=i    !x gets the position of the first taxa/node with a go-to empty
        if (x/=2*n-2) then !If such node is the last one, then the tree was already finished, if not we go on..
            a=0
            b=0
            i=1
            do               !This cycle works by walking the places of tree and each one of them is used as starting point to  
              !------------  !..start a go-to cycle until reaching the last node (2n-2)-th or the x-th node, once we have one of each
              j=i            !stored in positions a & b, we leave
              dist=0
              do 
                 dist=dist+unrootree(j,2)
                 if ((j==2*n-2).or.(j==x)) exit
                 j=NINT(unrootree(j,1))
              enddo  
              if (j==2*n-2) then
                  a=i
                  d1=dist
              elseif (j==x) then
                  b=i
                  d2=dist
              endif    
              !------------  
              if ((a/=0).and.(b/=0)) exit
              i=i+1     
            enddo    
            !a & b are taxa/nodes whose go-to path leads to the node (2n-2)-th and x-th respectively    
            unrootree(x,1)=2*n-2     
            unrootree(x,2)=dij(a,b)-d1-d2
        endif 
    
    
        return
    

   end subroutine NJ
   
   
   
   
   subroutine NJ2(n,dij,unrootree)
   !This subroutine estimate a phylogenetic tree by emplying the 
   !Neighbor-Joining method of Saitou& Nei (1987) but with the modified
   !algorithm of Studier & Keppler (1988). It can be consulted in
   !Felsestein's book "Inferring Phylogenies"
 
     implicit none

     !External
     integer, intent(in):: n                                 ! Number of taxa
     real(8), dimension(n,n), intent(in):: dij               ! Matrix of distances among taxa
     real(8), dimension(n+n-2,2), intent(out):: unrootree    ! The NJ tree: 1st col are go-to, 2nd branch lenght (the position in the array is the identifier)
     !Internals 
     integer:: i, j, k, l, r                                 !Counters
     integer:: m, ai, aj                                     !m is the in-cycle n, ai, aj store i, j positions of chosen taxa
     real(8), allocatable, dimension(:,:):: idij, Ndij       !dij for in-cycle use
     real(8), allocatable, dimension(:):: u
     real(8):: vi, vj, duj, mduj
     integer, allocatable, dimension(:):: taxa, ntaxa        !Keeps track of the taxa (including new-nodes-taxa) as labels of dij
     !--------------------------------------------------------
     !START
     allocate(idij(n+1,n+1),Ndij(n+1,n+1),u(n),taxa(n),ntaxa(n))
     !---------
     u=0
     idij=0
     do j=1,n,1
        do i=1,n,1    
           idij(i,j)=dij(i,j)
        enddo
     enddo
     unrootree=0
     taxa=(/(i,i=1,n)/)
     m=n
     do r=1,n-2,1   
        !----------------------------------------------------------------
        !(1) We compute ui for each taxon 
        do i=1,m,1
           u(i)=SUM(idij(:,i))/(m-2)
        enddo   
        !----------------------------------------------------------------
        !(2) We found the taxa pair that minimizes Dij-u(i)-u(j)
        do i=1,m-1,1
           do j=i+1,m,1                
              duj=idij(i,j)-u(i)-u(j)
              if ((duj<mduj).or.((i==1).and.(j==2))) then
                 mduj=duj    
                 ai=i
                 aj=j
              endif   
           enddo
        enddo    
        !----------------------------------------------------------------
        !(3) Join items i & j in new node x:
        !Compute branch length from x to i:
        vi=0.5*idij(ai,aj)+0.5*(u(ai)-u(aj))
        !Compute branch length from x to j:
        vj=0.5*idij(ai,aj)+0.5*(u(aj)-u(ai))
        !Make the changes
        unrootree(taxa(ai),1)=n+r     !n+r is the id of the new node
        unrootree(taxa(aj),1)=n+r    
        unrootree(taxa(ai),2)=MAX(vi,dble(0.0))  !MAX is to prevent negative branches
        unrootree(taxa(aj),2)=MAX(vj,dble(0.0))
        !----------------------------------------------------------------
        !(4) Compute the distance between new node and remaining taxa           
        do i=1,m,1
           if ((i/=ai).and.(i/=aj)) then  
              idij(i,m+1)=(idij(ai,i)+idij(aj,i)-idij(ai,aj))/2
           endif
        enddo
        do i=1,m,1           
           idij(m+1,i)=idij(i,m+1)
        enddo
        !----------------------------------------------------------------        
        !(5) Delete tips i & j and reeplace by new-node
        ntaxa=0
        j=1
        do i=1,m,1
           if ((i/=ai).and.(i/=aj)) then 
              ntaxa(j)=taxa(i)
              j=j+1
           endif
        enddo
        ntaxa(m-1)=n+r    !n+k is the id of the new node
        !Collapse the columns/rows of the chosen taxa
        Ndij=0           
        l=1
        do i=1,m+1,1
           if ((i/=ai).and.(i/=aj)) then 
              k=1 
              do j=1,m+1,1 
                 if ((j/=ai).and.(j/=aj)) then 
                    Ndij(k,l)=idij(i,j)
                    k=k+1
                 endif
              enddo
              l=l+1
           endif
        enddo                   
        !----------------------------------------------------------------        
        !(6) Finish it, if only 2-nodes remain
        m=m-1
        if (m==2) then
           unrootree(ntaxa(1),1)=ntaxa(2)
           unrootree(ntaxa(1),2)=Ndij(1,2)
           unrootree(ntaxa(2),1)=ntaxa(2)
           unrootree(ntaxa(2),2)=0              
        else              
           idij=0
           do i=1,m,1
              do j=1,m,1
                 idij(i,j)=Ndij(i,j)
              enddo
           enddo    
           taxa=ntaxa
        endif    
        !----------------------------------------------------------------
     enddo
     
     deallocate(idij,Ndij,u,taxa,ntaxa)  
     
      
     return    

   end subroutine NJ2
   
   
   
   
   subroutine Get_Neighbors(n,dij,a,b,L1x,L2x)
 
        implicit none

        !External
        integer, intent(in):: n                         ! Number of taxa
        real(8), dimension(n,n), intent(in):: dij          ! Matrix of distances among taxa
        integer, intent(out):: a                        ! One of the two neighbors found
        integer, intent(out):: b                        ! The other of the two neighbors found
        real(8), intent(out):: L1x                      ! Lenght of the first branch
        real(8), intent(out):: L2x                      ! Lenght of the second branch
        !Internals 
        integer:: i, j, z
        real(8):: lenght, L1, L2                        ! Minimum lenght of the trees 
        real(8), dimension(n*(n-1)/2):: L               ! Stores the lenghts all the alternative trees

        !--------------------------------------------------------
        !START
        z=1
        do i=1,n-1,1
           do j=i+1,n,1
              call tree_branch_lenghts(n,i,j,dij,L(z),L1,L2)
              if (z>1) then
                 if (L(z)<lenght) then
                    a=i
                    b=j
                    L1x=L1
                    L2x=L2
                    lenght=L(z)
                 endif 
              else
                 a=i
                 b=j
                 L1x=L1
                 L2x=L2
                 lenght=L(z)  
              endif   
              z=z+1 
           enddo
        enddo
     

        return

   end subroutine Get_Neighbors

    

   subroutine tree_branch_lenghts(n,x,y,dij,LT,L1x,L2x)


        implicit none


        !External
        integer, intent(in):: n                    ! Number of taxa
        integer, intent(in):: x                    ! First neighbor taxa
        integer, intent(in):: y                    ! Second neighbor taxa
        real(8), dimension(:,:), intent(in):: dij  ! Matrix of distances among taxa
        real(8), intent(out):: LT                  ! Lenght of the tree
        real(8), intent(out):: L1x                 ! Lenght of the first branch (taxa to node)
        real(8), intent(out):: L2x                 ! Lenght of the second branch (taxa to node)
        !Internals
        real(8):: Liy, D1z, D2z

        !-------------------------------------------
    
        LT = ( SUM(dij(x,:)) + SUM(dij(y,:)) - 2*dij(x,y) ) * (0.5/(real(n)-2))  !This is the sum of all paths x-z, y-z (excluding x-y)
        LT = LT + dij(x,y)/2
        Liy = SUM(dij)/2 - SUM(dij(x,:)) - SUM(dij(y,:)) + dij(x,y)  !This is the sum of all distances not including x or y
        LT = LT + Liy/(real(n)-2)
        D1z = ( SUM(dij(x,:)) - dij(x,y) ) / (real(n)-2) 
        D2z = ( SUM(dij(y,:)) - dij(x,y) ) / (real(n)-2)
        L1x = (dij(x,y)+D1z-D2z)/2
        L2x = (dij(x,y)+D2z-D1z)/2


        return
 
   end subroutine tree_branch_lenghts
   
   
   
   
   integer function iMINMAXLOC1(m,n,iarray,x1,x2,l1,l2) 
   
      implicit none
   
      !This function finds the location of the maximum/minimum value (depending if m=1/-1)) of   
      ! an array (iarray). We can define an interval of meaningful positions [l1,l2] and an
      ! interval of values [x1,x2]-including x1 & x2-, so if the value of the cell is outside
      ! [x1,x2] then the cell is ignored
   
      ! THE DIFFERENCE OF iMINMAX AND rMINMAX IS THAT THE ARRAY IS INTEGER IN iMINMAX AND REAL IN rMINMAX
      ! THE DIFFERENCE BETWEEN MINMAX1 AND MINMAX2 IS IN x1,x2 BEING EXCLUDED OR NOT FROM THE INTERVAL  
   
      integer, intent(in):: m                      
      integer, intent(in):: n
      integer, dimension(n), intent(in):: iarray
      integer, intent(in):: x1, x2
      integer, intent(in):: l1, l2
      
      integer:: i
      real(8):: ref
      logical:: first
      !start
      iMINMAXLOC1=0
      first=.true.      
      do i=l1,l2,1
         if ((iarray(i)>=x1).AND.(iarray(i)<=x2)) then 
            if (m==1) then 
               if (first) then
                  ref=iarray(i)   
                  iMINMAXLOC1=i   
                  first=.false.
               else   
                  if (iarray(i)>ref) then  
                     ref=iarray(i)   
                     iMINMAXLOC1=i                  
                  endif
               endif   
            elseif (m==-1) then
               if (first) then 
                  ref=iarray(i)   
                  iMINMAXLOC1=i  
                  first=.false.
               else    
                  if (iarray(i)<ref) then
                     ref=iarray(i)   
                     iMINMAXLOC1=i
                  endif 
               endif   
            endif    
         endif
      enddo

      
      return    
      
   end function iMINMAXLOC1
   
   
   
   integer function rMINMAXLOC1(m,n,iarray,x1,x2,l1,l2) 
   
      implicit none
   
      !This function finds the location of the maximum/minimum value (depending if m=1/-1)) of   
      ! an array (iarray). We can define an interval of meaningful positions [l1,l2] and an
      ! interval of values [x1,x2]-including x1 & x2-, so if the value of the cell is outside
      ! [x1,x2] then the cell is ignored
      
      ! THE DIFFERENCE OF iMINMAX AND rMINMAX IS THAT THE ARRAY IS INTEGER IN iMINMAX AND REAL IN rMINMAX
      ! THE DIFFERENCE BETWEEN MINMAX1 AND MINMAX2 IS IN x1,x2 BEING EXCLUDED OR NOT FROM THE INTERVAL  
     
      integer, intent(in):: m                      
      integer, intent(in):: n
      real(8), dimension(n), intent(in):: iarray
      real(8), intent(in):: x1, x2
      integer, intent(in):: l1, l2
      
      integer:: i
      real(8):: ref
      logical:: first
      !start
      rMINMAXLOC1=0
      first=.true.      
      do i=l1,l2,1
         if ((iarray(i)>=x1).AND.(iarray(i)<=x2)) then 
            if (m==1) then 
               if (first) then
                  ref=iarray(i)   
                  rMINMAXLOC1=i   
                  first=.false.
               else   
                  if (iarray(i)>ref) then  
                     ref=iarray(i)   
                     rMINMAXLOC1=i                  
                  endif
               endif   
            elseif (m==-1) then
               if (first) then 
                  ref=iarray(i)   
                  rMINMAXLOC1=i  
                  first=.false.
               else    
                  if (iarray(i)<ref) then
                     ref=iarray(i)   
                     rMINMAXLOC1=i
                  endif 
               endif   
            endif    
         endif
      enddo

      
      return    
      
   end function rMINMAXLOC1
   
   
   
   integer function rMINMAXLOC2(m,n,iarray,x1,x2,l1,l2) 
   
      implicit none
      
      ! This function finds the location of the maximum/minimum value (depending if m=1/-1)) of   
      ! an array (iarray) of size n. We can define a meaningful interval in the array [l1,l2] and an
      ! interval of values (x1,x2)-noninclusive x1 & x2-, so if the value of the cell is outside
      ! (x1,x2) then the cell is ignored
      
      ! THE DIFFERENCE OF iMINMAX AND rMINMAX IS THAT THE ARRAY IS INTEGER IN iMINMAX AND REAL IN rMINMAX
      ! THE DIFFERENCE BETWEEN MINMAX1 AND MINMAX2 IS IN x1,x2 BEING EXCLUDED OR NOT FROM THE INTERVAL  
   
      integer, intent(in):: m                      
      integer, intent(in):: n
      real(8), dimension(n), intent(in):: iarray
      real(8), intent(in):: x1, x2
      integer, intent(in):: l1, l2
      
      integer:: i
      real(8):: ref
      logical:: first
      !start
      rMINMAXLOC2=0
      first=.true.      
      do i=l1,l2,1
         if ((iarray(i)>x1).AND.(iarray(i)<x2)) then 
            if (m==1) then 
               if (first) then
                  ref=iarray(i)   
                  rMINMAXLOC2=i   
                  first=.false.
               else   
                  if (iarray(i)>ref) then  
                     ref=iarray(i)   
                     rMINMAXLOC2=i                  
                  endif
               endif   
            elseif (m==-1) then
               if (first) then 
                  ref=iarray(i)   
                  rMINMAXLOC2=i  
                  first=.false.
               else    
                  if (iarray(i)<ref) then
                     ref=iarray(i)   
                     rMINMAXLOC2=i
                  endif 
               endif   
            endif    
         endif
      enddo

      
      return    
      
   end function rMINMAXLOC2
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   real(8) function calendar_date(c14date,switch)
   
        implicit none
        
        real(8), intent(in):: c14date    !14C-age BP
        logical, intent(in):: switch
        !real, intent(out):: calendar_date
        
        integer:: i  
        real:: x1,x2,y1,y2
        !---------------------------------        
        if (switch) then
           i=NINT(c14date)
           Select Case (i)
           Case(:56)
               calendar_date=0             
           Case(57)
               calendar_date=19
           Case(58)
               calendar_date=21   
           Case(59:60)
               calendar_date=25
           Case(61:63)
               calendar_date=26             
           Case(64:66)
               calendar_date=27
           Case(67:69)
               calendar_date=28             
           Case(70)
               calendar_date=29
           Case(71)
               calendar_date=30             
           Case(72)
               calendar_date=54
           Case(73)
               calendar_date=55             
           Case(74)
               calendar_date=58
           Case(75)
               calendar_date=88     
           Case(76)
               calendar_date=89   
           Case(77)
               calendar_date=90   
           Case(78)
               calendar_date=91  
           Case(79)
               calendar_date=93   
           Case(80)
               calendar_date=142  
           Case(81)
               calendar_date=143   
           Case(82:83)
               calendar_date=144  
           Case(84:88)
               calendar_date=145   
           Case(89:92)
               calendar_date=146   
           Case(93:94)
               calendar_date=147   
           Case(95)
               calendar_date=148   
           Case(96)
               calendar_date=149   
           Case(97)
               calendar_date=151   
           Case(98)
               calendar_date=151   
           Case(99:101)
               calendar_date=152
           Case(102:107)
               calendar_date=153
           Case(108:111)
               calendar_date=152             
           Case(112:114)
               calendar_date=151
           Case(115:119)
               calendar_date=150             
           Case(120)
               calendar_date=149
           Case(121)
               calendar_date=148             
           Case(122)
               calendar_date=137 
           Case(123:125)
               calendar_date=133                 
           Case(126)
               calendar_date=132                
           Case(127:131)
               calendar_date=131
           Case(132)
               calendar_date=142               
           Case(133:143)
               calendar_date=143               
           Case(144:146)
               calendar_date=142             
           Case(147:152)
               calendar_date=143
           Case(153:154)
               calendar_date=142 
           Case(155)
               calendar_date=123
           Case(156:159)
               calendar_date=122
           Case(160)
               calendar_date=123
           Case(161)
               calendar_date=124
           Case(162)
               calendar_date=184    
           Case(163)
               calendar_date=207
           Case(164)
               calendar_date=208
           Case(165)
               calendar_date=209
           Case(166:168)
               calendar_date=210
           Case(169)
               calendar_date=211    
           Case(170:175)
               calendar_date=212
           Case(176:178)
               calendar_date=146
           Case(179)
               calendar_date=145 
           Case(180)
               calendar_date=146
           Case(181)
               calendar_date=145    
           Case(182:183)
               calendar_date=146
           Case(184)
               calendar_date=145
           Case(185:188)
               calendar_date=146
           Case(189:193)
               calendar_date=147
           Case(194:195)
               calendar_date=148    
           Case(196)
               calendar_date=219
           Case(197)
               calendar_date=220
           Case(198:199)
               calendar_date=221  
           Case(200)
               calendar_date=222
           Case(201:202)
               calendar_date=223    
           Case(203:205)
               calendar_date=224
           Case(206:207)
               calendar_date=225
           Case(208:210)
               calendar_date=226  
           Case(211:212)
               calendar_date=227
           Case(213:215)
               calendar_date=228    
           Case(216:217)
               calendar_date=229
           Case(218:219)
               calendar_date=230
           Case(220)
               calendar_date=231
           Case(221)
               calendar_date=232
           Case(222)
               calendar_date=292    
           Case(223)
               calendar_date=293
           Case(224:225)
               calendar_date=294
           Case(226:228)
               calendar_date=295
           Case(229:231)
               calendar_date=296
           Case(232:234)
               calendar_date=297    
           Case(235:237)
               calendar_date=298
           Case(238:240)
               calendar_date=299
           Case(241:244)
               calendar_date=300
           Case(245:248)
               calendar_date=301
           Case(249:251)
               calendar_date=302    
           Case(252:255)
               calendar_date=303
           Case(256:258)
               calendar_date=304
           Case(259:261)
               calendar_date=305 
           Case(262:264)
               calendar_date=306                 
           Case(265:269)
               calendar_date=307 
           Case(270:271)
               calendar_date=308     
           Case(272:274)
               calendar_date=309
           Case(275:278)
               calendar_date=310 
           Case(279:281)
               calendar_date=311 
           Case(282)
               calendar_date=312 
           Case(283)
               calendar_date=314 
           Case(284)
               calendar_date=359
           Case(285)
               calendar_date=361
           Case(286)
               calendar_date=362  
           Case(287)
               calendar_date=363 
           Case(288:290)
               calendar_date=364 
           Case(291:292)
               calendar_date=365 
           Case(293:294)
               calendar_date=366
           Case(295:296)
               calendar_date=367 
           Case(297:300)
               calendar_date=368 
           Case(301)
               calendar_date=369
           Case(302:304)
               calendar_date=370  
           Case(305)
               calendar_date=371 
           Case(306:308)
               calendar_date=372 
           Case(309:310)
               calendar_date=373 
           Case(311)
               calendar_date=374
           Case(312:313)
               calendar_date=375 
           Case(314:316)
               calendar_date=376 
           Case(317:320)
               calendar_date=377
           Case(321:324)
               calendar_date=378  
           Case(325:326)
               calendar_date=379 
           Case(327:328)
               calendar_date=380 
           Case(329:330)
               calendar_date=381 
           Case(331)
               calendar_date=382
           Case(332:333)
               calendar_date=383 
           Case(334)
               calendar_date=384
           Case(335)
               calendar_date=385
           Case(336)
               calendar_date=386  
           Case(337:338)
               calendar_date=387 
           Case(339)
               calendar_date=388 
           Case(340)
               calendar_date=389 
           Case(341:342)
               calendar_date=390
           Case(343)
               calendar_date=392 
           Case(344)
               calendar_date=393 
           Case(345)
               calendar_date=394
           Case(346:347)
               calendar_date=396  
           Case(348)
               calendar_date=397 
           Case(349)
               calendar_date=398 
           Case(350)
               calendar_date=399 
           Case(351)
               calendar_date=400
           Case(352)
               calendar_date=401 
           Case(353)
               calendar_date=402                
           Case(354)
               calendar_date=403
           Case(355:356)
               calendar_date=404 
           Case(357)
               calendar_date=406 
           Case(358)
               calendar_date=407 
           Case(359)
               calendar_date=408   
           Case(360)
               calendar_date=409
           Case(361)
               calendar_date=410 
           Case(362)
               calendar_date=412 
           Case(363)
               calendar_date=413 
           Case(364)
               calendar_date=414
           Case(365)
               calendar_date=415
           Case(366)
               calendar_date=416 
           Case(367)
               calendar_date=418 
           Case(368)
               calendar_date=460 
           Case(369)
               calendar_date=464   
           Case(370)
               calendar_date=466
           Case(371)
               calendar_date=469 
           Case(372)
               calendar_date=471 
           Case(373)
               calendar_date=473 
           Case(374)
               calendar_date=475 
           Case(375)
               calendar_date=476
           Case(376)
               calendar_date=477 
           Case(377)
               calendar_date=478 
           Case(378:379)
               calendar_date=479 
           Case(380)
               calendar_date=481   
           Case(381:382)
               calendar_date=482
           Case(383)
               calendar_date=483 
           Case(384:386)
               calendar_date=484 
           Case(387)
               calendar_date=486 
           Case(388)
               calendar_date=487  
           Case(389:390)
               calendar_date=488
           Case(391)
               calendar_date=489 
           Case(392)
               calendar_date=490 
           Case(393)
               calendar_date=491 
           Case(394:395)
               calendar_date=492   
           Case(396:397)
               calendar_date=493
           Case(398)
               calendar_date=494 
           Case(399:400)
               calendar_date=495 
           Case(401:402)
               calendar_date=496 
           Case(403:404)
               calendar_date=497 
           Case(405:406)
               calendar_date=498
           Case(407:408)
               calendar_date=499
           Case(409:410)
               calendar_date=500
           Case(411:413)
               calendar_date=501
           Case(414:416)
               calendar_date=502
           Case(417:418)
               calendar_date=503 
           Case(419:420)
               calendar_date=504
           Case(421:422)
               calendar_date=505
           Case(423:426)
               calendar_date=506
           Case(427:429)
               calendar_date=507
           Case(430:432)
               calendar_date=508
           Case(433:435)
               calendar_date=509
           Case(436:438)
               calendar_date=510
           Case(439:442)
               calendar_date=511
           Case(443:446)
               calendar_date=512
           Case(447:449)
               calendar_date=513
           Case(450:453)
               calendar_date=514
           Case(454:457)
               calendar_date=515
           Case(458:460)
               calendar_date=516
           Case(461:464)
               calendar_date=517
           Case(465:468)
               calendar_date=518
           Case(469:473)
               calendar_date=519
           Case(474:476)
               calendar_date=520
           Case(477:480)
               calendar_date=521
           Case(481:482)
               calendar_date=522
           Case(483:485)
               calendar_date=523
           Case(486:487)
               calendar_date=524
           Case(488:490)
               calendar_date=525
           Case(491:493)
               calendar_date=526
           Case(494:497)
               calendar_date=527  
           Case(498:499)
               calendar_date=528
           Case(500:502)
               calendar_date=529
           Case(503:505)
               calendar_date=530
           Case(506:510)
               y1=531
               y2=532
               x1=506
               x2=510
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(511:515)
               y1=533
               y2=534
               x1=511
               x2=515
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(516:520)
               y1=535
               y2=536
               x1=516
               x2=520
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(521:525)
               y1=537
               y2=538
               x1=521
               x2=525
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(526:530)
               y1=539
               y2=540
               x1=526
               x2=530
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(531:535)
               y1=541
               y2=543
               x1=531
               x2=535
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(536:540)
               y1=543
               y2=544
               x1=536
               x2=540
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(541:545)
               y1=545
               y2=546
               x1=541
               x2=545
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(546:550)
               y1=547
               y2=548
               x1=546
               x2=550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(551:555)
               y1=548
               y2=550
               x1=551
               x2=555
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(556:560)
               y1=550
               y2=581
               x1=556
               x2=560
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(561:565)
               y1=582
               y2=585
               x1=561
               x2=565
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(566:570)
               y1=586
               y2=588
               x1=566
               x2=570
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(571:575)
               y1=588
               y2=591
               x1=571
               x2=575
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(576:580)
               y1=591
               y2=593
               x1=576
               x2=580
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(581:585)
               y1=593
               y2=595
               x1=581
               x2=585
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(586:590)
               y1=595
               y2=597
               x1=586
               x2=590
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(591:595)
               y1=597
               y2=598
               x1=591
               x2=595
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(596:600)
               y1=599
               y2=600
               x1=596
               x2=600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(601:605)
               y1=600
               y2=602
               x1=601
               x2=605
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(606:610)
               y1=602
               y2=604
               x1=606
               x2=610
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(611:615)
               y1=604
               y2=605
               x1=611
               x2=615
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(616:620)
               y1=606
               y2=607
               x1=616
               x2=620
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(621:625)
               y1=608
               y2=609
               x1=621
               x2=625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(626:630)
               y1=609
               y2=611
               x1=626
               x2=630
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(631:635)
               y1=611
               y2=612
               x1=631
               x2=635
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(636:640)
               y1=613
               y2=614
               x1=636
               x2=640
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(641:645)
               y1=614
               y2=615
               x1=641
               x2=645
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(646:650)
               y1=615
               y2=617
               x1=646
               x2=650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(651:655)
               y1=617
               y2=618
               x1=651
               x2=655
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(656:660)
               y1=619
               y2=620
               x1=656
               x2=660
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(661:665)
               y1=620
               y2=622
               x1=661
               x2=665
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(666:670)
               y1=623
               y2=661
               x1=666
               x2=670
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(671:675)
               y1=662
               y2=664
               x1=671
               x2=675
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(676:680)
               y1=664
               y2=666
               x1=676
               x2=680
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(681:685)
               y1=666
               y2=667
               x1=681
               x2=685
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(686:690)
               y1=667
               y2=668
               x1=686
               x2=690
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(691:695)
               y1=668
               y2=669
               x1=691
               x2=695
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(696:700)
               y1=669
               y2=670
               x1=696
               x2=700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(701:705)
               y1=670
               y2=671
               x1=701
               x2=705
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(706:710)
               y1=672
               y2=673
               x1=706
               x2=710
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(711:715)
               y1=673
               y2=674
               x1=711
               x2=715
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(716:720)
               y1=674
               y2=675
               x1=716
               x2=720
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(721:725)
               y1=676
               y2=677
               x1=721
               x2=725
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(726:730)
               y1=677
               y2=678
               x1=726
               x2=730
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(731:735)
               y1=678
               y2=679
               x1=731
               x2=735
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(736:740)
               y1=680
               y2=681
               x1=736
               x2=740
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(741:745)
               y1=681
               y2=683
               x1=741
               x2=745
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(746:750)
               y1=683
               y2=684
               x1=746
               x2=750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(751:755)
               y1=685
               y2=686
               x1=751
               x2=755
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(756:760)
               y1=687
               y2=688
               x1=756
               x2=760
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(761:765)
               y1=688
               y2=690
               x1=761
               x2=765
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(766:770)
               y1=690
               y2=692
               x1=766
               x2=770
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(771:775)
               y1=693
               y2=695
               x1=771
               x2=775
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(776:780)
               y1=696
               y2=702
               x1=776
               x2=780
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(781:785)
               y1=704
               y2=707
               x1=781
               x2=785
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(786:790)
               y1=707
               y2=710
               x1=786
               x2=790
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(791:795)
               y1=710
               y2=713
               x1=791
               x2=795
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(796:800)
               y1=713
               y2=715
               x1=796
               x2=800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(801:805)
               y1=716
               y2=718
               x1=801
               x2=805
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(806:810)
               y1=718
               y2=720
               x1=806
               x2=810
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(811:815)
               y1=721
               y2=724
               x1=811
               x2=815
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(816:820)
               y1=724
               y2=729
               x1=816
               x2=820
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(821:825)
               y1=731
               y2=736
               x1=821
               x2=825
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(826:830)
               y1=736
               y2=739
               x1=826
               x2=830
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(831:835)
               y1=739
               y2=742
               x1=831
               x2=835
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(836:840)
               y1=743
               y2=747
               x1=836
               x2=840
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(841:845)
               y1=749
               y2=754
               x1=841
               x2=845
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(846:850)
               y1=755
               y2=758
               x1=846
               x2=850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(851:855)
               y1=759
               y2=762
               x1=851
               x2=855
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(856:860)
               y1=762
               y2=765
               x1=856
               x2=860
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(861:865)
               y1=766
               y2=768
               x1=861
               x2=865
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(866:870)
               y1=769
               y2=771
               x1=866
               x2=870
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(871:875)
               y1=772
               y2=775
               x1=871
               x2=875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(876:880)
               y1=776
               y2=779
               x1=876
               x2=880
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(881:885)
               y1=781
               y2=785
               x1=881
               x2=885
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(886:890)
               y1=787
               y2=831
               x1=886
               x2=890
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(891:895)
               y1=834
               y2=839
               x1=891
               x2=895
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(896:900)
               y1=840
               y2=843
               x1=896
               x2=900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(901:905)
               y1=843
               y2=845
               x1=901
               x2=905
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(906:910)
               y1=846
               y2=849
               x1=906
               x2=910
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(911:915)
               y1=850
               y2=852
               x1=911
               x2=915
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(916:920)
               y1=852
               y2=853
               x1=916
               x2=920
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(921:925)
               y1=854
               y2=855
               x1=921
               x2=925
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(926:930)
               y1=855
               y2=857
               x1=926
               x2=930
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(931:935)
               y1=858
               y2=861
               x1=931
               x2=935
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(936:940)
               y1=863
               y2=869
               x1=936
               x2=940
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(941:945)
               y1=870
               y2=874
               x1=941
               x2=945
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(946:950)
               y1=875
               y2=877
               x1=946
               x2=950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(951:955)
               y1=878
               y2=880
               x1=951
               x2=955
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(956:960)
               y1=881
               y2=883
               x1=956
               x2=960
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(961:965)
               y1=883
               y2=887
               x1=961
               x2=965
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(966:970)
               y1=889
               y2=921
               x1=966
               x2=970
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(971:975)
               y1=922
               y2=924
               x1=971
               x2=975
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(976:980)
               y1=924
               y2=925
               x1=976
               x2=980
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(981:985)
               y1=926
               y2=927
               x1=981
               x2=985
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(986:990)
               y1=927
               y2=928
               x1=986
               x2=990
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(991:995)
               y1=929
               y2=930
               x1=991
               x2=995
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(996:1000)
               y1=930
               y2=931
               x1=996
               x2=1000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1001:1010)
               y1=931
               y2=935
               x1=1001
               x2=1010
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(1011:1020)
               y1=935
               y2=939
               x1=1011
               x2=1020
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1021:1030)
               y1=940
               y2=945
               x1=1021
               x2=1030
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1031:1040)
               y1=946
               y2=950
               x1=1031
               x2=1040
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1041:1050)
               y1=950
               y2=954
               x1=1041
               x2=1050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1051:1060)
               y1=955
               y2=961
               x1=1051
               x2=1060
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1061:1070)
               y1=962
               y2=968
               x1=1061
               x2=1071
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1071:1080)
               y1=968
               y2=974
               x1=1071
               x2=1080
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1081:1090)
               y1=975
               y2=1007
               x1=1081
               x2=1090
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1091:1100)
               y1=1009
               y2=1015
               x1=1091
               x2=1100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(1101:1110)
               y1=1015
               y2=1020
               x1=1101
               x2=1110
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1111:1120)
               y1=1020
               y2=1025
               x1=1111
               x2=1120
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1121:1130)
               y1=1026
               y2=1033
               x1=1121
               x2=1130
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1131:1140)
               y1=1034
               y2=1039
               x1=1131
               x2=1140
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1141:1150)
               y1=1040
               y2=1046
               x1=1141
               x2=1150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1151:1160)
               y1=1048
               y2=1069
               x1=1151
               x2=1160
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1161:1170)
               y1=1069
               y2=1090
               x1=1161
               x2=1170
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1171:1180)
               y1=1093
               y2=1103
               x1=1171
               x2=1180
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1181:1190)
               y1=1104
               y2=1115
               x1=1181
               x2=1190
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1191:1200)
               y1=1116
               y2=1128
               x1=1191
               x2=1200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(1201:1210)
               y1=1129
               y2=1137
               x1=1201
               x2=1210
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1211:1220)
               y1=1138
               y2=1159
               x1=1211
               x2=1220
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1221:1230)
               y1=1160
               y2=1194
               x1=1221
               x2=1230
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1231:1240)
               y1=1197
               y2=1205
               x1=1231
               x2=1240
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1241:1250)
               y1=1205
               y2=1211
               x1=1241
               x2=1250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1251:1260)
               y1=1212
               y2=1219
               x1=1251
               x2=1260
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1261:1270)
               y1=1219
               y2=1224
               x1=1261
               x2=1270
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1271:1280)
               y1=1225
               y2=1229
               x1=1271
               x2=1280
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1281:1290)
               y1=1230
               y2=1235
               x1=1281
               x2=1290
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1291:1300)
               y1=1236
               y2=1271
               x1=1291
               x2=1300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
            Case(1301:1310)
               y1=1272
               y2=1277
               x1=1301
               x2=1310
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1311:1320)
               y1=1277
               y2=1282
               x1=1311
               x2=1320
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1321:1330)
               y1=1282
               y2=1286
               x1=1321
               x2=1330
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1331:1340)
               y1=1286
               y2=1289
               x1=1331
               x2=1340
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1341:1350)
               y1=1289
               y2=1292
               x1=1341
               x2=1350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1351:1360)
               y1=1293
               y2=1295
               x1=1351
               x2=1360
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1361:1370)
               y1=1296
               y2=1298
               x1=1361
               x2=1370
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1371:1380)
               y1=1298
               y2=1301
               x1=1371
               x2=1380
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1381:1390)
               y1=1301
               y2=1304
               x1=1381
               x2=1390
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1391:1400)
               y1=1304
               y2=1307
               x1=1391
               x2=1400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1401:1410)
               y1=1307
               y2=1311
               x1=1401
               x2=1410
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1411:1420)
               y1=1311
               y2=1318
               x1=1411
               x2=1420
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1421:1430)
               y1=1319
               y2=1325
               x1=1421
               x2=1430
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1431:1440)
               y1=1326
               y2=1331
               x1=1431
               x2=1440
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1441:1450)
               y1=1332
               y2=1338
               x1=1441
               x2=1450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1451:1460)
               y1=1339
               y2=1349
               x1=1451
               x2=1460
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1461:1470)
               y1=1350
               y2=1359
               x1=1461
               x2=1470
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1471:1480)
               y1=1360
               y2=1369
               x1=1471
               x2=1480
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1481:1490)
               y1=1369
               y2=1376
               x1=1481
               x2=1490
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1491:1500)
               y1=1377
               y2=1386
               x1=1491
               x2=1500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1501:1510)
               y1=1386
               y2=1393
               x1=1501
               x2=1510
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1511:1520)
               y1=1394
               y2=1400
               x1=1511
               x2=1520
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1521:1530)
               y1=1401
               y2=1406
               x1=1521
               x2=1530
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1531:1540)
               y1=1407
               y2=1412
               x1=1531
               x2=1540
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1541:1550)
               y1=1414
               y2=1454
               x1=1541
               x2=1550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1551:1560)
               y1=1455
               y2=1462
               x1=1551
               x2=1560
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1561:1570)
               y1=1462
               y2=1469
               x1=1561
               x2=1570
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1571:1580)
               y1=1469
               y2=1471
               x1=1571
               x2=1580
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1581:1590)
               y1=1471
               y2=1475
               x1=1581
               x2=1590
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1591:1600)
               y1=1476
               y2=1483
               x1=1591
               x2=1600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1601:1610)
               y1=1484
               y2=1527
               x1=1601
               x2=1610
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1611:1620)
               y1=1528
               y2=1533
               x1=1611
               x2=1620
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1621:1630)
               y1=1534
               y2=1538
               x1=1621
               x2=1630
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1631:1640)
               y1=1538
               y2=1542
               x1=1631
               x2=1640
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1641:1650)
               y1=1543
               y2=1548
               x1=1641
               x2=1650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1651:1660)
               y1=1549
               y2=1553
               x1=1652
               x2=1660
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1661:1670)
               y1=1554
               y2=1560
               x1=1661
               x2=1670
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1671:1680)
               y1=1561
               y2=1579
               x1=1671
               x2=1680
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1681:1690)
               y1=1580
               y2=1586
               x1=1681
               x2=1690
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1691:1700)
               y1=1587
               y2=1594
               x1=1691
               x2=1700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1701:1710)
               y1=1594
               y2=1633
               x1=1701
               x2=1710
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1711:1720)
               y1=1635
               y2=1652
               x1=1711
               x2=1720
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1721:1730)
               y1=1653
               y2=1658
               x1=1721
               x2=1730
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1731:1740)
               y1=1658
               y2=1662
               x1=1731
               x2=1740
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1741:1750)
               y1=1663
               y2=1667
               x1=1741
               x2=1750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1751:1760)
               y1=1668
               y2=1673
               x1=1751
               x2=1760
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1761:1770)
               y1=1673
               y2=1681
               x1=1761
               x2=1770
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1771:1780)
               y1=1683
               y2=1711
               x1=1771
               x2=1780
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1781:1790)
               y1=1712
               y2=1716
               x1=1781
               x2=1790
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1791:1800)
               y1=1717
               y2=1722
               x1=1791
               x2=1800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
            Case(1801:1810)
               y1=1723
               y2=1747
               x1=1801
               x2=1810
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1811:1820)
               y1=1749
               y2=1763
               x1=1811
               x2=1820
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1821:1830)
               y1=1764
               y2=1773
               x1=1821
               x2=1830
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1831:1840)
               y1=1773
               y2=1780
               x1=1831
               x2=1840
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1841:1850)
               y1=1780
               y2=1787
               x1=1841
               x2=1850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1851:1860)
               y1=1788
               y2=1818
               x1=1851
               x2=1860
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1861:1870)
               y1=1820
               y2=1833
               x1=1861
               x2=1870
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1871:1880)
               y1=1836
               y2=1843
               x1=1871
               x2=1880
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1881:1890)
               y1=1843
               y2=1848
               x1=1881
               x2=1890
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1891:1900)
               y1=1848
               y2=1853
               x1=1891
               x2=1900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1901:1910)
               y1=1854
               y2=1859
               x1=1901
               x2=1910
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1911:1920)
               y1=1860
               y2=1877
               x1=1911
               x2=1920
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(1921:1930)
               y1=1878
               y2=1884
               x1=1921
               x2=1930
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1931:1940)
               y1=1884
               y2=1890
               x1=1931
               x2=1940
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1941:1950)
               y1=1891
               y2=1903
               x1=1941
               x2=1950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1951:1960)
               y1=1904
               y2=1909
               x1=1951
               x2=1960
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1961:1970)
               y1=1910
               y2=1916
               x1=1961
               x2=1970
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1971:1980)
               y1=1917
               y2=1925
               x1=1971
               x2=1980
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1981:1990)
               y1=1926
               y2=1947
               x1=1981
               x2=1990
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(1991:2000)
               y1=1949
               y2=1960
               x1=1991
               x2=2000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2001:2010)
               y1=1961
               y2=1967
               x1=2001
               x2=2010
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2011:2020)
               y1=1968
               y2=1974
               x1=2011
               x2=2020
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(2021:2030)
               y1=1975
               y2=1986
               x1=2021
               x2=2030
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2031:2040)
               y1=1987
               y2=1998
               x1=2031
               x2=2040
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2041:2050)
               y1=1998
               y2=2015
               x1=2041
               x2=2050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2051:2060)
               y1=2016
               y2=2022
               x1=2051
               x2=2060
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2061:2070)
               y1=2023
               y2=2029
               x1=2061
               x2=2070
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2071:2080)
               y1=2029
               y2=2040
               x1=2071
               x2=2080
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2081:2090)
               y1=2043
               y2=2074
               x1=2081
               x2=2090
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2091:2100)
               y1=2076
               y2=2086
               x1=2091
               x2=2100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2101:2110)
               y1=2087
               y2=2094
               x1=2101
               x2=2110
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2111:2120)
               y1=2095
               y2=2101
               x1=2111
               x2=2120
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(2121:2130)
               y1=2103
               y2=2128
               x1=2121
               x2=2130
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2131:2140)
               y1=2129
               y2=2136
               x1=2131
               x2=2140
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2141:2150)
               y1=2137
               y2=2142
               x1=2141
               x2=2150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2151:2160)
               y1=2143
               y2=2213
               x1=2151
               x2=2160
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2161:2170)
               y1=2214
               y2=2221
               x1=2161
               x2=2170
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2171:2180)
               y1=2222
               y2=2226
               x1=2171
               x2=2180
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2181:2190)
               y1=2227
               y2=2231
               x1=2181
               x2=2190
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2191:2200)
               y1=2231
               y2=2234
               x1=2191
               x2=2200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2201:2210)
               y1=2234
               y2=2237
               x1=2201
               x2=2210
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2211:2220)
               y1=2237
               y2=2221
               x1=2211
               x2=2220
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2221:2230)
               y1=2221
               y2=2251
               x1=2221
               x2=2230
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2231:2240)
               y1=2253
               y2=2263
               x1=2231
               x2=2240
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2241:2250)
               y1=2264
               y2=2273
               x1=2241
               x2=2250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2251:2260)
               y1=2274
               y2=2328
               x1=2251
               x2=2260
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2261:2270)
               y1=2329
               y2=2333
               x1=2261
               x2=2270
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2271:2280)
               y1=2334
               y2=2337
               x1=2271
               x2=2280
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2281:2290)
               y1=2337
               y2=2341
               x1=2281
               x2=2290
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2291:2300)
               y1=2341
               y2=2343
               x1=2291
               x2=2300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2301:2310)
               y1=2343
               y2=2345
               x1=2301
               x2=2310
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2311:2320)
               y1=2345
               y2=2347
               x1=2311
               x2=2320
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(2321:2330)
               y1=2347
               y2=2349
               x1=2321
               x2=2330
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2331:2340)
               y1=2349
               y2=2351
               x1=2331
               x2=2340
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2341:2350)
               y1=2351
               y2=2353
               x1=2341
               x2=2350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2351:2360)
               y1=2353
               y2=2355
               x1=2351
               x2=2360
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2361:2370)
               y1=2355
               y2=2358
               x1=2361
               x2=2370
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2371:2380)
               y1=2358
               y2=2360
               x1=2371
               x2=2380
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2381:2390)
               y1=2360
               y2=2364
               x1=2381
               x2=2390
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2391:2400)
               y1=2386
               y2=2400
               x1=2391
               x2=2400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2401:2410)
               y1=2401
               y2=2410
               x1=2401
               x2=2410
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2411:2420)
               y1=2411
               y2=2418
               x1=2411
               x2=2420
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(2421:2430)
               y1=2419
               y2=2425
               x1=2421
               x2=2430
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2431:2440)
               y1=2426
               y2=2525
               x1=2431
               x2=2440
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2441:2450)
               y1=2527
               y2=2573
               x1=2441
               x2=2450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2451:2460)
               y1=2575
               y2=2586
               x1=2451
               x2=2460
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2461:2470)
               y1=2587
               y2=2595
               x1=2461
               x2=2470
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2471:2480)
               y1=2596
               y2=2600
               x1=2471
               x2=2480
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2481:2490)
               y1=2599
               y2=2583
               x1=2481
               x2=2490
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2491:2500)
               y1=2583
               y2=2592
               x1=2491
               x2=2500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(2501:2525)
               y1=2593
               y2=2673
               x1=2501
               x2=2525
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(2526:2550)
               y1=2673
               y2=2730
               x1=2526
               x2=2550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2551:2575)
               y1=2731
               y2=2741
               x1=2551
               x2=2575
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2576:2600)
               y1=2742
               y2=2749
               x1=2576
               x2=2600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2601:2625)
               y1=2749
               y2=2756
               x1=2601
               x2=2625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(2626:2650)
               y1=2756
               y2=2763
               x1=2626
               x2=2650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(2651:2675)
               y1=2763
               y2=2772
               x1=2651
               x2=2675
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2676:2700)
               y1=2773
               y2=2786
               x1=2676
               x2=2700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2701:2725)
               y1=2787
               y2=2821
               x1=2701
               x2=2725
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2726:2750)
               y1=2822
               y2=2839
               x1=2726
               x2=2750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(2751:2775)
               y1=2840
               y2=2870
               x1=2751
               x2=2775
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(2776:2800)
               y1=2870
               y2=2902
               x1=2776
               x2=2800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2801:2825)
               y1=2903
               y2=2926
               x1=2801
               x2=2825
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2826:2850)
               y1=2927
               y2=2964
               x1=2826
               x2=2850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2851:2875)
               y1=2965
               y2=2990
               x1=2851
               x2=2875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(2876:2900)
               y1=2992
               y2=3037
               x1=2876
               x2=2900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(2901:2925)
               y1=3038
               y2=3098
               x1=2901
               x2=2925
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2926:2950)
               y1=3101
               y2=3120
               x1=2926
               x2=2950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2951:2975)
               y1=3120
               y2=3182
               x1=2951
               x2=2975
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(2976:3000)
               y1=3183
               y2=3199
               x1=2976
               x2=3000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
          Case(3001:3025)
               y1=3200
               y2=3241
               x1=3001
               x2=3025
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3026:3050)
               y1=3242
               y2=3294
               x1=3026
               x2=3050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3051:3075)
               y1=3295
               y2=3308
               x1=3051
               x2=3075
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3076:3100)
               y1=3309
               y2=3352
               x1=3076
               x2=3100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3101:3125)
               y1=3352
               y2=3366
               x1=3101
               x2=3125
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(3126:3150)
               y1=3366
               y2=3379
               x1=3126
               x2=3150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3151:3175)
               y1=3380
               y2=3394
               x1=3151
               x2=3175
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3176:3200)
               y1=3395
               y2=3423
               x1=3176
               x2=3200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3201:3225)
               y1=3424
               y2=3441
               x1=3201
               x2=3225
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3226:3250)
               y1=3444
               y2=3465
               x1=3226
               x2=3250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(3251:3275)
               y1=3466
               y2=3511
               x1=3251
               x2=3275
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3276:3300)
               y1=3512
               y2=3525
               x1=3276
               x2=3300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3301:3325)
               y1=3525
               y2=3552
               x1=3301
               x2=3325
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3326:3350)
               y1=3569
               y2=3598
               x1=3326
               x2=3350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3351:3375)
               y1=3599
               y2=3622
               x1=3351
               x2=3375
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(3376:3400)
               y1=3623
               y2=3664
               x1=3376
               x2=3400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3401:3425)
               y1=3665
               y2=3675
               x1=3401
               x2=3425
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3426:3450)
               y1=3675
               y2=3706
               x1=3426
               x2=3450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3451:3475)
               y1=3706
               y2=3766
               x1=3451
               x2=3475
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3476:3500)
               y1=3767
               y2=3774
               x1=3476
               x2=3500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3501:3525)
               y1=3774
               y2=3811
               x1=3501
               x2=3525
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3526:3550)
               y1=3812
               y2=3853
               x1=3526
               x2=3550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3551:3575)
               y1=3854
               y2=3875
               x1=3551
               x2=3575
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3576:3600)
               y1=3876
               y2=3906
               x1=3576
               x2=3600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3601:3625)
               y1=3907
               y2=3944
               x1=3601
               x2=3625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(3626:3650)
               y1=3945
               y2=3962
               x1=3626
               x2=3650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3651:3675)
               y1=3963
               y2=4029
               x1=3651
               x2=3675
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3676:3700)
               y1=4030
               y2=4041
               x1=3676
               x2=3700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3701:3725)
               y1=4042
               y2=4056
               x1=3701
               x2=3725
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3726:3750)
               y1=4073
               y2=4122
               x1=3726
               x2=3750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(3751:3775)
               y1=4122
               y2=4132
               x1=3751
               x2=3775
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3776:3800)
               y1=4132
               y2=4196
               x1=3776
               x2=3800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3801:3825)
               y1=4196
               y2=4211
               x1=3801
               x2=3825
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3826:3850)
               y1=4212
               y2=4266
               x1=3826
               x2=3850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3851:3875)
               y1=4267
               y2=4333
               x1=3851
               x2=3875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(3876:3900)
               y1=4337
               y2=4359
               x1=3876
               x2=3900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(3901:3925)
               y1=4360
               y2=4415
               x1=3901
               x2=3925
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3926:3950)
               y1=4416
               y2=4422
               x1=3926
               x2=3950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3951:3975)
               y1=4423
               y2=4433
               x1=3951
               x2=3975
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(3976:4000)
               y1=4435
               y2=4476
               x1=3976
               x2=4000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(4001:4025)
               y1=4477
               y2=4486
               x1=4001
               x2=4025
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4026:4050)
               y1=4486
               y2=4498
               x1=4026
               x2=4050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4051:4075)
               y1=4499
               y2=4552
               x1=4051
               x2=4075
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4076:4100)
               y1=4553
               y2=4677
               x1=4076
               x2=4100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4101:4125)
               y1=4678
               y2=4699
               x1=4101
               x2=4125
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(4126:4150)
               y1=4701
               y2=4709
               x1=4126
               x2=4150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4151:4175)
               y1=4709
               y2=4765
               x1=4151
               x2=4175
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4176:4200)
               y1=4766
               y2=4786
               x1=4176
               x2=4200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4201:4225)
               y1=4787
               y2=4836
               x1=4201
               x2=4225
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4226:4250)
               y1=4837
               y2=4843
               x1=4226
               x2=4250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(4251:4275)
               y1=4844
               y2=4850
               x1=4251
               x2=4275
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4276:4300)
               y1=4850
               y2=4857
               x1=4276
               x2=4300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4301:4325)
               y1=4858
               y2=4865
               x1=4301
               x2=4325
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4326:4350)
               y1=4865
               y2=4908
               x1=4326
               x2=4350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4351:4375)
               y1=4909
               y2=4923
               x1=4351
               x2=4375
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(4376:4400)
               y1=4924
               y2=4978
               x1=4376
               x2=4400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4401:4425)
               y1=4986
               y2=5011
               x1=4401
               x2=4425
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4426:4450)
               y1=5012
               y2=5027
               x1=4426
               x2=4450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4451:4475)
               y1=5099
               y2=5184
               x1=4451
               x2=4475
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4476:4500)
               y1=5186
               y2=5175
               x1=4476
               x2=4500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4501:4525)
               y1=5174
               y2=5192
               x1=4501
               x2=4525
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4526:4550)
               y1=5193
               y2=5303
               x1=4526
               x2=4550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4551:4575)
               y1=5303
               y2=5312
               x1=4551
               x2=4575
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4576:4600)
               y1=5312
               y2=5317
               x1=4576
               x2=4600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4601:4625)
               y1=5318
               y2=5380
               x1=4601
               x2=4625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(4626:4650)
               y1=5380
               y2=5388
               x1=4626
               x2=4650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4651:4675)
               y1=5388
               y2=5416
               x1=4651
               x2=4675
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4676:4700)
               y1=5415
               y2=5402
               x1=4676
               x2=4700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4701:4725)
               y1=5402
               y2=5520
               x1=4701
               x2=4725
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4726:4750)
               y1=5521
               y2=5529
               x1=4726
               x2=4750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(4751:4775)
               y1=5529
               y2=5514
               x1=4751
               x2=4775
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4776:4800)
               y1=5513
               y2=5545
               x1=4776
               x2=4800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4801:4825)
               y1=5545
               y2=5594
               x1=4801
               x2=4825
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4826:4850)
               y1=5594
               y2=5598
               x1=4826
               x2=4850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4851:4875)
               y1=5598
               y2=5604
               x1=4851
               x2=4875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(4876:4900)
               y1=5605
               y2=5626
               x1=4876
               x2=4900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(4901:4925)
               y1=5627
               y2=5639
               x1=4901
               x2=4925
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4926:4950)
               y1=5639
               y2=5680
               x1=4926
               x2=4950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4951:4975)
               y1=5681
               y2=5695
               x1=4951
               x2=4975
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(4976:5000)
               y1=5696
               y2=5731
               x1=4976
               x2=5000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(5001:5050)
               y1=5732
               y2=5820
               x1=5001
               x2=5050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5051:5100)
               y1=5821
               y2=5845
               x1=5051
               x2=5100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5101:5150)
               y1=5845
               y2=5918
               x1=5101
               x2=5150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5151:5200)
               y1=5918
               y2=5959
               x1=5151
               x2=5200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5201:5250)
               y1=5959
               y2=5980
               x1=5201
               x2=5250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5251:5300)
               y1=5981
               y2=6097
               x1=5251
               x2=5300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5301:5350)
               y1=6097
               y2=6162
               x1=5301
               x2=5350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5351:5400)
               y1=6163
               y2=6240
               x1=5351
               x2=5400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5401:5450)
               y1=6240
               y2=6260
               x1=5401
               x2=5450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5451:5500)
               y1=6260
               y2=6300
               x1=5451
               x2=5500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(5501:5550)
               y1=6300
               y2=6349
               x1=5501
               x2=5550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5551:5600)
               y1=6349
               y2=6375
               x1=5551
               x2=5600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5601:5650)
               y1=6377
               y2=6431
               x1=5601
               x2=5650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5651:5700)
               y1=6431
               y2=6481
               x1=5651
               x2=5700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5701:5750)
               y1=6482
               y2=6535
               x1=5701
               x2=5750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5751:5800)
               y1=6537
               y2=6613
               x1=5751
               x2=5800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5801:5850)
               y1=6614
               y2=6668
               x1=5801
               x2=5850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5851:5900)
               y1=6668
               y2=6716
               x1=5851
               x2=5900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5901:5950)
               y1=6716
               y2=6772
               x1=5901
               x2=5950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(5951:6000)
               y1=6773
               y2=6837
               x1=5951
               x2=6000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1         
           Case(6001:6050)
               y1=6838
               y2=6917
               x1=6001
               x2=6050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6051:6100)
               y1=6918
               y2=6970
               x1=6051
               x2=6100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6101:6150)
               y1=6971
               y2=7080
               x1=6101
               x2=6150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6151:6200)
               y1=7081
               y2=7104
               x1=6151
               x2=6200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6201:6250)
               y1=7105
               y2=7205
               x1=6201
               x2=6250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6251:6300)
               y1=7205
               y2=7229
               x1=6251
               x2=6300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6301:6350)
               y1=7232
               y2=7276
               x1=6301
               x2=6350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6351:6400)
               y1=7277
               y2=7354
               x1=6351
               x2=6400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6401:6450)
               y1=7357
               y2=7382
               x1=6401
               x2=6450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6451:6500)
               y1=7383
               y2=7430
               x1=6451
               x2=6500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(6501:6550)
               y1=7430
               y2=7449
               x1=6501
               x2=6550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6551:6600)
               y1=7450
               y2=7492
               x1=6551
               x2=6600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6601:6650)
               y1=7493
               y2=7542
               x1=6601
               x2=6650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6651:6700)
               y1=7542
               y2=7580
               x1=6651
               x2=6700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6701:6750)
               y1=7580
               y2=7602
               x1=6701
               x2=6750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6751:6800)
               y1=7603
               y2=7646
               x1=6751
               x2=6800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6801:6850)
               y1=7646
               y2=7681
               x1=6801
               x2=6850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6851:6900)
               y1=7681
               y2=7714
               x1=6851
               x2=6900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6901:6950)
               y1=7715
               y2=7782
               x1=6901
               x2=6950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(6951:7000)
               y1=7783
               y2=7846
               x1=6951
               x2=7000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(7001:7050)
               y1=7847
               y2=7902
               x1=7001
               x2=7050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7051:7100)
               y1=7903
               y2=7949
               x1=7051
               x2=7100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7101:7150)
               y1=7949
               y2=7971
               x1=7101
               x2=7150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7151:7200)
               y1=7971
               y2=8008
               x1=7151
               x2=7200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7201:7250)
               y1=8009
               y2=8034
               x1=7201
               x2=7250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7251:7300)
               y1=8035
               y2=8112
               x1=7251
               x2=7300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7301:7350)
               y1=8113
               y2=8180
               x1=7301
               x2=7350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7351:7400)
               y1=8180
               y2=8238
               x1=7351
               x2=7400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7401:7450)
               y1=8238
               y2=8273
               x1=7401
               x2=7450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7451:7500)
               y1=8273
               y2=8354
               x1=7451
               x2=7500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(7501:7550)
               y1=8354
               y2=8383
               x1=7501
               x2=7550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7551:7600)
               y1=8384
               y2=8404
               x1=7551
               x2=7600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7601:7650)
               y1=8405
               y2=8424
               x1=7601
               x2=7650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7651:7700)
               y1=8424
               y2=8481
               x1=7651
               x2=7700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7701:7750)
               y1=8483
               y2=8550
               x1=7701
               x2=7750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7751:7800)
               y1=8551
               y2=8590
               x1=7751
               x2=7800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7801:7850)
               y1=8591
               y2=8620
               x1=7801
               x2=7850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7851:7900)
               y1=8621
               y2=8677
               x1=7851
               x2=7900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7901:7950)
               y1=8678
               y2=8854
               x1=7901
               x2=7950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(7951:8000)
               y1=8857
               y2=8900
               x1=7951
               x2=8000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(8001:8050)
               y1=8901
               y2=9004
               x1=8001
               x2=8050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8051:8100)
               y1=9004
               y2=9019
               x1=8051
               x2=8100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8101:8150)
               y1=9019
               y2=9063
               x1=8101
               x2=8150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8151:8200)
               y1=9064
               y2=9175
               x1=8151
               x2=8200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8201:8250)
               y1=9178
               y2=9218
               x1=8201
               x2=8250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8251:8300)
               y1=9219
               y2=9343
               x1=8251
               x2=8300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8301:8350)
               y1=9344
               y2=9383
               x1=8301
               x2=8350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8351:8400)
               y1=9384
               y2=9456
               x1=8351
               x2=8400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8401:8450)
               y1=9456
               y2=9484
               x1=8401
               x2=8450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8451:8500)
               y1=9484
               y2=9516
               x1=8451
               x2=8500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(8501:8550)
               y1=9517
               y2=9538
               x1=8501
               x2=8550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8551:8600)
               y1=9538
               y2=9547
               x1=8551
               x2=8600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8601:8650)
               y1=9548
               y2=9561
               x1=8601
               x2=8650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8651:8700)
               y1=9562
               y2=9644
               x1=8651
               x2=8700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8701:8750)
               y1=9646
               y2=9721
               x1=8701
               x2=8750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8751:8800)
               y1=9722
               y2=9837
               x1=8751
               x2=8800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8801:8850)
               y1=9838
               y2=10012
               x1=8801
               x2=8850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8851:8900)
               y1=10012
               y2=10052
               x1=8851
               x2=8900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8901:8950)
               y1=10053
               y2=10179
               x1=8901
               x2=8950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(8951:9000)
               y1=10180
               y2=10203
               x1=8951
               x2=9000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(9001:9050)
               y1=10203
               y2=10229
               x1=9001
               x2=9050
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9051:9100)
               y1=10230
               y2=10245
               x1=9051
               x2=9100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9101:9150)
               y1=10245
               y2=10262
               x1=9101
               x2=9150
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9151:9200)
               y1=10263
               y2=10343
               x1=9151
               x2=9200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9201:9250)
               y1=10344
               y2=10458
               x1=9201
               x2=9250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9251:9300)
               y1=10459
               y2=10525
               x1=9251
               x2=9300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9301:9350)
               y1=10527
               y2=10574
               x1=9301
               x2=9350
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9351:9400)
               y1=10575
               y2=10629
               x1=9351
               x2=9400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9401:9450)
               y1=10630
               y2=10700
               x1=9401
               x2=9450
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9451:9500)
               y1=10701
               y2=10741
               x1=9451
               x2=9500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(9501:9550)
               y1=10742
               y2=10919
               x1=9501
               x2=9550
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9551:9600)
               y1=10920
               y2=10692
               x1=9551
               x2=9600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9601:9650)
               y1=10973
               y2=11125
               x1=9601
               x2=9650
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9651:9700)
               y1=11126
               y2=11179
               x1=9651
               x2=9700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9701:9750)
               y1=11179
               y2=11203
               x1=9701
               x2=9750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9751:9800)
               y1=11204
               y2=11226
               x1=9751
               x2=9800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9801:9850)
               y1=11227
               y2=11245
               x1=9801
               x2=9850
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9851:9900)
               y1=11246
               y2=11267
               x1=9851
               x2=9900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9901:9950)
               y1=11267
               y2=11332
               x1=9901
               x2=9950
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(9951:10000)
               y1=11334
               y2=11482
               x1=9951
               x2=10000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(10001:10100)
               y1=11484
               y2=11685
               x1=10001
               x2=10100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(10101:10200)
               y1=11690
               y2=11907
               x1=10101
               x2=10200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(10201:10300)
               y1=11908
               y2=12163
               x1=10201
               x2=10300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(10301:10400)
               y1=12166
               y2=12337
               x1=10301
               x2=10400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(10401:10500)
               y1=12338
               y2=12457
               x1=10401
               x2=10500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(10501:10600)
               y1=12458
               y2=12632
               x1=10501
               x2=10600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(10601:10700)
               y1=12633
               y2=12699
               x1=10601
               x2=10700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(10701:10800)
               y1=12700
               y2=12765
               x1=10701
               x2=10800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(10801:10900)
               y1=12766
               y2=12847
               x1=10801
               x2=10900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(10901:11000)
               y1=12848
               y2=12910
               x1=10901
               x2=11000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(11001:11100)
               y1=12910
               y2=12985
               x1=11001
               x2=11100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(11101:11200)
               y1=12986
               y2=13107
               x1=11101
               x2=11200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(11201:11300)
               y1=13109
               y2=13195
               x1=11201
               x2=11300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(11301:11400)
               y1=13196
               y2=13291
               x1=11301
               x2=11400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(11401:11500)
               y1=13293
               y2=13407
               x1=11401
               x2=11500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(11501:11600)
               y1=13408
               y2=13478
               x1=11501
               x2=11600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(11601:11700)
               y1=13478
               y2=13588
               x1=11601
               x2=11700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(11701:11800)
               y1=13589
               y2=13697
               x1=11701
               x2=11800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(11801:11900)
               y1=13698
               y2=13791
               x1=11801
               x2=11900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(11901:12000)
               y1=13792
               y2=13971
               x1=11901
               x2=12000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(12001:12100)
               y1=13973
               y2=14110
               x1=12001
               x2=12100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(12101:12200)
               y1=14111
               y2=14232
               x1=12101
               x2=12200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(12201:12300)
               y1=14233
               y2=14417
               x1=12201
               x2=12300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(12301:12400)
               y1=14419
               y2=14618
               x1=12301
               x2=12400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(12401:12500)
               y1=14620
               y2=14821
               x1=12401
               x2=12500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(12501:12600)
               y1=14822
               y2=14958
               x1=12501
               x2=12600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(12601:12700)
               y1=14960
               y2=15100
               x1=12601
               x2=12700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(12701:12800)
               y1=15102
               y2=15294
               x1=12701
               x2=12800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(12801:12900)
               y1=15295
               y2=15587
               x1=12801
               x2=12900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(12901:13000)
               y1=15590
               y2=15854
               x1=12901
               x2=13000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(13001:13100)
               y1=15855
               y2=16006
               x1=13001
               x2=13100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(13101:13200)
               y1=16007
               y2=16128
               x1=13101
               x2=13200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(13201:13300)
               y1=16129
               y2=16233
               x1=13201
               x2=13300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(13301:13400)
               y1=16234
               y2=16355
               x1=13301
               x2=13400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(13401:13500)
               y1=16356
               y2=16472
               x1=13401
               x2=13500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(13501:13600)
               y1=16474
               y2=16690
               x1=13501
               x2=13600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(13601:13700)
               y1=16692
               y2=16855
               x1=13601
               x2=13700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(13701:13800)
               y1=16856
               y2=16990
               x1=13701
               x2=13800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(13801:13900)
               y1=16992
               y2=17145
               x1=13801
               x2=13900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(13901:14000)
               y1=17147
               y2=17245
               x1=13901
               x2=14000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(14001:14100)
               y1=17246
               y2=17330
               x1=14001
               x2=14100
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(14101:14200)
               y1=17331
               y2=17423
               x1=14101
               x2=14200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(14201:14300)
               y1=17424
               y2=17486
               x1=14201
               x2=14300
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(14301:14400)
               y1=17487
               y2=17549
               x1=14301
               x2=14400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(14401:14500)
               y1=17550
               y2=17642
               x1=14401
               x2=14500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(14501:14600)
               y1=17643
               y2=17891
               x1=14501
               x2=14600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(14601:14700)
               y1=17895
               y2=18035
               x1=14601
               x2=14700
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(14701:14800)
               y1=18036
               y2=18153
               x1=14701
               x2=14800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(14801:14900)
               y1=18154
               y2=18209
               x1=14801
               x2=14900
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(14901:15000)
               y1=18210
               y2=18252
               x1=14901
               x2=15000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(15001:15125)
               y1=18252
               y2=18301
               x1=15001
               x2=15125
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(15126:15250)
               y1=18302
               y2=18358
               x1=15126
               x2=15250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(15251:15375)
               y1=18359
               y2=18449
               x1=15251
               x2=15375
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(15376:15500)
               y1=18450
               y2=18729
               x1=15376
               x2=15500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(15501:15625)
               y1=18730
               y2=18895
               x1=15501
               x2=15625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(15626:15750)
               y1=18896
               y2=18988
               x1=15625
               x2=15750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(15751:15875)
               y1=18989
               y2=19064
               x1=15751
               x2=15875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(15876:16000)
               y1=19065
               y2=19144
               x1=15876
               x2=16000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(16001:16125)
               y1=19145
               y2=19236
               x1=16001
               x2=16125
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(16126:16250)
               y1=19237
               y2=19419
               x1=16126
               x2=16250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(16251:16375)
               y1=19420
               y2=19627
               x1=16251
               x2=16375
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(16376:16500)
               y1=19628
               y2=19815
               x1=16376
               x2=16500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(16501:16625)
               y1=19816
               y2=19904
               x1=16501
               x2=16625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(16626:16750)
               y1=19905
               y2=19980
               x1=16626
               x2=16750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(16751:16875)
               y1=19981
               y2=20109
               x1=16751
               x2=16875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(16876:17000)
               y1=20110
               y2=20272
               x1=16876
               x2=17000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(17001:17125)
               y1=20274
               y2=20524
               x1=17001
               x2=17125
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(17126:17250)
               y1=20526
               y2=20713
               x1=17126
               x2=17250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(17251:17375)
               y1=20714
               y2=20825
               x1=17251
               x2=17375
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(17376:17500)
               y1=20826
               y2=20922
               x1=17376
               x2=17500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(17501:17625)
               y1=20922
               y2=21028
               x1=17501
               x2=17625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(17626:17750)
               y1=21029
               y2=21194
               x1=17626
               x2=17750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(17751:17875)
               y1=21195
               y2=21452
               x1=17751
               x2=17875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(17876:18000)
               y1=21454
               y2=21713
               x1=17876
               x2=18000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(18001:18125)
               y1=21715
               y2=21868
               x1=18001
               x2=18125
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(18126:18250)
               y1=21869
               y2=21952
               x1=18126
               x2=18250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(18251:18375)
               y1=21953
               y2=22024
               x1=18251
               x2=18375
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(18376:18500)
               y1=22024
               y2=22114
               x1=18376
               x2=18500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(18501:18625)
               y1=22115
               y2=22250
               x1=18501
               x2=18625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(18626:18750)
               y1=22252
               y2=22467
               x1=18626
               x2=18750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(18751:18875)
               y1=22469
               y2=22710
               x1=18751
               x2=18875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(18876:19000)
               y1=22712
               y2=22866
               x1=18876
               x2=19000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19001:19125)
               y1=22866
               y2=22931
               x1=19001
               x2=19125
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19126:19250)
               y1=22932
               y2=23018
               x1=19126
               x2=19250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19251:19375)
               y1=23019
               y2=23133
               x1=19251
               x2=19375
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19376:19500)
               y1=23134
               y2=23282
               x1=19376
               x2=19500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19501:19625)
               y1=23283
               y2=23419
               x1=19501
               x2=19625
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19626:19750)
               y1=23420
               y2=23633
               x1=19626
               x2=19750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19751:19875)
               y1=23634
               y2=23800
               x1=19751
               x2=19875
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(19876:20000)
               y1=23801
               y2=23925
               x1=19876
               x2=20000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(20001:20200)
               y1=23926
               y2=24151
               x1=20001
               x2=20200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(20201:20400)
               y1=24152
               y2=24343
               x1=20201
               x2=20400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(20401:20600)
               y1=24344
               y2=24597
               x1=20401
               x2=20600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(20601:20800)
               y1=24598
               y2=24787
               x1=20601
               x2=20800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(20801:21000)
               y1=24788
               y2=25138
               x1=20801
               x2=21000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(21001:21200)
               y1=25139
               y2=25359
               x1=21001
               x2=21200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(21201:21400)
               y1=25359
               y2=25513
               x1=21201
               x2=21400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(21401:21600)
               y1=25514
               y2=25701
               x1=21401
               x2=21600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(21601:21800)
               y1=25703
               y2=26157
               x1=21601
               x2=21800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(21801:22000)
               y1=26159
               y2=26409
               x1=21801
               x2=22000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(22001:22200)
               y1=26410
               y2=26803
               x1=22001
               x2=22200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(22201:22400)
               y1=26804
               y2=27092
               x1=22201
               x2=22400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(22401:22600)
               y1=27095
               y2=27335
               x1=22401
               x2=22600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(22601:22800)
               y1=27336
               y2=27467
               x1=22601
               x2=22800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(22801:23000)
               y1=27468
               y2=27585
               x1=22801
               x2=23000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(23001:23200)
               y1=27585
               y2=27969
               x1=23001
               x2=23200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(23201:23400)
               y1=27970
               y2=28161
               x1=23201
               x2=23400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(23401:23600)
               y1=28162
               y2=28488
               x1=23401
               x2=23600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(23601:23800)
               y1=28490
               y2=28691
               x1=23601
               x2=23800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(23801:24000)
               y1=28692
               y2=28831
               x1=23801
               x2=24000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(24001:24200)
               y1=28832
               y2=28965
               x1=24001
               x2=24200
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(24201:24400)
               y1=28966
               y2=29166
               x1=24201
               x2=24400
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(24401:24600)
               y1=29168
               y2=29577
               x1=24401
               x2=24600
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(24601:24800)
               y1=29580
               y2=29848
               x1=24601
               x2=24800
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(24801:25000)
               y1=29849
               y2=29988
               x1=24801
               x2=25000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(25001:25250)
               y1=29989
               y2=30124
               x1=25001
               x2=25250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(25251:25500)
               y1=30125
               y2=30355
               x1=25251
               x2=25500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(25501:25750)
               y1=30358
               y2=30801
               x1=25501
               x2=25750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(25751:26000)
               y1=30802
               y2=30962
               x1=25751
               x2=26000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(26001:26250)
               y1=30962
               y2=31117
               x1=26001
               x2=26250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(26251:26500)
               y1=31118
               y2=31296
               x1=26251
               x2=26500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1     
           Case(26501:26750)
               y1=31297
               y2=31616
               x1=26501
               x2=26750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(26751:27000)
               y1=31616
               y2=31773
               x1=26751
               x2=27000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1                
           Case(27001:27250)
               y1=31773
               y2=31896
               x1=27001
               x2=27250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(27251:27500)
               y1=31896
               y2=32046
               x1=27251
               x2=27500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(27501:27750)
               y1=32046
               y2=32270
               x1=27501
               x2=27750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(27751:28000)
               y1=32271
               y2=32468
               x1=27751
               x2=28000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(28001:28250)
               y1=32469
               y2=32654
               x1=28001
               x2=28250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(28251:28500)
               y1=32654
               y2=32892
               x1=28251
               x2=28500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(28501:28750)
               y1=32893
               y2=33259
               x1=28501
               x2=28750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(28751:29000)
               y1=33260
               y2=33515
               x1=28751
               x2=29000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(29001:29250)
               y1=33516
               y2=33697
               x1=29001
               x2=29250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(29251:29500)
               y1=33697
               y2=33891
               x1=29251
               x2=29500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(29501:29750)
               y1=33892
               y2=34123
               x1=29501
               x2=29750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(29751:30000)
               y1=34124
               y2=34287
               x1=29751
               x2=30000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(30001:30250)
               y1=34288
               y2=34415
               x1=30001
               x2=30250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(30251:30500)
               y1=34416
               y2=34555
               x1=30251
               x2=30500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(30501:30750)
               y1=34555
               y2=34856
               x1=30501
               x2=30750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(30751:31000)
               y1=34857
               y2=35038
               x1=30751
               x2=31000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(31001:31250)
               y1=35039
               y2=35191
               x1=31001
               x2=31250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(31251:31500)
               y1=35192
               y2=35396
               x1=31251
               x2=31500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(31501:31750)
               y1=35397
               y2=35759
               x1=31501
               x2=31750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(31751:32000)
               y1=35760
               y2=36015
               x1=31751
               x2=32000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(32001:32250)
               y1=36016
               y2=36387
               x1=32001
               x2=32250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(32251:32500)
               y1=36390
               y2=37044
               x1=32251
               x2=32500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(32501:32750)
               y1=37045
               y2=37277
               x1=32501
               x2=32750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(32751:33000)
               y1=37278
               y2=37505
               x1=32751
               x2=33000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(33001:33250)
               y1=37506
               y2=37787
               x1=33001
               x2=33250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(33251:33500)
               y1=37788
               y2=39049
               x1=33251
               x2=33500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(33501:33750)
               y1=39054
               y2=39479
               x1=33501
               x2=33750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(33751:34000)
               y1=39480
               y2=39611
               x1=33751
               x2=34000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(34001:34250)
               y1=39612
               y2=39698
               x1=34001
               x2=34250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(34251:34500)
               y1=39698
               y2=39738
               x1=34251
               x2=34500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(34501:34750)
               y1=39737
               y2=39922
               x1=34501
               x2=34750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(34751:35000)
               y1=39923
               y2=40085
               x1=34751
               x2=35000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(35001:35250)
               y1=40085
               y2=40254
               x1=35001
               x2=35250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(35251:35500)
               y1=40255
               y2=40481
               x1=35251
               x2=35500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1    
           Case(35501:35750)
               y1=40483
               y2=41202
               x1=35501
               x2=35750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(35751:36000)
               y1=41203
               y2=41356
               x1=35751
               x2=36000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(36001:36250)
               y1=41356
               y2=41484
               x1=36001
               x2=36250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(36251:36500)
               y1=41485
               y2=41606
               x1=36251
               x2=36500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(36501:36750)
               y1=41606
               y2=41728
               x1=36501
               x2=36750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(36751:37000)
               y1=41728
               y2=41854
               x1=36751
               x2=37000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(37001:37250)
               y1=41855
               y2=41987
               x1=37001
               x2=37250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(37251:37500)
               y1=41988
               y2=42126
               x1=37251
               x2=37500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(37501:37750)
               y1=42126
               y2=42268
               x1=37501
               x2=37751
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(37751:38000)
               y1=42268
               y2=42425
               x1=37751
               x2=38000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(38001:38250)
               y1=42425
               y2=42612
               x1=38001
               x2=38250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(38251:38500)
               y1=42615
               y2=42895
               x1=38251
               x2=38500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1      
           Case(38501:38750)
               y1=42896
               y2=43111
               x1=38501
               x2=38750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(38751:39000)
               y1=43111
               y2=43256
               x1=38751
               x2=39000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(39001:39250)
               y1=43257
               y2=43380
               x1=39001
               x2=39250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(39251:39500)
               y1=43380
               y2=43496
               x1=39251
               x2=39500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(39501:39750)
               y1=43496
               y2=43612
               x1=39501
               x2=39750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(39751:40000)
               y1=43612
               y2=43738
               x1=39751
               x2=40000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1   
           Case(40001:40250)
               y1=43738
               y2=43887
               x1=40001
               x2=40250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(40251:40500)
               y1=43887
               y2=44072
               x1=40251
               x2=40500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(40501:40750)
               y1=44073
               y2=44289
               x1=40501
               x2=40750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(40751:41000)
               y1=44290
               y2=44553
               x1=40751
               x2=41000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(41001:41250)
               y1=44554
               y2=44804
               x1=41001
               x2=41250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(41251:41500)
               y1=44804
               y2=45005
               x1=41251
               x2=41500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(41501:41750)
               y1=45006
               y2=45189
               x1=41501
               x2=41750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(41751:42000)
               y1=45190
               y2=45381
               x1=41751
               x2=42000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(42001:42250)
               y1=45381
               y2=45586
               x1=42001
               x2=42250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(42251:42500)
               y1=45586
               y2=45830
               x1=42251
               x2=42500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(42501:42750)
               y1=45831
               y2=46129
               x1=42501
               x2=42750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(42751:43000)
               y1=46130
               y2=46632
               x1=42751
               x2=43000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1 
           Case(43001:43250)
               y1=46633
               y2=46844
               x1=43001
               x2=43250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(43251:43500)
               y1=46845
               y2=47045
               x1=43251
               x2=43500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(43501:43750)
               y1=47046
               y2=47242
               x1=43501
               x2=43750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(43751:44000)
               y1=47243
               y2=47438
               x1=43751
               x2=44000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(44001:44250)
               y1=47439
               y2=47637
               x1=44001
               x2=44250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(44251:44500)
               y1=47638
               y2=47834
               x1=44251
               x2=44500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(44501:44750)
               y1=47834
               y2=48035
               x1=44501
               x2=44750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(44751:45000)
               y1=48036
               y2=48246
               x1=44751
               x2=45000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1              
           Case(45001:45250)
               y1=48247
               y2=48470
               x1=45001
               x2=45500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(45251:45500)
               y1=48471
               y2=48710
               x1=45251
               x2=45500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(45501:45750)
               y1=48711
               y2=48969
               x1=45501
               x2=45750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(45751:46000)
               y1=48970
               y2=49248
               x1=45751
               x2=46000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(46001:46250)
               y1=49249
               y2=49546
               x1=46001
               x2=46250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(46251:46500)
               y1=49547
               y2=49863
               x1=46251
               x2=46500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(46501:46750)
               y1=49864
               y2=50209
               x1=46501
               x2=46750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(46751:47000)
               y1=50211
               y2=50576
               x1=46751
               x2=47000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(47001:47250)
               y1=50578
               y2=50939
               x1=47001
               x2=47250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(47251:47500)
               y1=50941
               y2=51282
               x1=47251
               x2=47500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(47501:47750)
               y1=51284
               y2=51592
               x1=47501
               x2=47750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(47751:48000)
               y1=51594
               y2=51886
               x1=47751
               x2=48000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(48001:48250)
               y1=51887
               y2=52175
               x1=48001
               x2=48250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(48251:48500)
               y1=52176
               y2=52461
               x1=48251
               x2=48500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(48501:48750)
               y1=52462
               y2=52745
               x1=48501
               x2=48750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(48751:49000)
               y1=52747
               y2=53035
               x1=48751
               x2=49000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1  
           Case(49001:49250)
               y1=53036
               y2=53340
               x1=49001
               x2=49250
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(49251:49500)
               y1=53341
               y2=53671
               x1=49251
               x2=49500
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(49501:49750)
               y1=53672
               y2=54028
               x1=49501
               x2=49750
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1
           Case(49751:50000)
               y1=54030
               y2=54426
               x1=49751
               x2=50000
               calendar_date=((real(i)-x1)*(y2-y1))/(x2-x1) + y1               
           Case(50001:)
               calendar_date=i                
           EndSelect    
        else
           calendar_date=c14date    
        endif
        
        return

    
    end function calendar_date





    
   
end module routines
