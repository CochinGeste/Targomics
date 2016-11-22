#!/usr/bin/env Rscript

#Sys.setlocale("LC_MESSAGES", 'en_GB.UTF-8')


ATCG_To_1234<-function(sequence_temp) ## OK 
{
  sequence_temp=gsub("A",as.numeric("1"),sequence_temp)
  sequence_temp=gsub("T",as.numeric("2"),sequence_temp) 
  sequence_temp=gsub("C",as.numeric("3"),sequence_temp)
  sequence_temp=gsub("G",as.numeric("4"),sequence_temp)
  sequence_temp = unlist(strsplit(sequence_temp,split = ""))
  return(sequence_temp)
}




creation_matrix_primer<- function(primer,Nombre_bases_demandees)
{
  taille_primer = length(primer)
  mini_primer1=primer[(taille_primer-Nombre_bases_demandees+1):taille_primer]
  
  mat_primer1=matrix(nrow=nombre_base_bornes -Nombre_bases_demandees +1 ,ncol=nombre_base_bornes ) 
  for (pos in 1: (nombre_base_bornes-Nombre_bases_demandees  +1)       )    
  {
    #pos = 1
    mat_primer1[pos,pos:(pos+length(mini_primer1)-1  )] = as.numeric(mini_primer1)
  }
  return(mat_primer1)
}








realignement<- function (Vecteur_position,taille_primer_temp,alignes_parfaits,Sequences_fastq_temp){
  # realignement(Vecteur_position,taille_primer_gauche,Position_alignes_parfait,Sequences_fastq_fwd) 
  #realignement(Vecteur_position,taille_primer_droit,Position_alignes_parfait,Sequences_fastq_fwd)
  #taille_primer_temp=taille_primer_droit
  #Sequences_fastq_temp=Sequences_fastq_revc
  #alignes_parfaits = Position_alignes_parfait
  position_ref = taille_primer_temp-ThreePrimeEnd_length+1
  alignes_parfaits=names(which(Vecteur_position==position_ref))
  Sequences_FullPrimer=Sequences_fastq_temp[alignes_parfaits,(taille_primer_temp+1):ncol(Sequences_fastq_temp)]
  
  Vecteur_position_temp_pas_alignes=Vecteur_position[which(Vecteur_position!=taille_primer_temp-ThreePrimeEnd_length+1)]
  
  Sequences_temp=matrix(nrow=length(Vecteur_position_temp_pas_alignes),ncol=ncol(Sequences_fastq_temp)-taille_primer_temp)
  rownames(Sequences_temp)=names(Vecteur_position_temp_pas_alignes)
  noms_pas_alignes=names(Vecteur_position_temp_pas_alignes)
  
  for (i in 1:length(Vecteur_position_temp_pas_alignes))
  {
    #i=1
    # i = 3493
    ## AS CHARACTER PRENDS PLUS DE TEMPS MAIS VRAIMENT PLUS 
    read_temp=as.numeric(Sequences_fastq_temp[as.numeric(noms_pas_alignes[i]),(as.numeric(Vecteur_position_temp_pas_alignes[i])+ThreePrimeEnd_length):ncol(Sequences_fastq_temp)])
    if (length(read_temp)>ncol(Sequences_fastq_temp)-taille_primer_temp)
    { 
      Sequences_temp[noms_pas_alignes[i],]=read_temp[1:(ncol(Sequences_fastq_temp)-taille_primer_temp)]
    } else {
      Sequences_temp[noms_pas_alignes[i],]=c(read_temp,rep(NA,ncol(Sequences_fastq_temp)-taille_primer_temp-length(read_temp)))
    }
  }
  colnames(Sequences_temp) = colnames(Sequences_FullPrimer)
  Sequences_P1Sens_Sens=rbind(Sequences_FullPrimer,Sequences_temp)
  
  
  return(Sequences_P1Sens_Sens)
}



###Matrix Substraction testing alignment of each read ("seq_read") on all possible positions of primer1 fragment ("mat_primer1"). Seuil is the tolerated number of mismatches (if more errors occur, NA is returned) 
alignement_seq_ref_contre_matrice<- function(seq_read,mat_primer1,seuil){
  sous= t(apply(mat_primer1,1,function(x) x=x-seq_read)) 
  sous[which(sous!=0)]=1  # Binary( 1 = miss match , 0 = match )
  SommeLignes=rowSums(1-sous,na.rm=T)
  vec= which(SommeLignes==max(SommeLignes))[1]
  if (SommeLignes[vec]>seuil)
  { 
    
    return(vec)
  } else {
    return(NA)
  }
  
}






complementaire <- function (seq_a_transformer){
  
  output=seq_a_transformer+1
  output[which(output==3)]=1
  output[which(output==5)]=3
  return(output)
}


ecriture_fasta  <- function (mat,file.out)
{
  #mat = Ilots_cpg_ACTG
  for (seq in 1:nrow(mat) )
  {
    seq_write=paste(mat[seq,],sep="")
    seq_write=gsub("NA","-",seq_write)
    write.fasta(seq_write, paste("Sequence_",seq,sep=""), nbchar = 60, file.out, open = "a")
  }
  
}

qc_output<-function(Sequence_temp) {  
  #Sequence_temp=Sequences_fastq
  noms_temp=rownames(Sequence_temp)
  vecteur_longueurs_temp=Vecteur_taille_total_reads_includedSmallReads[noms_temp]
  return(round(c(nrow(Sequence_temp),nrow(Sequence_temp)/length(Vecteur_taille_total_reads_includedSmallReads)*100,mean(vecteur_longueurs_temp),sd(vecteur_longueurs_temp))  ))
} 


compression_homopolymere <- function(matrice){
  #matrice = Sequences_fastq
  write.table(matrice,"temp2_toto",sep="",col.names = FALSE,row.names = FALSE)
  ligne_commande  = paste(path_script,"/compress160129.sh temp2_toto",sep="")
  #ligne_commande = "/home/mario/Documents/Mario/Projets/Bisufilte_seq/FastQ/compress.sh temp2"
  system(ligne_commande)
  return(lecture_homopolymere())
}

lecture_homopolymere <- function(){
  
  Sequences_fastq2 <- read.csv("temp2_toto",sep="\n",stringsAsFactors = TRUE ,colClasses = "vector" ,header= F) # Lecture du fichier_reads 
  #nchar(Sequences_fastq2)
  tailles=sort(apply(Sequences_fastq2,1,function(x) nchar(x))/2)
  
  #tailles=    apply(Sequences_fastq2  ,1,function(x) nchar(x))/2 
  ncolonnes=tailles[round(0.975*length(tailles))]
  chaine=paste(paste(rep(NA,ncolonnes),sep='',collapse=','),",",sep='')
  write.table(rbind(chaine,chaine,chaine,chaine,chaine,Sequences_fastq2),"temp2.csv",row.names = FALSE, col.names = FALSE,quote=FALSE )
  #a1<- read.table(input,colClasses = "character",header=FALSE )
  #Sequences_fastq <- read.fwf(input,widths=rep(1, max(nchar(a1$V1))) ,colClasses = "numeric",header=FALSE) # Lecture du fichier_reads 
  mat_return =  read.csv("temp2.csv",sep=",",header=FALSE,fill = T)
  mat_return = mat_return[6:dim(mat_return)[1],]
  system("rm temp2_toto")
  system("rm temp2.csv")
  #mat_return [,which(apply(mat_return,2,function(x) all(is.na(x))))] = NULL
  return(mat_return) 
}


wrapper_alignement = function(path_absolu,nom_fichier_tableau,nom_fichier_fastq,dossier_output,nom_tissu){
  dir.create(dossier_output,recursive = T, showWarnings = FALSE)
  path_script <<- paste(path_absolu,"Scripts/",sep="") 
  input_fastq = paste(path_absolu,"Input/",nom_fichier_fastq,sep="")
  input_tableau = paste(path_absolu,"Input/",nom_fichier_tableau,sep="")
  tableau_input= read.csv(input_tableau ,header=T)
  tableau_input$Sequence=sapply(tableau_input$Sequence,function(x) gsub("[\r\n]", "", x) ) ## Removing the "return" char
  ##Function
  
  ### Part 1 : Loading and reformating the input files
  ##Fastq sequences are transformed into a matrix of numbers, following the ATCG_to_1234 conversion
  
  ## launching a bash script "ATCG_to_1234.sh", transforming sequences (ex: "ATCTTTTG") within the fastq, into coma separated numbers (ex: "1,2,3,2,2,2,2,4,") in a text file
  
  bash_script=paste(path_script,"ATCG_to_1234.sh",sep="")
  system(paste(bash_script,input_fastq))
  fichier_csv= paste(input_fastq,"_numeric.tmp",sep="")
  Sequences_fastq2 <- read.table(fichier_csv,sep="\n",stringsAsFactors = TRUE ,colClasses = "vector" ) # Reading previously created # 
  
  
  ### Transforming the sequences of the fastq transformed file into a matrix of bases (bases are coded in 1,2,3 and 4) ##
  #finding the longest read (will determine the width of the array, after exclusion of the 5% extreme sizes)
  tailles=sort(apply(Sequences_fastq2,1,function(x) nchar(x))/2) 
  ncolonnes=tailles[round(0.975*length(tailles))] # we consider the size of the 97.5% longest read in case longer reads are odds
  chaine=paste(paste(rep(NA,ncolonnes),sep='',collapse=','),",",sep='') ## We create a temporary file "temp.csv", by adding, to the fastq transformed file, 5 lines with "ncolonnes" repeats of "NA,". Subsequently, R ("read.csv" function) will interprete the "temp.csv" file a table of "ncolonnes" columns.
  write.table(rbind(chaine,chaine,chaine,chaine,chaine,Sequences_fastq2),"temp.csv",row.names = FALSE, col.names = FALSE,quote=FALSE ) 
  
  Sequences_fastq =  read.csv("temp.csv",sep=",",header=FALSE,colClasses = "numeric" )  ## creating the matrix of bases
  Sequences_fastq = Sequences_fastq[6:dim(Sequences_fastq)[1],] 
  ## Removing the "NA," repeats (5 first lines)
  Sequences_fastq=as.matrix(Sequences_fastq) 
  
  Sequences_fastq_comp = compression_homopolymere(Sequences_fastq) ### Compression of the homopolymeres ### 
  
  Sequences_fastq = Sequences_fastq_comp 
  
  rownames(Sequences_fastq_comp) = seq(1,nrow(Sequences_fastq_comp))
  colnames(Sequences_fastq_comp) = paste("B",seq(1,ncol(Sequences_fastq_comp)),sep='')
  
  Nombre_total_reads=nrow(Sequences_fastq) # Total Number of reads
  
  Vecteur_taille_total_reads_includedSmallReads <<- ncol(Sequences_fastq) - ( apply(Sequences_fastq_comp ,1, function (x) length(which(is.na(x))) ) )+1 #vector of reads sizes
  
  Sequences_fastq_comp = Sequences_fastq_comp[,  1:max(Vecteur_taille_total_reads_includedSmallReads) ] ## Recadrage 
  
  # Removal of littles reads #

  Vecteur_taille_total_reads = Vecteur_taille_total_reads_includedSmallReads
  
  Sequences_fastq = Sequences_fastq_comp 
  
  #creating the reverse complementary strand of all reads ("Sequences_fastq_revc")
  Sequences_fastq_fwd <<- Sequences_fastq_comp 
  
  Sequences_fastq_fwd_save = Sequences_fastq_fwd
  
  Sequences_fastq_revc<<-t(apply(Sequences_fastq_fwd,1,function(x) x = c(rev(na.omit(x)),rep(NA,length(which(is.na(x))))))) 
  Sequences_fastq_revc<<-t(apply(Sequences_fastq_revc,1,function(x) x=complementaire(x)))
  
  Sequences_fastq_revc_save = Sequences_fastq_revc
  
  rownames(Sequences_fastq_revc)=rownames(Sequences_fastq_fwd)
  colnames(Sequences_fastq_revc)=colnames(Sequences_fastq_fwd)
  
  
  #initializing the global QC output matrix
  QC_total=matrix(nrow=0,ncol=9)
  colnames(QC_total)=c("Total_reads","P1Sens_sens","P1Antisens_sens","P2Sens_antisens","P2Antisens_antisens","Total_AlignedOnPrimers","Total_P1","Total_P2","Total_Aligned_OnRefSeq")
  
  
  file.remove(fichier_csv)
  
  ################
  ################
  ################
  ###############
  #################
  ################
  vecteur_deja_alignes = NULL
  
  ## For each region in "Reference_sequences_example1.csv" (main loop), alignment from "Example.fastq" is attempted ##
  for (u in 1:nrow(tableau_input)) 
  {
    #u=1
    print (tableau_input[u,1])
    ptm <- proc.time()
    nom_ech = as.character(tableau_input[u,1])
    
    seq_ref = tableau_input[u,2] #reading the reference sequence
    seq_ref=toupper(seq_ref)
    # Creating the methylated reference sequence after bisulfite: the CG remain CG; the other C become T
    seq_ref_temp1 = gsub("CG","X",  seq_ref )  
    seq_ref_temp2 = gsub("C","T",  seq_ref_temp1 )
    seq_ref_methyle = gsub("X","CG",  seq_ref_temp2 ) 
    taille_primer_gauche=tableau_input[u,3] #primers size
    taille_primer_droit=tableau_input[u,4]
    
    # Creating the unmethylated reference sequence after bisulfite: all the C become T
    seq_ref_unmethyle =  gsub("C","T",  seq_ref )
    taille_amplic <<- nchar(seq_ref_methyle) ##  Ref sequence size (including primers) 
    ## Transforming the bases in the reference sequences into numbers (A,T,C,G to 1,2,3,4)
    
    
    
    
    seq_ref_methyle_1234=as.numeric(ATCG_To_1234(seq_ref_methyle))
    seq_ref_unmethyle_1234=as.numeric(ATCG_To_1234(seq_ref_unmethyle)) 
    
    ### Compress of the homopolymeres ## AAAAAAA -> A etc
    seq_ref_methyle_1234=rbind(seq_ref_methyle_1234,rep(NA,length(seq_ref_methyle_1234))) ## Passage en matrix for using the compression_homopolymere
    seq_ref_methyle_1234_compress = compression_homopolymere(seq_ref_methyle_1234)[1,] ## Getting the 1st line 
    taille_amplic_compress <<- min(which(is.na(seq_ref_methyle_1234_compress))) - 1
    seq_ref_methyle_1234_compress = as.numeric(seq_ref_methyle_1234_compress [1:taille_amplic_compress] ) # Recadrage
    seq_ref_methyle_1234 = seq_ref_methyle_1234_compress
    rm(seq_ref_methyle_1234_compress)
    taille_amplic_initial = taille_amplic
    taille_amplic_total = taille_amplic
    taille_amplic <<- taille_amplic_compress
    
    seq_ref_unmethyle_1234=rbind( seq_ref_unmethyle_1234,rep(NA,length( seq_ref_unmethyle_1234))) ## Passage en matrix for using the compression_homopolymere
    seq_ref_unmethyle_1234_compress = compression_homopolymere( seq_ref_unmethyle_1234)[1,] ## Getting the 1st line 
    seq_ref_unmethyle_1234_compress = as.numeric(seq_ref_unmethyle_1234_compress[1:taille_amplic_compress])
    seq_ref_unmethyle_1234 =  seq_ref_unmethyle_1234_compress
    rm(seq_ref_unmethyle_1234_compress)
    
    ## Creation of compressed primers  ##
    
    primer1 = as.numeric( ATCG_To_1234(seq_ref_methyle)  [1:taille_primer_gauche] )
    primer1 = rbind(primer1,rep(NA,length(primer1)))
    primer1_compress = compression_homopolymere(primer1)[1,] 
    taille_primer_gauche = min(which(is.na(primer1_compress))) -1 
    primer1_compress  =  primer1_compress [ 1 : taille_primer_gauche]
    primer1 = as.numeric(primer1_compress)
    rm(primer1_compress)
    
    primer2 = as.numeric( ATCG_To_1234(seq_ref_methyle) [ ((taille_amplic_initial - taille_primer_droit)+1) :taille_amplic_initial ] )
    primer2= rbind(primer2,rep(NA,length(primer2)))
    primer2_compress = compression_homopolymere(primer2)[1,]
    taille_primer_droit = min(which(is.na(primer2_compress)))  -1 
    primer2_compress = primer2_compress[1: taille_primer_droit ]
    primer2 = as.numeric( primer2_compress)
    rm(primer2_compress)
    
    
    ## Removal of primers ## 
    
    seq_ref_coupee_meth=seq_ref_methyle_1234[(taille_primer_gauche+1) : (taille_amplic -taille_primer_droit) ]
    seq_ref_coupee_meth_rev_comp  = rev ( complementaire(seq_ref_methyle_1234)) 
    seq_ref_coupee_meth_rev_comp  = seq_ref_coupee_meth_rev_comp[(taille_primer_droit+1) : (taille_amplic -taille_primer_gauche) ]
    
    
    ### Genotypage of CG ###
    
    position_des_C=which(seq_ref_coupee_meth==3)
    position_des_C_plus1=which(seq_ref_coupee_meth[position_des_C+1]==4)
    pos_sens =position_des_C[position_des_C_plus1]
    pos=pos_sens
    pos_reverse =c(pos_sens+1) 
    pos_total=c(pos_sens,pos_reverse+1)    
    
    
    
    
    nombre_base_bornes <<- 40 ## Fraction of the length of the read used 
    
    ### Part 2 : Alignements of reads on primers 
    ####### Sens strand: 5' Primer1 XXXXXXXXXX Primer2_c  3'
    ####### Antisens strand: 5' Primer2 XXXXXXXX Primer1_c 3'
    
    
    ###Creating 2 matrixes (one for each primer), with "ThreePrimeEnd_length" bases of each primer.
    ##This matrix corresponds to all possible positions of the "ThreePrimeEnd_length" bases within "nombre_base_bornes" bases 
    ThreePrimeEnd_length <<- 5 #number of bases from primer1 to be used for the alignment (starting from 3' end) #Should not be higher than primer length
    
    ThreePrimeEnd_length = min(ThreePrimeEnd_length, taille_primer_droit,taille_primer_gauche)
    
    MinimalProportionOfMatches<<-0.9 # error rate for alignment on primer
    
    mat_primer1 = creation_matrix_primer(primer1,ThreePrimeEnd_length) ## creating a matrix with the 3' end of the primer of length "nomb
    colnames(mat_primer1)=colnames(Sequences_fastq_fwd)[1:nombre_base_bornes]
    primer2 =   rev ( complementaire(seq_ref_methyle_1234)) ##primer2 is the rev compl of the end of the reference seq
    primer2 = primer2[1:taille_primer_droit]
    mat_primer2 = creation_matrix_primer(primer2,ThreePrimeEnd_length)
    colnames(mat_primer2)=colnames(Sequences_fastq_revc)[1:nombre_base_bornes]
    
    
    
    ### Alignement on Primer1 of Sens strands (fwd reads) ; output: Sens reads ("Sequences_P1Sens_Sens")
    
    Vecteur_position  = apply(Sequences_fastq_fwd[,1:nombre_base_bornes],1,function(x)  alignement_seq_ref_contre_matrice(x,mat_primer1,ThreePrimeEnd_length*MinimalProportionOfMatches) ) # A quelle position ca s'aligne sur le read 
    Position_alignes_parfait = Vecteur_position[which(!is.na(Vecteur_position))] # Ceux qui sont alignés au début ou au milieu 
    if (length(Position_alignes_parfait)>0) {
      
      Position_Notaligned = names(Vecteur_position[which(is.na(Vecteur_position))])
      Sequences_P1Sens_sens= realignement(Vecteur_position,taille_primer_gauche,Position_alignes_parfait,Sequences_fastq_fwd)
      
      vecteur_deja_alignes = c(vecteur_deja_alignes,Position_alignes_parfait)
      #Sequences_fastq_fwd = Sequences_fastq_fwd[Position_Notaligned,]  #removing the reads successfully aligned
      #Sequences_fastq_revc = Sequences_fastq_revc[Position_Notaligned,] #removing the reads successfully aligned
    } else {
      Sequences_P1Sens_sens = matrix(nrow=0,ncol = ncol(Sequences_fastq_fwd) -taille_primer_gauche  )
    }
    

    ### Alignement on Primer2 of antisens strands (fwd reads); output: Antisens reads ("Sequences_P2Sens_Antisens")
    Vecteur_position = apply(Sequences_fastq_fwd[,1:nombre_base_bornes],1,function(x)  alignement_seq_ref_contre_matrice(x,mat_primer2,ThreePrimeEnd_length*MinimalProportionOfMatches) )
    
    Position_alignes_parfait =  Vecteur_position[which(!is.na(Vecteur_position))]

    
    if (length(Position_alignes_parfait)>0) {
      Position_Notaligned = names(Vecteur_position[which(is.na(Vecteur_position))])
      Sequences_P2Sens_antisens= realignement(Vecteur_position,taille_primer_droit,Position_alignes_parfait,Sequences_fastq_fwd)
      vecteur_deja_alignes = c(vecteur_deja_alignes,Position_alignes_parfait)
    }else {
      Sequences_P2Sens_antisens = matrix(nrow=0,ncol = ncol(Sequences_fastq_fwd) -taille_primer_droit  )
    }
    
 
    ### Alignement on Primer 2 of Sens strands (revc reads); output: Antisens reads ("Sequences_P2AntiSens_Antisens")
    Vecteur_position  = apply(Sequences_fastq_revc[,1:nombre_base_bornes],1,function(x)  alignement_seq_ref_contre_matrice(x,mat_primer2,ThreePrimeEnd_length*MinimalProportionOfMatches) )
    Position_alignes_parfait= Vecteur_position[which(!is.na(Vecteur_position))]

    #Position_alignes_parfait =Position_alignes_parfait[setdiff(names(Position_alignes_parfait),names(vecteur_deja_alignes) )]
    
    if (length(Position_alignes_parfait)>0) {
      Position_Notaligned = names(Vecteur_position[which(is.na(Vecteur_position))])
      Sequences_P2AntiSens_antisens= realignement(Vecteur_position,taille_primer_droit,Position_alignes_parfait,Sequences_fastq_revc)
      vecteur_deja_alignes = c(vecteur_deja_alignes,Position_alignes_parfait)
    }else {
      Sequences_P2AntiSens_antisens = matrix(nrow=0,ncol = ncol(Sequences_fastq_revc) -taille_primer_droit  )
    }
 
    
    ### Alignement on Primer1 of Antisens strands (revc reads; output: Sens reads ("Sequences_P1Antisens_Sens")
    Vecteur_position  = apply(Sequences_fastq_revc[,1:nombre_base_bornes],1,function(x)  alignement_seq_ref_contre_matrice(x,mat_primer1,ThreePrimeEnd_length*MinimalProportionOfMatches) )
    Position_alignes_parfait= Vecteur_position[which(!is.na(Vecteur_position))]

    #Position_alignes_parfait =Position_alignes_parfait[setdiff(names(Position_alignes_parfait),names(vecteur_deja_alignes) )]
    if (length(Position_alignes_parfait)>0) {
      Position_Notaligned = names(Vecteur_position[which(is.na(Vecteur_position))])
      Sequences_P1Antisens_sens= realignement(Vecteur_position,taille_primer_gauche,Position_alignes_parfait,Sequences_fastq_revc)
      vecteur_deja_alignes = c(vecteur_deja_alignes,Position_alignes_parfait)
    }else {
      Sequences_P1Antisens_sens = matrix(nrow=0,ncol = ncol(Sequences_fastq_revc) -taille_primer_gauche  )
    }
    
    
  

    
    #Merging all the reads aligned on P1 (all the reads are in Sens strand)
    Merge_sens = rbind(Sequences_P1Sens_sens,Sequences_P1Antisens_sens)
    #Merging all the reads aligned on P2 (all the reads are in AntiSens strand)
    Merge_antisens = rbind(Sequences_P2Sens_antisens,Sequences_P2AntiSens_antisens)

    ## Part 3 : Alignements of reads (successfully aligned on primers) on the reference sequence between the primers
    
    
    
    if (nrow(Merge_sens) > 2) {
      sous=t(apply(  Merge_sens[,1:length(seq_ref_coupee_meth)][,-pos]  ,1,function(x) x=x-seq_ref_coupee_meth[-pos] )) 

      sous[which(sous!=0)]=1  # Binary( 1 = miss match , 0 = match )
      SommeLignes=rowSums(1-sous,na.rm=T)
      n_na=apply(Merge_sens[,1:length(seq_ref_coupee_meth)],1,function(x) x=length(which(is.na(x))))
      t_ref=length(seq_ref_coupee_meth)*0.4
      indice_reads_trop_petit = names(which(n_na > t_ref))
      seuils=0.9*(length(seq_ref_coupee_meth[-pos_total]) -n_na)
      
      Position_aligned = names(which( SommeLignes-seuils > 0 ))
      Position_aligned = setdiff(Position_aligned,indice_reads_trop_petit) ## On retire les trop petits
      Position_Notaligned = names(which( SommeLignes-seuils <= 0 ))
      Merge_sens_alignedRef = Merge_sens[Position_aligned,1:length(seq_ref_coupee_meth)]
      Merge_sens=Merge_sens[Position_Notaligned, ]
     
    } else {
      Merge_sens_alignedRef= matrix(NA,nrow=0,ncol = length(seq_ref_coupee_meth ) ) 
    }
    
    
    
    
    if ( nrow(Merge_sens_alignedRef)==0 ) {
      Merge_sens_alignedRef= matrix(NA,nrow=0,ncol = length(seq_ref_coupee_meth ) )
      sens = Merge_sens_alignedRef
    } else {
      
      sens = Merge_sens_alignedRef[,1:ncol(Merge_sens_alignedRef)]
    }
    
    
    
    if (nrow(Merge_antisens) > 2 ) {
      sous=t(apply(  Merge_antisens[,1:length(seq_ref_coupee_meth_rev_comp)] [,-pos] ,1,function(x) x=x-seq_ref_coupee_meth_rev_comp[-pos]))
      sous[which(sous!=0)]=1  # Binary( 1 = miss match , 0 = match )
      SommeLignes=rowSums(1-sous,na.rm=T)
      n_na=apply(Merge_antisens[,1:length(seq_ref_coupee_meth_rev_comp)],1,function(x) x=length(which(is.na(x))))
      t_ref=length(seq_ref_coupee_meth)*0.4
      indice_reads_trop_petit = names(which(n_na > t_ref))
      seuils=0.9*(length(seq_ref_coupee_meth[-pos_total]) -n_na)
      Position_aligned = names(which( SommeLignes-seuils > 0 ))
      Position_aligned = setdiff(Position_aligned,indice_reads_trop_petit) ## On retire les trop petits
      Position_Notaligned = names(which( SommeLignes-seuils <= 0 ))
      Merge_antisens_alignedRef = Merge_antisens[Position_aligned,1:length(seq_ref_coupee_meth_rev_comp)]
      Merge_antisens=Merge_antisens[Position_Notaligned, ]
    } else {
      Merge_antisens_alignedRef = matrix(NA,nrow=0,ncol = length(seq_ref_coupee_meth_rev_comp))
    }
    
    
    if ( nrow(Merge_antisens_alignedRef)<2 ) {
      Merge_antisens_alignedRef = matrix(NA,nrow=0,ncol = length(seq_ref_coupee_meth_rev_comp))
      Merge_antisens_final_revc = Merge_antisens_alignedRef
      antisens = Merge_antisens_final_revc

    }else {
      
      Merge_antisens_final_revc=t(apply(Merge_antisens_alignedRef,1,function(x) x = c(rep(NA,length(which(is.na(x)))),rev(na.omit(x))))) 
      Merge_antisens_final_revc=t(apply(Merge_antisens_final_revc,1,function(x) x=complementaire(x)))
      antisens = cbind(    Merge_antisens_final_revc[,  1:( ncol(Merge_antisens_final_revc) ) ]   )
    } 
    
    
    
    colnames(antisens) = colnames(sens)
    
    Merge_alignedRef=rbind(sens,antisens)
    
    
    if ( nrow(Merge_alignedRef) == 0 )
    {
      Merge_alignedRef = matrix(NA,nrow = 2,ncol = length(seq_ref_coupee_meth_rev_comp)   )
    }
    
    ############################""
    ## Part 4 : Genotyping the CpG and tabular outputs
    
    ##creating a Fasta file with the aligned sequences ("Res/Alignement_IslandName.fasta")
    #converting back the sequence code from 1,2,3,4 to A,T,C,G
    
 
    nom_file_out=paste(dossier_output,"/Alignement_",  nom_ech    ,".fasta",sep="")
    write.table(Merge_alignedRef,"tempo_ACTG.tmp",row.names = F,col.names = F)
    bash_script2=paste(path_script,"1234_to_ATCG.sh",sep="")
    system(paste(bash_script2,"tempo_ACTG.tmp"))
    Merge_alignedRef_ATCG  <- read.table("tempo_ACTG.tmp",sep="\n" ,colClasses = "vector") # Reading previously created # 
    Merge_alignedRef_ATCG = t( apply(Merge_alignedRef_ATCG,1, function(x) x= unlist(       strsplit(x,split=',') )  ))
    colnames(Merge_alignedRef_ATCG) = colnames(Merge_alignedRef)
    rownames(Merge_alignedRef_ATCG) =rownames(Merge_alignedRef)
    system("rm tempo_ACTG.tmp")
    if (file.exists(nom_file_out)) file.remove(nom_file_out)
    
    for (seq in 1:nrow(Merge_alignedRef_ATCG) )
    {
      seq_write=paste(Merge_alignedRef_ATCG[seq,],sep="")
      seq_write=gsub("NA","-",seq_write)
      write.fasta(seq_write, paste("Sequence_",seq,sep=""), nbchar = 60, nom_file_out, open = "a")
    }

    if (length(pos) == 1) {
      nC_sens = which(Merge_alignedRef[,pos_sens]==3)
      nT_sens = which(Merge_alignedRef[,pos_sens]==2)
      nG_sens = which(Merge_alignedRef[,pos_sens]==4)
      nA_sens = which(Merge_alignedRef[,pos_sens]==1)
      
      nC_reverse = which(Merge_alignedRef[,pos_reverse]==3)
      nT_reverse = which(Merge_alignedRef[,pos_reverse]==2)
      nG_reverse = which(Merge_alignedRef[,pos_reverse]==4)
      nA_reverse = which(Merge_alignedRef[,pos_reverse]==1)
      
      
      if (length(nC)==0) nC =0
      if (length(nT)==0) nT =0
      if (length(nG)==0) nG =0
      if (length(nA)==0) nA =0
    } else {
      
      nC_sens = apply(Merge_alignedRef[,pos_sens],2, function(x) length(which(x==3)) )
      nT_sens = apply(Merge_alignedRef[,pos_sens],2, function(x) length(which(x==2)) )
      nG_sens = apply(Merge_alignedRef[,pos_sens],2, function(x) length(which(x==4)) )
      nA_sens = apply(Merge_alignedRef[,pos_sens],2, function(x) length(which(x==1)) )
      
      nC_reverse = apply(Merge_alignedRef[,pos_reverse],2, function(x) length(which(x==3)) )
      nT_reverse = apply(Merge_alignedRef[,pos_reverse],2, function(x) length(which(x==2)) )
      nG_reverse = apply(Merge_alignedRef[,pos_reverse],2, function(x) length(which(x==4)) )
      nA_reverse = apply(Merge_alignedRef[,pos_reverse],2, function(x) length(which(x==1)) )
      
      
      
    }
    

    
    
    tab_out_sens = cbind(pos_sens,nC_sens,nT_sens)
    tab_out_rev = cbind(pos_reverse,nG_reverse,nA_reverse)
    
    tab_out_final = rbind(tab_out_sens,tab_out_rev)
    tab_out_final2= cbind(   rep(nom_ech,length(pos_total))  ,    tab_out_final    )
    colnames(tab_out_final2)= c("Ilot_CpG","Index","N_Methyl","N_Unmethyl")
    
    tab_out_final2 = tab_out_final2[order(    as.numeric   (tab_out_final2[,"Index"])      ),]
    
    tab_out_final2[,"Index"]=as.numeric(tab_out_final2[,"Index"]) + taille_primer_gauche
    
    nom_tab_out_final2=paste(dossier_output,"/Coverage_CpG",nom_ech,"_",nom_tissu   ,".txt",sep="")
    
    
    
    
    
    write.table(tab_out_final2,file =nom_tab_out_final2,row.names=FALSE, col.names=FALSE, sep=" ")
    
    
    
    QC=matrix(nrow=6,ncol=9)
    
    ## Merging for stats ##
    
    noms_temp=c( rownames(Sequences_P2Sens_antisens) , rownames(Sequences_P2AntiSens_antisens) , rownames(Sequences_P1Sens_sens) , rownames(Sequences_P2AntiSens_antisens) )
    vecteur_longueurs_temp=Vecteur_taille_total_reads_includedSmallReads[noms_temp]
    res =(round(c(length(noms_temp),length(noms_temp)/length(Vecteur_taille_total_reads_includedSmallReads)*100,mean(as.numeric(vecteur_longueurs_temp)),sd(as.numeric(vecteur_longueurs_temp)))  ))
    
    QC=cbind(qc_output(Sequences_fastq),qc_output(Sequences_P1Sens_sens),qc_output(Sequences_P1Antisens_sens),qc_output(Sequences_P2Sens_antisens),qc_output(Sequences_P2AntiSens_antisens),res,qc_output(Merge_sens_alignedRef),qc_output(Merge_antisens_alignedRef),qc_output(Merge_alignedRef))
    QC=rbind(rep(NA,9),c( paste(nom_ech,nom_tissu,sep="") ,NA,NA,NA,NA,NA,NA,NA,NA),QC)
    colnames(QC)=c("Total_reads","P1Sens_sens","P1Antisens_sens","P2Sens_antisens","P2Antisens_antisens","Total_AlignedOnPrimers","Total_P1","Total_P2","Total_Aligned_OnRefSeq")
    rownames(QC)=c("","Island","N","Percent_ofTotalReads","Mean_ReadLength","Sd_ReadLength")  
    QC_total=rbind(QC_total,QC)
    print(proc.time() - ptm)
    
  }
  
  QC_total=as.matrix(QC_total)
  
  
  nom_fichier_QC=paste(dossier_output,"/zQC_project_comp.csv",sep="")
  
  write.csv(QC_total,paste(dossier_output,"/zQC_project_comp.csv",sep=""),na="")
  system("rm temp.csv")
  
}


##############################################################################################
##Main



library(seqinr) ## Loading librairy 


# ## Usage ##
## Script developped by Mario NEOU & Guillaume Assie ##
## Version 1.0 (January 2016) ###
## This script MUST BE RUN IN LINUX ## (Bash functions are included)


# Aim : This script allows the alignment of reads (Fastq file), generated by NGS, on short reference sequences (a few hundred reads) #
# It was specificly developped for aligning Bisulfite-treated DNA of CpG islands. 
# These regions are hard to align, with >50% of GC and many repeats, with almost no C (3-bases genetic code). Standard aligment programs have not been developped to align properly these types of sequences 
# This program then counts the methylated and the unmethylated alleles for each CPG 

## Inputs : For each sample:
## 1. A Fastq file contaning the reads to be aligned ("Reads_example1.fastq")
## 2. A  csv table with 4 columns containing ("Refence_sequences_example1.csv"):
##    - column 1 (string): the name of the regions
##    - column 2 (string): the reference sequence for alignement. This sequence must include the 2 primers . Methylated and non methylated versions of the sequence after bisulfite transformation will be automatically generated)
##    - column 3 (number): size of the 5' primer
##    - column 4 (number): Size of the 3' primer 


## Output:
## 1. One table in csv format reporting for each region : 
##    - in columns: each CpG
##    - in lines: the number of methylated Cytosin, the number of unmethylated cytosin, the total  number of cytosins, the % of methylated cytosins
## 2. One Fasta file for each region, containing all the reads that could be aligned, in positive strand (the 2 primers sequences are removed)
## 3. One QC table with:
##    - Alignment efficiencies for each region
##    - Size of aligned reads

## Dependance : Seqinr : https://cran.r-project.org/web/packages/seqinr/seqinr.pdf ## A retirer ## 

##Working directory structure
#path/     #contains scripts and temp files
#path/res  #contains the output files
#path/input #contains the input files


## Plan :

## Parameters ##to be specified by the user (path and input file names)
## Functions
## Main
## Part 1 : Loading and reformating the input files
## Part 2 : Alignements of reads on primers 
## Part 3 : Alignements of reads (successfully aligned on primers) on the reference sequence between the primers
## Part 4 : Genotyping the CpG and tabular outputs

#################################################

##parameters

working_directory="/home/simong/Bureau/MS_Version161018/Targomics_methyl_temp/" # Set your working directory, don't forget the "/" at the end
## working_directory="/home/AssieG/Bureau/MS_Version/Targomics/"



path_absolu <<- working_directory
nom_fichier_tableau ="Reference_sequences_example1.csv"

nom_fichier_fastq ="Exemple.fastq"
dossier_output = paste(path_absolu,"Exemple_output",sep="")
nom_tissu = "Exemple"
wrapper_alignement(path_absolu,nom_fichier_tableau,nom_fichier_fastq,dossier_output,nom_tissu)


















