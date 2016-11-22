
#Different steps:
#1. Loading / Normalisation / filtering NGS data (read counts for all amplicons, SNP read counts
#2. Calling homozygous deletions and high level amplification (by gene and by gene subregions)
#3. Calling gains and losses of one DNA copy by combining allelic ratios and read counts
#4. Plotting outputs

#Normalisation / filtering raw NGS data
#1.1 Opening coverage data
data=read.csv("Example_NReads_TotalBySample_161013.csv",sep=';') #beware of the sep=';': mandatory in France, but default should be ',' in all other countries!!
nbre_col=ncol(data)
#str(data) #data[1:2,]
rownames(data)=as.character(data[,1])
data_couv=data[,7:ncol(data)] #specify the columns containing the numbers of amplicons only and discard the columns with annotations
#str(data_couv)
#data_couv[1:2,]
echantillons=colnames(data_couv)


#opening amplicons annotations
annotations_Amplicons=read.csv("Example_IAD37558_384WellPlateDataSheet_130920.csv",sep=';')#beware of the sep=';'
#str(annotations_Amplicons)
rownames(annotations_Amplicons)=as.character(annotations_Amplicons$Amplicon_ID)
annotations_Amplicons=annotations_Amplicons[,c(1,8,10,11,14)] #please keep only the Amplicon name (rowname), the Library pool (1 or 2), the gene symbol, the Chr (in numeric value: X->23; Y->24), the amplicon start and end
colnames(annotations_Amplicons)=c("Library","gene_id","Chr","Start","End") #annotations_Amplicons[1:2,]

annotations_Amplicons=annotations_Amplicons[rownames(data),]  #sorting annotations according to data
#dim(annotations_Amplicons) #dim(data)

#1.2 Filtering unefficient probes
#hist(rowSums(data_couv)/ncol(data),breaks=100) #average number of X for each amplicon
cutoff=30 #by default, minimum of 30x is requiered for coverage.. 
ProbesToFlag=names(which(rowSums(data_couv)<cutoff*ncol(data)))
ProbesToKeep=setdiff(rownames(data_couv),ProbesToFlag)
data_couv=data_couv[ProbesToKeep,] 
annotations_Amplicons=annotations_Amplicons[ProbesToKeep,]

#1.3 Filtering samples (each of the 2 libraries for each sample must reach a minimal coverage)
librairie1=rownames(annotations_Amplicons)[which(annotations_Amplicons$Library==1)] #amplicon names from library1
librairie2=rownames(annotations_Amplicons)[which(annotations_Amplicons$Library==2)] #amplicon names from library1
genes=unique(as.character(annotations_Amplicons[,"gene_id"]))
SamplesToFlag=unique(c(names(which(colMeans(data_couv[librairie1,])<cutoff)),names(which(colMeans(data_couv[librairie2,])<cutoff))))
SamplesToKeep=setdiff(colnames(data_couv),SamplesToFlag)
data_couv=data_couv[,SamplesToKeep]
echantillons=colnames(data_couv)

#1.4 Equilibring the 2 librairies for each sample
ratio_librairies=colMeans(data_couv[librairie1,])/colMeans(data_couv[librairie2,])
data_couv_new=data_couv
if (mean(ratio_librairies)>1) {
    data_couv_new[librairie1,]=data_couv[librairie1,]
    data_couv_new[librairie2,]=t(apply(data_couv[librairie2,],1,function(x) x*ratio_librairies))
    } else {
data_couv_new[librairie1,]=t(apply(data_couv[librairie1,],1,function(x) x/ratio_librairies))
data_couv_new[librairie2,]=data_couv[librairie2,]
    }
data_couv=cbind(annotations_Amplicons,data_couv_new)

######################Setting up a baseline for each sample based on CN
##################################
#2. relative CN ratios between genes:
#algo: ratios are expressed relative to a baseline, generated from baseline genes
#baseline genes are identified as components of the largest set of genes with no singificant CN difference (after testing all possible pairs of genes). 

liste_genes=as.character(unique(data_couv$gene_id))
CNRelatif=matrix(nrow=length(liste_genes),ncol=length(echantillons))
dimnames(CNRelatif)=list(liste_genes,echantillons)
n_genes=length(liste_genes)
matrice_gene=matrix(nrow=n_genes,ncol=n_genes)
colnames(matrice_gene)=liste_genes
rownames(matrice_gene)=liste_genes
MeanSdRefEchantillons=matrix(nrow=length(echantillons),ncol=2)
colnames(MeanSdRefEchantillons)=c("MeanReads","SdReads")
rownames(MeanSdRefEchantillons)=echantillons
CNSd=CNRelatif
CNMean=CNRelatif
CNMedianRelatif=CNRelatif

for (echantillon in echantillons)
	{
	#echantillon=echantillons[1]
	for (gene in liste_genes)
        	{
        	pos_gene=which(data_couv[,"gene_id"]==gene)
	  	for (gene2 in liste_genes)
			{
			pos_gene2=which(data_couv[,"gene_id"]==gene2)
			ThereIsEnoughData=(length(unique(data_couv[pos_gene,echantillon]))>1)&(length(unique(data_couv[pos_gene2,echantillon]))>1)
			if (ThereIsEnoughData) {
				if (t.test(data_couv[pos_gene,echantillon],data_couv[pos_gene2,echantillon])$p.value<0.05)
					{
					matrice_gene[gene,gene2]=mean(data_couv[pos_gene,echantillon],na.rm=T)/mean(data_couv[pos_gene2,echantillon],na.rm=T)
					} else {
					matrice_gene[gene,gene2]=1
					}
				} else {
				matrice_gene[gene,gene2]=NA
				}
			}
	 	}
	#identify the maximal number of genes equal
	matrice_gene[which(matrice_gene!=1)]=0
	matrice_gene[which(is.na(matrice_gene))]=0
	baseline=which(rowSums(matrice_gene,na.rm=T)==max(rowSums(matrice_gene,na.rm=T),na.rm=T))[1]
	genes_baseline=names(which(matrice_gene[baseline,]==1))
	#determine the CN relative to baseline
	#filtering out outliers (0.25-0.75)
	pos_temp=which(data_couv[,"gene_id"]%in%genes_baseline)
	data_temp=data_couv[pos_temp,echantillon]
	data_temp=data_temp[which((data_temp>quantile(data_temp,0.25))&(data_temp<quantile(data_temp,0.75)))]
	MeanSdRefEchantillons[echantillon,"MeanReads"]=mean(data_temp,na.rm=T)
	MeanSdRefEchantillons[echantillon,"SdReads"]=sd(data_temp,na.rm=T)
	for (gene in liste_genes)
        	{
		pos_gene=which(data_couv[,"gene_id"]==gene)
		#filtering out outliers (0.1-0.9)
		pos_temp=which(findInterval(data_couv[pos_gene,echantillon],quantile(data_couv[pos_gene,echantillon],c(0.1,0.9),na.rm=T))==1)
		pos_temp=pos_gene[pos_temp]
		if (length(pos_temp)<2) pos_temp=pos_gene #protecting from error in case there are very few amplicons
		CNSd[gene,echantillon]=sd(data_couv[pos_temp,echantillon],na.rm=T)
		CNMean[gene,echantillon]=mean(data_couv[pos_temp,echantillon],na.rm=T)
		if (length(which(genes_baseline==gene))==0)
			{
	        	CNRelatif[gene,echantillon]=mean(data_couv[pos_gene,echantillon],na.rm=T)/mean(data_temp)
			} else {
			CNRelatif[gene,echantillon]=1
			}
		CNMedianRelatif[gene,echantillon]=median(data_couv[pos_gene,echantillon],na.rm=T)/mean(data_temp)
		}
	plot(1,1,type='n')
	text(1,1,paste(round(which(echantillons==echantillon)*100/length(echantillons)),'%'))
	}


##########################
#2. Detection of Homozygous deletions and amplicons: by gene, by subregions
#Algo: based only on read counts: 1. Whole gene analysis: Detection of genes with very low (or very high) read counts compared to others; 2. Sub-genic analysis: detection of gene regions with very low/high read counts

# Parameters
NConsecutiveAmplicons=9 ##9 is maximal value (>9: bug) ##number of amplicons defining the minimal gene sub-region
threshold_ratio=0.3 #threshold for calling a CN very low (relative to baseline CN)
threshold_HighAmplif=3 #threshold for calling a CN very high (relative to baseline)

##Whole gene analysis
GenesHomoDel=matrix(nrow=length(echantillons),ncol=length(genes))
rownames(GenesHomoDel)=echantillons
colnames(GenesHomoDel)=genes
GenesHomoDel[,]=0
GenesHighAmplif=GenesHomoDel

segments_temp_genes=matrix(nrow=0,ncol=4) #Col 4: 0, delhomo; 1: amplif

for (echantillon in echantillons)
	{
	#echantillon="ACC110"
	for (gene in genes)
		{
#		#gene="ZNRF3"
		pos_gene=which(data_couv[,"gene_id"]==gene)
		if (!is.na(CNSd[gene,echantillon]))
			{
			if (!is.na(MeanSdRefEchantillons[echantillon,"SdReads"]))
				{
				if ((CNRelatif[gene,echantillon]<threshold_ratio)&(CNSd[gene,echantillon]<MeanSdRefEchantillons[echantillon,"SdReads"]/1.5)) 
					{
					GenesHomoDel[echantillon,gene]=1
					segments_temp_genes=rbind(segments_temp_genes,c(min(pos_gene),max(pos_gene),echantillon,0))
					}
				}
			if (CNRelatif[gene,echantillon]>threshold_HighAmplif*1.5) 
				{
				GenesHighAmplif[echantillon,gene]=1
				segments_temp_genes=rbind(segments_temp_genes,c(min(pos_gene),max(pos_gene),echantillon,1))
				}
			}
		}
	}



#Sub-genic analysis (detection of small regions within genes with very low read counts)
segments_temp=matrix(nrow=0,ncol=3)
for (echantillon in echantillons)
	{
	#echantillon="ACC110"
	donnees=na.omit(data_couv[,echantillon])
	for (gene in genes)
		{
		#gene="ZNRF3"
		if (GenesHomoDel[echantillon,gene]==0)
			{
			pos_gene=which(data_couv[,"gene_id"]==gene)
			donnees=na.omit(data_couv[pos_gene,echantillon])
			pos_delhomo=which(donnees/MeanSdRefEchantillons[echantillon,"MeanReads"]<threshold_ratio)
			donnees_bin=rep(0,length(donnees))
			donnees_bin[pos_delhomo]=1
			if (length(pos_delhomo)>0)
				{
				#selecting regions with at least NConsecutiveAmplicons outliers
				donnees_bin_n=c(donnees_bin,rep(0,NConsecutiveAmplicons-1))
				for (n in 2:NConsecutiveAmplicons)
					{
					#n=3
					donnees_bin_n=donnees_bin_n+c(rep(0,n-1),donnees_bin,rep(0,NConsecutiveAmplicons-n))
					}	
				if (max(donnees_bin_n)>=NConsecutiveAmplicons)
					{
					string=paste(as.character(donnees_bin_n),collapse='')			
					pattern1=paste(as.character(c(NConsecutiveAmplicons,NConsecutiveAmplicons-1)),collapse='')
					pattern2=paste(as.character(c(NConsecutiveAmplicons-1,NConsecutiveAmplicons)),collapse='')
					pattern=paste(as.character(c(NConsecutiveAmplicons,NConsecutiveAmplicons)),collapse='')
					string=gsub(pattern1,pattern,string)
					string=gsub(pattern2,pattern,string)
					string=gsub(as.character(NConsecutiveAmplicons),"A",string)
					#positions=as.numeric(gregexpr(as.character(NConsecutiveAmplicons),string)[[1]])
					positions=as.numeric(gregexpr("A",string)[[1]])
					if (length(positions)>0)
						{
						#finding the break points
						string=gsub("[[:digit:]]","B",string)
						starts=pos_gene[as.numeric(gregexpr("BA",string)[[1]])-NConsecutiveAmplicons+3]
						ends=pos_gene[as.numeric(gregexpr("AB",string)[[1]])-1]
						results_temp=cbind(starts,ends,rep(echantillon,length(starts)))
						}
					for (i in 1:nrow(results_temp))
						{
						#i=4
						pos=as.numeric(results_temp[i,"starts"]):as.numeric(results_temp[i,"ends"])
						donnees=na.omit(data_couv[pos,echantillon])
						#removing outliers
						donnees=donnees[which((donnees>quantile(donnees,0.1))&(donnees<quantile(donnees,0.9)))]    ###|(mean(donnees)/MeanSdRefEchantillons[echantillon,"MeanReads"]<0.15)
						GenesHomoDel[echantillon,gene]=1 
						segments_temp=rbind(segments_temp,results_temp[i,])
						}
					}
				}
			}	
		}
	}


#Identification of genes with high level amplification (whole gene analysis only; no sub-genic analysis is implemented because of low biological meaning)
segments_HighAmplif_NGS=matrix(nrow=0,ncol=7)
colnames(segments_HighAmplif_NGS)=c("start","end","sample","chr","Pos_Start_37","Pos_End_37","gene")
pos_amplif=which(segments_temp_genes[,4]==1)
if (length(pos_amplif)>1) {
	segments_temp_high=segments_temp_genes[pos_amplif,1:3]
	segments_HighAmplif_NGS=as.data.frame(cbind(segments_temp_high,data_couv[as.numeric(segments_temp_high[,1]),"Chr"],data_couv[as.numeric(segments_temp_high[,1]),"Start"],data_couv[as.numeric(segments_temp_high[,1]),"End"],as.character(data_couv[as.numeric(segments_temp_high[,1]),"gene_id"])))
	colnames(segments_HighAmplif_NGS)=c("start","end","sample","chr","Pos_Start_37","Pos_End_37","gene")
	segments_HighAmplifHighAmplif_NGS=segments_HighAmplif_NGS[order(as.numeric(as.character(segments_HighAmplif_NGS[,"start"]))),]
	segments_HighAmplif_NGS=segments_HighAmplif_NGS[order(as.character(segments_HighAmplif_NGS[,"sample"])),]
	} 
if (length(pos_amplif)==1) {
	segments_temp_high=segments_temp_genes[pos_amplif,1:3]
	segments_HighAmplif_NGS=c(segments_temp_high,data_couv[as.numeric(segments_temp_high[1]),"Chr"],data_couv[as.numeric(segments_temp_high[1]),"Start"],data_couv[as.numeric(segments_temp_high[2]),"End"],as.character(data_couv[as.numeric(segments_temp_high[1]),"gene_id"]))
	names(segments_HighAmplif_NGS)=c("start","end","sample","chr","Pos_Start_37","Pos_End_37","gene")
	segments_HighAmplif_NGS=t(as.matrix(segments_HighAmplif_NGS))
	} 



#merging delhomo segments with entiere genes
segments_temp_genes_delhomo=segments_temp_genes[which(segments_temp_genes[,4]==0),1:3]
segments_temp=rbind(segments_temp_genes_delhomo,segments_temp)
segments_delhomo_NGS=as.data.frame(cbind(segments_temp,data_couv[as.numeric(segments_temp[,1]),"Chr"],data_couv[as.numeric(segments_temp[,1]),"Start"],data_couv[as.numeric(segments_temp[,2]),"End"],as.character(data_couv[as.numeric(segments_temp[,1]),"gene_id"])))
colnames(segments_delhomo_NGS)=c("start","end","sample","chr","Pos_Start_37","Pos_End_37","gene")
segments_delhomo_NGS=segments_delhomo_NGS[order(as.numeric(as.character(segments_delhomo_NGS[,"start"]))),]
segments_delhomo_NGS=segments_delhomo_NGS[order(as.character(segments_delhomo_NGS[,"sample"])),]


#Merging duplicates (e.g. segments detected both as whole gene and sub-gene regions)
liste=paste(segments_delhomo_NGS[,"sample"],segments_delhomo_NGS[,"gene"])
posToKeep=NULL
for (i in unique(liste))
	{
	#i=unique(liste)[1]
	pos=which(liste==i)
	posToKeep=c(posToKeep,pos[1])
	}
segments_delhomo_NGS=segments_delhomo_NGS[posToKeep,]




###3. Calling gains and losses of one DNA copy by combining allelic ratios and read counts

#3.1 reading VCF. Structure of file should respect :"chrom	start	end	ref	alt	genotype	Nature	N_ref	N_var	N_var_Sum	N_tot	rs	sample	mergeID	N_ref_L	N_var_L

data_vcf=read.csv("Example_VCF_SnpHet_TumorLeuco_160323.csv",as.is=T , sep = ",")
posToKeep=which(!is.na(data_vcf$ref))
data_vcf=data_vcf[posToKeep,]

#Preparing the list of segments (here 1 segment= 1 gene)
liste_segments_NGS=matrix(nrow=length(genes),ncol=3)
colnames(liste_segments_NGS)=c("Chr","Start","End")
rownames(liste_segments_NGS)=genes

for (gene in genes)
	{
	#gene="RPL22"
	pos_gene=which(data_couv[,2]==gene)
	liste_segments_NGS[gene,]=c(as.numeric(as.character(data_couv[min(pos_gene),3])),as.numeric(as.character(data_couv[min(pos_gene),4])),as.numeric(as.character(data_couv[max(pos_gene),5])))
	}

#adding gene names to vcf
gene_temp=rep(NA,nrow(data_vcf))
for (gene in genes)
	{
	#gene="TP53"
	Chr=liste_segments_NGS[gene,"Chr"]
	if (Chr==23) Chr="X"
	Debut=liste_segments_NGS[gene,"Start"]
	Fin=liste_segments_NGS[gene,"End"]
	pos_gene=which((data_vcf$chrom==paste("chr",Chr,sep=''))&(data_vcf$start>=Debut)&(data_vcf$end<=Fin))
	gene_temp[pos_gene]=gene
	}	
gene=gene_temp
data_vcf=cbind(data_vcf,gene)

#MBAF
MBAF=(abs(0.5-data_vcf$N_var/(data_vcf$N_var+data_vcf$N_ref))+0.5)
data_vcf=cbind(data_vcf,MBAF)

#checking identity of sample names between table_vcf AND table
data_vcf$sample 
echantillons
#echantillons_sorted=c("ACC106","ACC91","ACC100","ACC125","ACC72","ACC88","ACC63","ACC80","ACC57","ACC81","ACC84","ACC103","ACC87","ACC59","ACC65","ACC69","ACC74","ACC82","ACC93","ACC94","ACC120","ACC124","ACC101","ACC109","ACC86","ACC67","ACC73","ACC89","ACC97","ACC66","ACC79","ACC83","ACC110","ACC112" ,"ACC115" ,"ACC116","ACC122","ACC58","ACC61","ACC76","ACC99","ACC104","ACC107","ACC70","ACC75","ACC92","ACC98","ACC113","ACC132","ACC111","ACC95","ACC129","ACC130","ACC60","ACC64","ACC68","ACC71","ACC102")
echantillons_sorted=echantillons
#,"ACC118","ACC123","ACC105","ACC114","ACC121","ACC126","ACC128","ACC56","ACC108","ACC131","ACC90","ACC78","ACC85","ACC62","ACC127")
#samples to discard (bec in vcf but not in CN ; 2 reasons: data missing for 4 (ongoing); Failure of 1 library (4 samples)
samplesToDiscard=unique(data_vcf$sample[which(!data_vcf$sample%in%echantillons)])
posToKeep=which(!data_vcf$sample%in%samplesToDiscard)
data_vcf=data_vcf[posToKeep,]

###setting up baseline: a more accurate definition of baseline is possible, compared to the previous step based on CN. This time we use the allelic ratio to define 2 copies of DNA, assuming that it will correspond to CN associated with segments showing MBAF of 0.5 (normal heterozygous)
#algo: based on obvious hets (MBAF<0.6) and obvious LOH (MBAF>0.75)
#If there are hets: Fitting a grid diploid in the CN / BAF plot; Nota: testing hypothesis of tetraploid grid was abandonned bec CN not reliable enough.
#If there are no hets: fit a grid where LOH is deletion of 1 copy. Nota: testing hypothesis of CN LOH was abandonned bec CN not reliable enough.  
#Adjusting CN baseline to minimize the distance according to the best fit



seuil_Loh=0.75 #threshold beyond which MBAF is considered to correspond to AI/LOH 
seuil_noLoh=0.6 #threshold before which MBAF is considered heterozygous and allelic ratio is balanced.

CNMedianRelatif_BaselineOptimized=CNMedianRelatif


for (echantillon in echantillons)
	{
	#echantillon="T152"
	CN_MBAF_genes=matrix(nrow=length(genes),ncol=6)
	rownames(CN_MBAF_genes)=genes
	colnames(CN_MBAF_genes)=c("MBAF","CN","DistDipl","DistTetra","DistCnLohAndLoss","DistCnOnly")
	for (gene in genes)
		{
		#gene="BAI1"
		CN_MBAF_genes[gene,"CN"]=CNMedianRelatif[gene,echantillon]
		pos_gene_vcf=which((data_vcf$sample==echantillon)&(data_vcf$gene==gene))
		if (length(pos_gene_vcf)>0) {
			#case with germline het snps
			CN_MBAF_genes[gene,"MBAF"]=median(data_vcf$MBAF[pos_gene_vcf],na.rm=T)
			}
		}


	#Identification of SNPs with noLoh and with Loh
	noLoh=names(which(CN_MBAF_genes[,"MBAF"]<seuil_noLoh) )
	Loh=names(which(CN_MBAF_genes[,"MBAF"]>seuil_Loh))
	ThereAreSnpsWithNoLoh=(length(noLoh)>0) 
	if (ThereAreSnpsWithNoLoh) {
		#cases of hets with no LOH
		#correctif pour CN het: on essaie de fitter SNPs with no LOH à CN=1
		ratio_y1=1/CN_MBAF_genes[noLoh,"CN"]
		} else {
		#case with no SNP with LOH
		ratio_y1=NA
		}
	ThereAreSnpsInLoh=(length(Loh)>0)
	if (ThereAreSnpsInLoh) {
		x3=CN_MBAF_genes[Loh,"MBAF"]
		names(x3)=Loh
		ytheorique=1.5-x3
		ratio_y3=ytheorique/CN_MBAF_genes[Loh,"CN"]
		} else {
		#case with no LOH
		ratio_y3=NA
		}
	#patch: test de fonction de delta_y en fonction du correctif a
	#fonction_delta=NULL
	#echelle=seq(0,2,by=0.01)
	#for (a in echelle )
	#	{
	#	fonction_delta=c(fonction_delta,median(c(abs(1-CN_MBAF_genes[noLoh,"CN"]*a),abs(ytheorique-CN_MBAF_genes[Loh,"CN"]*a))))
		#fonction_delta=c(fonction_delta,sum(c(abs(1-CN_MBAF_genes[noLoh,"CN"]*a),abs(ytheorique-CN_MBAF_genes[Loh,"CN"]*a))))
	#	}
	#plot(echelle,fonction_delta,type='l')
	#names(fonction_delta)=echelle
	#which(fonction_delta==min(fonction_delta))
	#conclusion: with median (instead of sum), outliers are impacting less. Very close to median 
	#end of patch

	ratio_y=na.omit(c(ratio_y1,ratio_y3))
	correctif=median(ratio_y)
	CNMedianRelatif_BaselineOptimized[,echantillon]=CNMedianRelatif[,echantillon]*correctif
	}
		


#analysis (CN / MBAF)
#structure of outpout: 
#by gene CN,MBAF,Call,Call_distance,N_SNPhet,Namplicon
#by sample: coverage_status(excellent,good,fair,limit,unsufficient),library_balance,noise_status,
ChrAlterations=matrix(nrow=length(echantillons)*length(genes),ncol=9)
colnames(ChrAlterations)=c("Sample","Gene","CN","MBAF","Call","Call_distance","Proportion_of_cells","N_amplicons","N_snpHet")
i=0
PropTumorCells_all=rep(NA,length(echantillons))
names(PropTumorCells_all)=echantillons
for (echantillon in echantillons)
	{
	PropTumorCells=NULL
	#echantillon="ACC81"
	for (gene in genes)
		{
		i=i+1
		#gene="RNF43"
		ChrAlterations[i,"Sample"]=echantillon
		ChrAlterations[i,"Gene"]=gene
		y=CNMedianRelatif_BaselineOptimized[gene,echantillon]
		if (is.na(y)) y=CNMedianRelatif[gene,echantillon]
		ChrAlterations[i,"CN"]=y
		ChrAlterations[i,"N_amplicons"]=length(which(data_couv[,"gene_id"]==gene))
		pos_gene_vcf=which((data_vcf$sample==echantillon)&(data_vcf$gene==gene))
		if (length(pos_gene_vcf)>0) {
			#case with het snps
			x=median(data_vcf$MBAF[pos_gene_vcf],na.rm=T)
			ChrAlterations[i,"MBAF"]=x
			ChrAlterations[i,"N_snpHet"]=length(pos_gene_vcf)
			#identification of obvious LOH
			if (x>0.7)
				{
				#obvious LOH
				PropTemp=(x-0.5)*2
				PropTumorCells=c(PropTumorCells,PropTemp)
				ChrAlterations[i,"Proportion_of_cells"]=PropTemp
				y_DelExpected=1-0.5*PropTemp
				ChrAlterations[i,"Call"]="Chromosome Loss"
				ChrAlterations[i,"Call_distance"]=y-y_DelExpected
				}
			#obvious hets
			if (x<0.6) {
				dist_DiplHet=sqrt((x-0.5)^2+(y-1)^2)
				ChrAlterations[i,"Call"]="Diploid_het"
				ChrAlterations[i,"Call_distance"]=dist_DiplHet
				}
			#other cases (MBAF between 0.6 and 0.7) #discriminating between gain, CN LOH and deletions
			if ((x>=0.6)&(x<=0.7)) {
				PropTempDel=(x-0.5)*2
				PropTempGain=(x-0.5)*6   
				dist_Gain=abs(y-0.5*PropTempGain-1)
				dist_Loss=abs(y+0.5*PropTempDel-1)
				if (dist_Gain<dist_Loss) {
					ChrAlterations[i,"Call"]="Chromosome Gain"
					ChrAlterations[i,"Call_distance"]=dist_Gain
					#PropTumorCells=c(PropTumorCells,PropTempGain) not considered bec not precise
					ChrAlterations[i,"Proportion_of_cells"]=PropTempGain
					}
				if (dist_Loss<dist_Gain) {
					ChrAlterations[i,"Call"]="Chromosome Loss"
					ChrAlterations[i,"Call_distance"]=dist_Loss
					PropTumorCells=c(PropTumorCells,PropTempDel)
					ChrAlterations[i,"Proportion_of_cells"]=PropTempDel
					}
				}
			} else {
			#case with no het snp
			ChrAlterations[i,"N_snpHet"]=0
			#calling only based on CN, with 4 options: "4, 3, 2, 1 DNA copies
			cn=as.numeric(ChrAlterations[i,"CN"])
			if (cn>=2.25) ChrAlterations[i,"Call"]=">4 DNA copies?"
			if (cn<0.25) ChrAlterations[i,"Call"]="0 DNA copy?"
			if ((cn>=1.75)&(cn<2.25)) ChrAlterations[i,"Call"]="4 DNA copies?"
			if ((cn>=0.25)&(cn<0.75)) ChrAlterations[i,"Call"]="1 DNA copy?"
			if ((cn>=0.75)&(cn<1.25)) ChrAlterations[i,"Call"]="2 DNA copies?"
			if ((cn>=1.25)&(cn<1.75)) ChrAlterations[i,"Call"]="3 DNA copies?"
			}
		}
	if (length(PropTumorCells)>0) PropTumorCells_all[echantillon]=max(PropTumorCells,na.rm=T)
	}




#Adding a Chr, start and a end for each ChrAlteration
pos_temp=matrix(nrow=nrow(ChrAlterations),ncol=5)
colnames(pos_temp)=c("Chr","Start_37","End_37","Start_Ind","End_Ind")
for (gene in genes) {
	pos_gene=which(ChrAlterations[,"Gene"]==gene)
	pos_gene2=which(data_couv[,"gene_id"]==gene)
	start_ind=pos_gene2[1]
	end_ind=pos_gene2[length(pos_gene2)]
	pos_temp[pos_gene,1]=data_couv[start_ind,"Chr"]
	pos_temp[pos_gene,2]=data_couv[start_ind,"Start"]
	pos_temp[pos_gene,3]=data_couv[end_ind,"End"]
	pos_temp[pos_gene,4]=start_ind
	pos_temp[pos_gene,5]=end_ind
	}
ChrAlterations=cbind(ChrAlterations,pos_temp)

#Adding homozygous del
delhomo_temp=matrix(nrow=nrow(segments_delhomo_NGS),ncol=ncol(ChrAlterations))
colnames(delhomo_temp)=colnames(ChrAlterations)
delhomo_temp[,"Gene"]=as.character(segments_delhomo_NGS[,"gene"])
delhomo_temp[,"Sample"]=as.character(segments_delhomo_NGS[,"sample"])
for (i in 1:nrow(segments_delhomo_NGS))
	{
	gene=delhomo_temp[i,"Gene"]
	echantillon=delhomo_temp[i,"Sample"]
	delhomo_temp[i,"CN"]=CNMedianRelatif_BaselineOptimized[gene,echantillon]
	}
delhomo_temp[,"Call"]="Homozygous Deletion"
delhomo_temp[,"N_amplicons"]=as.numeric(as.character(segments_delhomo_NGS[,"end"]))-as.numeric(as.character(segments_delhomo_NGS[,"start"]))+1
delhomo_temp[,"Chr"]=as.character(segments_delhomo_NGS[,"chr"])
delhomo_temp[,"Start_37"]=as.character(segments_delhomo_NGS[,"Pos_Start_37"])
delhomo_temp[,"End_37"]=as.character(segments_delhomo_NGS[,"Pos_End_37"])
delhomo_temp[,"Start_Ind"]=as.character(segments_delhomo_NGS[,"start"])
delhomo_temp[,"End_Ind"]=as.character(segments_delhomo_NGS[,"end"])
ChrAlterations=rbind(ChrAlterations,delhomo_temp)

#Adding high-level amplicons
HighAmplif_temp=matrix(nrow=nrow(segments_HighAmplif_NGS),ncol=ncol(ChrAlterations))
colnames(HighAmplif_temp)=colnames(ChrAlterations)
HighAmplif_temp[,"Gene"]=as.character(segments_HighAmplif_NGS[,"gene"])
HighAmplif_temp[,"Sample"]=as.character(segments_HighAmplif_NGS[,"sample"])
for (i in 1:nrow(segments_HighAmplif_NGS))
	{
	gene=HighAmplif_temp[i,"Gene"]
	echantillon=HighAmplif_temp[i,"Sample"]
	HighAmplif_temp[i,"CN"]=CNMedianRelatif_BaselineOptimized[gene,echantillon]
	}
HighAmplif_temp[,"Call"]="High Level Amplification"
HighAmplif_temp[,"N_amplicons"]=as.numeric(as.character(segments_HighAmplif_NGS[,"end"]))-as.numeric(as.character(segments_HighAmplif_NGS[,"start"]))+1
HighAmplif_temp[,"Chr"]=as.character(segments_HighAmplif_NGS[,"chr"])
HighAmplif_temp[,"Start_37"]=as.character(segments_HighAmplif_NGS[,"Pos_Start_37"])
HighAmplif_temp[,"End_37"]=as.character(segments_HighAmplif_NGS[,"Pos_End_37"])
HighAmplif_temp[,"Start_Ind"]=as.character(segments_HighAmplif_NGS[,"start"])
HighAmplif_temp[,"End_Ind"]=as.character(segments_HighAmplif_NGS[,"end"])
ChrAlterations=rbind(ChrAlterations,HighAmplif_temp)




ChrAlterations=ChrAlterations[order(as.numeric(ChrAlterations[,"Start_37"])),]
ChrAlterations=ChrAlterations[order(as.numeric(ChrAlterations[,"Chr"])),]
ChrAlterations=ChrAlterations[order(ChrAlterations[,"Sample"]),]

#Suppressing duplicates delhomo et amplif
pos_DelHomo=which(ChrAlterations[,"Call"]=="Homozygous Deletion")
zero=paste(ChrAlterations[,"Sample"],ChrAlterations[,"Gene"])
delhomo=paste(ChrAlterations[pos_DelHomo,"Sample"],ChrAlterations[pos_DelHomo,"Gene"])
pos_ToFlag=which(zero%in%delhomo)
pos_ToFlag=setdiff(pos_ToFlag,pos_DelHomo)
pos_ToKeep=setdiff(seq(1,nrow(ChrAlterations)),pos_ToFlag)
ChrAlterations=ChrAlterations[pos_ToKeep,]

pos_HighAmplif=which(ChrAlterations[,"Call"]=="High Level Amplification")
HighAmplif=paste(ChrAlterations[pos_HighAmplif,"Sample"],ChrAlterations[pos_HighAmplif,"Gene"])
zero=paste(ChrAlterations[,"Sample"],ChrAlterations[,"Gene"])
pos_ToFlag=which(zero%in%HighAmplif)
pos_ToFlag=setdiff(pos_ToFlag,pos_HighAmplif)
pos_ToKeep=setdiff(seq(1,nrow(ChrAlterations)),pos_ToFlag)
ChrAlterations=ChrAlterations[pos_ToKeep,]

#write.csv(ChrAlterations,"zChrAlterations_Example_161012.csv")
	
#plotting
limits_genes=NULL
for (gene in genes)
	{
	limits_genes=c(limits_genes,max(which(data_couv[,"gene_id"]==gene)))
	}
names(limits_genes)=genes
limits_genes=c(0,limits_genes)
zoom=4 #zoom=1

cercle<-function(a,b)
	{#circle
	for (j in 1:25)
		{
		par(new=T)
		symbols(a,b, circles = 0.1-j/250, inches = FALSE,xlim=limx,ylim=limy,fg=gray(1-j/40),bg=gray(1-j/40),xaxt='n',yaxt='n',xlab='',ylab='')
		}
	}
ellipse<-function(xa,ya,xb,yb,pas,r)
	{#circle
	for (j in 1:25)
		{
		par(new=T)
		symbols(xa,ya, circles = r*(1-j/25/pas), inches = FALSE,xlim=limx,ylim=limy,fg=gray(1-j/40),bg=gray(1-j/40),xaxt='n',yaxt='n',xlab='',ylab='')
		rect(xa-r*(1-j/25/pas),ya,xb+r*(1-j/25/pas),yb,border=gray(1-j/40),col=gray(1-j/40))
		par(new=T)
		symbols(xb,yb, circles = r*(1-j/25/pas), inches = FALSE,xlim=limx,ylim=limy,fg=gray(1-j/40),bg=gray(1-j/40),xaxt='n',yaxt='n',xlab='',ylab='')
		}
	}



rectangle<-function(a,b,c,d)
	{
	for (j in 1:25)
		{
		rect(a+j/200,b+j/100,c-j/200,d-j/100, col=gray(1-j/40),border=gray(1-j/40))
		}
	}
polygone1<-function(xa,xb,xc,xd,ya,yb,yc,yd,pasX,pasY)
	{
	for (j in 1:25)
		{
		a=(yc-yd)/(yb-ya)
		polygon(c(xa+j/pasX,xb+j/pasX,xc-j/pasX,xd-j/pasX),c(ya+j/pasY/a*0.7,yb-j/pasY/a*2.5,yc-j/pasY*1.5,yd+j/pasY*2), col=gray(1-j/40),border=gray(1-j/40))
		}
	}
polygone2<-function(xa,xb,xc,xd,ya,yb,yc,yd,pasX,pasY)
	{
	for (j in 1:25)
		{
		a=(yc-yd)/(yb-ya)
		polygon(c(xa+j/pasX,xb+j/pasX,xc-j/pasX,xd-j/pasX),c(ya-j/pasY,yb-2*j/pasY,yc-j/pasY,yd+3*j/pasY), col=gray(1-j/40),border=gray(1-j/40))
		}
	}
triangle<-function(xa,xb,xc,ya,yb,yc)
	{
	for (j in 1:25)
		{
		polygon(c(xa+j/500,xb-0.5*j/500,xc-0.5*j/500),c(ya+3*j/500,yb-7*j/500,yc+3*j/500), col=gray(1-j/40),border=gray(1-j/40))
		}
	}

for (echantillon in echantillons)
#for (echantillon in echantillons_sorted) #echantillon=echantillons_sorted[1]
	{	
	png(paste("z",echantillon,"_MbafCn_Example_","161013.png",sep=''),height=2000,width=2000)
	par(mar=c(zoom*2,zoom*3,zoom,0),mgp=c(6,2,0))
	limy=c(-1400,1800)
	limx=c(0,max(limits_genes))
	#CN
	#echantillon="ACC91" #echantillon="ACC100" # "ACC106" "ACC91"  "ACC100" "ACC125" "ACC72"  "ACC88"  "ACC63"  "ACC80"  "ACC57"  "ACC81"  "ACC84"  "ACC103" "ACC87"  "ACC59"  "ACC65"  "ACC69"  "ACC74"  "ACC82"  "ACC93"  "ACC94"
	#echantillon="ACC106"
	
	donnees=na.omit(data_couv[,echantillon])
	convert_y=max(donnees,na.rm=T)
	donnees[which(donnees>convert_y)]=convert_y
	plot(donnees/convert_y*1000,xlim=limx,ylim=limy,pch=20,axes=F,xlab=paste('MBAF',paste(rep(' ',50),collapse='')),ylab=paste('CN (relative to ploidy)',paste(rep(' ',20),collapse=''),'N Reads',paste(rep(' ',31),collapse=''),'MBAF'),main=echantillon,cex=zoom,cex.main=zoom,cex.axis=zoom,cex.lab=zoom)
	abline(h=c(0,1000,1300,1800),col='grey') 
   	abline(h=MeanSdRefEchantillons[echantillon,"MeanReads"]*threshold_ratio,col='red')
	abline(h=MeanSdRefEchantillons[echantillon,"MeanReads"]*threshold_HighAmplif,col='red')
	abline(h=MeanSdRefEchantillons[echantillon,"MeanReads"],col='green')
	#looking for homozygous deletions and high level amplification
	pos=which(((ChrAlterations[,"Call"]=="Homozygous Deletion")|(ChrAlterations[,"Call"]=="High Level Amplification"))&(ChrAlterations[,"Sample"]==echantillon))
	if (length(pos)>0) {
		for (i in 1:length(pos))
			{
			pos_i=pos[i]
			range_i=as.numeric(ChrAlterations[pos_i,"Start_Ind"]):as.numeric(ChrAlterations[pos_i,"End_Ind"])
		    	par(new=T)
			couleur="red"
			plot(range_i,donnees[range_i]/convert_y*1000,xlim=limx,ylim=limy,col=couleur,pch=20,axes=F,xlab='',ylab='',cex=zoom)
			}
		}
	nbre_couleurs=15
	couleurs=rainbow(nbre_couleurs)
	rang_couleur=rep(sample(seq(1,nbre_couleurs),nbre_couleurs),10)
	i=0
	for (gene in genes)
		{
		i=i+1
		lines(rep(limits_genes[gene],2),c(1300,1800),col='grey')
		lines(rep(limits_genes[gene],2),c(0,1000),col='grey')
		text((limits_genes[i]+limits_genes[i+1])/2,1000,gene,srt=90,col=couleurs[rang_couleur[i]],adj=c(1,0.5),cex=zoom)
		}
	axis(2,c(0,500,1000),c(0,round(500*convert_y/1000),round(1000*convert_y/1000)),cex=zoom,cex.axis=zoom)


	#MBAFplot 
	
	axis(2,c(1300,1400,1500,1600,1700,1800),seq(0.5,1,by=0.1),cex=zoom,cex.axis=zoom)
	MBAF_all=NULL #by het SNP
	MBAF=NULL #by gene
	CN=NULL
	pos_het_ech=which((ChrAlterations[,"Sample"]==echantillon)&(!is.na(ChrAlterations[,"MBAF"])))
	if (length(pos_het_ech)>0) {
		for (gene in genes)
			{
			pos_MBAF=which((data_vcf[,"gene"]==gene)&(data_vcf[,"sample"]==echantillon))
			if (length(pos_MBAF)>0) {
				MBAF_temp=as.numeric(data_vcf[pos_MBAF,"MBAF"])
				names(MBAF_temp)=rep(gene,length(pos_MBAF))
				MBAF_all=c(MBAF_all,MBAF_temp)
				}
			}
		MBAF=as.numeric(ChrAlterations[pos_het_ech,"MBAF"])
		names(MBAF)=ChrAlterations[pos_het_ech,"Gene"]
		CN=as.numeric(ChrAlterations[pos_het_ech,"CN"])
		names(CN)=ChrAlterations[pos_het_ech,"Gene"]
		}
	#plotting MBAF
	for (gene in unique(names(MBAF_all)))
		{
		i=which(names(limits_genes)==gene)-1
		pos_gene=which(names(MBAF_all)==gene)
		xmini=limits_genes[i]
		xmaxi=limits_genes[i+1]
		if (length(pos_gene)==1) {
			x=(xmini+xmaxi)/2
			} else {
			delta=(xmaxi-xmini)/length(pos_gene)
			x=seq(xmini+delta/2,xmaxi-delta/2,by=delta)
			}
		par(new=T)
		plot(x,(MBAF_all[pos_gene]-0.5)*1000+1300,col=couleurs[rang_couleur[i]],pch=20,xlim=limx,ylim=limy,axes=F,xlab='',ylab='',cex=zoom)
		lines(c(min(x),max(x)),(rep(mean(MBAF_all[pos_gene],na.rm=T),2)-0.5)*1000+1300,col=couleurs[rang_couleur[i]])
		}

	#gap plot
	if (!is.null(CN)) {	
		limy=c(0,3*max(2,max(CN,na.rm=T)))
		} else {
		limy=c(0,6)
		}
	limx=c(0.4,1.5)
	par(new=T)
	plot(MBAF,CN,pch=20,xlim=limx,ylim=limy,col='grey',axes=F,xlab='',ylab='',cex=zoom)
	triangle(0.55,2/3,2/3,1.15,1.75,1.25)

	cercle(0.5,1)
	rect(0,0,0.5,max(limy)*5/12,border="white",col="white")

	for (j in seq(0.5,1,by=0.1))
		{
		lines(c(j,j),c(0.15,0.2),col='grey')
		text(j,0,paste((2*j-1)*100,"%",sep=""),pos=3,srt=90,col='grey',cex=zoom/2)
		}
	text(0.65,0.05,"%",font=3,col='grey',cex=zoom/2)
	text(0.75,0.05,"TUMOR",font=3,col='grey',cex=zoom/2)
	text(0.85,0.05,"CELLS",font=3,col='grey',cex=zoom/2)

	for (j in seq(0.5,2,by=0.5))
		{
		lines(c(0.4,0.42),c(j,j),col='grey')
		text(0.42,j,j*2,col='grey',pos=4,cex=zoom/2)
		lines(c(0.45,0.47),c(j,j),col='grey')
		}
	text(0.43,1.75,srt=270,"DNA",col='grey',font=3,cex=zoom/2)
	text(0.43,1.25,srt=270,"COPY",col='grey',font=3,cex=zoom/2)
	text(0.43,0.75,srt=270,"NUMBER",col='grey',font=3,cex=zoom/2)
	
	text(0.55,1.05,"heterozygous",col="white",font=3,cex=zoom/2)
	text(0.55,0.95,"diploid",col="white",font=3,cex=zoom/2)



	text(0.64,1.4,"Gain",col="white",font=3,cex=zoom/2)

	polygone1(0.6,0.6,1,1,0.85,1.05,1.25,0.25,150,90)
	text(0.83,0.9,"Loss",col="white",font=3,cex=zoom/2)
	


	for (gene in unique(names(MBAF)))
		{
		i=which(names(limits_genes)==gene)-1
		pos_gene=which(names(MBAF)==gene)
		par(new=T)
		plot(MBAF[pos_gene],CN[pos_gene],col=couleurs[rang_couleur[i]],pch=20,xlim=limx,ylim=limy,axes=F,xlab='',ylab='',cex=zoom)
#		par(new=T)
#		plot(MBAF[pos_gene],CN[pos_gene],col=couleurs[rang_couleur[i]],type='l',xlim=limx,ylim=limy,axes=F,xlab='',ylab='',cex=zoom)
		pitch=20
		n_pitch=1
		call=ChrAlterations[which((ChrAlterations[,"Sample"]==echantillon)&(ChrAlterations[,"Gene"]==gene)),"Call"]
		if ((call=="Diploid_het")|(call=="Copy-neutral LOH")|(call=="2 DNA copies?")) pitch=23
		if ((call=="Chromosome Gain")|(call=="3 DNA copies?")) pitch=24
		if ((call=="Homozygous Deletion")|(call=="Chromosome Loss")|(call=="1 DNA copy?")) pitch=25
		if ((call=="Tetraploid_het")|(call=="4 DNA copies?")|(call==">4 DNA copies?")) {
			pitch=24
			n_pitch=2
			}
		if ((call== "Homozygous Deletion")|(call=="0 DNA copy?")) {
			pitch=25
			n_pitch=2
			}
		if (n_pitch==1) points(MBAF[pos_gene],CN[pos_gene],pch=pitch,bg=couleurs[rang_couleur[i]],cex=2*zoom,col=couleurs[rang_couleur[i]])
		if (n_pitch==2)
			{
			points(MBAF[pos_gene],CN[pos_gene]-0.05,pch=pitch,bg=couleurs[rang_couleur[i]],cex=2*zoom,col=couleurs[rang_couleur[i]])
			points(MBAF[pos_gene],CN[pos_gene]+0.05,pch=pitch,bg=couleurs[rang_couleur[i]],cex=2*zoom,col=couleurs[rang_couleur[i]])
			}
 		y_sd=CNSd[gene,echantillon]/MeanSdRefEchantillons[echantillon,"MeanReads"]
		lines(rep(MBAF[pos_gene],2),c(max(0,mean(CN[pos_gene])-y_sd/2),min(max(CN),CN[pos_gene]+y_sd/2)),col=couleurs[rang_couleur[i]])
		pos_gene2=which(names(MBAF_all)==gene)
		x_sd=sd(MBAF_all[pos_gene2])/mean(MBAF_all[pos_gene2])
		lines(c(max(0,MBAF[pos_gene]-x_sd/2),min(max(MBAF),MBAF[pos_gene]+x_sd/2)),rep(CN[pos_gene],2),col=couleurs[rang_couleur[i]])
		}
	axis(1,c(0.5,0.6,0.7,0.8,0.9,1),c(0.5,0.6,0.7,0.8,0.9,1),cex=zoom,cex.axis=zoom)
	axis(2,c(0.5,1,1.5,2),c(0.5,1,1.5,2),cex=zoom,cex.axis=zoom)
	#metrix
	text(1.05,2,paste("% of tumor cells:",round(PropTumorCells_all[echantillon]*100),"%"),pos=4,cex=zoom)
	text(1.05,1.5,paste("Coverage mean:",round(median(CNMean[,echantillon],na.rm=T))),pos=4,cex=zoom)
	text(1.05,1,paste("Coverage coef of variation: ",round(median(CNSd[,echantillon]/CNMean[,echantillon],na.rm=T)*100),"%",sep=''),pos=4,cex=zoom)
	text(1.05,0.5,paste("Germline Heterozygous SNPs:",length(MBAF_all)),pos=4,cex=zoom)
	text(1.05,0,paste("Library ratios:",round(ratio_librairies[echantillon],1)),pos=4,cex=zoom)
	dev.off()
	}
