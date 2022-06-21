
# Import libraries
library(shiny); library(shinythemes); library(data.table); library(RCurl);library(DT);library(shinyjs)
require(stringdist); require(irr); require(xlsx);require(limma);require(pbapply); require(parallelDist); library(mgcv); require(psych); require(RColorBrewer);library(dplyr)
#source("../local_functions.R")
#library("devtools"); reload(pkg = "D:/AAAA/Projects/ROR", quiet = FALSE); library("ROR");

#Read data
# metadata=read.delim(paste0("../fast0.0/data/metadata_2022_05_18.tsv"))
# metadata=data.frame(metadata,row.names = "Accession.ID")
# sampleID=rownames(metadata)
# variant_pred=substr(metadata$Variant,5,11); variant_pred=strsplit2(variant_pred," ")[,1]
# variant_pred[variant_pred==""]="Zother"
# metadata=data.frame(metadata, variant_pred)
# save(metadata, file="data/pred_data_2022_05_18.Rdata")

# Loading data
# load(file="../fast0.0/data/training ID.Rdata")
# load(file="../fast0.0/result/full training result.Rdata")
load(file="data/pred_data_2022_05_18.Rdata")
load(file="D:/AAAA/Projects/COVID19-viral genomics/Fast/data/ref seq.Rdata"); rownames(ref_AA)=sub("-","_", rownames(ref_AA))
load(file="../fast0.0/workingfiles/hap cores.Rdata")

#save(ex_data, file="data/exampledata.Rdata")
#load(file="data/exampledata.Rdata"); pred_data=ex_data
#write.csv(pred_data, file="data/example.csv")
# AIS model
xN_limit=10000
# variant_prop is modifiable per population
# variant_corehap  is chosen by version
AIS<-function(pred_data, variant_prop, variant_corehap, m_parameters){
	prob_threshold=m_parameters["prob_threshold"]   # default=0.99
	min_prob=m_parameters["min_prob"]               # default=0.1
	pred_AA=pred_data$AA.Substitutions
	pred_AA=sub("\\(","", pred_data$AA.Substitutions); pred_AA=sub("\\)","", pred_AA)
	predID=rownames(pred_data); names(pred_AA)=predID
	variant_pred=pred_data$variant_pred
	if (is.null(variant_pred)) variant_pred=rep("na", length(pred_AA)); names(variant_pred)=predID
	variant_list=names(variant_corehap); variantN=length(variant_list); variant_list_s=setdiff(variant_list,"Zother")
	minMT=NULL; for (vj in variant_list) minMT=c(minMT, min(10, floor(0.5*length(variant_corehap[[vj]]$genes))))
	names(minMT)=variant_list

	pred_prob=lapply(pred_AA, function(x,y){
		tmpy=y[,"AA"]
		tmp=strsplit(x,",")[[1]]; tmp=sub("del","d", tmp); tmp=sub("stop","s", tmp)
		tmp1=nchar(tmp); tmpx=substr(tmp, 1, tmp1-1);
		tmpMut=substr(tmp, tmp1, tmp1);names(tmpMut)=tmpx
		tmpy[tmpx]=tmpMut[tmpx]
		
		jointp=NULL;tmpV=NULL; for (vj in variant_list){
			local_hap=tmpy[variant_corehap[[vj]]$genes]; 
			local_hap=paste(local_hap, collapse = "")
			hap_list=rownames(variant_corehap[[vj]]$freq)[variant_corehap[[vj]]$mutN>=floor(0.5*minMT[vj])]
			tmp=is.element(local_hap, hap_list)
			local_p=ifelse(tmp,variant_corehap[[vj]]$freq[local_hap,2],0)
			jointp=c(jointp, local_p)
			if (local_p>0) tmpV=c(tmpV,variant_corehap[[vj]]$freq[local_hap,1]) else tmpV=c(tmpV,0)
		}
		names(jointp)=variant_list; names(tmpV)=variant_list
		if (sum(jointp[variant_list_s])>0) jointp["Zother"]=0    # favor known variants
		tmpV=sum(tmpV["Zother"])
		if (tmpV<=1) {
			jointp=NULL; for (vj in variant_list){
				local_hap=tmpy[variant_corehap[[vj]]$genes];
				local_hap=paste(local_hap, collapse = "")
				lib_freq=cbind(variant_corehap[[vj]]$freq, MT=variant_corehap[[vj]]$mutN);
				lib_freq=lib_freq[lib_freq[,3]>=minMT[vj],]
				lib_dist=stringdist(rownames(lib_freq), local_hap, method="lv")     #Levenshtein distance
				if (sum(lib_dist<=2)==0) local_p=0 else local_p=max(lib_freq[lib_dist<=2,2])      # assumption: differ by 2 AA
				jointp=c(jointp, local_p)
			}
			names(jointp)=variant_list
			if (sum(jointp[variant_list_s])>0) jointp["Zother"]=0    # favor known variants
		}
		jointp=jointp*variant_prop
		jointp}, ref_AA)
	pred_prob=matrix(unlist(pred_prob), ncol=length(variant_list), byrow=T)
	colnames(pred_prob)=variant_list; rownames(pred_prob)=predID
	tmp=rowSums(pred_prob);  tmp0=min(tmp[tmp>0])*min_prob                          # assumption
	pred_prob[tmp==0,"Zother"]=tmp0
	
	pred_prob1=pred_prob;tmp=rowSums(pred_prob); for (j in colnames(pred_prob)) pred_prob1[,j]=pred_prob[,j]/tmp
	
	pred_variant=apply(pred_prob1, 1, function(x, var_names) {
		tmp=max(x); if (tmp<prob_threshold) output="zMixture" else output=var_names[x==tmp][1]   #assumption
		c(output, tmp)}, colnames(pred_prob1))
	tmp=matrix(pred_variant,ncol=2, byrow=T)
	pred_variant=data.frame(pred=tmp[,1], postp=as.numeric(tmp[,2]))
	rownames(pred_variant)=predID
	
	# post-prediction modification
	tmp_data=data.frame(pred_AA,pred_prob1[,variant_list_s])
	ppred=apply(tmp_data,1, function(x,y){
		target_variant=x[-1];
		target_variant=names(target_variant)[target_variant>=1-prob_threshold]; other_variant=setdiff(variant_list, target_variant)  # assumption
		if (length(target_variant)==0) ppred_variant="Zother" else
			if (length(target_variant)==1) ppred_variant=target_variant else {
				if (length(target_variant)==2) {
					core1=variant_corehap[[target_variant[1]]]$genes; core2=variant_corehap[[target_variant[2]]]$genes
					tmp=setdiff(core1, core2); core2=setdiff(core2, core1); core1=tmp
					tmp=strsplit2(x[1],","); tmp=sub("del","d", tmp); tmp=sub("stop","s", tmp)
					tmp1=nchar(tmp); tmp2=substr(tmp, 1, tmp1-1);
					tmpMut=substr(tmp, tmp1, tmp1);names(tmpMut)=tmp2
					tmp_core1=intersect(core1, tmp2);
					tmpIND1=sum(ref_AA[tmp_core1,"AA"]!=tmpMut[tmp_core1])
					tmp_core2=intersect(core2, tmp2);
					tmpIND2=sum(ref_AA[tmp_core2,"AA"]!=tmpMut[tmp_core2])
					if (tmpIND1==0 & tmpIND2==0) stop("check") else
						if (tmpIND1>0 & tmpIND2==0) ppred_variant=target_variant[1] else
							if (tmpIND1==0 & tmpIND2>0) ppred_variant=target_variant[2] else
								ppred_variant=paste(target_variant,collapse = "-")
				} else {if (length(target_variant)==3) {
					core1=variant_corehap[[target_variant[1]]]$genes;
					core2=variant_corehap[[target_variant[2]]]$genes
					core3=variant_corehap[[target_variant[3]]]$genes
					tmp1=setdiff(core1, c(core2,core3)); tmp2=setdiff(core2,c(core1,core3));
					tmp3=setdiff(core3,c(core1,core2));core1=tmp1;core2=tmp2;core3=tmp3;
					tmp=strsplit2(x[1],","); tmp=sub("del","d", tmp); tmp=sub("stop","s", tmp)
					tmp1=nchar(tmp); tmp2=substr(tmp, 1, tmp1-1);
					tmpMut=substr(tmp, tmp1, tmp1);names(tmpMut)=tmp2
					tmp=intersect(core1, tmp2);tmpIND1=sum(ref_AA[tmp,"AA"]!=tmpMut[tmp])
					tmp=intersect(core2, tmp2);tmpIND2=sum(ref_AA[tmp,"AA"]!=tmpMut[tmp])
					tmp=intersect(core3, tmp2);tmpIND3=sum(ref_AA[tmp,"AA"]!=tmpMut[tmp])
					tmp=as.numeric(c(tmpIND1, tmpIND2, tmpIND3)>0)
					if (sum(tmp==c(0,0,0))==3) stop("check") else
						if (sum(tmp==c(1,0,0))==3) ppred_variant=target_variant[1] else
							if (sum(tmp==c(0,1,0))==3) ppred_variant=target_variant[2] else
								if (sum(tmp==c(0,0,1))==3) ppred_variant=target_variant[3] else
									if (sum(tmp==c(1,1,0))==3) ppred_variant=paste(target_variant[1:2],collapse = "-") else
										if (sum(tmp==c(1,0,1))==3) ppred_variant=paste(target_variant[c(1,3)],collapse = "-") else
											if (sum(tmp==c(0,1,1))==3) ppred_variant=paste(target_variant[2:3],collapse = "-") else
												if (sum(tmp==c(1,1,1))==3) ppred_variant=paste(target_variant,collapse = "-")
				} else {print((target_variant))
					stop("Check for possible >four variant mixtures")}
				}
			}
			return(ppred_variant)}, ref_AA)
	pred_variant[,1]=ifelse(pred_variant[,1]=="Zother","UP", pred_variant[,1])
	pred_variant[,1]=ifelse(pred_variant[,1]=="Zmixture","MV", pred_variant[,1])
	ppred=ifelse(ppred=="Zother","UP", ppred)
	tmp=colnames(pred_prob1); tmp[tmp=="Zother"]="Others"; colnames(pred_prob1)=tmp
	
	# concordances
	variant_pred=ifelse(variant_pred=="Zother","UC",variant_pred)
	kappa_freq=table(pred_variant[,1], variant_pred)
	tmpc=colSums(kappa_freq); tmpr=rowSums(kappa_freq)
	kappa_freq=cbind(kappa_freq, total=tmpr); kappa_freq=rbind(kappa_freq, total=c(tmpc,sum(tmpc)))
	pkappa_freq=table(ppred, variant_pred)
	tmpc=colSums(pkappa_freq); tmpr=rowSums(pkappa_freq)
	pkappa_freq=cbind(pkappa_freq, total=tmpr); pkappa_freq=rbind(pkappa_freq, total=c(tmpc,sum(tmpc)))
	
	# output
	summary=data.frame(GISAID=variant_pred, AIS=pred_variant[,"pred"], post.AIS=ppred, prob= round(100*pred_variant[,"postp"],0))
	if (!is.null(pred_data$Lineage)) summary=data.frame(summary,  lineage=pred_data$Lineage)
	if (!is.null(pred_data$Pango.lineage)) summary=data.frame(summary,  lineage=pred_data$Pango.lineage)
	if (!is.null(pred_data$Clade)) summary=data.frame(summary,  clade=pred_data$Clade)
	if (!is.null(pred_data$Collection.date)) summary=data.frame(summary,  clade=pred_data$Collection.date)
	if (!is.null(pred_data$Location)) summary=data.frame(summary,  clade=pred_data$Location)

	result=NULL; result$summary=summary; 
	result$probability=data.frame(AIS=summary[,"AIS"], signif(pred_prob1,4)); result$inputdata=pred_data
	result$kappa=kappa_freq;result$kappa_pp=pkappa_freq
	return(result)
}


