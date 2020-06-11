
options <- commandArgs(trailingOnly=T)
if(length(options) < 1) stop("Invalid argument number\n\nRscript NMF_analysis.R -m ERP004264_PD.table \n\nor if you need more analysis option \n Rscript NMF_analysis.R -h\n\n")

library("argparser")

m_par <- arg_parser('NMF analysis R scripts', hide.opts = TRUE)
m_par <- add_argument(m_par, "--np", help="Setting the number of processors", default="p2")
m_par <- add_argument(m_par, "--nr", help="Setting the number of NMF calculation", default=50)
m_par <- add_argument(m_par, "--f", help="Setting the max number of feature", default=25)
m_par <- add_argument(m_par, "-m", help="input matrix")

argv <- parse_args(m_par)
max_f_n <- argv$f
nrun_set <- argv$nr
proc_n_set <- argv$np
data_file = argv$m

require(pheatmap)
require(RColorBrewer)
require(gplots)
require(preprocessCore)
require(NMF)
require(dendextend)
require(dplyr)
require(tidyr)
require(reshape2)
require(ggplot2)

# range of k values in NMF
k_range <- c(2:max_f_n)
adj_idx = k_range[1] - 1	# to start matrix index from 1 

# Clustering Purity calculation
ClusterPurity <- function(clusters, classes) {
	      sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

micro_abu <- read.table(data_file,stringsAsFactors=FALSE)	# row (micro-organism) x col (samples)
species_info <- micro_abu[,1]
micro_abu <- micro_abu[,-1]
disease_class <- micro_abu[1,]
micro_abu <- micro_abu[-1,]
micro_abu<-as.matrix(sapply(micro_abu, as.numeric))
rownames(micro_abu) <- species_info[-1]
colnames(micro_abu) <- disease_class

# check distribution...
micro_dist_per_sample <- c()
for(i in c(1:ncol(micro_abu))){
	micro_dist_per_sample <- c(micro_dist_per_sample, length( which(micro_abu[,i] > 0) ))  
}

write.table(micro_abu, file="micro_count_mat_pre.txt", quote=F, col.names=T,row.names=T)
quantile(micro_dist_per_sample, c(0.01,.05, .1, .15, .95, .99))

col_rem_idx <- c()
for(i in c(1:ncol(micro_abu)) ){	# for each sample
	if( length(which(micro_abu[,i] > 0)) < length(micro_abu[,i]) * 0.005){ #  below 0.5% 		
		col_rem_idx <- c(col_rem_idx,i)
	}
}

if(length(col_rem_idx)>1){ 
	micro_abu <- micro_abu[,-col_rem_idx]
}	

row_rem_idx <- c()
for(i in c(1:nrow(micro_abu)) ){	# for each microbiome
	if( length(which(micro_abu[i,] > 0)) < 2 ){	# at least 1 sequences and 1 samples	
		row_rem_idx <- c(row_rem_idx,i)
	}
}
if(length(row_rem_idx)!=0){ 
	micro_abu <- micro_abu[-row_rem_idx,]
}


#remove all row data sum is 0
all_0_row_id <- which(apply(micro_abu, 1, sum)==0)

if (length(all_0_row_id) != 0){
	all_0_row_id <- which(apply(micro_abu, 1, sum)==0)
	micro_abu <- micro_abu[-all_0_row_id,]
}


disease_class <- colnames(micro_abu)
species_name <- rownames(micro_abu)

write.table(species_info,file="microbe_list_original.txt",quote=F,row.names=F,col.names=F)
write.table(species_name,file="microbe_list_trimmed.txt",quote=F,row.names=F,col.names=F)
write.table(disease_class,file="sample_lable_list_trimmed.txt",quote=F,row.names=F,col.names=F)
write.table(micro_abu, file="micro_count_mat.txt", quote=F, col.names=T,row.names=T)

dim(micro_abu)

#quantile normalization
norm_micro_abu <- normalize.quantiles(micro_abu)

nmf_model_list <- c()
nmf_coef_list <- c()

purity_value_mat <- matrix(0,length(k_range),3)
pdf(file="coefmatrix.pdf")

for(i in k_range){
	res  <-  nmf(norm_micro_abu, i, "nsNMF",.opt=proc_n_set, nrun=nrun_set)    # single run, 18 body sites...2 institutes..
	nmf_model_list <- c(nmf_model_list,res)
	#raw_coefmat_gen
	coef_matrix <- scoef(res)
	colnames(coef_matrix) <- disease_class
	rownames(coef_matrix) <- paste("K",1:nrow(coef_matrix),sep="")
	write.table(coef_matrix,file=paste("raw_coefmat_",i,".txt", sep=""), quote=F, col.names=T,row.names=T)
	out <- pheatmap(coef_matrix,
	show_rownames=T,cluster_cols=T, cluster_rows=T, scale="row", clustering_distance_rows="euclidean",
	clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE, cex=0.7)
	res_coef <- coefmap(res,track=NA)
	nmf_coef_list <- c(nmf_coef_list,res_coef$Colv) 
	#purity calculation
	cut_idx <- cutree(res_coef$Colv,k=i)
	purity_disease <- ClusterPurity(cut_idx,unlist(disease_class))
	cat(paste("nmf",i,purity_disease,"\n"))
	purity_value_mat[i-adj_idx,1] <- i 
	purity_value_mat[i-adj_idx,2] <- purity_disease 
}
dev.off()

#basismap_coefmat_ generate
pdf(file="basismatrix.pdf")
for(K in k_range){
	#K=5  # need to decide K-value first
	basis_mat <- (basis(nmf_model_list[[K - adj_idx]]))   # K x 4471
	rownames(basis_mat) <- species_name
	colnames(basis_mat) <- paste("K",1:ncol(basis_mat),sep="")
	write.table(basis_mat,file=paste("raw_bmap_mat_",K+1-1,".txt", sep=""), quote=F, col.names=T,row.names=T)
	out <- pheatmap(basis_mat,
	show_rownames=T,cluster_cols=T, cluster_rows=T, scale="row", clustering_distance_rows="euclidean",
	clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE, cex=0.7)
}
dev.off()

#control for clusteringpurty Raw count hclust result
t_hc <- hclust(dist(t(norm_micro_abu)))
for(i in k_range){
	t_memb <- cutree(t_hc, k = i)
	purity_disease_c <- ClusterPurity(t_memb,unlist(disease_class))
	cat(paste("hclust",i,purity_disease_c,"\n"))
	purity_value_mat[i-adj_idx,3] <- purity_disease_c 
}
write.table(purity_value_mat,file="Clustering_purity_value.txt",quote=F,col.names=F,row.names=F)
# write to output file : causion, it will append the results into the same files again and again, delete all extracted...txt file first. 
for(K in k_range){
	for(i in (1:K)){
		nmf_m <- nmf_model_list[[K-1]]
        bmap <- basis(nmf_m)
        cat(paste(i,paste(rownames(micro_abu)[bmap[,i]>0],collapse=" ")),"\n",file=paste("extracted_microb_feat_by_",K,".txt",sep=''),append=TRUE)
    }
}

# extractFeatures and make kk_species_list.txt
col_r1 <- c(2:max_f_n)
for (k in col_r1){
	nmf_m <- nmf_model_list[[k-1]]
    col_r <- c(1:k)
    bmap <- basis(nmf_m)
    for (i in col_r){
		efx <- bmap[,i]>0
        sp_l <- rownames(micro_abu)[efx]
        dt_l <- bmap[efx , i ]
        write(paste(sp_l,":",dt_l,"','"), file="MMF_microbe_list.txt", append=TRUE)
        write(paste('@@',k,'-',i,'\n'), file="MMF_microbe_list.txt", append=TRUE)
    }
}

#connecting K - disease
for(K in k_range){
      coef_mat <- (scoef(nmf_model_list[[K - adj_idx]]))   # K x 4471
      colnames(coef_mat) <- disease_class
      rownames(coef_mat) <- paste("K",1:nrow(coef_mat),sep="")
      group_coefmat <- melt(coef_mat) %>% group_by(Var1, Var2) %>%  summarize(mean_contribution = mean(value))
      group_coefmat_unmelt <- as.data.frame(group_coefmat %>% spread("Var2","mean_contribution") %>% arrange(Var1))
      group_coefmat_unmelt_matrix <- group_coefmat_unmelt[,-1]
      rownames(group_coefmat_unmelt_matrix) <- group_coefmat_unmelt[,1]
      write.table(group_coefmat_unmelt_matrix,file=paste("group_coefmat_",K+1-1,".txt", sep=""), quote=F, col.names=T,row.names=T)
      tif_name <- paste("disease_k_connect_by_",K,".tiff",sep="")
 	  tiff(tif_name, width = 6, height = 6, units = 'in', res = 300)
	  aheatmap(group_coefmat_unmelt_matrix)
 	  dev.off()
}
