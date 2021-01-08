require(tidyverse)
require(cmapR)
require(factoextra)
require(FactoMineR)
require(ICC)
require(qvalue)
require(RColorBrewer)
require(gridExtra)
require(ggrepel)

###1. GCT Parsing, Normalization, and Filtering function ------------------------------------

sm_to_gct_lfq <- function(x,md,col.suffix="totalIntensity",join.by="sample"){
	#Function to convert a spectrum mill LFQ ratio report into a gct file
	#Args:
	#	x: a spectrum mill lfq report stored as data frame 
	#	species: species to be selected from the table. If NA, all species are allowed
	#	col.suffix: the suffix of the columns that contain the data. Spetrum Mill
	#	provides total intensities (defualt), peptide counts, and spectra counts.
	#	join.by: Specifies the column that is common between the 
	#Returns:
	#	A gct object

	mat <- x %>% select(ends_with(col.suffix)) %>% as.matrix
	cn <- data.frame(col = colnames(mat)) %>% separate("col",c("sample",NA)," ") 
	
	rownames(mat) <- x$id
	rdesc <- select(x,-ends_with(col.suffix)) 
	cdesc <- full_join(md,cn,by=join.by) %>% rename(id = sample)

	colnames(mat) <- cdesc$id
	if(ncol(mat)!=nrow(cdesc)|nrow(mat)!=nrow(rdesc)){
		stop("Column or row annotations do not match 
		     the matrix size.")
	}
	
	return(new("GCT", mat=mat, rdesc=rdesc, cdesc=cdesc))
}


sm_to_gct <- function(x,md,cols=NA,join.by="tmtlabel"){
	#Function to convert a spectrum mill tmt ratio report into a gct file
	#Args:
	#	x: a spectrum mill tmt ratio report stored as data frame 
	#	species: species to be selected from the table. If NA, all species are allowed
	#	cols: the columns that contain the data (as names or as column indices). If NA,
	#	then use the columns that contain ":", which are the ratio columns
	#	join.by: Specifies the column that is common between the 
	#Returns:
	#	A gct object
	if(is.na(cols[1])){
		mat <- select(x,contains(":")) %>% as.matrix
		cn <- parse_sm_colnames_ratios(colnames(mat))
	} else {
		mat <- x[,cols] %>% as.matrix
		cn <- parse_sm_colnames(colnames(mat))
	}
	
	rownames(mat) <- x$id
	rdesc <- select(x,-contains(",")) 
	cdesc <- right_join(md,cn,by=join.by) %>% rename(id = sample)

	#colnames(cdesc)[which(colnames(cdesc)=="sample")] <- "id"
	colnames(mat) <- cdesc$id
	if(ncol(mat)!=nrow(cdesc)|nrow(mat)!=nrow(rdesc)){
		stop("Column or row annotations do not match 
		     the matrix size.")
	}
	return(new("GCT", mat=mat, rdesc=rdesc, cdesc=cdesc))
}


parse_sm_colnames <- function(x){
	#Function to parse the column names for TMT intensity datasets obtained from the process report of spectrum mill
	#Args:
	#	x: the column names for tmt ratios obtained by spectrum mill process report
	#Returns:
	#	A data frame with the plex, tmt label used, and sample name as columns
	
	data.frame(col = x) %>% separate("col",c("plex","tmtlabel","sample"),", ") 
	
}

parse_sm_colnames_ratios <- function(x){
	#Function to parse the column names for TMT ratio datasets obtained from the process report of spectrum mill
	#Args:
	#	x: the column names for tmt ratios obtained by spectrum mill process report
	#Returns:
	#	A data frame with the plex, tmt label used, and sample name as columns
	data.frame(col = x) %>% separate("col",c("plex","tmtlabel","sample"),",") %>% 
		transmute(plex = plex,
			  tmtlabel = str_match(tmtlabel," (.*):")[,2],
			  sample = str_match(sample," (.*):")[,2])

}

gct_ratio <- function(x,denom){
	#Convert data to log2 ratios by dividing to the median of samples selected as denominators
	#Args:
	#	x: a GCT object
	#	denom: Character vector indicating the sample ids to be used as denominators
	
	if(length(denom)>1){
		d <- apply(x@mat[,denom],1,mean,na.rm=T)
	} else {
		d <- x@mat[,denom]
	}
	
	
	#x <- subset_gct(x,cid = setdiff(colnames(x@mat),denom)) 
	x@mat <- log2(x@mat/d)
	return(x)
}

median_mad_norm <- function(x,mad=T){
	#Perform median normalization
	#Args:
	#	x: a GCT object
	#	mad: logic indicating to use mad normalization
	#Returns:
	#	A normalized GCT object
	
	
	if(mad){
		x@mat <- scale(x@mat,center = apply(x@mat,2,median,na.rm=T),
			       scale = apply(x@mat,2,mad,na.rm=T))
		x@mat <- x@mat * (attr(x@mat,"scaled:scale") %>% mean)
	} else{
		x@mat <- scale(x@mat,center = apply(x@mat,2,median,na.rm=T),
			       scale = F)

	}
	return(x)
}

remove_all_na <- function(x){
	#Remove proteins that have all missing values
	#Args:
	#	x: a GCT object
	#Returns:
	#	A GCT object without proteins that have all missing values
	x <- subset_gct(x,rid = which(rowSums(is.na(x@mat)) < ncol(x@mat))) 
	return(x)
}

remove_na <- function(x,pct){
	#Remove proteins that have more than the set percentage of missing values
	#Args:
	#	x: a GCT object
	#	pct: the percent cutoff
	#Returns:
	#	A GCT object without proteins that have any missing values
	
	x <- subset_gct(x,rid = which(rowSums(is.na(x@mat)) <= pct*ncol(x@mat))) 
	return(x)
}

filter_gct <- function(x,column, value){
	#Filter GCT object based on row metadata
	#Args:
	#	x: a GCT object
	#	column: name of metadata column to use for filtering
	#	value: The value that should be matched for the filtering
	#Returns:
	#	A GCT object with species filtered
	x <- subset_gct(x,rid = which(unlist(x@rdesc[,column]) %in% value)) 
	return(x)
}

impute_quantile <- function(x,quant,cols=NA){
	#Impute values in the matrix by taking the set quantile. 
	#Imputation is performed per sample (column)
	#Args:
	#	x: a GCT object
	#	quant: the quantile used for imputation
	#	cols: the columns to be imputed or NA for all
	list.imputed <- rowSums(is.na(x@mat))>0
	
	impute_q <- function(y,quant){
		return.y <- y
		return.y[is.na(y)] <- quantile(y,quant,na.rm=T)
		return(return.y)
	}
	
	if(is.na(cols)){
		x@mat <- apply(x@mat,2,impute_q,quant=quant)
	} else {
		x@mat[,cols] <- apply(x@mat[,cols],2,impute_q,quant=quant)
	}
	
	x@rdesc <- data.frame(x@rdesc,has.imputed = list.imputed)
	
	return(x)
}

normalize_0to1 <- function(x){
  #Normalize values in a numeric vector from 0 to 1
  return((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm = T)))
}


###2. Format and manipulate data-------------------------------------------

make_results <- function(x,contrast.matrix,fit.eb,one.sided=F,lower=F,qval =T){
	#Calculate q values and generates a result table
	#Args:
	#	x: a GCT object
	#	contrast.matrix: contrast matrix used to calculate contrasts
	#	fit.eb: The fitted linear model with empirical bayes calculation for adj. pvalue
	#Returns:
	#	A dataframe with the original expression values, pvalues, and qvalues
	res <- as.data.frame(x@mat) %>% rownames_to_column(var="id")
	
	for(i in colnames(contrast.matrix)){
		
		stats <- topTable(fit.eb, number = nrow(x@mat), sort.by = "none", coef = i) %>%
			select(logFC, P.Value, adj.P.Val) 
		
		if(one.sided){
			stats$P.Value <- limma.one.sided(fit.eb)[,i]
		}
		
		if(qval){
			stats <- stats %>%
				mutate(Q.Value = qvalue(P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues) %>%
				mutate_at(vars(), signif, 4)
		} else {
			stats <- stats %>%
				mutate(Q.Value = adj.P.Val) %>%
				mutate_at(vars(), signif, 4)
		}
		
		colnames(stats) <- paste0(i,".",colnames(stats))
		res <- cbind(res,stats)
		print(i)
		#summary(qvalue(fit.eb$p.value[,i]))
		
	}
	
	return(res)
}

make_nested_results <- function(x, contrast.matrix, fit.eb,one.sided=F,lower=F,qval=T){
	#Calculate q values and generates a result table in nested format
	#Args:
	#	x: a GCT object
	#	contrast.matrix: contrast matrix used to calculate contrasts
	#	fit.eb: The fitted linear model with empirical bayes calculation for adj. pvalue
	#Returns:
	#	A nested dataframe for each contrast with the original expression values, 
	#       pvalues, and qvalues
	
	res <- as.data.frame(x@mat) %>% rownames_to_column(var="id")
	
	#Holder for list of data to be place in the nested dataframe
	data <- list()
	
	for(i in colnames(contrast.matrix)){
		
		stats <- topTable(fit.eb, number = nrow(x@mat), sort.by = "none", coef = i) %>%
			select(logFC, P.Value, adj.P.Val)
		
		if(one.sided){
			stats$P.Value <- limma.one.sided(fit.eb)[,i]
		}
		
		if(qval){
			stats <- stats %>% mutate(Q.Value = qvalue(P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues) %>%
				mutate_at(vars(), signif, 4)
		} else {
			stats <- stats %>% mutate(Q.Value = adj.P.Val)	
		}
		
		data <- c(data,list(cbind(res,stats)))
		
		print(i)
		#summary(qvalue(fit.eb$p.value[,i]))
	}
	return(tibble(contrasts = colnames(contrast.matrix),
		      data = data))
}

make_results_ftest <- function(x,fit.eb,coef=NULL){
	#Calculate q values and generates a result table
	#Args:
	#	x: a GCT object
	#	coef: lsit of coefficients to perform the F-test (or empty for all)
	#	fit.eb: The fitted linear model with empirical bayes calculation for adj. pvalue
	#Returns:
	#	A dataframe with the original expression values, pvalues, and qvalues
	res <- as.data.frame(x@mat) %>% rownames_to_column(var="id")
	
	stats <- topTable(fit.eb, number = nrow(x), sort.by = "none", coef = coef) %>%
			mutate(Q.Value = qvalue(P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues) %>%
		mutate_at(vars(), signif, 4)
	res <- cbind(res,stats)
	#summary(qvalue(stats$P.Value))
	
	
	
	return(res)
}

map.accessions <- function(x,table){
	# Separates a string accessions delimited by pipes and then tries to find a mapped value 
	# in a mapping table. Column 1 is assumed to be From and column 2 to be To
	# Args:
	#   x: A character string of accession numbers separated by a pipe
	#   table: The list of accessions to match to
	# Returns:
	#   True if a match is found
	x.split <- unlist(strsplit(as.character(x),split="\\|"))
	x.split <- sub("-.*","",x.split)
	
	map.single <- function(m){
		to.return <- table %>% filter_at(1,all_vars(. == m))
		if(nrow(to.return) > 0){
			return(as.character(to.return[,2]))
		} else {
			return("N/A")
		}
	}
	
	return(paste0(unlist(map(x.split,map.single)),collapse = "|"))
}

#' Internal algorithm: Make limma test one-sided
#' @param fit Result of "lmFit" and "eBayes" functions in "limma" package.
#' @param lower Shall one-sided p-value indicated down-regultation?
limma.one.sided = function (fit,lower=FALSE) {
	se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
	df.total <- fit$df.prior + fit$df.residual
	pt(fit$t, df=df.total, lower.tail=lower)
}

get_significant_summary <- function(results,qtreshold=0.05,logFC = 0){
	#Creates a table that summarizes the number of significantly regulated genes
	#Args:
	#	results: A nested result table as produced by make_nested_results
	#	qtreshold: q-value treshold to count a significant differential abundance
	output <- results %>% mutate(
		total = map_dbl(data,function(x){sum(!is.na(x$Q.Value<qtreshold))}),
		significant = map_dbl(data,function(x){sum(x$Q.Value<qtreshold & abs(x$logFC)>logFC,na.rm=T)}),
		up = map_dbl(data,function(x){sum(x$logFC>logFC & x$Q.Value<qtreshold,na.rm=T)}),
		down = map_dbl(data,function(x){sum(x$logFC< (-1*logFC) & x$Q.Value<qtreshold,na.rm=T)})
	) %>% select(-data)
	return(output)
}

###3. Produce plots-----------------------------------------------------------------

pca_plots <- function(pca,cdesc,axes = c(1,2),title="PCA",ellipses = T,...){
	#Generate pca plots colored with specific metadata
	#Args:
	#	pca: a pca object as generated by PCA
	#	cdesc: a data.frame indicating the features to color
	#	axes: the pca dimensions to plot
	#	title: the desired title for the plot
	#	...: Other arguments passed to fviz_pca
	#Returns:
	#	A list of pca plots
	output <- list()
	for(i in colnames(cdesc)){
		g <- fviz_pca_ind(pca,col.ind = as.factor(cdesc[,i]),
				  mean.point=FALSE,geom="point",axes=axes,
				  addEllipses = ellipses,...)+
			labs(title=title,subtitle = paste0("Labeled by ",i))+
			theme_bw()
		output <- c(output,list(g))
	}
	names(output) <- colnames(cdesc)
	return(output)
}

pca_correlations <- function(pca,cdesc,components = c(1:5)){
	#Calculates the intragroup correlation coefficient between each principal
	#component and the sample annotation
	#Args:
	#	pca: a pca object as generated by PCA
	#	cdesc: a data.frame indicating the features to color
	#	components: the components to be correlated. defaults to the first 10
	#
	
	#Obtain the principal component coordinates
	p <- data.frame(pca$ind$coord[,components])
	data <- data.frame(dims=character(),
			   correlation=numeric(),
			   experimental.factor = character() )
	
	for(i in colnames(cdesc)){
		#Create the class column to be used for correlation with the PC coordinates
		class = data.frame(matrix(rep(as.factor(cdesc[,i]),ncol(p)),ncol = ncol(p)))
		
		#Calculate intragroup correlations
		y <- map2(class,p,function(class,p){
			return(ICCest(class,p)$ICC)
		}) %>% unlist
		
		dimensions <- as.numeric(gsub("Dim\\.","",colnames(p)))
		
		data <- rbind(data,
			      data.frame(dims = dimensions,
			      	   correlation = y,
			      	   experimental.factor=rep(i,ncol(p))))
		
	}
	ggplot(data=data,aes(as.factor(dims),correlation,group=experimental.factor,color=experimental.factor)) +
		geom_line() + geom_point() + labs(x="Component (% Variance Explained)") + theme_bw()+
		scale_x_discrete(labels = paste0("PC",data$dims, " (",round(pca$eig[,2],1),"%)"))+
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_protein <- function(x,id,cdesc.group="treatment",cdesc.fill=NA){
	#Plot the values of a gct object for a specific proteins
	#Args:
	#	x: a GCT object
	#	id: id of the protein to be plotted
	#	cdesc.value: The cdesc column to be used for the grouping
	#	cdesc.color: The fill color to use for additioanl separation
	#	             If NA, no separation occurs. 
	data <- melt.gct(x) %>% filter(id.x==id) 
	g<- ggplot(data, aes_string(x=cdesc.group,y="value"))+labs(title=id)+theme_bw()
	if(!is.na(cdesc.fill)){
		g<- g+geom_boxplot(aes_string(fill=cdesc.fill),outlier.shape=NA)+
			geom_point(position=position_dodge(width=0.75),aes_string(group=cdesc.fill))	
	} else{
		g <- g<- g+geom_boxplot(outlier.shape=NA)+
			geom_point(position=position_dodge(width=0.75))
	}
	return(g)
}

plot_volcano <- function(contrasts,results,qtreshold=0.05,lab.sig=F,qval=T){
	#Generates a volcano plot using p-values and logFC ratios
	#Args:
	#	contrasts: Name of the contrast being plotted
	#	results: a results table with logFC, q values, and p values
	#	qtreshold: The q-value treshold to be used. Default is 0.05.
	#	lab.sig: Whether to add labels to the significant data points	
  if(!qval){
    results <- results %>% mutate(Q.Value = adj.P.Val)
  }
	results <- results[!is.na(results$Q.Value),]
	
	g <- ggplot(results,aes(logFC,-log10(P.Value),color = Q.Value < qtreshold)) + 
		geom_point(size=2,alpha=0.7) + 
		scale_color_manual(values= c("Grey",brewer.pal(n=3,name = "Set2")[1]),
				   labels=c("Non-significant","Significant"),
				   name = paste0("q value < ",qtreshold))+
		labs(x= paste0("Log2(",regmatches(contrasts,regexpr("^[[:alnum:]]+",contrasts))," over ",regmatches(contrasts,regexpr("[[:alnum:]]+$",contrasts)),")"), 
		     y="-Log10(p-value)",
		     title = paste0("Comparison ",contrasts),
		     subtitle = paste0(sum(results$Q.Value<qtreshold)," significant out of ",nrow(results),
		     		  ". ",sum(results$logFC>0&results$Q.Value<qtreshold), " increasing and ",
		     		  sum(results$logFC<0&results$Q.Value<qtreshold)," decreasing"))+
		theme_bw()
	
	
	if(min(results$Q.Value) < qtreshold){
		p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% 
			    	arrange(min))[[1,"P.Value"]]
		
		g<- g+	annotate("text",x=max(results$logFC),y=-log10(p.tresh),
				label=paste0("q value < ",qtreshold),vjust = -1,hjust="inward")+
			geom_hline(yintercept = -log10(p.tresh))
	}
	
	if(lab.sig){
		labels <- results$id
		labels[results$Q.Value>=qtreshold] <- ""
		g <- g+geom_text_repel(aes(label=labels),color="gray25")
	}
	
	
	return(g)
}

plot_volcano_color <- function(contrasts,results,color.list=NA,qtreshold=0.05,color.name="",
			       colors=c("Grey",brewer.pal(n=3,name = "Set2")[2]),
			       xlims = NA, ylims = NA,labels=character(),alpha.list = 0.7){
	#Generates a volcano plot using p-values and logFC ratios
	#Args:
	#	contrasts: Name of the contrast being plotted
	#	results: a results table with logFC, q values, and p values
	#	color.list: A vector with color categories to be plotted
	#	qtreshold: The q-value treshold to be used. Default is 0.05.
	to.keep <- !is.na(results$Q.Value)
	results <- results[to.keep,]
	p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% arrange(min))[[1,"P.Value"]]
	
	g <- ggplot(results,aes(logFC,-log10(P.Value),color = id %in% color.list)) + 
		geom_point(size=2,alpha=0.3) + 
		scale_color_manual(values=colors,labels=c("False","True"), name = color.name)+
		geom_hline(yintercept = -log10(p.tresh))+
		labs(x= paste0("Log2(",regmatches(contrasts,regexpr("^[[:alnum:]]+",contrasts)),
			       " over ",regmatches(contrasts
			       		    ,regexpr("[[:alnum:]]+$",contrasts)),")"), 
		     y="-Log10(p-value)",
		     title = paste0("Comparison ",contrasts),
		     subtitle = paste0(sum(results$Q.Value<qtreshold)," significant proteins out of ",nrow(results)))+
		theme_bw()+
		annotate("text",x=max(results$logFC),y=-log10(p.tresh),
			 label=paste0("q value < ",qtreshold),vjust = -1,hjust="inward")	
	
	
	if(!is.na(xlims[1])){
		g <- g+scale_x_continuous(limits=xlims)
	}
	
	if(!is.na(ylims[1])){
		g <- g+scale_y_continuous(limits=ylims)
	}
	
	if(length(labels)>0){
		labels <- labels[to.keep]
		g<-g+geom_text_repel(aes(label= as.character(labels)),nudge_y = 1)
	}
	
	
	return(g)
}

plot_volcano_color2 <- function(contrasts,results,color.list="Gray",qtreshold=0.05,color.name="",
				colors=c("Grey",brewer.pal(n=3,name = "Set2")[2]),
				labels=character(), label.size = 5,
				alpha.list = 0.7,
				shape.list="", shape.name="",shapes=c(16,4,15,17),
				xlims = NA, ylims = NA){
	#Generates a volcano plot using p-values and logFC ratios
	#Args:
	#	contrasts: Name of the contrast being plotted
	#	results: a results table with logFC, q values, and p values
	#	color.list: A vector with color categories to be plotted
	#	qtreshold: The q-value treshold to be used. Default is 0.05.
	#	labels: A character vector with equal rows as results 
	#results <- results[!is.na(results$adj.P.Val),]
	p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% arrange(min))[[1,"P.Value"]]
	
	g <- ggplot(results,aes(logFC,-log10(P.Value),color = color.list)) + 
		geom_point(size=2,aes(shape=shape.list,alpha = alpha.list)) + 
		scale_shape_manual(values=shapes,name=shape.name)+
		scale_color_manual(values=colors, name = color.name)+
		geom_hline(yintercept = -log10(p.tresh))+
		labs(x= paste0("Log2(",regmatches(contrasts,regexpr("^[[:alnum:]]+",contrasts)),
			       " over ",regmatches(contrasts
			       		    ,regexpr("[[:alnum:]]+$",contrasts)),")"), 
		     y="-Log10(p-value)",
		     title = paste0("Comparison ",contrasts),
		     subtitle = paste0(sum(results$Q.Value<qtreshold,na.rm=T)," significant proteins out of ",nrow(results)))+
		theme_bw()+
		guides(alpha=F)+
		annotate("text",x=max(results$logFC),y=-log10(p.tresh),
			 label=paste0("q value < ",qtreshold),vjust = -1,hjust="inward")
	if(length(labels)>0){
		g<-g+geom_text_repel(aes(label= as.character(labels)),size=label.size)
	}
	if(!is.na(xlims[1])){
		g <- g+scale_x_continuous(limits=xlims)
	}
	
	if(!is.na(ylims[1])){
		g <- g+scale_y_continuous(limits=ylims)
	}
	return(g)
}

plot_time_clust <- function(results,cont, clust, k){
	#Plot temporal time traces clustered according to an hclust object
	# and a set number of clusters (k)
	#Args:
	#	results: A results data frame
	#	cont: The column names for the contrasts to be plotted
	#	clust: A vector indicating clusters (as produced by cuttree or pheatmap)
	#	k: THe number of clusters to make
	
	results <- cbind(results,cluster = clust)
	
	
	line_plot <- function(cluster,data,cont,y.limits){
		data.l <- gather(data,"timepoint","value",cont)
		data.l$timepoint <- factor(data.l$timepoint,levels = cont)
		g <- ggplot(data.l,aes(x=timepoint,y=value,group=id))+geom_line(alpha=0.5)+
			theme_bw() +labs(title=paste0("Cluster ",cluster),
					 x="Time point",
					 y="log2(Ratio to control)",
					 subtitle= paste0(nrow(data)," proteins or peptides"))+
			scale_y_continuous(limits=y.limits)+
			geom_hline(aes(yintercept=0),color="red")
		
		return(g)
	}
	
	y.limits <- c( results %>% select(cont) %>% unlist() %>% min() ,
		       results %>% select(cont) %>% unlist() %>% max() )
	
	# y.limits <- c( results %>% select(cont) %>% unlist() %>% quantile(0.01) ,
	# 	       results %>% select(cont) %>% unlist() %>% quantile(0.99) )
	
	print(y.limits)
	gl <- results %>% group_by(cluster) %>% nest %>%
		pmap(line_plot,cont,y.limits)
	
	
	#grid.arrange(grobs=gl,ncol=2)
	return(gl)
}


plot_correlations <- function(data,col.x,col.y,color.list="Gray",color.name="",
			      colors=c("Grey",brewer.pal(n=3,name = "Set2")[2]),
			      labels=character(),alpha.list = 0.7,shape.list=16,
			      shapes=c(4,15,16,17)){
	#Plot correlations between two numerical vectors with options to add colors or labels
	#Args:
	#	data:	Data frame containing columns to be plotted
	#	col.x,col.y: Numerical vectors
	#	color.list: A vector with color categories to be plotted
	#	qtreshold: The q-value treshold to be used. Default is 0.05.
	#	labels: A character vector with equal rows as results 
	
	g<- ggplot(data, aes_string(x=col.x,y=col.y))+
		geom_point(alpha=alpha.list,aes(color = color.list,shape=shape.list))+
		scale_color_manual(values=colors, name = color.name)+
		theme_bw()
	
	if(length(shape.list>1)){
		g <- g+scale_shape_manual(values=shapes)
	}
	
	if(length(labels)>0){
		g<-g+geom_text_repel(aes(label= as.character(labels)))
	}
	return(g)
}
