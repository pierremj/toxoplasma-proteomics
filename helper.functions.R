require(tidyverse)
require(cmapR)
require(factoextra)
require(FactoMineR)
require(ICC)
require(qvalue)
require(RColorBrewer)
require(gridExtra)

###1. GCT Parsing, Normalization, and Filtering function ------------------------------------

sm_to_gct <- function(x,md){
	#Function to convert a spectrum mill tmt ratio report into a gct file
	#Args:
	#	x: a spectrum mill tmt ratio report stored as data frame 
	#	species: species to be selected from the table. If NA, all species are allowed
	#Returns:
	#	A gct object
	mat <- select(x,contains(":")) %>% as.matrix
	rownames(mat) <- x$id
	cn <- parse_sm_colnames(colnames(mat))
	colnames(mat) <- cn$sample
	rdesc <- select(x,-contains(" ")) 
	cdesc <- full_join(md,cn,by="sample") %>% rename(id = sample)
	if(ncol(mat)!=nrow(cdesc)|nrow(mat)!=nrow(rdesc)){
		stop("Column or row annotations do not match 
		     the matrix size.")
	}
	return(new("GCT", mat=mat, rdesc=rdesc, cdesc=cdesc))
}


parse_sm_colnames <- function(x){
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
	x <- subset.gct(x,rid = which(rowSums(is.na(x@mat)) < ncol(x@mat))) 
	return(x)
}

remove_na <- function(x,pct){
	#Remove proteins that have more than the set percentage of missing values
	#Args:
	#	x: a GCT object
	#	pct: the percent cutoff
	#Returns:
	#	A GCT object without proteins that have any missing values

	x <- subset.gct(x,rid = which(rowSums(is.na(x@mat)) <= pct*ncol(x))) 
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
	x <- subset.gct(x,rid = which(x@rdesc[,column] == value)) 
	return(x)
}


###2. Format and manipulate data-------------------------------------------

make_results <- function(x,contrast.matrix,fit.eb){
	#Calculate q values and generates a result table
	#Args:
	#	x: a GCT object
	#	contrast.matrix: contrast matrix used to calculate contrasts
	#	fit.eb: The fitted linear model with empirical bayes calculation for adj. pvalue
	#Returns:
	#	A dataframe with the original expression values, pvalues, and qvalues
	res <- as.data.frame(x@mat) %>% rownames_to_column(var="id")
	
	for(i in colnames(contrast.matrix)){
		
		stats <- topTable(fit.eb, number = nrow(x), sort.by = "none", coef = i) %>%
			select(logFC, P.Value, adj.P.Val) %>%
			mutate(Q.Value = qvalue(P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues) %>%
			mutate_at(vars(), signif, 4)
		colnames(stats) <- paste0(i,".",colnames(stats))
		res <- cbind(res,stats)
		print(i)
		summary(qvalue(fit.eb$p.value[,i]))
		
	}
	
	return(res)
}

make_nested_results <- function(x, contrast.matrix, fit.eb){
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
		
		stats <- topTable(fit.eb, number = nrow(x), sort.by = "none", coef = i) %>%
			select(logFC, P.Value, adj.P.Val) %>%
			mutate(Q.Value = qvalue(P.Value, fdr.level=0.05, pi0.method="bootstrap")$qvalues) %>%
			mutate_at(vars(), signif, 4)
		data <- c(data,list(cbind(res,stats)))

		print(i)
		summary(qvalue(fit.eb$p.value[,i]))
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
	summary(qvalue(stats$P.Value))
	
	
	
	return(res)
}

###3. Produce plots-----------------------------------------------------------------

pca_plots <- function(pca,cdesc,axes = c(1,2),title="PCA",...){
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
			     addEllipses = T,...)+
			labs(title=title,subtitle = paste0("Labeled by ",i))+
			theme_bw()
		output <- c(output,list(g))
	}
	names(output) <- colnames(cdesc)
	return(output)
}

pca_correlations <- function(pca,cdesc,components = c(1:20)){
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

plot_volcano <- function(contrasts,results,qtreshold=0.05){
	#Generates a volcano plot using p-values and logFC ratios
	#Args:
	#	contrasts: Name of the contrast being plotted
	#	results: a results table with logFC, q values, and p values
	#	qtreshold: The q-value treshold to be used. Default is 0.05.
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
	
	
	return(g)
}

plot_volcano_color <- function(contrasts,results,color.list=NA,qtreshold=0.05,color.name="",colors=c("Grey",brewer.pal(n=3,name = "Set2")[2])){
	#Generates a volcano plot using p-values and logFC ratios
	#Args:
	#	contrasts: Name of the contrast being plotted
	#	results: a results table with logFC, q values, and p values
	#	color.list: A vector with color categories to be plotted
	#	qtreshold: The q-value treshold to be used. Default is 0.05.
	results <- results[!is.na(results$Q.Value),]
	p.tresh <- (results %>% mutate(min = abs(qtreshold-results$Q.Value)) %>% arrange(min))[[1,"P.Value"]]
	
	g <- ggplot(results,aes(logFC,-log10(P.Value),color = id %in% color.list)) + 
		geom_point(size=2,alpha=0.7) + 
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
