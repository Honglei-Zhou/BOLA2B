library(gplots)
library(Rtsne)
library(png)
library(sva)
library(plyr)
library(cowplot)
library(edgeR)
library(limma)
library(Glimma)
library(ggplot2)
library(RColorBrewer)
library(ggfortify)
library(RUVcorr)
library(reshape)
library(data.table)
library(EDASeq)
library(biomaRt)
library(ggvenn)
library(multiMiR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(idmap3)
library(oligo)
library(survival)
library(survminer)
library(dplyr)
library(stringr)
library(htmlTable)
library(GGally)
library(network)
library(purrr)
library(tidyverse)
library(rsample)
library(pROC)
library(kableExtra)

# survival analysis
# survival.df : contains survival data
# my_df: contains expression data of BOLA2B in different cancers

covariates <- c('BOLA2B')

result <- NULL

for(sample in tcga_ids){

	my_df.sub <- my_df[which(my_df$TCGA_ID == sample),]

	my_means <- colMedians(as.matrix(my_df.sub[,covariates]))
	for(i in c(1:length(covariates))){
	  v <- my_means[i]
	  if(v == 0){
	    high <- my_df.sub[,covariates[i]] > v
	  }else{
	    high <- my_df.sub[,covariates[i]] >= v
	  }
	  my_df.sub[,covariates[i]][which(high==TRUE)] <- 'high'
	  my_df.sub[,covariates[i]][which(high==FALSE)] <- 'low'
	}
	  

	univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS_time, OS)~', x)))
	model_names <- names(univ_formulas)
	# print(model_names)
	i <- 1
	univ_models <- list()
	for(m in univ_formulas){
	  skip_to_next <- FALSE
	  p <- tryCatch({
	      print(m)    
	      print(model_names[i])
	      fit <- coxph(m, data = my_df.sub)
	      univ_models[[model_names[i]]] <- fit
	      i <- i+1 
	      
	    },error = function(e){
	      i <- i+1
	      print(sample)
	      print('\n')
	      skip_to_next <- TRUE
	      }, 
	    warning = function(cond) {
	      message(cond)
	      print(sample)
	      print('\n')
	      i <- i+1
	      skip_to_next <- TRUE
	    })
	  # print(p)
	  if(skip_to_next == TRUE){next}

	} 
	j <- 1
	km_models <- list()
	for(m in univ_formulas){
	  skip_to_next <- FALSE
	  p <- tryCatch({
	      print(m)    
	      print(model_names[j])
	      fit <- surv_fit(m, data = my_df.sub)
	      km_models[[model_names[j]]] <- fit
	      j <- j+1 
	      
	    },error = function(e){
	      i <- i+1
	      print(sample)
	      print('\n')
	      skip_to_next <- TRUE
	      }, 
	    warning = function(cond) {
	      message(cond)
	      print(sample)
	      print('\n')
	      j <- j+1
	      skip_to_next <- TRUE
	    })
	  # print(p)
	  if(skip_to_next == TRUE){next}
	} 

	print('Start to analyze coxph')
	univ_results <- lapply(univ_models,
	                       function(x){ 
	                        fit.coxph <- x
	                        x <- summary(x)
	                        print(x)
	                        p.value<-signif(x$coefficients[,'Pr(>|z|)'], digits=2)
	                        wald.test<-signif(x$wald["test"], digits=2)
	                        beta<-signif(x$coef[1], digits=2);#coeficient beta
	                        HR <-signif(x$coef[2], digits=2);#exp(beta)
	                        total <- x$n
	                        nevent <- x$nevent
	                        HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	                        HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	                        HR.confint <- paste0(HR, " (", 
	                                   HR.confint.lower, "-", HR.confint.upper, ")")
	                        file_name <- rownames(x$conf.int)[1]
	                        file_name <- str_split(file_name, 'l')[[1]][1]
	                      
	                        expr_gene <- table(my_df.sub[,file_name])
	                        expr_gene <- as.data.frame(expr_gene)
	                        tmp<-c(file_name, expr_gene[1,2],expr_gene[2,2],total, nevent, beta, HR, HR.confint.lower, HR.confint.upper, HR.confint, wald.test, p.value, sample)
	                   
	                        names(tmp)<-c("Symbol",levels(expr_gene[1,1])[expr_gene[1,1]],levels(expr_gene[2,1])[expr_gene[2,1]],"Total", "Event", "beta", "HR", "HR_low", "HR_high","HR (95% CI for HR)", "wald.test", 
	                                                "p.value", "TCGA_ID")
	                        return(tmp)
	                         })
	univ_km <- lapply(km_models, function(x){
	  my_fit <- x
	  # print(fit)
	  x <- summary(x)
	  print(x$table)
	  file_name <- rownames(x$table)[1]
	  file_name <- str_split(file_name, '=')[[1]][1]
	  print(file_name)
	  p.value <- surv_pvalue(my_fit)$pval
	  # p.value<-signif(x$coefficients[,'Pr(>|z|)'], digits=2)
	  print(p.value)
	  if(!is.na(p.value) && p.value < 0.05){
	  	# save KM plot
	    png(paste('result/TCGA.', sample,'.', file_name,'.KM.png', sep=""), height=2000, width=2000, res=300)
	    print(my_fit)
	    sample_s <- str_split(sample,'_')[[1]][[2]]
	    p <- ggsurvplot(
	          my_fit,                     
	          data = my_df.sub,            
	          risk.table = TRUE,       
	          risk.table.col = file_name, 
	          pval = TRUE,         
	          conf.int = TRUE,     
	          palette = "Dark2",
	          xlab = "Time in days", 
	          ggtheme = theme_bw(), 

	          conf.int.style = "step",
	          legend.title = file_name,
	          legend.labs = c("High", "Low"))
	    print(p)
	    # ggsave(paste('result/TCGA.', sample,'.', file_name,'.KM.png', sep=""),height=3,width=3)
	    dev.off()
	  }
	  return(file_name)
	  })
	tmp <- t(as.data.frame(univ_results, check.names = FALSE))
	print(tmp)
	tmp <- as.data.frame(tmp)

	if(dim(tmp)[1] == 0) next

	if (is.null(result)){
	  result <- tmp
	}else{
	  result <- rbind(result, tmp)
	}

	rownames(result) <- c(1:dim(result)[1])
}

# save HR data
write.csv(result, 'result/TCGA.TCGA_PAAD.mRNA.HR.csv')

# Make HR plot
result <- result[which(result$p.value < 0.05),]

for(sample in unique(result$TCGA_ID)){
  my_plot_data <- result[which(result$TCGA_ID == sample),]

  my_plot_data$Index <- as.numeric(c(1:dim(my_plot_data)[1]))
  my_plot_data$HR <- as.numeric(my_plot_data$HR)
  my_plot_data$HR_low <- as.numeric(my_plot_data$HR_low)
  my_plot_data$HR_high <- as.numeric(my_plot_data$HR_high)


  ############################################
  ### CUSTOMIZE APPEARANCE WITH THESE     ####
  ############################################
  blankRows<-2    # blank rows under boxplot
  titleSize<-4
  dataSize<-4
  boxColor<-"pink"
  ############################################
  ############################################

  ## BASIC THEMES (SO TO PLOT BLANK GRID)
  theme_grid <- theme(
    axis.line = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks.length = unit(0.0001, "mm"),
    axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
    legend.position = "none", 
    panel.background = element_rect(fill = "transparent"), 
    panel.border = element_blank(), 
    panel.grid.major = element_line(colour="grey"), 
    panel.grid.minor = element_line(colour="grey"), 
    panel.margin = unit(c(-0.1,-0.1,-0.1,-0.1), "mm"), 
    plot.margin = unit(c(5,0,5,0.01), "mm")
  )

  theme_bare <- theme_grid +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )
  group_low <- my_plot_data[,c(1,3,7,8,9,10,12)]
  group_high <- my_plot_data[,c(1,2,7,8,9,10,12)]
  colnames(group_low) <- c('Gene', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  colnames(group_high) <- c('Gene', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  group_low$Group <- 'Low'
  group_high$Group <- 'High'
  group_high[,c(3,4,5)] <- 1
  group_high[,c(6,7)] <- NA
  group_data <- rbind(group_low, group_high)
  group_data <- group_data[order(group_data$Gene),]
  print(group_data)
  rownames(group_data) <- c(1:dim(group_data)[1])
  group_data$ID <- c(1:dim(group_data)[1])
  hazard_data<-expand.grid(ID=1:nrow(group_data),HR=1)
  hazard_data$HR <- group_data$HR
  hazard_data<-rbind(hazard_data,ddply(group_data,.(Gene),summarise,ID=max(ID)+0.1,HR=NA)[,2:3])

  hazard_data<-rbind(hazard_data,data.frame(ID=c(0,-1:(-2-blankRows),max(group_data$ID)+1,max(group_data$ID)+2),HR=NA))

  hazard_data$HR_low <- NA
  hazard_data$HR_high <- NA

  hazard_data[1:dim(group_data)[1],]$HR_low <- group_data$HR_low
  hazard_data[1:dim(group_data)[1],]$HR_high <- group_data$HR_high

  hr_labels <- group_data[,c('ID', 'HR (95% CI for HR)')]
  colnames(hr_labels) <- c("ID", "lab")

  upper <- max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)
  lower <- min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)
  lower <- min(0.4, lower)
  upper <- max(2.8, upper)

  upper_int <- round(((upper - 1)/2),2)
  lower_int <- round(((1 - lower)/2),2)

  scale_data <- data.frame(ID=0,HR=c((1-2*lower_int-0.05),(1-lower_int), 1, (1 + upper_int), (1+2*upper_int+0.2)))

  group_mirna<-ddply(group_data,.(Gene),summarise,y=max(ID)+0.1)


  hl_rows <- data.frame(ID=(1:floor(length(unique(hazard_data$ID[which(hazard_data$ID>0)]))/2))*2,col="lightgrey")
  hl_rows$ID <- hl_rows$ID+blankRows+1

  hl_rect <- function(col="white",alpha=0.5){
    rectGrob(   x = 0, y = 0, width = 1, height = 1, just = c("left","bottom"), gp=gpar(alpha=alpha, fill=col))
  }

  ## DATA FOR TEXT LABELS
  md_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=1)),
                        y=c(1),
                        lab=c("Hazard Ratio"),drop=FALSE)

  rt_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                        y=c(1,4),
                        lab=c("Hazard Ratio\n(95% CI)","P Value"))

  lf_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                       y=c(0.5,4),
                       lab=c("Gene","No. of\nPatients"))

  legend_labels <- data.frame(x=c(rep(1,times=2)),
                       y=c((1-lower_int),(1+upper_int)),
                       lab=c("Low Better","High Better"))

  haz <- ggplot(hazard_data,aes(factor(ID),HR))+ labs(x=NULL, y=NULL) 

  ## MIDDLE PANEL WITH LOG SCALE
  middle_panel <- haz +
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # geom_segment(aes(x = 2, y = 1, xend = 1.5, yend = 1)) + 
    geom_hline(aes(yintercept=1),linetype=2, size=0.5)+
    geom_point() + 
    geom_errorbar(data=hazard_data, aes(ymin=HR_low, ymax=HR_high), width=.3)+
    # geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+ 
    scale_y_log10() + 
    # scale_y_continuous() +
    coord_flip() +
    geom_text(data=scale_data,aes(3,HR,label=HR),vjust=0.5, size=dataSize) +
    geom_text(data=md_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=hr_labels,aes(factor(ID),0.4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=legend_labels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) +
    geom_point(data=scale_data,aes(2.5,HR),shape=3,size=3) + 
    geom_point(aes(2,6),shape=3,alpha=0,vjust=0) + 
    # geom_errorbar(width=.08)
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = max(8,max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)))) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare


  ## RIGHT PANEL WITH LOG SCALE
  rightPanel <- haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # scale_y_log10() +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=hr_labels,aes(factor(ID),3,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=rt_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare

  ## LEFT PANEL WITH NORMAL SCALE
  leftPanel<-haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=group_mirna,aes(factor(y),0.5,label=Gene, fontface="bold"),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),1,label=Group),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=NoP),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=lf_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, hjust=0, size=4, size=titleSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    theme_bare

  ## PLOT THEM BOTH IN A GRID SO THEY MATCH UP
  height <- 4000 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  resolution <- 400 + 100 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  print(hazard_data)
  # png(paste('result/TCGA.', sample,'.DEG.mRNA.all.HR.png', sep=""), height=height, width=6000, res=resolution)
  p <- ggarrange(leftPanel,middle_panel,rightPanel, widths=c(3,4,3), ncol=3, nrow=1)
  ggsave(paste('result/TCGA.', sample,'.DEG.mRNA.all.HR.png', sep=""),height=5,width=5)
  # print(p)
  # dev.off() 
}


out <- result[,c(1,2,3,10,12),drop=FALSE]

colnames(out) <- c('Gene', 'High', 'Low', 'HR (95% CI for HR)', 'p.value')
out$p.value <- as.numeric(out$p.value)
if(dim(out[which((out$p.value * 10) < 0.01),])[1] > 0){
    out[which(out$p.value < 0.001),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value < 0.001),]$`HR (95% CI for HR)`, ' ***')

  }
if(dim(out[which(out$p.value > 0.001 & out$p.value < 0.01),])[1] > 0){
  out[which(out$p.value > 0.001 & out$p.value < 0.01),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value > 0.001 & out$p.value < 0.01),]$`HR (95% CI for HR)`, ' **')
}
if(dim(out[which(out$p.value > 0.01 & out$p.value < 0.05),])[1] > 0){
  out[which(out$p.value > 0.01 & out$p.value < 0.05),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value > 0.01 & out$p.value < 0.05),]$`HR (95% CI for HR)`, ' *')
}

out$p.value <- factor(out$p.value)

out$High <- paste0(out$High, '(Ref)')
# out$Gene <- paste0(result$TCGA_ID, '_', out$Gene)
out <- as.matrix(out)
rownames(out) <- out[,1]
cgroup <- c('NoP', 'Univariate Model')
rgroup <- unique(result$TCGA_ID)
n.rgroup <- as.numeric(as.matrix(table(result$TCGA_ID))[,1])

out <- out[,-1,drop=FALSE]


colnames(out) <- c('High', 'Low', 'HR<br />(95% CI for HR)', '<i>p</i> value')
colnames(out) <- sprintf('<b>%s</b>', colnames(out))

out %>% 
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7"),
                    css.cell = "padding-left: .5em; padding-right: .5em; line-height: 1.8;color:black",
                    align="r") %>% 
  htmlTable(rowlabel = 'Tumor/Gene', 
            rgroup=rgroup,
            n.rgroup=n.rgroup,
            ctable = TRUE, align = 'cccc',
            ## number of columns that each cgroup label spans:
            n.cgroup = c(2, 2), cgroup = cgroup,
            caption = "<font size=3 color=black>Table: Hazard ratios and <i>p</i> values of cox model with DEG vs OS in COAD</font>", 
            tfoot = '<font size=1><sup>&dagger;</sup>* <i>p</i><0.05,** <i>p</i><0.01, *** <i>p</i><0.001.</font>') %>% 
  save_kable(file = "result/TCGA.table.png",zoom = 4)


