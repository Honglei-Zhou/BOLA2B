library(IOBR)
library(EPIC)
library(estimate) 
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(psych)
library(maftools)

# my_df: contains gene expression matrix

timer<-deconvo_tme(eset = my_df, method = "timer", group_list = rep("stad",dim(my_df)[2]))

#> # A tibble: 6 x 7
#>   ID           B_cell_TIMER T_cell_CD4_TIMER T_cell_CD8_TIMER Neutrophil_TIMER
#>   <chr>               <dbl>            <dbl>            <dbl>            <dbl>
#> 1 TCGA-B7-5818       0.0954            0.130            0.191            0.114
#> 2 TCGA-BR-4187       0.0983            0.131            0.202            0.125
#> 3 TCGA-BR-4201       0.0996            0.122            0.200            0.127
#> 4 TCGA-BR-4253       0.101             0.131            0.240            0.129
#> 5 TCGA-BR-4256       0.0945            0.133            0.213            0.137
#> 6 TCGA-BR-4257       0.0907            0.126            0.199            0.121
#> # … with 2 more variables: Macrophage_TIMER <dbl>, DC_TIMER <dbl>

timer_cor_result = NULL
for(cell_type in colnames(timer)){
	cor.result <- cor.test(timer[,cell_type],my_df[,'BOLA2B '],method='pearson')
	tmp <- data.frame('r'=t$estimate,'pval'=t$p.value,'cell_type'=cell_type)
	if(is.null(timer_cor_result)){
		timer_cor_result <- tmp
	}else{
		timer_cor_result <- rbind(timer_cor_result,tmp)
	}
}

write.csv(timer_cor_result,'timer.csv')

# ESTIMATE
estimate<-deconvo_tme(eset = my_df, method = "estimate")
#> # A tibble: 6 x 5
#>   ID      StromalScore_est… ImmuneScore_esti… ESTIMATEScore_es… TumorPurity_est…
#>   <chr>               <dbl>             <dbl>             <dbl>            <dbl>
#> 1 TCGA-B…             -111.             1762.             1651.            0.662
#> 2 TCGA-B…             2054.             2208.             4262.            0.334
#> 3 TCGA-B…             1411.             1770.             3181.            0.479
#> 4 TCGA-B…              483.             2905.             3389.            0.451
#> 5 TCGA-B…             1659.             2541.             4200.            0.342
#> 6 TCGA-B…              831.             1722.             2553.            0.557


df <- cbind(estimate[,'ESTIMATEScore_estimate',drop=FALSE],my_df[,'BOLA2B',drop=FALSE])

sp <- ggscatter(df, x = "BOLA2B", y = "ESTIMATEScore_estimate",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson", label.x = 3, label.y = 30)#> `geom_smooth()` using formula 'y ~ x'


# TMB
tmb_score <- tmb(maf, captureSize = 50, logScale = TRUE)


cor.result <- cor.test(tmb_score[,'total_perMB'],my_df[,'BOLA2B '],method='pearson')
tmb_cor_result <- data.frame('r'=t$estimate,'pval'=t$p.value,'TCGA_ID'=tcga_id)

write.csv(tmb_cor_result,'tmb.csv')


