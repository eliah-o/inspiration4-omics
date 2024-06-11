library(ggplot2)

library(tidyverse)
library(dplyr)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ggsignif)
library(tidyverse)
library(RColorBrewer)


library("pheatmap")
library(ggplotify)
library(patchwork)
library(grid)

library(dplyr)
library(data.table)

all.list <- c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B',
             'NK', 'CD14_Mono', 'CD16_Mono', 'DC', 'other')

for (i in all.list){
    x <- read.csv(paste0("path_to_gene_expression_file.separated.by.timepoint.csv"))
    
    assign(paste0(i, ".gene.tp"), x)
}

for (i in all.list){
    x <- read.csv(paste0("path_to_gene_expression_file.separated.by.orig.ident.csv"))
    
    assign(paste0(i, ".gene.orig.ident"), x)
}

for (i in all.list){
    x <- read.csv(paste0("path_to_gene_expression_file.separated.by.timepoint.ver2.csv"))
    
    assign(paste0(i, ".gene.TP"), x)
}

for (i in all.list){
    for (j in c('C001', 'C002', 'C003', 'C004')){
        x <- read.csv(paste0("path_to_gene_expression_file.separated.by.timepoint.ver2.csv"))
    
        assign(paste0(i, '.', j, ".gene.TP.ID"), x)
}}



for (i in all.list){
    for (j in c('C001', 'C002', 'C003', 'C004')){
        x <- read.csv(paste0("../../Differential_expression/with_filter/celltype/pval/update/expression/average.", i, '.', j, ".peak.expression.TP.csv"))
    
        assign(paste0(i, '.', j, ".peak.TP.ID"), x)
}}



for (i in all.list){
    x <- read.csv(paste0("../../Differential_expression/with_filter/celltype/pval/update/expression/average.", i, ".peak.expression.tp.csv"))
    
    assign(paste0(i, ".peak.tp"), x)
}

for (i in all.list){
    x <- read.csv(paste0("../../Differential_expression/with_filter/celltype/pval/update/expression/average.", i, ".peak.expression.orig.ident.csv"))
    
    assign(paste0(i, ".peak.orig.ident"), x)
}

for (i in all.list){
    x <- read.csv(paste0("../../Differential_expression/with_filter/celltype/pval/update/expression/average.", i, ".peak.expression.TP.csv"))
    
    assign(paste0(i, ".peak.TP"), x)
}





for (i in all.list[! all.list %in% c('CD16_Mono', 'DC')]){
    x1 <- get(paste0(i, ".gene.tp"))
    x2 <- get(paste0(i, ".gene.orig.ident"))
    x3 <- get(paste0(i, ".peak.tp"))
    x4 <- get(paste0(i, ".peak.orig.ident"))
    
    x5 <- get(paste0(i, ".gene.TP"))
    x6 <- get(paste0(i, ".peak.TP"))
    
    y1 <- x1 %>% select(-c('X', 'gene'))
    y2 <- x2 %>% select(-c('X', 'gene'))
    y3 <- x3 %>% select(-c('X', 'gene'))
    y4 <- x4 %>% select(-c('X', 'gene'))
    
    y5 <- x5 %>% select(-c('X', 'gene'))
    y6 <- x6 %>% select(-c('X', 'gene'))
    
    colnames(y1) <- c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82')
    colnames(y2) <- colnames(pbmc.gene.orig.ident)[2:25]
    colnames(y3) <- c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82')
    colnames(y4) <- colnames(pbmc.peak.orig.ident)[2:25]
    
    colnames(y5) <- c('Pre-flight', 'R+1', 'R+45', 'R+82')
    colnames(y6) <- c('Pre-flight', 'R+1', 'R+45', 'R+82')
    
    assign(paste0(i, '.gene.tp.2'), y1)
    assign(paste0(i, '.gene.orig.ident.2'), y2)
    assign(paste0(i, '.peak.tp.2'), y3)
    assign(paste0(i, '.peak.orig.ident.2'), y4)
    
    assign(paste0(i, '.gene.TP.2'), y5)
    assign(paste0(i, '.peak.TP.2'), y6)
    
    
    z1 <- prcomp(y1, scale.=TRUE)
    z2 <- prcomp(y2, scale.=TRUE)
    z3 <- prcomp(y3, scale.=TRUE)
    z4 <- prcomp(y4, scale.=TRUE)
    
    z5 <- prcomp(y5, scale.=TRUE)
    z6 <- prcomp(y6, scale.=TRUE)
    
    assign(paste0(i, ".pr.gene.tp"), z1)
    assign(paste0(i, ".pr.gene.orig.ident"), z2)
    assign(paste0(i, ".pr.peak.tp"), z3)
    assign(paste0(i, ".pr.peak.orig.ident"), z4)
    assign(paste0(i, ".pr.gene.TP"), z5)
    assign(paste0(i, ".pr.peak.TP"), z6)
    
    
    
    
    df1 <- as.data.frame(z1$rotation)
    df2 <- as.data.frame(z2$rotation)
    df3 <- as.data.frame(z3$rotation)
    df4 <- as.data.frame(z4$rotation)
    
    df5 <- as.data.frame(z5$rotation)
    df6 <- as.data.frame(z6$rotation)
    
    
    DF1 <- data.frame(Timepoint <- rownames(df1), PC1 <- df1$PC1, PC2 <- df1$PC2)
    colnames(DF1) <- c("Timepoint", 'PC1', 'PC2')
    DF2 <- data.frame(Timepoint <- rownames(df2), PC1 <- df2$PC1, PC2 <- df2$PC2)
    colnames(DF2) <- c("Timepoint", 'PC1', 'PC2')
    DF3 <- data.frame(Timepoint <- rownames(df3), PC1 <- df3$PC1, PC2 <- df3$PC2)
    colnames(DF3) <- c("Timepoint", 'PC1', 'PC2')
    DF4 <- data.frame(Timepoint <- rownames(df4), PC1 <- df4$PC1, PC2 <- df4$PC2)
    colnames(DF4) <- c("Timepoint", 'PC1', 'PC2')
    
    
    DF5 <- data.frame(Timepoint <- rownames(df5), PC1 <- df5$PC1, PC2 <- df5$PC2)
    colnames(DF5) <- c("Timepoint", 'PC1', 'PC2')
    DF6 <- data.frame(Timepoint <- rownames(df6), PC1 <- df6$PC1, PC2 <- df6$PC2)
    colnames(DF6) <- c("Timepoint", 'PC1', 'PC2')
    
    v1 <- data.frame(sd = z1$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v2 <- data.frame(sd = z2$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v3 <- data.frame(sd = z3$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v4 <- data.frame(sd = z4$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    
    v5 <- data.frame(sd = z5$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v6 <- data.frame(sd = z6$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    
    
    assign(paste0(i, ".df.gene.tp"), DF1)
    assign(paste0(i, ".df.gene.orig.ident"), DF2)
    assign(paste0(i, ".df.peak.tp"), DF3)
    assign(paste0(i, ".df.peak.orig.ident"), DF4)
    
    assign(paste0(i, ".df.gene.TP"), DF5)
    assign(paste0(i, ".df.peak.TP"), DF6)
    
    assign(paste0(i, ".gene.tp.pct"), v1)
    assign(paste0(i, ".gene.orig.ident.pct"), v2)
    assign(paste0(i, ".peak.tp.pct"), v3)
    assign(paste0(i, ".peak.orig.ident.pct"), v4)
    
    assign(paste0(i, ".gene.TP.pct"), v5)
    assign(paste0(i, ".peak.TP.pct"), v6)
    

    
}

for (i in c('CD16_Mono', 'DC')){
    x1 <- get(paste0(i, ".gene.tp"))
    x2 <- get(paste0(i, ".gene.orig.ident"))
    x3 <- get(paste0(i, ".peak.tp"))
    x4 <- get(paste0(i, ".peak.orig.ident"))
    
    x5 <- get(paste0(i, ".gene.TP"))
    x6 <- get(paste0(i, ".peak.TP"))
    
    y1 <- x1 %>% select(-c('X', 'gene'))
    y2 <- x2 %>% select(-c('X', 'gene'))
    y3 <- x3 %>% select(-c('X', 'gene'))
    y4 <- x4 %>% select(-c('X', 'gene'))
    
    y5 <- x5 %>% select(-c('X', 'gene'))
    y6 <- x6 %>% select(-c('X', 'gene'))
    
    colnames(y1) <- c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82')
    colnames(y2) <- colnames(DC.gene.orig.ident)[2:24]
    colnames(y3) <- c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82')
    colnames(y4) <- colnames(DC.peak.orig.ident)[2:24]
    
    colnames(y5) <- c('Pre-flight', 'R+1', 'R+45', 'R+82')
    colnames(y6) <- c('Pre-flight', 'R+1', 'R+45', 'R+82')
    
    assign(paste0(i, '.gene.tp.2'), y1)
    assign(paste0(i, '.gene.orig.ident.2'), y2)
    assign(paste0(i, '.peak.tp.2'), y3)
    assign(paste0(i, '.peak.orig.ident.2'), y4)
    
    assign(paste0(i, '.gene.TP.2'), y5)
    assign(paste0(i, '.peak.TP.2'), y6)
    
    
    z1 <- prcomp(y1, scale.=TRUE)
    z2 <- prcomp(y2, scale.=TRUE)
    z3 <- prcomp(y3, scale.=TRUE)
    z4 <- prcomp(y4, scale.=TRUE)
    
    z5 <- prcomp(y5, scale.=TRUE)
    z6 <- prcomp(y6, scale.=TRUE)
    
    assign(paste0(i, ".pr.gene.tp"), z1)
    assign(paste0(i, ".pr.gene.orig.ident"), z2)
    assign(paste0(i, ".pr.peak.tp"), z3)
    assign(paste0(i, ".pr.peak.orig.ident"), z4)
    assign(paste0(i, ".pr.gene.TP"), z5)
    assign(paste0(i, ".pr.peak.TP"), z6)
    
    
    
    
    df1 <- as.data.frame(z1$rotation)
    df2 <- as.data.frame(z2$rotation)
    df3 <- as.data.frame(z3$rotation)
    df4 <- as.data.frame(z4$rotation)
    
    df5 <- as.data.frame(z5$rotation)
    df6 <- as.data.frame(z6$rotation)
    
    
    DF1 <- data.frame(Timepoint <- rownames(df1), PC1 <- df1$PC1, PC2 <- df1$PC2)
    colnames(DF1) <- c("Timepoint", 'PC1', 'PC2')
    DF2 <- data.frame(Timepoint <- rownames(df2), PC1 <- df2$PC1, PC2 <- df2$PC2)
    colnames(DF2) <- c("Timepoint", 'PC1', 'PC2')
    DF3 <- data.frame(Timepoint <- rownames(df3), PC1 <- df3$PC1, PC2 <- df3$PC2)
    colnames(DF3) <- c("Timepoint", 'PC1', 'PC2')
    DF4 <- data.frame(Timepoint <- rownames(df4), PC1 <- df4$PC1, PC2 <- df4$PC2)
    colnames(DF4) <- c("Timepoint", 'PC1', 'PC2')
    
    
    DF5 <- data.frame(Timepoint <- rownames(df5), PC1 <- df5$PC1, PC2 <- df5$PC2)
    colnames(DF5) <- c("Timepoint", 'PC1', 'PC2')
    DF6 <- data.frame(Timepoint <- rownames(df6), PC1 <- df6$PC1, PC2 <- df6$PC2)
    colnames(DF6) <- c("Timepoint", 'PC1', 'PC2')
    
    v1 <- data.frame(sd = z1$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v2 <- data.frame(sd = z2$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v3 <- data.frame(sd = z3$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v4 <- data.frame(sd = z4$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    
    v5 <- data.frame(sd = z5$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    v6 <- data.frame(sd = z6$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))
    
    
    assign(paste0(i, ".df.gene.tp"), DF1)
    assign(paste0(i, ".df.gene.orig.ident"), DF2)
    assign(paste0(i, ".df.peak.tp"), DF3)
    assign(paste0(i, ".df.peak.orig.ident"), DF4)
    
    assign(paste0(i, ".df.gene.TP"), DF5)
    assign(paste0(i, ".df.peak.TP"), DF6)
    
    assign(paste0(i, ".gene.tp.pct"), v1)
    assign(paste0(i, ".gene.orig.ident.pct"), v2)
    assign(paste0(i, ".peak.tp.pct"), v3)
    assign(paste0(i, ".peak.orig.ident.pct"), v4)
    
    assign(paste0(i, ".gene.TP.pct"), v5)
    assign(paste0(i, ".peak.TP.pct"), v6)
    

    
}





for (i in all.list){
    
    x1 <- get(paste0(i, '.C001.gene.TP.ID'))
    x2 <- get(paste0(i, '.C002.gene.TP.ID'))
    x3 <- get(paste0(i, '.C003.gene.TP.ID'))
    x4 <- get(paste0(i, '.C004.gene.TP.ID'))
    
    colnames(x1) <- c('X', 'Pre.flight - C001', 'R.1 - C001', 'R.45 - C001', 'R.82 - C001', 'gene')
    colnames(x2) <- c('X', 'Pre.flight - C002', 'R.1 - C002', 'R.45 - C002', 'R.82 - C002', 'gene')
    colnames(x3) <- c('X', 'Pre.flight - C003', 'R.1 - C003', 'R.45 - C003', 'R.82 - C003', 'gene')
    colnames(x4) <- c('X', 'Pre.flight - C004', 'R.1 - C004', 'R.45 - C004', 'R.82 - C004', 'gene')
    
    y <- cbind(x1,x2,x3,x4)
    
    z <- y %>% select('X', 'Pre.flight - C001', 'R.1 - C001', 'R.45 - C001', 'R.82 - C001',
                                   'Pre.flight - C002', 'R.1 - C002', 'R.45 - C002', 'R.82 - C002',
                                   'Pre.flight - C003', 'R.1 - C003', 'R.45 - C003', 'R.82 - C003',
                                   'Pre.flight - C004', 'R.1 - C004', 'R.45 - C004', 'R.82 - C004', 'gene')
    
    
    
    
    assign(paste0(i, ".gene.TP.ID"), z)
}

for (i in all.list){
    
    x1 <- get(paste0(i, '.C001.peak.TP.ID'))
    x2 <- get(paste0(i, '.C002.peak.TP.ID'))
    x3 <- get(paste0(i, '.C003.peak.TP.ID'))
    x4 <- get(paste0(i, '.C004.peak.TP.ID'))
    
    colnames(x1) <- c('X', 'Pre.flight - C001', 'R.1 - C001', 'R.45 - C001', 'R.82 - C001', 'gene')
    colnames(x2) <- c('X', 'Pre.flight - C002', 'R.1 - C002', 'R.45 - C002', 'R.82 - C002', 'gene')
    colnames(x3) <- c('X', 'Pre.flight - C003', 'R.1 - C003', 'R.45 - C003', 'R.82 - C003', 'gene')
    colnames(x4) <- c('X', 'Pre.flight - C004', 'R.1 - C004', 'R.45 - C004', 'R.82 - C004', 'gene')
    
    y <- cbind(x1,x2,x3,x4)
    
    z <- y %>% select('X', 'Pre.flight - C001', 'R.1 - C001', 'R.45 - C001', 'R.82 - C001',
                                   'Pre.flight - C002', 'R.1 - C002', 'R.45 - C002', 'R.82 - C002',
                                   'Pre.flight - C003', 'R.1 - C003', 'R.45 - C003', 'R.82 - C003',
                                   'Pre.flight - C004', 'R.1 - C004', 'R.45 - C004', 'R.82 - C004', 'gene')
    
    
    
    
    assign(paste0(i, ".peak.TP.ID"), z)
}







for (i in all.list){
    x1 <- get(paste0(i, ".gene.TP.ID"))
    
    y1 <- x1 %>% select(-c('X', 'gene'))
    
    colnames(y1) <- c('Pre-flight', 'R+1', 'R+45', 'R+82')
    
    assign(paste0(i, '.gene.TP.ID.2'), y1)
    
    
    z1 <- prcomp(y1, scale.=TRUE)
    
    assign(paste0(i, ".pr.gene.TP.ID"), z1)
    
    
    
    
    df1 <- as.data.frame(z1$rotation)
    
    
    DF1 <- data.frame(Timepoint <- rownames(df1), PC1 <- df1$PC1, PC2 <- df1$PC2)
    colnames(DF1) <- c("Timepoint", 'PC1', 'PC2')
    
    DF1$timepoint <- c('Pre-flight', 'R+1', 'R+45', 'R+82',
                         'Pre-flight', 'R+1', 'R+45', 'R+82',
                         'Pre-flight', 'R+1', 'R+45', 'R+82',
                         'Pre-flight', 'R+1', 'R+45', 'R+82')
    
    DF1$ID <- c('C001', 'C001', 'C001', 'C001',
                         'C002', 'C002', 'C002', 'C002',
                         'C003', 'C003', 'C003', 'C003',
                         'C004', 'C004', 'C004', 'C004')
    
    
    v1 <- data.frame(sd = z1$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))

    
    
    assign(paste0(i, ".df.gene.TP.ID"), DF1)
    
    assign(paste0(i, ".gene.TP.ID.pct"), v1)
    

    
}

for (i in all.list){
    x1 <- get(paste0(i, ".peak.TP.ID"))
    
    y1 <- x1 %>% select(-c('X', 'gene'))
    
    colnames(y1) <- c('Pre-flight', 'R+1', 'R+45', 'R+82')
    
    assign(paste0(i, '.peak.TP.ID.2'), y1)
    
    
    z1 <- prcomp(y1, scale.=TRUE)
    
    assign(paste0(i, ".pr.peak.TP.ID"), z1)
    
    
    
    
    df1 <- as.data.frame(z1$rotation)
    
    
    DF1 <- data.frame(Timepoint <- rownames(df1), PC1 <- df1$PC1, PC2 <- df1$PC2)
    colnames(DF1) <- c("Timepoint", 'PC1', 'PC2')
    
    DF1$timepoint <- c('Pre-flight', 'R+1', 'R+45', 'R+82',
                         'Pre-flight', 'R+1', 'R+45', 'R+82',
                         'Pre-flight', 'R+1', 'R+45', 'R+82',
                         'Pre-flight', 'R+1', 'R+45', 'R+82')
    
    DF1$ID <- c('C001', 'C001', 'C001', 'C001',
                         'C002', 'C002', 'C002', 'C002',
                         'C003', 'C003', 'C003', 'C003',
                         'C004', 'C004', 'C004', 'C004')
    
    
    v1 <- data.frame(sd = z1$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd)))

    
    
    assign(paste0(i, ".df.peak.TP.ID"), DF1)
    
    assign(paste0(i, ".peak.TP.ID.pct"), v1)
    

    
}





for (i in all.list){
    DF1 <- get(paste0(i, ".df.gene.TP.ID"))
    DF2 <- get(paste0(i, ".df.peak.TP.ID"))
    
    DF1$celltype <- i
    DF2$celltype <- i
    
    assign(paste0(i, '.gene'), DF1)
    assign(paste0(i, '.peak'), DF2)

    
    
}



df.gex <- rbind(pbmc.gene, CD4_T.gene, CD8_T.gene, other_T.gene,
               B.gene, NK.gene, CD14_Mono.gene, CD16_Mono.gene,
               DC.gene)

df.ATAC <- rbind(pbmc.peak, CD4_T.peak, CD8_T.peak, other_T.peak,
               B.peak, NK.peak, CD14_Mono.peak, CD16_Mono.peak,
               DC.peak)




options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list){
    DF1 <- get(paste0(i, ".df.gene.TP.ID"))
    v1 <- get(paste0(i, ".gene.TP.ID.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = ID, color = timepoint), size = 6)+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
    print(a)
    
    
    
    
}

options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list){
    DF1 <- get(paste0(i, ".df.peak.TP.ID"))
    v1 <- get(paste0(i, ".peak.TP.ID.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = ID, color = timepoint), size = 6)+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
    print(a)
    
    
    
    
    
}







options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list){
    DF1 <- get(paste0(i, ".df.gene.tp"))
    v1 <- get(paste0(i, ".gene.tp.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = Timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 2)+
  geom_line()+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
    print(a)
    
    
     
    
}

options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list){
    DF1 <- get(paste0(i, ".df.peak.tp"))
    v1 <- get(paste0(i, ".peak.tp.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = Timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 2)+
  geom_line()+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
   
}



options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list){
    DF1 <- get(paste0(i, ".df.gene.TP"))
    v1 <- get(paste0(i, ".gene.TP.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = Timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 2)+
  geom_line()+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
    
    
}

options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list){
    DF1 <- get(paste0(i, ".df.peak.TP"))
    v1 <- get(paste0(i, ".peak.TP.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = Timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 2)+
  geom_line()+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", 
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
   
    
}





options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list[! all.list %in% c('CD16_Mono', 'DC')]){
    DF1 <- get(paste0(i, ".df.gene.orig.ident"))
    v1 <- get(paste0(i, ".gene.orig.ident.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    DF1$timepoint <- c(rep(c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82'), 4))
    DF1$timepoint <- factor(DF1$timepoint,
                           levels=unique(DF1$timepoint))
    
    
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 1)+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
   
    
}





options(repr.plot.width=10, repr.plot.height=8)

for (i in c('CD16_Mono', 'DC')){
    DF1 <- get(paste0(i, ".df.gene.orig.ident"))
    v1 <- get(paste0(i, ".gene.orig.ident.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    DF1$timepoint <- c('L-92', 'L-44', 'R+1', 'R+45', 'R+82',
                      'L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82',
                      'L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82',
                      'L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82')
    DF1$timepoint <- factor(DF1$timepoint,
                           levels=c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82'))
    
    
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 1)+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
    
    
}





options(repr.plot.width=10, repr.plot.height=8)

for (i in all.list[! all.list %in% c('CD16_Mono', 'DC')]){
    DF1 <- get(paste0(i, ".df.peak.orig.ident"))
    v1 <- get(paste0(i, ".peak.orig.ident.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    DF1$timepoint <- c(rep(c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82'), 4))
    DF1$timepoint <- factor(DF1$timepoint,
                           levels=unique(DF1$timepoint))
    
    
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 1)+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
   
    
}

options(repr.plot.width=10, repr.plot.height=8)

for (i in c('CD16_Mono', 'DC')){
    DF1 <- get(paste0(i, ".df.peak.orig.ident"))
    v1 <- get(paste0(i, ".peak.orig.ident.pct"))
    
    DF1$Timepoint <- factor(DF1$Timepoint,
                           levels=unique(DF1$Timepoint))
    
    DF1$timepoint <- c('L-92', 'L-44', 'R+1', 'R+45', 'R+82',
                      'L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82',
                      'L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82',
                      'L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82')
    DF1$timepoint <- factor(DF1$timepoint,
                           levels=c('L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82'))
    
    
    
    
    
    
    a <- ggplot(DF1, aes(x = PC1, fill = timepoint, y = PC2)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 1)+
    theme_bw()+
        labs(x = paste0("PC1 (", round(v1$pct[1],1), "%)"), y = paste0("PC2 (", round(v1$pct[2],1), "%)"), fill = 'Timepoint')+ 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(size = 24, vjust = 0.5, hjust = 0, angle = -90), 
        axis.text.y = element_text(size=24))+ 
  scale_fill_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))+ 
  scale_color_manual(values=c("#68BBDA", "#377EB8", "#3E4D8A",
                            '#E41A1C', '#D9742C', '#ECB03B'))
    
      
    
}






