#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#April 1, 20201

# Performs statistical analysis of CETSA mass spec data via area under the curve comparison
# This script calculates:
# i. area under the curve via trapezoid method
# ii. Statistics via Welch Two Sample t-test
# iii. Qvalues using bioconductor's qvalue package (not working)

# OUTPUT: 
# i. processed excel workbook with multiple sheets
# ii. plots of melting curves on the replicate level
# iii. plots of melting curves on the protein level - averaged across replicates

### Begin Script ###

#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2021_0330_BRC4/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2021_0330_BRC4/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES

packages = c("stringr", "hablar", 
             "readxl", "writexl", "openxlsx",
             "dplyr", "tidyr")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# to generate qvalues
BiocManager::install("qvalue")
suppressPackageStartupMessages(library(qvalue))
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# load data

# .tsv BGS Report of CETSA data on the peptide level # for very large files like this, make sure it is saved in either .tsv or .csv file format to prevent data loss
df.input <- read.csv("//bigrock/GibsonLab/users/birgit/Barecia/BRC4_CDK/21_0325_BRC4_V1/Spectronaut/20210325_174612_210325_BRC4_all files_V1/20210401_164439_210325_BRC4_all_files_V1_BGS_V2_Report.tsv", 
                     sep="\t", stringsAsFactors = FALSE)

# check out the unique files names
df.input$R.FileName %>% unique()
# how many files names?
df.input$R.FileName %>% unique() %>% length()

# candidates file (which contains the number of unique peptides for each protein)
# only need to use candidates file if the BGS Report does not contain number of unique peptides information (Christina will include this variable in future BGS Reports)
df.candidates <- read.csv("//bigrock/GibsonLab/users/birgit/Barecia/BRC4_CDK/21_0325_BRC4_V1/Spectronaut/20210325_174612_210325_BRC4_all files_V1/Candidates.tsv", sep="\t", stringsAsFactors = FALSE)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# remove one peptide wonders

# one peptide wonders
wonder.prots <- df.candidates %>%
  filter(X..Unique.Total.Peptides == 1) %>%
  pull(UniProtIds) %>%
  unique()

# remove one peptide wonders from input data
df.input <- df.input %>%
  filter(! PG.ProteinAccessions %in% wonder.prots) # exclude one peptide wonders

# remove one peptide wonders from candidates data
df.candidates <- df.candidates %>%
  filter(X..Unique.Total.Peptides > 1) # exclude one peptide wonders
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# prepare data

# split treatment group and temperature from R.Condition variable
df <- df.input %>%
  separate(col = R.Condition, into = c("Treatment", "Temperature"), sep="-") %>% # treatment and temperature separated by dash
  mutate(Temperature = as.numeric(Temperature)) %>% # convert temperature variable to numeric
  group_by(PG.ProteinAccessions, PG.Genes, R.Replicate, Treatment, Temperature) %>% # sum FG.Quantity across peptides
  summarise(Total.FG.Quantity = sum(FG.Quantity)) %>% 
  ungroup() 
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# check out data a bit before proceeding with analysis

# how many temperature points and proteins were measured for each condition and replicate?
df %>%
  group_by(R.Replicate, Treatment) %>%
  mutate(Number.of.Temperature.Points = unique(Temperature) %>% length()) %>%
  mutate(Number.of.Proteins = unique(PG.ProteinAccessions) %>% length()) %>%
  select(R.Replicate, Treatment, Number.of.Temperature.Points, Number.of.Proteins) %>%
  unique()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# calculate area under the curve (AUC) 

# first calcualte areas via trapezoid method
df.trapezoid <- df %>%
  group_by(PG.ProteinAccessions, PG.Genes, R.Replicate, Treatment) %>%
  mutate(TEMPERATURE_INITIAL = Temperature) %>% # define initial temp
  mutate(TEMPERATURE_FINAL = lead(Temperature)) %>% # define final temp
  mutate(FG.QUANTITY_INITIAL = Total.FG.Quantity) %>% # define initial base
  mutate(FG.QUANTITY_FINAL = lead(Total.FG.Quantity)) %>% # define final base
  mutate(FG.QUANTITY_INITIAL_TEMP = Total.FG.Quantity[1]) %>% # define initial quantity at the initial temperature
  filter(!is.na(TEMPERATURE_FINAL)) %>% # remove last row with final temperature as it is no longer needed
  mutate(TEMPERATURE_CHANGE = TEMPERATURE_FINAL - TEMPERATURE_INITIAL) %>% # calculate change in temperature (which we can think about as the trapezoid height)
  mutate(TRAPEZOID_AREA = (FG.QUANTITY_FINAL + FG.QUANTITY_INITIAL)*TEMPERATURE_CHANGE/2) %>% # calculate trapezoid area which is = height*(base1+base2)/2
  mutate(TRAPEZOID_AREA_NORM = TRAPEZOID_AREA/FG.QUANTITY_INITIAL_TEMP) %>% # calculate normalized trapezoid area
  select(-Temperature, -Total.FG.Quantity) # drop these variables

# calculate total area under curve
df.auc <- df.trapezoid %>%
  group_by(PG.ProteinAccessions, PG.Genes, R.Replicate, Treatment) %>%
  summarise(AREA_UNDER_CURVE_NORMALIZED = sum(TRAPEZOID_AREA_NORM)) 
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# statistics
# perform Welch Two Sample t-test comparing normalized areas between treatment groups for each protein

# define vectors for loop
proteins <- df.auc$PG.ProteinAccessions %>% unique()
treatments <- df.auc$Treatment %>% unique()

# initialize output dataframe
df.statistics <- as.data.frame(matrix(NA, nrow=length(proteins), ncol=6))
names(df.statistics) <- c("ProteinAccessions", "Genes", "Normalized_Area_Treatment1", "Normalized_Area_Treatment2", "Pvalue", "T_Statistic")
names(df.statistics)[3:4] <- paste("Normalized_Area", treatments, sep="_") # customize these variable names using treatment names

for(i in seq_along(proteins)){
  # subset data by protein
  df.loop <- df.auc %>%
    filter(PG.ProteinAccessions==proteins[i])
  
  # split data by treatment
  loop.split <- split(df.loop, df.loop$Treatment)
  
  # grab areas 
  areas_x <- loop.split[[1]] %>% pull(AREA_UNDER_CURVE_NORMALIZED)
  areas_y <- loop.split[[2]] %>% pull(AREA_UNDER_CURVE_NORMALIZED)
  
  # perform Welch Two Sample t-test
  ttest <- t.test(areas_x, areas_y)
  
  # grab statistical results
  pvalue <- ttest$p.value
  t_stat <- ttest$statistic
  
  # write out results
  df.statistics[i, "ProteinAccessions"] <- df.loop %>% pull(PG.ProteinAccessions) %>% unique()
  df.statistics[i, "Genes"] <- df.loop %>% pull(PG.Genes) %>% unique()
  df.statistics[i, "Normalized_Area_CTL"] <- areas_x %>% mean()
  df.statistics[i, "Normalized_Area_NMN"] <- areas_y %>% mean()
  df.statistics[i, "Pvalue"] <- pvalue
  df.statistics[i, "T_Statistic"] <- t_stat
}

# arrange by Pvalue
df.statistics <- df.statistics %>%
  arrange(Pvalue)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# example using broom::tidy to run t.test on a lot dataset ... look into this?
# set.seed(111)
# exampleDF <- data.frame(s1 = rnorm(10), s2 = rnorm(10), s3 = rnorm(10), s4 = rnorm(10), s5 = rnorm(10))
# group = factor(rep(c("A","B"),c(3,2)))
# res = apply(exampleDF,1,function(i)tidy(t.test(i ~ group))) # broom::tidy
# res = do.call(rbind,res)
# res$adjP = p.adjust(res$p.value,"BH")
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# write out statistics, area under curve and total.fg.quantity tables
write.xlsx(list(df.statistics, df.auc, df), "Table_statistics_results.xlsx")
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# plot significant proteins - individual replicates

# filter for significant proteins P<0.05
df.sig <- df.statistics %>%
  filter(Pvalue<0.05)

# define significant proteins for plotting
proteins <- df.sig$ProteinAccessions %>% unique()

# grab data to plot from df dataframe
df.plot <- df %>%
  filter(PG.ProteinAccessions %in% proteins)

# define temperatures
temps <- df.plot$Temperature %>% unique()

# define replicates
reps <- df.plot$R.Replicate %>% unique()

#inititate PDF
pdf(file="CETSA_AUC_Significant_Replicates.pdf")
par(mfrow=c(2,3)) # two rows by three columns

for(i in seq_along(proteins)){
  for(j in seq_along(reps)){
    
    # define data for this iteration
    df.loop <- df.plot %>%
      filter(PG.ProteinAccessions==proteins[i]) %>%
      filter(R.Replicate==reps[j])
      
    # define x-axis values which is temperature
    x <- df.loop %>% pull(Temperature) %>% unique()
    
    # split data by treatment group
    df.first <- split(df.loop, df.loop$Treatment)[[1]] %>% mutate(COLOR = "red") # add color for plotting
    df.second <- split(df.loop, df.loop$Treatment)[[2]] %>% mutate(COLOR = "blue") # add color for plotting
    
    # define y-data by treatment group
    y_first <- df.first %>% pull(Total.FG.Quantity)
    y_second <- df.second %>% pull(Total.FG.Quantity)
    
    # normalize y-data
    y_first_n <- y_first/y_first[1]
    y_second_n <- y_second/y_second[1]
    
    # create title
    gene <- df.loop %>% pull(PG.Genes) %>% unique() # grab gene name
    plot.title <- paste0(gene, " Replicate = ", reps[j]) # combine gene and replicate number
    
    # plot
    plot(0, xlim=c(min(x), max(x)), ylim=c(0,max(c(y_first_n, y_second_n))*1.4), xlab=expression("Temperature " (degree*C)), ylab="Relative Abundance", main = plot.title)
    
    # plot data points
    points(x, y_first_n, col=df.first$COLOR %>% unique())
    points(x, y_second_n, col=df.second$COLOR %>% unique())
    
    #plot line connecting the data points
    lines(x, y_first_n, col=df.first$COLOR %>% unique())
    lines(x, y_second_n, col=df.second$COLOR %>% unique())
    
    # legend
    legend("topleft", legend=unique(df.loop$Treatment), col=c(df.first$COLOR %>% unique(), df.second$COLOR %>% unique()), lty=1, cex=0.8)
    
    # #flag the significant data points
    # shift <- max(c(y_tn,y_cn))*0.15 #to shift flags up the same amount
    # points(x[y_q], y_tn[y_q]+shift, pch=8) #adds asterisks (flags)
    
  } # end for replicate loop
} # end for protein loop
graphics.off()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# plot significant proteins - averaged across replicates

# filter for significant proteins P<0.05
df.sig <- df.statistics %>%
  filter(Pvalue<0.05)

# define significant proteins for plotting
proteins <- df.sig$ProteinAccessions %>% unique()

# grab data to plot from df dataframe
df.plot <- df %>%
  filter(PG.ProteinAccessions %in% proteins)

# define temperatures
temps <- df.plot$Temperature %>% unique()

# inititate PDF
pdf(file="CETSA_AUC_Significant_Average.pdf")
par(mfrow=c(2,3)) # two rows by three columns

for(i in seq_along(proteins)){
  # define data for this iteration
  df.loop <- df.plot %>%
    filter(PG.ProteinAccessions==proteins[i]) %>%
    group_by(PG.ProteinAccessions, PG.Genes, Treatment, Temperature) %>%
    summarise(AVERAGE_TOTAL.FG.QUANTITY = mean(Total.FG.Quantity)) # calculate average across replicates
  
  # define x-axis values which is temperature
  x <- df.loop %>% pull(Temperature) %>% unique()
  
  # split data by treatment group
  df.first <- split(df.loop, df.loop$Treatment)[[1]] %>% mutate(COLOR = "red") # add color for plotting
  df.second <- split(df.loop, df.loop$Treatment)[[2]] %>% mutate(COLOR = "blue") # add color for plotting
  
  # define y-data by treatment group
  y_first <- df.first %>% pull(AVERAGE_TOTAL.FG.QUANTITY)
  y_second <- df.second %>% pull(AVERAGE_TOTAL.FG.QUANTITY)
  
  # normalize y-data
  y_first_n <- y_first/y_first[1]
  y_second_n <- y_second/y_second[1]
  
  # create title
  plot.title <- df.loop %>% pull(PG.Genes) %>% unique() # grab gene name
  
  # plot
  plot(0, xlim=c(min(x), max(x)), ylim=c(0,max(c(y_first_n, y_second_n))*1.4), xlab=expression("Temperature " (degree*C)), ylab="Relative Abundance", main = plot.title)
  
  # plot data points
  # points(x, y_first_n, col="red")
  # points(x, y_second_n, col="blue")
  points(x, y_first_n, col=df.first$COLOR %>% unique())
  points(x, y_second_n, col=df.second$COLOR %>% unique())
  
  #plot line connecting the data points
  # lines(x, y_first_n, col="red")
  # lines(x, y_second_n, col="blue")
  lines(x, y_first_n, col=df.first$COLOR %>% unique())
  lines(x, y_second_n, col=df.second$COLOR %>% unique())
  
  # legend
  # legend("topleft", legend=unique(df.loop$Treatment), col=c("red", "blue"), lty=1, cex=0.8)
  legend("topleft", legend=unique(df.loop$Treatment), col=c(df.first$COLOR %>% unique(), df.second$COLOR %>% unique()), lty=1, cex=0.8)
  
  # legend pvalue
  # grab pvalue
  pvalue <- df.sig %>%
    filter(ProteinAccessions==proteins[i]) %>%
    pull(Pvalue) %>%
    unique()
  pvalue.title <- paste0("Pvalue = ", pvalue %>% round(digits=3)) # make title for pvalue, round to 3 decimals
  legend("topright", legend=pvalue.title, col="black", cex=0.8)
  
  # #flag the significant data points
  # shift <- max(c(y_tn,y_cn))*0.15 #to shift flags up the same amount
  # points(x[y_q], y_tn[y_q]+shift, pch=8) #adds asterisks (flags)
  
} # end for protein loop
graphics.off()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# what is up with replicate 3?
df.test <- df %>%
  filter(R.Replicate==3)

df.test %>% pull(Temperature) %>% unique()
# replicate 3 only has measurements at three temperatures
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Qvalue 

qobj <- qvalue(p = df.statistics$Pvalue)
pi0 <- qobj$pi0

#adjust pvalues
padj <- p.adjust(p = df.statistics$Pvalue, method=c("BH"), n = length(df.statistics$Pvalue))
hist(padj, breaks=100)
unique(padj)
#------------------------------------------------------------------------------------


# END 
