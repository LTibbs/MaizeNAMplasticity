# Laura Tibbs Cortes
# Aug 21, 2020

# based on CERIS-JGRA code: https://github.com/jmyu/CERIS_JGRA
# as well as SimpleM https://github.com/LTibbs/SimpleM/
# and PET extraction code https://github.com/LTibbs/PET_extraction 

# set working directory with setwd()

#Start off printing the version information so you know what versions were used to generate the results
log_file <- paste('NAM_CERES.run.', format(Sys.time(), "%Y-%m-%d.%R:%S"), '.log',  sep = '');
sink(log_file);

#print start time:
print(paste0("Start time: ", Sys.time()))

# load libraries
library(R.utils)
library(rnoaa);
library(dplyr);
library(RCurl);
library(htmltab);
library(data.table);
# library(corrgram);
library(geosphere);
library(lubridate);
# library(animation);
library("colorspace");
library(RColorBrewer);
library(tidyverse)

source('sub_funcs_20200707_LTC.R');

print(sessionInfo())

# set working directories:
Dir <- ""
sp_dir <- ""
exp_dir <- "NAM_BLUE/"

# set filtering thresholds
line.outlier.filter <- 3 # remove lines before CERIS if there are LESS THAN 3; keep if there are AT LEAST 3 
filter.less.than <- 4 # remove data from predictions if there are LESS THAN X obs of the genotype for this trait, keep if there are AT LEAST (>=) X 
min.window.size <- 6 # set minimum window size: don't even look at windows less than X days when you subtract start from end day (windows are >= X). # Using window size with i < j-5, which gives min window size of 6, from Multicrop paper methods
last.FT <- 46 # set last day to check for FT traits (last day before any flower)
last.harvest <- last.FT + 60 # set last day to check for harvest traits
max.window <- 60 # don't consider any windows larger than 60 days 

# set up which traits must be decided by FT and which by harvest, based on when they were measured in original papers
FT.traits <- c("EH", "LL", "LW", "PH", "TPBN", "TL",
               "ULA", "ASI", "DTA", "DTS")
harvest.traits <- c("CD", "CL", "CM", "EM", "ERN",
                    "KN", "KPR", "TKW", "T20KW")

# Should we include precipitation-based indices in the search?
incl.precip="precip" # or "no.precip"

# Set up GDD calculations
## for maize 
t_base <- 50; t_max1 <- 86; t_max2 <- 1000; Haun_threshold <- 0;
#  this is for GDD. If the temp is below 50, it's set to 50 for GDD. If the temp is above 86, it's set to 86.
# Haun threshold is used in Barley etc for GDD and is also the reason for t_max2; not used here (see https://ndawn.ndsu.nodak.edu/help-barley-growing-degree-days.html)

# read in environmental metadata:
env_meta_file <- paste("11", 'Env_meta_table', sep = ''); 

env_meta_info_0 <- read.table(env_meta_file, header = T, sep = "\t", stringsAsFactors = F);

# set start and end years of the experiments
exp_s_year <- min(env_meta_info_0$TrialYear) 
exp_e_year <- max(env_meta_info_0$TrialYear) + 1

searching_daps <- 150 # how many days after planting do we care about?

# make name for ptt ptr file (holds env index values) out of : directory, number of environments, DAP to search up until
ptt_ptr_file <- paste(exp_dir, nrow(env_meta_info_0), 'Envs_envParas_DAP',  searching_daps, sep = '');

# read in the environmental data: 
# this was made from data from ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/by_year/ and https://earlywarning.usgs.gov/fews/product/81
# based on https://github.com/jmyu/CERIS_JGRA and https://github.com/LTibbs/PET_extraction
PTT_PTR <- read.table(ptt_ptr_file, header = T , sep = "\t"); 

# if want to exclude precip-based env variables:
if(incl.precip=="no.precip") {
  PTT_PTR <- PTT_PTR %>%
    select(-PRECIP, -PET, -H20.balance)
}

exp_traits_file <- paste(exp_dir, 'traits_ori_NAM', sep = ''); 
exp_fnd_file <- paste(exp_dir, 'traits_ori_FND', sep='')

pop.table <- fread("NAM.population.list.txt", sep="\t", data.table=F)


# Data pre-processing ---------------------------------------------------------

# remove outliers
  data.join <- fread() %>%
    filter(!is.na(ril_code)) %>%
    filter(ril_code != "") # remove entries that we don't know the genotype of
  # Now, using Aaron's Nature Plants paper on NAM as guide:
  # ". Outliers were removed as follows...
  # the interquartile ranges (IQRs) were calculated for each RIL
  # across environments and for each environment across RILs within a phenotype.
  # Any trait measurements of RILs that were more than 1.5 times larger or smaller
  # than either of the IQRs was removed."
  outlier.detect <- data.join %>%
    gather(key="Trait", value = "Value",
           -c(env_code, pop_code, Entry_id,
           ril_code, field, rep, pblock, block, Rangeposition, Rowposition)) %>%
    group_by(env_code, Trait) %>%
    mutate(Q1=quantile(Value, na.rm=T)[[2]],
           Q3=quantile(Value, na.rm=T)[[4]]) %>%
    distinct() %>%
    mutate(IQR=Q3-Q1) %>%
    mutate(lo.bound=Q1-1.5*IQR) %>%
    mutate(hi.bound=Q3+1.5*IQR) %>%
    mutate(to.keep.env=!(Value<lo.bound | Value > hi.bound))
  outlier.detect <- outlier.detect %>%
    ungroup() %>%
    group_by(ril_code, Trait) %>%
    mutate(ril.Q1=quantile(Value, na.rm=T)[[2]],
           ril.Q3=quantile(Value, na.rm=T)[[4]]) %>%
    distinct() %>%
    mutate(ril.IQR=ril.Q3-ril.Q1) %>%
    mutate(ril.lo.bound=ril.Q1-1.5*ril.IQR) %>%
    mutate(ril.hi.bound=ril.Q3+1.5*ril.IQR) %>%
    mutate(to.keep.ril=!(Value<ril.lo.bound | Value > ril.hi.bound)) %>%
    ungroup()
  outlier.detect <- outlier.detect %>% # keep observations that do NOT have to.keep==F
    filter(to.keep.env == T) %>% # (I'm ok with keeping to.keep==NA because those are just missing anyway)
    filter(to.keep.ril == T) %>%
    select(env_code, pop_code, Entry_id,
           field, rep, pblock, block, Rangeposition, Rowposition,
           ril_code, Trait, Value) %>%
    distinct()%>%
    pivot_wider(names_from=Trait, values_from=Value, values_fill=NA)

  # output data after outlier removal
  fwrite(outlier.detect, "traits_ori_NAM_no.outliers", sep="\t")
  
# Calculate BLUEs - - - - - - - - - -
# following Hung et al Heredity BUT
# excluding "env" and "family" from the model because
# I'm interested in those things and don't want to average them out!
library(lme4)
BLUE.data  <- fread("NAM/traits_ori_NAM_no.outliers") %>% 
  mutate(rep=paste(field, rep, sep="_"), # make rep etc unique within each field so no need to nest terms in model (for efficiency)
         pblock=paste(rep, pblock, sep="_"),
         block=paste(pblock, block, sep="_"),
         Rowposition=paste(rep, Rowposition, sep="_"),
         Rangeposition=paste(rep, Rangeposition, sep="_"))

# Ideal BLUE model: other things are just a simplification of this
# lmer(data=data[[n]], Value ~ Grp + Grp:ril_code +
#                             (1|field) + (1|field:rep) + (1|field:rep:pblock) +
#                             (1|field:rep:pblock:block) +(1|field:rep:Rowposition)+
#                             (1|field:rep:Rangeposition)

# make things factors:
BLUE.data$field <- factor(BLUE.data$field)
BLUE.data$rep <- factor(BLUE.data$rep)
BLUE.data$pblock <- factor(BLUE.data$pblock)
BLUE.data$ril_code <- factor(BLUE.data$ril_code)
BLUE.data$block <- factor(BLUE.data$block)
BLUE.data$Rangeposition <- factor(BLUE.data$Rangeposition)
BLUE.data$Rowposition <- factor(BLUE.data$Rowposition)

# make function to pull BLUEs (efficiently):
pheno.efficient <- function(predf, original_data, BLUEs.from.file=F) {
if(BLUEs.from.file) { # if reading in BLUEs from a file from a previous run, don't need to make the tibble
  df <- predf
} else {df <- tibble(Gen=names(predf), Estimate=predf)}
df$Gen <- gsub("ril_code", "", df$Gen) # remove extra characters

# pull intercept value
df.b1 <- df$Estimate[df$Gen=="(Intercept)"]

# calculate BLUEs:
df <- mutate(df, BLUE=ifelse(Gen=="(Intercept)", Estimate, # if intercept, pheno=estimate
                             Estimate + df.b1)) # otherwise, add estimate + intercept

# but the very last Gen (missing from named Gens) will have BLUE=intercept
thing <- original_data %>%
  filter(!is.na(Value))
stopifnot(length(setdiff(thing$ril_code, df$Gen))==1)
df$Gen <- gsub("\\(Intercept\\)", setdiff(thing$ril_code, df$Gen), df$Gen) # set the intercept as the missing genotype's value

return(df)

}

# Calculate the BLUEs:
R.BLUEs  <- vector("list", length(unique(BLUE.data$env_code))*length(c(FT.traits, harvest.traits)))
data <- vector("list", length(unique(BLUE.data$env_code))*length(c(FT.traits, harvest.traits)))

n <- 1
for(my.trait in c(FT.traits, harvest.traits)) {
# for(my.trait in c(FT.traits, harvest.traits)[1:2]) {
  print(paste0("beginning trait ", my.trait))
  print(Sys.time())

    for(my.env in unique(BLUE.data$env_code)) {
# for(my.env in unique(BLUE.data$env_code)[1:2]) { # make mini
  print(paste0("beginning env ", my.env))
  print(Sys.time())


  # filter out trait and env of interest
  data[[n]] <- BLUE.data %>%
    filter(env_code==my.env) %>%
    # filter(pblock %in% c("1_1_1", "1_1_2")) %>% # make mini
    select(env_code, pop_code, Entry_id, field, rep, pblock, block, Rangeposition, Rowposition, ril_code, Value=(my.trait))


  # skip if no data for trait in env
  if(nrow(data[[n]] %>% filter(!(is.na(Value))))==0) {print(paste("No data for", my.trait, "in", my.env))
    next}
  # comment the following out if only want to read in the fixef to feed into function:

  # Some envs have multiple fields, others multiple reps
  if(length(unique(data[[n]]$field))==1) {
    if(length(unique(data[[n]]$rep))==1) { # 1 rep, 1 field
      BLUE.model <- lmer(data=data[[n]], Value ~ ril_code +
                           (1|pblock) +
                           (1|block) +(1|Rowposition)+
                           (1|Rangeposition),
                         na.action = "na.exclude", control=lmerControl(calc.derivs=FALSE))
    } else { # >1 rep, 1 field
      BLUE.model <- lmer(data=data[[n]], Value ~ ril_code +
                           (1|rep) + (1|pblock) +
                           (1|block) +(1|Rowposition)+
                           (1|Rangeposition),
                         na.action = "na.exclude", control=lmerControl(calc.derivs=FALSE))
    }

  } else {
    # if(length(unique(data[[n]]$rep))==1) { # 1 rep, >1 field. It will look like >1 rep because "field" is now included in "rep" but it's really just one. No environments exist with really >1 field AND >1 rep
    BLUE.model <- lmer(data=data[[n]], Value ~ ril_code +
                         (1|field) + (1|pblock) +
                         (1|block) +(1|Rowposition)+
                         (1|Rangeposition),
                       na.action = "na.exclude", control=lmerControl(calc.derivs=FALSE))
    # } else { # full model when rep>1 AND field >1
    # warning(print("no environment has >1 rep AND >1 field!"))
    # BLUE.model <- lmer(data=data[[n]], Value ~ Grp + Grp:ril_code +
    #                      (1|field) + (1|field:rep) + (1|field:rep:pblock) +
    #                      (1|field:rep:pblock:block) +(1|field:rep:Rowposition)+
    #                      (1|field:rep:Rangeposition),
    #                    na.action = "na.exclude")
    # }
  }


  # pull BLUEs
  R.BLUEs.temp <- fixef(BLUE.model)
  # output in case this fails:
  fwrite(tibble(Gen=names(R.BLUEs.temp), Estimate=R.BLUEs.temp),
         paste0("BLUEs.", my.env, ".", my.trait, ".csv"))

  # read in BLUEs:
  R.BLUEs.temp <- fread(paste0("BLUEs.", my.env, ".", my.trait, ".csv"))

  R.BLUEs[[n]] <- pheno.efficient(predf=R.BLUEs.temp, original_data=data[[n]], BLUEs.from.file=T)
  R.BLUEs[[n]] <- R.BLUEs[[n]] %>%
    mutate(env_code=my.env) %>%
    mutate(Trait=my.trait)

  n <- n+1
}

}
out.BLUEs <- do.call("rbind", R.BLUEs) %>% select(-Estimate)# rbind resuls together
out.wide <- out.BLUEs %>%
  pivot_wider(names_from = Trait, values_from="BLUE") %>%
  rename(ril_code=Gen) %>%
  left_join(BLUE.data %>% select(env_code, pop_code, Entry_id, ril_code) %>% distinct,
          by = c("ril_code", "env_code"))%>%
  select(colnames(mean.data))

# output BLUEs
fwrite(out.wide, "NAM/traits_BLUE_full", sep="\t") # output ALL BLUEs, including diversity panel
fwrite(out.wide %>% filter(pop_code != 27 | ril_code %in% c("B73","HP301", pop.table$parent)), exp_traits_file, sep="\t") # output BLUEs for RILs and parents
fwrite(out.wide %>% filter(ril_code %in% c("B73","HP301", pop.table$parent)), exp_fnd_file, sep="\t") # output only parents


# Read in processed data --------------------------------------------------

    # read in the trait data:
    exp_traits <- read.table(exp_traits_file, sep = "\t", header = T, stringsAsFactors = F) %>% # read in multi-environment phenotype
      filter(pop_code !=17) %>% # remove IBM, it's not actually in NAM and has no genotype
      filter(ril_code != "Mo17") # remove IBM parent
    total_n <- length(unique(exp_traits$ril_code)) # find total number of lines
    
    fnd_traits <-  read.table(exp_fnd_file, sep="\t", header=T, stringsAsFactors = F)
    
    # every time -- make exp_traits match genotype
    mytaxa <- fread("NAM.GD.hidensity.matchY.txt", select=1)
    exp_traits <- exp_traits %>%
      filter(ril_code %in% mytaxa$Taxa)
    stopifnot(setdiff(mytaxa$Taxa, exp_traits$ril_code)==0)
    stopifnot(setdiff(exp_traits$ril_code,mytaxa$Taxa)==0)
    exp_traits$ril_code <- factor(exp_traits$ril_code,levels=mytaxa$Taxa)
    exp_traits <- exp_traits %>%
      arrange(ril_code)
    exp_traits$ril_code <- as.character(exp_traits$ril_code)
    stopifnot(all.equal(unique(exp_traits$ril_code), mytaxa$Taxa))
    
    all_env_codes <- sort(unique(exp_traits$env_code)); # pull all environment codes
    env_cols <- colorspace::rainbow_hcl(length(all_env_codes), c = 80, l = 60, start = 0, end = 300, fixup = TRUE, alpha = 0.75); # make colors
    
    # read in geno with ril names in order:
    taxa <- mytaxa$Taxa # for the folds from CERIS to work, taxa needs to contain EVERY taxa, not just those with phenotypes
    
    # read in data formatted for FR-gBLUP
    fr.exp_traits <- fread("NAM_BLUE/traits_ori_NAM_FR", data.table = F)
    
# CERES for whole NAM population: -----------------------------------------
    # based on CERIS-JGRA code: https://github.com/jmyu/CERIS_JGRA
    
    print(paste0("CERES for the whole population, with parameters: line.outlier.filter=", line.outlier.filter,
                 ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
                 ", last.harvest=", last.harvest, ", max.window.size=", max.window))
    output.results <- vector("list", ncol(exp_traits)-2)
    for (trait_num in c(5:ncol(exp_traits))) {
      # for (trait_num in c(23)) {
      trait <- colnames(exp_traits)[trait_num];
      current.data <- exp_traits[,trait_num,]
      exp_trait_dir <- paste(exp_dir, trait, '/', sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir)};
      exp_pop_dir <- paste0(exp_trait_dir, "whole_pop", "/") # make directory
      if(incl.precip=="no.precip") {exp_pop_dir <- paste0(exp_trait_dir, "whole_pop_no_precip", "/")} # make no precip directory
      if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir)};# make directory
      
      # set indices for which column name is which
      lInd <- which(colnames(exp_traits) == 'ril_code'); # assumint that each ril is essentially the "line" from before
      eInd <- which(colnames(exp_traits) == 'env_code');
      tInd <- which(colnames(exp_traits) == trait);
      
      exp_trait <- exp_traits[,c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
      colnames(exp_trait)[ncol(exp_trait)] <- 'Yobs'; # rename phenotype data column as Yobs
      
      # remove missing values--negative is ok though
      exp_trait <- exp_trait[!is.na(exp_trait$Yobs),];
      # exp_trait <- exp_trait[exp_trait$Yobs > 0,]
      
      # How many observations are there per pedigree?
      obs_lne_n <- aggregate(Yobs ~ ril_code, data = exp_trait, length);
      # hist(obs_lne_n$Yobs)
      
      # find lines with only one or two observations and remove them:
      line_outlier <- obs_lne_n$ril_code[obs_lne_n$Yobs < line.outlier.filter]
      exp_trait <- exp_trait[!(exp_trait$ril_code %in% line_outlier),]
      
      line_codes <- unique(exp_trait$ril_code); # pull unique line codes
      
      # find mean value for each environment:
      env_mean_trait_0 <- na.omit(aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code),
                                            mean, na.rm = T));
      colnames(env_mean_trait_0)[2] <- 'meanY';
      env_mean_trait <- env_mean_trait_0[order(env_mean_trait_0$meanY),]; # sort by mean phenotype value
      
      #### Do you want to leave one environment out while calculating correlations?
      LOO_env <- 0; ### for maize diversity panel, it is 0
      
      # set param for graphing:
      col_wdw <- 25;
      col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1) # this function makes a diverging color palette
      
      # Perform exhaustive search in WHOLE NAM pop for optimal window and environmental parameter to use
      Exhaustive_search_full(env_mean_trait, PTT_PTR, searching_daps, exp_pop_dir,
                             current.data,
                             trait,
                             searching_daps, searching_daps, LOO_env, min.window.size = min.window.size)
      
      # plot the exhaustive search triangles:
      pop_cor_file <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', LOO_env, 'LOO_cor_whole_pop', sep = '');
      pop_pval_file <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', LOO_env, 'LOO_P_whole_pop', sep = '');
      
      exhaustive_plot(pop_cor_file, pop_pval_file, PTT_PTR, searching_daps, searching_daps,
                      current.data, type="r.only")
      exhaustive_plot(pop_cor_file, pop_pval_file, PTT_PTR, searching_daps, searching_daps,
                      current.data, type="p.only")
      exhaustive_plot(pop_cor_file, pop_pval_file, PTT_PTR, searching_daps, searching_daps,
                      current.data, type="both")
      
      # filter windows by established criteria:
      search.results <- fread(paste(exp_pop_dir, trait, '_', nrow(env_mean_trait),
                                    'Envs_PTTPTR_', LOO_env, 'LOO_cor_whole_pop', sep = ''))%>%
        select(-starts_with("nR")) # don't need the negative version
      search.results <- search.results %>%
        tidyr::gather(key="Parameter", value="Corr", -Day_x, -Day_y, -window)
      
      search.results <- search.results %>%
        filter(window >= min.window.size) %>% # i < j-5
        filter(window <= max.window)
      
      if(trait %in% FT.traits) {
        search.results <- search.results %>%
          filter(Day_y <= last.FT) # last possible day for FT traits = 55
      } else if (trait %in% harvest.traits) {
        search.results <- search.results %>%
          filter(Day_y <= last.harvest) # last possible day for harvest traits
      } else (warning("trait not in FT or harvest traits"))
      
      # Only now choose the remaining window with the best corr:
      search.results  <- search.results %>%
        filter(abs(Corr)==max(abs(Corr), na.rm = T))
      ties <- nrow(search.results) # record if there are ties
      # save.search.results <- search.results # output all results
      fwrite(x = search.results, file=paste0(exp_pop_dir,trait, "_optimalparameters_", LOO_env, "LOO_",
                                             searching_daps, "DAP.txt"), sep="\t")
      # if there are ties, just pick the top one in the table.
      search.results <- search.results[1,]%>%
        mutate(trait=trait)
      output.results[[trait_num-2]] <- search.results
      print(output.results[[trait_num-2]])
      
      maxR_dap1 <- search.results$Day_x;
      maxR_dap2 <- search.results$Day_y;
      kPara_Name <- search.results$Parameter;
      kPara_Name <- gsub("R_", "", kPara_Name)
      PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); ## DL -> 5, GDD -> 6; PTT -> 7; PTR -> 8; PTD1 -> 8, PTD2 -> 9, PTS -> 10
      
      
      # need to SAVE the optimal window chosen, the associated correlation, and whether there were ties, for each trait.
      
      # Plot the environmental parameter means vs trait mean for each environment
      Plot_Trait_mean_envParas(env_mean_trait, PTT_PTR, maxR_dap1, maxR_dap2, trait, exp_pop_dir, env_cols);
      
      # Make output file containing the slopes and intercepts from the linear models for each line across environments
      Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, exp_trait, line_codes, exp_pop_dir, PTT_PTR_ind, filter.less.than = filter.less.than);
      
      #### LOOCV function
      # aka: leave one out cross validation:
      obs_prd_file <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, maxR_dap1, '_', maxR_dap2, sep = '');
      LOO_pdf_file <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, maxR_dap1, '_', maxR_dap2, '.png', sep = '');
      
      # Output data file with results of leave-one-out cross validation using both environ mean and param as predictors
      # if (!file.exists(obs_prd_file)) {
      if(length(unique(exp_trait$env_code))>=filter.less.than){ # only do this if there's enough environments to actually work
        prdM <- LOOCV(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, PTT_PTR_ind, exp_trait,  obs_prd_file, filter.less.than=filter.less.than, line_codes)
      }
      
      # Plot the observed vs predicted results using several prediction methods (see subfunction notes)
      Plot_prediction_result(obs_prd_file, all_env_code, prdM, kPara_Name, LOO_pdf_file, trait=trait);
    }
    
    # output the auto-chosen optimal parameter and window for each trait
    if (incl.precip=="no.precip") {
      fwrite(do.call("rbind", output.results), paste0(exp_dir, "CERES.auto.optimal.parameters.whole_pop.no_precip.txt"), sep="\t")
    } else {
      fwrite(do.call("rbind", output.results), paste0(exp_dir, "CERES.auto.optimal.parameters.whole_pop.txt"), sep="\t")
    }
    

#  FR-gBLUP 1->3 --------------------------------------------------
    # Goal: 1 ->3 prediction using windows from each training set
    # So:
    #  training set = randomly-selected 50% of genotypes from each family (all environments)
    #  testing set = the other 50% of genotypes 
    # 1. Get optimal window from CERIS, 
    # 2. Use asreml to get slope and intercept parameters for ALL lines (training and testing),
    #      which will give us the final pheno
    # 4. Repeat 50x
    
    source("genomic_prediction_functions.R")
    
    # This is designed for 2 folds-- split 50/50 training/testing
    nfolds <- 2
    rep.num <- c(1:10) # what rep to run -- eventually, want 1:50
    predict.type <- "within.pop" # set to within.pop or across.pop to predict either within each NAM pop or across the whole thing
    CERIS.done <- TRUE # use to determine if CERIS already run to prevent re-running it
    
    # make the folders for the .pop version (equivalent to .v4 in rrBLUP)
    exp_prepop_dir <- paste0(exp_dir, "FR.TEUG.1to3.pop/"); if (!dir.exists(exp_prepop_dir))  { dir.create(exp_prepop_dir, recursive=T)};
    exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop/"); if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir, recursive = T)};
    if(incl.precip=="no.precip") {exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop_no_precip", "/")} # make no precip directory
    if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir)};# make directory
    if(!dir.exists(paste0(exp_pop_dir, "output"))) {dir.create(paste0(exp_pop_dir, "output")) }
    
    #### Do you want to leave one environment out while calculating correlations?
    LOO_env <- 0; ### leave as 0 here because I already manually left one out
    
    # set param for graphing:
    col_wdw <- 25;
    col_palette <- colorspace::diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1) # this function makes a diverging color palette
    
    print(paste0("FR-gBLUP 1 to 3 for the whole population, rrBLUP within families, with parameters: line.outlier.filter=", line.outlier.filter,
                 ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
                 ", last.harvest=", last.harvest, ", max.window.size=", max.window))
    
    # First time only: make grm files (sparse matrix of lower triangle of kinship):
    # read in kinship:
    # kinship was calculated in GAPIT (maf filter 0.01) based on public genotype data from 
    # ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023.vcf.gz (https://datacommons.cyverse.org/browse/iplant/home/shared/panzea/genotypes/GBS/v27)
    preK <- fread("NAM.GAPIT.Kin.VanRaden.csv", 
                  data.table=F, stringsAsFactors = F) %>%
      rename(taxa=Taxa)
    for(p in unique(pop.table$number)) {
      if(p==17) {next} # pop 17 has no genotype data
    # prove it
        # temp <- exp_traits %>% filter(pop_code==17)
        # sum(temp$ril_code %in% taxa)
      pop.taxa <- unique(fr.exp_traits %>% filter(pop_code==p) %>% select(ril_code))
      popK <- preK[,-1][preK$taxa %in% pop.taxa$ril_code, preK$taxa %in% pop.taxa$ril_code] 
      stopifnot(dim(popK)==nrow(pop.taxa))

      sa.K <- matrix(ncol=3, nrow=(.5*nrow(popK)^2)+nrow(popK))
      m <- 1 # make iterator (row in new sparse matrix)
      for (i in 1:nrow(popK)) { # go through all of K
        for(j in 1:ncol(popK)) {
          if(isZero(popK[i,j])) {next}
          if(j>i) {next}
          # get the info for the sparse matrix:
          sa.K[m,1] <- i
          sa.K[m,2] <- j
          sa.K[m,3] <- popK[i,j]
          if(is.na(sa.K[m,3])) {break}
          m <- m+1
        }
        print(i)
      }
      sa.K <- sa.K[rowSums(is.na(sa.K)) != ncol(sa.K),]
      fwrite(sa.K,paste0("asreml.pop",p,".grm"), col.names =F, sep="\t")
      rm(sa.K, popK, pop.taxa)
    }
    
    # run CERIS first:
if(!CERIS.done) {
  # read in geno with ril names in order:
  geno <- fread(file="NAM.GD.hidensity.matchY.txt", select=1)
  taxa <- geno$Taxa
  rm(geno)
  
  
  lInd <- which(colnames(exp_traits) == 'ril_code'); 
  eInd <- which(colnames(exp_traits) == 'env_code');
  exp_prepop_dir <- paste0(exp_dir, "TEUG.1to3.v4/"); if (!dir.exists(exp_prepop_dir))  { dir.create(exp_prepop_dir, recursive=T)};
  exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop/"); if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir, recursive = T)};
  if(incl.precip=="no.precip") {exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop_no_precip", "/")} # make no precip directory
  if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir)};# make directory
  
  #### Do you want to leave one environment out while calculating correlations?
  LOO_env <- 0; ### leave as 0 here because I already manually left one out
  
  # set param for graphing:
  col_wdw <- 25;
  col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1) # this function makes a diverging color palette
  
  print(paste0("1 to 3 for the whole population, rrBLUP within families, with parameters: line.outlier.filter=", line.outlier.filter,
               ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
               ", last.harvest=", last.harvest, ", max.window.size=", max.window))
  
  for (i in rep.num) {
    # first, make folds:
    
    # make random folds for this rep, across all subpops:
    taxa.table <- as.tibble(taxa)
    colnames(taxa.table) <- c("value")
    taxa.table$obs.id <- 1:length(taxa)
    taxa.table <- left_join(taxa.table, exp_traits, by=c("value"="ril_code")) %>%
      select(value, obs.id, pop_code) %>%
      distinct
    #make the folds--
    # for any number of folds:
    myF <- vector("list", nfolds)
    # preF <- vector("list", nfolds)
    for(subpop in unique(exp_traits$pop_code)) {
      # get the obs id's for this pop:
      pop.ids <- taxa.table %>%
        filter(pop_code==subpop) %>%
        select(obs.id)
      
      preF <- fold.maker.2(pop.ids$obs.id, nfolds)
      for(j in c(1:nfolds)) {
        myF[[j]] <- c(myF[[j]], preF[[j]])
      }
    }
    
    output.results <- vector("list", 
                             length(c(FT.traits, harvest.traits))*nrow(pop.table))
    n <- 1  
    # get optimal windows and slope + intercept:
    for(trait in c(FT.traits, harvest.traits)) {
      tInd <- which(colnames(exp_traits) == trait);
      
      for(fold in c(1:nfolds)) {
        # pull training data:
        exp_trait_trn <- exp_traits[!(exp_traits$ril_code %in% taxa[myF[[fold]]]),c(lInd, eInd, tInd)]
        colnames(exp_trait_trn)[ncol(exp_trait_trn)] <- 'Yobs'; # rename phenotype data column as Yobs
        # and pull the training data as a vector only:
        current.data <- exp_trait_trn[,c(colnames(exp_trait_trn)=="Yobs")]
        
        # remove missing values
        exp_trait_trn <- exp_trait_trn[!is.na(exp_trait_trn$Yobs),];
        # How many observations are there per pedigree?
        obs_lne_n <- aggregate(Yobs ~ ril_code, data = exp_trait_trn, length);
        # find lines with only one or two observations and remove them:
        line_outlier <- obs_lne_n$ril_code[obs_lne_n$Yobs < line.outlier.filter]
        exp_trait_trn <- exp_trait_trn[!(exp_trait_trn$ril_code %in% line_outlier),]
        line_codes <- unique(exp_trait_trn$ril_code); # pull unique line codes
        # also pull testing data:
        exp_trait_test <- exp_traits[(exp_traits$ril_code %in% taxa[myF[[fold]]]),c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
        exp_trait_test <- exp_trait_test[!(exp_trait_test$ril_code %in% line_outlier),] # also remove the outlier lines from your testing set to not mess up the training means
        colnames(exp_trait_test)[ncol(exp_trait_test)] <- 'Yobs'; # rename phenotype data column as Yobs
        
        
        # find mean value for each environment:
        env_mean_trait_0_trn <- na.omit(aggregate(x = exp_trait_trn$Yobs, by = list(env_code = exp_trait_trn$env_code),
                                                  mean, na.rm = T));
        colnames(env_mean_trait_0_trn)[2] <- 'meanY';
        env_mean_trait_trn <- env_mean_trait_0_trn[order(env_mean_trait_0_trn$meanY),]; # sort by mean phenotype value
        
        env_mean_trait_0_test <- na.omit(aggregate(x = exp_trait_test$Yobs, by = list(env_code = exp_trait_test$env_code),
                                                   mean, na.rm = T));
        colnames(env_mean_trait_0_test)[2] <- 'meanY';
        env_mean_trait_test <- env_mean_trait_0_test
        
        # Perform exhaustive search in training set from WHOLE NAM pop, with 1 environment left out, for optimal window and environmental parameter to use
        Exhaustive_search_full(env_mean_trait_trn, PTT_PTR, searching_daps, exp_pop_dir,
                               current.data,
                               trait,
                               searching_daps, searching_daps, LOO_env, min.window.size = min.window.size)
        
        # need to read in default output and rename them:
        pop_cor_file.old <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_cor_whole_pop', sep = '');
        pop_pval_file.old <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_P_whole_pop', sep = '');
        
        pop_cor <- fread(pop_cor_file.old)
        pop_pval <- fread(pop_pval_file.old) 
        
        pop_cor_file <- paste(exp_pop_dir, trait, '_TEUG.1to3.fold', fold, ".rep", i, '_cor_whole_pop', sep = '')
        fwrite(pop_cor, pop_cor_file, sep="\t")
        
        pop_pval_file <-paste(exp_pop_dir, trait, '_TEUG.1to3.fold.', fold, ".rep", i, '_P_whole_pop', sep = '')
        fwrite(pop_pval, pop_pval_file, sep="\t")
        
        # tidy up:
        file.remove(pop_cor_file.old, pop_pval_file.old) 
        rm(pop_cor, pop_pval)
        
        # filter windows by established criteria:
        search.results <- fread(pop_cor_file)%>%
          select(-starts_with("nR")) # don't need the negative version
        search.results <- search.results %>%
          tidyr::gather(key="Parameter", value="Corr", -Day_x, -Day_y, -window)
        
        search.results <- search.results %>%
          filter(window >= min.window.size) %>% # i < j-5
          filter(window <= max.window)
        
        if(trait %in% FT.traits) {
          search.results <- search.results %>%
            filter(Day_y <= last.FT) # last possible day for FT traits = 55
        } else if (trait %in% harvest.traits) {
          search.results <- search.results %>%
            filter(Day_y <= last.harvest) # last possible day for harvest traits
        } else (warning("trait not in FT or harvest traits"))
        
        # Only now choose the remaining window with the best corr:
        search.results <- search.results <- search.results %>%
          filter(abs(Corr)==max(abs(Corr), na.rm = T))
        
        # if there are ties, just pick the top one in the table.
        search.results <- search.results[1,]%>%
          mutate(trait=trait)
        
        # output results:
        output.results[[n]] <- cbind(rep=i, fold=fold, search.results)
        print(output.results[[n]])
        
        maxR_dap1 <- search.results$Day_x;
        maxR_dap2 <- search.results$Day_y;
        kPara_Name <- search.results$Parameter;
        kPara_Name <- gsub("R_", "", kPara_Name)
        PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); ## DL -> 5, GDD -> 6; PTT -> 7; PTR -> 8; PTD1 -> 8, PTD2 -> 9, PTS -> 10
        
        # Make output file containing the slopes and intercepts from the linear models for each line across environments
        Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait_trn, PTT_PTR, exp_trait_trn, line_codes, exp_pop_dir, PTT_PTR_ind, 
                        filter.less.than); # do filtering for number of observations here
        
        # read file in and rename it:
        slope_file.old <-paste0(exp_pop_dir,"Intcp_Slope")
        pop_slope <- fread(slope_file.old)
        slope_file <- paste0(exp_pop_dir, trait, "_TEUG.1to3.fold", fold, ".rep", i, "_Intcp_Slope")
        fwrite(pop_slope,slope_file, sep="\t")
        
        # tidy up:
        file.remove(slope_file.old) 
        rm(pop_slope)
        rm(exp_trait_trn, exp_trait_test)
        
        n <- n+1 # increase iterator
        
      }
    }
    # output windows:
    opt.windows <- do.call("rbind", output.results)
    fwrite(opt.windows, paste0(exp_pop_dir, "TEUG.1to3.auto.optimal.parameters.whole_pop.rep", i, ".csv"))
    
  } 
  
}
    
    # Now prep for and run FR-gBLUP in ASREML
    
    # make the folders for the .pop version (equivalent to .v4 in rrBLUP)
    exp_prepop_dir <- paste0(exp_dir, "FR.TEUG.1to3.pop/"); if (!dir.exists(exp_prepop_dir))  { dir.create(exp_prepop_dir, recursive=T)};
    exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop/"); if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir, recursive = T)};
    if(incl.precip=="no.precip") {exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop_no_precip", "/")} # make no precip directory
    if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir)};# make directory
    if(!dir.exists(paste0(exp_pop_dir, "output"))) {dir.create(paste0(exp_pop_dir, "output")) }
    
    for (i in rep.num){
      if(CERIS.done) { # if CERIS-JGRA already run for 1->3, need to read in: optimal window, which line in which fold
        # first, read in folds:
        myF <-vector("list", nfolds)
        for (j in c(1:nfolds)) {
          fold.long <- vector()
          for(my.trait in c(FT.traits, harvest.traits)){ # need to loop through traits because some had missing data
            fold <- fread(paste0(exp_dir, "TEUG.1to3.v4/whole_pop/",my.trait,"_TEUG.1to3.fold",j,".rep", i, "_Intcp_Slope")) %>%
              select(line_codes)
            fold.long <- unique(c(fold.long, fold$line_codes)) 
          }
          myF[[j]] <- which(taxa %in% fold.long)
        }
        stopifnot(sum(taxa[setdiff(1:length(taxa), c(myF[[1]], myF[[2]]))] %in% exp_traits$ril_code)==0) # check that everything with a phenotype is in here
        stopifnot(length(intersect(myF[[1]], myF[[2]]))==0) # check that no lines are in both folds
        
        
      }  
      # Now, read in optimal windows:
      opt.windows <- fread(paste0(exp_dir, "TEUG.1to3.v4/whole_pop/TEUG.1to3.auto.optimal.parameters.whole_pop.rep", i, ".csv"))
      if(incl.precip=="no.precip") {    
        opt.windows <- fread(paste0(exp_dir, "TEUG.1to3.v4/whole_pop_no_precip/TEUG.1to3.auto.optimal.parameters.whole_pop.rep", i, ".csv"))
      }
      
      # prep the input files for asreml:
      
      # get the EC (environmental covariate) value for each trait x fold combination:
      EC.table <- tibble(rep=NA, fold=NA, trait=NA, env_code=NA, EC=NA)
      for(my.trait in c(FT.traits, harvest.traits)) {
        for(my.fold in 1:nfolds) {
          current.window <- opt.windows %>% filter(fold==my.fold, trait==my.trait)
          maxR_dap1 <- current.window$Day_x;
          maxR_dap2 <- current.window$Day_y;
          kPara_Name <- current.window$Parameter;
          kPara_Name <- gsub("R_", "", kPara_Name)
          
          current.EC <- PTT_PTR %>% 
            group_by(env_code) %>% 
            mutate(day.num=1:150) %>% 
            filter(day.num %in% c(maxR_dap1:maxR_dap2)) %>%
            select(env_code, kPara_Name, day.num) %>%
            mutate(EC=mean(get(kPara_Name))) %>%
            mutate(trait=my.trait,
                   rep=i,
                   fold=my.fold
            ) %>%
            select(rep, fold, trait, env_code, EC) %>%distinct
          current.EC$EC <- scale(current.EC$EC,center = T,scale = T) # scale the environmental covariate as required for FR-gBLUP
          
          EC.table <- rbind(EC.table, current.EC)
          rm(current.EC)
        }
      }
      EC.table <- EC.table %>% filter(!is.na(rep))
      wide.EC <- EC.table %>%
        group_by(trait) %>%
        mutate(row = row_number()) %>%
        pivot_wider(names_from=trait, values_from=EC) %>%
        select(-row) %>%
        rename_with(.cols=-c(rep, fold, env_code), ~ paste0(.x, ".EC"))
      
      # add the EC to the phenotype data:
      
      # make asreml input files, with different training sets set to NAs:
        for(my.pop in unique(fr.exp_traits$pop_code)) {
          if(is.na(my.pop)) {next}
          for(my.fold in 1:nfolds) {
            pheno.fold <- fr.exp_traits %>%
              rename(Genotype=ril_code) %>%
              filter(pop_code ==my.pop) %>%
              mutate(
                ASI=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, ASI),
                CD=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, CD),
                CL=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, CL),
                CM=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, CM),
                DTA=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, DTA),
                DTS=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, DTS),
                EH=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, EH),
                EM=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, EM),
                ERN=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, ERN),
                KN=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, KN),
                KPR=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, KPR),
                TKW=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, TKW),
                LL=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, LL),
                LW=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, LW),
                PH=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, PH),
                T20KW=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, T20KW),
                TPBN=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, TPBN),
                TL=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, TL),
                ULA=ifelse(Genotype %in% taxa[myF[[my.fold]]], NA, ULA)) %>%
              mutate(rep=i) %>%
              mutate(fold=my.fold) %>%
              select(Genotype, env_code, pop_code, rep, fold, everything(), -Entry_id)
            for(my.trait in c(FT.traits, harvest.traits)) { # filter out lines without enough observations
              temp.pheno.fold <- pheno.fold[,c(1,which(colnames(pheno.fold)==my.trait))]
              colnames(temp.pheno.fold)[2] <- "Yobs" 
              obs_lne_n <- aggregate(Yobs~Genotype, data=temp.pheno.fold, length)
              line_outlier <- obs_lne_n$Genotype[obs_lne_n$Yobs < (filter.less.than)] # want to make sure have AT LEAST filter.less.than-1 training obs, so need at least filter.less.than total obs
              
              temp.pheno.fold$Yobs[temp.pheno.fold$Genotype %in% line_outlier] <- NA
              pheno.fold[,which(colnames(pheno.fold)==my.trait)] <- temp.pheno.fold$Yobs
            }
            pheno.fold <- left_join(pheno.fold, wide.EC)
            colnames(pheno.fold) <- gsub("\\.", "_", colnames(pheno.fold))
            
            fwrite(pheno.fold, paste0(exp_pop_dir, paste0("asreml_pop",my.pop, "_fold", my.fold, "_rep", i, ".csv")))
            rm(pheno.fold)
          }
          
        }

    }
        
    # Run ASREML:
    # ASREML code is:
    # ../../../asreml.pop$J.grm !NSD 
    # asreml_pop$J_fold$K_rep$1.csv !SKIP 1 !MAXIT 10000 # should be 10000 # run job as: -R asreml_pop_rep.as repnum
    # !ROWFAC env_code !COLFAC Genotype
    # # !AISINGULARITIES
    # $I !SIGMAP ~ mu $I_EC mv !r  str(Genotype Genotype.$I_EC us(2).grm(Genotype))    
    # residual idv(env_code).id(Genotype)
    
    # Compile ASREML results - - -
    
      n <- 1
      result.id <- vector()
      for(my.trait in c(FT.traits, harvest.traits)){
        for(my.pop in unique(fr.exp_traits$pop_code)) {
          result.id[n] <- paste(my.trait, my.pop, sep="_")
          print(paste(my.trait, my.pop, sep=";"))
          n <- n+1
        }
        
      }

      sln <- fread(paste0(temp.dir, "asreml_full.pop.sln"), data.table=F) %>%
        mutate(Set = result.id[cumsum(Model_Term == "Model_Term")+1]) %>% # set the cycle name used for these table rows
        filter(Model_Term != "Model_Term") %>%
        mutate(tempcol=str_split(Set, '_')) %>% # split the string
        rowwise() %>%
        mutate(trait=unlist(tempcol)[1], pop=unlist(tempcol)[2]) %>%
        dplyr::select(-tempcol) %>%
        ungroup() 
      
      thing <- sln %>%
        filter(str_detect(Model_Term, '_EC')) %>%
        filter(!(str_detect(Model_Term, "Genotype"))) %>%
        mutate(tempcol=str_split(Model_Term, '_')) %>% # split the string
        rowwise() %>%
        mutate(check.trait=unlist(tempcol)[1]) %>%
        dplyr::select(-tempcol) %>%
        ungroup()
      stopifnot(all.equal(thing$trait, thing$check.trait))
      
      sln.full <- sln
      sln.full$pop <- as.numeric(sln.full$pop)
      sln.full$Effect <- as.numeric(sln.full$Effect)
      
      # calculate predictions:
      for (my.trait in c(FT.traits, harvest.traits)){
        pre.combine <- vector("list", length(unique(fr.exp_traits$pop_code)))
        n <- 1
        
        for (my.pop in sort(unique(fr.exp_traits$pop_code))){
          if (is.na(my.pop)) {next}
          
          # pull current data from sln file
          current.sln <- sln.full %>%
            filter(Set %in% sln.full$Set[which(sln.full$Model_Term==paste0(my.trait,"_EC"))]) %>%
            filter(pop==my.pop) %>%
            distinct
          if(nrow(current.sln)==0) {next} # skip if no data
          
          # pull input data for this pop
          current.pop <- fread(paste0(temp.dir, "asreml_pop",my.pop,".csv"), data.table=F) %>%
            dplyr::select(Genotype, env_code, paste0(my.trait), paste0(my.trait, "_EC"))%>%
            # filter(is.na(UQ(rlang::sym(my.trait)))) %>% # commenting this out because need to predict ALL of them, not just the unknowns
            # filter(env_code==my.env) %>% # only need the testing environ
            distinct
          
          current.pop$y.hat <- rep(NA, nrow(current.pop))
          current.pop$intcp.hat <- rep(NA, nrow(current.pop))
          current.pop$slope.hat <- rep(NA, nrow(current.pop))
          
          
          for (j in 1:nrow(current.pop)) {
            # calculate the prediction for each line
            current.genotype <- current.pop$Genotype[j]
            current.env <- current.pop$env_code[j]
            current.EC <- current.pop[j,4][[1]] # this only works because I select the column order above!
            current.intcp <- current.sln$Effect[current.sln$Level==current.genotype]
            current.slope <- current.sln$Effect[current.sln$Level==paste0(current.genotype, ".001")]
            
            current.pop$y.hat[j] <- current.sln$Effect[current.sln$Model_Term=="mu"]+ #mu
              current.sln$Effect[current.sln$Model_Term==paste0(my.trait, "_EC")]*current.EC + # + BXe
              current.intcp + # + uG
              current.slope*current.EC # + bG*xE
            current.pop$intcp.hat[j] <- current.intcp
            current.pop$slope.hat[j] <- current.slope
          }
          pre.combine[[n]] <- current.pop
          print(paste(my.trait, " pop ", my.pop, ": cor ", cor(pre.combine[[n]][,7], pre.combine[[n]]$y.hat, use="complete")))
          
          rm(current.pop)
          n <- n+1
          
        }
        
        combined <- do.call(what = "rbind", pre.combine) %>%
          mutate(obs.id=NA) %>%
          select(env_code, ril_code=Genotype, paste0(my.trait), obs.id, intcp.hat, slope.hat, y.hat)
        # plot(combined$PH, combined$y.hat)
        # cor(combined$PH, combined$y.hat, use="complete")
        
        kPara_Name <- "Env Param"
        trait <- my.trait
        
        Plot_TEUG_result(obs_prd_file = combined, all_env_codes, kPara_Name, trait=trait, 
                         forecast_png_file = paste0(temp.dir, "output/FR_", trait, "_wholepop_obs.prd.png"),
                         path = F, save.output = T)
        fwrite(combined, paste0(temp.dir, "output/FR_", trait, "_wholepop_obs.prd.txt"))
        
        
      }
    
      
      # FR-gBLUP 1->2 -----------------------------------------------------------
      
      # training set = n-1 environments (all genotypes)
      # testing set = the other environment (all genotypes, predict families one at a time)
      
      # use training set to find optimal window and to get slope and intercept
      # parameters for each genotype in CERIS. Then, use that slope and intercept, plus
      # the environmental index of the other environment, to get predicted pheno of
      # testing set
      
      # Repeat until each environment has been in the testing set
      #set current trait: 
      trait <- "DTA" # provide desired trait for CERIS; can also make loop
      
      # run CERIS to get optimal windows
      temp.dir <- paste0(exp_dir, "UETG.1to2/")
      if (!dir.exists(temp.dir))  { dir.create(temp.dir)};
      exp_trait_dir <- paste0(temp.dir, trait, '/', sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir)};
      exp_pop_dir <- paste0(exp_trait_dir, "whole_pop", "/"); 
      if(incl.precip=="no.precip") {exp_pop_dir <- paste0(exp_trait_dir, "whole_pop_no_precip", "/")} # make no precip directory
      if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir)};# make directory
      
      print(paste0("1 to 2 for ", trait, " for the whole population, with parameters: line.outlier.filter=", line.outlier.filter,
                   ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
                   ", last.harvest=", last.harvest, ", max.window.size=", max.window))
      
      # set indices for which column name is which
      lInd <- which(colnames(exp_traits) == 'ril_code'); # assumint that each ril is essentially the "line" from before
      eInd <- which(colnames(exp_traits) == 'env_code');
      tInd <- which(colnames(exp_traits) == trait);
      
      output.results <- vector("list", length(all_env_codes))
      UETG.results <- vector("list", length(all_env_codes))
      for (env_num in seq_along(all_env_codes)) {
        # for (env_num in 1:2) {
        current.env <- all_env_codes[env_num] # set current environment
        
        # skip if testing data is blank:
        if (sum(!(is.na(exp_traits[exp_traits$env_code==current.env,c(tInd)])))==0) {next}
        
        # pull training data
        exp_trait_trn <- exp_traits[exp_traits$env_code!=current.env,c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
        colnames(exp_trait_trn)[ncol(exp_trait_trn)] <- 'Yobs'; # rename phenotype data column as Yobs
        
        
        # remove missing values--negative is ok though
        exp_trait_trn <- exp_trait_trn[!is.na(exp_trait_trn$Yobs),];
        # exp_trait <- exp_trait[exp_trait$Yobs > 0,]
        
        # How many observations are there per pedigree?
        obs_lne_n <- aggregate(Yobs ~ ril_code, data = exp_trait_trn, length);
        # hist(obs_lne_n$Yobs)
        
        # find lines with only one or two observations and remove them:
        line_outlier <- obs_lne_n$ril_code[obs_lne_n$Yobs < line.outlier.filter]
        exp_trait_trn <- exp_trait_trn[!(exp_trait_trn$ril_code %in% line_outlier),]
        
        line_codes <- unique(exp_trait_trn$ril_code); # pull unique line codes
        
        # also pull testing data:
        exp_trait_test <- exp_traits[exp_traits$env_code==current.env,c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
        exp_trait_test <- exp_trait_test[!(exp_trait_test$ril_code %in% line_outlier),] # also remove the outlier lines from your testing set to not mess up the training means
        colnames(exp_trait_test)[ncol(exp_trait_test)] <- 'Yobs'; # rename phenotype data column as Yobs
        
        # and pull the training data as a vector only:
        current.data <- exp_trait_trn[,c(colnames(exp_trait_trn)=="Yobs")]
        
        # find mean value for each environment:
        env_mean_trait_0_trn <- na.omit(aggregate(x = exp_trait_trn$Yobs, by = list(env_code = exp_trait_trn$env_code),
                                                  mean, na.rm = T));
        colnames(env_mean_trait_0_trn)[2] <- 'meanY';
        env_mean_trait_trn <- env_mean_trait_0_trn[order(env_mean_trait_0_trn$meanY),]; # sort by mean phenotype value
        
        env_mean_trait_0_test <- na.omit(aggregate(x = exp_trait_test$Yobs, by = list(env_code = exp_trait_test$env_code),
                                                   mean, na.rm = T));
        colnames(env_mean_trait_0_test)[2] <- 'meanY';
        env_mean_trait_test <- env_mean_trait_0_test
        
        #### Do you want to leave one environment out while calculating correlations?
        LOO_env <- 0; ### leave as 0 here because I already manually left one out
        
        # set param for graphing:
        col_wdw <- 25;
        col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1) # this function makes a diverging color palette
        
        # Perform exhaustive search in WHOLE NAM pop, with 1 environment left out, for optimal window and environmental parameter to use
        Exhaustive_search_full(env_mean_trait_trn, PTT_PTR, searching_daps, exp_pop_dir,
                               current.data,
                               trait,
                               searching_daps, searching_daps, LOO_env, min.window.size = min.window.size)
        
        # need to read in default output and rename them:
        pop_cor_file.old <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_cor_whole_pop', sep = '');
        pop_pval_file.old <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_P_whole_pop', sep = '');
        
        pop_cor <- fread(pop_cor_file.old)
        pop_pval <- fread(pop_pval_file.old) 
        
        pop_cor_file <- paste(exp_pop_dir, trait, '_UETG.1to2.whole.no', current.env, '_cor_whole_pop', sep = '')
        fwrite(pop_cor, pop_cor_file, sep="\t")
        
        pop_pval_file <- paste(exp_pop_dir, trait, '_UETG.1to2.whole.no', current.env, '_P_whole_pop', sep = '')
        fwrite(pop_pval, pop_pval_file, sep="\t")
        
        # tidy up:
        file.remove(pop_cor_file.old, pop_pval_file.old) 
        rm(pop_cor, pop_pval)
        
        # filter windows by established criteria:
        search.results <- fread(pop_cor_file)%>%
          select(-starts_with("nR")) # don't need the negative version
        search.results <- search.results %>%
          tidyr::gather(key="Parameter", value="Corr", -Day_x, -Day_y, -window)
        
        search.results <- search.results %>%
          filter(window >= min.window.size) %>% # i < j-5
          filter(window <= max.window)
        
        if(trait %in% FT.traits) {
          search.results <- search.results %>%
            filter(Day_y <= last.FT) # last possible day for FT traits = 55
        } else if (trait %in% harvest.traits) {
          search.results <- search.results %>%
            filter(Day_y <= last.harvest) # last possible day for harvest traits
        } else (warning("trait not in FT or harvest traits"))
        
        # Only now choose the remaining window with the best corr:
        search.results <- search.results <- search.results %>%
          filter(abs(Corr)==max(abs(Corr), na.rm = T))
        
        # if there are ties, just pick the top one in the table.
        search.results <- search.results[1,]%>%
          mutate(trait=trait)
        output.results[[env_num]] <- cbind(testing.env=current.env, search.results)
        print(output.results[[env_num]])
        
        maxR_dap1 <- search.results$Day_x;
        maxR_dap2 <- search.results$Day_y;
        kPara_Name <- search.results$Parameter;
        kPara_Name <- gsub("R_", "", kPara_Name)
        PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); ## DL -> 5, GDD -> 6; PTT -> 7; PTR -> 8; PTD1 -> 8, PTD2 -> 9, PTS -> 10
        
        # Make output file containing the slopes and intercepts from the linear models for each line across environments
        Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait_trn, PTT_PTR, 
                        exp_trait_trn, line_codes, exp_pop_dir, PTT_PTR_ind,filter.less.than-1);
        
        
        # read file in and rename it:
        slope_file.old <-paste0(exp_pop_dir,"Intcp_Slope")
        pop_slope <- fread(slope_file.old)
        slope_file <- paste0(exp_pop_dir, trait, "_UETG.1to2.whole.no", current.env, "_Intcp_Slope")
        fwrite(pop_slope,slope_file, sep="\t")
        
        # tidy up:
        file.remove(slope_file.old) 
        rm(pop_slope)
        
        exp_trait <- rbind(exp_trait_trn, exp_trait_test) # combine training and testing data
        
        if(max(obs_lne_n$Yobs, na.rm = T)>= (filter.less.than-1) & sum(!(is.na(exp_trait_test$Yobs)))>=1) { # only run if enough obs:
          # Use slopes and intercepts to predict new environment:
          function.output <- UETG.1to2.function(maxR_dap1, maxR_dap2, PTT_PTR_ind, env_mean_trait_trn, env_mean_trait_test,
                                                PTT_PTR, exp_trait, test_env=current.env, filter.less.than)
          # and save the output:
          prdM <- function.output[[2]]
          UETG.results[[env_num]] <- function.output[[1]]
        }
        
        rm(search.results, exp_trait_test, exp_trait_trn, kPara_Name)
      }
      
      # Now, run ASREML
      
      CERIS.done <- TRUE # use to determine if CERIS already run before FRgBLUP
      
      temp.dir <- paste(exp_dir, "FR.UETG.1to2.pop", sep=""); if (!dir.exists(temp.dir))  { dir.create(temp.dir)};
      if(incl.precip=="no.precip") {
        temp.dir <- paste(exp_dir, "FR.UETG.1to2.pop_no_precip", sep=""); if (!dir.exists(temp.dir))  { dir.create(temp.dir)};
      }
      
      print(paste0("FR-gBLUP 1 to 2 for the whole population, with parameters: line.outlier.filter=", line.outlier.filter,
                   ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
                   ", last.harvest=", last.harvest, ", max.window.size=", max.window))
      
      #read in optimal windows:
      window.list <- vector("list", length=length(c(harvest.traits, FT.traits)))
      for (j in 1:length(c(harvest.traits, FT.traits))) {
        my.trait <- c(harvest.traits, FT.traits)[j]
        if(incl.precip=="no.precip") {
          window.list[[j]] <- fread(paste0(exp_dir, "UETG.1to2/", my.trait, "/whole_pop_no_precip/", my.trait, "_UETG.1to2.auto.optimal.parameters.whole_pop.txt"))
        } else{
          window.list[[j]] <- fread(paste0(exp_dir, "UETG.1to2/", my.trait, "/whole_pop/", my.trait, "_UETG.1to2.auto.optimal.parameters.whole_pop.txt"))
        }
      }
      opt.windows <- do.call("rbind", window.list)
      
      # prep the input files for asreml:
      
      # get the EC (environmental covariate) value for each trait x left-out env combination:
      EC.table <- tibble(trait=NA, env_code=NA, testing.env=NA, EC=NA)
      
      for(my.trait in c(FT.traits, harvest.traits)) {
        for(my.env in env_meta_info_0$env_code) {
          current.window <- opt.windows %>% filter(testing.env==my.env, trait==my.trait)
          if(nrow(current.window)==0) {
            print(paste0("No window found for ", my.trait, " in training env ", my.env))
            next
          }
          
          maxR_dap1 <- current.window$Day_x;
          maxR_dap2 <- current.window$Day_y;
          kPara_Name <- current.window$Parameter;
          kPara_Name <- gsub("R_", "", kPara_Name)
          
          current.EC <- PTT_PTR %>% 
            group_by(env_code) %>% 
            mutate(day.num=1:150) %>% 
            filter(day.num %in% c(maxR_dap1:maxR_dap2)) %>%
            select(env_code, kPara_Name, day.num) %>%
            mutate(EC=mean(get(kPara_Name))) %>%
            mutate(trait=my.trait,
                   testing.env=my.env) %>%
            select(trait, env_code, testing.env, EC) %>%distinct
          current.EC$EC <- scale(current.EC$EC,center = T,scale = T) # scale the environmental covariate as required for FR-gBLUP
          
          EC.table <- rbind(EC.table, current.EC)
          rm(current.EC)
        }
      }
      EC.table <- EC.table %>% filter(!is.na(trait))
      wide.EC <- EC.table %>%
        pivot_wider(names_from=trait, values_from=EC)%>%
        rename_with(.cols=-c(env_code, testing.env), ~ paste0(.x, ".EC"))
      
      # add EC to the phenotype data:
      pheno.full <- left_join(fr.exp_traits, wide.EC, by="env_code") %>%
        select(-Entry_id) %>%
        arrange(env_code, ril_code)
      
      # make asreml input files, with different training sets set to NAs:
      
      for(my.env in env_meta_info_0$env_code) {
        for(my.pop in unique(fr.exp_traits$pop_code)){
          pheno.part <- pheno.full %>%
            filter(testing.env==my.env) %>%
            filter(pop_code==my.pop) %>%
            mutate(
              ASI=ifelse(env_code==my.env, NA, ASI),
              CD=ifelse(env_code==my.env, NA, CD),
              CL=ifelse(env_code==my.env, NA, CL),
              CM=ifelse(env_code==my.env, NA, CM),
              DTA=ifelse(env_code==my.env, NA, DTA),
              DTS=ifelse(env_code==my.env, NA, DTS),
              EH=ifelse(env_code==my.env, NA, EH),
              EM=ifelse(env_code==my.env, NA, EM),
              ERN=ifelse(env_code==my.env, NA, ERN),
              KN=ifelse(env_code==my.env, NA, KN),
              KPR=ifelse(env_code==my.env, NA, KPR),
              TKW=ifelse(env_code==my.env, NA, TKW),
              LL=ifelse(env_code==my.env, NA, LL),
              LW=ifelse(env_code==my.env, NA, LW),
              PH=ifelse(env_code==my.env, NA, PH),
              T20KW=ifelse(env_code==my.env, NA, T20KW),
              TPBN=ifelse(env_code==my.env, NA, TPBN),
              TL=ifelse(env_code==my.env, NA, TL),
              ULA=ifelse(env_code==my.env, NA, ULA)) %>%
            select(Genotype=ril_code, env_code, pop_code, testing.env, everything()) 
          for(my.trait in c(FT.traits, harvest.traits)) { # filter out lines without enough observations
            temp.pheno.part <- pheno.part[,c(1,which(colnames(pheno.part)==my.trait))]
            colnames(temp.pheno.part)[2] <- "Yobs" 
            obs_lne_n <- aggregate(Yobs~Genotype, data=temp.pheno.part, length)
            line_outlier <- obs_lne_n$Genotype[obs_lne_n$Yobs < (filter.less.than)] # want to make sure have AT LEAST filter.less.than-1 training obs, so need at least filter.less.than total obs
            
            temp.pheno.part$Yobs[temp.pheno.part$Genotype %in% line_outlier] <- NA
            pheno.part[,which(colnames(pheno.part)==my.trait)] <- temp.pheno.part$Yobs
          }
          
          colnames(pheno.part) <- gsub("\\.", "_", colnames(pheno.part))
          
          fwrite(pheno.part, paste0(temp.dir, paste0("/asreml_1to2_no",my.env,"_pop", my.pop, ".csv")))
          rm(pheno.part)
        }
        
      }
      
      
      # compile FRgBLUP1->2 --- 
      
      # get the CYCLE arguments for asreml_rep.as:
      n <- 1
      result.id <- vector()
      for(my.trait in c(FT.traits, harvest.traits)) {
        for(my.env in env_meta_info_0$env_code) {
          current.window <- opt.windows %>% filter(testing.env==my.env, trait==my.trait)
          if(nrow(current.window)==0) {
            # print(paste0("No window found for ", my.trait, " in training env ", my.env))
            next
          }
          for(my.pop in unique(fr.exp_traits$pop_code)) {
            result.id[n] <- paste(my.trait, my.env, my.pop, sep="_")
            # print(paste(my.trait, my.env, my.pop, sep=";"))
            n <- n+1
            
          }
          
        }
      }
      # fwrite(as.tibble(result.id), "output_1to2_pop.txt")
      
      
      sln <- fread(paste0(temp.dir, "/asreml_1to2.sln"), data.table=F) %>%
        mutate(Set = result.id[cumsum(Model_Term == "Model_Term")+1]) %>% # set the cycle name used for these table rows
        filter(Model_Term != "Model_Term") %>%
        mutate(tempcol=str_split(Set, '_')) %>% # split the string
        rowwise() %>%
        mutate(trait=unlist(tempcol)[1], testing.env=unlist(tempcol)[2], pop=unlist(tempcol)[3]) %>%
        dplyr::select(-tempcol) %>%
        ungroup() 
      
      thing <- sln %>%
        filter(str_detect(Model_Term, '_EC')) %>%
        filter(!(str_detect(Model_Term, "Genotype"))) %>%
        mutate(tempcol=str_split(Model_Term, '_')) %>% # split the string
        rowwise() %>%
        mutate(check.trait=unlist(tempcol)[1]) %>%
        dplyr::select(-tempcol) %>%
        ungroup()
      stopifnot(all.equal(thing$trait, thing$check.trait))
      
      # tidy
      sln.full <- sln
      sln.full$pop <- as.numeric(sln.full$pop)
      sln.full$Effect <- as.numeric(sln.full$Effect)
      
      # calculate predictions:
      
      for (my.trait in c(FT.traits, harvest.traits)){
        pre.combine <- vector("list", length(unique(fr.exp_traits$pop_code))*nfolds*nrow(env_meta_info_0))
        n <- 1
        for(my.env in env_meta_info_0$env_code) {
          for (my.pop in sort(unique(fr.exp_traits$pop_code))){
            if (is.na(my.pop)) {next}
            
            # pull current data from sln file
            current.sln <- sln.full %>%
              filter(Set %in% sln.full$Set[which(sln.full$Model_Term==paste0(my.trait,"_EC"))]) %>%
              filter(pop==my.pop) %>%
              filter(testing.env==my.env) %>%
              distinct
            if(nrow(current.sln)==0) {next} # skip if no data
            
            # pull input data for this pop
            current.pop <- fread(paste0(temp.dir, "/asreml_1to2_no", my.env, "_pop",my.pop,".csv"), data.table=F) %>%
              dplyr::select(Genotype, env_code, paste0(my.trait), paste0(my.trait, "_EC"))%>%
              filter(is.na(UQ(rlang::sym(my.trait)))) %>% # only need to predict the ones that weren't known
              filter(env_code==my.env) %>% # only need the testing environ
              distinct
            
            current.pop$y.hat <- rep(NA, nrow(current.pop))
            current.pop$intcp.hat <- rep(NA, nrow(current.pop))
            current.pop$slope.hat <- rep(NA, nrow(current.pop))
            
            # and we'll need it for the opposite fold too; need to get it from some other testing environ too!
            other.pop <- fread(paste0(temp.dir, "/asreml_1to2_no",
                                      env_meta_info_0$env_code[!(env_meta_info_0$env_code == my.env)][1],
                                      "_pop",my.pop, ".csv"), data.table=F) %>%
              dplyr::select(Genotype, env_code, paste0(my.trait), paste0(my.trait, "_EC")) %>%
              filter(env_code==my.env) # only need the testing environ
            
            
            for (j in 1:nrow(current.pop)) {
              # calculate the prediction for each line
              current.genotype <- current.pop$Genotype[j]
              current.env <- current.pop$env_code[j]
              current.EC <- current.pop[j,4][[1]] # this only works because I select the column order above!
              current.intcp <- current.sln$Effect[current.sln$Level==current.genotype]
              current.slope <- current.sln$Effect[current.sln$Level==paste0(current.genotype, ".001")]
              
              current.pop$y.hat[j] <- current.sln$Effect[current.sln$Model_Term=="mu"]+ #mu
                current.sln$Effect[current.sln$Model_Term==paste0(my.trait, "_EC")]*current.EC + # + BXe
                current.intcp + # + uG
                current.slope*current.EC # + bG*xE
              # sln1b$Effect[sln1b$Level==current.env] # take out + uE
              current.pop$intcp.hat[j] <- current.intcp
              current.pop$slope.hat[j] <- current.slope
            }
            pre.combine[[n]] <- left_join(current.pop %>% select(-paste0(my.trait)),
                                          other.pop %>% select(-paste0(my.trait, "_EC")))
            # print(paste("Fold ", my.fold, ", pop ", my.pop, ": cor ", cor(pre.combine[[n]][,7], pre.combine[[n]]$y.hat, use="complete")))
            
            rm(current.pop)
            n <- n+1
            
          }
        }
        
        combined <- do.call(what = "rbind", pre.combine) %>%
          mutate(obs.id=NA) %>%
          select(env_code, ril_code=Genotype, paste0(my.trait), obs.id, intcp.hat, slope.hat, y.hat)
        
        # plot(combined$EH, combined$y.hat)
        # cor(combined$EH, combined$y.hat, use="complete")
        kPara_Name <- "Env Param"
        trait <- my.trait
        
        Plot_TEUG_result(obs_prd_file = combined, all_env_codes, kPara_Name, trait=trait, 
                         forecast_png_file = paste0(temp.dir, "/output/FR_", trait, "_UETG.1to2.pop_obs.prd.png"),
                         # forecast_png_file = paste0("NAM/FR.TEUG.1to3/no_covar", "FR_", trait, "_rep",i, "_TEUG.1to3_obs.prd.png"),
                         path = F, save.output = T)
        fwrite(combined, paste0(temp.dir, "/output/FR_", trait, "_UETG.1to2.pop_obs.prd.txt"))
      }
      
      
# FR-gBLUP 1->4 -----------------------------------------------------------
      
      
      # Goal: 1 ->4 prediction using windows from each training set
      
      # CERIS training set: n-1 environment, 50% of genotypes from each family
      # testing set: other 1 environment, other 50% of all genotypes
       # 1.	Use CERIS to find optimal window for training set, plus Slope and Intcp
      # 2. Use asreml to predict slope and intercept for other genotypes in testing environ
      # 3. Use environmental indices plus estimated slope and intercept to predict phenotype in testing set.
      # 4. Repeat until each environment has been in the testing set.
      # 5. Repeat 50x for each environment left out.
      
      source("genomic_prediction_functions.R")
      
      # This is designed for 2 folds-- split 50/50 training/testing
      nfolds <- 2
      rep.num <- c(12) # what rep to run -- eventually, want 1:50
      CERIS.done <- TRUE # use to determine if CERIS already run before FRgBLUP

      #### Do you want to leave one environment out while calculating correlations?
      LOO_env <- 0; ### leave as 0 here because I already manually left one out
      
      # set param for graphing:
      col_wdw <- 25;
      col_palette <- colorspace::diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1) # this function makes a diverging color palette
      
      # Run CERIS itself- - - -
      # read in geno with ril names in order:
      geno <- fread(file="NAM.GD.hidensity.matchY.txt", select=1)
      taxa <- geno$Taxa
      rm(geno)
      
      lInd <- which(colnames(exp_traits) == 'ril_code'); 
      eInd <- which(colnames(exp_traits) == 'env_code');
      exp_prepop_dir <- paste0(exp_dir, "UEUG.1to4.v4/"); if (!dir.exists(exp_prepop_dir))  { dir.create(exp_prepop_dir, recursive=T)};
      # exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop/"); if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir, recursive = T)};
      if(incl.precip=="no.precip") {exp_prepop_dir <- paste0(exp_dir, "UEUG.1to4.v4_no_precip", "/")} # make no precip directory
      if (!dir.exists(exp_prepop_dir))  { dir.create(exp_prepop_dir, recursive=T)};
      
      #### Do you want to leave one environment out while calculating correlations?
      LOO_env <- 0; ### leave as 0 here because I already manually left one out
      
      # set param for graphing:
      col_wdw <- 25;
      col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1) # this function makes a diverging color palette
      
      print(paste0("1 to 4 v4 for the whole population, rrBLUP within families, with parameters: line.outlier.filter=", line.outlier.filter,
                   ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
                   ", last.harvest=", last.harvest, ", max.window.size=", max.window))
      
      
      
      for (i in rep.num) {
        
        # first, make folds:
        # make random folds for this rep, across all subpops:
        taxa.table <- as.tibble(taxa)
        colnames(taxa.table) <- c("value")
        taxa.table$obs.id <- 1:length(taxa)
        taxa.table <- left_join(taxa.table, exp_traits, by=c("value"="ril_code")) %>%
          select(value, obs.id, pop_code) %>%
          distinct
        #make the folds--
        # for any number of folds:
        myF <- vector("list", nfolds)
        # preF <- vector("list", nfolds)
        for(subpop in unique(exp_traits$pop_code)) {
          # for(subpop in pop.table$number[c(1,6)]) {
          # get the obs id's for this pop:
          pop.ids <- taxa.table %>%
            filter(pop_code==subpop) %>%
            select(obs.id)
          
          preF <- fold.maker.2(pop.ids$obs.id, nfolds)
          # preF[[1]] <- pop.ids$obs.id[1:90]
          # preF[[2]] <- pop.ids$obs.id[-c(1:90)]
          for(j in c(1:nfolds)) {
            myF[[j]] <- c(myF[[j]], preF[[j]])
          }
        }
        
        # get optimal windows and slope + intercept, for each trait x fold x environment:
        output.results <- vector("list", 
                                 length(c(FT.traits, harvest.traits))*nrow(pop.table)*length(all_env_codes))
        n <- 1
        for(trait in c(FT.traits, harvest.traits)) {
          # for(trait in c("ERN")) {
          tInd <- which(colnames(exp_traits) == trait);
          exp_pop_dir <- paste0(exp_prepop_dir, trait, "/"); if (!(dir.exists(exp_pop_dir))) { dir.create(exp_pop_dir, recursive=T)};
          if(!dir.exists(paste0(exp_pop_dir, "/whole_pop/"))) {dir.create(paste0(exp_pop_dir, "/whole_pop/"), recursive=T)}
          current.trait <- trait
          
          for(env_num in seq_along(all_env_codes)) {
            current.env <- all_env_codes[env_num] # set current environment
            
            for(fold in c(1:nfolds)) {
              
              # pull training data
              exp_trait_trn <- exp_traits[exp_traits$env_code!=current.env,c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
              exp_trait_trn <- exp_trait_trn[!(exp_trait_trn$ril_code %in% taxa[myF[[fold]]]),]
              colnames(exp_trait_trn)[ncol(exp_trait_trn)] <- 'Yobs'; # rename phenotype data column as Yobs
              
              # also pull testing data:
              exp_trait_test <- exp_traits[exp_traits$env_code==current.env,c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
              exp_trait_test <- exp_trait_test[(exp_trait_test$ril_code %in% taxa[myF[[fold]]]),]; # filter by fold
              colnames(exp_trait_test)[ncol(exp_trait_test)] <- 'Yobs'; # rename phenotype data column as Yobs
              
              # skip this environment if there is no testing data in this fold:
              if(sum(!(is.na(exp_trait_test$Yobs)))==0) {next}
              
              # remove missing values--negative is ok though
              exp_trait_trn <- exp_trait_trn[!is.na(exp_trait_trn$Yobs),];
              # How many observations are there per pedigree?
              obs_lne_n <- aggregate(Yobs ~ ril_code, data = exp_trait_trn, length);
              # find lines with only one or two observations and remove them:
              line_outlier <- obs_lne_n$ril_code[obs_lne_n$Yobs < line.outlier.filter]
              exp_trait_trn <- exp_trait_trn[!(exp_trait_trn$ril_code %in% line_outlier),]
              line_codes <- unique(exp_trait_trn$ril_code); # pull unique line codes
              
              # and pull the training data as a vector only:
              current.data <- exp_trait_trn[,c(colnames(exp_trait_trn)=="Yobs")]
              
              # find mean value for each environment: 
              env_mean_trait_0_trn <- na.omit(aggregate(x = exp_trait_trn$Yobs, by = list(env_code = exp_trait_trn$env_code),
                                                        mean, na.rm = T));
              colnames(env_mean_trait_0_trn)[2] <- 'meanY';
              env_mean_trait_trn <- env_mean_trait_0_trn[order(env_mean_trait_0_trn$meanY),]; # sort by mean phenotype value
              
              # don't need this section for 1->4:
              # env_mean_trait_0_test <- na.omit(aggregate(x = exp_trait_test$Yobs, by = list(env_code = exp_trait_test$env_code),
              #                                       mean, na.rm = T));
              # colnames(env_mean_trait_0_test)[2] <- 'meanY';
              # env_mean_trait_test <- env_mean_trait_0_test
              
              # Perform exhaustive search in training set from WHOLE NAM pop, with 1 environment left out, for optimal window and environmental parameter to use
              Exhaustive_search_full(env_mean_trait_trn, PTT_PTR, searching_daps, 
                                     exp_pop_dir=paste0(exp_pop_dir,"whole_pop/"),
                                     current.data,
                                     trait,
                                     searching_daps, searching_daps, LOO_env, min.window.size = min.window.size)
              
              # need to read in default output and rename them:
              pop_cor_file.old <- paste(exp_pop_dir,'whole_pop/', trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_cor_whole_pop', sep = '');
              pop_pval_file.old <- paste(exp_pop_dir, 'whole_pop/',trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_P_whole_pop', sep = '');
              
              pop_cor <- fread(pop_cor_file.old)
              pop_pval <- fread(pop_pval_file.old) 
              
              pop_cor_file <- paste(exp_pop_dir, 'whole_pop/', trait, '_UEUG.1to4.no', current.env,'.fold', fold, ".rep", i, '_cor_whole_pop', sep = '')
              fwrite(pop_cor, pop_cor_file, sep="\t")
              
              pop_pval_file <-paste(exp_pop_dir, 'whole_pop/', trait, '_UEUG.1to4.no',current.env,'.fold', fold, ".rep", i, '_P_whole_pop', sep = '')
              fwrite(pop_pval, pop_pval_file, sep="\t")
              
              # tidy up:
              file.remove(pop_cor_file.old, pop_pval_file.old) 
              rm(pop_cor, pop_pval)
              
              # filter windows by established criteria:
              search.results <- fread(pop_cor_file)%>%
                select(-starts_with("nR")) # don't need the negative version
              search.results <- search.results %>%
                tidyr::gather(key="Parameter", value="Corr", -Day_x, -Day_y, -window)
              
              search.results <- search.results %>%
                filter(window >= min.window.size) %>% # i < j-5
                filter(window <= max.window)
              
              if(trait %in% FT.traits) {
                search.results <- search.results %>%
                  filter(Day_y <= last.FT) # last possible day for FT traits = 55
              } else if (trait %in% harvest.traits) {
                search.results <- search.results %>%
                  filter(Day_y <= last.harvest) # last possible day for harvest traits
              } else (warning("trait not in FT or harvest traits"))
              
              # Only now choose the remaining window with the best corr:
              search.results <- search.results <- search.results %>%
                filter(abs(Corr)==max(abs(Corr), na.rm = T))
              
              # if there are ties, just pick the top one in the table.
              search.results <- search.results[1,]%>%
                mutate(trait=trait)
              
              # output results:
              output.results[[n]] <- cbind(testing.env=current.env,rep=i, fold=fold, search.results)
              print(output.results[[n]])
              
              maxR_dap1 <- search.results$Day_x;
              maxR_dap2 <- search.results$Day_y;
              kPara_Name <- search.results$Parameter;
              kPara_Name <- gsub("R_", "", kPara_Name)
              PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); ## DL -> 5, GDD -> 6; PTT -> 7; PTR -> 8; PTD1 -> 8, PTD2 -> 9, PTS -> 10
              
              # Make output file containing the slopes and intercepts from the linear models for each line across environments
              Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait_trn, PTT_PTR, exp_trait_trn, line_codes, paste0(exp_pop_dir,"whole_pop/"), PTT_PTR_ind, 
                              filter.less.than-1); # do filtering for number of observations here
              # read file in and rename it:
              slope_file.old <-paste0(exp_pop_dir,"whole_pop/Intcp_Slope")
              pop_slope <- fread(slope_file.old)
              slope_file <- paste0(exp_pop_dir,"whole_pop/", trait, "_UEUG.1to4.no", current.env,".fold", fold, ".rep", i, "_Intcp_Slope")
              fwrite(pop_slope,slope_file, sep="\t")
              
              # tidy up:
              file.remove(slope_file.old) 
              rm(pop_slope)
              
              n <- n+1 # increase iterator
            }
            
          }
        }
        # output windows: 
        opt.windows <- do.call("rbind", output.results)
        if (incl.precip=="no.precip"){
          fwrite(opt.windows, paste0("NAM_BLUE/UEUG.1to4.v4_no_precip/UEUG.1to4_auto.optimal.parameters.whole_pop.rep",i, ".csv"))
          
        } else {
          fwrite(opt.windows, paste0("NAM_BLUE/UEUG.1to4.v4/UEUG.1to4_auto.optimal.parameters.whole_pop.rep",i, ".csv"))
          
        }
      }
      
      # read in geno with ril names in order:
      geno <- fread(file="NAM.GD.hidensity.matchY.txt", select=1)
      taxa <- geno$Taxa
      rm(geno)
      
      lInd <- which(colnames(exp_traits) == 'ril_code'); 
      eInd <- which(colnames(exp_traits) == 'env_code');
      exp_prepop_dir <- paste0(exp_dir, "UEUG.1to4.v4/"); if (!dir.exists(exp_prepop_dir))  { dir.create(exp_prepop_dir, recursive=T)};
      # exp_pop_dir <- paste0(exp_prepop_dir, "whole_pop/"); if (!dir.exists(exp_pop_dir))  { dir.create(exp_pop_dir, recursive = T)};
      if(incl.precip=="no.precip") {exp_prepop_dir <- paste0(exp_dir, "UEUG.1to4.v4_no_precip", "/")} # make no precip directory
      if (!dir.exists(exp_prepop_dir))  { dir.create(exp_prepop_dir, recursive=T)};
      
      #### Do you want to leave one environment out while calculating correlations?
      LOO_env <- 0; ### leave as 0 here because I already manually left one out
      
      # set param for graphing:
      col_wdw <- 25;
      col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1) # this function makes a diverging color palette
      
      print(paste0("1 to 4 v4 for the whole population, rrBLUP within families, with parameters: line.outlier.filter=", line.outlier.filter,
                   ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
                   ", last.harvest=", last.harvest, ", max.window.size=", max.window))
      
      
      
      for (i in rep.num) {
        
        # first, make folds:
        # make random folds for this rep, across all subpops:
        taxa.table <- as.tibble(taxa)
        colnames(taxa.table) <- c("value")
        taxa.table$obs.id <- 1:length(taxa)
        taxa.table <- left_join(taxa.table, exp_traits, by=c("value"="ril_code")) %>%
          select(value, obs.id, pop_code) %>%
          distinct
        #make the folds--
        # for any number of folds:
        myF <- vector("list", nfolds)
        # preF <- vector("list", nfolds)
        for(subpop in unique(exp_traits$pop_code)) {
          # for(subpop in pop.table$number[c(1,6)]) {
          # get the obs id's for this pop:
          pop.ids <- taxa.table %>%
            filter(pop_code==subpop) %>%
            select(obs.id)
          
          preF <- fold.maker.2(pop.ids$obs.id, nfolds)
          # preF[[1]] <- pop.ids$obs.id[1:90]
          # preF[[2]] <- pop.ids$obs.id[-c(1:90)]
          for(j in c(1:nfolds)) {
            myF[[j]] <- c(myF[[j]], preF[[j]])
          }
        }
        
        # get optimal windows and slope + intercept, for each trait x fold x environment:
        output.results <- vector("list", 
                                 length(c(FT.traits, harvest.traits))*nrow(pop.table)*length(all_env_codes))
        n <- 1
        for(trait in c(FT.traits, harvest.traits)) {
          # for(trait in c("ERN")) {
          tInd <- which(colnames(exp_traits) == trait);
          exp_pop_dir <- paste0(exp_prepop_dir, trait, "/"); if (!(dir.exists(exp_pop_dir))) { dir.create(exp_pop_dir, recursive=T)};
          if(!dir.exists(paste0(exp_pop_dir, "/whole_pop/"))) {dir.create(paste0(exp_pop_dir, "/whole_pop/"), recursive=T)}
          current.trait <- trait
          
          for(env_num in seq_along(all_env_codes)) {
            current.env <- all_env_codes[env_num] # set current environment
            
            for(fold in c(1:nfolds)) {
              
              # pull training data
              exp_trait_trn <- exp_traits[exp_traits$env_code!=current.env,c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
              exp_trait_trn <- exp_trait_trn[!(exp_trait_trn$ril_code %in% taxa[myF[[fold]]]),]
              colnames(exp_trait_trn)[ncol(exp_trait_trn)] <- 'Yobs'; # rename phenotype data column as Yobs
              
              # also pull testing data:
              exp_trait_test <- exp_traits[exp_traits$env_code==current.env,c(lInd, eInd, tInd)]; ### make sure the colname order is line env trait
              exp_trait_test <- exp_trait_test[(exp_trait_test$ril_code %in% taxa[myF[[fold]]]),]; # filter by fold
              colnames(exp_trait_test)[ncol(exp_trait_test)] <- 'Yobs'; # rename phenotype data column as Yobs
              
              # skip this environment if there is no testing data in this fold:
              if(sum(!(is.na(exp_trait_test$Yobs)))==0) {next}
              
              # remove missing values--negative is ok though
              exp_trait_trn <- exp_trait_trn[!is.na(exp_trait_trn$Yobs),];
              # How many observations are there per pedigree?
              obs_lne_n <- aggregate(Yobs ~ ril_code, data = exp_trait_trn, length);
              # find lines with only one or two observations and remove them:
              line_outlier <- obs_lne_n$ril_code[obs_lne_n$Yobs < line.outlier.filter]
              exp_trait_trn <- exp_trait_trn[!(exp_trait_trn$ril_code %in% line_outlier),]
              line_codes <- unique(exp_trait_trn$ril_code); # pull unique line codes
              
              # and pull the training data as a vector only:
              current.data <- exp_trait_trn[,c(colnames(exp_trait_trn)=="Yobs")]
              
              # find mean value for each environment: 
              env_mean_trait_0_trn <- na.omit(aggregate(x = exp_trait_trn$Yobs, by = list(env_code = exp_trait_trn$env_code),
                                                        mean, na.rm = T));
              colnames(env_mean_trait_0_trn)[2] <- 'meanY';
              env_mean_trait_trn <- env_mean_trait_0_trn[order(env_mean_trait_0_trn$meanY),]; # sort by mean phenotype value
              
              # don't need this section for 1->4:
              # env_mean_trait_0_test <- na.omit(aggregate(x = exp_trait_test$Yobs, by = list(env_code = exp_trait_test$env_code),
              #                                       mean, na.rm = T));
              # colnames(env_mean_trait_0_test)[2] <- 'meanY';
              # env_mean_trait_test <- env_mean_trait_0_test
              
              # Perform exhaustive search in training set from WHOLE NAM pop, with 1 environment left out, for optimal window and environmental parameter to use
              Exhaustive_search_full(env_mean_trait_trn, PTT_PTR, searching_daps, 
                                     exp_pop_dir=paste0(exp_pop_dir,"whole_pop/"),
                                     current.data,
                                     trait,
                                     searching_daps, searching_daps, LOO_env, min.window.size = min.window.size)
              
              # need to read in default output and rename them:
              pop_cor_file.old <- paste(exp_pop_dir,'whole_pop/', trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_cor_whole_pop', sep = '');
              pop_pval_file.old <- paste(exp_pop_dir, 'whole_pop/',trait, '_', nrow(env_mean_trait_trn), 'Envs_PTTPTR_', LOO_env, 'LOO_P_whole_pop', sep = '');
              
              pop_cor <- fread(pop_cor_file.old)
              pop_pval <- fread(pop_pval_file.old) 
              
              pop_cor_file <- paste(exp_pop_dir, 'whole_pop/', trait, '_UEUG.1to4.no', current.env,'.fold', fold, ".rep", i, '_cor_whole_pop', sep = '')
              fwrite(pop_cor, pop_cor_file, sep="\t")
              
              pop_pval_file <-paste(exp_pop_dir, 'whole_pop/', trait, '_UEUG.1to4.no',current.env,'.fold', fold, ".rep", i, '_P_whole_pop', sep = '')
              fwrite(pop_pval, pop_pval_file, sep="\t")
              
              # tidy up:
              file.remove(pop_cor_file.old, pop_pval_file.old) 
              rm(pop_cor, pop_pval)
              
              # filter windows by established criteria:
              search.results <- fread(pop_cor_file)%>%
                select(-starts_with("nR")) # don't need the negative version
              search.results <- search.results %>%
                tidyr::gather(key="Parameter", value="Corr", -Day_x, -Day_y, -window)
              
              search.results <- search.results %>%
                filter(window >= min.window.size) %>% # i < j-5
                filter(window <= max.window)
              
              if(trait %in% FT.traits) {
                search.results <- search.results %>%
                  filter(Day_y <= last.FT) # last possible day for FT traits = 55
              } else if (trait %in% harvest.traits) {
                search.results <- search.results %>%
                  filter(Day_y <= last.harvest) # last possible day for harvest traits
              } else (warning("trait not in FT or harvest traits"))
              
              # Only now choose the remaining window with the best corr:
              search.results <- search.results <- search.results %>%
                filter(abs(Corr)==max(abs(Corr), na.rm = T))
              
              # if there are ties, just pick the top one in the table.
              search.results <- search.results[1,]%>%
                mutate(trait=trait)
              
              # output results:
              output.results[[n]] <- cbind(testing.env=current.env,rep=i, fold=fold, search.results)
              print(output.results[[n]])
              
              maxR_dap1 <- search.results$Day_x;
              maxR_dap2 <- search.results$Day_y;
              kPara_Name <- search.results$Parameter;
              kPara_Name <- gsub("R_", "", kPara_Name)
              PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); ## DL -> 5, GDD -> 6; PTT -> 7; PTR -> 8; PTD1 -> 8, PTD2 -> 9, PTS -> 10
              
              # Make output file containing the slopes and intercepts from the linear models for each line across environments
              Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait_trn, PTT_PTR, exp_trait_trn, line_codes, paste0(exp_pop_dir,"whole_pop/"), PTT_PTR_ind, 
                              filter.less.than-1); # do filtering for number of observations here
              # read file in and rename it:
              slope_file.old <-paste0(exp_pop_dir,"whole_pop/Intcp_Slope")
              pop_slope <- fread(slope_file.old)
              slope_file <- paste0(exp_pop_dir,"whole_pop/", trait, "_UEUG.1to4.no", current.env,".fold", fold, ".rep", i, "_Intcp_Slope")
              fwrite(pop_slope,slope_file, sep="\t")
              
              # tidy up:
              file.remove(slope_file.old) 
              rm(pop_slope)
              
              n <- n+1 # increase iterator
            }
            
          }
        }
        # output windows: 
        opt.windows <- do.call("rbind", output.results)
        if (incl.precip=="no.precip"){
          fwrite(opt.windows, paste0("NAM_BLUE/UEUG.1to4.v4_no_precip/UEUG.1to4_auto.optimal.parameters.whole_pop.rep",i, ".csv"))
          
        } else {
          fwrite(opt.windows, paste0("NAM_BLUE/UEUG.1to4.v4/UEUG.1to4_auto.optimal.parameters.whole_pop.rep",i, ".csv"))
          
        }
        
      
      
      if(predict.type=="within.pop") {
        if(incl.precip=="no.precip") {
          temp.dir <- paste(exp_dir, "FR.UEUG.1to4.pop_no_precip", sep=""); if (!dir.exists(temp.dir))  { dir.create(temp.dir)};
        } else {
          temp.dir <- paste(exp_dir, "FR.UEUG.1to4.pop", sep=""); if (!dir.exists(temp.dir))  { dir.create(temp.dir)};
        }
      } else if (predict.type=="across.pop") {
        if(incl.precip=="no.precip") {
          temp.dir <- paste(exp_dir, "FR.UEUG.1to4.whole_no_precip", sep=""); if (!dir.exists(temp.dir))  { dir.create(temp.dir)};
        } else {
          temp.dir <- paste(exp_dir, "FR.UEUG.1to4.whole", sep=""); if (!dir.exists(temp.dir))  { dir.create(temp.dir)};
        }
      }

      print(paste0("FR-gBLUP 1 to 4 for the whole population, with parameters: line.outlier.filter=", line.outlier.filter,
                   ", filter.less.than=",filter.less.than,", min.window.size=", min.window.size, ", last.FT=",last.FT,
                   ", last.harvest=", last.harvest, ", max.window.size=", max.window))
      for (i in rep.num){
        if(CERIS.done) { # if CERIS-JGRA already run for 1->4, need to read in: optimal window, which line in which fold
          # first, read in folds:
          myF <-vector("list", nfolds)
          for (j in c(1:nfolds)) {
            fold.long <- vector()
            for(my.trait in c(FT.traits, harvest.traits)){ # need to loop through traits because some had missing data
              for(my.env in env_meta_info_0$env_code) {
                my.file <- paste0(exp_dir, "UEUG.1to4.v4/",my.trait,"/whole_pop/", my.trait, "_UEUG.1to4.no", my.env, ".fold",j,".rep", i, "_Intcp_Slope")
                
                if(file.exists(my.file)) {
                  fold <- fread(my.file) %>%
                    select(line_codes)
                  fold.long <- unique(c(fold.long, fold$line_codes))
                }
              }
            }
            myF[[j]] <- which(taxa %in% fold.long)
          }
          stopifnot(sum(taxa[setdiff(1:length(taxa), c(myF[[1]], myF[[2]]))] %in% exp_traits$ril_code)==0) # check that everything with a phenotype is in here
          stopifnot(length(intersect(myF[[1]], myF[[2]]))==0) # check that no lines are in both folds
          
          
        }   
        # Now, read in optimal windows:
        opt.windows <- fread(paste0(exp_dir, "UEUG.1to4.v4/UEUG.1to4_auto.optimal.parameters.whole_pop.rep", i, ".csv"))
        if(incl.precip=="no.precip") {
          opt.windows <- fread(paste0(exp_dir, "UEUG.1to4.v4_no_precip/UEUG.1to4_auto.optimal.parameters.whole_pop.rep", i, ".csv"))
        }
        # prep the input files for asreml:
        # get the EC (environmental covariate) value for each trait x fold x testing env combination:
        # EC.table <- tibble(rep=NA, fold=NA, trait=NA, env_code=NA, testing.env=NA, EC=NA)
        # sink("outfile.txt")
        for(my.trait in c(FT.traits, harvest.traits)) {
          # for(my.trait in c(FT.traits)[1:2]) {
          for(my.env in env_meta_info_0$env_code) {
            for(my.fold in 1:nfolds) {
              current.window <- opt.windows %>% filter(testing.env==my.env, trait==my.trait, fold==my.fold)
              if(nrow(current.window)==0) {
                # print(paste0("No window found for ", my.trait, " in training env ", my.env))
                next
              }
              # if(my.fold==1) { # only print them for one fold
              #   for(my.pop in unique(fr.exp_traits$pop_code)) {
              #     print(paste(my.trait, my.pop, my.env, sep="_"))
              #   }
              # }
              
              
              maxR_dap1 <- current.window$Day_x;
              maxR_dap2 <- current.window$Day_y;
              kPara_Name <- current.window$Parameter;
              kPara_Name <- gsub("R_", "", kPara_Name)
              
              current.EC <- PTT_PTR %>%
                group_by(env_code) %>%
                mutate(day.num=1:150) %>%
                filter(day.num %in% c(maxR_dap1:maxR_dap2)) %>%
                select(env_code, kPara_Name, day.num) %>%
                mutate(EC=mean(get(kPara_Name))) %>%
                mutate(trait=my.trait,
                       testing.env=my.env,
                       fold=my.fold,
                       rep=i) %>%
                select(rep, fold, trait, env_code, testing.env, EC) %>%distinct
              current.EC$EC <- scale(current.EC$EC,center = T,scale = T) # scale the environmental covariate as required for FR-gBLUP
              
              if(exists("EC.table")) {
                EC.table <- rbind(EC.table, current.EC)
              } else {
                EC.table <- current.EC
              }
              rm(current.EC)
            }
            
          }
        }
        # sink()
        EC.table <- EC.table %>% filter(!is.na(rep))
         wide.EC <- EC.table %>%
          spread(key=trait, value=EC)
        colnames(wide.EC)[-c(1:4)] <- paste0(  colnames(wide.EC)[-c(1:4)], ".EC")
        
        # add EC to the phenotype data:
        pheno.full <- left_join(fr.exp_traits, wide.EC, by="env_code") %>%
          select(-Entry_id)
        
        # make asreml input files, with different training sets set to NAs:
        for(my.env in env_meta_info_0$env_code) {
          pheno.env <- pheno.full %>%
            filter(testing.env==my.env) %>%
            mutate(
              ASI=ifelse(env_code==my.env, NA, ASI),
              CD=ifelse(env_code==my.env, NA, CD),
              CL=ifelse(env_code==my.env, NA, CL),
              CM=ifelse(env_code==my.env, NA, CM),
              DTA=ifelse(env_code==my.env, NA, DTA),
              DTS=ifelse(env_code==my.env, NA, DTS),
              EH=ifelse(env_code==my.env, NA, EH),
              EM=ifelse(env_code==my.env, NA, EM),
              ERN=ifelse(env_code==my.env, NA, ERN),
              KN=ifelse(env_code==my.env, NA, KN),
              KPR=ifelse(env_code==my.env, NA, KPR),
              TKW=ifelse(env_code==my.env, NA, TKW),
              LL=ifelse(env_code==my.env, NA, LL),
              LW=ifelse(env_code==my.env, NA, LW),
              PH=ifelse(env_code==my.env, NA, PH),
              T20KW=ifelse(env_code==my.env, NA, T20KW),
              TPBN=ifelse(env_code==my.env, NA, TPBN),
              TL=ifelse(env_code==my.env, NA, TL),
              ULA=ifelse(env_code==my.env, NA, ULA)) %>%
            select(Genotype=ril_code, env_code, pop_code, rep, fold, testing.env, everything())
          
          colnames(pheno.env) <- gsub("\\.", "_", colnames(pheno.env))
          
          for(my.pop in unique(pheno.env$pop_code)) {
            pheno.pop <- pheno.env %>%
              filter(pop_code==my.pop)
            for(my.fold in 1:nfolds) {
              fold.taxa <- taxa[myF[[my.fold]]]
              pheno.part <- pheno.pop %>%
                filter(fold==my.fold) %>%
                mutate(
                  ASI=ifelse(Genotype %in% fold.taxa, NA, ASI),
                  CD=ifelse(Genotype %in% fold.taxa, NA, CD),
                  CL=ifelse(Genotype %in% fold.taxa, NA, CL),
                  CM=ifelse(Genotype %in% fold.taxa, NA, CM),
                  DTA=ifelse(Genotype %in% fold.taxa, NA, DTA),
                  DTS=ifelse(Genotype %in% fold.taxa, NA, DTS),
                  EH=ifelse(Genotype %in% fold.taxa, NA, EH),
                  EM=ifelse(Genotype %in% fold.taxa, NA, EM),
                  ERN=ifelse(Genotype %in% fold.taxa, NA, ERN),
                  KN=ifelse(Genotype %in% fold.taxa, NA, KN),
                  KPR=ifelse(Genotype %in% fold.taxa, NA, KPR),
                  TKW=ifelse(Genotype %in% fold.taxa, NA, TKW),
                  LL=ifelse(Genotype %in% fold.taxa, NA, LL),
                  LW=ifelse(Genotype %in% fold.taxa, NA, LW),
                  PH=ifelse(Genotype %in% fold.taxa, NA, PH),
                  T20KW=ifelse(Genotype %in% fold.taxa, NA, T20KW),
                  TPBN=ifelse(Genotype %in% fold.taxa, NA, TPBN),
                  TL=ifelse(Genotype %in% fold.taxa, NA, TL),
                  ULA=ifelse(Genotype %in% fold.taxa, NA, ULA)) 
              for(my.trait in c(FT.traits, harvest.traits)) { # filter out lines without enough observations
                temp.pheno.part <- pheno.part[,c(1,which(colnames(pheno.part)==my.trait))]
                colnames(temp.pheno.part)[2] <- "Yobs" 
                obs_lne_n <- aggregate(Yobs~Genotype, data=temp.pheno.part, length)
                line_outlier <- obs_lne_n$Genotype[obs_lne_n$Yobs < (filter.less.than)] # want to make sure have AT LEAST filter.less.than-1 training obs, so need at least filter.less.than total obs
                
                temp.pheno.part$Yobs[temp.pheno.part$Genotype %in% line_outlier] <- NA
                pheno.part[,which(colnames(pheno.part)==my.trait)] <- temp.pheno.part$Yobs
              }
              pheno.part <- pheno.part %>% 
                select(Genotype,env_code,pop_code,rep,fold,testing_env,
                       ASI,CD,CL,CM,DTA,DTS,EH,EM,ERN,KN,KPR,TKW,LL,LW,
                       PH,T20KW,TPBN,TL,ULA,EH_EC,LL_EC,LW_EC,PH_EC,TPBN_EC,
                       TL_EC,ULA_EC,ASI_EC,DTA_EC,DTS_EC,CD_EC,CL_EC,CM_EC,
                       EM_EC,ERN_EC,KN_EC,KPR_EC,TKW_EC,T20KW_EC)
              fwrite(pheno.part, paste0(temp.dir, paste0("/asreml_1to4_no",my.env,"_pop", my.pop, "_fold", my.fold, "_rep", i, ".csv")))
              rm(pheno.part)
            }
          }
          
          
        }
        rm(EC.table)
      }
      
      # Compile FR 1->4 -- -- -- - --- -
      
      for (i in rep.num) {
        if(CERIS.done) {
          # first, read in folds:
          myF <-vector("list", nfolds)
          for (j in c(1:nfolds)) {
            fold.long <- vector()
            for(my.trait in c(FT.traits, harvest.traits)){ # need to loop through traits because some had missing data
              for(my.env in env_meta_info_0$env_code) {
                my.file <- paste0(exp_dir, "UEUG.1to4.v4/",my.trait,"/whole_pop/", my.trait, "_UEUG.1to4.no", my.env, ".fold",j,".rep", i, "_Intcp_Slope")
                
                if(file.exists(my.file)) {
                  fold <- fread(my.file) %>%
                    select(line_codes)
                  fold.long <- unique(c(fold.long, fold$line_codes))
                }
              }
            }
            myF[[j]] <- which(taxa %in% fold.long)
          }
          stopifnot(sum(taxa[setdiff(1:length(taxa), c(myF[[1]], myF[[2]]))] %in% exp_traits$ril_code)==0) # check that everything with a phenotype is in here
          stopifnot(length(intersect(myF[[1]], myF[[2]]))==0) # check that no lines are in both folds
          
        } 
        # Now, read in optimal windows:
        opt.windows <- fread(paste0(exp_dir, "UEUG.1to4.v4/UEUG.1to4_auto.optimal.parameters.whole_pop.rep", 1, ".csv"))
        
        # if predicting within each pop one at a time:
          # get the CYCLE arguments for asreml_rep.as:
          n <- 1
          result.id <- vector()
          for(my.trait in c(FT.traits, harvest.traits)) {
            for(my.env in env_meta_info_0$env_code) {
              for(my.fold in 1:nfolds) {
                current.window <- opt.windows %>% filter(testing.env==my.env, trait==my.trait, fold==my.fold)
                if(nrow(current.window)==0) {
                  # print(paste0("No window found for ", my.trait, " in training env ", my.env, " fold", my.fold))
                  next
                }
                if(my.fold==1) { # only print them for one fold
                  for(my.pop in unique(fr.exp_traits$pop_code)) {
                    result.id[n] <- paste(my.trait, my.pop, my.env, sep="_")
                    n <- n+1
                  }
                }
                
              }
              
            }
          }
        
        
         # Need to read in each fold separately:
        sln.1 <- fread(paste0(temp.dir, "/asreml_1to4_fold1_full_rep", i, ".sln"), data.table=F) %>%
          # sln.1 <- fread(paste0(temp.dir, "/asreml_1to4_fold1_full_rep", i, "_begin1.sln"), data.table=F) %>%
          mutate(Set = result.id[cumsum(Model_Term == "Model_Term")+1]) %>% # set the cycle name used for these table rows
          filter(Model_Term != "Model_Term") %>%
          mutate(tempcol=str_split(Set, '_')) %>% # split the string
          rowwise() %>%
          mutate(trait=unlist(tempcol)[1], pop=unlist(tempcol)[2], testing.env=unlist(tempcol)[3]) %>%
          dplyr::select(-tempcol) %>%
          ungroup() %>%
          mutate(fold=1)
        
        # double check:
        thing <- sln.1 %>%
          filter(str_detect(Model_Term, '_EC')) %>%
          filter(!(str_detect(Model_Term, "Genotype"))) %>%
          mutate(tempcol=str_split(Model_Term, '_')) %>% # split the string
          rowwise() %>%
          mutate(check.trait=unlist(tempcol)[1]) %>%
          dplyr::select(-tempcol) %>%
          ungroup()
        stopifnot(all.equal(thing$trait, thing$check.trait))
        tail(unique(sln.1$Set))
        
        sln.2 <- fread(paste0(temp.dir, "/asreml_1to4_fold2_full_rep", i, ".sln"), data.table=F) %>%
          # sln.2 <- fread(paste0(temp.dir, "/asreml_1to4_fold2_full_rep", i, "_begin.sln"), data.table=F) %>%
          mutate(Set = result.id[cumsum(Model_Term == "Model_Term")+1]) %>% # set the cycle name used for these table rows
          filter(Model_Term != "Model_Term") %>%
          mutate(tempcol=str_split(Set, '_')) %>% # split the string
          rowwise() %>%
          mutate(trait=unlist(tempcol)[1], pop=unlist(tempcol)[2], testing.env=unlist(tempcol)[3]) %>%
          dplyr::select(-tempcol) %>%
          ungroup() %>%
          mutate(fold=2)
        
        # double check:
        thing <- sln.2 %>%
          filter(str_detect(Model_Term, '_EC')) %>%
          filter(!(str_detect(Model_Term, "Genotype"))) %>%
          mutate(tempcol=str_split(Model_Term, '_')) %>% # split the string
          rowwise() %>%
          mutate(check.trait=unlist(tempcol)[1]) %>%
          dplyr::select(-tempcol) %>%
          ungroup()
        stopifnot(all.equal(thing$trait, thing$check.trait))
        tail(unique(sln.2$Set))
        
        # combine the folds, tidy formatting:
        sln.full <- rbind(sln.1, sln.2)
        sln.full$pop <- as.numeric(sln.full$pop)
        sln.full$fold <- as.numeric(sln.full$fold)
        sln.full$Effect <- as.numeric(sln.full$Effect)
        
        # calculate predictions:
        for (my.trait in c(FT.traits, harvest.traits)){
          # for (my.trait in c(FT.traits, harvest.traits)[-c(15:19)]){
          pre.combine <- vector("list", length(unique(fr.exp_traits$pop_code))*nfolds*nrow(env_meta_info_0))
          n <- 1
          for(my.env in env_meta_info_0$env_code) {
            for (my.pop in sort(unique(fr.exp_traits$pop_code))){
              if (is.na(my.pop)) {next}
              for(my.fold in 1:nfolds) {
                # pull current data from sln file
                current.sln <- sln.full %>%
                  filter(Set %in% sln.full$Set[which(sln.full$Model_Term==paste0(my.trait,"_EC"))]) %>%
                  filter(pop==my.pop) %>%
                  filter(fold==my.fold) %>%
                  filter(testing.env==my.env) %>%
                  distinct
                if(nrow(current.sln)==0) {next} # skip if no data
                
                # pull input data for this fold
                current.fold <- fread(paste0(temp.dir, "/asreml_1to4_no", my.env, "_pop",my.pop, "_fold",my.fold,"_rep", i,".csv"), data.table=F) %>%
                  dplyr::select(Genotype, env_code, paste0(my.trait), paste0(my.trait, "_EC"))%>%
                  filter(is.na(UQ(rlang::sym(my.trait)))) %>% # only need to predict the ones that weren't known
                  filter(env_code==my.env) %>% # only need the testing environ
                  filter((Genotype %in% taxa[myF[[my.fold]]])) %>% # only need the testing genotypes
                  distinct
                
                current.fold$y.hat <- rep(NA, nrow(current.fold))
                current.fold$intcp.hat <- rep(NA, nrow(current.fold))
                current.fold$slope.hat <- rep(NA, nrow(current.fold))
                
                # and we'll need it for the opposite fold too; need to get it from some other testing environ too!
                other.fold <- fread(paste0(temp.dir, "/asreml_1to4_no",
                                           env_meta_info_0$env_code[!(env_meta_info_0$env_code == my.env)][1],
                                           "_pop",my.pop, "_fold",c(1:nfolds)[-my.fold],"_rep", i,".csv"), data.table=F) %>%
                  dplyr::select(Genotype, env_code, paste0(my.trait), paste0(my.trait, "_EC")) %>%
                  filter((Genotype %in% taxa[myF[[my.fold]]])) %>% # only need testing genotypes
                  filter(env_code==my.env) # only need the testing environ
                
                
                for (j in 1:nrow(current.fold)) {
                  # calculate the prediction for each line
                  current.genotype <- current.fold$Genotype[j]
                  current.env <- current.fold$env_code[j]
                  current.EC <- current.fold[j,4][[1]] # this only works because I select the column order above!
                  current.intcp <- current.sln$Effect[current.sln$Level==current.genotype]
                  current.slope <- current.sln$Effect[current.sln$Level==paste0(current.genotype, ".001")]
                  
                  current.fold$y.hat[j] <- current.sln$Effect[current.sln$Model_Term=="mu"]+ #mu
                    current.sln$Effect[current.sln$Model_Term==paste0(my.trait, "_EC")]*current.EC + # + BXe
                    current.intcp + # + uG
                    current.slope*current.EC # + bG*xE
                  # sln1b$Effect[sln1b$Level==current.env] # take out + uE
                  current.fold$intcp.hat[j] <- current.intcp
                  current.fold$slope.hat[j] <- current.slope
                }
                # pre.combine[[n]] <- left_join(current.fold %>% select(-paste0(my.trait)),
                #                               other.fold %>% select(-paste0(my.trait, "_EC")))
                # print(paste("Fold ", my.fold, ", pop ", my.pop, ": cor ", cor(pre.combine[[n]][,7], pre.combine[[n]]$y.hat, use="complete")))
                
                # to use on Nova:
                my.EC <- paste0(my.trait, "_EC")
                pre.combine[[n]] <- left_join(current.fold %>% select(-my.trait),
                                              other.fold %>% select(-my.EC))
                
                rm(current.fold)
                n <- n+1
              }
            }
          }
          
          combined <- do.call(what = "rbind", pre.combine) %>%
            mutate(obs.id=NA) %>%
            select(env_code, ril_code=Genotype, paste0(my.trait), obs.id, intcp.hat, slope.hat, y.hat)
          
          # plot(combined$EH, combined$y.hat)
          # cor(combined$EH, combined$y.hat, use="complete")
          kPara_Name <- "Env Param"
          trait <- my.trait
          
          Plot_TEUG_result(obs_prd_file = combined, all_env_codes, kPara_Name, trait=trait, 
                           forecast_png_file = paste0(temp.dir, "/output/FR_", trait, "_rep",i, "_UEUG.1to4.pop_obs.prd.png"),
                           # forecast_png_file = paste0("NAM/FR.TEUG.1to3/no_covar", "FR_", trait, "_rep",i, "_TEUG.1to3_obs.prd.png"),
                           path = F, save.output = T)
          fwrite(combined, paste0(temp.dir, "/output/FR_", trait, "_rep",i, "_UEUG.1to4.pop_obs.prd.txt"))
          
        }
        
      }
      
 
    
    
# GCTA GWAS --------------------------------------------------------------------

# Use GCTA to get high-density GWAS results

# on server, run the following: 
# NOTE: the following lines are to be run on the command line, NOT R
# Then, the results can be used in R
module load gcta/1.91.2beta-q2ognx7
module load plink

# merge the SV and SNP data
plink --merge-list sv_and_snp_mb --memory 40000 --out merge # memory says to limit to 40GB, sv_and_snp_mb has the list of all files to combine (SNP files plus SV ones)

# make files to run SV+SNP GWAS:
# make GRM (K), took ~11 hr
nohup gcta64 --bfile merge --make-grm --maf 0.01 --out sv.snp.grm &
  
  # make it sparse: took 2 seconds
  gcta64 --grm sv.snp.grm --make-bK-sparse 0.05 --out sparse.grm

# need to make maf-filtered version: took ~25 min
# nohup gcta64 --bfile merge --maf 0.01 --out merge.maf &
plink --bfile merge --memory 40000 --maf 0.01 --make-bed --out merge.maf01

# and run GWAS:
gcta.format <- fread("merge.maf01.fam") %>% select(V1, V2)
gcta.env <- left_join(gcta.format, fr.exp_traits, by=c("V1"="ril_code"))
env.list <- sort(unique(fr.exp_traits$env_code))
for(current.trait in c(FT.traits, harvest.traits)) {
  
  # read in phenotypes
  myY <- fread(paste0(exp_dir, current.trait, "/whole_pop/Intcp_Slope"), header = T, stringsAsFactors = F, data.table=F)
  
  # format for GCTA GWAS
  gcta.y <- left_join(gcta.format, 
                      myY %>% 
                        select(line_codes, Intcp_para_adj, Slope_para) %>% 
                        distinct,
                      by=c("V1"="line_codes"))
  
  # write slope pheno file
  fwrite(gcta.y %>% select(V1,V2,slope.hat=Slope_para),
         paste0(exp_dir, "GCTA_fixed/", current.trait, "_slope_fixed.phen"),
         sep="\t", col.names = F, na="NA", quote=F)
  
  # write intcp pheno file
  fwrite(gcta.y %>% select(V1,V2,intcp.hat=Intcp_para_adj),
         paste0(exp_dir, "GCTA_fixed/", current.trait, "_intcp_fixed.phen"),
         sep="\t", col.names = F, na="NA", quote = F)
  
  for(env_num in seq_along(env.list)){ # make phenotype in each environment
    my.env <- env.list[env_num]
    fwrite(gcta.env %>% filter(env_code==my.env) %>% select(V1, V2, current.trait),
           paste0(exp_dir, "GCTA_fixed/", current.trait, "_BLUE_", env_num, ".phen"),
           sep="\t", col.names=F, na="NA", quote=F)
  }
}

# Calculate SimpleM threshold for high-density:
# use https://github.com/LTibbs/SimpleM/blob/main/simpleM_efficient.R
# Result was:
threshold <- 4.750729e-07

# Run GCTA GWAS: on server (this is command line code, NOT R code)
# need to run for each trait; example (DTA) shown here
# This uses slurm array numbers to iterate through environments
module load gcta/1.94.0b
gcta64 --mlma --bfile merge.maf01 --grm sv.snp.grm --pheno DTA_intcp.phen --out geno_assoc_DTA_intcp --thread-num 28
gcta64 --mlma --bfile merge.maf01 --grm sv.snp.grm --pheno DTA_BLUE_$SLURM_ARRAY_TASK_ID.phen --out geno_assoc_DTA_BLUE_$SLURM_ARRAY_TASK_ID --thread-num 28

# Make Manhattan plots and P value plots for 20M snps/svs
for(my.trait in c(FT.traits, harvest.traits)) {
  
  untar(paste0("GCTA_fixed/results_tar/",my.trait,".mlma.fixed.tar.gz"),
        files=paste0("geno_assoc_", my.trait, "_slope_fixed.mlma"),
        list=F, exdir="GCTA_fixed/results_tar/temp")
  slope.snp <- fread(paste0("/GCTA_fixed/results_tar/temp/geno_assoc_",
                            my.trait, "_slope_fixed.mlma"),
                     select=c(1:3,9))
  
  # untar intcp, then read in:
  untar(paste0("GCTA_fixed/results_tar/",my.trait,".mlma.fixed.tar.gz"),
        files=paste0("geno_assoc_", my.trait, "_intcp_fixed.mlma"),
        list=F, exdir="GCTA_fixed/results_tar/temp")
  int.snp <- fread(paste0("GCTA_fixed/results_tar/temp/geno_assoc_",
                          my.trait, "_intcp_fixed.mlma"),
                   select=c(1:3,9))
  min.p <- min(c(int.snp$p, slope.snp$p), na.rm=T)
  
  png(paste0('Figures/', my.trait, '_fixed_int_slope_fastman.png'), width =6.5, height = 8, units = 'in', res = 300)
  layout(matrix(c(1:2),2,1))
  fastman(int.snp,
          suggestiveline=F,
          genomewideline=-log10(threshold),
          chr="Chr", bp = "bp", p = "p",
          maxP=-log10(min.p)
  )
  fastman(slope.snp,
          suggestiveline=F,
          genomewideline=-log10(threshold),
          chr="Chr", bp = "bp", p = "p",
          maxP=-log10(min.p)
  )
  dev.off()
  
  # Make plots of intercept vs slope P values
  stopifnot(all.equal(slope.snp$SNP, int.snp$SNP))
  pdat <- tibble(trait=my.trait,
                 neg.log.int.p=-log10(int.snp$p),
                 neg.log.slope.p=-log10(slope.snp$p))
  max.p <- max(c(pdat$neg.log.int.p, pdat$neg.log.slope.p),
               na.rm = T)
  pplot <- ggplot(pdat) +
    scattermore::geom_scattermore(aes(x=neg.log.int.p, y=neg.log.slope.p),
                                  pointsize=2.5)+
    facet_wrap(.~trait) +
    geom_abline(slope=1, intercept=0) +
    geom_vline(xintercept=-log10(threshold), color="red") +
    geom_hline(yintercept=-log10(threshold), color="red") +
    xlab("Intercept (-log P)")+
    ylab("Slope (-log P)")+
    xlim(0,max.p)+
    ylim(0,max.p)
  ggsave(paste0("Figures/",my.trait,"pplot.jpg"), plot=pplot, device="jpeg", width=2.5, height=2.5, units="in", dpi=500)
  
  
  rm(slope.snp, int.snp, min.p)
}

# PCA ---------------------------------------------------------------------

      # env mean:
      # For PC of the traits, we want to have the environments as column names plus have a "Trait" column
      pc.mean <- fr.exp_traits %>%
        select(-c(pop_code, Entry_id)) %>%
        pivot_longer(-c(env_code, ril_code), names_to="Trait", values_to="BLUE") %>%
        group_by(env_code, Trait) %>%
        mutate(env.mean=mean(BLUE, na.rm = T)) %>%
        select(-ril_code, -BLUE) %>%
        distinct
      pc.mean.traits <- pc.mean %>% # for the PCs of the traits using the env mean:
        pivot_wider(names_from=Trait, values_from=env.mean)
      plot.pc.mean.traits <- prcomp(na.omit(pc.mean.traits[,-1]), center=T, scale=T)$rotation
      
      
