#!/usr/bin/env Rscript
library(optparse)
#library(here)
library(this.path)
src_path <- this.dir(verbose = getOption("verbose"))
# get arguments
option_list <- list(
    make_option(c("-i", "--infile"), type="character", default=NULL, 
            help="[REQUIRED] COI input file path", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
            help="output file name [default= %default]", metavar="character"),
    make_option("--maxCOI", type="integer", default=25, 
            help="max COI threshold [default= %default]", metavar="number"), 
    make_option("--maxMissTol", type="integer", default=20,
            help="maximum tolerance for missing data for individuals/sample and site[default= %default]",
            metavar="number"),
    make_option("--totalRun", type="integer", default=100000,
            help="Total 'production' iteration steps [default= %default]",
            metavar="number"),
    make_option("--totalBurnIn", type="integer", default=1000,
            help="Total burn in steps [default= %default]",
            metavar="number"),
    make_option("--M0", type="integer", default=20,
            help="Initial value for COI [default= %default]",
            metavar="number"),
    make_option("--e1", type="double", default=0.05,
            help="probability of errouniously call homozygous [default= %default]",
            ),
    make_option("--e2", type="double", default=0.05,
            help="probability of errouniously call heterozygous [default= %default]",
            metavar="number"),
    make_option("--estimate_e", type="logical", default=FALSE,
            help="Estimate e1 and e2 [default= %default]"),
    make_option("--outdir", type="character", default=getwd(),
            help=" path to write output file [default= current working directory]",
            metavar="character"),
    make_option("--src_path", type="character", default=src_path,
            help="source code path [default=location of the repository]",
            metavar="character"),
    make_option("--outPrefix", type="character", default="COIout",
            help="Prefix for output files [default='COIout']",
            metavar="character"),
    make_option("--seed", default=NULL,
            help="set seed for random number generator [default= NULL]")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# TODO implement arguments sanity checks

# check if required was provided
if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$seed) == FALSE){
  # set random number seed
  set.seed(opt$seed)
}


print(" ---| THERealMcCOIL.R |---")
print(" PARAMETERS: ")
print("  -- paths --")
print(paste("  infile      :", opt$infile))
print(paste("  outdir      :", opt$outdir))
print(paste("  src_path    :", opt$src_path))
print(paste("  outPrefix   :", opt$outPrefix))
print(" -- algorithm settings --")
print(paste("  maxCOI      :", opt$maxCOI))
print(paste("  maxMissTol  :", opt$maxMissTol))
print(paste("  totalRun    :", opt$totalRun))
print(paste("  totalBurnIn :", opt$totalBurnIn))
print(paste("  M0          :", opt$M0))
print(paste("  e1          :", opt$e1))
print(paste("  e2          :", opt$e2))
print(paste("  estimate_e  :", opt$estimate_e))
print(paste("  seed        :", opt$seed))
print(" -------------------------")

# set err_method
# if estimate e values, then set err_metho to 3
if (opt$estimate_e == TRUE){
    err_method <- 3
}  
# if use e values provided, then set err_method to 1
if (opt$estimate_e == FALSE){
    err_method <- 1
}

source(paste(opt$src_path, "categorical_method/McCOIL_categorical.R", sep="/"))

#example dataset
print("@ loading data...")
data0= read.table(opt$infile, head=T)
data=data0[,-1]
rownames(data)=data0[,1]
print("@ Run McCOIl...")
McCOIL_categorical(data,
                maxCOI=opt$maxCOI,
                threshold_ind=opt$maxMissTol,
                threshold_site=opt$maxMissTol, 
                totalrun=opt$totalRun,
                burnin=opt$totalBurnIn, 
                M0=opt$M0,
                e1=opt$e1,
                e2=opt$e2,
                err_method=err_method,
                path=opt$outdir,
                output= opt$outPrefix,
                src_path=opt$src_path)
print(":: DONE ::")