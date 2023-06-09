
for (package in c('optparse', 'rmarkdown','kableExtra','knitr')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos = "http://cran.us.r-project.org")
    library(package, character.only=T)
  }
}

#Parse arguments
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix", metavar="character"),
  make_option(c("-m", "--mutations"), type="character", default=NULL, 
              help="mutations table", metavar="character"),
  make_option(c("-c", "--cnv_plot"), type="character", default=NULL,
              help="cnv plot", metavar="character"),
  make_option(c("-r", "--rf_details"), type="character", default=NULL,
              help="RF details tsv", metavar="character"),
  make_option(c("-v", "--votes"), type="character", default=NULL,
              help="votes file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-n", "--patient"), type="character", default=NULL, 
              help="patient", metavar="character"),
  make_option(c("-e", "--coverage"), type="character", default=NULL,
              help="coverage summary", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample", metavar="character"),
  make_option(c("-t", "--mgmt"), type="character", default=NULL,
              help="mgmt prediction", metavar="character"),
  ## following option added by GF
  make_option(c("-u", "--report_UKHD"), type="character", default=NULL,
              help="report_UKHD R markdown doc", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
prefix <- opt$prefix
mutations <- opt$mutations
rf_details <- opt$rf_details
cnv_plot <- opt$cnv_plot
votes <- opt$votes
coverage <- opt$coverage
patient <- opt$patient
mgmt<- opt$mgmt
sample <- opt$sample
# add option for report_UKHD above
report_UKHD <- opt$report_UKHD

#render("/Rapid_CNS2_report_UKHD.Rmd", 
render(report_UKHD, 
       output_format = "html_document", 
       output_dir=opt$output_dir,
       output_file=paste0(prefix,"_Rapid-CNS2_report.html")
)
