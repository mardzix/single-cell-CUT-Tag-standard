cat("*** Loading libraries*** \n")
library(argparse)
library(yaml)
library(ggplot2)
library(funr)
# library(rtracklayer)
# library(GenomicRanges)


# Source aux functions
source(paste0(dirname(funr::sys.script()),"/func.R"))
set.seed(1234)

########### Arguments parser

parser <- ArgumentParser()

parser$add_argument("-s", "--sample", type="character", default='foo', 
                    help="sample name [as in config file key]")

parser$add_argument("-a", "--antibody", type="character", default='foo',
                    help="antibody name [as in config file key]")

parser$add_argument("-o", "--out_prefix", type="character", default="/",
                    help="folder for the output in clustering_snakemake folder")

parser$add_argument("--bcd_all", type="character",
                    help="Path to barcodes summary statistics per cell")

parser$add_argument("--bcd_peak", type="character",
                    help="Path to barcodes in peaks summary statistics per cell")

parser$add_argument("--metadata", type="character",
                    help="Path to the cellranger metadata singlecell.csv file")

parser$add_argument("--fragments", type="character",
                    help="Path to the cellranger fragments.tsv.gz file")

parser$add_argument("--min_reads", type="double", default='3.0',
                    help="Minimum number of reads per cell")

parser$add_argument("--max_reads", type="double", default='5.5',
                    help="Maximum number of reads per cell")

parser$add_argument("--peak_fraction_min", type="double", default='0.2',
                    help="Minimum fraction of reads within peak")

parser$add_argument("--peak_fraction_max", type="double", default='1',
                    help="Minimum fraction of reads within peak")



args      <- parser$parse_args()
print(args)

cutoff_reads_min            = args$min_reads
cutoff_reads_max            = args$max_reads
cutoff_peak_percentage_low  = args$peak_fraction_min
cutoff_peak_percentage_high = args$peak_fraction_max

# Load ensembl annotation
genebodyandpromoter.coords.flat <- load_ensembl_annot()

########## Filter the barcodes
cat("*** Reading barcode statistics files \n")

all_barcodes_file  <- args$bcd_all
peak_barcodes_file <- args$bcd_peak
metadata_file      <- args$metadata

metadata = read.csv(metadata_file, header = 1)
metadata = metadata[2:nrow(metadata),]
metadata$logUMI = log10(metadata$passed_filters + 1)
metadata$promoter_ratio = (metadata$promoter_region_fragments+1) / (metadata$passed_filters + 1)
metadata$peak_region_ratio = (metadata$peak_region_fragments+1) / (metadata$passed_filters + 1)

# Fix metadata cell barcodes
metadata$barcode <- paste0(args$sample,"_",metadata$barcode)

# Read barcode statistics files
all_barcodes <- read.table(file=all_barcodes_file)
peak_barcodes <- read.table(file=peak_barcodes_file)
bcd <- merge(all_barcodes,peak_barcodes,by="V2")  
colnames(bcd) <- c("barcode","all_unique_MB","peak_MB")
bcd$peak_ratio_MB <- bcd$peak_MB/bcd$all_unique_MB
bcd$sample <- args$sample


# Merge 10x metadata with barcode statistics
metadata <- merge(metadata,bcd,by='barcode')


metadata$is__cell_barcode <- as.factor(metadata$is__cell_barcode)

print(head(metadata))
################ MB filtering
cat("*** Filtering cells \n")
metadata[,"passedMB"] <- FALSE
metadata[metadata$all_unique_MB > 10^cutoff_reads_min &
         metadata$all_unique_MB < 10^cutoff_reads_max &
         metadata$peak_ratio_MB > cutoff_peak_percentage_low &
         metadata$peak_ratio_MB < cutoff_peak_percentage_high,"passedMB"] <- TRUE


################ Cell picking scatterplot nreads ~ percent in peaks
dir.create(args$out_prefix,recursive=TRUE)

p1 <- ggplot(data = metadata,aes(x=log10(all_unique_MB),y=peak_ratio_MB,col=is__cell_barcode)) + 
  geom_point(size=0.1) + 
  scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(as.numeric(as.character(metadata$is__cell_barcode)))),"FALSE")) +
  theme(legend.position="bottom",text=element_text(size=26)) 

p2 <- ggplot(data = metadata,aes(x=log10(passed_filters),y=peak_region_fragments/passed_filters,col=is__cell_barcode)) + 
  geom_point(size=0.1) +
  scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(as.numeric(as.character(metadata$is__cell_barcode)))),NA)) + 
  theme(legend.position="bottom",text=element_text(size=26))

ggsave(plot = p1+p2,
       filename=paste0(args$out_prefix,'cells_10x.png'),width = 20,height = 10,units = 'in')


p1 <- ggplot(data = metadata) + 
      geom_point(aes(x=log10(all_unique_MB),y=peak_ratio_MB,col=passedMB),size=0.1) + 
      geom_hline(yintercept = c(cutoff_peak_percentage_high,cutoff_peak_percentage_low)) + 
      geom_vline(xintercept = c(cutoff_reads_min,cutoff_reads_max)) + 
      scale_color_manual(values=c("black","gold"),labels=c(paste("TRUE",sum(metadata$passedMB)),"FALSE")) +
      theme(legend.position="bottom",text=element_text(size=26))
  
p2 <- ggplot(data = metadata) + 
      geom_point(aes(x=log10(passed_filters),y=peak_region_fragments/passed_filters,col=passedMB),size=0.1) +
      scale_color_manual(values=c("black","gold")) + 
      theme(legend.position="bottom",text=element_text(size=26))

ggsave(plot = p1+p2,
       filename=paste0(args$out_prefix,'cells_picked.png'),width = 20,height = 10,units = 'in')


################# Export bw selected / unselected
cat("*** Reading fragments file \n")

fragments         <- args$fragments
fragments_gr      <- rtracklayer::import(fragments,format = "bed")


cat("*** Exorting merged bw files \n")
barcode_pass   <- metadata$barcode[metadata$passedMB]
barcode_nopass <- metadata$barcode[!metadata$passedMB]
  
fragments.pass   <- fragments_gr[fragments_gr$name %in% barcode_pass]
fragments.nopass <- fragments_gr[fragments_gr$name %in% barcode_nopass]
  
cat("*** Calculating coverage \n")
coverage.pass   <- GenomicRanges::coverage(fragments.pass)
coverage.nopass <- GenomicRanges::coverage(fragments.nopass)
  
cat("*** Normalizing \n")
coverage.pass <- coverage.pass/length(fragments.pass)
coverage.nopass <- coverage.nopass/length(fragments.nopass)
  
cat("*** Exporting \n")
rtracklayer::export(object=coverage.pass,  con = paste0(args$out_prefix,'cells_picked.bw'))
rtracklayer::export(object=coverage.nopass,con = paste0(args$out_prefix,'cells_not_picked.bw'))

write.csv(x = metadata,file = paste0(args$out_prefix,'metadata.csv'))