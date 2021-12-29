library(tidyverse)
library(deepSNV)
library(VariantAnnotation)
library(optparse)

do_shearwater_full = function(counts, samples, muts_2_call, cutoff = 0.05){
    pvals = betabinLRT(counts, rho=1e-4, maxtruncate = 1, truncate = 1)$pvals
    qvals = p.adjust(pvals, method="BH")
    dim(qvals) = dim(pvals)
    vcf = qvals2Vcf(qvals, counts, muts_2_call, samples = samples, mvcf = TRUE, cutoff = cutoff)
    gr = granges(vcf)

    #Ensure that the alt and ref allele are not the same. (This sometimes happens.)
    alts = gr$ALT %>%
        unlist() %>%
        as.vector()
    if (length(alts) != length(gr)){
        stop("Some variants seem to have multiple alternative alleles.", call. = F)
    }
    refs = gr$REF %>%
        as.vector()
    mut_f = refs != alts
    vcf = vcf[mut_f]

    #Remove variants with no correct ref allele
    gr = granges(vcf)
    vcf = vcf[gr$REF != "-"]

    #Set infinite gqs to 99
    gq = geno(vcf)$GQ
    gq[is.infinite(gq)] = 99
    geno(vcf)$GQ = gq

    #vcf = vcf[isSNV(vcf)]
    return(vcf)
}

option_list = list(
  make_option(c("-q", "--QUAL"), type="integer", default=10,
              help="Minimum base quality to count"),
  make_option(c("-m", "--MQ"), type="integer", default=0,
              help="The output path")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

q = opt$QUAL
mq = opt$MQ

#Set output directory
out_dir = str_c("/hpc/pmc_vanboxtel/processed/VAN4509/targetcapture/", "q", q, "_mq", mq)
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)
if (!dir.exists("counts")){
    dir.create("counts")
}


#Read bed regions
muts_2_call = read_tsv("/hpc/pmc_vanboxtel/processed/VAN4509/targetcapture/muts.bed", col_types = "cii", col_names = c("Chrom", "Start", "End"), comment = "#") %>%
    makeGRangesFromDataFrame(starts.in.df.are.0based = T) %>%
    sort()

#Determine bam files
samples = c("N01P0ISRB", "N01SKINRB", "N01SPLEENRB", "NR1P0ISRB", "NR1SKINRB", "NR2P0ISRB", "NR2SKINRB", "OS1P0ISRB", "OS1SKINRB")
bam_fnames = str_c("/hpc/pmc_vanboxtel/processed/VAN4509/", samples, "/mapping/", samples, "_dedup.bam")

counts = loadAllData(bam_fnames, region = muts_2_call, q = q, mq = mq)
saveRDS(counts, "counts.rds")

vcf = do_shearwater_full(counts, samples, muts_2_call)


gr = rowRanges(vcf)
names(gr) = str_c("id_", seq_along(vcf))
rowRanges(vcf) = gr
writeVcf(vcf, "shear.vcf")
#The written out vcf does not contain chromosome lengths in its header.
#These need to be manually added.

vcf_snv = vcf[isSNV(vcf)]
writeVcf(vcf_snv, "shear_snv.vcf")
