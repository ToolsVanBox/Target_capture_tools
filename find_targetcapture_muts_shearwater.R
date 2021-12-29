library(tidyverse)
library(VariantAnnotation)
library(deepSNV)
library(broom)
library(ggbeeswarm)
library(ape)
library(nlme)
library(performance)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/basic_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_targeted/filter_exons_shear_functions.R")


#Perform multiple analyses in one go.
do_whole_shear_analysis = function(shear_vcf, prev_gr_l, fetuses, dir = "."){
    
    if(!dir.exists(dir)){
        dir.create(dir)
    }
    old_dir = setwd(dir)
    on.exit(setwd(old_dir), add = T)
    
    writeVcf(shear_vcf, "shearwater.vcf")
    
    #Look how many overlaps there are between the targeted and the original sequencing.
    overlaps_tb = create_overlaps_tb(fetuses, shear_vcf, prev_gr_l, overview_samples) %>% 
        write_tsv("overlaps_withknown.txt")
    
    #Compare the vafs between bulks
    purrr::walk(fetuses, function(fetus) analyze_vafs(fetus, shear_vcf, prev_vcfs, overview_samples))
    
    #Compare the mut counts between disomy and trisomy.
    compare_mutcounts_divstri(fetuses)
    
    invisible(0)
}

#Look how many overlaps there are between the targeted sequencing and the original
create_overlaps_tb = function(fetuses, shear_vcf, prev_gr_l, overview_samples){
    overlaps_tb = purrr::map(fetuses, count_overlaps_fetus, shear_vcf, prev_gr_l, overview_samples) %>% 
        set_names(fetuses) %>% 
        bind_rows(.id = "fetus")
    return(overlaps_tb)
}

#Count the number of overlaps between the targeted sequencing and the original for one fetus.
count_overlaps_fetus = function(fetus, shear_vcf, prev_gr_l, overview_samples){
    
    #Get capture mutations and existing mutations belonging to the fetus.
    nr_captured = nrow(shear_vcf)
    shear_vcf = filter_fetus_shearvcf(shear_vcf, fetus)
    nr_captured_fetus = nrow(shear_vcf)
    gr = granges(shear_vcf)
    prev_gr = prev_gr_l[[fetus]]
    nr_total = length(prev_gr)
    
    #Overlap mutations with existing muts.
    hits = findOverlaps_muts(gr, prev_gr)
    nr_overlap = length(hits)
    gr_occur = gr[queryHits(hits)]
    
    #find sample overviews
    overview_fetus = overview_samples %>% 
        dplyr::filter(fetus == !!fetus)
    
    #Count overlaps with different vcf sets.
    vcfs_shared_fnames = c(overview_fetus$shared_vcf[[1]], overview_fetus$indel_shared[[1]])
    nr_shared = count_overlaps(vcfs_shared_fnames, gr_occur)
    
    vcfs_unique_fnames = c(overview_fetus$unique_vcf, overview_fetus$indel_uniq)
    nr_uniq = count_overlaps(vcfs_unique_fnames, gr_occur)
    
    vcfs_bulk_fnames = c(overview_fetus$uniq_vcf_inbulk, overview_fetus$indel_uniq_inbulk, overview_fetus$indel_shared_inbulk[[1]], overview_fetus$shared_inbulk_vcf[[1]])
    nr_bulk = count_overlaps(vcfs_bulk_fnames, gr_occur)
    
    vcfs_notbulk_fnames = c(overview_fetus$uniq_vcf_notbulk, overview_fetus$indel_uniq_notbulk, overview_fetus$indel_shared_notbulk[[1]], overview_fetus$shared_notbulk_vcf[[1]])
    nr_notbulk = count_overlaps(vcfs_notbulk_fnames, gr_occur)
    
    tb = tibble("nr_total_old" = nr_total, "nr_captured" = nr_captured, "nr_captured_fetus" = nr_captured_fetus, "nr_overlap" = nr_overlap, "nr_shared" = nr_shared, "nr_uniq" = nr_uniq, "nr_bulk" = nr_bulk, "nr_notbulk" = nr_notbulk)
    return(tb)
}

#Find overlapping mutations between two granges objects.
findOverlaps_muts = function(gr, gr2){
    hits = findOverlaps(gr, gr2, type = "equal")
    same_alt = unlist(gr[queryHits(hits)]$ALT == gr2[subjectHits(hits)]$ALT)
    hits = hits[same_alt]
    return(hits)
}

count_overlaps = function(vcf_fnames, gr_occur){
    gr_shared = read_vcfs_as_gr(vcf_fnames)
    nr_overlap = findOverlaps_muts(gr_occur, gr_shared) %>% 
        length()
    return(nr_overlap)
}



read_vcfs_as_gr = function(vcf_fnames){
    vcf_fnames = vcf_fnames[file.exists(vcf_fnames)]
    vcfs = purrr::map(vcf_fnames, readVcf, genome = "hg19")
    gr = purrr::map(vcfs, granges) %>% 
        do.call(c, .)
    return(gr)
}

remove_nonpresent_variants = function(vcf){
    gt = geno(vcf)$GT
    gt = merge_gt(gt)
    gt_f = str_detect(gt, "0/1|1/1")
    vcf_present = vcf[gt_f]
    return(vcf_present)
}

#Merge the genotype columns into a single vector.
merge_gt = function(gt){
    gt = gt %>% 
        as.data.frame %>%
        as_tibble() %>% 
        unite(., "gt", 1:ncol(.), sep = ";") %>% 
        pull(gt)
    return(gt)
}



#Function to analyze the vafs per fetus.
analyze_vafs = function(fetus, shear_vcf, prev_vcfs, overview_samples){
    
    if(!dir.exists(fetus)){
        dir.create(fetus)
    }
    old_dir = setwd(fetus)
    on.exit(setwd(old_dir), add = T)
    
    #Get data from fetus
    shear_vcf = filter_fetus_shearvcf(shear_vcf, fetus)
    gr = granges(shear_vcf)
    prev_vcf = prev_vcfs[[fetus]]
    prev_gr = granges(prev_vcf)
    
    #Overlap mutations with existing muts.
    hits = findOverlaps_muts(gr, prev_gr)
    vcf_occur = shear_vcf[queryHits(hits)]
    prev_vcf = prev_vcf[subjectHits(hits)]
    writeVcf(prev_vcf, "wgs_captured.vcf")
    writeVcf(vcf_occur, str_c(fetus, "_shear.vcf"))
    
    #Get gts from prev_vcf. (These show in what branches the muts were originally.)
    overview_fetus = overview_samples %>% 
        dplyr::filter(fetus == !!fetus)
    gt = get_gt_prev_merge(prev_vcf, overview_fetus)
    
    #Calculate whether muts have the same vaf beween bulks
    bulks_all = samples(header(vcf_occur))
    bulks = str_subset(bulks_all, fetus)
    ref = get_ref_shearvcf(vcf_occur)[,bulks, drop = F]
    alt = get_alt_shearvcf(vcf_occur)[,bulks, drop = F]
    
    combine_ref_alt_col = function(ref_col, alt_col){
        purrr::map2(ref_col, alt_col, function(x, y) as.integer(c(x, y)))
    }
    ad = purrr::map(seq(1, ncol(ref)), function(i) combine_ref_alt_col(ref[,i], alt[,i])) %>% 
        do.call(cbind, .)
    colnames(ad) = bulks
    same_vaf_tb = purrr::map(seq(1, nrow(ad)), function(i) test_same_vaf_bulks(ad[i,])) %>% 
        bind_rows(.id = "id") %>% 
        write_tsv("vaf_betweenbulks.txt")
    
    
    #Calculate whether vafs are the same within a branch. This is done per bulk.
    same_vaf_muts_tb = test_same_vaf_muts(ad, gt) %>% 
        write_tsv("vaf_betweenmuts_inbranch.txt")
    
    
    #Calculate the vafs and plot them per branch.
    vaf_gt_tb = make_vaf_gt_tb(vcf_occur, gt, bulks)
    write_tsv(vaf_gt_tb, "vaf_gt_tb.txt")
    vaf_fig = plot_vaf_bulks(vaf_gt_tb)
    ggsave("vafs.pdf", vaf_fig, width = 40)
    
    #Plot tree with only captured mutations:
    plot_tree(prev_vcf, overview_fetus)
    
    invisible(0)
}

#Function to plot a tree for a fetus.
plot_tree = function(vcf, overview_fetus){
    
    gt = get_gt_prev(vcf, overview_fetus)
    co_occur_m = cbind(gt, "root" = c(rep(0, nrow(gt))))
    max_muts = colSums(co_occur_m) %>% max()
    
    #Create tree
    tree = co_occur_m %>% t() %>% dist.gene() %>% nj() #neighbour joining tree construction
    rooted_tree = root(tree, outgroup = "root", resolve.root = T) %>% drop.tip("root", trim.internal = T)
    
    #Color edges based on tissue
    #color_map = tibble("origin" = c("Liver", "SI", "BM", "Intestine"), "color" = c("Red", "Darkblue", "Darkred", "Blue"))
    #colors_tb = left_join(overview_fetus[,c("origin", "sample")], color_map, by = "origin")
    #colors_tb_ordered = enframe(clone_names, name = "nr", value = "sample") %>% inner_join(colors_tb, by = "sample") #Switch color ordering from the order in overview_samples to the order of cols in the vcf.
    #color_index = match(colors_tb_ordered$nr, rooted_tree$edge[,2])
    #tree_edge_cols = rep("black", nrow(rooted_tree$edge))
    #tree_edge_cols[color_index] = colors_tb_ordered$color
    
    #Plot tree
    pdf("tree.pdf")
    my_tree_plot = plot(rooted_tree, use.edge.length = T, label.offset = 1, edge.width = 2)
    axis(1,at =seq(from =0, to = max_muts, by = 1),cex.axis = 1)
    #nodelabels() # we id node label
    edgelabels(rooted_tree$edge.length, bg="black", col="white", font=2)
    dev.off()
    
    invisible(0)
}

#Function to test whether mutations have the same vaf between bulks
test_same_vaf_bulks = function(ads_mut){
    ref = purrr::map_int(ads_mut, 1)
    alt = purrr::map_int(ads_mut, 2)
    mat = rbind(ref, alt)
    sample_names = str_c(names(ads_mut), collapse = ";")
    prop_res = prop.test(mat[2,], n = colSums(mat))
    tb = broom::tidy(prop_res) %>% 
        tidyr::unite(!!sym(sample_names), contains("estimate"), sep = ";")
    return(tb)
}

#Function to test whether mutations wihtin a tree branch have the same vaf.
#This is tested for all branches in all bulks
test_same_vaf_muts = function(ad, gt){
    ad_branch_l = purrr::map(levels(gt), function(gt_level) ad[gt_level == gt, , drop = F])
    gt_nr_tb = tibble(gt = levels(gt), "nr_gt" = seq_along(gt))
    
    tb = purrr::map(ad_branch_l, test_same_vaf_mut) %>% 
        set_names(levels(gt)) %>% 
        bind_rows(.id = "gt") %>% 
        dplyr::filter(str_detect(method, "Pearson's")) %>% 
        left_join(gt_nr_tb, by = "gt")
    return(tb)
}

#Helper function to test whether mutations within a single branch have the same vaf within the bulks
test_same_vaf_mut = function(ad_branch){
    tb = purrr::map(colnames(ad_branch), function(bulk_name) test_same_vaf_mut_bulk(ad_branch[,bulk_name])) %>% 
        set_names(colnames(ad_branch)) %>% 
        bind_rows(.id = "bulk")
    return(tb)
}

#Helper function to test whether mutations within a single branch have the same vaf in a specific bulk
test_same_vaf_mut_bulk = function(ad_mut_branch){
    mat = do.call(rbind, ad_mut_branch)
    chisq_res = chisq.test(mat, simulate.p.value = T)
    tb = broom::tidy(chisq_res) %>% 
        dplyr::mutate(nr_muts = nrow(mat))
    
    # prop_res = prop.test(mat)
    # tb = broom::tidy(prop_res) %>% 
    #     tidyr::unite("estimate", contains("estimate"), sep = ";")
    return(tb)
}

#Make a table containing the vafs and the tree branches.
make_vaf_gt_tb = function(vcf, gt, bulks){
    vafs = geno(vcf)$VF[,bulks, drop = F]
    
    vaf_gt_tb = vafs %>%
        as.data.frame() %>% 
        dplyr::mutate(id = row_number(), gt = gt) %>%
        dplyr::arrange(desc(!!sym(bulks[1]))) %>% 
        dplyr::mutate(gt = factor(gt, levels = unique(gt))) %>%
        dplyr::group_by(gt) %>% 
        dplyr::mutate(vaf_order = row_number()) %>%
        dplyr::ungroup() %>% 
        gather(key = "sample", value = "vaf", -id, -vaf_order, -gt)
    return(vaf_gt_tb)
}


#Get the genotypes from the previously called mutations. (These show the tree branches.)
get_gt_prev = function(vcf, overview_fetus){
    bulks = c(overview_fetus$bulk[1], overview_fetus$other_bulks[1])
    gt = geno(vcf)$GT
    bulkcols = which(colnames(gt) %in% bulks)
    gt = gt[,-bulkcols, drop = F]
    gt[gt == "1/1" | gt == "0/1"] = 1 #Makes sure homozygous mutations on the x chromosome are placed in the same combi as muts on the autosomes.
    gt[gt == "0/0"] = 0
    gt = apply(gt, 2, as.integer)
    return(gt)
}
get_gt_prev_merge = function(vcf, overview_fetus){
    gt = get_gt_prev(vcf, overview_fetus)
    gt = gt %>% 
        merge_gt() %>% 
        factor()
    return(gt)
}

#Plot the vaf for all bulks, per tree branch.
plot_vaf_bulks = function(vaf_gt_tb){
    vaf_fig = ggplot(vaf_gt_tb, aes(x = vaf_order, y = vaf, colour = sample)) +
        geom_point() +
        geom_line() +
        facet_grid(. ~ gt, scales = "free_x") +
        labs(x = "") +
        theme_freek() +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank())
    return(vaf_fig)
}


get_ref_shearvcf = function(shear_vcf){
    alt = get_alt_shearvcf(shear_vcf)
    gen = geno(shear_vcf)
    dp = gen$FD + gen$BD
    ref = dp - alt
    return(ref)
}
get_alt_shearvcf = function(shear_vcf){
    gen = geno(shear_vcf)
    alt = gen$FW + gen$BW
    return(alt)
}

#Count the number of mutations in a prev vcf.
count_prev_muts = function(fetus){
    vcf_fname = str_c(fetus, "/wgs_captured.vcf")
    prev_vcf = readVcf(vcf_fname)
    overview_fetus = overview_samples %>% 
        dplyr::filter(fetus == !!fetus)
    
    gt = get_gt_prev(prev_vcf, overview_fetus)
    count_tb = gt %>% 
        colSums() %>% 
        enframe(value = "nrmuts", name = "sample") %>% 
        dplyr::mutate(fetus = fetus,
                      celltype = ifelse(str_detect(sample, "SI|Si"), "SI", "HSPC"))
    return(count_tb)
}

get_fixed = function(m, name = NA){
    lme_sum = summary(m)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    ci = intervals(m, which = "fixed")$fixed %>% as_tibble()
    vals = cbind(vals, ci[,c("lower", "upper")])
    call = lme_sum$call %>% as.character() %>% str_c(collapse = " ")
    invisible(capture.output(r2s <- suppressWarnings(performance::r2(m))))
    if (is_na(r2s)){
        r2s = list("R2_conditional" = NA ,"R2_marginal" = NA)
    }
    vals = vals %>% 
        dplyr::mutate("r2_cond" = r2s$R2_conditional, "r2_marg" = r2s$R2_marginal, logLik = m$logLik, BIC = lme_sum$BIC, AIC = lme_sum$AIC, call = call, name = name) %>% 
        as_tibble()
    return(vals)
}


compare_mutcounts_divstri = function(fetuses){
    count_tb = map(fetuses, count_prev_muts) %>% 
        bind_rows() %>% 
        dplyr::mutate(trisomy = ifelse(fetus %in% c("N01", "OS1"), "T21", "D21"),
                      trisomy = factor(trisomy, levels = c("D21", "T21")),
                      celltype = factor(celltype, levels = c("HSPC", "SI")))
    
    m1 = lme(nrmuts ~ trisomy*celltype, random = ~ 1 | fetus, data = count_tb, weights = varIdent(form = ~1|trisomy))
    m1_tb = get_fixed(m1)
    pval = m1_tb %>%
        dplyr::filter(variable == "trisomyT21") %>% 
        pull(`p-value`)
    
    no_age_fig = ggplot(count_tb, aes(y = nrmuts, x = trisomy, colour = fetus, shape = celltype)) +
        geom_boxplot(outlier.shape = NA, colour = "black", shape = "") +
        geom_quasirandom(size = 3) +
        labs(title = "DivsTri no age", x = "", y = "Nr. of substitutions", colour = "Fetus", shape = "Celltype") +
        annotate("text", x = 1, y = 23, label = paste0("P (D21 vs T21): ", round(pval, 3)), size = 6) +
        theme_freek()
    ggsave("mut_counts.pdf", no_age_fig, useDingbats = F)
}




fetuses = c("NR1", "NR2", "OS1", "N01")
setwd("~/surfdrive/Shared/Projects/Freek/Freek_trees_targeted/shear_strict/")

#Read previously called mutations
overview_samples = read_tsv("~/hpc/pmc_vanboxtel/projects/Freek_trees/overview_samples.txt")
prev_vcf_fnames = str_c("~/hpc/pmc_vanboxtel/projects/Freek_trees/", fetuses, "/", fetuses, "_complete.vcf")
prev_vcfs = purrr::map(prev_vcf_fnames, readVcf, genome = "hg19")
names(prev_vcfs) = fetuses
prev_gr_l = purrr::map(prev_vcfs, granges)
prev_gr_l = purrr::map2(prev_gr_l, fetuses, function(gr, fetus) {gr$fetus = fetus; return(gr)})


#Set sample names
#samples = c("N01P0ISRB", "N01SKINRB", "N01SPLEENRB", "NR1P0ISRB", "NR1SKINRB", "NR2P0ISRB", "NR2SKINRB", "OS1P0ISRB", "OS1SKINRB")
#muts_2_call = read_tsv("muts.bed", col_types = "cii", col_names = c("Chrom", "Start", "End"), comment = "#") %>%
#muts_2_call = read_tsv("muts.bed", col_types = "cii") %>%
#    makeGRangesFromDataFrame(starts.in.df.are.0based = T) %>% 
#    sort()

#Read shearwater vcf
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
vcf_fname = "~/surfdrive/Shared/Projects/Freek/Freek_trees_targeted/shear_strict/hpc_output/shear.vcf"
vcf = readVcf(vcf_fname)
shear_vcf = vcf[isSNV(vcf)]
indel_vcf = vcf[isIndel(vcf)]


#Remove variants near indels
shear_vcf = remove_variants_near_indels(shear_vcf, indel_vcf)

#Filter on quality
#shear_vcf = filter_quality(shear_vcf)

#Filter on highvaf others
#shear_vcf = filter_highvaf_others(shear_vcf)

writeVcf(shear_vcf, "filtered_shear.vcf")

#Do analyses
do_whole_shear_analysis(shear_vcf, prev_gr_l, fetuses)

#Filter the shear vcf for mutations that are consistently present in the bulks.
shear_vcf_consistent = filter_consistent_shearvcf(shear_vcf)
do_whole_shear_analysis(shear_vcf_consistent, fetuses, dir = "consistent_in_bulk_shear")



#Try different cutoffs.
cutoffs = c(0.01, 0.05, 0.1, 0.2, 0.5)
shear_vcf_l = purrr::map(cutoffs, function(cutoff) do_shearwater_full(counts, samples, muts_2_call, cutoff))
overlaps_cutoffs_tb = purrr::map(shear_vcf_l, function(shear_vcf) create_overlaps_tb(fetuses, shear_vcf, prev_gr_l, overview_samples)) %>% 
    set_names(cutoffs) %>% 
    bind_rows(.id = "cutoff")
write_tsv(overlaps_cutoffs_tb, "overlaps_cutoffs.txt")

#Plot how the different cutoffs affect the overlaps
overlaps_cutoffs_tb %>%
    dplyr::select(cutoff, fetus, nr_overlap)
overlap_ratio_fig = ggplot(overlaps_cutoffs_tb, aes(x = cutoff, fill = fetus, y = nr_overlap/nr_captured_fetus)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_freek()
ggsave("overlap_ratio.pdf", overlap_ratio_fig)

overlap_fig = ggplot(overlaps_cutoffs_tb, aes(x = cutoff, fill = fetus, y = nr_overlap)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_freek()
ggsave("overlap.pdf", overlap_fig)

#Cluster samples based on the raw counts
cluster_counts = function(counts, samples){
    counts_2d = purrr::map(seq(1, dim(counts)[3]), function(i) counts[,,i]) %>% 
        do.call(cbind, .)
    rownames(counts_2d) = samples
    counts_cluster = counts_2d %>% 
        t() %>% 
        scale() %>% 
        t() %>% 
        dist() %>% 
        hclust() 
    return(counts_cluster)
}

pdf("counts_cluster.pdf")
plot(cluster_counts(counts, samples))
dev.off()

pvals = betabinLRT(counts, rho=1e-4, maxtruncate = 1, truncate = 1)$pvals
qvals = p.adjust(pvals, method="BH")
dim(qvals) = dim(pvals)

mut_qvals = purrr::map_dbl(seq(1, dim(counts)[2]), function(i) min(qvals[,i,], na.rm = T))
mut_f = mut_qvals < 0.05

counts_mut = counts[,mut_f,]
counts_nonmut = counts[,!mut_f,]

pdf("counts_cluster_mut.pdf")
plot(cluster_counts(counts_mut, samples))
dev.off()

pdf("counts_cluster_nonmut.pdf")
plot(cluster_counts(counts_nonmut, samples))
dev.off()



#Look at the vafs between matching and not matching mutations
get_vaf_overlapvsnot = function(bulk, shear_vcf, prev_gr_l){
    #Get capture mutations and existing mutations belonging to the fetus.
    shear_vcf = filter_fetus_shearvcf(shear_vcf, bulk)
    gr = granges(shear_vcf)
    fetus = str_sub(bulk, 1, 3)
    prev_gr = prev_gr_l[[fetus]]
    
    #Overlap mutations with existing muts.
    hits = findOverlaps_muts(gr, prev_gr)
    
    #Determine vafs for overlapping and non overlapping mutation calls
    vaf_match = geno(shear_vcf[queryHits(hits)])$VF[,bulk]
    vaf_nomatch = geno(shear_vcf[-queryHits(hits)])$VF[,bulk]
    vaf_match = tibble("bulk" = bulk, "vaf" = vaf_match, type = "match")
    vaf_nomatch = tibble("bulk" = bulk, "vaf" = vaf_nomatch, type = "nomatch")
    vaf = rbind(vaf_match, vaf_nomatch)
    return(vaf)
}

vafs = purrr::map(samples, function(bulk) get_vaf_overlapvsnot(bulk, shear_vcf, prev_gr_l)) %>% 
    bind_rows()

vaf_fig = ggplot(vafs, aes(x = type, y = vaf, colour = bulk)) +
    geom_quasirandom() +
    theme_freek()
ggsave("vaf_matchvsnomatch_noindels.pdf", vaf_fig)



# get_shearwater_vcf = function(sample, samples, counts, muts_2_call){
#     samples_sub = select_samples(sample, samples)
#     counts_sub = counts[samples %in% samples_sub,,]
#     vcf_sample = do_shearwater(counts_sub, sample, samples_sub, muts_2_call)
#     return(vcf_sample)
# }
# select_samples = function(sample, samples){
#     fetus = str_sub(sample, 1, 3)
#     samples_other_fetus = samples %>% str_subset(fetus, negate = T)
#     samples_keep = c(sample, samples_other_fetus)
#     return(samples_keep)
# }
# 
# do_shearwater = function(counts_sub, sample, samples_sub, muts_2_call){
#     pvals <- betabinLRT(counts_sub, rho=1e-4, maxtruncate = 1, truncate = 1)$pvals
#     qvals <- p.adjust(pvals, method="BH")
#     dim(qvals) = dim(pvals)
#     vcf = qvals2Vcf(qvals, counts_sub, muts_2_call, samples = samples_sub, mvcf = TRUE, cutoff = 0.05)
#     vcf_sample = vcf[geno(vcf)$GT[,sample] == 1]
#     return(vcf_sample)
# }
# 
# shear_vcf_l = purrr::map(samples, function(sample) get_shearwater_vcf(sample, samples, counts, muts_2_call)) %>%
#     set_names(samples)

# calc_ref_alt_dp = function(shear_vcf){
#     gen = geno(shear_vcf)
#     alt = gen$FW + gen$BW
#     dp = gen$FD + gen$BD
#     ref = dp - alt
#     geno(shear_vcf)$alt = alt
#     geno(shear_vcf)$ref = ref
#     }
#
# merge_vcfs = function(fetus, shear_vcf_l){
#     fetus_samples = str_subset(names(shear_vcf_l), fetus)
#     vcf_l = shear_vcf_l[fetus_samples]
#     gr_l = purrr::map(vcf_l, get_shear_gr)
#     findOverlaps_muts(gr_l[[1]], gr_l[[2]])
#
# }
#
# merge_2_grs = function(gr1, gr2){
#     hits = findOverlaps_muts(gr1, gr2)
#     gr_total = c(gr1, gr2) %>%
#         unique()
#     gr_total$ref_dp = NA
#
# }
#
# get_shear_gr = function(shear_vcf){
#     gr = granges(shear_vcf)
#     gr$ref_dp = get_ref_shearvcf(shear_vcf)
#     gr$alt_dp = get_alt_shearvcf(shear_vcf)
#     gr$vaf = geno(shear_vcf)$VF[,1]
#     return(gr)
# }

# quick_count_overlaps_fetus = function(fetus, shear_vcf, prev_gr_l, overview_samples){
#     
#     #Get capture mutations and existing mutations belonging to the fetus.
#     nr_captured = nrow(shear_vcf)
#     shear_vcf = filter_fetus_shearvcf(shear_vcf, fetus)
#     nr_captured_fetus = nrow(shear_vcf)
#     gr = granges(shear_vcf)
#     prev_gr = prev_gr_l[[fetus]]
#     nr_total = length(prev_gr)
#     
#     #Overlap mutations with existing muts.
#     hits = findOverlaps_muts(gr, prev_gr)
#     nr_overlap = length(hits)
#     return(nr_overlap)
# }









# gr_n01 = granges(vcfML_n01)
# 
# hits = findOverlaps_muts(gr_n01, prev_n01)
# vcfML_n01 = vcfML_n01[queryHits(hits)]
# geno(vcfML_n01)$VF
# 
# #bayesian
# bf = bbb(counts, model = "AND", truncate = 0.05)
# vcf = bf2Vcf(bf, counts, muts_2_call, cutoff = 0.05, samples = samples, prior = 0.5, mvcf = T)
# 
# #ML
# pvals <- betabinLRT(counts, rho=1e-4, maxtruncate = 1, truncate = 1)$pvals
# qvals <- p.adjust(pvals, method="BH")
# dim(qvals) = dim(pvals)
# vcfML = qvals2Vcf(qvals, counts, muts_2_call, samples = samples, mvcf = TRUE)
# vcfML = vcfML[geno(vcfML)$GT[,1] == 1]
# gr_fullcounts = granges(vcfML)
# 
# findOverlaps_muts(gr_n01, prev_n01)
# findOverlaps_muts(gr_fullcounts, prev_n01)
