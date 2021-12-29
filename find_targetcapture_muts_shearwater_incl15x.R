library(tidyverse)
library(VariantAnnotation)
library(deepSNV)
library(broom)
library(ggbeeswarm)
library(ape)
library(nlme)
library(performance)
source("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/basic_functions.R")
source("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Boxtel_General/Scripts/pmc_vanboxtel/Freek_targeted/filter_exons_shear_functions.R")


#Perform multiple analyses in one go.
do_whole_shear_analysis = function(shear_vcf, prev_vcfs, tree_tables_l, fetuses, dir = "."){
    
    if(!dir.exists(dir)){
        dir.create(dir)
    }
    old_dir = setwd(dir)
    on.exit(setwd(old_dir), add = T)
    
    writeVcf(shear_vcf, "shearwater.vcf")
    
    #Compare the vafs between bulks
    purrr::walk(fetuses, function(fetus) analyze_vafs(fetus, shear_vcf, prev_vcfs, tree_tables_l))

    invisible(0)
}




#Find overlapping mutations between two granges objects.
findOverlaps_muts = function(gr, gr2){
    hits = findOverlaps(gr, gr2, type = "equal")
    same_alt = unlist(gr[queryHits(hits)]$ALT == gr2[subjectHits(hits)]$ALT)
    hits = hits[same_alt,]
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
analyze_vafs = function(fetus, shear_vcf, prev_vcfs, tree_tables_l){
    
    if(!dir.exists(fetus)){
        dir.create(fetus)
    }
    old_dir = setwd(fetus)
    on.exit(setwd(old_dir), add = T)
    
    #Get data from fetus
    shear_vcf = filter_fetus_shearvcf(shear_vcf, fetus)
    gr = granges(shear_vcf)
    
    
    #Get gt
    tree_table = tree_tables_l[[fetus]]
    shear_tb = data.frame("chr" = as.vector(seqnames(gr)), "HG38" = start(gr), "i" = seq_len(length(gr)))
    gt_join = inner_join(shear_tb, tree_table, by = c("chr", "HG38"))
    gt = gt_join %>% 
        dplyr::select(-chr, -HG38, -i) %>% 
        as.matrix()
    
    vcf_occur = shear_vcf[gt_join$i,]
 
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
    
    

    
    gt_merged = gt %>% 
        merge_gt() %>% 
        factor()
    
    #Calculate whether vafs are the same within a branch. This is done per bulk.
    same_vaf_muts_tb = test_same_vaf_muts(ad, gt_merged) %>% 
        write_tsv("vaf_betweenmuts_inbranch.txt")
    
    #Calculate the vafs and plot them per branch.
    vaf_gt_tb = make_vaf_gt_tb(vcf_occur, gt_merged, bulks)
    write_tsv(vaf_gt_tb, "vaf_gt_tb.txt")
    vaf_fig = plot_vaf_bulks(vaf_gt_tb)
    ggsave("vafs.pdf", vaf_fig, width = 40)
    
    #Plot tree with only captured mutations:
    plot_tree(gt)
    
    invisible(0)
}

#Function to plot a tree for a fetus.
plot_tree = function(gt){
    
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
    pdf("tree.pdf", width = 12)
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
    gt = geno(vcf)$GT
    gt[gt == "1/1" | gt == "0/1"] = 1 #Makes sure homozygous mutations on the x chromosome are placed in the same combi as muts on the autosomes.
    gt[gt == "0/0" | gt == "."] = 0
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
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(), text = element_text(size = 7))
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
setwd("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/Freek_trees_targeted/shear_strict_incl15x_samples/")

#Read previously called mutations
overview_samples = read_tsv("overview_samples.txt")
prev_vcf_fnames = str_c("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/Freek_trees_targeted/shear_strict_incl15x_samples/wgs_vcfs/", fetuses, "_only_complete.vcf")
prev_vcfs = purrr::map(prev_vcf_fnames, readVcf, genome = "hg38")
names(prev_vcfs) = fetuses
# prev_gr_l = purrr::map(prev_vcfs, granges)
# prev_gr_l = purrr::map2(prev_gr_l, fetuses, function(gr, fetus) {gr$fetus = fetus; return(gr)})




tree_table_fnames = str_c("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/Freek_trees_targeted/shear_strict_incl15x_samples/wgs_vcfs/", fetuses, "_tree_table_hg38.txt")

read_tree_tables = function(tree_table_fname){
    tree_table = read.delim(tree_table_fname, sep = " ")
    tree_table$chr = str_remove(rownames(tree_table), ":.*")
    return(tree_table)
}

tree_tables_l = purrr::map(tree_table_fnames, read_tree_tables)
names(tree_tables_l) = fetuses

#Read shearwater vcf
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
vcf_fnames = list.files("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/Freek_trees_targeted/shear_strict_incl15x_samples/shear_vcfs/", 
                        full.names = TRUE, pattern = ".vcf$")

shear_vcf = purrr::map(vcf_fnames, readVcf) %>% 
    do.call(rbind, .) %>% 
    sort()

writeVcf(shear_vcf, "filtered_shear.vcf")

do_whole_shear_analysis(shear_vcf, prev_vcfs, tree_tables_l, fetuses, dir = "consistent_in_bulk_shear")


