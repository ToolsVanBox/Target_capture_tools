add_filler_filter_header = function(vcf){
    DF_ori = fixed(header(vcf))[[1]]
    DF = DF_ori[1, , drop = F]
    rownames(DF) = "Filler"
    DF$Description = "This is a filler, to prevent igv complaints"
    DF_new = rbind(DF_ori, DF)
    fixed(header(vcf))[[1]] = DF_new
    return(vcf)
}


#Filter for variants occurring in both bulks. (excluding spleen)
filter_consistent_shearvcf = function(shear_vcf){
    gt = geno(shear_vcf)$GT
    gt[gt == "."] = 0
    gt = apply(gt, 2, as.numeric)
    spleen_f = str_detect(colnames(gt), "SPLEEN", negate = T)
    gt = gt[,spleen_f, drop = F]
    fetuses = str_sub(colnames(gt), 1, 3) %>% 
        unique()
    
    #Check per fetus whether the mutations are consistenlty present in the bulk
    nr_consistent_fetuses = purrr::map(fetuses, function(fetus) filter_consistent_fetus_shearvcf(fetus, gt)) %>% 
        bind_cols() %>% 
        rowSums()
    filtered_vcf = shear_vcf[nr_consistent_fetuses >= 1]
    return(filtered_vcf)
}
#Helper function of filter_consistent_shearvcf
filter_consistent_fetus_shearvcf = function(fetus, gt){
    gt = gt[,str_detect(colnames(gt), fetus)]
    nrbulks_found = rowSums(gt)
    consistent_infetus_f = nrbulks_found == 2
    return(consistent_infetus_f)
}

#Filter for variants that occur in only a single tissue, instead of in both tissues.
filter_singletissue_shearvcf = function(shear_vcf){
    gt = get_gt(shear_vcf)
    fetuses = str_sub(colnames(gt), 1, 3) %>% 
        unique()
    
    #Check per fetus whether the mutations are consistenlty present in the bulk
    nr_consistent_fetuses = purrr::map(fetuses, function(fetus) filter_singletissue_fetus_shearvcf(fetus, gt)) %>% 
        bind_cols() %>% 
        rowSums()
    filtered_vcf = shear_vcf[nr_consistent_fetuses >= 1]
    return(filtered_vcf)
}
#Helper function of filter_singletissue_shearvcf
filter_singletissue_fetus_shearvcf = function(fetus, gt){
    gt = gt[,str_detect(colnames(gt), fetus)]
    nrbulks_found = rowSums(gt)
    singletissue_f = nrbulks_found == 1
    return(singletissue_f)
}

#Filter for variants that have the correct REF
filter_correct_ref = function(vcf, ref_genome){
    gr = granges(vcf)
    seqlevelsStyle(gr) = "UCSC"
    ref_seq = getSeq(get(ref_genome), seqnames(gr), start(gr), end(gr))
    ref_f = ref_seq == gr$REF
    vcf = vcf[ref_f]
    return(vcf)
}

#Filter for non germline variants by comparing against the original non-targeted vcf.
filter_not_germline = function(shear_vcf, prev_vcfs, fetuses, type = c("all", "indel")){
    somatic_vcf = purrr::map(fetuses, filter_not_germline_fetus, shear_vcf, prev_vcfs, type) %>% 
        do.call(rbind, .) %>% 
        sort() %>% 
        unique()
    return(somatic_vcf)
}

#Helper function of filter_not_germline that works on a single fetus.
#Filter for non germline variants by comparing against the original non-targeted vcf.
filter_not_germline_fetus = function(fetus, shear_vcf, prev_vcfs, type = c("all", "indel")){
    
    type = match.arg(type)
   
    #Get shear gr
    fetus_shear_vcf = filter_fetus_shearvcf(shear_vcf, fetus)
    gr = granges(fetus_shear_vcf)
    
    #Get wgs vcf and get germline variants
    prev_vcf = prev_vcfs[[fetus]]
    gt_prev = get_gt(prev_vcf)
    germ_f = rowSums(gt_prev) >= ncol(gt_prev) - 1
    prev_vcf = prev_vcf[germ_f]
    
    if (type == "indel"){
        prev_vcf = prev_vcf[isIndel(prev_vcf)]
    }
    
    prev_gr = granges(prev_vcf)
    
    if (nrow(prev_vcf) == 0){
        return(fetus_shear_vcf)
    }
    
    #Overlap to find variants that were already found with the wgs
    if (type == "all"){
    hits = findOverlaps_muts(prev_gr, gr)
    } else if (type == "indel"){
        hits = find_nearby_indels(prev_gr, gr)
    }
    germ_i = subjectHits(hits)
    
    if (length(germ_i)){
        somatic_fetus_shear_vcf = fetus_shear_vcf[-germ_i]
    } else{
        somatic_fetus_shear_vcf = fetus_shear_vcf
    }
    return(somatic_fetus_shear_vcf)
}

#Helper function of filter_not_germline_fetus. Filters for variants in a specific fetus.
filter_fetus_shearvcf = function(shear_vcf, fetus){
    gt = get_gt(shear_vcf)
    sample_f = str_detect(colnames(gt), fetus)
    gt = gt[,sample_f, drop = F]
    mut_f = rowSums(gt) > 0
    fetus_vcf = shear_vcf[mut_f]
    return(fetus_vcf)
}

#Function to find variants present in two grs
findOverlaps_muts = function(gr, gr2){
    hits = findOverlaps(gr, gr2, type = "equal")
    same_alt = any(gr[queryHits(hits)]$ALT == gr2[subjectHits(hits)]$ALT)
    hits = hits[same_alt]
    return(hits)
}

#Function to find variants near an indel
find_nearby_indels = function(indel_gr, gr){
    indel_gr = GenomicRanges::reduce(indel_gr)
    hits = findOverlaps(gr, indel_gr, maxgap = 10)
    return(hits)
}


#Function to get th gt of a vcf
get_gt = function(vcf){
    gt = geno(vcf)$GT
    #gt[gt == "1/1" | gt == "0/1"] = 1 #Makes sure homozygous mutations on the x chromosome are placed in the same combi as muts on the autosomes.
    gt[gt == "0/0" | gt == "."] = 0
    gt[gt != 0] = 1
    
    #Transform to integer. Also works for single row gt.
    col_names = colnames(gt)
    row_names = rownames(gt)
    gt = purrr::map(seq_len(nrow(gt)), function(i) as.integer(gt[i,])) %>% 
        do.call(rbind, .)
    colnames(gt) = col_names
    rownames(gt) = row_names
    return(gt)
}

#Functions to filter out variants occuring in multiple fetuses
filter_shear_multiple_fetuses = function(shear_vcf, fetuses){
    gt = get_gt(shear_vcf)
    pos_i = purrr::map(fetuses, get_fetus_i, gt) %>% 
        do.call(c, .)
    multiple_fetus_i = unique(pos_i[duplicated(pos_i)])
    unique_fetus_shear_vcf = shear_vcf[-multiple_fetus_i]
    return(unique_fetus_shear_vcf)
}
#Helper function of filter_shear_multiple_fetuses
get_fetus_i = function(fetus, gt){
    sample_f = str_detect(colnames(gt), fetus)
    gt = gt[, sample_f, drop = F]
    present_i = which(rowSums(gt) > 0)
    return(present_i)
}

#Function to create some graphs and a list of genes.
create_graphs  = function(vcf, out_dir, ref_genome){
    
    #Go to out dir
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = T)
    
    #Look at genes.
    ann = info(vcf)$ANN %>%
        as.list()
    ann[elementNROWS(ann) == 0] = ""
    ann = purrr::map_chr(ann, str_c, collapse = ";")
    ann_l = str_split(ann, pattern = "\\|")
    gene = purrr::map_chr(ann_l, 4)
    
    uniq_gene = unique(gene)
    gene_tb = tibble("gene" = uniq_gene) %>% 
        write_tsv("uniq_genes.txt")
    
    #Create overview table
    type = purrr::map_chr(ann_l, 2)
    effect_size = purrr::map_chr(ann_l, 3)
    cosm_ids = str_remove(names(vcf), "id_[0-9]*;?")
    gt = get_gt(vcf)
    called_samples = purrr::map_chr(seq(1, nrow(gt)), ~get_samples_gt(gt[.x,]))
    
    
    gr = granges(vcf)
    gr$ALT = gr$ALT %>% 
        unlist() %>% 
        as.vector()
    gr_tb = gr %>% 
        as.data.frame() %>% 
        as_tibble() %>% 
        dplyr::select(chromosome = seqnames, start, end, REF, ALT) %>% 
        dplyr::mutate(gene = gene, type = type, effect_size = effect_size, cosm_ids = cosm_ids, called_samples = called_samples)
    
    vafs = geno(vcf)$VF
    colnames(vafs) = str_c("vaf_", colnames(vafs))
    colnames(gt) = str_c("gt_", colnames(gt))
    vaf_gt_tb = cbind(vafs, gt) %>% 
        as.data.frame() %>% 
        as_tibble()
    
    overview_tb = bind_cols(gr_tb, vaf_gt_tb)
    write_tsv(overview_tb, "overview.txt")
    
    #Look at nr muts
    nr_muts_tb = vcf %>%
        get_gt() %>% 
        colSums() %>% 
        enframe(name = "sample", value = "nr_muts")
    
    nr_muts_fig = ggplot(nr_muts_tb, aes(x = sample, y = nr_muts)) +
        geom_bar(stat = "identity") +
        labs(y = "nr variants") +
        theme_freek() +
        theme(axis.text.x = element_text(angle = 90))
    ggsave("nr_muts.pdf", nr_muts_fig, useDingbats = F)
    
    #Look at vaf
    vaf_tb = geno(vcf)$VF %>% 
        as.data.frame() %>% 
        as_tibble() %>% 
        gather(key = "bulk", value = "vaf")
    
    vaf_fig = ggplot(vaf_tb, aes(x = bulk, y = vaf)) +
        geom_quasirandom(size = 0.5) +
        theme_freek() +
        theme(axis.text.x = element_text(angle = 90))
    ggsave("vafs.pdf", vaf_fig, useDingbats = F)
    
    
    seqlevelsStyle(gr) = "UCSC"
    type_occur = mut_type_occurrences(GRangesList(gr), ref_genome)
    spec_fig = plot_spectrum(type_occur)
    ggsave("spectrum.pdf", spec_fig, useDingbats = F)
    
    invisible(0)
    
    
}

#Helper function of create_graphs
get_samples_gt = function(gt_mut){
    gt_f = gt_mut == 1
    present_samples = names(gt_mut)[gt_f]
    present_samples = str_c(present_samples, collapse = ";")
    return(present_samples)
}

#Filter for variants in a specific tissue
filter_tissue = function(vcf, tissue){
    gt = get_gt(vcf)
    tissue_cols = str_detect(colnames(gt), tissue)
    gt = gt[,tissue_cols, drop = F]
    vcf = vcf[rowSums(gt)>=1,]
    return(vcf)
}

#Filter on several quality parameters
filter_quality = function(vcf){
    
    #Set filters
    gt = get_gt(vcf)
    gq_f = geno(vcf)$GQ >= 20
    back_alt_f = geno(vcf)$BW >= 4
    forw_alt_f = geno(vcf)$FW >= 4
    
    #Combine filters
    combi_f = gq_f & back_alt_f & forw_alt_f & gt
    combi_f = rowSums(combi_f) >= 1
    
    #Filter and return
    vcf = vcf[combi_f]
    return(vcf)
}

#Filter highvaf others
filter_highvaf_others = function(vcf){
    gt = get_gt(vcf)
    vaf = geno(vcf)$VF
    vaf_f = rowSums(gt != 1 & vaf > 0.1) >= 1
    vcf = vcf[!vaf_f]
    return(vcf)
}

#Remove variants near indels
remove_variants_near_indels = function(vcf, indel_vcf){
    indel_gr = GenomicRanges::reduce(granges(indel_vcf))
    hits = findOverlaps(vcf, indel_gr, maxgap = 10)
    
    if (length(hits)){
        vcf = vcf[-queryHits(hits)]
    }
    return(vcf)
}
