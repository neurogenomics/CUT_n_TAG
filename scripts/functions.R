
GITHUB.list_files <- function(creator="RajLabMSSM",
                              repo="Fine_Mapping_Shiny",
                              query=NULL,
                              return_download_api=T){
  repo_api <- file.path("https://api.github.com/repos",creator,repo,
                        "git/trees/master?recursive=1")
  req <- httr::GET(repo_api)
  httr::stop_for_status(req)
  filelist <- unlist(lapply(httr::content(req)$tree, "[", "path"), use.names = F)
  print(paste(length(filelist),"files found in GitHub repo:", file.path(creator,repo)))
  if(!is.null(query)){
    # query_string <- "*Nalls23andMe_2019.*UKB.multi_finemap.csv.gz"
    bool <- grepl(query, filelist)
    filelist <- filelist[bool]
    print(paste(length(filelist),"files found matching query."))
  }
  if(return_download_api){
    filelist <- file.path("https://github.com",creator,repo,"raw/master",filelist)
  }
  return(filelist)
}


GITHUB.download_files <- function(filelist,
                                  download_dir="./",
                                  overwrite=F,
                                  nThread=parallel::detectCores()){
  local_files <- parallel::mclapply(filelist, function(x){
    print(paste("Downloading",x))
    destfile <-  gsub("https://github.com/*.*/raw/master/www/data",
                      download_dir,x)
    dir.create(dirname(destfile), showWarnings = F, recursive = T)
    if(!file.exists(destfile) & overwrite==F) download.file(url = x, destfile=destfile)
    return(destfile)
  }, mc.cores = nThread) %>% unlist()
}


GITHUB.make_data_dict <- function(named_lists){
  data_dict <- list()
  for(x in names(named_list)){
    dat <- named_list[[x]]
    data_dict[[x]] <- as.list(setNames(dat, basename(dirname(dirname(dat))) ) )
  }
  # Add the locus dir as a bonus
  data_dict$locus_dir <- as.list(setNames(dirname(dirname(dat)), basename(dirname(dirname(dat))) ))
  # Dataset dir
  data_dict$dataset_dir <- dirname(dirname(dirname(dat)))[1]
  # Dataset name
  data_dict$dataset_type <- basename(dirname(dirname(dirname(dirname(dat)))))[1]
  # data
  data_dict$dataset <- basename(dirname(dirname(dirname(dat))))[1]
  return(data_dict)
}




GITHUB.find_pages <- function(creator="RajLabMSSM",
                              repo="Fine_Mapping",
                              local_repo=NULL,
                              return_table=T,
                              save_path=NULL,
                              dir_start=2){
  if(is.null(local_repo)){
    filelist <- GITHUB.list_files(creator=creator,
                                  repo=repo,
                                  query="*.*.html",
                                  return_download_api=F)
  } else {
    filelist <- gsub("^[.][/]","",list.files(path = local_repo, pattern="*.*\\.html", full.names=T, recursive=T))
  }
  gh_pages_url <- file.path(paste0("https://",creator,".github.io"),repo)
  gh_pages_links <- file.path(gh_pages_url, filelist)
  if(return_table){
    links_df <- data.frame(creator=creator,
                           repo=repo,
                           dir1=unlist(lapply(filelist, function(x){strsplit(x,"/")[[1]][dir_start]} )),
                           dir2=unlist(lapply(filelist, function(x){strsplit(x,"/")[[1]][dir_start+1]} )),
                           dir3=unlist(lapply(filelist, function(x){strsplit(x,"/")[[1]][dir_start+2]} )), 
                           dir4=unlist(lapply(filelist, function(x){strsplit(x,"/")[[1]][dir_start+3]} )), 
                           url=gh_pages_links,
                           link=paste0("<a href='",gh_pages_links,"' target='blank'>",filelist,"</a>"))
    if(!is.null(save_path)){
      print(paste("Writing links to ==>",save_path))
      data.table::fwrite(links_df,save_path, sep=",")
      write.table(paste(pages[["link"]],collapse="\n\n"), "results.md")
    }
    return(links_df)
  } else {
    return(gh_pages_links)
  }
}



message_parallel <- function(...){
  system(sprintf('echo "%s"', paste0(..., collapse="")))
}



# Extension of echolocatoR::import.bw.filt
import.bw.parallel <- function(bw.file,
                               gr.dat,
                               parallel_chrom=T,
                               nThread=parallel::detectCores()-1,
                               bw.file_format="UCSC"){
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.dat) <- bw.file_format)
  if(parallel_chrom){ 
    # It's muuuch faster to just iterate over each chromosome, 
    ## Get all ranges within min/max POS, 
    #and then filter to just the relevant ranges within your ranges.
    if(nThread>1) print(paste("Parallelizing queries across",nThread,"cores"))
    chroms <- unique(GenomicRanges::seqnames(gr.dat))
    gr.list <- parallel::mclapply(chroms, function(chr){
      message_parallel(paste("Querying chrom =",chr)) 
      gr.chr <-  gr.dat[GenomicRanges::seqnames(gr.dat)==chr,] 
      bw.chr <- rtracklayer::import.bw(con = bw.file,
                                       selection = gr.chr) 
      hits <- GenomicRanges::findOverlaps(query = gr.chr, subject = bw.chr)
      bw.chr <- bw.chr[S4Vectors::subjectHits(hits), ] 
      return(bw.chr)
    }, mc.cores = nThread) %>%
      GenomicRanges::GRangesList(compress = F) 
    gr.bind <- GenomicRanges::sort(unlist(gr.list))
    return(gr.bind)
  
  } else {
    # Otherwise, just use the score for the exact values 
    # Tends to be farrrr slower
    bw.filt <- rtracklayer::import.bw(con = bw.file,
                                      selection = gr.dat)
    # bw.dat <- rtracklayer::BigWigSelection(ranges = gr.dat,  colnames = "score")
    return(bw.filt)
  } 
}

 

prepare_tagMatrix <- function(peaks_list,
                              nThread=parallel::detectCores()-1){ 
  TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  tagMatrix_list <- parallel::mclapply(names(peaks_list), function(x){
    message_parallel(x) 
    gr <- peaks_list[[x]]
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC")
    promoter <- ChIPseeker::getPromoters(TxDb=TxDb,
                             upstream=3000, 
                             downstream=3000)
    tagMatrix <- ChIPseeker::getTagMatrix(peak = gr,
                              windows=promoter)
    return(tagMatrix)
  }, mc.cores = nThread) %>% `names<-`(names(peaks_list))
  return(tagMatrix_list)
}




prepare_annotatePeak <- function(peaks_list,
                                 nThread=parallel::detectCores()-1){
  TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene 
  parallel::mclapply(names(peaks_list), function(x){  
    message_parallel(x)
    gr <- peaks_list[[x]] 
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC") 
    peakAnno <- ChIPseeker::annotatePeak(peak = gr,
                                         tssRegion=c(-3000, 3000),
                                         TxDb=TxDb,
                                         annoDb="org.Hs.eg.db")
    return(peakAnno)
  }, mc.cores = nThread) %>% `names<-`(names(peaks_list))
}



peaks_plot <- function(CT.bw_filt,
                       geom="histogram",
                       show_plot=T){ 
  gg_hist <- ggbio::autoplot(CT.bw_filt, 
                             geom=geom, 
                             aes(fill=score),
                             adjust=.2,
                             binwidth=1000
  ) +
    theme_classic() 
  
  max_height <- echolocatoR::PLOT.get_max_histogram_height(gg=gg_hist)
  rect_height <- max_height / if(geom=="density") 8 else 10
  CT_narrowPeaks_filt$y <- 0 - rect_height
  
  gg <- gg_hist +
    ggbio::geom_rect(CT_narrowPeaks_filt,
                     stat="identity", 
                     rect.height= rect_height,
                     aes(y=y, fill=peak_score, label=Annotation),
                     alpha = .7, inherit.aes = F,
                     color="transparent",
                     hjust=1) +
    theme_classic() +
    ggbio::zoom_in(fac = 1/20) +
    ggbio::scale_x_sequnit(unit = "Mb")
  if(show_plot) print(gg)
  return(gg)
}



compare_peak_overlap <- function(gr.query, 
                                 gr.subject){
  # gr.query <- CT.narrowPeaks; gr.subject <- ENCODE.broadPeaks
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.query) <- "UCSC")
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr.subject) <- "UCSC")
  
  hits <- GenomicRanges::findOverlaps(query = gr.query, 
                                      subject = gr.subject)
  query_hits <- gr.query[S4Vectors::queryHits(hits)]
  print(
    paste0(length(query_hits)," / ",length(gr.query),
           " (",round(length(query_hits)/length(gr.query)*100,2),"%)",
           " of query peaks overlap with subject peaks."
    )
  )
  subject_hits <- gr.subject[S4Vectors::subjectHits(hits)]
  print(
    paste0(length(subject_hits)," / ",length(gr.subject),
           " (",round(length(subject_hits)/length(gr.subject)*100,2),"%)",
           " of subject peaks overlap with query peaks."
    )
  )
}





extract_gene_lists <- function(annotatePeak_list,
                               use_seq2gene=F){
  TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_lists <- lapply(names(annotatePeak_list), function(x,
                                                          .use_seq2gene=use_seq2gene){ 
    if(.use_seq2gene){
      genes <- ChIPseeker::seq2gene(annotatePeak_list[[x]]@anno, 
                                    tssRegion = c(-1000, 1000), 
                                    flankDistance = 3000,
                                    TxDb=TxDb) 
    } else{ genes <- annotatePeak_list[[x]]@anno$geneId } 
    return(genes)
  }) %>% `names<-`(names(annotatePeak_list))
  return(gene_lists)
}



enrich_clusterProfiler <- function(annotatePeak_list,
                                   fun="enrichKEGG",
                                   pvalueCutoff  = 0.05,
                                   pAdjustMethod = "BH",
                                   use_seq2gene=F,
                                   show_plot=T){ 
  gene_lists <- extract_gene_lists(annotatePeak_list, 
                                   use_seq2gene = use_seq2gene) 
  pathways <- clusterProfiler::compareCluster(geneCluster   = gene_lists,
                                              fun = fun,
                                              pvalueCutoff  = pvalueCutoff,
                                              pAdjustMethod = pAdjustMethod)
  .plot <- clusterProfiler::dotplot(pathways, 
                                    showCategory = 15, 
                                    title = paste(if(use_seq2gene)"seq2gene +",
                                                  gsub("^enrich","",fun),"pathway enrichment"))
  if(show_plot)print(.plot)
  return(list(pathways=pathways,
              dotplot=.plot ) )
}




enrich_ReactomePA <- function(annotatePeak_list,
                              pvalueCutoff  = 0.05,
                              pAdjustMethod = "BH",
                              use_seq2gene=F,
                              show_plot=T){ 
  # ReactomePA::enrichPathway can only take one gene list at once for you have to iterate
  res <- lapply(names(annotatePeak_list), function(x,
                                                   .use_seq2gene=use_seq2gene,
                                                   .pvalueCutoff=pvalueCutoff,
                                                   .pAdjustMethod=pAdjustMethod){ 
    print(x)  
    genes <- extract_gene_lists(annotatePeak_list[x], 
                                use_seq2gene = .use_seq2gene)[[1]] 
    pathways <- ReactomePA::enrichPathway(genes, 
                                          pvalueCutoff=.pvalueCutoff,
                                          pAdjustMethod=.pAdjustMethod)
    .plot <- ReactomePA::dotplot(pathways,
                                 title=paste(x,":\n",
                                             if(.use_seq2gene)"seq2gene +",
                                             "ReactomePA pathway enrichment"))
    if(show_plot) print(.plot) 
    return(list(pathways=pathways, 
                dotplot=.plot) )
  }) %>% `names<-`(names(annotatePeak_list))
  return(res)
}



gene_vennplot <- function(annotatePeak_list,
                          use_seq2gene=F,
                          show_plot=T){
  TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_lists <- extract_gene_lists(annotatePeak_list,
                                   use_seq2gene=use_seq2gene)
  # Plot
  vp <- ChIPseeker::vennplot(gene_lists)
  if(show_plot)print(vp)
  # Return
  return(list(gene_lists=gene_lists,
              vennplot=vp))
}


