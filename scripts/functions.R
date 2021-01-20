
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
                               nThread=parallel::detectCores()){
  if(parallel_chrom){ 
    # It's muuuch faster to just iterate over each chromosome, 
    ## Get all ranges within min/max POS, 
    #and then filter to just the relevant ranges within your ranges.
    if(nThread>1) print(paste("Parallelizing queries across",nThread,"cores"))
    chroms <- unique(GenomicRanges::seqnames(gr.dat))
    gr.list <- parallel::mclapply(chroms, function(chrom){
      message_parallel(paste("Querying chrom =",chrom))
      gr.chr <-  gr.dat[GenomicRanges::seqnames(gr.dat)==chrom,] 
      bw.chr <- rtracklayer::import.bw(con = bw.file,
                                       selection = gr.chr) 
      hits <- GenomicRanges::findOverlaps(query = gr.chr, subject = bw.chr)
      bw.chr <- bw.chr[S4Vectors::subjectHits(hits), ] 
      return(bw.chr)
    }, mc.cores = nThread) %>%
      GenomicRanges::GRangesList(compress = F)
    GenomicRanges::sort(gr.bind)
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



#### 

prepare_tagMatrix <- function(peaks_list,
                              nThread=parallel::detectCores()-1){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(ChIPseeker)
  TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
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
  parallel::mclapply(names(peaks_list), function(x){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    message_parallel(x)
    gr <- peaks_list[[x]] 
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC") 
    peakAnno <- ChIPseeker::annotatePeak(peak = gr,
                                         tssRegion=c(-3000, 3000),
                                         TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                                         annoDb="org.Hs.eg.db")
    return(peakAnno)
  }, mc.cores = nThread) %>% `names<-`(names(peaks_list))
}





