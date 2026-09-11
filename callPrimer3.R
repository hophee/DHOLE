#' call primer3 for a given set of DNAstringSet object
#' Forked from https://gist.github.com/al2na/8540391
#' TODO: Add support for target amplicon region (Maybe as [] in the fasta input)
#' @param seq: DNA template as a character string (required)
#' @param fw_primer: optional forward (left) primer, if provided Primer3 will assess it and will try to find a suitable reverse primer. Default: NULL
#' @param rv_primer: optional reverse (right) primer (must be reverse-complemented to the template), if provided Primer3 will assess it and will try to find a suitable forward primer. Default: NULL
#' @param size_range: a string with space separated list of desired amplicon size ranges. Default: '151-500'
#' @param Tm: range of melting temprature parameters as a numerical vector containing (min,optimal,max). default: c(55,57,58)
#' @param name: name of the amplicon in 'chr_start_end' format
#' @param sequence_target: a string containing space separated list of target pairs: 'starting_position,target_length starting_position2,target_length2'. Default: NULL
#' @param SSR_exclude_region: Usually same as target, make sure the SSR is not considered as an internal oligo. Default: FALSE
#' @param report: path to output file for Primer3 report (use R valid path). Default: NULL
#' @param seq_oh_left: a string of nucleotides to be added to the 5' end of the left primer (requires Primer3 >2.5). Default: NULL
#' @param seq_oh_right: a string of nucleotides to be added to the 5' end of the right primer (requires Primer3 >2.5). Default: NULL
#' @param primer_num: number of primers to return (an integer). Default: 4  
#' @param p3_input_file: path to primer3 input file that will be created (use R valid path). Default: temporary file
#' @param p3_output_file: path to primer3 output file that will be created (use R valid path). Default: temporary file
#' @param primer3: path to primer3 location (use R valid path)
#' @param thermo.param: path to thermodynamic parameters folder (use R valid path)
#' @param settings: path to text file with saved parameters (use R valid path)
#' @author Ido Bar forked from Altuna Akalin which modified Arnaud Krebs' original function
#' @example
#'
#' primers=mapply(callPrimer3,seq=sample.seq,
#'                                   name=names(sample.seq),
#'                      MoreArgs=list(primer3=primer3,size_range=size_range,Tm=Tm,
#'                                    thermo.param=thermoParam,settings=settings)
#'                      ,SIMPLIFY = FALSE)
#'
#'
callPrimer3 <- function(seq,size_range='151-500',Tm=c(55,57,58), Tm_diff=5, name,sequence_target=NULL, fw_primer=NULL, rv_primer=NULL,
                        SSR_exclude_region=FALSE,report=NULL,seq_oh_left=NULL, seq_oh_right=NULL, primer_num=4,
                        p3_input_file=tempfile(),
                        p3_output_file=tempfile(),
                        primer3=path.expand("~/OneDrive - Griffith University/Research/Methods/Bioinformatics/primer3-2.6.1/src/primer3_core.exe"),
                        thermo.param=path.expand("~/OneDrive - Griffith University/Research/Methods/Bioinformatics/primer3-2.6.1/src/primer3_config/"),
                        settings=path.expand("~/OneDrive - Griffith University/Research/Methods/Bioinformatics/primer3-2.6.1/settings_files/primer3_melt_mama_settings.txt")){
  
  stderr_file <- if (is.null(report)) tempfile() else paste0(report, ".stderr.log")
  on.exit(unlink(c(p3_input_file, p3_output_file,
                  if (is.null(report)) stderr_file)), add = TRUE)
  primer3_error <- function(message) {
    stop(structure(list(message = message, call = NULL, stage = "primer3"),
                   class = c("primer3_error", "error", "condition")))
  }
  run_primer3 <- function(args) {
    output <- tryCatch(
      suppressWarnings(system2(primer3, shQuote(args), stdout = TRUE,
                               stderr = stderr_file)),
      error = function(e) primer3_error(conditionMessage(e))
    )
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
      detail <- if (file.exists(stderr_file)) paste(readLines(stderr_file, warn = FALSE), collapse = "\n") else ""
      primer3_error(sprintf("Primer3 exit status %d: %s", status, detail))
    }
    output
  }
  tryCatch({
  seq_target_str <- if (is.null(sequence_target)) NULL else sprintf("SEQUENCE_TARGET=%s",sequence_target)
  exclude_region_str <- if (SSR_exclude_region & !is.null(sequence_target)) sprintf("SEQUENCE_INTERNAL_EXCLUDED_REGION=%s\n", sequence_target) else NULL
  # identify primer3 version used
  pr3_ver <-  as.numeric(sub(".+release (\\d+\\.\\d+)\\.*(\\d*)$", "\\1\\2", 
                             run_primer3("--about")))
  if (length(pr3_ver) != 1L || is.na(pr3_ver)) primer3_error("Некорректный ответ Primer3 --about")
  sequence_overhung_left <- if (pr3_ver>2.5 & !is.null(seq_oh_left)) sprintf("SEQUENCE_OVERHANG_LEFT=%s",seq_oh_left) else  NULL
  sequence_overhung_right <- if (pr3_ver>2.5 & !is.null(seq_oh_right)) sprintf("SEQUENCE_OVERHANG_RIGHT=%s",seq_oh_right) else NULL
  fw_primer_str <- if (is.null(fw_primer)) NULL else sprintf("SEQUENCE_PRIMER=%s",fw_primer)
  rv_primer_str <- if (is.null(rv_primer)) NULL else sprintf("SEQUENCE_PRIMER_REVCOMP=%s",rv_primer)
  
  #print(excluded.regions)
  # make primer 3 input file
  cat(sprintf("SEQUENCE_ID=%s",name  ),
            sprintf("SEQUENCE_TEMPLATE=%s",as.character(seq)),
            seq_target_str,  
            fw_primer_str, 
            rv_primer_str,
            sequence_overhung_left,
            sequence_overhung_right,
            exclude_region_str,
            "PRIMER_TASK=pick_detection_primers",
            "PRIMER_PICK_LEFT_PRIMER=1",
            "PRIMER_PICK_RIGHT_PRIMER=1",
            "PRIMER_PICK_INTERNAL_OLIGO=0",
            "PRIMER_EXPLAIN_FLAG=1",
            sprintf("PRIMER_NUM_RETURN=%s",primer_num),
            sprintf("PRIMER_PAIR_MAX_DIFF_TM=%s" ,Tm_diff),
            sprintf("PRIMER_MIN_TM=%s" ,Tm[1]),
            sprintf("PRIMER_OPT_TM=%s" ,Tm[2]),
            sprintf("PRIMER_MAX_TM=%s" ,Tm[3]),
            sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s" ,size_range),
            sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s" ,thermo.param),
            "=", sep="\n", file = p3_input_file
  )
  run_primer3(c(paste0("--p3_settings_file=", settings),
                paste0("--output=", p3_output_file), p3_input_file))
  if (!file.exists(p3_output_file) || file.size(p3_output_file) == 0L) {
    primer3_error("Primer3 output is missing or empty")
  }
  lines <- readLines(p3_output_file, warn = FALSE)
  errors <- sub("^PRIMER_ERROR=", "", grep("^PRIMER_ERROR=.+", lines, value = TRUE))
  if (length(errors)) primer3_error(paste(errors, collapse = "; "))
  out <- read.delim(p3_output_file, sep = "=", header = FALSE, fill = TRUE,
                    quote = "", comment.char = "", colClasses = "character")
  returned.primers <- as.integer(out[out[, 1] == "PRIMER_PAIR_NUM_RETURNED", 2])
  if (length(returned.primers) != 1L || is.na(returned.primers) || returned.primers < 0L) {
    primer3_error("Primer3 output lacks a valid PRIMER_PAIR_NUM_RETURNED")
  }
  if (returned.primers == 0L) return(data.frame())
  if (returned.primers>0){
    designed.primers=data.frame()
    for (i in seq(0,returned.primers-1,1)){
      #IMPORT SEQUENCES
      id=sprintf(  'PRIMER_LEFT_%i_SEQUENCE',i)
      PRIMER_LEFT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      id=sprintf(  'PRIMER_RIGHT_%i_SEQUENCE',i)
      PRIMER_RIGHT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      
      #IMPORT PRIMING POSITIONS
      id=sprintf(  'PRIMER_LEFT_%i',i)
      PRIMER_LEFT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #PRIMER_LEFT_LEN=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      id=sprintf(  'PRIMER_RIGHT_%i',i)
      PRIMER_RIGHT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #IMPORT Tm
      id=sprintf(  'PRIMER_LEFT_%i_TM',i)
      PRIMER_LEFT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      id=sprintf(  'PRIMER_RIGHT_%i_TM',i)
      PRIMER_RIGHT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      
      res=out[grep(i,out[,1]),]
      extra.inf=t(res)[2,,drop=FALSE]
      colnames(extra.inf)=sub( paste("_",i,sep=""),"",res[,1])
      extra.inf=extra.inf[,-c(4:9),drop=FALSE] # remove redundant columns
      extra.inf=apply(extra.inf,2,as.numeric)
      #Aggegate in a dataframe
      primer.info=data.frame(PRIMER_ID=name, PRIMER_NUM=i+1, 
                             PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE,
                             PRIMER_LEFT_TM, PRIMER_RIGHT_TM,
                             PRIMER_LEFT_pos=PRIMER_LEFT[1],
                             PRIMER_LEFT_len=PRIMER_LEFT[2],
                             PRIMER_RIGHT_pos=PRIMER_RIGHT[1],
                             PRIMER_RIGHT_len=PRIMER_RIGHT[2],
                             t(data.frame(extra.inf))
                             
      )
      rownames(primer.info)=NULL
      designed.primers=rbind(designed.primers, primer.info)
      #print(primer.info)
    }
    
    
    
    #colnames(designed.primers)=c('PrimerID',
    #                             'Fwseq','Rvseq',
    #                             'FwTm','RvTm',
    #                             'FwPos','Fwlen',
    #                             'RvPos','Rvlen',
    #                             'fragLen' )
    if (!is.null(report)) {
      run_primer3(c(paste0("--p3_settings_file=", settings), "--format_output",
                    paste0("--output=", report), p3_input_file))
    }
  }
  return(designed.primers)
  }, error = function(e) {
    if (inherits(e, "primer3_error")) stop(e)
    primer3_error(conditionMessage(e))
  })
}


#' convert primer list to genomic intervals
#'
#' function returns a GRanges or data frame object from a list of primers designed
#' by \code{designPrimers} function after calculating genomic location of the amplicon
#' targeted by the primers.
#'
#' @param primers a list of primers returned by \code{designPrimers} function
#' @param as.data.frame logical indicating if a data frame should be returned
#'        instead of \code{GRanges} object. Default=FALSE
#'
#' @examples
#'  data(bisPrimers)
#'  # remove data with no primers found
#'  bisPrimers=bisPrimers[!is.na(bisPrimers)]
#'  gr.pr=primers2ranges(bisPrimers) # convert primers to GRanges
#'
#' @seealso \code{\link{filterPrimers}}, \code{\link{designPrimers}}
#'
#' @export
#' @docType methods
primers2ranges<-function(primers,as.data.frame=FALSE){
  if (!require(GenomicRanges)) stop("GenomicRanges not installed, use BiocManager::install('GenomicRanges') and try again.")
  
  if(any(is.na(primers))){
    warning( "There are targets without primers\nfiltering those before conversion")
    primers=primers[ !is.na(primers) ]
  }
  df=do.call("rbind",primers) # get primers to a df
  locs=gsub("\\|\\.+","",rownames(df)) # get the coordinates from list ids
  temp=do.call("rbind",strsplit(locs,"_")) #
  
  start=as.numeric(temp[,2])
  chr=as.character(temp[,1])
  
  
  
  amp.start= start + as.numeric(df$PRIMER_LEFT_pos)
  amp.end  = start + as.numeric(df$PRIMER_RIGHT_pos)
  res=data.frame(chr=chr,start=amp.start,end=amp.end,df)
  #saveRDS(res,file="/work2/gschub/altuna/projects/DMR_alignments/all.designed.primers.to.amps.rds")
  if(as.data.frame) return(res)
  gr=GRanges(seqnames=res[,1],ranges=IRanges(res[,2],res[,3]) )
  values(gr)=res[,-c(1,2,3)]
  return(gr)
}

