## -------------------------------------------------------------------------------
## contig check functions:
## - count valid reads
## - degenerated symbols
## - (then apply reads function) unsupposed mutation
## - (then apply reads function) silent mutation
## - recap 2nd basecall
##
## input:
##
## output:
##
## -------------------------------------------------------------------------------
checkContig <- function(contig,  output_dir){

    log_info(sprintf('%s started.', 'checkContig()'))

    ## get contigName    
    ctn <- unique(contig$contigName)
    
    ## dataframe for return
    ## For each contig
    tibble(contigName=ctn) %>%
        mutate(result='') %>%
        mutate(n_degenerated_symbol=NA) %>%
        mutate(n_degenerated_symbol_msg=NA) %>%
        mutate(n_unmatch=NA) %>%
        mutate(n_unmatch_msg=NA) %>%
        mutate(n_unmatch_silent=NA) %>%
        mutate(n_unmatch_silent_msg=NA) %>%    
        mutate(n_2nd_basecall=NA) %>%
        mutate(n_2nd_basecall_msg=NA) %>%
        mutate(n_2nd_basecall_N=NA) %>%
        mutate(n_2nd_basecall_N_msg=NA) %>%
        mutate(SequenceLength=NA) %>%
        mutate(Sequence='-') -> re
    
    if (CHECK_N_MUTATION==TRUE){
        re %>%
            mutate(n_unmatch_N=NA) %>%
            mutate(n_unmatch_N_msg=NA) %>%
            mutate(n_unmatch_N_silent=NA) %>%
            mutate(n_unmatch_N_silent_msg=NA) -> re
    }
    
    ## create output_dir
    outdir_contigs <- sprintf('%s/contigs', output_dir)
    dir.create(outdir_contigs, showWarnings=FALSE)

    ## check insertion/deletion
    contig %>%
        filter(is.na(n_insertion)) %>%
        filter(is.na(n_deletion)) -> contig.valid

    if (nrow(contig.valid)==0){
        re$result <- 'Invalid reads.'
        return(re)
    }
    ## assemble reads
    contig.valid %>%
        select(readName, Sequence) %>%
        rename(Name=readName) %>%
        dataframe2ds() -> aln

    ## get consensus
    con <- ConsensusSequence(aln)
    names(con) <- ctn
    consensussequence <- as.character(con)
    re$Sequence <- consensussequence

    ## check Length
    re$SequenceLength <- nchar(consensussequence) - count_dash(consensussequence)
    ## re$SequenceLength <- nchar(df.c$Sequence[2]) - count_dash(df.c$Sequence[2])
    log_info(sprintf('SequenceLength: %s', re$SequenceLength))
    
    if (re$SequenceLength < nchar(as.character(TEMPLATE))){
        re$result <- 'Contig length is shorter than template.'

        ## if length is shorter than threshold, skip sequence check
        if (re$SequenceLength < minContigLength){
            return(re)
        }
    } else {
        re$result <- 'OK'        
    }
    
    ## check degenerated symbols
    l <- checkDegeneratedSymbol(consensussequence)
    re$result <- str_c(re$result, ', ', l$msg)
    re$n_degenerated_symbol <- l$n_degenerated_symbol
    re$n_degenerated_symbol_msg <- l$n_degenerated_symbol_msg

    ## re-assesment(consensus) unsupposed mutations
    con.aln <- c(TEMPLATE, con)
    names(con.aln) <- c('template', 'primary')
    df.c <- dnastringset2df(con.aln)
    writeXStringSet(con.aln, sprintf('%s/%s_contig.fa', outdir_contigs, ctn))
    
    ## check unmatch
    re <- checkUnmatchMContig(df.c, re, outdir_contigs)

    ## check 2nd basecall
    re <- recap2ndBaseCall(contig.valid, re)

    return(re)
}



recap2ndBaseCall <- function(contig.valid, re){
    
    log_info(sprintf('%s started.', 'recap2ndBaseCall()'))
    
    ## re-assesment 2nd_basecall
    contig.valid %>%
        filter(!is.na(n_2nd_basecall_msg)) -> contig_2b
    
    contig.valid %>%
        filter(!is.na(n_2nd_basecall_N_msg)) -> contig_2b_N
    
    n_2nd_basecall <-nrow(contig_2b)     
    n_2nd_basecall_N <-nrow(contig_2b_N)

    log_info(sprintf('n_2nd_basecall: %s', n_2nd_basecall))
    log_info(sprintf('n_2nd_basecall_N: %s', n_2nd_basecall_N))
    
    if (n_2nd_basecall>0){
        n_2nd_basecall_msg <- unique(str_split(str_c(contig_2b$n_2nd_basecall_msg, collapse=','), pattern=',')[[1]])
    }
    if (n_2nd_basecall_N>0){
        n_2nd_basecall_N_msg <- unique(str_split(str_c(contig_2b_N$n_2nd_basecall_N_msg, collapse=','), pattern=',')[[1]])
    }

    ## Send message (where different peak detected)
    if (n_2nd_basecall_N > 0 & n_2nd_basecall > 0){
        n_2nd_basecall_N<- length(n_2nd_basecall_N_msg)
        n_2nd_basecall <- length(n_2nd_basecall_msg)        
        
        if (allow2ndBasecallNotN == TRUE){
            msg <-  'Secondary basecall detected(N).'
        } else {
            msg <- 'Secondary basecall detected(N/not N).'
        }
    } else if (n_2nd_basecall_N ==0  & n_2nd_basecall > 0){
        n_2nd_basecall <- length(n_2nd_basecall_msg)
        n_2nd_basecall_N <- NA
        n_2nd_basecall_N_msg <- NA
        if (allow2ndBasecallNotN == TRUE){
            msg <- 'OK'
        } else {
            msg <- 'Secondary basecall detected(not N).'
        }
    } else if (n_2nd_basecall_N > 0 & n_2nd_basecall == 0){
        n_2nd_basecall_N <- length(n_2nd_basecall_N_msg)        
        n_2nd_basecall <- NA
        n_2nd_basecall_msg <- NA
        msg <- 'Secondary basecall detected(N).'
    } else {
        n_2nd_basecall <- NA
        n_2nd_basecall_N <- NA
        n_2nd_basecall_msg <- NA
        n_2nd_basecall_N_msg <- NA        
        msg <- 'OK'
    }

    log_info(sprintf('n_2nd_basecall split: %s', n_2nd_basecall))
    log_info(sprintf('n_2nd_basecall_N split: %s', n_2nd_basecall_N))
        
    n_2nd_basecall_msg <- str_c(n_2nd_basecall_msg, collapse=',')
    n_2nd_basecall_N_msg <- str_c(n_2nd_basecall_N_msg, collapse=',')
    
    re$result <- str_c(re$result, ', ', msg)
    re$n_2nd_basecall     <- n_2nd_basecall
    re$n_2nd_basecall_msg <- n_2nd_basecall_msg
    re$n_2nd_basecall_N   <- n_2nd_basecall_N
    re$n_2nd_basecall_N_msg <- n_2nd_basecall_N_msg
    
    return(re)
}
