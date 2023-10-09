## -------------------------------------------------------------------------------
## checkReads()
##
## check indel, unsupposed/silent mutation, secondary base call
## --------------------------------------
##
## input: contigList
## # A tibble: 2 × 4
##   Path                                          File             Stem  Direction
##   <chr>                                         <chr>            <chr> <chr>    
## 1 ./data/JC00030617test//Dar0106pro_premix.ab1  Dar0106pro_prem… Dar0… forward  
## 2 ./data/JC00030617test//Dar0106term_premix.ab1 Dar0106term_pre… Dar0… reverse  
##
## input: reads
## # A tibble: 2 × 29
##   contigName readName          direction msg   checkedPrimaryS… checkedSecondar…
##   <chr>      <chr>             <chr>     <chr> <chr>            <chr>           
## 1 Dar0106    Dar0106pro_premi… forward   Inde… ---------------… ---------------…
## 2 Dar0106    Dar0106term_prem… reverse   Inde… ATGAAATACCTATTG… ATGAAATACCTATTG…
## # … with 23 more variables: trimmedPrimarySequence <list>,
##
## input: output_dir
## [1] "./result/JC00030617test/"
## --------------------------------------
i <- 1

checkReads <- function(contigList, reads, output_dir){

    log_info(sprintf('%s started.', 'checkReads()'))

    ## create output_dir
    outdir_reads <- sprintf('%s/reads', output_dir)
    dir.create(outdir_reads, showWarnings=FALSE)

    ## check each reads in same contigList. 
    for (i in (1: nrow(contigList))){
        
        p <- contigList$path[i]
        ctn <- contigList$contigName[i]

        ## all reads are treated as 'forward'
        read <- mySangerRead(p)

        ## read information, initialize, update direction using contig info
        gr <- getRead(read, ctn, read5primeOffset, read3primeOffset) %>%
            mutate(result='') %>%
            mutate(direction=contigList$direction[i]) %>%
            mutate(n_insertion=NA) %>%
            mutate(n_deletion=NA) %>%
            mutate(n_unmatch=NA) %>%
            mutate(n_unmatch_msg=NA) %>%
            mutate(n_unmatch_silent=NA) %>%
            mutate(n_unmatch_silent_msg=NA) %>%
            mutate(checkedPrimarySequeneLength=NA) %>%
            mutate(checkedPrimarySequence='-') %>%
            mutate(checkedSecondarySequence='-') %>%
            mutate(n_2nd_basecall=NA) %>%
            mutate(n_2nd_basecall_msg=NA) %>%
            mutate(n_2nd_basecall_N=NA) %>%
            mutate(n_2nd_basecall_N_msg=NA) 

        ## get trimmed sequence
        trimmed <- c(DNAStringSet(gr$trimmedPrimarySequence[[1]]),
                     DNAStringSet(gr$trimmedSecondarySequence[[1]]))
        
        ## align to template
        if (gr$direction == 'forward'){
            aln <- c(TEMPLATE, trimmed)
        } else if (gr$direction == 'reverse'){
            aln <- c(TEMPLATE, reverseComplement(trimmed))
        } else {
            log_error(sprintf('Read direction error: %s', gr$readName))
            stop('Quit.')
        }

        ## Align and logging
        names(aln) <- c('template', 'primary', 'secondary')
        aln <- AlignSeqs(aln, processors=7, verbose=FALSE)
        out_path <- sprintf('%s/%s_0_rawseq.fa', outdir_reads, contigList$file[i])
        writeXStringSet(aln, out_path)
        
        ## Crop and logging
        df <- cropSeq(dnastringset2df(aln))
        out_path <- sprintf('%s/%s_1_cropped.fa', outdir_reads, contigList$file[i])
        writeXStringSet(dataframe2ds(df), out_path)

        ## minRead check
        gr$checkedPrimarySequenceLength <- (nchar(df$Sequence[2]) - count_dash(df$Sequence[2]))
        if (gr$checkedPrimarySequenceLength < minReadLength){
            gr$result <- 'Read length is shorter than minReadLength'
            reads <- bind_rows(reads, gr)
            next
        }
        
        ## Insertion check
        n_insertion <- checkInsertion(df)
        log_info(sprintf('n_insertion: %s', n_insertion))
        if (n_insertion > 0) {
            ## store information
            msg <- 'Insertion detected.'
            gr$n_insertion <- n_insertion
            gr$result <- msg
            reads <- bind_rows(reads, gr)
            next
        } 
        
        ## Deletion check
        n_deletion <- checkDeletion(df)
        log_info(sprintf('n_deletion: %s', n_deletion))        
        if (n_deletion>0) {
            ## store information
            gr$n_deletion <- n_deletion
            msg = 'Deletion detected.'
            gr$result <- msg
            reads <- bind_rows(reads, gr)
            next
        }
        ## check umnatch
        gr <- checkUnmatchM(df, gr, outdir_reads)

        
        ## Contamination check
        l <- check2ndBaseCall(df)
        
        gr$result <- str_c(gr$result, ', ', l$msg)
        gr$checkedPrimarySequence <- df$Sequence[2]
        gr$checkedSecondarySequence <- df$Sequence[3]
        gr$n_2nd_basecall <- l$n_2nd_basecall
        gr$n_2nd_basecall_msg <- l$n_2nd_basecall_msg
        gr$n_2nd_basecall_N <- l$n_2nd_basecall_N
        gr$n_2nd_basecall_N_msg <- l$n_2nd_basecall_N_msg
        log_info(sprintf('n_2nd_basecall: %s', l$n_2nd_basecall))
        log_info(sprintf('n_2nd_basecall_N: %s', l$n_2nd_basecall_N))                 
        ## Concat
        reads <- bind_rows(reads, gr)
    }
    
    return(reads)
}


## -------------------------------------------------------------------------------

