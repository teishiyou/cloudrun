## -------------------------------------------------------------------------------
## reads check functions:
## - insertion
## - deletion
## - secondary base call
## - mutation
##   - unmatch
##   - silent
## - mutation_N
##   - unamtch
##   - silent
##
## -------------------------------------------------------------------------------
## df[1]: templte
## df[2]: query: primary sequence
## df[3]: query: secondary sequence
## 
##   Name      Sequence                                                            
##   <chr>     <chr>                                                               
## 1 template  ATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCG…
## 2 primary   -------------------------------------------TCGCGGCCCAGCCGGCCATGGCCG…
## 3 secondary -------------------------------------------TCGCGGCCCAGCCGGCCATGGCCG…

## ---------------------------
## cropSeq()
##
## crop sequence dataframe as same length as template
## ---------------------------
cropSeq <- function(df){

    ## get sequence ends
    l <- getSequenceEnds(df$Sequence[1])
    
    ## crop start from end
    aln <- subseq(dataframe2ds(df), start=l$start, end=l$end)

    ## convert to dataframe
    df <- dnastringset2df(aln)

    return(df)
}

## ---------------------------
## checkInsertion()
##
## check sequence contains insertion or not.
## return number of inserted nt
## ---------------------------

checkInsertion <- function(df){
    log_info(sprintf('%s started.', 'checkInsertion()'))
    insert_pos <- str_locate_all(df$Sequence[1], pattern='[A|T|G|C|N]-+[A|T|G|C|N]')[[1]]

    if (length(insert_pos) == 0 ){
        return(0) # no insertion
    } else {
        n_insertion <- count_inside_dash(df$Sequence[1])
        return(n_insertion)
    }
}

## ---------------------------
## checkDeletion()
##
## check sequence contains deletion or not
## return number of deleted nt
## ---------------------------
checkDeletion <- function(df){
    log_info(sprintf('%s started.', 'checkDeletion()'))
    delete_pos <- str_locate_all(df$Sequence[2], pattern='[A|T|G|C|N]-+[A|T|G|C|N]')[[1]]

    if (length(delete_pos) == 0 ){
        return(0) # no deletion
    } else {
        n_deletion <- count_inside_dash(df$Sequence[2])
        return(n_deletion) 
    }
}

## ---------------------------
## transposeSeqDf()
##
## transpose df
## ---------------------------
transposeSeqDf <- function(df){
    ## transpose
    ## TODO: should be more easier way
    nseq <- nchar(df$Sequence[1])    
    
    df %>%
        separate(col=Sequence, sep='', into=(str_c('s', sprintf('%s', seq(0, nseq))))) %>%
        select(-c(s0)) -> df.tmp

    df.tmp %>%
        gather(key=Pos, value=Sequence, 2:ncol(df.tmp)) %>%
        mutate(Pos=as.numeric(str_replace(Pos, 's', ''))) %>%
        spread(key=names(df.tmp)[1], value=Sequence) -> df.t

    return(df.t)
}

## ---------------------------
## createDf()
##
## createDf for analysis
## ---------------------------
createDf <- function(df){
    
    ## --------------------------------- nucleotide
    df %>%
        transposeSeqDf() %>%
        rename(Template=template) %>%
        rename(Primary=primary) %>%
        rename(Secondary=secondary) %>%        
        mutate(Posa=map(Pos, posn2posa)) %>%
        unnest(cols=c(Posa)) %>%
        mutate(PosString=str_c(Template, Pos, Primary)) %>%        
        mutate(CompareNt=case_when(
                   Template==Primary ~ 'identical',
                   TRUE ~ 'different'
               )) -> df.nt

    ## --------------------------------- aa
    tibble(Posa=seq(1, nchar(df$Sequence[1])/3),
           TemplateTriplet=sequencent2triplet(df$Sequence[1]),
           PrimaryTriplet=sequencent2triplet(df$Sequence[2])) %>%
        mutate(TemplateAA=map(TemplateTriplet, translateA.seq)) %>%
        mutate(PrimaryAA=map(PrimaryTriplet, translateA.seq)) %>%
        unnest(cols=c(TemplateAA, PrimaryAA)) %>%
        mutate(PosStringAA=str_c(TemplateAA, Posa, PrimaryAA)) -> df.aa
    
    df.nt %>%
        left_join(df.aa, by='Posa') %>%
        mutate(PosStringFat=str_c(PosString, '(', PosStringAA, ')')) -> df.nt
    
    return(df.nt)
}

## ---------------------------
## createDf()
##
## createDf contig for analysis
## ---------------------------

createDfContig <- function(df){
    log_info(sprintf('%s started.', 'createDfContig()'))
    
    ## --------------------------------- nucleotide
    df %>%
        transposeSeqDf() %>%
        rename(Template=template) %>%
        rename(Primary=primary) %>%
        mutate(Posa=map(Pos, posn2posa)) %>%
        unnest(cols=c(Posa)) %>%
        mutate(PosString=str_c(Template, Pos, Primary)) %>%        
        mutate(CompareNt=case_when(
                   Template==Primary ~ 'identical',
                   TRUE ~ 'different'
               )) -> df.nt

    ## --------------------------------- aa
    tibble(Posa=seq(1, nchar(df$Sequence[1])/3),
           TemplateTriplet=sequencent2triplet(df$Sequence[1]),
           PrimaryTriplet=sequencent2triplet(df$Sequence[2])) %>%
        mutate(TemplateAA=map(TemplateTriplet, translateA.seq)) %>%
        mutate(PrimaryAA=map(PrimaryTriplet, translateA.seq)) %>%
        unnest(cols=c(TemplateAA, PrimaryAA)) %>%
        mutate(PosStringAA=str_c(TemplateAA, Posa, PrimaryAA)) -> df.aa
    
    df.nt %>%
        left_join(df.aa, by='Posa') %>%
        mutate(PosStringFat=str_c(PosString, '(', PosStringAA, ')')) -> df.nt
    
    return(df.nt)
}

## ---------------------------
## checkUnmatchMContig()
##
## check following:
##   unmatch
##   unmatch_silent
##   if codon specified: 
##       unmatch_N
##       unmatch_N_silent
## ---------------------------

checkUnmatchMContig <- function(df, gr, outdir_contigs){
    log_info(sprintf('%s started.', 'checkUnmatchMContig()'))

    dfa <- createDfContig(df)

    write_csv(dfa, sprintf('%s/%s.csv', outdir_contigs, gr$contigName))
    
    ## check except NNN 
    l <- checkUnmatchMutation(dfa)
    
    gr$result <- str_c(gr$result, ', ', l$msg)
    gr$n_unmatch <- l$n_unmatch
    gr$n_unmatch_msg <- l$n_unmatch_msg
    gr$n_unmatch_silent <- l$n_silent
    gr$n_unmatch_silent_msg <- l$n_silent_msg
    log_info(sprintf('n_unmatch: %s', l$n_unmatch))
    log_info(sprintf('n_unmatch_silent: %s', l$n_silent))
    
    ## check NNN
    if (CHECK_N_MUTATION==TRUE){
        l <- checkUnmatchMutation_N(dfa)
        gr$result <- str_c(gr$result, ', ', l$msg)
        gr$n_unmatch_N <- l$n_unmatch
        gr$n_unmatch_N_msg <- l$n_unmatch_msg
        gr$n_unmatch_N_silent <- l$n_silent
        gr$n_unmatch_N_silent_msg <- l$n_silent_msg
        log_info(sprintf('n_unmatch_N: %s', l$n_unmatch))
        log_info(sprintf('n_unmatch_N_silent: %s', l$n_silent))
    }

    return(gr)
}


## ---------------------------
## checkUnmatchM()
##
## check following:
##   unmatch
##   unmatch_silent
##   if codon specified: 
##       unmatch_N
##       unmatch_N_silent
## ---------------------------


checkUnmatchM <- function(df, gr, outdir_reads){
     log_info(sprintf('%s started.', 'checkUnmatchM()'))

     dfa <- createDf(df)

     write_csv(dfa, sprintf('%s/%s.csv', outdir_reads, gr$readName))
     
     ## check except NNN 
     l <- checkUnmatchMutation(dfa)

     gr$result <- l$msg
     gr$n_unmatch <- l$n_unmatch
     gr$n_unmatch_msg <- l$n_unmatch_msg
     gr$n_unmatch_silent <- l$n_silent
     gr$n_unmatch_silent_msg <- l$n_silent_msg
     log_info(sprintf('n_unmatch: %s', l$n_unmatch))
     log_info(sprintf('n_unmatch_silent: %s', l$n_silent))
     
     ## check NNN
     if (CHECK_N_MUTATION==TRUE){
         l <- checkUnmatchMutation_N(dfa)
         
         gr$result <- str_c(gr$result, ', ', l$msg)
         gr$n_unmatch_N <- l$n_unmatch
         gr$n_unmatch_N_msg <- l$n_unmatch_msg
         gr$n_unmatch_N_silent <- l$n_silent
         gr$n_unmatch_N_silent_msg <- l$n_silent_msg
         log_info(sprintf('n_unmatch_N: %s', l$n_unmatch))
         log_info(sprintf('n_unmatch_N_silent: %s', l$n_silent))
     }
     
     return(gr)
}

## ---------------------------
## checkUnmatchMutation()
##
## check whether the mutaion is unsupposed or not, then silent or notificationItem
##
## N='N'     : check inside  N
## N='not_N' : check outside N
##
## dfa:
## # A tibble: 200 × 12
##       Pos Primary Secondary Template  Posa PosString CompareNt TemplateTriplet PrimaryTriplet TemplateAA PrimaryAA PosStringAA
##     <dbl> <chr>   <chr>     <chr>    <dbl> <chr>     <chr>     <chr>           <chr>          <chr>      <chr>     <chr>      
##   1     1 -       -         A            1 A1-       different ATG             ---            M          -         M1-        
##   2     2 -       -         T            1 T2-       different ATG             ---            M          -         M1-        
##   3     3 -       -         G            1 G3-       different ATG             ---            M          -         M1-        
##   .     . .       .         .            .    .              .   .               .            .          .           .
## 135   135 A       A         A           45 A135A     identical GCA             GCA            A          A         A45A       
## 136   136 G       G         N           46 N136G     different NNN             GCG            -          A         -46A       
## 137   137 C       C         N           46 N137C     different NNN             GCG            -          A         -46A       
## 138   138 G       G         N           46 N138G     different NNN             GCG            -          A         -46A       
## 139   139 A       A         A           47 A139A     identical AGC             AGC            S          S         S47S       
## 140   140 G       G         G           47 G140G     identical AGC             AGC            S          S         S47S       
## 141   141 C       C         C           47 C141C     identical AGC             AGC            S          S         S47S       
## 142   142 G       G         G           48 G142G     identical GGT             GGT            G          G         G48G       
## 143   143 G       G         G           48 G143G     identical GGT             GGT            G          G         G48G       
## 144   144 T       T         T           48 T144T     identical GGT             GGT            G          G         G48G       
## 145   145 T       T         T           49 T145T     identical TTT             TTT            F          F         F49F       
## 146   146 T       T         T           49 T146T     identical TTT             TTT            F          F         F49F       
## 147   147 T       T         T           49 T147T     identical TTT             TTT            F          F         F49F       
##   .     . .       .         .            .    .              .   .               .            .          .           .
## ---------------------------

checkUnmatchMutation_N <- function(dfa){

    log_info(sprintf('%s started, N=%s', 'checkUnmatchMutation_N()', 'N'))

    ## Posa to exclude from check (PrimaryTriplet='---)
    dfa %>%
        ## select(-c(Pos, Primary, Secondary, Template, PosString, PosStringFat)) %>%
        select(Posa, CompareNt, TemplateTriplet, PrimaryTriplet, TemplateAA, PrimaryAA, PosStringAA) %>%
        filter(str_detect(PrimaryTriplet, '-')) %>%
        filter(str_detect(TemplateTriplet, 'N')) %>%
        distinct() %>%
        .$Posa -> exPos
    ## exPos: (ex)128, # of residue of AA

    
    ## all NNN are contained in PrimarySequence
    if (length(exPos)==0){
        exN <- 0
    ## some NNN are not contained in Primary Sequence
    } else {
        DESIGN %>%
            select(N, Posa) %>%
            distinct() %>%
            filter(Posa %in% exPos) %>%
            .$N -> exN
        ## exN: (ex)11, pos of mutated AA (1 ~ SEQLEN).
        ## when Mseq is 'QSAFAS', N is c(1, 2, 3, 4, 5, 6)
    }
    
    ## posN to compare
    posNComp <- setdiff(seq(1, max(DESIGN$N)), exN)
    
    ## filter 'NNN' Sequence
    dfa %>%
        ## select(-c(Pos, Primary, Secondary, Template, PosString, PosStringFat)) %>%
        select(Posa, CompareNt, TemplateTriplet, PrimaryTriplet, TemplateAA, PrimaryAA, PosStringAA) %>%
        filter(!str_detect(PrimaryTriplet, '-')) %>%
        filter(TemplateTriplet=='NNN') %>%
        distinct() %>%
        left_join(DESIGN, by='Posa', multiple='all') -> dfa.m

    ## compare PrimaryTriplet and Triplet(Designed)
    dfa.m %>%
        select(Posa, PrimaryTriplet, PrimaryAA, PosStringAA, N, Design, Triplet, AA) %>%
        mutate(CompareTriplet=case_when(
                   PrimaryTriplet==Triplet ~ 'identical',
                   TRUE ~ 'different')) %>%
        filter(CompareTriplet=='identical') -> a

    if (nrow(a)==length(unique(dfa.m$Posa))){
        n_unmatch <- 0
    } else {
        ## unmatch position N(1~length(MUTATION AA))
        unmatch.posN <- setdiff( posNComp, unique(a$N))
        n_unmatch <- length(unmatch.posN)
        log_info(sprintf('unmatch.posN: %s', str_c(unmatch.posN, collapse=',')))
    }
    
    if (n_unmatch==0){
        ret <- list(msg='OK',
                    n_unmatch=NA,
                    n_unmatch_msg=NA,
                    n_silent=NA,
                    n_silent_msg=NA)
        return(ret)
    }
    
    ## unmatch detected, then check codon are designed or not
    dfa %>%
        filter(Primary!='-') %>%
        filter(CompareNt=='different') %>%
        filter(Template=='N') %>%
        left_join(DESIGN, by='Posa', multiple='all') %>%
        filter(N %in% unmatch.posN) -> aa

## aa
## # A tibble: 6 × 18
##     Pos Primary Secondary Template  Posa PosString CompareNt TemplateTriplet PrimaryTriplet TemplateAA PrimaryAA PosStringAA PosStringFat     N Design   Triplet AA    Qc   
##   <dbl> <chr>   <chr>     <chr>    <dbl> <chr>     <chr>     <chr>           <chr>          <chr>      <chr>     <chr>       <chr>        <dbl> <chr>    <chr>   <chr> <chr>
## 1   484 C       C         N          162 N484C     different NNN             CAC            -          H         -162H       N484C(-162H)     7 WT       GAA     E     OK   
## 2   484 C       C         N          162 N484C     different NNN             CAC            -          H         -162H       N484C(-162H)     7 Design_1 CAG     Q     OK   
## 3   485 A       A         N          162 N485A     different NNN             CAC            -          H         -162H       N485A(-162H)     7 WT       GAA     E     OK   
## 4   485 A       A         N          162 N485A     different NNN             CAC            -          H         -162H       N485A(-162H)     7 Design_1 CAG     Q     OK   
## 5   486 C       C         N          162 N486C     different NNN             CAC            -          H         -162H       N486C(-162H)     7 WT       GAA     E     OK   
## 6   486 C       C         N          162 N486C     different NNN             CAC            -          H         -162H       N486C(-162H)     7 Design_1 CAG     Q     OK   

    ## create n_unmatch, n_silent msg
    n_unmatch_msg <- ''
    n_silent_msg  <- ''
    n_silent <- 0
    j <- 0

    for (posN in unmatch.posN) {
        aa %>%
            filter(N==posN) -> tmpaa
        
        designed_aa <- unique(tmpaa$AA)
        observed_aa <- unique(tmpaa$PrimaryAA)
        observed_triplet <- unique(tmpaa$PrimaryTriplet)

        ## unmatch
        n_unmatch_tmp <- str_c(posN, observed_aa, '(', observed_triplet, ')')
        
        if (j == 0){
            n_unmatch_msg <- str_c(n_unmatch_tmp)
        } else {
            log_info(sprintf('j: %s', j))
            ## log_info(n_unmatch_msg)
            ## log_info(n_unmatch_tmp)
            n_unmatch_msg <- str_c(n_unmatch_msg, ',', n_unmatch_tmp)
            ## log_info(n_unmatch_msg)
        }
        
        ## silent
        if (sum(observed_aa ==  designed_aa) > 0) {
            n_silent <- n_silent + 1
            n_silent_tmp <- str_c(posN, observed_aa, '(', observed_triplet, ')')
            n_silent_msg <- str_c(n_silent_msg, n_silent_tmp, collapse=',')
        }
        j <- j + 1
    }

    ## create return msg
    if (n_silent == 0){
        msg <- 'Unsupposed mutation detected(N).'
    } else if (n_unmatch == n_silent){
        msg <- 'Silent mutation detected(N).'
    } else {
        msg <- 'Unsupposed mutation detected(N)., Silent mutation detected(N).'
    }
    
    ret <- list(msg=msg,
                n_unmatch=n_unmatch,
                n_unmatch_msg=n_unmatch_msg,
                n_silent=n_silent,
                n_silent_msg=n_silent_msg)
    return(ret)
}

checkUnmatchMutation <- function(dfa){

    log_info(sprintf('%s started.', 'checkUnmatchMutation()'))

    dfa %>%
        filter(Primary!='-') %>%
        filter(Template != Primary) %>%
        filter(Template!='N') -> diff
    
## # A tibble: 1 × 13
##     Pos Primary Secondary Template  Posa PosString CompareNt TemplateTriplet PrimaryTriplet TemplateAA PrimaryAA PosStringAA PosStringFat
##   <dbl> <chr>   <chr>     <chr>    <dbl> <chr>     <chr>     <chr>           <chr>          <chr>      <chr>     <chr>       <chr>       
## 1   325 T       T         C          109 C325T     different CGT             TGT            R          C         R109C       C325T(R109C)

    n_unmatch <- nrow(diff)
    
    ## no unmatch
    if (n_unmatch == 0 ){
        n_silent <- NA
        ret <- list(msg='OK',
                    n_unmatch=NA,
                    n_unmatch_msg=NA,
                    n_silent=NA,
                    n_silent_msg=NA)
        
    ## unmatch detected.
    } else {
        ## number of unsupposed mutation
        diff %>%
            .$PosStringFat -> unmatch.pos

        ## check silent
        diff %>%
            filter(TemplateAA==PrimaryAA) -> diff.silent

        n_silent <- nrow(diff.silent)
        diff.silent %>%
            .$PosStringFat -> unmatch.silent.pos

        ## set msg
        if (n_silent == 0){
            msg <- 'Unsupposed mutation detected.'
        } else if (n_unmatch == n_silent){
            msg <- 'Silent mutation detected.'
        } else {
            msg <- 'Unsupposed mutation detected., Silent mutation detected.'
        }
        
        ## Send message (where different peak detected)
        pos.string <- str_c(unmatch.pos, collapse=',')
        pos.silent.string <- NA
        if (n_silent > 0) {
            pos.silent.string <- str_c(unmatch.silent.pos, collapse=',')    
        }
        
        ret <- list(msg=msg,
                    n_unmatch=n_unmatch,
                    n_unmatch_msg=pos.string,
                    n_silent=n_silent,
                    n_silent_msg=pos.silent.string)
    }
    return(ret) 
}


## ---------------------------
## check2ndBaseCall()
##
## check secondary base call
## ---------------------------
check2ndBaseCall <- function(df){

    log_info(sprintf('%s started.', 'check2ndBaseCall()'))

    nseq <- nchar(df$Sequence[1])

    ## separate, transpose
    df %>%
        transposeSeqDf() %>%
        mutate(PosString=str_c(primary, Pos, secondary)) %>%
        mutate(compare=case_when(
                   primary==secondary ~ 'identical',
                   TRUE ~ 'different'
               )) -> df.t

    ## differnt sequence
    df.t %>%
        filter(compare=='different') -> diff        
    
    ## no contamination
    if (nrow(diff) == 0 ){
        ret <- list(msg='OK',
                    n_2nd_basecall=NA,
                    n_2nd_basecall_msg=NA,
                    n_2nd_basecall_N= NA,
                    n_2nd_basecall_N_msg=NA)
    ## secondary base called.
    } else {
        ## Position primary seq and secondary seq is different
        df.t %>%
            filter(compare =='different') %>%
            .$PosString -> different.pos
        
        ## Template sequence where primary and secondary are different 
        df.t %>%
            filter(compare =='different') %>%
            .$template -> different.template

        ## Different positions where template is 'N': mutation permitted
        n_2nd_basecall_N_msg <- different.pos[different.template == 'N']
        n_2nd_basecall_N <- length(n_2nd_basecall_N_msg)
        if (n_2nd_basecall_N == 0){
            n_2nd_basecall_N_msg <- NA
        }
        
        ## Different positions where template is '[A|T|C|G]': mutation not permitted
        n_2nd_basecall_msg <- different.pos[different.template != 'N']
        n_2nd_basecall <- length(n_2nd_basecall_msg)
        if (n_2nd_basecall == 0){
            n_2nd_basecall_msg <- NA
        }

        ## Send message (where different peak detected)
        if (n_2nd_basecall_N > 0 & n_2nd_basecall > 0){
            if (allow2ndBasecallNotN == TRUE){
                msg <- 'Secondary basecall detected(N).'
            } else {
                msg <- 'Secondary basecall detected(N/not N).'
            }

        } else if (n_2nd_basecall_N ==0  & n_2nd_basecall > 0){
            if (allow2ndBasecallNotN == TRUE){
                msg <- 'OK'
            } else {
                msg <- 'Secondary basecall detected(not N).'
            }
        } else if (n_2nd_basecall_N > 0 & n_2nd_basecall == 0){
            msg <- 'Secondary basecall detected(N).'
        }
        ret <- list(msg=msg,
                    n_2nd_basecall=n_2nd_basecall,
                    n_2nd_basecall_msg=str_c(n_2nd_basecall_msg, collapse=','),                    
                    n_2nd_basecall_N=n_2nd_basecall_N,
                    n_2nd_basecall_N_msg=str_c(n_2nd_basecall_N_msg, collapse=','))
    }
    return(ret)
}

## -------------------------------------------------------------------------------



checkDegeneratedSymbol <- function(consensussequence){
    log_info('checkDegenatedSymbol()')

    DEGENERATED_SYMBOLS <- names(IUPAC_CODE_MAP)[5: length(IUPAC_CODE_MAP)]
    ##  [1] "M" "R" "W" "S" "Y" "K" "V" "H" "D" "B" "N"

    tibble(Sequence=str_split(consensussequence, pattern='')[[1]]) %>%
        tibble::rowid_to_column() %>%
        rename(Pos=rowid) %>%
        mutate(PosString=str_c(Pos, Sequence)) %>%
        filter(Sequence %in% DEGENERATED_SYMBOLS) -> deg

    msg <- ''
    ## no degenerated symbols
    if (nrow(deg) == 0){
        ret <- list(msg='OK',
                    n_degenerated_symbol=NA,
                    n_degenerated_symbol_msg=NA)
    ## degenerated symbols detected
    } else {
        ## detected symbols
        deg %>%
            .$PosString %>%
            str_c(collapse=',') -> posstring
        
        ## Send message
        msg <- str_c(msg, posstring)
        
        ret <- list(msg='Degenerated symbols detected.',
                    n_degenerated_symbol=nrow(deg),
                    n_degenerated_symbol_msg=msg)
    }
    return(ret)
}

