## *************************************************************
## Functions
## *************************************************************


## split position of aa(space splitted string) to integer vector
splitPosa <- function(posa){
    posa <- str_replace(posa, pattern='^\\s+', '')
    posa <- str_replace(posa, pattern='\\s+$', '') 
    posa <- str_replace_all(posa, pattern='\\s+', ' ')
    posa <- as.integer(unlist(str_split(posa, pattern=' ')))   ## Postions in amino acids
    return(posa)
}

## letter
translateA <- function(ltr){
    ltr <- unname(GENETIC_CODE[ltr])
    if (is.na(ltr)){ltr <- '-' }
    return(ltr)
}

## sequence
translateA.seq <- function(seq){
    seq <- rev(rev(str_split(gsub('(.{3})', '\\1 ', seq), ' ')[[1]])[-1])
    seq <- str_c(unlist(map(seq, translateA)), collapse = '')
    return(seq)
}


## Dataframe to AAStringSet
dataframe2as <- function(df){
    as <- AAStringSet(df$Sequence)
    names(as) <- df$Name
    return(as)
}

## Dataframe to DNAStringSet
dataframe2ds <- function(df){
    ds <- DNAStringSet(df$Sequence)
    names(ds) <- df$Name
    return(ds)
}

dnastringset2df <- function(ds){
    name <- names(ds)
    sequence <- paste(ds)
    df <- as_tibble(data.frame(Name=name, Sequence=sequence))

    return(df)
    
}

## ------------------------------------
## aa to triplet
## ------------------------------------
## https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=83333
## DFCODON <- read_csv('/home/docker/Lab/codon//codon.83333/codon.ecoli.most.frequent.csv')
aa2triplet <- function(ltr){
    trp <- DFCODON$triplet[DFCODON$aa == ltr]  ## DFCODON: global variable
    if (length(trp)==0){
        trp <- '---'
    }
    return(trp)
}

aa2triplet.seq <- function(aa){
    seq <- str_split(aa, '')[[1]]
    seq <- str_c(unlist(map(seq, aa2triplet)), collapse='')
    return(seq)
}
## ------------------------------------

str_split.my <- function(string){
    str_split(string, pattern='')[[1]]
}


## from
##  [1] "E" "A" "L" "E" "L" "R" "A" "S" "E" "V" "Q" "R" "D" "K" "I" "A" "F" "T"
## to
##  [1] "GAA" "GCA" "CTG" "GAA" "CTG" "CGT" "GCA" "AGC" "GAA" "GTT" "CAG" "CGT"
## [13] "GAT" "AAA" "ATT" "GCA" "TTT" "ACC"
splitaa2triplet <- function(spltaa){
    unlist(map(spltaa, aa2triplet))
}


## mutated aa letter -> full aa sequence
## posa : mutated position (ex: c(23, 31, 35, 39, 43, 67, 82, 99, 101, 107, 122, 162, 167, 169, 173, 195, 226, 244))
## mseq : mutation aa seq (ex: 'EALELRASEVQRDKIAFT')
## org.a: original aa full sequence
subLetterA <- function(posa, mseq, org.a){
    mltr <- str_split(mseq, '')[[1]]
    seq <- str_split(org.a, '')[[1]]
    seq[posa] <- mltr
    seq <- str_c(seq, collapse = '')
    return(seq)
}


## sequencent: 'AAATTTCCCGGGAAATTT'
## rtn       : "AAA" "TTT" "CCC" "GGG" "AAA" "TTT"
sequencent2triplet <- function(sequencent){
    if (nchar(sequencent) %% 3 == 0){
        rtn <- rev(rev(str_split(gsub('(.{3})', '\\1 ', sequencent), ' ')[[1]])[-1])
        return(rtn)
    } else {
        print(sprintf('sequencent: %s', sequencent))
        stop('nchar of sequencent should be multiple of 3.')
        
    }
}

## mutated nt letter -> full nt sequence
## Caution!! it's not nt position but aa position.
## posa   : mutated position aa
##           ex.):[1]  23  31  35  39  43  67  82  99 101 107 122 162 167 169 173 195 226 244    
## spltnt : splitted nt seq. length is identical to posa.
##           ex.) [1] "CAG" "GCG" "GCG" "GAA" "ATT" "CGC" "GCG" "AGC" "GAA" "GAT" "GGC" "GGC"
##               [13] "GAT" "AAA" "TTT" "GGC" "TAT" "ACC"
## org    : original(wt) nt full sequence
subLetter <- function(posa, spltnt, org){
    ## seq <- rev(rev(str_split(gsub('(.{3})', '\\1 ', org), ' ')[[1]])[-1]) ## [1] "Add" "Ajj" "Dkk" "Nll"
    seq <- sequencent2triplet(org)
    ##seq[posa] <- spltnt

    ## Compare with wt aa
    mseq.s <- unlist(map(spltnt, translateA))
    wtseq.s <- str_split(translateA.seq(org), '')[[1]][posa]
    idx <- (mseq.s == wtseq.s) ## FALSE: sholud be substisuted, TRUE: wt original nt sequence should be retained.

    ## posa.sub: position should be substituted.
    ## ex).
    ## posa    :[1]  23  31  35  39  43  67  82  99 101 107 122 162 167 169 173 195 226 244    
    ## wtseq.s :[1] "E" "A" "L" "E" "L" "R" "A" "S" "E" "V" "Q" "R" "D" "K" "I" "A" "F" "T"
    ## mseq.s  :[1] "Q" "A" "A" "E" "I" "R" "A" "S" "E" "D" "G" "G" "D" "K" "F" "G" "Y" "T"
    ## idx     :[1]  F   T   F   T   F   T   T   T   T   F   F   F   T   T   F   F   F   T
    ## posa.sub:[1]  23      35      43                 107 122 162         173 195 226 
    posa.sub <- posa[!idx] 
    seq[posa.sub] <- spltnt[!idx]
    
    seq <- str_c(seq, collapse='')
    return(seq)
}

### ----------------------------------------------------------------------------
### count dash / plus
###   - dash means sequence contains deletion
###   - plus means sequence is generated from reads contains insertion/deletion
### ----------------------------------------------------------------------------
count_inside_dash <- function(Sequence){
    ## this function counts '-' in sequences
    num <- sum(nchar(str_extract_all(Sequence, '(?<=[A|T|G|C|N])-+(?=[A|T|G|C|N])')[[1]]))
    return(num)
}

count_dash <- function(Sequence){
    ## this function counts '-' in sequences
    num <- sum(nchar(str_extract_all(Sequence, '-')[[1]]))
    return(num)
}

count_plus <- function(Sequence){
    ## this function counts '-' in sequences
    num <- length(str_extract_all(Sequence, '\\+')[[1]])
    return(num)
}

### ----------------------------------------------------------------------------
### cliporf
### ----------------------------------------------------------------------------
cliporf <- function(sequence, orfstart, orfend){
    sequence <- DNAStringSet(sequence)

    ## clip orf
    match <- vmatchPattern(as.character(orfstart), sequence)
    s <- match[[1]]@start
    match <- vmatchPattern(as.character(orfend), sequence)
    e <- match[[1]]@start + match[[1]]@width - 1

    ## when multiple ends hits 
    if (length(e)>1){
        e <- e[1]
    } else if (length(e)==0) {
        return(-1)
    }
     
    sequence <- subseq(sequence, s, e)

    return(sequence)
}

### ----------------------------------------------------------------------------
### generate ApEformat File from fasta
### ----------------------------------------------------------------------------
genApeFormat <- function(query){

    ## Specify ape template file
    template <- read_lines('./ApE_Template.fa')
    
    ## remove dash if exists
    sequence <- str_replace_all(as.character(query), '-', '')

    ## split into 60 characters/line
    each60 <- str_split(gsub('(.{60})', '\\1 ', sequence), pattern=' ')[[1]]

    ## split each lines into 10 characters, then remove last space
    sequence_lines <- str_replace(gsub('(.{10})', '\\1 ', each60), ' $', '')

    ## make header
    n <- length(sequence_lines)
    header <- (seq(1, n) -1) * 60 + 1
    header <- str_pad(header, 9, 'left')

    ## concat header and sequence
    sequence_lines <- paste(header, sequence_lines)

    ## get orf start and end pos
    start <- vmatchPattern(as.character(palorfstart), DNAStringSet(sequence))[[1]]@start
    match <- vmatchPattern(as.character(palorfend), DNAStringSet(sequence))[[1]]
    end <- match@start + match@width - 1

    ## replace name
    res <- str_replace_all(template, 'RK003_CL_A1_50', names(query))

    ## replace length
    len <- sprintf('%s bp', nchar(sequence))
    res <- str_replace(res, '2235 bp', len)

    ## replace date
    today <- str_to_upper(format(Sys.Date(), format='%d-%b-%Y'))
    res <- str_replace(res, '05-AUG-2022', today)

    ## replace orf position
    res <- str_replace(res, '53..2218', str_c(start, '..', end))

    ## replace sequence body
    res <- c(res[1: 18], sequence_lines, res[57])

    return(res)
}



getSequenceStart <- function(sequence){

    first_letter <- str_split(sequence, pattern='')[[1]][1]

    ## 1st letter is [A|T|G|C]
    if (first_letter != '-') {
        start <- 1

    ## 1st letter is '-'
    } else {
        ## get template start position(find 1st letter position)    
        start <- str_locate(sequence, pattern='-[A|T|G|C]')[[1]] + 1
    }

    ## when no match
    if (is.na(start)){ start <- 1 }

    return(start)
    
}

getSequenceEnds <- function(sequence){

    ## get start pos
    start <- getSequenceStart(sequence)

    ## get end pos
    s <- rev(str_split(sequence, '')[[1]])          ## reverse sequence(array)
    s <- str_c(s, collapse = '')                    ## reverse sequence(one-string)
    
    end <- getSequenceStart(s)
    end <- nchar(s) - end + 1
    
    ends <- list(start=start, end=end)
    
    return(ends)
}

## get 3 bases from specified posa
get3bases <- function(sequence, pos){
    s <- str_split(sequence, '')[[1]][pos: (pos+2)]
    s <- str_c(s, collapse='')
    return(s)
} 

get3basesA <- function(sequence, posa){

    ## print(sprintf('sequence: %s', sequence))
    ## print(sprintf('posa: %s', posa))
    
    if (nchar(sequence) %% 3 != 0){
        return(sequence)
    }
    
    s <- sequencent2triplet(sequence)

    ## print(sprintf('s[posa]: %s', s[posa]))    
    ## s <- str_c(s[posa], collapse='')
    return(s[posa])
} 


## --------------------------------------------------------
## TODO:
## should be more tidy coding.
## --------------------------------------------------------
## getSingleMseq <- function(df, posa){
##     for (i in 1: length(posa)) {
##         cn  <- paste0('Md', i)
##         cna <- paste0('Ma', i)
##         df <- df %>%
##             mutate(!!sym(cn) := pmap(list(sequence=Sequence, posa=posa[i]), get3basesA)) %>%
##             mutate(!!sym(cna) := map(!!sym(cn), translateA))
##     }
##     return(df)
## }
getSingleMseq <- function(df, posa){
    SEQLEN <- length(posa)

    df %>%
        separate(Mseq,
                 into=as.character(sprintf('Ma%s', seq(1, SEQLEN))),
                 sep=seq(1, SEQLEN), remove=FALSE) %>%
        ## following code is not appropriate because it does not use raw sequence
        ## mutate(Mseqnt = map(Mseq, aa2triplet.seq)) %>% 
        mutate(Mseqnt=pmap(list(sequence=Sequence, posa=list(posa)), getMseqNt)) %>%
        separate(Mseqnt,
                 into=as.character(sprintf('Md%s', seq(1, SEQLEN))),
                 sep=seq(1, SEQLEN)*3, remove=FALSE) -> df
    return(df)
}


getMseq <- function(sequence, posa){

    if ((nchar(sequence) %% 3) != 0){
        return(sequence)
    }
    
    tibble(Triplet=sequencent2triplet(sequence)) %>%
    mutate(Mseq=map(Triplet, translateA.seq)) %>%
    unnest(cols=c(Mseq)) %>%
    .$Mseq -> Mseq

    Mseq <- str_c(Mseq[posa], collapse='')

    return(Mseq)

}


getMseqNt <- function(sequence, posa){

    if ((nchar(sequence) %% 3) != 0){
        return(sequence)
    }
    
    tibble(Triplet=sequencent2triplet(sequence)) %>%
    .$Triplet -> Triplet

    MseqNt <- str_c(Triplet[posa], collapse='')

    return(MseqNt)

}

## --------------------------------------------------------------
## template: fasta contains NNN where mutation designged residue
## --------------------------------------------------------------
getPosaFromTemplate <- function(template){

    tibble(Triplet=sequencent2triplet(as.character(template))) %>%
        tibble::rowid_to_column() %>%
        rename(Posa=rowid) %>%
        filter(Triplet=='NNN') %>%
        .$Posa -> posa

    return(posa)
    
}

## getRegion (start:end)in sequence
getRegion <- function(sequence, start, end){
    rgn <- str_split(sequence, pattern='')[[1]][start:end]
    rgn <- str_c(rgn, collapse = '')
    return(rgn)
}




posn2posa <- function(n){
    remainder <- n %% 3
    quotient  <- n %/% 3
    
    if (remainder== 0){
        return(quotient )
    } else {
        return(quotient+1)
    }
}
