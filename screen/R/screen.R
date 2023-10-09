## -----------------------------------------------------------------
## read_design
##
## read codon design file and check format
## 
## design_path     : design path. should be .xlsx
## design_sheet    : sheet name of .xlsx
## output_dir      : design output directory
## -----------------------------------------------------------------
read_design <- function(design_path, design_sheet, output_dir){
    design <- read_excel(design_path)

    ## check design index, stop if N is not sequential
    n_max <- max(design$N)
    idx <- seq(1, n_max)
    if (sum(design$N == idx) != n_max){
        stop('Design check failed. Check your design: unappropriate column N.')
    }
    
    design %>%
        pivot_longer(cols=!N, names_to='Design', values_to='Triplet') %>%
        filter(!is.na(Triplet)) %>%
        mutate(Triplet=map(Triplet, toupper)) %>%
        mutate(AA=map(Triplet, translateA.seq)) %>%
        unnest(cols=c(Triplet, AA)) %>%
        mutate(Qc = case_when(
                   str_detect(AA, '-') ~ 'NG',
                   str_detect(AA, '\\*') ~ 'NG',
                   nchar(Triplet) != 3 ~ 'NG',
                   TRUE ~ 'OK'))  -> design
    
    design_name <- str_extract(design_path, pattern='(?<=/)[^/]+$')
    design_out_path <- sprintf('%s/%s', output_dir, design_name)
    design_out_path <- str_replace(design_out_path, '\\.xlsx', '\\_p.xlsx')
    write.xlsx(design, design_out_path)

    ## check design
    if (sum(design$Qc!='OK') > 0){
        stop(sprintf('Design check failed. Check your design: %s', design_out_path))
    }
    
    return(design)
}


## -----------------------------------------------------------------
## read_flist
##
## input_dir      : data input directory. ex: './data/JC00030617test/'
## stem           : clone name.           ex: 'Dar[0-9]+'
## suffix         : file suffix           ex: '.ab1$'
## forward_string : forward_string.       ex: c('pro', 'promotor')
## reverse_string : reverse_string        ex: c('term', 'terminator')
## 
## -----------------------------------------------------------------
read_flist <- function(input_dir, stem, suffix, forward_string, reverse_string){
    tibble(path=list.files(input_dir, full.names=TRUE, pattern=suffix, recursive=TRUE)) %>%
        mutate(file=str_extract(path, '[^/]+$')) %>%
        mutate(contigName=str_extract(path, stem)) %>%
        mutate(direction=case_when(
                   str_detect(file, paste(forward_string, collapse='|')) ~ 'forward',
                   str_detect(file, paste(reverse_string, collapse='|')) ~ 'reverse',
                   TRUE ~ 'unknown'
               )) -> fList

    return(fList)
}

## -----------------------------------------------------------------
## read_flist2
##
## stem is depreciated.
## string that precedes 'direction strings' is extracted as contigNames
## 
## -----------------------------------------------------------------
read_flist2 <- function(input_dir, stem, suffix, forward_string, reverse_string){

    direction_strings <- c(forward_string, reverse_string)
    
    extractDirection <- function(x){
        x <- str_extract_all(x, direction_strings)
        x <- x[!sapply(x, is.character0)][[1]]
        return(x)
    }
    
    is.character0 <- function(x){
        if (length(x) == 0){
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    tibble(path=list.files(input_dir, full.names=TRUE, pattern=suffix, recursive=TRUE)) %>%
        mutate(file=str_extract(path, '[^/]+$')) %>%
        filter(str_detect(file, str_c(direction_strings, collapse='|'))) %>%
        mutate(direction=map(file, extractDirection)) %>%
        unnest(cols=c(direction)) %>%
        mutate(contigName=str_extract(file, paste0('^.+(?=[_|-]', direction, ')'))) %>%
        mutate(direction=case_when(
                   str_detect(file, paste(forward_string, collapse='|')) ~ 'forward',
                   str_detect(file, paste(reverse_string, collapse='|')) ~ 'reverse',
                   TRUE ~ 'unknown'
               )) -> fList

    return(fList)
}


## ----------------------------------------------------------------------------
## rerun_rule()
##
## exclude samples when sample is re-runned
## ----------------------------------------------------------------------------
rerun_rule <- function(fList){

    fList %>%
        mutate(stem=str_extract(file, '.+(?=\\.ab1$)')) %>%
        mutate(stem=str_replace(stem, '_1', '')) -> fList
    
    fList %>%
        group_by(stem) %>%
        tally() -> n
    
    fList %>%
        left_join(n, by='stem') %>%
        mutate(rerun=case_when(
                   str_detect(file, '_1\\.ab1$') ~ 'Yes',
                   TRUE ~ 'No'
               )) %>%
        filter( n==1 | (n ==2 & rerun == 'Yes' )) -> fList
    
    return(fList)
}


## ----------------------------------------------------------------------------
## mySangerRead()
##
## call SangerRead()
## ----------------------------------------------------------------------------
mySangerRead <- function(readFileName){

    readf <- SangerRead(readFeature           = 'Forward Read',
                        readFileName          = readFileName,
                        geneticCode           = GENETIC_CODE,
                        TrimmingMethod        = TrimmingMethod,
                        M1TrimmingCutoff      = M1TrimmingCutoff,
                        M2CutoffQualityScore  = M2CutoffQualityScore,
                        M2SlidingWindowSize   = M2SlidingWindowSize,
                        baseNumPerRow         = baseNumPerRow,
                        heightPerRow          = heightPerRow,
                        signalRatioCutoff     = signalRatioCutoff,
                        showTrimmed           = TRUE)

    return(readf)
}



## ----------------------------------------------------------------------------
## trimSeq()
## 
## re-trimming with specified offset
## ----------------------------------------------------------------------------
trimSeq <- function(primarySeq, trimmedStartPos, trimmedFinishPos, read5primeOffset, read3primeOffset){
    s <- trimmedStartPos + read5primeOffset
    e <- trimmedFinishPos - read3primeOffset

    if (s >= nchar(primarySeq) | e<=0 | s > e){
        return(primarySeq)
    }

    trimmed <- subseq(primarySeq, s, e)
    return(trimmed)
}

## ----------------------------------------------------------------------------
## getRead()
##
## get read properties as dataframe from sanger read object
## ----------------------------------------------------------------------------
getRead <- function(read, contigName, read5primeOffset, read3primeOffset){
    ## read              : sangerRead objet
    ## contigName        : contigname
    ## read5primeOffset  : 5' offset from trim_startPos
    ## read3primeOffset  : 3' offset from trim_endtPos    
    
    df <- tibble(
        readName                 = read@objectResults@readResultTable$readName,
        readPath                 = read@readFileName,
        trimmedSeqLength         = read@QualityReport@trimmedSeqLength,
        trimmedStartPos          = read@QualityReport@trimmedStartPos,
        trimmedFinishPos         = read@QualityReport@trimmedFinishPos,
        rawMeanQualtiyScore      = read@QualityReport@rawMeanQualityScore,
        trimmedMeanQualityScore  = read@QualityReport@trimmedMeanQualityScore,
        rawMinQualtiyScore       = read@QualityReport@rawMinQualityScore,
        trimmedMinQualityScore   = read@QualityReport@trimmedMinQualityScore,
        remainingRatio           = read@QualityReport@remainingRatio,
        qualityBaseScores        = list(read@QualityReport@qualityBaseScores),
        qualityPhredScores       = list(read@QualityReport@qualityPhredScores),
        rawSeqLength             = read@QualityReport@rawSeqLength,        
        primarySeqRaw            = list(read@primarySeqRaw),
        primarySeq               = list(read@primarySeq),
        secondarySeq             = list(read@secondarySeq),        
        TrimmingMethod           = read@QualityReport@TrimmingMethod,
        M1TrimmingCutoff         = read@QualityReport@M1TrimmingCutoff,
        M2CutoffQualityScore     = read@QualityReport@M2CutoffQualityScore,
        M2SlidingWindowSize      = read@QualityReport@M2SlidingWindowSize,
        read5primeOffset         = read5primeOffset,
        read3primeOffset         = read3primeOffset        
    )

    df %>%
        mutate(contigName               = contigName) %>%
        mutate(direction                = str_replace(read@objectResults@readResultTable$direction, ' Read', '')) %>%
        mutate(trimmedPrimarySequence   = pmap(list(primarySeq,   trimmedStartPos, trimmedFinishPos, read5primeOffset, read3primeOffset), trimSeq)) %>%
        mutate(trimmedSecondarySequence = pmap(list(secondarySeq, trimmedStartPos, trimmedFinishPos, read5primeOffset, read3primeOffset), trimSeq)) %>%
        mutate(rawSequence              = list(read@primarySeqRaw)) %>%
        select(contigName, readName, trimmedPrimarySequence, trimmedSecondarySequence, trimmedSeqLength, everything()) -> df
        
    return(df)
}

