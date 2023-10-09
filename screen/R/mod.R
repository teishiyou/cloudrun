M2inside_calculate_trimming <-function(qualityPhredScores,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize) {
    rawSeqLength <- length(qualityPhredScores)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    if (M2SlidingWindowSize > 40 || M2SlidingWindowSize < 0 ||
        M2SlidingWindowSize%%1!=0 ||
        M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
        M2CutoffQualityScore%%1!=0) {
        trimmedStartPos = NULL
        trimmedFinishPos = NULL
    } else if (rawSeqLength <= M2SlidingWindowSize){
        trimmedStartPos = NULL
        trimmedFinishPos = NULL
    } else {
        ### ------------------------------------------------------------------------
        ### Find the trimming start point
        ###  First window that the average score is bigger than threshold score.
        ###  (Whole window will be kept)
        ### ------------------------------------------------------------------------
        totalThresholdScore <- M2CutoffQualityScore * M2SlidingWindowSize
        trimmedStartPos = 1
        trimmedFinishPos = 2
        for (i in seq_len((rawSeqLength-M2SlidingWindowSize+1))) {
            totalScore <-
                sum(qualityPhredScores[i:(i+M2SlidingWindowSize-1)])
            if (totalScore > totalThresholdScore) {
                trimmedStartPos = i
                break
            }
        }
        qualityPhredScoresRev <- rev(qualityPhredScores)
        for (i in seq_len((rawSeqLength-M2SlidingWindowSize+1))) {
            totalScore <-
                sum(qualityPhredScoresRev[i:(i+M2SlidingWindowSize-1)])
            if (totalScore > totalThresholdScore) {
                trimmedFinishPos = i
                break
            }
        }
        trimmedFinishPos <- length(qualityPhredScoresRev) - trimmedFinishPos + 1
        # for (i in (trimmedStartPos+M2SlidingWindowSize-1):(rawSeqLength-M2SlidingWindowSize+1)) {
        #     totalScore <-
        #         sum(qualityPhredScores[i:(i+M2SlidingWindowSize-1)])
        #     if (totalScore < totalThresholdScore) {
        #         # Keep all base pairs in the previous window.
        #         trimmedFinishPos = i + M2SlidingWindowSize - 2
        #         break
        #     }
        # }
        if (trimmedStartPos == (rawSeqLength-M2SlidingWindowSize+1) ||
            trimmedStartPos == trimmedFinishPos) {
            trimmedStartPos = 1
            trimmedFinishPos = 2
        }
        trimmedSeqLength = trimmedFinishPos - trimmedStartPos
        trimmedQualityPhredScore <-
            qualityPhredScores[(trimmedStartPos+1):trimmedFinishPos]
        trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
        trimmedMinQualityScore <- min(trimmedQualityPhredScore)
        remainingRatio = trimmedSeqLength / rawSeqLength
    }
    return(list("rawSeqLength" = rawSeqLength,
                "rawMeanQualityScore" = rawMeanQualityScore,
                "rawMinQualityScore" = rawMinQualityScore,
                "trimmedStartPos" = trimmedStartPos,
                "trimmedFinishPos" = trimmedFinishPos,
                "trimmedSeqLength" = trimmedSeqLength,
                "trimmedMeanQualityScore" = trimmedMeanQualityScore,
                "trimmedMinQualityScore" = trimmedMinQualityScore,
                "remainingRatio" = remainingRatio))
}
