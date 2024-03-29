## ============================================================================
## The motifProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome getBSgenome
#' @importFrom GenomicRanges start

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Transfer the input motif
.transMotif <- function(motif)
{
    ## IUPAC codes
    IUPAC <- data.frame(IUPAC=c("A", "C", "G", "T", "R", "Y", "S", "W",
                            "K", "M", "B", "D", "H", "V", "N"),
                        code=c("[A]", "[C]", "[G]", "[T]", "[AG]", "[CT]",
                            "[GC]", "[AT]", "[GT]", "[AC]", "[CGT]",
                            "[AGT]", "[ACT]", "[ACG]", "[ATCG]"))

    ## Transfer the input motif
    splitmotif <- strsplit(motif,"") %>% as.data.frame()
    splitmotif$code <- IUPAC$code[match(splitmotif[,1], IUPAC$IUPAC)]
    output <- paste0(splitmotif$code, collapse = '')
    output <- paste0("^", output)
    return(output)
}

## Calculation as fraction

.getFractionPlot <- function(flankingSeq, l, checkmotif, flanking, motif, title)
{
    number <- c()
    for (j in seq_len(l)){
        sub <- substr(flankingSeq,j,j+nchar(motif)-1)
        P <- length(sub[grep(checkmotif, sub)])/
            length(sub)
        number <- c(number,P)
    }

    ## Making dataframe for the plot
    df <- data.frame(Position=seq(from=-(flanking+nchar(motif)-1),
        to=flanking), Fraction=number)

    df$Position <- factor(df$Position,
        levels=seq(from=-(flanking+nchar(motif)-1), to=flanking))
    df <- df[-seq_len(nchar(motif)-1),]

    p1 <- ggplot(df, aes(x=Position, y=Fraction)) +
        geom_bar(stat = "identity", fill = "Orange") + ggtitle(title) +
        theme_bw() + ylab(paste0("Fraction of ", motif)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

    out <- list(Numbers=df, Plot=p1)

    return(out)
}

## Making the plot with number
.getNumberPlot <- function(flankingSeq, l, checkmotif, flanking, motif, title)
{
    number <- c()
    for (j in seq_len(l)){
        sub <- substr(flankingSeq,j,j+nchar(motif)-1)
        P <- length(sub[grep(checkmotif, sub)])
        number <- c(number,P)
    }

    ## Making dataframe for the plot
    df <- data.frame(Position = seq(from=-(flanking+nchar(motif)-1),
                                    to=flanking),
        Number = number)

    df$Position <- factor(df$Position,
        levels=seq(from=-(flanking+nchar(motif)-1), to=flanking))
    df <- df[-seq_len(nchar(motif)-1),]

    p1 <- ggplot(df, aes(x=Position, y=Number)) +
        geom_bar(stat = "identity", fill = "Orange") + ggtitle(title) +
        theme_bw() + ylab(paste0("Number of ", motif)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

    out <- list(Numbers=df, Plot=p1)
    return(out)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "motifProfile" methods for GRanges objects.
##

#' @title motifProfile for the GRanges objects
#'
#' @description An function to plot the frequency or fraction of the interested
#'              motif around the center of input peaks.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param motif A character string which use the IUPAC nucleotide code, e.g.
#'              DRACH, TTAGGG.
#' @param genome The name of the full genome sequences package in the
#'               Bioconductor, e.g. "BSgenome.Mmusculus.UCSC.mm10". You should
#'               install the package before running this function.
#' @param fraction A logical vector (TRUE or FALSE) that the result should be
#'                 presented in fraction or number.
#' @param title The main title for the output meta gene profile plot.
#' @param flanking The size of the flanking windows that you would like to
#'                 check. Flanking=5 will give you the result of the 10+1nt
#'                 windows around the center of peaks.
#'
#' @return A list object, the list 1 contains the information of the
#'         frequency of specified motif around the center of peaks. The list 2
#'         includes the plot of motifProfile.
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' ## Please make sure that the correct BSgenome package have installed before
#' ## running motifProfile. For example,library("BSgenome.Mmusculus.UCSC.mm10")
#' ## would be required for the mouse data.
#'
#' output <- motifProfile(test,
#'   motif = "DRACH",
#'   genome = "BSgenome.Mmusculus.UCSC.mm10",
#'   flanking = 20
#' )
#' @export
#'

motifProfile <- function(object, motif=NA, genome=NA, fraction=TRUE,
    title="Motif Profile", flanking=10)
    {
        if (!is.character(motif))
            stop("The 'motif' should be a character string.")
        if (missing(motif))
            stop("The parameter 'motif' is missing.")
        if (missing(object))
            stop("The input object is missing.")
        if (missing(genome))
            stop("The parameter genome is missing.")
        if (nchar(motif) > 1+flanking)
            stop("The motif is longer than the flanking region. Please
            set a bigger number to the parameter flanking.")
        else {
            object <- .centerPeaks(object)

            ## Transfer input into uppercase
            motif <- toupper(motif)
            splitString <- strsplit(motif,"") %>% as.data.frame()

            ## Make sure the motif is using the IUPAC code
            IUPAC <- c("A", "T", "C", "G", "R", "Y", "S", "W", "K",
                "M", "B", "D", "H", "V", "N")
            if (sum(!splitString[,1] %in% IUPAC) != 0)
                stop("The 'motif' should be the IUPAC nucleotide code.")

            ## Transfer the IUPAC code
            checkmotif <- .transMotif(motif)
            genome <- getBSgenome(genome)
            flankingSeq <- Biostrings::getSeq(genome,
                object+flanking + nchar(motif) - 1)

            ## Calculate the number of iteration for the string splite
            l <- 1+((flanking+nchar(motif)-1)*2)-nchar(motif)+1

            ## Calculate the result and generate plot
            if (fraction == TRUE) {
                output <- .getFractionPlot(flankingSeq, l, checkmotif,
                    flanking, motif, title)
            }
            if (fraction == FALSE) {
                output <- .getNumberPlot(flankingSeq, l, checkmotif,
                    flanking, motif, title)
            }

            return(output)
        }
    }
