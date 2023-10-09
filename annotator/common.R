require(Rsamtools)
library(parallel)
library(randomForest)
library(Biostrings)
library(seqinr)
library(TFBSTools)
library(GenomicRanges)
library(IRanges)
require("JASPAR2014")
library(parallel)
library(randomForest)



#library(BSgenome)



locateMotif <- function(seq,  min.score="50%" )
{
    motif_Jaspar_id = getOption("jasparID")
    pfm     <- getMatrixByID(JASPAR2014, ID = motif_Jaspar_id)
    pwm     <- toPWM(pfm) # convert to position weighted matrix
    subject <- DNAString(seq)
    siteset <- searchSeq(pwm, subject, seqname="seq1", min.score=min.score, strand="*")
    siteset =  as(siteset, "data.frame")
    cols = c("start", "end", "absScore", "strand")
    siteset[which.max(siteset$absScore), cols, drop =F]

}

rangeOverlap <- function(tab1, tab2)
{
    query <- with(tab1, GenomicRanges::GRanges(seqnames = paste0("chr", V1),ranges=IRanges::IRanges(start=V2,end=V3))  )
    subject <- with(tab2, GenomicRanges::GRanges(seqnames = paste0("chr", V1),ranges=IRanges::IRanges(start=V2,end=V3))  )
    idx = GenomicRanges::findOverlaps(query, subject, type="any", ignore.strand=TRUE , select="first")
    ifOverlapped = !is.na(idx)
    ifOverlapped

}

getJaspar <- function(dnaSeqs) {
    ncpus = getOption("ncpus") 
    if(ncpus >1)
    { 
        res     = mclapply(dnaSeqs, function(baseStr){
            locateMotif(baseStr)
    }, mc.cores = ncpus)
    } else
    {

        res     = lapply(dnaSeqs, function(baseStr){
            locateMotif(baseStr)
        })
    }

    res=do.call(rbind.data.frame, res)
    res
}

getDeepBind <- function(dnaSeqs, tool) {
    ncpus = getOption("ncpus") 
    deepBindIDFile = getOption("deepBindID") 
    if(ncpus > 1)
    {

        res     = mclapply(dnaSeqs, function(baseStr){
            cmd = sprintf("%s %s < <(echo '%s')",tool, deepBindIDFile, baseStr)
            cmd = paste("/bin/bash -c", shQuote(cmd))
            res = read.table(pipe(cmd),header =F,skip =1)
            res = as.numeric(res)
            dim(res) = c(4, length(res)/4) # score, coreLoc, size, strand
            res = t(res)
            res[,2] = res[,2] + res[,3]/2
            res
        },mc.cores=ncpus)
    } else
    {
        res    = lapply(dnaSeqs, function(baseStr){
            cmd = sprintf("%s %s < <(echo '%s')",tool, deepBindIDFile, baseStr)
            cmd = paste("/bin/bash -c", shQuote(cmd))
            res = read.table(pipe(cmd),header =F,skip =1)
            res = as.numeric(res)
            dim(res) = c(4, length(res)/4) # score, coreLoc, size, strand
            res = t(res)
            res[,2] = res[,2] + res[,3]/2
            res
        })

    }
    do.call(rbind, res)
}

center <- function(tab) {
    mid = with(tab, (V2+V3)/2)
    tab$V2 = mid -400; tab$V3 = mid +400
    tab
}


getDistToTAD <- function(tab, file, span = 10000)
{
    hiC_TAD_F   =  file
    tab_TAD=  read.table(hiC_TAD_F, header =F, as.is = T)
    colnames(tab_TAD)[1:3] = c("V1", "V2", "V3")
    tab_TAD_leftBound =  data.frame(V1 = tab_TAD[[1]],
                            V2=tab_TAD[[2]]-span,
                            V3=tab_TAD[[2]]+span, stringsAsFactors =F)
    tab_TAD_rightBound =  data.frame(V1 = tab_TAD[[1]],
                            V2=tab_TAD[[3]]-span,
                            V3=tab_TAD[[3]]+span, stringsAsFactors =F)
    map_chr = c(1:24); names(map_chr) = c(1:22,"X", "Y")
    toSimplePos <- function(pos) {
        simplePos = pos[1] * 1e6 + (pos[2] + pos[3])/2 / 1000
        simplePos
    }
    simplePos1 = apply(with(tab, data.frame(map_chr [ V1 ], V2, V3)), 1, toSimplePos)
    simplePos2 = apply(with(tab_TAD_leftBound, data.frame(map_chr [ V1 ], V2, V3)), 1, toSimplePos)
    simplePos3 = apply(with(tab_TAD_rightBound, data.frame(map_chr [ V1 ], V2, V3)), 1, toSimplePos)
    findClosest <- function(v, ll, strand = c()) {
        idx = which.min(abs(v-ll))
        if(length(strand) ==0) strand = rep(1, length(ll))
        (v - ll[idx]) *strand[idx]
    }
    dist_to_LBoundary = sapply (simplePos1, findClosest, ll = simplePos2)
    dist_to_RBoundary = sapply (simplePos1, findClosest, ll = simplePos3)
    cbind(TAD_left_bound = dist_to_LBoundary, TAD_right_bound = dist_to_RBoundary)
}
getGC <- function(fullTab,ref)
{
    offset = 80
    outSeqFile = sprintf("offset.%d.fastq", offset)
    V2 = with(fullTab,V2 + J_start)
    V3 = with(fullTab,V2 + J_end)
    grtest <- GRanges(seqnames =paste0("chr",fullTab$V1), ranges=IRanges(start=V2-offset,end=V2+20+offset))
    myset = readDNAStringSet(ref)
    mySeq = getSeq(myset, grtest)
    names(mySeq) = paste(fullTab$V1, V2,V3, sep = "-")
    writeXStringSet(mySeq, filepath=outSeqFile, compress=FALSE, compression_level=NA, format="fasta")
    rm(mySeq); rm(myset)
    getGC <- function(file) {
        seqs <- read.fasta(file = file)
        gc = lapply(seqs,GC) # calling the GC funciton from seqinr Package
        unlist(gc)
    }
    getGC(outSeqFile)
}


count_cpg <- function(dna_seq) {
  # Convert the sequence to uppercase to simplify pattern matching
  dna_seq <- toupper(dna_seq)
  # Count the number of CpG islands in the sequence
  num_cpg <- length(gregexpr("CG", dna_seq)[[1]])
  # Return the number of CpG islands
  return(num_cpg)
}

hasCpG <- function(fullTab, ref)
{
    offset = 80
    outSeqFile = sprintf("offset.%d.fastq", offset)
    V2 = with(fullTab,V2 + J_start)
    V3 = with(fullTab,V2 + J_end)
    grtest <- GRanges(seqnames =paste0("chr",fullTab$V1), ranges=IRanges(start=V2,end=V3))
    myset = readDNAStringSet(ref)
    mySeq = getSeq(myset, grtest)

    nCpG = sapply(mySeq, count_cpg)
    nCpG

}






getUpMotif <-function(fullTab, ref, offset = 20, ncpus=4)
{
    locateMotif <- function(seq, min.score="50%", pwm,strand )
    {
        subject <- DNAString(seq)
        siteset <- searchSeq(pwm, subject, seqname="seq1", min.score=min.score, strand=strand)
        siteset =  as(siteset, "data.frame")
        if(nrow(siteset) == 0 ) return (  data.frame(start=NA,end=NA,absScore=NA,strand = NA) )
        cols = c("start", "end", "absScore", "strand")
        siteset[which.max(siteset$absScore), cols, drop =F]
    }

    getJaspar <- function(dnaSeqs,strand) {
        motif = getOption("jaspar_upMotif")
        pfm = readJASPARMatrix(motif, matrixClass="PFM")[[1]]
        pwm     <- toPWM(pfm)
        res     = mclapply(dnaSeqs, function(baseStr){
            locateMotif(baseStr, pwm=pwm,strand=strand)
        }, mc.cores=ncpus)
        res=do.call(rbind.data.frame, res)
        res
    }

    fullTab$upMotif_score = NA; fullTab$upMotif_start = NA
    fullTab$upMotif_end   = NA; fullTab$upMotif_strand = NA
    tab_orig = fullTab
    for ( targetStrand in c( "+", "-") )
    {
        tab = tab_orig
        tab = subset(tab, J_strand == targetStrand)
        V2  = with(tab,V2 + J_start); V3  = with(tab,V2 + J_end)
        if(targetStrand == "+")
            grtest <- GRanges(seqnames =paste0("chr",tab$V1), ranges=IRanges(start=V2-offset,end=V2), strand =tab[, "J_strand"])
        if(targetStrand == "-")
            grtest <- GRanges(seqnames =paste0("chr",tab$V1), ranges=IRanges(start=V3,end=V3+offset), strand =tab[, "J_strand"])
        dnaSeqs    = scanFa(ref,grtest)
        dnaSeqs    = as.character(dnaSeqs)
        tab_jaspar = getJaspar(dnaSeqs, strand = targetStrand)
        tab_orig[tab_orig$J_strand== targetStrand, ]$upMotif_score = tab_jaspar$absScore
        tab_orig[tab_orig$J_strand== targetStrand, ]$upMotif_start = tab_jaspar$start
        tab_orig[tab_orig$J_strand== targetStrand, ]$upMotif_end   = tab_jaspar$end
        tab_orig[tab_orig$J_strand== targetStrand, ]$upMotif_strand   = tab_jaspar$strand
    }
    tab_orig
}

