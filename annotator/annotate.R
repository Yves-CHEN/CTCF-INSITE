require("JASPAR2014")
source("./common.R")
library(parallel)
library(randomForest)





#---------------------------------------------------------
#  MAIN
#---------------------------------------------------------

args    = commandArgs(T)

if(length(args)<3) stop("[error] 3 arguments are required: macsF, ncpus, genomeRef")
macsF   = args[1]
ncpus   = as.numeric(args[2])
genomeRef = (args[3])


  options("genomeRef" = "~/work/Data/hg19.fa",
        "deepBindPath" = "./resource/deepbind/deepbind",
        "offsetValue" = 80,
        "minScore" = "50%",
        "jasparID" = "MA0139.1",
        "jaspar_upMotif" = "./resource/upstreamCTCF.jaspar",
        "deepBindID" = "./resource/CTCF.motif.ids",
        "macsFPath" = "input/example.narrowPeaks.gz",
        "model"  = "./tools/RF.model.Oct2022.dat",
        "consitutiveCTCF_file" = "./resource/annotation/hg19/40.celltype.chipseq/constitutive.CTCF.txt",
        "genomicAnnotation_file" = "./resource/annotation/hg19/refGeneExtent.hg19.bed.gz",
        "TAD" = "./resource/annotation/hg19/LNCaP.TAD.40kb.bed.gz",
        "ncpus" = ncpus)




deepBind               =  getOption ("deepBindPath")
ctcfInsite_file        =  getOption ("model")
consitutiveCTCF_file   = getOption ("consitutiveCTCF_file")
genomicAnnotation_file = getOption ("genomicAnnotation_file")
TAD_file               = getOption ("TAD")

resultDir = "results"
dir.create(resultDir)



tab_ctcfCore = read.table(macsF, header =F, as.is = T)
# check if macsF chrID ,which is V1 column,  is in chr1 or 1. If it does not have chr prefix, then add chr
#tab_ctcfCore = addChr(tab_ctcfCore)

tab_ctcfCore = tab_ctcfCore[,c(1,2,3,7)]

sel = tab_ctcfCore$V1 %in% paste0("chr", c(1:22, 'X', 'Y'))
tab_ctcfCore = tab_ctcfCore[sel,]

tab_ctcfCore$len = with(tab_ctcfCore,V3-V2) 
tab_ctcfCore = subset(tab_ctcfCore, len > 80)
len = tab_ctcfCore$len

# center segment >1000bp as deepBind has a limite of 1000
tab_ctcfCore[len>=999,] = center(tab_ctcfCore[len>=999,])
tab_ctcfCore$len = with(tab_ctcfCore,V3-V2) 
tab_ctcfCore$V1 = substring(tab_ctcfCore$V1,4)




gr      = with(tab_ctcfCore, GenomicRanges::GRanges(seqnames = paste0("chr", V1),ranges=IRanges::IRanges(start=V2,end=V3))  )
dnaSeqs = scanFa(genomeRef, gr)
dnaSeqs = as.character(dnaSeqs)
tab_jaspar    = getJaspar(dnaSeqs)
colnames(tab_jaspar) = paste0("J_", colnames(tab_jaspar))
tab_deepBind  = getDeepBind(dnaSeqs, tool = deepBind)
colnames(tab_deepBind  ) = c("chip.score", "D_loc", "D_motif_len","D_strand")




# --------------------------------------
#    constitutive.CTCF.txt: is a set of commonly found CTCFs across 40 cell types. Refer to
#        Maurano, "Role of DNA methylation in modulating transcription factor occupancy."
#         Cell reports 12.7 (2015): 1184-1195
# --------------------------------------
print("----------------------------")
tab_common_CTCF  = read.table(consitutiveCTCF_file, header = T, as.is = T)
colnames(tab_common_CTCF  )[1:3] = c("V1", "V2", "V3")
tab_common_CTCF$V1 = (substring(tab_common_CTCF$V1,4))
ifConstitutive = rangeOverlap(tab_ctcfCore,tab_common_CTCF)

#---------------------------------------------------
#   Exon, 5-UTR, promotor, intron
#---------------------------------------------------
checkDomainOverlap <- function(tab_ctcf, tab_annot, domain)
{
    tab_annot= subset(tab_annot, V7 == domain)
    rangeOverlap(tab_ctcf, tab_annot)
}
tab_annot = read.table(genomicAnnotation_file,header =F, as.is =T)
tab_annot$V1 = substring(tab_annot$V1, 4)
domains = unique (tab_annot$V7)

domainLabel = lapply(domains, function(domain){
    ifOverlap = checkDomainOverlap(tab_ctcfCore, tab_annot, domain)
    res = data.frame( ifOverlap)
    colnames(res)[1] = domain
    res
})

domainLabel = do.call("cbind", domainLabel)

#---------------------------------------------------
#   TAD
#---------------------------------------------------
distToTAD = getDistToTAD (tab_ctcfCore, TAD_file)
fullTab = data.frame(tab_ctcfCore,tab_deepBind, distToTAD,  ifConstitutive,  tab_jaspar, domainLabel, stringsAsFactors =T)
colnames(fullTab)[4] = "macsFoldEnrich_rep1"


# --------------------------------------
#  Upstream Motif
# --------------------------------------
print("Get Upstream Motif")
fullTab = getUpMotif(fullTab, ref=genomeRef )


# --------------------------------------
#  GC content
# --------------------------------------
fullTab$gc = getGC(fullTab, ref=genomeRef )
fullTab$ifHasCPG = hasCpG(fullTab, ref = genomeRef) > 0


# ------------------------------------------------
#  predict persistence binding  
# ------------------------------------------------
# load(ctcfInsite_file)
# probPersistent_EnrichScaled <- predict(models[["scale"]],newdata=fullTab,type='prob')[,2]
# probPersistent <- predict(models[["noScale"]],newdata=fullTab,type='prob')[,2]
# fullTab = data.frame(fullTab, probPersistent_EnrichScaled, probPersistent,stringsAsFactors=F)
write.table(file = sprintf("%s/%s.fullAnnot.txt", resultDir, basename(macsF)), fullTab, quote=F, row.names =F)



