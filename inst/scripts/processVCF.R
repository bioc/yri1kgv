Require = function(...) suppressPackageStartupMessages(require(..., quietly=TRUE))

dochrGr = function(chrn, ids, which, genome="hg19") {
Require(VariantAnnotation)
Require(snpStats)
 allvc = dir("/proj/rerefs/reref00/1000Genomes_Phase1_v3/ALL/", patt="ALL.*gz$",
  full=TRUE)
 f_ind = grep(paste0(chrn, "\\."), allvc)
 cpath = allvc[f_ind]
 stopifnot(file.exists(cpath))
 stopifnot(length(cpath)==1)
 vp = ScanVcfParam( info=NA, geno="GT", fixed="ALT", samples=ids,
        which=which )
 tmp = readVcf( cpath, genome=genome, param=vp )
 sz = prod(dim(tmp))
 if (sz == 0) return(NULL)
 genotypeToSnpMatrix(tmp)$genotypes
}

fixdupcol = function(x) {
 allc = unlist(nl <- lapply(x, colnames))
 bad = which(duplicated(allc))  # assumed very rare
 if (length(bad) == 0) return(x)
 message("NOTE: duplicated SNP id detected.")
 badids = allc[bad]
 lkc = sapply(x, is.null)
 if (any(lkc)) x = x[-which(lkc)]
 colcounts = sapply(x, dim)[2,]  # will fail if NULL passed
 ends = cumsum(colcounts)
 sts = c(1, ends[-length(ends)]+1)
 for (i in 1:length(ends)) {
   for (j in 1:length(bad)) {
     if (bad[j] >= sts[i] & bad[j] <= ends[i]) {
          kill = which(colnames(x[[i]]) == badids[j])
          x[[i]] = x[[i]][,-kill]
          }
     }
  }
 x
}
  
# badids = allc[bad]
# hasdup = sapply
 

dochr = function(chrn, ids, genome="hg19", size=5e6) {
 Require(piFDR2)
Require(GenomicRanges)
 Require(foreach)
 Require(doParallel)
 tt = tiling(chrn, size)
 seqlevels(tt) = gsub("chr", "", seqlevels(tt))
 ans <- foreach( x=1:length(tt) ) %dopar%  dochrGr( chrn, ids, tt[x], genome=genome ) 
 ans = fixdupcol(ans)
 ans = do.call(cbind, ans) # may lose class
 if (!is(ans, "SnpMatrix")) return(new("SnpMatrix", ans))
 ans
}

savechr = function(chrn, ids, genome="hg19", size=5e6) {
 Require(chan1kg)
 Require(GenomicRanges)
 obn = chrn
 assign(obn, dochr(chrn=chrn, ids=ids, genome=genome, size=size))
 save(list=obn, file=paste0(obn, ".rda"))
}

library(parallel)
library(foreach)
library(doParallel)
registerDoParallel(cores=2)
#registerDoSEQ()
#load("okids.rda")
library(yri1kgv)
data(eset)
sn = colnames(exprs(ex))
chrs = paste("chr", 4:1, sep="")
#debug(dochrGr)
lapply(chrs, function(x) {gc(); savechr(x, sn)})
 
 
#
#library(parallel)
#options(mc.cores=10)
#date()
#til = tiling(chrn, 1000000)
#GenomicRanges::seqlevels(til) = gsub("chr", "", chrn) 
#inds = 1:length(til)
#gr = floor(inds/10)
#sinds = split(inds, gr)
#nspl = length(sinds)
#obs = paste0("allm_", chrn, "_", 1:nspl)
#for (i in 1:nspl) {
#   tmp = mclapply(sinds[[i]], function(x) vcf2sm.gr1(cpath, til[x], idsToKeep=ids))
##   tmp = lapply(sinds[[i]], function(x) vcf2sm.gr1(cpath, til[x], idsToKeep=ids))
#   assign(obs[i], chk <- try(combSM(tmp)))
#   if (!inherits(chk, "try-error"))
#   save(list=obs[i], file=paste0(obs[i], ".rda"))
#   }
#}
