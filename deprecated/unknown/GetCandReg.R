##############################################################################
GetCandReg <- function(highobj, annot, traits)
{
  ## currently this only gets max over genome; want max by chr, yes?
  traits <- unique(traits)
  n <- length(traits)
  out <- data.frame(matrix(NA, n, 6))
  names(out) <- c("gene", "phys.chr", "phys.pos", "peak.chr", "peak.pos",
                  "peak.lod")
  
  ## Get annotation of trait name and physical chromosome and position (in cM).
  m <- match(traits, annot[,1])
  out[, 1:3] <- annot[m, c(1,3,5)]
  
  ## Get LOD peak information
  m <- !is.na(out[,3])
  
  pheno.cols <- match(traits[m], highobj$names)
  if(any(is.na(pheno.cols)))
    stop("some traits do not have scans")
  
  chr.pos <- highobj$chr.pos
  highlod <- highobj$highlod[highobj$highlod$phenos %in% pheno.cols,]
  tmp <- cumsum(table(highlod[,"phenos"]))
  tmp <- c(0, tmp[-length(tmp)])
  peak.index <- tmp + tapply(highlod$lod, highlod$phenos, which.max)
  
  ## now relate to lod, chr, pos, but get order right with traits
  m <- match(highobj$names[highlod[peak.index, "phenos"]], as.character(out[,1]))
  if(any(is.na(m)))
    stop("cannot match highlod with pheno names")
  
  out[m, 6] <- highlod[peak.index, "lod"]
  out[m, 4] <- chr.pos[highlod[peak.index, "row"], "chr"]
  out[m, 5] <- chr.pos[highlod[peak.index, "row"], "pos"]
  
  out[!is.na(out[,4]),, drop = FALSE]
}
