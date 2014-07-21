### R code from vignette source 'postgwas.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: postgwas.Rnw:22-23
###################################################
library(xtable)


###################################################
### code chunk number 2: init
###################################################
library(postgwas)


###################################################
### code chunk number 3: postgwas.Rnw:90-92
###################################################
load("bufferHS.RData")
setPostgwasBuffer(bufferHS)


###################################################
### code chunk number 4: postgwas.Rnw:110-116
###################################################
manhattanplot(
    gwas.dataset = "whrTrunc.txt.remapped.xz",
    highlight.text = NULL,
    use.buffer = TRUE,
    toFile = NULL
)


###################################################
### code chunk number 5: postgwas.Rnw:148-159
###################################################
jpeg(filename = "manh.jpeg", width = 1600, height = 400, pointsize = 25, quality = 95)
manhattanplot( "whrTrunc.txt.remapped.xz",
    highlight.logp = 7.3,
    highlight.text = "SNP",
    highlight.cex = 0.8,
    highlight.col = "cyan",
    highlight.fontface = "plain",
    use.buffer = TRUE,
    toFile = NULL
)
dev.off()


###################################################
### code chunk number 6: manh3
###################################################
manhattanplot(
    "whrTrunc.txt.remapped.xz",
    highlight.logp = c(6, 8, 10),
    highlight.fontface = c("italic", "italic", "bold"),
    highlight.cex = c(0.6, 0.6, 1),
    highlight.text = c("SNP", "genes", "genes"),
    highlight.win = c(75000, 200000, 200000),
    ticks.y = TRUE,
    plot.title = "WHR dataset plus gene annotation",
    use.buffer = TRUE,
    reduce.dataset = 5,
    toFile = NULL
)


###################################################
### code chunk number 7: reg1
###################################################
snps <- data.frame(SNP = c("rs2745353", "rs10923712", "rs4846567"))
regionalplot(
    snps = snps,
    gwas.datasets = c("heightTrunc.txt.remapped.xz", "whrTrunc.txt.remapped.xz"),
    max.logp = 15,
    ld.options = NULL,
    out.format = NULL,
    use.buffer = TRUE
)


###################################################
### code chunk number 8: postgwas.Rnw:292-294
###################################################
load("LDrs4846567region.RData")
setPostgwasBuffer(ld.regionalplot = LDrs4846567)


###################################################
### code chunk number 9: reg2
###################################################
regionalplot(
  snps = data.frame(SNP = "rs4846567"),
  gwas.datasets = "whrTrunc.txt.remapped.xz",
  ld.options = list(gts.source = 2),
  out.format = list(file = "pdf", panels.per.page = 3),
  max.logp = 15,
  use.buffer = TRUE
)


###################################################
### code chunk number 10: postgwas.Rnw:339-343
###################################################
if(library("Rsamtools", logical = TRUE)) {
  bgzip("reseq.vcf", overwrite = TRUE)
  indexTabix("reseq.vcf.bgz", format = "vcf4")
}


###################################################
### code chunk number 11: postgwas.Rnw:392-394
###################################################
chrmap <- data.frame(CHR = 1:23, CHR.VCF = paste("chr", 1:23, sep = ""))
head(chrmap)


###################################################
### code chunk number 12: postgwas.Rnw:429-431
###################################################
load("LDrs10923712region.RData")
setPostgwasBuffer(ld.regionalplot = LDrs10923712)


###################################################
### code chunk number 13: reg3
###################################################
regionalplot(
  snps = data.frame(SNP = "rs10923712"),
  gwas.datasets = "whrTrunc.txt.remapped.xz",
  window.size = 600000,
  ld.options = list(gts.source = 2, max.snps.per.window = 100),
  out.format = list(file = "pdf", panels.per.page = 3),
  max.logp = 15,
  var.options = list(
      vcf = list(
          file = "reseq.vcf.bgz",
          remap.positions = FALSE,
          chrom.map = chrmap
      ),
      vcf.info.colorize = c(EFF = "SYNONYMOUS", EFF = "NON_SYNONYMOUS")
  ),
  use.buffer = TRUE
)


###################################################
### code chunk number 14: s2g
###################################################
snps <- data.frame(SNP = c("rs10923712", "rs4846567"))
s2g <- as.matrix(snp2gene.prox(snps, use.buffer = TRUE))
xtable(
    s2g,
    caption = "SNP to gene annotation by proximity. The closest gene
in each direction (covering, up- and downstream genes, including overlapping genes)
is given for each query SNP.",
    label = "tab:s2g",
    align = rep(">{\\scriptsize\\sffamily}c", ncol(s2g) +1)
)


###################################################
### code chunk number 15: s2gLD
###################################################
s2g <- snp2gene.LD(snps, use.buffer = TRUE, gts.source = "s2g.ped.xz")
xtable(
    as.matrix(s2g[!is.na(s2g$geneid), ]),
    caption = "SNP to gene annotation by LD (showing only genes with entrez ID).
The column ld.max lists the highest $r^2$ value between the query SNP and all SNPs within the gene in question,
ld.mean the mean of all $r^2$ values of the query SNP with SNPs in that gene and ld.sdev its standard deviation, respectively.",
    label = "tab:s2gLD",
    align = rep(">{\\scriptsize\\sffamily}c", ncol(s2g) +1)
)


###################################################
### code chunk number 16: whr22Prep
###################################################
whr <- read.table("whrTrunc.txt.remapped.xz", header = T)
whr22 <- whr[whr$CHR == "22", ]


###################################################
### code chunk number 17: whr22Prep2
###################################################
whr22 <- snp2gene.prox(whr22, level = 0, use.buffer = TRUE)
whr22 <- whr22[whr22$direction == "cover", ]


###################################################
### code chunk number 18: gates
###################################################
whr22.gp <- gene2p(
    whr22,
    method = GATES,
    gts.source = "whr22.ped.xz"
)


###################################################
### code chunk number 19: postgwas.Rnw:607-617
###################################################
whr22.gp <- whr22.gp[, !grepl("original", colnames(whr22.gp))] # remove old chr and bp cols
whr22.gp$gene.p <- signif(as.numeric(as.vector(whr22.gp$gene.p)), 3)
xtable(
    as.matrix(head(whr22.gp)),
    caption = "SNP to gene annotation including a column \\emph{gene.p} for gene-wise aggregated p-values.
Only the first six lines of the original result table are shown. Calculations were based on intragenic SNPs on
chromosome 22 of the \\emph{WHR} dataset.",
    label = "tab:g2p",
    align = rep(">{\\scriptsize\\sffamily}c", ncol(whr22.gp) +1)
)


###################################################
### code chunk number 20: SpD
###################################################
whr22.gp.SpD <- gene2p(
    whr22,
    method = SpD,
    gts.source = "whr22.ped.xz"
)
whr22.gp.SpD <- whr22.gp.SpD[order(whr22.gp.SpD$SNP), ]
whr22.gp  <- whr22.gp[order(whr22.gp.SpD$SNP), ]


###################################################
### code chunk number 21: postgwas.Rnw:636-637
###################################################
mean(log10(whr22.gp$gene.p) - log10(whr22.gp.SpD$gene.p), na.rm = TRUE)


###################################################
### code chunk number 22: GOenrich
###################################################
enrich.res <- gwasGOenrich(gwas = whr22.gp, ontology = "CC", pruneTermsBySize = 8)


###################################################
### code chunk number 23: postgwas.Rnw:665-671
###################################################
xtable(
    enrich.res,
    caption = "Result of applying the gene set enrichment analysis functions of the topGO package to the \\emph{WHR} chromosome 22 dataset , using the \\emph{gwasGOenrich} function. ",
    label = "tab:GOenrich",
    align = rep(">{\\scriptsize\\sffamily}c", ncol(enrich.res) +1)
)


###################################################
### code chunk number 24: net1prep
###################################################
snpsW <- removeNeighborSnps(whr[whr$P < 10^-6, ])
genesW.closest <- snp2gene.prox(snpsW, level = 0, use.buffer = TRUE)


###################################################
### code chunk number 25: postgwas.Rnw:756-757
###################################################
network.data <- network.data.net1


###################################################
### code chunk number 26: net1
###################################################
network.igraph <- gwas2network(
  gwas.mapped.genes = genesW.closest,
  network = network.data,
  max.communities = 0,
  vertexcolor.GO.regex = list(red = "metaboli"),
  min.transparency = 0.4,
  max.transparency = 1,
  use.buffer = TRUE
)


###################################################
### code chunk number 27: postgwas.Rnw:805-810
###################################################
gwas2network.plot(
    network.igraph,
    filename = "exampleNet.jpeg",
    device = jpeg, res = 100, pointsize = 16
)


###################################################
### code chunk number 28: net2prep
###################################################
genesW.LD <- snp2gene.LD(snpsW, use.buffer = TRUE, gts.source = "whr.ped.xz")
genesW <- rbind(
  genesW.closest[, c("SNP", "P", "geneid")],
  genesW.LD[genesW.LD$ld.max > 0.6 | genesW.LD$ld.mean > 0.1, c("SNP", "P", "geneid")]
)


###################################################
### code chunk number 29: postgwas.Rnw:868-869
###################################################
network.data <- network.data.net2


###################################################
### code chunk number 30: net2
###################################################
network.igraph <- gwas2network(
  gwas.mapped.genes = genesW,
  network = network.data,
  max.communities = 0,
  vertexcolor.GO.regex = list(
                           red = "metaboli",
                           yellow = "(LDL|HDL|lipoprotein)",
                           blue = "(insulin|diabet|glucose)"
                         ),
  min.transparency = 0.4,
  max.transparency = 1,
  use.buffer = TRUE
)


###################################################
### code chunk number 31: postgwas.Rnw:888-889
###################################################
gwas2network.plot(network.igraph, file = "exampleNet2.pdf")


###################################################
### code chunk number 32: GenABEL
###################################################
data(srdta)
gwas <- ccfast("bt", srdta)


###################################################
### code chunk number 33: GenABEL2
###################################################
regionalplot(
  snps = data.frame(SNP = "rs1020"),
  gwas.datasets = gwas,
  ld.options = list(gts.source = srdta@gtdata),
  plot.genes = FALSE
)


###################################################
### code chunk number 34: GenABEL3
###################################################
gwas.custom <- data.frame(
    SNP = snpnames(gwas),
    P = gwas[, "Pc1df"],
    BP = gwas[, "Position"],
    CHR = chromosome(gwas)
)


###################################################
### code chunk number 35: postgwas.Rnw:1142-1144
###################################################
load("bufferMM.RData")
setPostgwasBuffer(bufferMM)


###################################################
### code chunk number 36: postgwas.Rnw:1169-1176
###################################################
manhattanplot(
    "mouse_lmm.assoc.txt.xz",
    biomart.config = biomartConfigs$mmusculus,
    reduce.dataset = FALSE,
    use.buffer = TRUE,
    toFile = FALSE
)


###################################################
### code chunk number 37: postgwas.Rnw:1194-1196
###################################################
mm <- postgwas:::readGWASdatasets("mouse_lmm.assoc.txt.xz")
snps.mm <- removeNeighborSnps(mm[mm$P < 5*10^-8, ], maxdist = 100000)


###################################################
### code chunk number 38: s2gMouse
###################################################
genes.mm <- snp2gene.prox(
    snps.mm,
    level = 0,
    use.buffer = TRUE,
    biomart.config = biomartConfigs$mmusculus
)
genes.mm.colsel <- as.matrix(genes.mm[, c("geneid", "genename", "start", "end", "CHR", "BP", "SNP",  "P", "direction")])
xtable(
    genes.mm.colsel,
    caption = "SNP to gene annotation by proximity for mouse data",
    label = "tab:s2gMM",
    align = rep(">{\\scriptsize\\sffamily}c", ncol(genes.mm.colsel) +1)
)


###################################################
### code chunk number 39: regMouse
###################################################
regionalplot(
    snps = snps.mm,
    gwas.datasets = c("mouse_lmm.assoc.txt.xz"),
    window.size = 350000,
    biomart.config = biomartConfigs$mmusculus,
    ld.options = NULL,
    use.buffer = TRUE
)


###################################################
### code chunk number 40: postgwas.Rnw:1264-1280
###################################################
if(library(org.Mm.eg.db, logical.return = TRUE)) {
  network.data <- getInteractions.GO(
      genes.mm$geneid,
      GOpackagename = "org.Mm.eg.db",
      similarity = "hausdorff"
  )
  network.igraph <- gwas2network(
      genes.mm,
      network.data,
      max.communities = 0,
      vertexcolor.GO.overrep = "org.Mm.eg.db",
      biomart.config = biomartConfigs$mmusculus,
      use.buffer = TRUE
  )
  gwas2network.plot(network.igraph, file = "exampleNetMouse.pdf")
}


###################################################
### code chunk number 41: postgwas.Rnw:1313-1314
###################################################
buffers <- getPostgwasBuffer()


###################################################
### code chunk number 42: postgwas.Rnw:1320-1321
###################################################
names(buffers)


###################################################
### code chunk number 43: postgwas.Rnw:1327-1329
###################################################
snpbuffer <- buffers[["snps"]]
head(snpbuffer)


###################################################
### code chunk number 44: postgwas.Rnw:1336-1338
###################################################
snpbuffer[1, "refsnp_id"] <- "alternativeID"
setPostgwasBuffer(snps = snpbuffer)


###################################################
### code chunk number 45: postgwas.Rnw:1345-1346
###################################################
setPostgwasBuffer(uselist = buffers)


###################################################
### code chunk number 46: postgwas.Rnw:1421-1422
###################################################
head(read.table("whrTrunc.txt.remapped.xz", header = TRUE))


