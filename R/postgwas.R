postgwas <- function(
              gwas.dataset, 
              suggestive.p = 1e-06, 
              genomewide.p = 5e-08, 
              gts.source = NULL, 
              biomart.config = biomartConfigs$hsapiens, 
			        GOpackagename = "org.Hs.eg.db"
            ) {

  clearPostgwasBuffer()
  
  if(!is.null(GOpackagename)) {
    if(!(GOpackagename %in% installed.packages()[, "Package"])) {
     message(paste("Annotation Package ", GOpackagename, " from bioconductor is not installed. Network will not be colorized by overrepresented GO terms."))
      vertexcolor.GO.overrep <- NULL
    } else {
      vertexcolor.GO.overrep <- GOpackagename
    }
  }
  
  
  message("\n\n****************** PREPARING SNPS ******************\n")
  
  snps <- readGWASdatasets(gwas.dataset, assert.unique.snps.per.dataset = TRUE)
  
  # check for rs-IDs
  if(sum(grepl("rs\\d+", snps$SNP)) < nrow(snps) / 2)
    warning("Most SNPs are not named as dbSNP rs-ids. Make sure that the supplied IDs match the biomart config.\n")
  
  
  message("\n\n****************** SNP TO GENE MAPPING ******************\n")
  
  message("\nAnnotating genes (proximity)...")
  tryCatch(
    {
      s2g <- snp2gene.prox(
        snps[snps$P <= suggestive.p, ], 
        level = 0, 
        biomart.config = biomart.config, 
        use.buffer = TRUE
        );
      s2g <- s2g[, c("SNP", "CHR", "BP", "P", "geneid", "genename", "start", "end", "direction")]
    },
    error = function(e) message(paste("snp2gene.prox: An error occured, proceeding with the next function.  ", e))
  )
  
  tryCatch(
    if(!is.null(gts.source)) {
      message("\nAnnotating genes (LD)...")
      s2g.LD <- snp2gene.LD(
        snps[snps$P <= suggestive.p, ], 
        gts.source = gts.source, 
        biomart.config = biomart.config, 
        use.buffer = TRUE
      )
      s2g.LD <- s2g.LD[s2g.LD$ld.max > 0.3, c("SNP", "CHR", "BP", "P", "geneid", "genename", "start", "end", "ld.max", "ld.mean", "ld.sdev")]
      s2g.LD$direction <- "LD"
      if(exists("s2g")) {
        s2g$ld.max <- NA
        s2g$ld.mean <- NA
        s2g$ld.sdev <- NA
        s2g <- rbind(s2g, s2g.LD)
      } else {
        s2g <- s2g.LD
      }
    },
    error = function(e) message(paste("snp2gene.LD: An error occured, proceeding with the next function.  ", e))
  )
  
  if(exists("s2g")) 
    write.table(s2g, paste("postgwas_snp2gene_p", suggestive.p, ".csv", sep = ""), row.names = FALSE, sep = "\t")
  
  
  
  message("\n\n****************** MANHATTAN PLOT ******************\n")
  
  tryCatch(
    manhattanplot(
      gwas.dataset = gwas.dataset, 
      highlight.logp = -log10(genomewide.p),
      reduce.dataset = 5, 
      biomart.config = biomart.config
    ), 
    error = function(e) message(paste("manhattanplot: An error occured, proceeding with the next function.  ", e))
  )
  
  message("\n\n****************** REGIONAL PLOTS ******************\n")
  
  tryCatch(
    {
      message("Select representative SNPs (lowest P) within 50 kb windows for regionalplots...");
      snps <- removeNeighborSnps(snps[snps$P <= suggestive.p, c("SNP", "P")], maxdist = 200000, biomart.config = biomart.config, use.buffer = TRUE);
      if(nrow(snps) <= 0)
        stop("The dataset does not contain SNPs with suggestive significance.\n");
      
      if(!is.null(gts.source)) {
        ld.options <- list(
            gts.source = gts.source, 
            max.snps.per.window = ceiling(300 / log10(nrow(snps) +1)), 
            rsquare.min = 0.2, 
            show.rsquare.text = FALSE
        )
      } else  {
        ld.options <- NULL
      };
      
      setPostgwasBuffer(snps = NULL);
      message(paste("Regionalplots of", nrow(snps), "loci..."));
      regionalplot(
        snps = snps, 
        gwas.datasets = gwas.dataset, 
        ld.options = ld.options, 
        biomart.config = biomart.config, 
        max.logp = -log10(genomewide.p) * 2.5, 
        use.buffer = TRUE
      )
    },
    error = function(e) message(paste("regionalplot: An error occured, proceeding with the next function.  ", e))
  )
  
  
  if(exists("s2g")) {
  
    s2g <- s2g[, c("SNP", "geneid", "genename", "P")]
    s2g <- s2g[order(s2g$P), ]
    s2g <- s2g[!duplicated(s2g$genename), ] # is sorted: keeps better p-values
    s2g <- vectorElements(s2g)
    
    if(nrow(s2g) <= 0)
      stop("No SNPs were mapped to genes.\n")
    
    message("\n\n****************** NETWORK ANALYSIS ******************\n")
    message("\nDownloading network data...")
    tryCatch(
      gwas2network(
        gwas.mapped.genes = s2g, 
        network = getInteractions.path(filter.ids = unique(s2g$geneid), toFile = "postgwas.interaction.download.path"),  
        prune = "gwasonly", 
        vertexcolor.GO.overrep = vertexcolor.GO.overrep, 
        biomart.config = biomart.config, 
        use.buffer = TRUE
      ),
      error = function(e) message(paste("gwas2network: An error occured.  ", e))
    )
  } 
  
  message("\n\n****************** THANKS FOR USING POSTGWAS ******************\n")
  message("Your output files have been produced.")
  message("Try also out the individual postgwas functions - ")
  message("they offer additional features and customization possibilities.\n")
  
}
