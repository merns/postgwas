# % File man/getAbstracts.Rd
# \name{getAbstracts}
# \alias{getAbstracts}
# \title{Retrive publication abstracts naming trait keywords and gene symbols associated to SNPs}
# 
# \description{
#   Writes files of PubMed abstracts relating each SNP with the given phenotype terms.
# }
# 
# \usage{
#   getAbstracts(
#     snps, 
#     pheno.terms = c("heparin", "thrombocytopenia"), 
#     genes.tag = "[TI]", 
#     pheno.tag = "[TIAB]"
#     )
# }
# 
# \arguments{
#   \item{snps}{A data frame containing columns 'SNP' and 'genename'. Each row forms a SNP to gene association that is considered for abstract retrieval. An appropriate dataframe is returned by the \code{\link{snp2gene}} function. 
#   }
# }
# 
# \details{
#   Literature search tools need to expand gene names with their synonyms and previous symbols so that all naming conventions of genes that might be used in all kinds of articles are covered. 
#   with rs-IDs \code{myconfig$gene$attr$name <- 'uniprot_swissprot_accession'}
#   
#   Filenames are [snps$SNP]_abstracts.txt
#   The usage of wildcards in pheno.terms is not recommended. 
#   Mesh term expansion is disabled for terms without wildcards (exact match).
# }
# 
# \value{
#   A
# }
# 
# \seealso{
#   \code{\link{snp2gene}}
# }
# 
# \examples{
#   
#   \dontshow{
#     ## preload data for offline usage
#     postgwas.buffer.snps <- read.table(
#       system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE
#       )
#     postgwas.buffer.genes <- read.table(
#       system.file("extdata", "postgwas.buffer.genes", package = "postgwas"), header = TRUE
#       )
#     
#     snps <- data.frame(SNP = c("rs172154", "rs759704"))
#     snps.mapped <- snp2gene.prox(snps, use.buffer = TRUE)
#     
#     postgwas:::getAbstracts(
#       snps.mapped = snps, 
#       pheno.terms = c("cancer", "tumor", "tumorigenesis")
#       )
#     
#   }
#   
# }




# getAbstracts <- function(
#                   snps.mapped, 
#                   pheno.terms, 
#                   genes.tag = "[TIAB]", 
#                   pheno.tag = "[TI]", 
#                   max.abstracts.per.snp = 100,
#                   hgnc.expansion = F
#                 ) {
#   #
#   # Value:
#   #  Beside writing files with SNP and abstract information, returns a list of PubMed query strings for each SNP which can be submitted in Entrez.
# 
#   # expand with hgnc symbols for human data
#   if(hgnc.expansion) {
# #   if(is.null(hgnc.file)) {
# #     # see http://www.genenames.org/cgi-bin/hgnc_downloads.cgi for HGNC download (and DB column names) documentation
# #     # we download only approved entries with a mapped ensembl id present
# #     hgnc.url <- url("http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?title=HGNC+output+data&submit=submit&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&col=md_ensembl_id&status=Approved&status_opt=2&level=pri&=on&where=md_ensembl_id+!%3D+%22%22&order_by=gd_app_sym_sort&limit=&format=text&.cgifields=&.cgifields=level&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag")
# #     hgnc.dl <- readLines(hgnc.url, warn = FALSE)
# #     close(hgnc.url)
# #     hgnc.list <- strsplit(hgnc.dl, "\t", fixed = TRUE)
# #     hgnc.df <- as.data.frame(list2df(hgnc.list))
# #     colnames(hgnc.df) <- hgnc.list[[1]]
# #     hgnc.df <- hgnc.df[-1, ]
# #   } else {
# #     hgnc <- read.table(hgnc_file, sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
# #     hgnc <- hgnc[hgnc$Status == "Approved", ]
# #     hgnc <- hgnc[, c("Ensembl ID (mapped data supplied by Ensembl)", "Approved Symbol", "Approved Name", "Previous Symbols", "Aliases")]
# #   }
#     
#     syms.snp <- merge(snps, hgnc, by.x = "SNP", by.y = "")
#     
#         syms.snp.string <- paste(
#                               syms.snp$"Approved Symbol", 
#                               syms.snp$"Previous Symbols", 
#                               syms.snp$Aliases, 
# #                              syms.snp$"Approved Name", 
#                               sep = ",", 
#                               collapse = ","
#                             )
#         syms.snp.vec <- strsplit(syms.snp.string, ",")[[1]]
#         syms.snp.vec <- syms.snp.vec[syms.snp.vec != ""]
# 
#   }
# 
# 
#   return(
#     sapply(
#       levels(factor(as.vector(snps.mapped$SNP))), 
#       function(snp) {
#         # all gene symbols for this snp
#         syms.snp <- snps.mapped[snps.mapped$SNP == snp, "genename"]
#         # concat gene symbols to string separated by 'genes.tag' and OR
#         syms.snp.string <- paste(
#                               syms.snp$"Approved Symbol", 
#                               syms.snp$"Previous Symbols", 
#                               syms.snp$Aliases, 
# #                              syms.snp$"Approved Name", 
#                               sep = ",", 
#                               collapse = ","
#                             )
#         syms.snp.vec <- strsplit(syms.snp.string, ",")[[1]]
#         syms.snp.vec <- syms.snp.vec[syms.snp.vec != ""]
#         syms.snp.vec <- gsub("[[:space:]]", "", syms.snp.vec)
#         syms <- paste(paste(syms.snp.vec, collapse = paste(genes.tag, "OR ")), genes.tag)
#         # concat with pheno.terms to pubmed query
#         query <- paste("(", syms, ") AND (", paste(paste(pheno.terms, collapse = paste(pheno.tag, "OR ")), pheno.tag), ")" )
#         # retrieve abstracts and write to file
#         pubmed.url <- url(gsub(" ", "+", paste("http://www.ncbi.nlm.nih.gov/pubmed?term=", query, "&dispmax=", max.abstracts.per.snp,"&report=medline&format=text", sep = "")))
#         pubmed.dl <- readLines(pubmed.url, warn = FALSE)
#         close(pubmed.url)
#         # reformat so that titles and abstracts span only one line
#         for(i in length(pubmed.dl):2) {
#           # when current line does not start with a medline tag, combine with the one before
#           if(!grepl("^[[:upper:]]+\\s*-", pubmed.dl[i])) {
#             pubmed.dl[i-1] <- gsub("\\s+", " ", paste(pubmed.dl[i-1], pubmed.dl[i]))
#             pubmed.dl <- pubmed.dl[-i]
#           }
#         }
#         # select only TI, TIAB, PMID lines and replace line breaks by html <p> tag
#         pubmed.res <- paste(grep("(PMID\\s*-)|(TI\\s*-)|(AB\\s*-)", pubmed.dl, value = TRUE), collapse = "<p>")
#         # highlight pheno terms and gene symbols with html tags
#         for(sym in syms.snp.vec)
#           pubmed.res <- gsub(
#                           gsub("[[:punct:]]", "", sym), 
#                           paste("<font color=\"red\">", sym, "</font>"), 
#                           pubmed.res, 
#                           ignore.case = TRUE
#                         )
#         for(term in pheno.terms)
#           pubmed.res <- gsub(
#                           term, 
#                           paste("<font color=\"blue\">", gsub("[[:punct:]]", "", term), "</font>"), 
#                           pubmed.res, 
#                           ignore.case = TRUE
#                         )
#         outfile <- file(paste(snp, "_abstracts.html", sep = ""), open = "w+")
#         writeLines(paste("PubMed query:", query, "<p>"), con = outfile)
#         writeLines(pubmed.res, con = outfile)
#         close(outfile)
#         return(query)
#       }, 
#       USE.NAMES = TRUE,
#       simplify = FALSE
#     )
#   )
# }
