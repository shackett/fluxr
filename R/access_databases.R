#' Associate genes with KEGG
#'
#' Find KEGG pathway IDs and names associated with an organisms' genes
#'
#' @param organism_code a three letter code
#' @param classification type of kegg classification to use: either module or pathway
#' @param optional character vector of gene names to return
#'
#' @return tibble of KEGG gene names and their corresponding pathway codes and names
#'
#' @export
#'
#' @examples
#' associated_genes_with_kegg("sce")
associated_genes_with_kegg <- function(organism_code, classification = "module", genes = NULL){
  stopifnot(is.character(genes) | is.null(genes))
  stopifnot(is.character(organism_code))
  stopifnot(length(organism_code) == 1)
  stopifnot(nchar(organism_code) == 3)

  if(classification == "module"){

  gene_class_ids <- suppressMessages(readr::read_tsv(paste0("http://rest.kegg.jp/link/", organism_code, "/module"), col_names = c("kegg_class_id", "kegg_gene_id")))
  gene_class_codes <- suppressMessages(readr::read_tsv(paste0("http://rest.kegg.jp/list/module/", organism_code), col_names = c("kegg_class_id", "kegg_class_name")))

  }else if(classification == "pathway"){

  gene_class_ids <- suppressMessages(readr::read_tsv(paste0("http://rest.kegg.jp/link/", organism_code, "/pathway"), col_names = c("kegg_class_id", "kegg_gene_id")))
  gene_class_codes <- suppressMessages(readr::read_tsv(paste0("http://rest.kegg.jp/list/pathway/", organism_code), col_names = c("kegg_class_id", "kegg_class_name")))

  }else{stop('"classification" type not supported')}

  gene_classes <- gene_class_ids %>% dplyr::inner_join(gene_class_codes, by = "kegg_class_id") %>%
    dplyr::mutate(gene = sub('^[a-z]+:', '', kegg_gene_id))

  if(!is.null(genes)){
    gene_classes <- gene_classes %>%
      dplyr::filter(gene %in% genes)
  }

  gene_classes
}
