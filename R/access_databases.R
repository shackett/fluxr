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
    dplyr::mutate(systematic_gene = sub('^[a-z]+:', '', kegg_gene_id))

  if(!is.null(genes)){
    gene_classes <- gene_classes %>%
      dplyr::filter(gene %in% genes)
  }

  gene_classes
}


#' Extended Gene Annotation
#'
#' @param dataset a biomaRt dataset contained within listDatasets(useMart('ensembl'))
#'
#' @return a tibble combining multiple IDs for a gene
#' @export
#'
#' @examples
#' extended_gene_annotation("scerevisiae_gene_ensembl")
extended_gene_annotation <- function(dataset){

  ensembl = biomaRt::useMart("ensembl", dataset = dataset)

  extended_annotations <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description'), mart = ensembl)
  colnames(extended_annotations) <- c("systematic_gene", "common_gene", "gene_description")

  extended_annotations
}


#'
download_BRENDA_regulators <- function(){
# all EC numbers  "http://www.brenda-enzymes.org/all_enzymes.php"
}

