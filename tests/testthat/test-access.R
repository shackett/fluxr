test_that("\code{associated_genes_with_kegg} finds yeast genes", {

  associated_df <- associated_genes_with_kegg(organism_code = "sce", classification = "module")

  expect_true("tbl_df" %in% class(associated_df))
  expect_gt(nrow(associated_df), 600)
  expect_equal(colnames(associated_df), c("kegg_class_id", "kegg_gene_id", "kegg_class_name", "gene"))

})
