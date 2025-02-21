library(PheWAS)

genes = c("YLPM1", "GIGYF1")

for (gene in genes) {
  PROJECT_DIR <- "/Users/deeprobanerjee/Documents/bmi_project/BMI_monogenic"
  filename <- paste0(PROJECT_DIR, "/data/phewas/", gene, "/phewas_meta.csv.gz")
  phe_df <- read.csv(filename, colClasses = c("phenotype" = "character"))
  phe_df["p"] <- phe_df["p_value"]
  
  # Save the plot to a PDF
  save_file <- paste0(PROJECT_DIR, "/data/phewas/", gene, "/phewas_manhattan.pdf")
  pdf(save_file, width = 6.5, height = 3.5)  # Adjust size as needed
  print(phewasManhattan(
    phe_df, annotate.phenotype.description = FALSE, 
    point.size=1, title=paste0("Phewas: ", gene), size.x.labels=6, size.y.labels=6))
  dev.off()  # Close the PDF device
}

