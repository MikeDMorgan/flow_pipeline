estimate_heritability <- function(MZ_data, DZ_data){
  # subset MZ and DZ datasets by markers, calculate correlations (plot?)
  DZ.list_of_subset = list()
  MZ.list_of_subset = list()
  h2_list = list()

  for(i in seq_len(length(levels(DZ_data$marker)))){
    mark = levels(DZ.data$marker)[i]
    DZ.subset = subset(DZ.data, marker == mark)
    DZ.list_of_subset[[mark]] = DZ.subset
  
    MZ.subset = subset(MZ_data, marker == mark)
    MZ.list_of_subset[[mark]] = MZ.subset
  
    # standardize values to zero mean and unit variance, within twin class
    DZtwin1_z = (DZ.subset$twin1 - mean(DZ.subset$twin1))/sd(DZ.subset$twin1)
    DZtwin2_z = (DZ.subset$twin2 - mean(DZ.subset$twin2))/sd(DZ.subset$twin2)
    
    DZ.cor_z = cor(DZtwin1_z, DZtwin2_z)
    
    MZtwin1_z = (MZ.subset$twin1 - mean(MZ.subset$twin1))/sd(MZ.subset$twin1)
    MZtwin2_z = (MZ.subset$twin2 - mean(MZ.subset$twin2))/sd(MZ.subset$twin2)
    
    #MZ.cor_z = cor(MZtwin1_z, MZtwin2_z)
    DZ.cor = cor(DZ.subset$twin1, DZ.subset$twin2)
    MZ.cor = cor(MZ.subset$twin1, MZ.subset$twin2)

    h2 = 2 * (MZ.cor - DZ.cor)
    #h2_z = 2 * (MZ.cor_z - DZ.cor_z)
    h2_list[[mark]] = h2
    #h2_list[[mark]] = h2_z
  }
  return(h2_list)
}

MZ.data <- read.table("/ifs/projects/proj052/flow_processing_tables/MZ_twins_data.tsv", h=T, sep="\t")
DZ.data <- read.table("/ifs/projects/proj052/flow_processing_tables/DZ_twins_data.tsv", h=T, sep="\t")
h2_list <- estimate_heritability(MZ.data, DZ.data)
print(unlist(h2_list))