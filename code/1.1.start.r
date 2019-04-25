#### Libraries ####
library(taxonomizr)### convert accession id to tax id and to names ###
dir.create('input/database/ac2tax')
prepareDatabase('input/database/ac2tax/accessionTaxa.sql',types = c("nucl_gb", "nucl_wgs","prot"))
#file.remove(list.files('./input/database/ac2tax','accession2taxid.gz$'))
#file.remove(list.files('./input/database/ac2tax','dmp$'))
#### End ###
