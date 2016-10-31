#####################################################
# Unserialize arrays and return them as a dataframe #
# with useful metadata added                        #
#####################################################
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(data.table))

unserializeArray <- function(rdb, statement, gate, panel, table){
    rconn <- dbConnect(SQLite(), rdb)
    r.df <- dbGetQuery(rconn, statement)

    # fields should include `array`
    r.df$real.array <- lapply(unlist(r.df$array),
                              FUN=function(x){unserialize(charToRaw(x))})
    batch_list <- add_batch(r.df)

    # drop the raw array
    out_df <- r.df[,c(1:4, 6, 7)]
   
   # use rbindlist as its quicker than rbind
   arrays <- data.frame(rbindlist(r.df$real.array))

   # add the metadata
   ids <- unlist(strsplit(r.df$rows, fixed=T, split="/"))
   arrays$twin.id <- ids
   arrays$batch <- unlist(batch_list)
   arrays$panel <- panel
   arrays$gate <- gate
   cols <- colnames(arrays)[!is.na(colnames(arrays))]
   all_cols <- cols[!grepl("NA", cols)]
   subarrays <- arrays[all_cols]

   return(subarrays)
}