#programa per fer correr LULU corregint el bug que Mahe ha trobat al fer el MUMU (veure https://github.com/tobiasgf/lulu/issues/8)
#el bug fa que nomes junti MOTUs quan la cooccurrencia es 1 encara que li hagis dit 0.95
#llegeix una matriu de OTUs, i calcula les distancies entre les sequencies


#per instal.lar LULU
#library(devtools)
#install_github("tobiasgf/lulu")

#library(lulu)
library(ape)
library(dplyr)

#funcio LULU modificada

lulu2<-function(otutable, matchlist, minimum_ratio_type = "min", 
          minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95) 
{
  require(dplyr)
  start.time <- Sys.time()
  colnames(matchlist) <- c("OTUid", "hit", "match")
  matchlist = matchlist[which(matchlist$hit != "*"), 
  ]
  matchlist = matchlist[which(matchlist$hit != matchlist$OTUid), 
  ]
  statistics_table <- otutable[, 0]
  statistics_table$total <- rowSums(otutable)
  statistics_table$spread <- rowSums(otutable > 0)
  statistics_table <- statistics_table[with(statistics_table, 
                                            order(spread, total, decreasing = TRUE)), ]
  otutable <- otutable[match(row.names(statistics_table), row.names(otutable)), 
  ]
  statistics_table$parent_id <- "NA"
 
 #linia corregida
 #log_con <- file(paste0("lulu.log_", format(start.time,"%Y%m%d_%H%M%S")), open = "a")

 log_con <- file(paste0("lulu_bug_corrected.log_", format(start.time,"%Y%m%d_%H%M%S")), open = "a")
  for (line in seq(1:nrow(statistics_table))) {
    print(paste0("progress: ", round(((line/nrow(statistics_table)) * 
                                        100), 0), "%"))
    potential_parent_id <- row.names(otutable)[line]
    cat(paste0("\n", "####processing: ", potential_parent_id, 
               " #####"), file = log_con)
    daughter_samples <- otutable[line, ]
    hits <- matchlist[which(matchlist$OTUid == potential_parent_id & 
                              matchlist$match > minimum_match), "hit"]
    cat(paste0("\n", "---hits: ", hits), file = log_con)
    last_relevant_entry <- sum(statistics_table$spread >= 
                                 statistics_table$spread[line])
    potential_parents <- which(row.names(otutable)[1:last_relevant_entry] %in% 
                                 hits)
    cat(paste0("\n", "---potential parent: ", 
               row.names(statistics_table)[potential_parents]), 
        file = log_con)
    success <- FALSE
    if (length(potential_parents) > 0) {
      for (line2 in potential_parents) {
        cat(paste0("\n", "------checking: ", 
                   row.names(statistics_table)[line2]), file = log_con)
        if (!success) {
          relative_cooccurence <- sum((daughter_samples[otutable[line2, 
          ] > 0]) > 0)/sum(daughter_samples > 0)
          cat(paste0("\n", "------relative cooccurence: ", 
                     relative_cooccurence), file = log_con)
          if (relative_cooccurence >= minimum_relative_cooccurence) {
            cat(paste0(" which is sufficient!"), 
                file = log_con)
            if (minimum_ratio_type == "avg") {
              # relative_abundance <- mean(otutable[line2, 
              # ][daughter_samples > 0]/daughter_samples[daughter_samples > 
              #                                            0])
              # cat(paste0("\n", "------mean avg abundance: ", 
              #            relative_abundance), file = log_con)
              #linies canviades a
              daughter_sampless<-daughter_samples[,which(otutable[line2,]>0)]
              relative_abundance <- mean(otutable[line2, which(otutable[line2,]>0)
              ][daughter_sampless > 0]/daughter_sampless[daughter_sampless > 
                                                         0])
              cat(paste0("\n", "------mean avg abundance: ", 
                         relative_abundance), file = log_con)
            }
            else {
              # relative_abundance <- min(otutable[line2, 
              # ][daughter_samples > 0]/daughter_samples[daughter_samples > 
              #                                            0])
              # cat(paste0("\n", "------min avg abundance: ", 
              #            relative_abundance), file = log_con)
              
              #linies canviades a
              daughter_sampless<-daughter_samples[,which(otutable[line2,]>0)]
              relative_abundance <- min(otutable[line2,which(otutable[line2,]>0)
              ][daughter_sampless > 0]/daughter_sampless[daughter_sampless >
                                                         0])
              cat(paste0("\n", "------min avg abundance: ",
                         relative_abundance), file = log_con)
              
            }
            if (relative_abundance > minimum_ratio) {
              cat(paste0(" which is OK!"), file = log_con)
              if (line2 < line) {
                statistics_table$parent_id[line] <- statistics_table[row.names(otutable)[line2], 
                                                                     "parent_id"]
                cat(paste0("\n", "SETTING ", 
                           potential_parent_id, " to be an ERROR of ", 
                           (statistics_table[row.names(otutable)[line2], 
                                             "parent_id"]), "\n"), 
                    file = log_con)
              }
              else {
                statistics_table$parent_id[line] <- row.names(otutable)[line2]
                cat(paste0("\n", "SETTING ", 
                           potential_parent_id, " to be an ERROR of ", 
                           (row.names(otutable)[line2]), "\n"), 
                    file = log_con)
              }
              success <- TRUE
            }
          }
        }
      }
    }
    if (!success) {
      statistics_table$parent_id[line] <- row.names(statistics_table)[line]
      cat(paste0("\n", "No parent found!", 
                 "\n"), file = log_con)
    }
  }
  close(log_con)
  total_abundances <- rowSums(otutable)
  curation_table <- cbind(nOTUid = statistics_table$parent_id, 
                          otutable)
  statistics_table$curated <- "merged"
  curate_index <- row.names(statistics_table) == statistics_table$parent_id
  statistics_table$curated[curate_index] <- "parent"
  statistics_table <- transform(statistics_table, rank = ave(total, 
                                                             FUN = function(x) rank(-x, ties.method = "first")))
  curation_table <- as.data.frame(curation_table %>% group_by(nOTUid) %>% 
                                    summarise_all(list(sum)))
  row.names(curation_table) <- as.character(curation_table$nOTUid)
  curation_table <- curation_table[, -1]
  curated_otus <- names(table(statistics_table$parent_id))
  curated_count <- length(curated_otus)
  discarded_otus <- setdiff(row.names(statistics_table), curated_otus)
  discarded_count <- length(discarded_otus)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  result <- list(curated_table = curation_table, curated_count = curated_count, 
                 curated_otus = curated_otus, discarded_count = discarded_count, 
                 discarded_otus = discarded_otus, runtime = time.taken, 
                 minimum_match = minimum_match, minimum_relative_cooccurence = minimum_relative_cooccurence, 
                 otu_map = statistics_table, original_table = otutable)
  return(result)
}

#main program





setwd("/piec2/j.zarcero/EXPERIMENTS/HARS/RBVC/MJOLNIR3/DEMULTIPLEXED_FILES/")
list.files()
data<-read.table("./RBVC_FRIGGA.tsv",sep="\t",head=T,stringsAsFactors=F)

otutab<-data[,16:266]
row.names(otutab)<-data$id

#generem fitxer seqs
rr<-(cbind(as.character(data$id),as.character(data$sequence)))
r<-strsplit(rr[,2],"")
names(r) <- rr[,1]
m<-as.DNAbin(r)
#calculem distancies entre seqs
d<-100-(dist.dna(m,model="N")*100/313)

#saveRDS(d,"./distance")
#d<-readRDS("distance")


nMOTUs<-length(m)
query<-vector("character",(nMOTUs*nMOTUs-nMOTUs)/2)
matching<-vector("character",(nMOTUs*nMOTUs-nMOTUs)/2)
index_query<-vector("integer",(nMOTUs*nMOTUs-nMOTUs)/2)
index_matching<-vector("integer",(nMOTUs*nMOTUs-nMOTUs)/2)
contador<-0

for (i in 1:(length(m)-1))
{
  if (i/1000-ceiling(i/1000)==0) message(Sys.time()," processant MOTU ",i," de ",length(m))
  for (j in (i+1):length(m))
    {
    contador<-contador+1
    index_query[contador]<-i
    index_matching[contador]<-j
    }
}


query<-names(r)[index_query]
matching<-names(r)[index_matching]

# aixo es massa lent
# query<-vector()
# matching<-vector()
# 
# 
# for (i in 1:(length(m)-1))
# {
#   if (i/100-ceiling(i/100)==0) message(Sys.time()," processant MOTU ",i," de ",length(m))
#   for (j in (i+1):length(m))
#   {
#     query<-append(query,names(r)[i])
#     matching<-append(matching,names(r)[j])
#   }
# }


#ULL que les columnes 1 i 2 de matchlist no poden ser factors
matchlist<-data.frame(cbind(query,matching,"similarity"=d),stringsAsFactors=F)



#saveRDS(matchlist,"./matchlist")
#matchlist<-readRDS(matchlist)

curated_result <- lulu2(otutab, matchlist)

saveRDS(curated_result,"./curated_result_bug_corrected")

message("he acabat")
