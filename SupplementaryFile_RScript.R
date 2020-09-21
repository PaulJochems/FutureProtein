

rm(list = ls())  # remove any variables in R's memory 
#pas dit aan naar de output omgeving
setwd("C:/Users/Willem/Desktop")

library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(vegan);
library(MASS);
library(mclust);
library(tidyverse);
library(cluster);
library(factoextra);


#we create our own function called matching, we use this to match peptides which are similar and whitin the boundaries. 
# id multiple lines are matched, we calcualte the neirest using the euclidean distance. 
matching <- function(data_referention, data_proteins, mz_threshold, rt_threshold) {
  
  
  #we use brute force to calculate the matches of the protein to the referentions. IT is possible that one protein matches to multiple referentions. 
  #in thath case the protein is sotred multiple times. After the matching is done we find the neares point (if multiple matches) by calculating the
  #distance in the 2D space, using pythagoras.
  
  
  #we create an oversized matrix to store the results, aka the matched lines
  m.result = matrix(nrow = nrow(data_proteins) * 100, 
                    ncol = 8,
                    dimnames =  list(paste("ID", 1:(nrow(data_proteins)*100)), 
                                     list("ID", "RT", "mz", "RT_ref", "mz_ref", "RT_delta", "mz_delta", "match_ID")))
  
  #we want to count the lines we add to the matrix 
  line_index = 1
  
  #now loop over every line of the proteins data, and try every line of the referention dat to match it to.
  #this kind of brute force is not optimally for computing fast and efficient. However the script only needs to run once.
  for (i in 1:nrow(data_proteins)){
    #we store the relevant data of the protein we try to match
    charge = data_proteins$charge[i] 
    mz =  data_proteins$mz[i]
    RT = data_proteins$retentiontime[i]
    ID = data_proteins$ID[i]
    for (t in 1:nrow(data_referention)){
      #we store the relevant data of the proteins of the referention dataset.
      charge_ref = data_referention$charge[t] 
      mz_ref =  data_referention$mz[t]
      RT_ref = data_referention$retentiontime[t]
      ID_ref = data_referention$ID[t]
      #if we have a match, which states equal charge, and Mz and RT within the states limits, we add a line to the result matrix.
      if(charge_ref == charge & abs(mz_ref - mz) <= mz_threshold & abs(RT_ref - RT) <= rt_threshold){
        m.result[line_index,] = c(ID,RT,mz,mz_ref,RT_ref, abs(RT_ref - RT), abs(mz_ref - mz), ID_ref)
        line_index = line_index + 1
      }
    }
  }

  
  # we store the result in a dataframe..
  result = as.data.frame(m.result[1:line_index -1,])
  #we want to calculate the distance, but we want to weigh both factors equally. Therefore we choose the correct for the ratio of the limits.
  result$mz_delta = result$mz_delta * (rt_threshold/mz_threshold)
  #we use pythagoras to calculate the absolute distance 
  result$distance = sqrt(result$mz_delta^2 + result$RT_delta^2 )
  
  #we now select only the line with the lowest distance, which is the best match for the proeitn to the referention protein.
  result = result %>% 
    group_by(ID) %>% 
    slice(which.min(distance))
  
  #We store the result in the data proteins dataset.
  data_proteins = left_join(data_proteins, result[c("ID", "match_ID")], by = "ID")
  #we return the dataset. 
  return(data_proteins)
}



#we load the requiring data into the R global environment
data = read.csv(file = "C:/Users/Willem/Desktop/output_proteins/inputs/Overview_Proteomic_Ion.csv", sep = ";")
peptides = read.csv(file = "c:/Users/Willem/Desktop/output_proteins/inputs/verdeling_groepen.csv", sep = ";")

#we define the relevant parameters, the threshold is defined at 100.000
treshold_RA = 100000
RT_limit = 0.5
MZ_limit = 0.003

old_length = nrow(data)
## we remove the values with a lower RA than stated in thetreshold RA
data = data[data$raw.abundance >= treshold_RA,]
print(old_length - nrow(data))
old_length = nrow(data)

#we set ID's on the dataset
data$ID = seq(1:nrow(data))
#we already set the type f protein based on the mass, we use the petides table to do so.
data$type = NA
for (i in 1:nrow(peptides)){
  data$type = ifelse(data$mass<peptides$max[i] & data$mass>peptides$min[i], as.character(peptides$Naam[i]), data$type)
}


####################################### STEP 1 - Correct for digestive enzymes ###########################################

#First we correct for digestive enzyme's measured in Blanco
blanco_protein = "Blanco"

#we create a table with oly the referention and a table with only the rest of the proteins.
data_referention = data[data$Protein == blanco_protein,]
#we want to remove the BLANCO dataset from the proteins, it is of no use anymore.
data_proteins  = data[!data$Protein == blanco_protein,]

#we call the function matching, which returns a dataframe of data_proteins with the added match_ID (the ID's of banco wich has been matched)
data_match_blanco = matching(data_referention = data_referention, data_proteins = data_proteins, mz_threshold = MZ_limit, rt_threshold = RT_limit)

#we store the new data in data, because we want to lose the Blanco protein anyways, we only want to keep lines that are not matched.
data = data_match_blanco[is.na(data_match_blanco$match_ID),]
#and we remove the match_ID, because it holds no value anymore (we already delteted all the matches).
data$match_ID = NULL
print(old_length - nrow(data))



####################################### STEP 2 - Match the Proteins to WPC     ###########################################
referention_protein = "WPC"


#we re-ID the dataset, based on the referention protein first, and then mass. 
data$Protein = as.character(data$Protein)
data = data[order(ifelse(data$Protein == referention_protein, "aaaa", data$Protein), data$mass),]
data$ID = seq(from= 1, to = nrow(data))

#We select the protein WPC to be in data_referention
data_referention = data[data$Protein == referention_protein,]
#we select the other proteins in the data_proteins
data_proteins  = data[!data$Protein == referention_protein,]

#we start the matching proces, resulting in the best match for every ine in proteins, if it is wihtin the limits.
data_proteins = matching(data_referention = data_referention, data_proteins = data_proteins, rt_threshold = RT_limit, mz_threshold = MZ_limit )


########################################## STEP 3 save the WPC results, check the some ID's on matches #############################

#we write the data WPC to the harddisk.
output_map = "/output_final"
paste_wd = paste(getwd(), output_map, sep = "")
dir.create(paste_wd)
setwd(paste_wd)
write.csv(data_referention, "WPC.csv")

#we know that the match_Id is the ID in WPC
data_referention$match_ID = data_referention$ID 
#we want to do math with the matches, so WPS gets a value of -1 and the proteins a value of 1.
data_referention$value = -1
data_proteins$value = 1

#we want to get some information on certain interesting peptide sequences
peptide_sequence = as.data.frame(x = c(33,37,43,50,128))
colnames(peptide_sequence) = "ID"
peptide_sequence$`name sequence` = c("LK","DK", "VR", "MK", "MAPK")
sequence_proteins = inner_join(data_proteins, peptide_sequence, by = c("match_ID" = "ID"))

#we join the data using innre join, resulting in a filter on only ID's representing in both tables.
wpc = inner_join(data_referention, peptide_sequence, by = c("ID" = "ID"))
sequence_proteins = bind_rows(sequence_proteins, wpc)
sequence_proteins = sequence_proteins[,c("match_ID", "Protein", "mz", "charge", "retentiontime", "raw.abundance", "mass", "type", "name sequence")]
colnames(sequence_proteins)[1] = "ID"
write.csv(sequence_proteins, "Result_peptide_sequence.csv")


################################################ STEP 4 presenting the results #####################################

#we create a for loop in which we present all results for every protein (except for WPC)

for (protein in unique(data_proteins$Protein)){
  
  
  #####################################STEP 4a - plot the tilemap ###############################################################
  
  #create heatmap raster from WPC
  #join data from WPC and the protein
  data_plot = union(data_referention, data_proteins[data_proteins$Protein == protein,])
  #if the union results in a NA line, we want to remove is
  data_plot = data_plot[!is.na(data_plot$Protein),]
  
  #We want to give the non mathced lines a match_ID, zo they show up in the plot
  data_plot$match_ID = ifelse(is.na(data_plot$match_ID), data_plot$ID, data_plot$match_ID)
  
  #Two lines can be matched to the same match ID and therefore should be presented in the same ID.
  #we solve this to calculate the min and max value, we calculate the mean mass
  data_heatmap <- data_plot %>%
               group_by(match_ID) %>%
               summarise(WPC = min(value), protein = max(value), mass = mean(mass))
  
  #We want to show the data ordered by mass, we reorder the frame by the mass and the match_ID (if they have the same mass).
  data_heatmap = data_heatmap[order(round(data_heatmap$mass,4), data_heatmap$match_ID),]
  data_heatmap$ID = seq(from = 0, to = nrow(data_heatmap)-1)
  #we calculate three options: Overlap (both proteins), WPC, or the protein we are plotting
  data_heatmap$Value = as.factor(ifelse(data_heatmap$protein > 0 & data_heatmap$WPC < 0, "Overlap", 
                               ifelse(data_heatmap$WPC < 0, "WPC" , protein)))
  #we again add the type of peptide
  data_heatmap$type = ""
  for (i in 1:nrow(peptides)){
    data_heatmap$type = ifelse(data_heatmap$mass<peptides$max[i] & data_heatmap$mass>peptides$min[i], as.character(peptides$Naam[i]), data_heatmap$type)
  }
  
  #we manually draw the tilemap, so we set the basic raster size on 20
  raster_size = 20
  
  #we manually  calculate the x and y for every line of protein
  #we build from (1,1). We build from the ground up, so first x then y.
  data_heatmap$y = ceiling((data_heatmap$ID + 1 ) / raster_size)
  data_heatmap$x = (data_heatmap$ID) %%raster_size
  data_heatmap$label = data_heatmap$match_ID
  #we log the number of proteins and the number of referentions in the dataset
  number_protein = length(data_heatmap$label[data_heatmap$Value == protein]) 
  number_ref = max(data_referention$match_ID)
  #we log the numbers in a way that we take the ID's from WPC and start counting form there.
  #we only do this for the unmatched proteins. The matched proteins get the WPC ID
  #Note that we create new ID's here. We show these ID's in the data and in the plots, called label.
  #so onl proteins who are unmatched receive a label.
  data_heatmap$label[data_heatmap$Value == protein] = seq(from = number_ref + 1,to = number_protein+number_ref)
  #we adjust the x and y a bit so it looks better.
  data_heatmap$y = data_heatmap$y + 0.5
  data_heatmap$x = data_heatmap$x + 0.5
  
  
  ##we now manually calculate the locations and size of the black lines, for a little make up on the plots. This calculation is based on 
  #the x and y position of the peptides. We want to show the boundaries of the peptides in the tilemap.We create a table with these boundaries. 
  lijnen_y = data_heatmap %>%
           group_by(type) %>%
           summarise(y = max(y) + 0.5)
  lijnen_y$x_begin = 0
  lijnen_y$x_eind = raster_size
  #we calculate the positions
  for (type in lijnen_y$type){
    y = lijnen_y$y[lijnen_y$type == type] - 0.5
    x_max = max(data_heatmap$x[data_heatmap$type == type & data_heatmap$y == y]) + 0.5
    lijnen_y$x_eind[lijnen_y$type == type] = x_max 
  }
  
  #we use this color palette...
  pal <- c("#D3F5FF","#B89E86", "#82A7B3" )
  names(pal) = c("Overlap", referention_protein, protein)
  
  #and we plot the heatmap
  heatmap = ggplot(data_heatmap, aes(x, y)) + 
    geom_tile(aes(fill = Value)) +
    geom_text(aes(label=label)) +
    scale_fill_manual(values = pal) +
    geom_segment(x= 0, y= 1, xend=0, yend= max(lijnen_y$y), size =1, colour = "gray26" ) +
    geom_segment(x= 0, y= 1, xend= raster_size, yend= 1,size =1, colour = "gray26" ) +
    geom_segment(x= raster_size, y= 1, xend=raster_size, yend= max(lijnen_y$y) - 1,size =1, colour = "gray26" ) +
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y =element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_blank()) 
  
  #and we plot the extra lines
  lijnen_y = lijnen_y[order(lijnen_y$y),]  
  for (i in 1:nrow(lijnen_y)) {
    heatmap <- heatmap + geom_segment(x=lijnen_y$x_begin[i], y= lijnen_y$y[i], xend=lijnen_y$x_eind[i], yend= lijnen_y$y[i], size =1, colour = "gray26" )
    heatmap <- heatmap + geom_segment(x=lijnen_y$x_eind[i], y= lijnen_y$y[i] - 1, xend  = raster_size, yend= lijnen_y$y[i] -1, size =1, colour = "gray26" )
    heatmap <- heatmap + geom_segment(x=lijnen_y$x_eind[i], y= lijnen_y$y[i], xend  = lijnen_y$x_eind[i], yend= lijnen_y$y[i] -1, size =1, colour = "gray26" )
    y_label = ifelse(i == 1, lijnen_y$y[i] - ((lijnen_y$y[i] - 1)/2), 
                     ifelse(i == nrow(lijnen_y), max(lijnen_y$y) - (max(lijnen_y$y) - lijnen_y$y[i-1])/2, 
                            lijnen_y$y[i]- (lijnen_y$y[i] - lijnen_y$y[i-1])/2))
    heatmap <- heatmap + annotate("text", x=-.5, y =y_label, label = lijnen_y$type[i], hjust = 1)
    heatmap <- heatmap + annotate("text", x=-3, y =y_label, label = "   ", hjust = 0)
  }
  
  #finally we save our heatmap.
  name = paste(protein, "_heatmap.jpeg", sep = "")
  ggsave(name, plot = heatmap, width = 10, height = 10)
  
  
  ##################################### STEP 4b - pot the barchart ###############################################################
  #we only want to plot the data which has label Overlap
  data_both = data_plot[data_plot$match_ID %in% data_heatmap$match_ID[data_heatmap$Value == "Overlap"],]
  #We want the ID's corresponding with the tilemaps, we don't need no labels because labels are only necessary in unmatched data
  data_both$ID = as.factor(data_both$match_ID)
  
  #If two pieces of protein are matched to the same piece op WPC, we take the sum of the RA of those both pieces to show in the barchart.
  data_both = data_both %>%
    group_by(Protein, ID) %>%
    summarise(raw_abundance = sum(raw.abundance), Mass = mean(mass))
  
  #we want mass as a factor.
  data_both$Mass = as.factor(data_both$Mass)
  
  #plot de barchart
  barchart = ggplot(data = data_both, aes(x = ID, log10(raw_abundance)))+
    geom_bar(aes(fill = Protein), stat="identity", position = "dodge")  +
    coord_cartesian(ylim=c(log10(treshold_RA), max(log10(data_both$raw_abundance)))) +
    theme(panel.background = element_rect(fill = "white"),
          axis.text.y = element_text(size=25),
          axis.title.y =element_text(size=25),
          axis.text.x = element_text(angle = 45, size=14 ,hjust = 1),
          axis.title.x = element_text(size=25),
          legend.position = "none")+
    scale_fill_manual(values = pal) +
    ylab(label = "10log Peptide Abundance")
  
  #we saven the barchart
  name = paste(protein, "_barchart.jpeg", sep = "")
  ggsave(name, plot = barchart, width = 20, height = 20)
  
  #we save the data from the plots. We only save the data which is present whitin this protein
  data_output  = left_join(data_proteins[data_proteins$Protein == protein,], data_heatmap[c("match_ID", "label")], by = c("ID" = "match_ID"))
  
  #If the label is NA, a match with WPC is confirmed..
  data_output$Match = ifelse(is.na(data_output$label), "Y", "N")
  #we want to use the IDs which are shown in the plots, so if the data has a label (unmatched), we use the label as ID.
  #otherwise we use the match ID. Note that one ID can occur multiple times if two peptides within one protein are being matched to 
  #the same piece of WPC. 
  data_output$ID = ifelse(is.na(data_output$label), data_output$match_ID, data_output$label)
  #we order the columns and select only the ones we want to print
  data_output = data_output[,c(7,1,2,3,4,5,6,8,12)]
  #We order the rows based on the IDs
  data_output = data_output[order(data_output$ID),]
  
  #save the data as CSV file
  name = paste(protein, ".csv", sep = "")
  write.csv(data_output, name, row.names = F)
}
#show the location in which all the plots were saved
print(paste("geprint in map:" , getwd()))


################################################### TIME for CLUSTERING ###################################################



setwd("C:/Users/Willem/Desktop/output_final/");


######################################################### STEP 5-  PROTEOMICS CLUSTERING #############################################################

#We want all data in one table for the proteomics analysis
data_total = union_all(data_referention, data_proteins)
#create a Y if a piece of protein is matched to the referention protein (this case = WPC)
data_total$match =  ifelse(data_total$Protein == "WPC", "Y", ifelse(is.na(data_total$match_ID), "N", "Y"))

#create the total table for every type of peptide (in long format). We use integrers to round the rsults of RA.. (otherwise the clustering will get messy)
#if we look at the height of the RA this is not a problem >10000
proteomics = data_total %>%
  group_by(Protein, type) %>%
  summarise(RA_total = as.integer(sum(raw.abundance)))
#convert long format to pivot
proteomics = dcast(proteomics, Protein ~ type, fun.aggregate = sum, value.var = "RA_total")
#we want to set the names for the RA_totals
names(proteomics) = paste(names(proteomics), "_RA_total", sep = "")
names(proteomics)[1] = "Protein"


#to calculate the overlap we first calculate the subtotals for matched and non-matched
overlap = data_total %>%
  group_by(Protein, match) %>%
  summarise(number_peptides = n(), abundance_peptides = sum(raw.abundance))
#now we create a total column to calculate the fraction.
overlap = overlap %>% 
  group_by(Protein) %>% 
  mutate(Total_number_peptides = sum(number_peptides), Total_abundance_peptides = sum(abundance_peptides) )
#we don't need the Not matched column any more
overlap = overlap[overlap$match == "Y",]

#we calculate the percentages
overlap$Match_number_percentage = 100 * overlap$number_peptides/overlap$Total_number_peptides
overlap$Match_abundance_percentage = 100 * overlap$abundance_peptides/overlap$Total_abundance_peptides
# we join the relevant information to the totals. 
overlap = overlap[c("Protein", "Match_number_percentage", "Match_abundance_percentage", "Total_number_peptides", "Total_abundance_peptides")]
proteomics= left_join(proteomics, overlap, by = "Protein")

#Rownames for proteitn to match the clustering
rownames(proteomics) = proteomics$Protein
proteomics$Protein = NULL

#we want to measure the fractions of certain pepties within the protein
proteomics$dipeptide_RA_total =  proteomics$dipeptide_RA_total / proteomics$Total_abundance_peptides 
proteomics$`di/tripeptide_RA_total` = proteomics$`di/tripeptide_RA_total` / proteomics$Total_abundance_peptides 
proteomics$`tri/oligopeptide_RA_total` = proteomics$`tri/oligopeptide_RA_total` / proteomics$Total_abundance_peptides 
proteomics$oligopeptide_RA_total =  proteomics$oligopeptide_RA_total / proteomics$Total_abundance_peptides 
#set the correct column order and names
proteomics = proteomics[,c(2,1,4,3,8,7,5,6)]
colnames(proteomics) = c("fraction dipeptide", "fraction di/tripeptide", "fraction tri/oligopeptide", "fraction oligopeptide", "total peptide abundance", 
                         "total unique peptides",  "relative peptide overlap WPC", "peptide abundance overlap WPC")

#we save proteomics file for the supplementaries..
write.csv(proteomics, "supplementary_proteomics.csv")

#we scale the data around zero, to measure every coumn equally
proteomics.std = apply(proteomics, 2, scale, center=TRUE, scale=TRUE)

#we set the rownames 
rownames(proteomics.std)<-rownames(proteomics)


#clustering and some meta information, set the correct column names
#make the clustering
mod2 =  Mclust(proteomics.std)
mod2$classification
#plot the BIC
plot(mod2$BIC)

#we build the grid plot to determine the best 2D representation.
plot(mod2,
     fillEllipses = T, CEX = 0.75, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "blue", "brown", "grey", "pink", "blue"),
     symbols = c(15))

#we build the plot, we use the lower plot to determine the best 2D representation of the clustering (in de dimens parameter)
plot(mod2, dimens=c(7,5),
     fillEllipses = T, CEX = 1, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "blue", "brown", "grey", "pink"),
     symbols = c(15))
text(mod2$data[,7], mod2$data[,5], labels = rownames(proteomics.std), pos = 1)

#we built another 2D graph
#we build the plot, we use the lower plot to determine the best 2D representation of the clustering (in de dimens parameter)
plot(mod2, dimens=c(8,4),
     fillEllipses = T, CEX = 1, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "blue", "brown", "grey", "pink"),
     symbols = c(15)) 
text(mod2$data[,8], mod2$data[,4], labels = rownames(proteomics.std), pos = 1)
#we saven the scaled clusterdata
proteomics.std  = as.data.frame(proteomics.std)
proteomics.std$clustergroup = mod2$classification
write.csv(proteomics.std, "Proteomics_cluster_data.csv")




################################################# STEP 6 - in vitro BIOLOGICAL DATA #######################################################################

#load the data
dat = as.data.frame(fread("c:/Users/Willem/Desktop/output_proteins/inputs/Data_with_interpolation.csv")[,1:8])
names(dat)[1]<-"Protein";

print(as.data.frame(dat[1:10,]));## Show the first 10 records in the data
print(dim(dat)); ## print out the number of records x number of variables


#calcualte the means per protein
dat.list<-split(dat,list(dat$Protein));
Bio_means<-matrix(data=rep(NA,length(dat.list)*(ncol(dat)-1)),ncol=ncol(dat)-1);
rownames(Bio_means)<-names(dat.list)
colnames(Bio_means)<-names(dat)[-1];

for(i in 1:length(dat.list)){
  curDat<-dat.list[[i]];
  Bio_means[i,]<-apply(curDat[,-1],2,mean,na.rm=TRUE);
}

#scale the data around zero
Bio_means.std<-apply(Bio_means, 2, scale, center=TRUE, scale=TRUE); # this is a nice shortcut
rownames(Bio_means.std)<-rownames(Bio_means);

#drop the nan columns
dropme<-which(is.nan(Bio_means.std[1,]));
if(length(dropme)){
  Bio_means.std<-Bio_means.std[,-dropme];
}


#Cluster the data
cluster_bio<- Mclust(Bio_means.std)
cluster_bio$classification
cluster_bio$BIC

#plot the BIC
plot(cluster_bio$BIC)

#plot some 2D representations
plot(cluster_bio, dimens=c(5,7),
     fillEllipses = T, CEX = 1, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "grey", "blue"),
     symbols = c(15,16,17,18))
text(cluster_bio$data[,5], cluster_bio$data[,7], labels = rownames(Bio_means.std), pos = 1)

#plot some 2D representations
plot(cluster_bio, dimens=c(4,6),
     fillEllipses = T, CEX = 1, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "grey", "blue"),
     symbols = c(15,16,17,18))
text(cluster_bio$data[,4], cluster_bio$data[,6], labels = rownames(Bio_means.std), pos = 1)


#plot the grid
plot(cluster_bio,
     fillEllipses = T, CEX = 0.75, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "grey", "blue"),
     symbols = c(15,16,17,18))


#print the data
biological_data.std = as.data.frame(Bio_means.std)
biological_data.std$classification = cluster_bio$classification
write.csv(biological_data.std, "Cluster_data_biological_efficacy.csv")


######################################## STEP 7 - plot the combined proteomics and biological data #######################################

#we want to combine both sets, but the proteomic set withous the peptide_types....
proteomics.std$Protein = rownames(proteomics.std)
biological_data.std$Protein = rownames(biological_data.std)

#for the combined plot we lose the information about the peptides
proteomics.std$`fraction dipeptide` = NULL
proteomics.std$`fraction tri/oligopeptide` = NULL
proteomics.std$`fraction oligopeptide` = NULL
proteomics.std$`fraction di/tripeptide` = NULL

#we join the two datasets together
tot.std = left_join(proteomics.std, biological_data.std, by = "Protein")

#we only wnat quantative data columns, so no groups or classifications or proteins
tot.std$clustergroup = NULL
tot.std$classification = NULL
rownames(tot.std) = tot.std$Protein
tot.std$Protein = NULL

#perform the clustering
total_clust<- Mclust(tot.std)
total_clust$classification

  
#plot the BIC
plot(total_clust$BIC)

#plot the grid to view the results
plot(total_clust,
     fillEllipses = T, CEX = 0.75, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "grey", "blue"),
     symbols = c(15,16,17,18))

#plot some 2D representations
plot(total_clust, dimens=c(1,5),
     fillEllipses = T, CEX = 1, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "grey", "blue"),
     symbols = c(15,16,17,18))
text(total_clust$data[,1], total_clust$data[,5], labels = rownames(tot.std), pos = 1)


#plot some 2D representations
plot(total_clust, dimens=c(4,8),
     fillEllipses = T, CEX = 1, what = "classification", colors = c("gold2", "cyan3", "firebrick3", "green4", "grey", "blue"),
     symbols = c(15,16,17,18))
text(total_clust$data[,4], total_clust$data[,8], labels = rownames(tot.std), pos = 1)


#save the data as results
tot.std = as.data.frame(tot.std)
tot.std$classification = total_clust$classification
write.csv(tot.std, "Total_cluster.csv")



