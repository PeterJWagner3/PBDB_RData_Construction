# load source files & data sets ####
common_r_source <- "/Users/peterjwagner/Documents/R_Projects/Common_R_Source_Files";
data_for_r <- "/Users/peterjwagner/Documents/R_Projects/Data_for_R";
this_folder <- "/Users/peterjwagner/Documents/R_Projects/PBDB_RData_Construction";
setwd(this_folder);
source(paste(common_r_source,"/Chronos.r",sep=""));
source(paste(common_r_source,"/Data_Downloading_v4.r",sep=""));
source(paste(common_r_source,"/Occurrence_Data_Routines.r",sep=""));
source(paste(common_r_source,"/Stratigraphy.r",sep=""));
source(paste(common_r_source,"/Update_Paleobiology_Database_RData.r",sep=""));
source(paste(common_r_source,"/Wagner_Kluges.r",sep=""));

load(paste(data_for_r,"/Gradstein_2020_Augmented.RData",sep=""));
load(paste(data_for_r,"/PaleoDB_Edits.RData",sep=""));
load(paste(data_for_r,"/Rock_Unit_Database.Rdata",sep=""));
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
options(warn=1);

pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
pbdb_opinions <- pbdb_data_list$pbdb_opinions;

# Execute program ####
taxon <- "Acanthoceratidae";
subtaxon_rank <- "subfamily";
related_taxa <- pbdb_taxonomy[pbdb_taxonomy$parent_no==pbdb_taxonomy$parent_no[pbdb_taxonomy$taxon_name==taxon],];
cc <- match(subtaxon_rank,taxonomic_rank);
cg <- match("genus",taxonomic_rank);
for (i in cc:cg)	related_taxa <- rbind(related_taxa,pbdb_taxonomy[pbdb_taxonomy$parent_no %in% related_taxa$accepted_no,]);
related_taxa <- related_taxa[!related_taxa$accepted_rank %in% c("species","subspecies"),]
if (!subtaxon_rank %in% colnames(related_taxa))
	related_taxa <- add_higher_taxon_to_pbdb_taxonomy(taxon_rank=subtaxon_rank,paleodb_taxonomy=related_taxa);
related_taxon_opinions <- pbdb_opinions[pbdb_opinions$orig_no %in% related_taxa$taxon_no | pbdb_opinions$child_spelling_no %in% related_taxa$taxon_no,];

related_genera <- related_taxa[related_taxa$family %in% taxon & related_taxa$taxon_rank %in% c("genus","subgenus"),];
related_genus_opinions <- pbdb_opinions[pbdb_opinions$orig_no %in% related_genera$taxon_no | pbdb_opinions$child_spelling_no %in% related_genera$taxon_no,];
taxon_members <- related_taxa[related_taxa$family %in% taxon,];
taxon_members <- unique(taxon_members);
#acanthoceratid_all <- add_higher_taxon_to_pbdb_taxonomy(taxon_rank="subfamily",paleodb_taxonomy=acanthoceratid_all);

treatise <- related_genus_opinions[related_genus_opinions$reference_no==28273,]
relv_treatise <- rbind(related_genera[related_genera$taxon_no %in% treatise$orig_no,],
					   related_genera[related_genera$taxon_no %in% treatise$child_spelling_no,],
					   related_genera[related_genera$orig_no %in% treatise$orig_no,],
					   related_genera[related_genera$orig_no %in% treatise$child_spelling_no,],
					   related_genera[related_genera$accepted_no %in% treatise$orig_no,],
					   related_genera[related_genera$accepted_noorig_no %in% treatise$child_spelling_no,]);
relv_treatise <- unique(relv_treatise);
pbdb_opinions$reference_no[pbdb_opinions$taxon_name=="Buchiceras"][3] <- 28273;
related_genera <- taxon_members[taxon_members$taxon_rank %in% c("genus","subgenus"),];
related_genera <- related_genera[related_genera$difference=="",];
ngenera <- nrow(related_genera);
subtaxon_history <- data.frame(genus=as.character(related_genera$taxon_name),current=as.character(rep("",ngenera)),treatise=as.character(rep("",ngenera)));
for (i in 1:nrow(related_genera))	{
	subtaxon_history$current[i] <- related_genera$subfamily[i]
	parent <- relv_treatise$parent_name[relv_treatise$orig_no %in% related_genera$orig_no[i]][1];
	if (!is.na(parent))	{
		while (taxon_members$taxon_rank[taxon_members$taxon_name %in% parent] %in% c("genus","tribe"))	{
			parent <- relv_treatise$parent_name[relv_treatise$orig_no %in% related_genera$orig_no[match(parent,related_genera$taxon_name)]];
			}
		if (taxon_members$taxon_rank[taxon_members$taxon_name %in% parent]=="subfamily")	{
			subtaxon_history$treatise[i] <- parent;
			}
		}
	}

write.csv(subtaxon_history,"Treatise_to_Today.csv",row.names = F)


mertz_taxa <- c("Acanthoceras","Cunningtoniceras","Nigericeras","Morrowites","Pseudaspidoceras","Mammites","Dunveganoceras","Spathites","Calycoceras","Sumitomoceras","Conlinoceras","Romaniceras","Kamerunoceras","Paraconlinoceras","Euomphaloceras","Tarrantoceras","Eucalycoceras","Pseudocalycoceras","Plesiacanthoceras","Watinoceras","Metoicoceras","Neocardioceras")
write.csv(cbind(mertz_taxa,pbdb_taxonomy$parent_name[match(mertz_taxa,pbdb_taxonomy$taxon_name)]),"Mertz_Tree.csv",row.names=FALSE)
{}
