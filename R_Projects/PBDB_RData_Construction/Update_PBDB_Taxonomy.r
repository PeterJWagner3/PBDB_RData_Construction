#install.packages("taskscheduleR")
library(cronR);

# load source files & data sets ####
data_for_r <- getwd();
data_for_r <- strsplit(data_for_r,"/")[[1]];
data_for_r <- paste(data_for_r[1:match("R_Projects",data_for_r)],collapse="/");
common_r_source <- paste(data_for_r,"/Common_R_Source_Files",sep="");
data_for_r <- paste(data_for_r,"/Data_for_R",sep="");
source(paste(common_r_source,"/Chronos.r",sep=""));
source(paste(common_r_source,"/Data_Downloading_v4.r",sep=""));
source(paste(common_r_source,"/Occurrence_Data_Routines.r",sep=""));
source(paste(common_r_source,"/Stratigraphy.r",sep=""));
source(paste(common_r_source,"/Update_Paleobiology_Database_RData.r",sep=""));
source(paste(common_r_source,"/Wagner_Kluges.r",sep=""));

load(paste(data_for_r,"/Gradstein_2020_Augmented.RData",sep=""));
load(paste(data_for_r,"/PaleoDB_Edits.RData",sep=""));
load(paste(data_for_r,"/Rock_Unit_Database.Rdata",sep=""));
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));

# Update Taxonomy for one taxon ####
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
load(paste(data_for_r,"/Gradstein_2020_Augmented.RData",sep=""));
load(paste(data_for_r,"/PaleoDB_Edits.RData",sep=""));
load(paste(data_for_r,"/Rock_Unit_Database.Rdata",sep=""));
pbdb_taxonomy_edits <- paleodb_fixes$paleodb_taxonomy_edits;
taxon <- "Cephalopoda";
bug_taxa <- c("Deuteropoda","Lobopoda","Onychophora","Radiodonta","Artiopoda");
for (bt in 1:length(bug_taxa))	{
#	options(warn=0);
#	load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon=bug_taxa[bt],pbdb_data_list,pbdb_taxonomy_edits);
#	pbdb_data_list$pbdb_taxonomy[pbdb_data_list$pbdb_taxonomy$taxon_no>1000000,];
	save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
#	options(warn=1);
	}
mollusc_taxa <- c("Bivalvia","Cephalopoda","Gastropoda","Polyplacophora","Rostroconchia","Helcionelloida","Tergomya","Maikhanellidae");
for (mt in 1:length(mollusc_taxa))	{
#	options(warn=0);
#	load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon=mollusc_taxa[mt],pbdb_data_list,pbdb_taxonomy_edits);
	save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
#	options(warn=1);
	}
ct <- 1;
chordate_taxa <- c("Cephalochordata","Conodontophorida","Craniata","Metaspriggina","Olfactores","Tunicata","Undichna","Acrania","Agnatha","Anaspidomorphi","Cephalaspidomorpha","Craniota","Cyclostomi","Eotetrapoda","Galeaspida","Herpetichthyes","Myzichthyes","Osteostracomorphii","Ostracophori","Pachycardia","Petromyzontida","Petromyzontomorphi","Chondrichthiomorphi","Chondrichthyes","Eugnathostomata","Ichthyes","Placodermiomorphi","Plectrodus","Qilinyu","Thelodonti","Actinopterygii","Coregonidae","Ganoidei","Hypostomata","Iniomi","Teleostomata","Sarcopterygii");
while (ct <= length(chordate_taxa))	{
#	options(warn=0);
#	load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon=chordate_taxa[ct],pbdb_data_list,pbdb_taxonomy_edits);
	save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
	ct <- ct+1;
#	options(warn=1);
	}
options(warn=1);
#pbdb_data_list$pbdb_finds <- pbdb_data_list$pbdb_finds[!pbdb_data_list$pbdb_finds$occurrence_no %in% 1461783,]

load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
pbdb_data_list_smol$pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
pbdb_data_list_smol$pbdb_finds <- pbdb_data_list$pbdb_finds;
save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));

#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
# Update Taxonomy for basic taxonomic groups ####
higher_taxa_to_update <- c("Ichnofossils","Problematica","Cyanobacteria","Retaria","Rhodophyta",
						   "Arboreomorpha","Rangeomorpha","Erniettomorpha","Tetraradialomorpha","Bilaterialomorpha",
						   "Porifera","Cnidaria","Lophophorata","Mollusca",
						   "Panarthropoda","Chordata","Hemichordata","Echinodermata");
#for (ht in 6:7)	{
#ht <- match("Artiopoda",higher_taxa_to_update);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
for (ht in 14:length(higher_taxa_to_update))	{
#for (ht in 16:16)	{
#	load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
#	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon="Chordata",pbdb_data_list,pbdb_taxonomy_edits);
	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon=higher_taxa_to_update[ht],pbdb_data_list,pbdb_taxonomy_edits);
#	pbdb_data_list$pbdb_taxonomy[pbdb_data_list$pbdb_taxonomy$taxon_name=="Erniettamorpha",]
#	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon="Marrellomorpha",pbdb_data_list);
	load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
	pbdb_data_list_smol$pbdb_finds <- pbdb_data_list$pbdb_finds;
	pbdb_data_list_smol$pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
	save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));
	}

orphan_taxonomy_info <- unorphan_species(pbdb_data_list$pbdb_taxonomy);
orphan_taxonomy <- orphan_taxonomy_info$all_orphan_taxonomy;
pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% orphan_taxonomy$taxon_no,] <- orphan_taxonomy;
pbdb_data_list_smol$pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));
pbdb_data_list$lost_taxon_numbers <- perm_missing;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));


#pbdb_finds <- expello_na_from_matrix(pbdb_finds,replacement="");

# add missing taxa to taxonomy table #
data_for_r <- getwd();
data_for_r <- strsplit(data_for_r,"/")[[1]];
data_for_r <- paste(data_for_r[1:match("R_Projects",data_for_r)],collapse="/");
data_for_r <- paste(data_for_r,"/Data_for_R",sep="");
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
latest_taxon_no <- max(pbdb_taxonomy$taxon_no);
missing_taxon_nos <- (1:latest_taxon_no)[!(1:latest_taxon_no) %in% pbdb_taxonomy$taxon_no];
permantently_missing <- read.csv("Lost_Taxon_Nos.csv",header=T);
perm_missing <- unique(permantently_missing[,1]);
missing_taxon_nos <- missing_taxon_nos[!missing_taxon_nos %in% perm_missing];
mtaxa <- length(missing_taxon_nos);
mt <- 0;
markers <- round(length(missing_taxon_nos)*(1:99)/100,0)
#mt <- 1001;
while (mt < mtaxa)	{
	mt <- mt+1;
	if (mt %in% markers)	{
		print(paste(match(mt,markers),"% done",sep=""));
		pbdb_taxonomy <- put_pbdb_dataframes_into_proper_type(pbdb_taxonomy);
		pbdb_taxonomy <- clean_the_bastards(pbdb_taxonomy);
		pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
		save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
		}
	add_this_taxon <- accersi_taxonomic_data_for_one_taxon_no(taxon_id=missing_taxon_nos[mt]);
	add_this_taxon <- clean_the_bastards(pbdb_data=add_this_taxon);
	pbdb_taxonomy <- rbind(pbdb_taxonomy,add_this_taxon);
	still_missing <- missing_taxon_nos[1:mt][!missing_taxon_nos[1:mt] %in% pbdb_taxonomy$taxon_no];
	perm_missing <- unique(c(perm_missing,still_missing));
	write.csv(perm_missing,"Lost_Taxon_Nos.csv",row.names = F);
	}

pbdb_taxonomy <- put_pbdb_dataframes_into_proper_type(pbdb_taxonomy);
pbdb_taxonomy <- clean_the_bastards(pbdb_taxonomy);
pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

saveRDS(pbdb_data_list,paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
#dummy <- file.choose()
#save(pbdb_data_list,file=dummy);
#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

# quick fixes
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
pbdb_data_list$pbdb_finds$class_no[pbdb_data_list$pbdb_finds$order %in% "Rhombifera"] <- pbdb_data_list$pbdb_finds$class_no[pbdb_data_list$pbdb_finds$class=="Rhombifera"][1]
pbdb_data_list$pbdb_finds$class[pbdb_data_list$pbdb_finds$order %in% "Rhombifera"] <- "Rhombifera";
pbdb_data_list$pbdb_finds$order_no[pbdb_data_list$pbdb_finds$order %in% "Rhombifera"] <- 0;
pbdb_data_list$pbdb_finds$order[pbdb_data_list$pbdb_finds$order %in% "Rhombifera"] <- "NO_ORDER_SPECIFIED";
pbdb_data_list$pbdb_taxonomy$class_no[pbdb_data_list$pbdb_taxonomy$order %in% "Rhombifera"] <- pbdb_data_list$pbdb_taxonomy$class_no[match("Rhombifera",pbdb_data_list$pbdb_taxonomy$class)]
pbdb_data_list$pbdb_taxonomy$class[pbdb_data_list$pbdb_taxonomy$order %in% "Rhombifera"] <- "Rhombifera"
pbdb_data_list$pbdb_taxonomy$order_no[pbdb_data_list$pbdb_taxonomy$order %in% "Rhombifera"] <- 0;
pbdb_data_list$pbdb_taxonomy$order_[pbdb_data_list$pbdb_taxonomy$order %in% "Rhombifera"] <- "NO_ORDER_SPECIFIED";
pbdb_data_list$pbdb_taxonomy$flags[pbdb_data_list$pbdb_taxonomy$orig_no==30998] <- "V";
pbdb_data_list$pbdb_taxonomy$orig_no[pbdb_data_list$pbdb_taxonomy$taxon_no==30938] <- 30998;

load(file.choose());
pbdb_data_list_smol$pbdb_finds$accepted_name[pbdb_data_list_smol$pbdb_finds$occurrence_no==149150]
pbdb_data_list_smol$pbdb_finds$accepted_name <- gsub("Echinosphaerites \\(Echinosphaerites\\)","Echinosphaerites",pbdb_data_list_smol$pbdb_finds$accepted_name);
pbdb_data_list_smol$pbdb_finds$accepted_name[pbdb_data_list_smol$pbdb_finds$occurrence_no==149150]

pbdb_data_list$pbdb_finds$accepted_name <- pbdb_data_list_smol$pbdb_finds$accepted_name;
pbdb_data_list$pbdb_finds$accepted_name[pbdb_data_list$pbdb_finds$occurrence_no==149150]
#pbdb_data_list$pbdb_finds$accepted_name <- gsub("Echinosphaerites (Echinosphaerites)","Echinosphaerites",pbdb_data_list$pbdb_finds$accepted_name)
pbdb_data_list$pbdb_finds$accepted_name <- gsub("Echinosphaerites \\(Echinosphaerites\\)","Echinosphaerites",pbdb_data_list$pbdb_finds$accepted_name)

pbdb_data_list$pbdb_finds$accepted_name <- gsub("Echinosphaerites (Echinosphaerites)","Echinosphaerites",pbdb_data_list$pbdb_finds$accepted_name);
pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$accepted_name_org %in% "Echinosphaerites (Echinosphaerites)",]
pbdb_data_list$pbdb_finds$genus[pbdb_data_list$pbdb_finds$genus %in% "Echinosphaerites (Echinosphaerites)"] <- "Echinosphaerites";

pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$occurrence_no==149150,]
pbdb_data_list_smol <- list(pbdb_data_list$pbdb_sites_refined,pbdb_data_list$pbdb_finds,pbdb_data_list$pbdb_taxonomy,pbdb_data_list$time_scale[pbdb_data_list$time_scale$scale=="Stage Slice",]);
names(pbdb_data_list_smol) <- c("pbdb_sites","pbdb_finds","pbdb_taxonomy","time_scale");
save(pbdb_data_list_smol,file=paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
{}

pbdb_finds <- pbdb_data_list$pbdb_finds;
pbdb_data_list$pbdb_finds$accepted_name <- gsub("\\( ","\\(",pbdb_data_list$pbdb_finds$accepted_name)
taxon_name <- pbdb_data_list$pbdb_finds$accepted_name;
subgenus_finds <- pbapply::pbsapply(taxon_name,is.species.in.subgenus);
unique_subgeneraed_species <- sort(unique(taxon_name[subgenus_finds]));
species_name <- unique_subgeneraed_species;
#species_name <- "Acadoparadoxides (Acadoparadoxides) pentagonalis"
run_markers <- round((length(unique_subgeneraed_species)*(1:99)/100),0);
for (uss in 1:length(unique_subgeneraed_species))	{
	if (uss %in% run_markers)	print(paste(100*round((uss/length(unique_subgeneraed_species)),2),"% done",sep=""));
	pbdb_finds <- standardize_genus_subgenus_combos(species_name=unique_subgeneraed_species[uss],pbdb_finds);
	}
#pbapply::pbsapply(species_name,standardize_genus_subgenus_combos,pbdb_finds)

pbdb_data_list_smol$pbdb_taxonomy[pbdb_data_list_smol$pbdb_taxonomy$taxon_name %in% "Keratosa",]
pbdb_data_list_smol$pbdb_taxonomy$accepted_rank[pbdb_data_list_smol$pbdb_taxonomy$taxon_name %in% "Keratosa"] <- "subclass";
pbdb_taxonomy$accepted_no[pbdb_taxonomy$taxon_name %in% "Keratosa"] <- 463040;
pbdb_taxonomy$flags[pbdb_taxonomy$taxon_no==273399] <- "V";

load("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData"); # refined PBDB
pbdb_data_list_smol$pbdb_taxonomy <- pbdb_taxonomy;
save(pbdb_data_list_smol,file=paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));

# Fix Dates ####
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
nsites <- nrow(pbdb_data_list$pbdb_sites);
modified <- created <- rep(convert_pbdb_date(pbdb_data_list$pbdb_sites$created[1]),nsites);
for (pd in 1:nsites)	{
	if (pd %% 5000==0)	print(pd);
	created[pd] <- convert_pbdb_date(pbdb_data_list$pbdb_sites$created[pd]);
	modified[pd] <- convert_pbdb_date(pbdb_data_list$pbdb_sites$modified[pd]);
	}
cn <- match("created",colnames(pbdb_data_list$pbdb_sites));
pbdb_data_list$pbdb_sites$created <- pbdb_data_list$pbdb_sites$modified <- NULL;
#tibble::add_column(paleodb_finds, subgenus=as.character(subgenus), .after = match("genus_no",colnames(paleodb_finds)));
pbdb_data_list$pbdb_sites <- tibble::add_column(pbdb_data_list$pbdb_sites, modified=modified, .after = match("modifier",colnames(pbdb_data_list$pbdb_sites)));
pbdb_data_list$pbdb_sites <- tibble::add_column(pbdb_data_list$pbdb_sites, created=created, .after = match("modifier",colnames(pbdb_data_list$pbdb_sites)));

nfinds <- nrow(pbdb_data_list$pbdb_finds);
#effed <- (1:nfinds)[created==""];
#for (ef in 1:length(created==""))	created[effed[ef]] <- created[effed[ef]-1];
modified <- created <- rep(convert_pbdb_date(pbdb_data_list$pbdb_finds$created[1]),nfinds);
for (pd in 1:nfinds)	{
	if (pd %% 5000==0)	print(pd);
	if (pbdb_data_list$pbdb_finds$created[pd]=="" || is.na(pbdb_data_list$pbdb_finds$created[pd]))	{
		pbdb_data_list$pbdb_finds$created[pd] <- pbdb_data_list$pbdb_finds$created[pd-1];
		} else	{
		created[pd] <- convert_pbdb_date(pbdb_data_list$pbdb_finds$created[pd]);
		}
	if (pbdb_data_list$pbdb_finds$created[pd]=="" || is.na(pbdb_data_list$pbdb_finds$created[pd]))	{
		pbdb_data_list$pbdb_finds$modified[pd] <- pbdb_data_list$pbdb_finds$created[pd];
		} else	{
		modified[pd] <- convert_pbdb_date(pbdb_data_list$pbdb_finds$modified[pd]);
		}
	}
pbdb_data_list$pbdb_finds$created <- pbdb_data_list$pbdb_finds$modified <- NULL;
pbdb_data_list$pbdb_finds <- tibble::add_column(pbdb_data_list$pbdb_finds, modified=modified, .after = match("abund_unit",colnames(pbdb_data_list$pbdb_finds)));
pbdb_data_list$pbdb_finds <- tibble::add_column(pbdb_data_list$pbdb_finds, created=created, .after = match("abund_unit",colnames(pbdb_data_list$pbdb_finds)));
pbdb_data_list$pbdb_finds$modified[is.na(pbdb_data_list$pbdb_finds$modified)] <- pbdb_data_list$pbdb_finds$created[is.na(pbdb_data_list$pbdb_finds$created)]
effed <- (1:nfinds)[is.na(pbdb_data_list$pbdb_finds$created)];
for (ef in 1:length(effed))	{
	pbdb_data_list$pbdb_finds$created[effed[ef]] <- pbdb_data_list$pbdb_finds$created[effed[ef]-1];
	if (is.na(pbdb_data_list$pbdb_finds$modified[effed[ef]]))
		pbdb_data_list$pbdb_finds$modified[effed[ef]] <- pbdb_data_list$pbdb_finds$created[effed[ef]];
	}

nsites <- nrow(pbdb_data_list$pbdb_sites_refined);
pbdb_data_list$pbdb_sites_refined$created <- as.Date.character(pbdb_data_list$pbdb_sites$created[match(pbdb_data_list$pbdb_sites_refined$collection_no,pbdb_data_list$pbdb_sites$collection_no)]);
pbdb_data_list$pbdb_sites_refined$modified <- as.Date.character(pbdb_data_list$pbdb_sites$modified[match(pbdb_data_list$pbdb_sites_refined$collection_no,pbdb_data_list$pbdb_sites$collection_no)]);

created_effed <- (1:nsites)[pbdb_sites$created==""];
created_effed <- c(created_effed,(1:nsites)[!is.na(as.numeric(pbdb_sites$created))]);
created_effed <- c(created_effed,(1:nsites)[is.na(pbdb_sites$created)]);
created_good <- (1:nsites)[!(1:nsites) %in% created_effed];
for (ce in 1:length(created_effed))	{
	i <- sum(created_good<created_effed[ce]);
	new_date <- convert_pbdb_date(strsplit(as.character(pbdb_sites$created[created_good[i]])," ")[[1]][1]);
	pbdb_sites$created[created_effed[ce]] <- as.character(as.Date(strsplit(as.character(new_date)," ")[[1]][1]));
	}
pbdb_date <- pbdb_sites$created[1];
convert_pbdb_date(pbdb_sites$created[1])
new_date <- as.Date.character(rep("",nsites));
for (i in 1:nsites)	new_date[i] <- convert_pbdb_date(pbdb_sites$created[i]);
pbdb_sites$created <- new_date;

pbdb_sites$modified[pbdb_sites$modified==""] <- pbdb_sites$created[pbdb_sites$modified==""];
modified_effed <- (1:nsites)[is.na(pbdb_sites$modified)];
pbdb_sites$modified[modified_effed] <- pbdb_sites$created[modified_effed];
new_date <- as.Date.character(rep("",nsites));
for (i in 1:nsites)	new_date[i] <- convert_pbdb_date(pbdb_sites$modified[i]);
pbdb_sites$modified <- new_date;

pbdb_sites_refined$created <- as.Date.character(pbdb_sites$created[match(pbdb_sites_refined$collection_no,pbdb_sites$collection_no)]);
pbdb_sites_refined$modified <- as.Date.character(pbdb_sites$modified[match(pbdb_sites_refined$collection_no,pbdb_sites$collection_no)]);

new_date <- strsplit(as.character(new_date)," ")[[1]][1]
pbdb_sites_refined[pbdb_sites_refined$collection_no==117142,colnames(pbdb_sites_refined) %in% colnames(pbdb_sites)] <- pbdb_sites[pbdb_sites$collection_no==117142,];
pbdb_sites_refined$rock_no[pbdb_sites_refined$collection_no==117142] <- 11729;
pbdb_sites_refined$rock_no_sr[pbdb_sites_refined$collection_no==117142] <- 11729;
pbdb_sites_refined$formation_no[pbdb_sites_refined$collection_no==117142] <- 11729;
pbdb_sites_refined$rock_unit_senior[pbdb_sites_refined$collection_no==117142] <- "E-Lert";
pbdb_sites_refined$ma_lb[pbdb_sites_refined$collection_no==117142] <- 276.6;
pbdb_sites_refined$ma_ub[pbdb_sites_refined$collection_no==117142] <- 274.35;
pbdb_sites_refined$interval_lb[pbdb_sites_refined$collection_no==117142] <- "Ku2";
pbdb_sites_refined$pbdb_rock_no[pbdb_sites_refined$collection_no==117142] <- 5859;
pbdb_sites_refined$pbdb_rock_no_sr[pbdb_sites_refined$collection_no==117142] <- 5859;
pbdb_sites_refined$pbdb_formation_no[pbdb_sites_refined$collection_no==117142] <- 5859;
pbdb_sites_refined[pbdb_sites_refined$collection_no==117142,];

# update rock ages ####
rock_database <- rock_unit_data$rock_unit_database;
pbdb_sites_refined <- redate_paleodb_collections_with_time_scale(paleodb_collections=pbdb_data_list$pbdb_sites_refined,time_scale=gradstein_2020_emended$time_scale,zone_database=gradstein_2020_emended$zones,stratchron="International");
pbdb_sites_refined <- redate_paleodb_collections_with_zones(paleodb_collections = pbdb_sites_refined,zone_database = gradstein_2020_emended$zones,time_scale=gradstein_2020_emended$time_scale);
pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
pbdb_data_list <- fix_sites_based_on_rock_data_only(pbdb_data_list,rock_database);

pbdb_sites_refined[pbdb_sites_refined$collection_no==158966,]
rock_unit_data$rock_unit_database[rock_unit_data$rock_unit_database$rock_no==20981,]

pbdb_sites_refined[pbdb_sites_refined$collection_no %in% 158966,]
pbdb_sites_refined$formation[pbdb_sites_refined$collection_no %in% 158966:158981]
pbdb_sites_refined$rock_no[pbdb_sites_refined$collection_no %in% 158966:158981] <- 2234;
pbdb_sites_refined$rock_no_sr[pbdb_sites_refined$collection_no %in% 158966:158981] <- 2234;
pbdb_sites_refined$formation_no[pbdb_sites_refined$collection_no %in% 158966:158981] <- 2234;

pbdb_sites_refined[pbdb_sites_refined$collection_no==170479,];
rock_unit_data$rock_unit_database$ma_ub[rock_unit_data$rock_unit_database$rock_no_sr %in% 4614] <- 431.75;
pbdb_sites_refined[pbdb_sites_refined$collection_no==174335,];
rock_unit_data$rock_unit_database[rock_unit_data$rock_unit_database$rock_no_sr %in% 19567,];
rock_unit_data$rock_unit_database$ma_ub[rock_unit_data$rock_unit_database$rock_no_sr %in% 19567] <- 95.9

pbdb_sites_refined[pbdb_sites_refined$collection_no %in% 226884:226887,]
pbdb_sites_refined$zone[pbdb_sites_refined$collection_no %in% 226884:226887] <- "Inoceramus kamuy";
pbdb_sites_refined$zone_species[pbdb_sites_refined$collection_no %in% 226884:226887] <- "kamuy";
rock_unit_data$rock_unit_database[rock_unit_data$rock_unit_database$rock_no_sr %in% 19736,];

#pbdb_sites_refined$zone[pbdb_sites_refined$collection_no %in% 158966:158981]
#pbdb_sites_refined$ma_lb[pbdb_sites_refined$collection_no %in% 158966:158981]
#pbdb_sites_refined$ma_ub[pbdb_sites_refined$collection_no %in% 158966:158981]
#time_scale <- gradstein_2012_emended$time_scale;
#finest_chronostrat <- subset(gradstein_2012_emended$time_scale,gradstein_2012_emended$time_scale$scale=="Stage Slice");
#finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
#zone_database <- gradstein_2012_emended$zones;

#pbdb_data_list$pbdb_sites_refined[pbdb_data_list$pbdb_sites_refined$collection_no %in% c(155912,155918),]
#pbdb_data_list$pbdb_sites_refined$formation_no[pbdb_data_list$pbdb_sites_refined$collection_no %in% c(155912,155918)] <- 3461;
#pbdb_data_list$pbdb_sites_refined$rock_unit_senior[pbdb_data_list$pbdb_sites_refined$collection_no %in% c(155912,155918)] <- "Dundee Flagstone";
pbdb_data_list$pbdb_taxonomy[pbdb_data_list$pbdb_taxonomy$taxon_name=="Synphoria",]

all_trilobites <- pbdb_data_list$pbdb_taxonomy[pbdb_data_list$pbdb_taxonomy$class=="Trilobita",];
trilobite_genera <- all_trilobites[all_trilobites$accepted_rank %in% c("genus","subgenus"),];
trilobite_genera <- trilobite_genera[trilobite_genera$accepted_name==trilobite_genera$taxon_name,];
nrow(trilobite_genera)

lost_genera <- data.frame(readxl::read_xlsx("Genera_to_Add.xlsx"));
lost_genera$genus
lgenera <- nrow(lost_genera);
l_genera <- lost_genera$genus;
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
for (lg in 1:lgenera)	{
#	lost_genera$genus[lg]
	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon=l_genera[lg],pbdb_data_list,pbdb_taxonomy_edits);
	save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
	}
which(pbdb_data_list$pbdb_opinions=="Aldaniopirifer",arr.ind=T);
max(pbdb_data_list$pbdb_opinions$modified)
max(pbdb_data_list$pbdb_opinions$orig_no)
options(warn=1);



{}
