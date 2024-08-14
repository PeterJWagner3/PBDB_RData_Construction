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
do_not_forget_taxa <- c();  # must be c() for nothing; "" will be read as taxon blank.
#Sys.Date()
do_not_forget_start <- "";
do_not_forget_end <- "";
x <- c("Proterozoic","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene")

# parameters for main function ####
#do_not_forget_start <- "Neogene";
#do_not_forget_end <- "Neogene";
#date()
# Get Daily Download ####
if (do_not_forget_start!="")	{
	pbdb_data_list <- update_pbdb_RData(pbdb_data_list,rock_unit_data,gradstein_2020_emended,paleodb_fixes,do_not_forget_taxa=do_not_forget_taxa,do_not_forget_start=do_not_forget_start,do_not_forget_end=do_not_forget_end);
	} else	{
	pbdb_data_list <- update_pbdb_RData(pbdb_data_list,rock_unit_data,gradstein_2020_emended,paleodb_fixes,do_not_forget_taxa=do_not_forget_taxa);
	}

# Update Ages again (just in case) ####
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
finest_chronostrat <- pbdb_data_list$time_scale;
pbdb_sites_refined <- pbdb_data_list$pbdb_sites_refined;
age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat);
a <- match(pbdb_sites_refined$interval_lb,finest_chronostrat$interval);
z <- match(pbdb_sites_refined$interval_ub,finest_chronostrat$interval);
fix_these <- (1:nrow(pbdb_sites_refined))[a>z];
write.csv(pbdb_sites_refined[fix_these,],"effed_up_collection_dates.csv")
i <- 0;
while (i < length(fix_these))	{
	i <- i+1;
	zz <- pbdb_sites_refined$interval_lb[fix_these[i]];
	aa <- pbdb_sites_refined$interval_ub[fix_these[i]];
	pbdb_sites_refined$interval_lb[fix_these[i]] <- aa;
	pbdb_sites_refined$interval_ub[fix_these[i]] <- zz;
	pbdb_sites_refined$ma_lb[fix_these[i]] <- (finest_chronostrat$ma_lb[finest_chronostrat$interval==aa]+finest_chronostrat$ma_ub[finest_chronostrat$interval==aa])/2;
	pbdb_sites_refined$ma_ub[fix_these[i]] <- (finest_chronostrat$ma_lb[finest_chronostrat$interval==zz]+finest_chronostrat$ma_ub[finest_chronostrat$interval==zz])/2;
#	pbdb_sites_refined$ma_lb[fix_these[i]] <- mean(as.numeric(finest_chronostrat[match(pbdb_sites_refined$interval_lb[fix_these[i]],finest_chronostrat$interval),c("ma_lb","ma_ub")]));
#	pbdb_sites_refined$ma_ub[fix_these[i]] <- mean(as.numeric(finest_chronostrat[match(pbdb_sites_refined$interval_ub[fix_these[i]],finest_chronostrat$interval),c("ma_lb","ma_ub")]));
	}
age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat);
pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
save(pbdb_data_list,file=paste(data_for_r,"/Paleobiology_Database.RData",sep=""));

# Fix Accepted Names ####
print("Fixing accepted names");
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
load(paste(data_for_r,"/PaleoDB_Edits.RData",sep=""));
options(warn=1);
pbdb_finds <- pbdb_data_list$pbdb_finds;
genus_finds <- pbdb_finds[pbdb_finds$accepted_rank %in% c("genus","subgenus"),];
genus_finds <- genus_finds[genus_finds$identified_rank %in% c("species","subspecies"),];
genus_finds$identified_name <- gsub("\\)","\\) ",genus_finds$identified_name);
genus_finds$identified_name <- gsub("  "," ",genus_finds$identified_name);
print("separating unentered formal species from informal species names");
xx <- pbapply::pbsapply(genus_finds$identified_name,identify_informal_species);
genus_finds$difference[xx] <- "informal species";
informal_finds <- genus_finds[xx,];
formal_finds <- genus_finds[!xx,];
#genus_finds <- genus_finds[!pbapply::pbsapply(genus_finds$identified_name,identify_informal_species),];
#occidere_indeterminate_species(genus_finds$identified_name[68])
for (i in 1:length(taxon_qualifiers))
	formal_finds$identified_name <- gsub(taxon_qualifiers[i],"",formal_finds$identified_name);
#genus_finds$identified_name <- pbapply::pbsapply(taxon_qualifiers,gsub,"",genus_finds$identified_name)
#genus_finds[genus_finds$occurrence_no==281754,]
#genus_finds$identified_name[1:100]
taxon_name <- formal_finds$identified_name;
species_epithet <- pbapply::pbsapply(taxon_name,divido_species_epithets);

good_epithets <- gsub("\\.","",species_epithet)==species_epithet;
bad_epithets <- gsub("\\.","",species_epithet)!=species_epithet;
formal_finds$difference[bad_epithets] <- "particular informal species";
#species_epithet <- species_epithet[good_epithets];
formal_finds$accepted_name <- paste(formal_finds$genus,species_epithet);
#pbdb_finds[pbdb_finds$occurrence_no==formal_finds$occurrence_no[967],];
formal_finds$accepted_rank <- "species";
#genus_finds <- rbind(informal_finds,formal_finds);
#genus_finds <- genus_finds[order(genus_finds$collection_no,genus_finds$occurrence_no),];
pbdb_finds <- pbdb_finds[!pbdb_finds$occurrence_no %in% informal_finds$occurrence_no,];
pbdb_finds <- pbdb_finds[!pbdb_finds$occurrence_no %in% formal_finds$occurrence_no,];
pbdb_finds <- rbind(pbdb_finds,informal_finds,formal_finds);
pbdb_finds <- pbdb_finds[order(pbdb_finds$collection_no,pbdb_finds$occurrence_no),];
#nrow(pbdb_finds)
#good_epithets <- gsub("\\.","",species_epithet)==species_epithet;
#good_epithets[pbdb_finds$accepted_rank %in% c("species","subspecies")] <- F;
#genus_finds$accepted_name[good_epithets] <- paste(genus_finds$genus[good_epithets],species_epithet[good_epithets]);
#pbdb_finds$accepted_name[good_epithets] <- "species";
#cbind(pbdb_finds$accepted_rank[pbdb_finds$genus %in% "Pseudobakewellia"],pbdb_finds$accepted_name[pbdb_finds$genus %in% "Pseudobakewellia"]);
pbdb_data_list$pbdb_finds <- pbdb_finds;
save(pbdb_data_list,file=paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
#load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
pbdb_data_list_smol$pbdb_finds <- pbdb_finds;
save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));

# Update Paleogeography ####
geog_models <- c("GOLONKA","MULLER2022","MERDITH2021","MULLER2019","MULLER2016","MATTHEWS2016_mantle_ref","MATTHEWS2016_pmag_ref","RODINIA2013","SETON2012","PALEOMAP");
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
load(paste(data_for_r,"/Paleobiology_Database_Extended_Paleogeographic_Data.RData",sep=""));
options(warn=1);

pbdb_sites_refined <- pbdb_data_list$pbdb_sites_refined;
print(paste("Getting paleogeography models from PBDB Information"));
pbdb_download <- accersi_paleogeographic_data_for_all_collections();
p_s_r <- (1:nrow(pbdb_sites_refined))[pbdb_sites_refined$collection_no %in% pbdb_download$collection_no];
np_s_r <- (1:nrow(pbdb_sites_refined))[!pbdb_sites_refined$collection_no %in% pbdb_download$collection_no];
coll_id <- pbdb_sites_refined$collection_no[np_s_r];
a <- 1;
while (a <= length(np_s_r))	{
	pbdb_download <- rbind(pbdb_download,accersi_paleogeographic_data_for_one_collection_no(coll_id[a]));
	a <- a+1;
	}
#pbdb_sites_refined$collection_no[is.na(match(pbdb_sites_refined$collection_no,pbdb_download$collection_no))]
pbdb_download <- unique(pbdb_download);
pbdb_download <- pbdb_download[order(pbdb_download$collection_no),];
#pbdb_sites_refined$geoplate[pbdb_sites_refined$collection_no %in% pbdb_download$collection_no[pbdb_download$formation %in% c("Bostow","Bostów","Bóstow")]] <- 305;
#pbdb_sites_refined[pbdb_sites_refined$rock_unit_senior=="Bostow",]
pbdb_sites_refined <- pbdb_sites_refined[order(pbdb_sites_refined$collection_no),];
pbdb_sites_refined$geoplate <- as.numeric(pbdb_download$geoplate[match(pbdb_sites_refined$collection_no,pbdb_download$collection_no)]);

#pbdb_sites_refined$geoplate[pbdb_sites_refined$rock_unit_senior=="Bostow"] <- 305;
pbdb_sites_refined$paleolng[pbdb_sites_refined$collection_no %in% pbdb_download$collection_no[pbdb_download$geoplate==302 & pbdb_download$geoplate2==305]] <- pbdb_download$paleolat2[pbdb_download$geoplate==302 & pbdb_download$geoplate2==305];
pbdb_sites_refined$paleolng[pbdb_sites_refined$collection_no %in% pbdb_download$collection_no[pbdb_download$geoplate==302 & pbdb_download$geoplate2==305]] <- pbdb_download$paleolng2[pbdb_download$geoplate==302 & pbdb_download$geoplate2==305];
pbdb_sites_refined$paleomodel[pbdb_sites_refined$collection_no %in% pbdb_download$collection_no[pbdb_download$geoplate==302 & pbdb_download$geoplate2==305]] <- "scotese";
pbdb_sites_refined$geoplate[pbdb_sites_refined$collection_no %in% pbdb_download$collection_no[pbdb_download$geoplate==302 & pbdb_download$geoplate2==305]] <- 305;
sum(is.na(match(pbdb_sites_refined$collection_no,pbdb_download$collection_no)))
pbdb_sites_refined$geoplate <- as.numeric(pbdb_sites_refined$geoplate);
pbdb_sites_refined$geoplate <- as.numeric(pbdb_download$geoplate[match(pbdb_sites_refined$collection_no,pbdb_download$collection_no)]);
pbdb_sites_refined$paleolat <- as.numeric(pbdb_download$paleolat[match(pbdb_sites_refined$collection_no,pbdb_download$collection_no)]);
pbdb_sites_refined$paleolng <- as.numeric(pbdb_download$paleolng[match(pbdb_sites_refined$collection_no,pbdb_download$collection_no)]);
pbdb_sites_refined$created <- pbdb_download$created[match(pbdb_sites_refined$collection_no,pbdb_download$collection_no)];
pbdb_sites_refined$modified <- pbdb_download$modified[match(pbdb_sites_refined$collection_no,pbdb_download$collection_no)];

pbdb_sites <- pbdb_data_list$pbdb_sites;
pbdb_sites <- pbdb_sites[order(pbdb_sites$collection_no),];
a_in_b <- (1:nrow(pbdb_sites))[pbdb_sites$collection_no %in% pbdb_download$collection_no];
b_in_a <- (1:nrow(pbdb_download))[pbdb_download$collection_no %in% pbdb_sites$collection_no];
pbdb_sites$geoplate[a_in_b] <- as.numeric(pbdb_download$geoplate[b_in_a]);
fix_me <- pbdb_sites$collection_no[a_in_b][is.na(pbdb_sites$geoplate[a_in_b])];
a_in_b <- match(fix_me,pbdb_sites$collection_no);
b_in_a <- match(fix_me,pbdb_download$collection_no);
pbdb_sites$geoplate[a_in_b] <- as.numeric(pbdb_download$geoplate2[b_in_a]);
fix_me <- pbdb_sites$collection_no[a_in_b][is.na(pbdb_sites$geoplate[a_in_b])];
a_in_b <- match(fix_me,pbdb_sites$collection_no);
b_in_a <- match(fix_me,pbdb_download$collection_no);
pbdb_sites$geoplate[a_in_b] <- as.numeric(pbdb_download$geoplate3[b_in_a]);
pbdb_data_list$pbdb_sites <- pbdb_sites;
#hist(pbdb_sites_refined$geoplate,breaks=c(0,unique(sort(pbdb_sites_refined$geoplate))))
#added_paleogeography <- read.csv("Added_Paleogeography.csv");
mcoll_id <- pbdb_sites_refined$collection_no[!pbdb_sites_refined$collection_no %in% added_paleogeography$collection_no];
mcoll <- length(mcoll_id);

if (mcoll>0)	{
	added <- added_paleogeography[1:mcoll,];
	added$collection_no <- mcoll_id;
	old_sites <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% mcoll_id,];
	for (m in 1:3)	{
		old_sites <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% mcoll_id,];
		paleomodel <- toupper(added_paleogeography[,paste("paleomodel",m,sep="")][1]);
		print(paste("Getting paleogeographic estimates from",paleomodel,"model"));
		if (paleomodel=="GOLONKA")
			old_sites <- old_sites[((old_sites$ma_lb+old_sites$ma_ub)/2)<540,];
		if (paleomodel %in% geog_models)	{
			paleolat <- paste("paleolat",m,sep="");
			paleolng <- paste("paleolng",m,sep="");
			added[,c(paleolat,paleolng)] <- 0.0;
			old_sites <- mundify_paleolatitude_and_paleolongitude(old_sites=old_sites,model=paleomodel);
			if (nrow(old_sites)>0)	{
				old_sites <- old_sites[order(old_sites$collection_no),];
				added[added$collection_no %in% old_sites$collection_no,c(paleolat,paleolng)] <- old_sites[,c("paleolat","paleolng")];
				}
			}
		}

	colnames(pbdb_download)[colnames(pbdb_download)=="paleomodel"] <- "paleomodel1";
	colnames(pbdb_download)[colnames(pbdb_download)=="paleolat"] <- "paleolat1";
	colnames(pbdb_download)[colnames(pbdb_download)=="paleolng"] <- "paleolng1";
	for (m in 4:6)	{
		paleolat <- paste("paleolat",m-3,sep="");
		paleolng <- paste("paleolng",m-3,sep="");
		added <- added[order(added$collection_no),];
		old_sites <- old_sites[order(old_sites$collection_no),];
		added[,c(paleolng,paleolat)] <- 0;

		relv_coll <- added$collection_no[added$collection_no %in% old_sites$collection_no];
		coll_relv <- old_sites$collection_no[old_sites$collection_no %in% added$collection_no]
		r_c <- (1:nrow(added))[added$collection_no %in% relv_coll];
		c_r <- (1:nrow(old_sites))[old_sites$collection_no %in% coll_relv];
		added[r_c,paleolng] <- pbdb_download[match(old_sites$collection_no[c_r],pbdb_download$collection_no),paleolng];
		added[r_c,paleolat] <- pbdb_download[match(old_sites$collection_no[c_r],pbdb_download$collection_no),paleolat];
		}
	#added_paleogeography$modifed <- Sys.Date()-23;
	added$modified <- Sys.Date();
	colnames(added) <- colnames(added_paleogeography);
	added_paleogeography <- rbind(added_paleogeography,added);
	}
#ncol(added_paleogeography)
#ncol(added)
added_paleogeography <- added_paleogeography[order(added_paleogeography$collection_no),];
write.csv(added_paleogeography,"Added_Paleogeography.csv",row.names=FALSE);
save(added_paleogeography,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database_Extended_Paleogeographic_Data.RData",sep=""));
#pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

nsites <- nrow(pbdb_sites_refined);
plateless <- sort(unique(c((1:nsites)[is.na(pbdb_sites_refined$geoplate)],(1:nsites)[!pbdb_sites_refined$geoplate %in% 100:1000])));
length(plateless)
pbdb_sites_refined$geoplate[plateless] <- as.numeric(pbdb_download$geoplate[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
if (is.null(pbdb_download$paleolat))	{
	pbdb_sites_refined$paleolat[plateless] <- as.numeric(pbdb_download$paleolat1[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
	pbdb_sites_refined$paleolng[plateless] <- as.numeric(pbdb_download$paleolng1[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
	} else	{
	pbdb_sites_refined$paleolat[plateless] <- as.numeric(pbdb_download$paleolat[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
	pbdb_sites_refined$paleolng[plateless] <- as.numeric(pbdb_download$paleolng[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
	}
plateless <- (1:nsites)[is.na(pbdb_sites_refined$geoplate)];
length(plateless)
pbdb_sites_refined$geoplate[plateless] <- as.numeric(pbdb_download$geoplate2[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
pbdb_sites_refined$paleolat[plateless] <- as.numeric(pbdb_download$paleolat2[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
pbdb_sites_refined$paleolng[plateless] <- as.numeric(pbdb_download$paleolng2[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
plateless <- (1:nsites)[is.na(pbdb_sites_refined$geoplate)];
length(plateless)
pbdb_sites_refined$geoplate[plateless] <- as.numeric(pbdb_download$geoplate3[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
pbdb_sites_refined$paleolat[plateless] <- as.numeric(pbdb_download$paleolat3[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
pbdb_sites_refined$paleolng[plateless] <- as.numeric(pbdb_download$paleolng3[match(pbdb_sites_refined$collection_no[plateless],pbdb_download$collection_no)]);
plateless <- (1:nsites)[is.na(pbdb_sites_refined$geoplate)];
length(plateless);
#pbdb_sites_refined$collection_no[plateless];
#pbdb_data_list$pbdb_sites_refined[pbdb_data_list$pbdb_sites_refined$collection_no %in% 9992,]
pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;

pbdb_finds <- pbdb_data_list$pbdb_finds;
nfinds <- nrow(pbdb_finds);
gfinds <- (1:nfinds)[pbdb_finds$collection_no %in% pbdb_download$collection_no]
pbdb_finds$max_ma[gfinds] <- pbdb_download$max_ma[match(pbdb_finds$collection_no[gfinds],pbdb_download$collection_no)];
pbdb_finds$min_ma[gfinds] <- pbdb_download$min_ma[match(pbdb_finds$collection_no[gfinds],pbdb_download$collection_no)];
pbdb_data_list$pbdb_finds <- pbdb_finds;
#pbdb_sites_refined$rock_no[pbdb_sites_refined$collection_no==231439] <- pbdb_sites_refined$rock_no_sr[pbdb_sites_refined$collection_no==231439] <- 14215;
#pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

#load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
pbdb_data_list_smol$pbdb_sites <- pbdb_data_list$pbdb_sites_refined;
pbdb_data_list_smol$pbdb_finds <- pbdb_data_list$pbdb_finds;
save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));

#max(pbdb_sites_refined$created)
beepr::beep("wilhelm");

# Update Taxonomy ####
higher_taxa_to_update <- c("Ichnofossils","Problematica","Cyanobacteria","Retaria","Rhodophyta",
						   "Arboreomorpha","Rangeomorpha","Erniettomorpha","Tetraradialomorpha","Bilaterialomorpha",
						   "Porifera","Cnidaria","Lophophorata",
						   "Hemichordata","Echinodermata","Mollusca","Panarthropoda","Chordata");
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
pbdb_taxonomy_edits <- paleodb_fixes$paleodb_taxonomy_edits;
#length(pbdb_taxonomy_edits$taxon_no); length(unique(pbdb_taxonomy_edits$taxon_no))
load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
options(warn=1);
ht <- 1;
while (ht <= length(higher_taxa_to_update))	{
	#while (ht <= match("Bilaterialomorpha",higher_taxa_to_update))	{
#for (ht in 16:16)	{
#	load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
#	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon="Chordata",pbdb_data_list,pbdb_taxonomy_edits);
#	pbdb_data_list$pbdb_finds$identified_name[pbdb_data_list$pbdb_finds$occurrence_no %in% c(1592363,1658884)] <- "Lancastria cf. plana Geyer and Peel 2011";
	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon=higher_taxa_to_update[ht],pbdb_data_list,pbdb_taxonomy_edits);
#	pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$occurrence_no %in% c(1592363,1658884),]
#	save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
#	pbdb_data_list$pbdb_taxonomy[pbdb_data_list$pbdb_taxonomy$taxon_name=="Erniettamorpha",]
#	pbdb_data_list <- update_taxonomic_data_in_pbdb_RData(taxon="Marrellomorpha",pbdb_data_list);
	pbdb_data_list_smol$pbdb_finds <- pbdb_data_list$pbdb_finds;
	pbdb_data_list_smol$pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
	save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));
	ht <- ht+1;
	}
#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
#pbdb_data_list$pbdb_taxonomy[nrow(pbdb_data_list$pbdb_taxonomy),];
#write.csv(pbdb_taxonomy[pbdb_taxonomy$taxon_name %in% c("Retillamina","Retilamina"),],"Retilamina.csv",row.names = F,fileEncoding = "UTF-8");
# pick up here if there was a timeout while updating taxonomy
options(warn=0);
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
options(warn=1);
pbdb_finds <- pbdb_data_list$pbdb_finds;
pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
pbdb_finds$genus[pbdb_finds$subgenus_no>0] <- pbdb_taxonomy$accepted_name[match(pbdb_finds$subgenus_no[pbdb_finds$subgenus_no>0],pbdb_taxonomy$taxon_no)];
pbdb_finds$accepted_name <- pbapply::pbsapply(pbdb_finds$accepted_name,mundify_taxon_names);
pbdb_data_list$pbdb_finds <- pbdb_finds;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

#load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
pbdb_data_list_smol$pbdb_finds <- pbdb_data_list$pbdb_finds;
pbdb_data_list_smol$pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));

beepr::beep("wilhelm");

{}
