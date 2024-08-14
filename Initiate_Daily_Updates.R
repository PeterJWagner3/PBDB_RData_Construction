#install.packages("taskscheduleR")
require(cronR);
# list the contenst of a crontab
cron_ls();
# list the full path of the rscript
path <- "/Users/peterjwagner/Documents/R_Projects/PBDB_RData_Construction/Daily_PBDB_Update.r";
# create the command to execute it
cmd <- cron_rscript(path);
# add the command & specify when
cron_add(command=cmd,frequency="daily",at="1:30 AM");

zeit <- format(as.POSIXct(Sys.time()), format = "%H:%M");
undone <- T;
counter <- 0;
path <- "/Users/peterjwagner/Documents/R_Projects/PBDB_RData_Construction/Daily_PBDB_Update.r";
while (!zeit %in% c("01:30","1:30"))	{
	counter <- counter+1;
	alt_zeit <- zeit;
	zeit <- format(as.POSIXct(Sys.time()), format = "%H:%M");
	if (zeit!=alt_zeit)	{
		print(zeit);
		if (alt_zeit<"01:30" && zeit>"01:30")	zeit <- "01:30";
		}
	}
path <- "/Users/peterjwagner/Documents/R_Projects/PBDB_RData_Construction/Daily_PBDB_Update.r";
source(path);

#load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));

write.csv(pbdb_taxonomy[pbdb_taxonomy$taxon_name %in% pbdb_taxonomy$accepted_name[pbdb_taxonomy$taxon_name %in% c("Pseudocryphaeus quarterspinosus","Minicryphaeus quaterspinosus")],],"Minicryphaeus_quaterspinosus.csv",row.names = FALSE)
pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$occurrence_no %in% c(1592363,1658884),]

accersi_all_occurrences_data <- function()	{
#paste("https://paleobiodb.org/data1.2/occs/single.txt?id=",occ_no,"&show=refattr,classext,rem,entname,abund,ecospace,crmod&limit=all",sep="");
#finds <- read.csv(http,header=TRUE,stringsAsFactors=FALSE,encoding="UTF-8");
http <- "https://paleobiodb.org/data1.2/occs/list.csv?base_name=Life,Ichnofossils";
finds <- read.csv(http,header=TRUE,stringsAsFactors=FALSE,encoding="UTF-8");
}

all_finds <- read.csv(file.choose(),header=TRUE,stringsAsFactors=FALSE,encoding="UTF-8");
pbdb_finds <- pbdb_data_list$pbdb_finds;
missing_finds <- all_finds$occurrence_no[!all_finds$occurrence_no %in% pbdb_finds$occurrence_no];
lost_finds <- pbdb_finds[1,];
lost_finds <- lost_finds[lost_finds$occurrence_no<1,];
#nn <- pbapply::pbsapply(missing_finds,get_occurrences_from_occurrence_number);
mfinds <- length(missing_finds);
triggers <- ceiling(((1:199)*mfinds)/200);
for (i in 1:mfinds)	{
	if (i %in% triggers)	print(paste(100*match(i,triggers)/(length(triggers)+1),"% done",sep=""))
	xx <- get_occurrences_from_occurrence_number(missing_finds[i])
	xx$accepted_name_orig <- xx$accepted_name;
	lost_finds <- rbind(lost_finds,xx)
#	xx$identified_name <- mundify_taxon_names(xx$identified_name)
	}
pbdb_finds <- rbind(pbdb_finds,lost_finds);
pbdb_finds <- pbdb_finds[order(pbdb_finds$occurrence_no),];
pbdb_data_list$pbdb_finds <- pbdb_finds;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

mundify_taxon_names("Lancastria cf. plana Geyer and Peel 2011 informal")
nn <- update_taxonomy("Macrobiotus")
nn$phylum <- "Tardigrada";
nn$is_extant <- "extant";
nn$phylum_no <-67137;
pbdb_data_list$pbdb_taxonomy <- rbind(pbdb_data_list$pbdb_taxonomy,nn);
pbdb_data_list$pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy[order(pbdb_data_list$pbdb_taxonomy$taxon_no),];
pbdb_data_list$pbdb_sites_refined$ma_lb[pbdb_data_list$pbdb_sites_refined$collection_no==159584] <- 516.3;
pbdb_data_list$pbdb_sites_refined$ma_ub[pbdb_data_list$pbdb_sites_refined$collection_no==159584] <- 515.5;
pbdb_data_list$pbdb_sites_refined$interval_lb[pbdb_data_list$pbdb_sites_refined$collection_no==159584] <- "Cm32";
pbdb_data_list$pbdb_sites_refined$interval_ub[pbdb_data_list$pbdb_sites_refined$collection_no==159584] <- "Cm33";
#write.csv(pbdb_taxonomy[pbdb_taxonomy$taxon_name %in% c("Gamadiscus","Gammadiscus","Cyrtodiscus"),],"Gamadiscus.csv",row.names=FALSE,fileEncoding = "UTF-8");
#pbdb_sites <- pbdb_data_list$pbdb_sites_refined;
#pbdb_sites_orig <- pbdb_data_list$pbdb_sites;
#nsites <- nrow(pbdb_sites);
#effed_sites <- (1:nsites)[pbdb_sites$geoplate!=floor(pbdb_sites$geoplate)];
#ffd_sites <- (1:nsites)[pbdb_sites_orig$geoplate!=floor(pbdb_sites_orig$geoplate)];

#fix_these <- pbdb_sites_orig[ffd_sites,];
#fix_these <- fix_these[!is.na(fix_these$lat),];
#fix_these <- mundify_paleolatitude_and_paleolongitude(old_sites=fix_these);

#pbdb_data_list$pbdb_sites$paleolat[ffd_sites] <- fix_these$paleolat;
#pbdb_data_list$pbdb_sites$paleolng[ffd_sites] <- fix_these$paleolng;
#pbdb_data_list$pbdb_sites_refined$paleolat[ffd_sites] <- fix_these$paleolat;
#pbdb_data_list$pbdb_sites_refined$paleolng[ffd_sites] <- fix_these$paleolng;
#for (i in 1:length(ffd_sites))	{
#	pbdb_data_list$pbdb_sites$paleolat[ffd_sites[i]] <- as.numeric(fix_these$paleolat[i]);
#	pbdb_data_list$pbdb_sites$paleolng[ffd_sites[i]] <- as.numeric(fix_these$paleolng[i]);
#	pbdb_data_list$pbdb_sites_refined$paleolat[ffd_sites[i]] <- as.numeric(fix_these$paleolat[i]);
#	pbdb_data_list$pbdb_sites_refined$paleolng[ffd_sites[i]] <- as.numeric(fix_these$paleolng[i]);
#	}
#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

#for (i in 1:nrow(fix_these))	{
#	if (i %in% ceiling(0.01*nrow(fix_these)*1:99))
#		print(paste(match(i,ceiling(0.01*nrow(fix_these)*1:99)),"% done",sep=""));
#	fix_these$geoplate[i] <- accersi_geoplate_data(as.integer(fix_these$collection_no[i]));
#	}
#fix_these$geoplate <- as.numeric(fix_these$geoplate);
#for (i in 1:length(ffd_sites))	{
#	pbdb_data_list$pbdb_sites$geoplate[ffd_sites[i]] <- fix_these$geoplate[i];
#	pbdb_data_list$pbdb_sites_refined$geoplate[ffd_sites[i]] <- fix_these$geoplate[i];
#	}
#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

#pbdb_data_list$pbdb_sites$geoplate[is.na(pbdb_data_list$pbdb_sites$geoplate)] <- 0;
#pbdb_data_list$pbdb_sites_refined$geoplate[is.na(pbdb_data_list$pbdb_sites_refined$geoplate)] <- 0;
#still_effed <- (1:nsites)[pbdb_data_list$pbdb_sites_refined$geoplate<100];
#fix_these2 <- pbdb_data_list$pbdb_sites_refined[still_effed,];
#fix_these2 <- fix_these2[!is.na(fix_these2$lat),];
#good_still <- pbdb_data_list$pbdb_sites_refined[!(1:nsites) %in% still_effed,];
#for (i in 1:nrow(fix_these2))	{
#	j <- match(fix_these2$collection_no[i],pbdb_data_list$pbdb_sites_refined$collection_no);
#	k <- match(fix_these2$collection_no[i],pbdb_data_list$pbdb_sites$collection_no);
##	pbdb_data_list$pbdb_sites_refined[j,]
#	if (i %in% ceiling(0.01*nrow(fix_these2)*1:99))
#		print(paste(match(i,ceiling(0.01*nrow(fix_these2)*1:99)),"% done",sep=""));
#	dlat <- abs(good_still$lat-fix_these2$lat[i]);
#	dlng <- abs(good_still$lng-fix_these2$lng[i]);
#	crude_dist <- (dlat^2+dlng^2)^0.5;
#	fix_these2$geoplate[i] <- good_still$geoplate[match(min(crude_dist),crude_dist)]
#	pbdb_data_list$pbdb_sites_refined$geoplate[j] <- fix_these2$geoplate[i];
#	pbdb_data_list$pbdb_sites$geoplate[k] <- fix_these2$geoplate[i];
##	accersi_geoplate_data(coll_id=as.integer(fix_these2$collection_no[i]));
#	}
#fix_these2 <- mundify_paleolatitude_and_paleolongitude(old_sites=fix_these2);

#pbdb_data_list$pbdb_sites$paleolat[match(fix_these2$collection_no,pbdb_data_list$pbdb_sites$collection_no)] <- fix_these2$paleolat;
#pbdb_data_list$pbdb_sites$paleolng[match(fix_these2$collection_no,pbdb_data_list$pbdb_sites$collection_no)] <- fix_these2$paleolng;
#pbdb_data_list$pbdb_sites$geoplate[match(fix_these2$collection_no,pbdb_data_list$pbdb_sites$collection_no)] <- fix_these2$geoplate;
for (i in 1:nrow(fix_these2))	{
	j <- match(fix_these2$collection_no[i],pbdb_data_list$pbdb_sites_refined$collection_no);
	pbdb_data_list$pbdb_sites_refined$paleolat[j] <- fix_these2$paleolat[i];
	pbdb_data_list$pbdb_sites_refined$paleolng[j] <- fix_these2$paleolng[i];
#	pbdb_data_list$pbdb_sites_refined$geoplate[j] <- fix_these2$geoplate[i];
	}
#pbdb_data_list$pbdb_sites_refined$paleolat[match(fix_these2$collection_no,pbdb_data_list$pbdb_sites_refined$collection_no)] <- fix_these2$paleolat;
#pbdb_data_list$pbdb_sites_refined$paleolng[match(fix_these2$collection_no,pbdb_data_list$pbdb_sites_refined$collection_no)] <- fix_these2$paleolng;
#pbdb_data_list$pbdb_sites_refined$geoplate[match(fix_these2$collection_no,pbdb_data_list$pbdb_sites_refined$collection_no)] <- fix_these2$geoplate;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));

coll_id <- fix_these$collection_no;
dummys <- pbapply::pbsapply(coll_id,accersi_geoplate_data);
pbdb_data_list$pbdb_sites$geoplate[ffd_sites] <- dummys$geoplate;
pbdb_data_list$pbdb_sites_refined$geoplate[ffd_sites] <- dummys$geoplate;

dummy <- pbapply::pbsapply(coll_id,accersi_data_for_one_collection);

write.csv(pbdb_taxonomy[pbdb_taxonomy$accepted_name %in% c("Glyptosphaerites","Sphaeronites","Sphaeronites (Sphaeronites)","Sphaeronites pomum","Sphaeronites (Sphaeronites) pomum"),],"Sphaeronites_pomum.csv",row.names=F)

write.csv(pbdb_taxonomy[pbdb_taxonomy$taxon_name=="Sphaeronites leuchtenbergi",],"leuchtenbergi.csv",row.names = F,fileEncoding = "UTF8");
pbdb_taxonomy[pbdb_taxonomy$accepted_name %in% "Glyptosphaerites leuchtenbergi",];
pbdb_taxonomy[pbdb_taxonomy$accepted_name %in% "Sphaeronites pomum",];
pbdb_finds[pbdb_finds$accepted_name %in% "Sphaeronites pomum",];
pbdb_finds[pbdb_finds$accepted_name %in% "Glyptosphaerites leuchtenbergi",];

load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
load(paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
pbdb_finds <- pbdb_data_list$pbdb_finds;
pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;

gen_data <- (1:nrow(pbdb_taxonomy))[pbdb_taxonomy$genus_no>0];
pbdb_taxonomy$family_no <- as.numeric(pbdb_taxonomy$family_no);
pbdb_taxonomy$family[gen_data] <- pbdb_taxonomy$family[match(as.numeric(pbdb_taxonomy$genus_no[gen_data]),pbdb_taxonomy$taxon_no)];
pbdb_taxonomy$family_no[gen_data] <- as.numeric(pbdb_taxonomy$family_no[match(as.numeric(pbdb_taxonomy$genus_no[gen_data]),pbdb_taxonomy$taxon_no)]);

pbdb_taxonomy$family_no[is.na(pbdb_taxonomy$family_no)] <- 0;
fam_data <- (1:nrow(pbdb_taxonomy))[pbdb_taxonomy$family_no>0];
fam_data <- fam_data[!is.na(match(as.numeric(pbdb_taxonomy$family_no[fam_data]),pbdb_taxonomy$taxon_no))];
pbdb_taxonomy$class[fam_data] <- pbdb_taxonomy$class[match(as.numeric(pbdb_taxonomy$family_no[fam_data]),pbdb_taxonomy$taxon_no)];
pbdb_taxonomy$class_no[fam_data] <- as.numeric(pbdb_taxonomy$class_no[match(as.numeric(pbdb_taxonomy$family_no[fam_data]),pbdb_taxonomy$taxon_no)]);
pbdb_taxonomy$order[fam_data] <- pbdb_taxonomy$order[match(as.numeric(pbdb_taxonomy$family_no[fam_data]),pbdb_taxonomy$taxon_no)];
pbdb_taxonomy$order_no[fam_data] <- as.numeric(pbdb_taxonomy$order_no[match(as.numeric(pbdb_taxonomy$family_no[fam_data]),pbdb_taxonomy$taxon_no)]);
pbdb_taxonomy$family[fam_data] <- pbdb_taxonomy$family[match(as.numeric(pbdb_taxonomy$family_no[fam_data]),pbdb_taxonomy$taxon_no)];
pbdb_taxonomy$family_no[fam_data] <- as.numeric(pbdb_taxonomy$family_no[match(as.numeric(pbdb_taxonomy$family_no[fam_data]),pbdb_taxonomy$taxon_no)]);

pbdb_finds$genus[pbdb_finds$subgenus_no>0] <- pbdb_taxonomy$accepted_name[match(pbdb_finds$subgenus_no[pbdb_finds$subgenus_no>0],pbdb_taxonomy$taxon_no)];
gfinds <- (1:nrow(pbdb_finds))[pbdb_finds$genus_no>0];
pbdb_finds$genus[gfinds] <- pbdb_taxonomy$accepted_name[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$genus_no[gfinds] <- pbdb_taxonomy$accepted_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$family[gfinds] <- pbdb_taxonomy$family[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$family_no[gfinds] <- pbdb_taxonomy$family_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$order[gfinds] <- pbdb_taxonomy$order[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$order_no[gfinds] <- pbdb_taxonomy$order_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];

pbdb_data_list$pbdb_finds <- pbdb_finds;
pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
pbdb_data_list_smol$pbdb_finds <- pbdb_data_list$pbdb_finds;
pbdb_data_list_smol$pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
save(pbdb_data_list_smol,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData",sep=""));

higher_taxa_to_update <- c("Ichnofossils","Problematica","Cyanobacteria","Retaria","Rhodophyta",
						   "Arboreomorpha","Rangeomorpha","Erniettomorpha","Tetraradialomorpha","Bilaterialomorpha",
						   "Porifera","Cnidaria","Lophophorata","Mollusca",
						   "Panarthropoda","Chordata","Hemichordata","Echinodermata");
load(paste(data_for_r,"/Paleobiology_Database.RData",sep=""));
for (ht in 11:length(higher_taxa_to_update))	{
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


while (undone)	{
	alt_zeit <- zeit;
	zeit <- format(as.POSIXct(Sys.time()), format = "%H:%M");
	counter <- counter+1;
	if (zeit>alt_zeit)	print(zeit);
	if (zeit %in% c("02:30","2:30"))	{
		source(path);
		undone <- F;
		}
	}

pbdb_finds <- pbdb_data_list$pbdb_finds
pbdb_finds[pbdb_finds$identified_name %in% "Xenodiscus whiteanus",]
pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy
pbdb_taxonomy[pbdb_taxonomy$taxon_name=="Xenodiscus whiteanus",]

{}
