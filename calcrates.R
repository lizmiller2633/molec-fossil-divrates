# Methods and Functions are camelCase. Variables and Data Structures are PascalCase
# Fields generally follow snake_case for better SQL compatibility
# Dependency functions are not embedded in master functions
# []-notation is used wherever possible, and $-notation is avoided.
# []-notation is slower, but more explicit and works for atomic vectors

#############################################################################################################
############################################## CONFIGURATION, SCRIPT ########################################
#############################################################################################################    
# Increase the timeout time and change the fancyquote settings
options(timeout=600, "useFancyQuotes"=FALSE)

# Load or install sf package
if (suppressWarnings(require("velociraptr"))==FALSE) {
        install.packages("velociraptr",repos="http://cran.cnr.berkeley.edu/");
        library("velociraptr");
        }

#############################################################################################################
########################################### FUNCTIONS, mole-fossil-rates ####################################
#############################################################################################################
# Assign the Footeian status of in-bin occurrences	
assignStatus<-function(Data,Intervals=Epochs) {
	AgeMatrix<-t(velociraptr::presenceMatrix(Data,"early_interval"))
	# Subset the data to only include intervals in both arguments
	Intervals<-subset(Intervals,Intervals[,"name"]%in%colnames(AgeMatrix)==TRUE)
	AgeMatrix<-AgeMatrix[,as.character(Intervals[,"name"])]
	for (i in 1:nrow(AgeMatrix)) {
		AgeMatrix[i,which(AgeMatrix[i,]==1)[1]]<-2
		AgeMatrix[i,rev(which(AgeMatrix[i,]>0))[1]]<-3
		if (length(which(AgeMatrix[i,]==2))==0) {
			AgeMatrix[i,which(AgeMatrix[i,]==3)]<-4
			}
		}
	return(AgeMatrix)
	}
	
# Calculat p and q
calcRates<-function(AgeMatrix,Rate="Origination") {
	Rates<-switch(Rate,
		"Origination"=apply(AgeMatrix,2,function(x) length(which(x==1))/length(which(x==3 | x==1))),
		"Extinction"=apply(AgeMatrix,2,function(x) length(which(x==1))/length(which(x==2 | x==1)))
		)
	return(-log(Rates))
	}
########################################### SCRIPT, mole-fossil-rates #######################################
# Download from the API
CanonicalTaxa<-c("Bivalvia","Gastropoda","Anthozoa","Brachiopoda","Trilobita","Bryozoa","Nautiloidea","Ammonoidea","Crinoidea","Blastoidea","Edrioasteroidea")
CanonicalPBDB<-downloadPBDB(CanonicalTaxa,"Cambrian","Pleistocene")

# Clean and simplify the data
CanonicalPBDB<-velociraptr::cleanRank(CanonicalPBDB,"genus")
# Remove blank geoplate occurrences
CanonicalPBDB<-subset(CanonicalPBDB,is.na(CanonicalPBDB[,"geoplate"])!=TRUE)
# Remove unnecessary columns for increased performance
CanonicalPBDB<-CanonicalPBDB[,c("genus","family","class","early_interval","late_interval","max_ma","min_ma","paleolat","paleolng","geoplate")]

# Download ages timescale
Ages<-velociraptr::downloadTime("international%20ages")

# Sort and constrain ages and epochs
CanonicalStages<-velociraptr::constrainAges(CanonicalPBDB,Ages)
# Remove intervals with too few occurrences
CanonicalStages<-subset(CanonicalStages,CanonicalStages[,"early_interval"]%in%names(which(table(CanonicalStages[,"early_interval"])>=1000))==TRUE)

# For each stage, calculate whether a taxon is within-bin, through-bin, start-bin, or end-bin category of Foote
BinStatus<-assignStatus(CanonicalStages,Ages)

# Calculate the origination and extinction rate, respecitvely
GlobalQ<-calcRates(BinStatus,"Extinction")
GlobalP<-calcRates(BinStatus,"Origination")

# Clean up NaN and Inf for youngest and oldest interval
GlobalQ[is.infinite(GlobalQ) | is.nan(GlobalQ)] <- NA
GlobalP[is.infinite(GlobalP) | is.nan(GlobalP)] <-NA

# Merge Rates with Ages
GlobalQ<-transform(merge(as.data.frame(GlobalQ),Ages,by="row.names",all=FALSE),row.names=Row.names,Row.names=NULL)
GlobalP<-transform(merge(as.data.frame(GlobalP),Ages,by="row.names",all=FALSE),row.names=Row.names,Row.names=NULL)

# Convert to Rates (i.e., by time)
GlobalQ[,"GlobalQ"]<-GlobalQ[,"GlobalQ"]/(GlobalQ[,"b_age"]-GlobalQ[,"t_age"])
GlobalP[,"GlobalP"]<-GlobalP[,"GlobalP"]/(GlobalP[,"b_age"]-GlobalP[,"t_age"])

# Reorder the dataset
GlobalQ<-GlobalQ[order(GlobalQ[,"Midpoint"]),]
GlobalP<-GlobalP[order(GlobalP[,"Midpoint"]),]
