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
############################################## FUNCTIONS, FOSSIL-RATES ######################################
#############################################################################################################
# Assign the Footeian status of in-bin occurrences	
assignStatus<-function(AgeMatrix) {
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
		    
############################################### SCRIPT, FOSSIL-RATES ########################################
# Download fossil data from the pbdb API
CanonicalTaxa<-c("Bivalvia","Gastropoda","Anthozoa","Brachiopoda","Trilobita","Bryozoa","Nautiloidea","Ammonoidea","Crinoidea","Blastoidea","Edrioasteroidea")
CanonicalPBDB<-velociraptr::downloadPBDB(CanonicalTaxa,"Cambrian","Pleistocene")

# Clean and simplify the data
CanonicalPBDB<-velociraptr::cleanTaxonomy(CanonicalPBDB,"genus")
# Remove occurrences with invalid paleocoordinates
CanonicalPBDB<-subset(CanonicalPBDB,is.na(CanonicalPBDB[,"paleolat"])!=TRUE)
# Remove unnecessary columns for increased performance
CanonicalPBDB<-CanonicalPBDB[,c("genus","family","class","early_interval","late_interval","max_ma","min_ma","paleolat","paleolng","geoplate")]

# Download ages timescale
Ages<-velociraptr::downloadTime("international%20ages")

# Sort and constrain ages
CanonicalStages<-velociraptr::constrainAges(CanonicalPBDB,Ages)
# Remove intervals with too few occurrences
CanonicalStages<-subset(CanonicalStages,CanonicalStages[,"early_interval"]%in%names(which(table(CanonicalStages[,"early_interval"])>=1000))==TRUE)
# Convert into presence matrix
CanonicalStages<-velociraptr::presenceMatrix(CanonicalStages,Rows="genus",Columns="early_interval")
# Sort the order of the intervals to be from youngest-oldest/left-right
CanonicalStages<-CanonicalStages[,rownames(subset(Ages,Ages[,"name"]%in%colnames(CanonicalStages)))]				    
				    
# For each stage, calculate whether a taxon is within-bin, through-bin, start-bin, or end-bin category of Foote
FooteFossils<-assignStatus(CanonicalStages)

# Calculate the origination and extinction rate, respecitvely
FossilQ<-calcRates(FooteFossils,"Extinction")
FossilP<-calcRates(FooteFossils,"Origination")

# Clean up NaN and Inf for youngest and oldest interval
FossilQ[is.infinite(FossilQ) | is.nan(FossilQ)] <-NA
FossilP[is.infinite(FossilP) | is.nan(FossilP)] <-NA

# Merge Rates with Ages
FossilQ<-transform(merge(as.data.frame(FossilQ),Ages,by="row.names",all=FALSE),row.names=Row.names,Row.names=NULL)
FossilP<-transform(merge(as.data.frame(FossilP),Ages,by="row.names",all=FALSE),row.names=Row.names,Row.names=NULL)

# Convert to Rates (i.e., by time)
FossilQ[,"GlobalQ"]<-FossilQ[,"GlobalQ"]/(FossilQ[,"b_age"]-FossilQ[,"t_age"])
FossilP[,"GlobalP"]<-FossilP[,"GlobalP"]/(FossilP[,"b_age"]-FossilP[,"t_age"])

# Reorder the dataset
FossilQ<-FossilQ[order(FossilQ[,"Midpoint"]),]
FossilP<-FossilP[order(FossilP[,"Midpoint"]),]

#############################################################################################################
############################################# FUNCTIONS, TRUNCATION-RATES ###################################
#############################################################################################################
# This is a hasty rewrite of some older	for analyzing macrostrat data. This simply converts the range-through of
# Macrostrat "sections" - i.e., "gap-bound-packages" - as our proxy of truncation rates
# Sections by Time matrix	
timeSections<-function(Data) {
	Data<-Data[,c("section_id","b_age","t_age")]
	Data<-na.omit(Data)
	Data[,"b_age"]<-ceiling(Data[,"b_age"])
	Data[,"t_age"]<-floor(Data[,"t_age"])	
	Data<-Data[!duplicated(Data),] # Why didn't I use unique( )? We will never know
	FinalMatrix<-matrix(0,nrow=nrow(Data),ncol=ceiling(max(Data[,"b_age"])))
	colnames(FinalMatrix)<-1:ncol(FinalMatrix)
	rownames(FinalMatrix)<-Data[,"section_id"]
	for (i in 1:nrow(Data)) {
		FinalMatrix[i,Data[i,"t_age"]:Data[i,"b_age"]]<-1
		}
	return(FinalMatrix)
	}

# Sum range through
sumbt<-function(x,y=Areas) {
	Weights<-y[which(y>0)]/sum(y[which(y>0)])
	bt<-Weights[which(x==1)]
	return(sum(bt))
	}

# sum origination denominator
sumFt<-function(x,y=Areas) {
	Weights<-y[which(y>0)]/sum(y[which(y>0)])
	Ft<-Weights[which(x==1 | x==3)]
	return(sum(Ft))
	}

# sum extinction denominator
sumbL<-function(x,y=Areas) {
	Weights<-y[which(y>0)]/sum(y[which(y>0)])
	bL<-Weights[which(x==2 | x==1)]
	return(sum(bL))
	}

# Calculat p and q
weightedRates<-function(AgeMatrix,Rate="Origination",Areas) {
	Rates<-switch(Rate,
		"Origination"=apply(AgeMatrix,2,function(col) sumbt(col,Areas)/sumFt(col,Areas)),
		"Extinction"=apply(AgeMatrix,2,function(col) sumbt(col,Areas)/sumbL(col,Areas))
		)
	return(-log(Rates))
	}			    
				    
############################################# SCRIPT, TRUNCATION-RATES ######################################
# Download north american sections
CanonicalSections<-read.csv("https://macrostrat.org/api/sections?format=csv&project_id=1",stringsAsFactors=FALSE)

# Assign each section to the 4 foote categories - e.g., range-through, singleton, etc.
FooteSections<-assignStatus(timeSections(CanonicalSections))
# Extract column areas for each section
Areas<-CanonicalSections[match(rownames(FooteSections),CanonicalSections[,"section_id"]),"col_area"]				    

# Calculate the origination and extinction rate, respecitvely
SectionsQ<-weightedRates(FooteSections,"Extinction",Areas)
SectionsP<-weightedRates(FooteSections,"Origination",Areas)

# Calculate the origination and extinction rate, respecitvely
# SectionsQ<-calcRates(FooteSections,"Extinction")
# SectionsP<-calcRates(FooteSections,"Origination")

# Clean up NaN and Inf for youngest and oldest interval
SectionsQ[is.infinite(SectionsQ) | is.nan(SectionsQ)] <- NA
SectionsP[is.infinite(SectionsP) | is.nan(SectionsP)] <- NA

# Reformat into the tab-separated values file required by PyRate
TruncateQ<-cbind(time=2:541,truncation=SectionsQ[2:541])
TruncateP<-cbind(time=2:541,truncation=SectionsP[2:541])

