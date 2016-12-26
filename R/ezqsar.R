#' ezqsar_f
#'
#' This function can easily create a MLR-QSAR model from a proper set of compounds. 
#'
#' It takes a structure file and an activity file and will give a cross-validated QSAR model as well as external test set prediction results. A table of calculated descriptors and a plot showing observed vs predicted activity of train and test sets will be generated in the working directory.    
#'
#' @param SDFfile It is a sdf file that includes structures of the all molecules (2D or 3D)
#' @param activityfile It is a csv file that contains activity data for the molecules. It should contain names of the molecules same as they are specified in a field of the SDFfile  and  the header of the column should be "name". The ranking of the compound also should be same in activityfile and SDFfile. The second column of the IC50 CSV file should contain reported activities of the molecules and the header of the column should be "IC50".
#' @param propertyfield It is a name of a field in the SDFfile that contains names of the molecules same as in the first coloumn of the activityfile (default is "title")
#' @param propertyfield_newset It is a name of a field in the newdataset SDFfile that contains names of the molecules (default is "title")
#' @param propertyfield_newset2 It is a name of a field in the newdataset2 SDFfile that contains names of the molecules (default is "title")
#' @param Nofdesc It is maximum Number of descriptors that can be present in the developed MLR model (3-6 is recomended, default is 6)
#' @param correlation It indicates level of correlation defined to omit highly correlated descriptors (default is 1)
#' @param partition It indicates the partition of the train set (default is 0.8)
#' @param des_sel_meth It defines a method for variable (descriptor) selection before MLR (it could be "exhaustive","backward", "forward" (default) or "seqrep")
#' @param testset If it is equal to zero (default) the test set will be selected according to the IC50 values and partition parameter otherwise a vector containing the row 
#' numbers of test set can be provided here (eg. c(5, 8, 12, 21)).
#' @param newdataset It is an optional sdf file that includes new set molecules to predict their activity by the developed QSAR model (2D or 3D)
#' @param newdataset2 It is a secondary optional sdf file
#' @param activity If it is equal to zero (default) it will be assumed that the reported activities are expressed in -log IC50 otherwise it will be considred #' as original IC50 values.
#' @param Cutoff It is a cut off value for reporting a descriptor value of a molecule as an outlier in AD_outlier_train, AD_outlier_test, AD_outlier_newset  and AD_outlier_newset2 tables  (default is 3). 
#' @return An object with several attributes 
#' @export
#' @examples
#'file1<-system.file("extdata", "molecules-3d.sdf", package = "ezqsar")
#'file2<-system.file("extdata", "IC50.csv", package = "ezqsar")
#'file3<-system.file("extdata", "newset-3d.sdf", package = "ezqsar")
#'model<-ezqsar_f(SDFfile=file1, activityfile=file2, newdataset=file3, testset=c(4,6,12,22))
#'attributes (model)
#'print (model$Q2)
#'print (model$R2)
#'print (model$test)
#'print (model$R2_pred)
#'print (model$Tanimoto_test_sum)
#'print (model$AD_outlier_test)
#'print (model$newset)
#'print (model$Tanimoto_newset_sum)


ezqsar_f<- function (SDFfile, activityfile,propertyfield= "title", propertyfield_newset= "title", propertyfield_newset2= "title",
Nofdesc=6 , correlation=1, partition=0.8, des_sel_meth="forward",testset=0, newdataset=0, newdataset2=0, activity=0 ,Cutoff=3)
{

centerscale="F"#"TRUE"
na.action<-stats::na.omit#This will double garantee omition of NAs  
#This is a function for descriptor generation
des_generation<- function(filename="descriptors.csv",  SDF)
{
	mols <- rcdk::load.molecules( SDF )#reading input molecules
	#Descriptor calculation
	dc <- rcdk::get.desc.categories()
	dc
	dA <- rcdk::get.desc.names(dc[1])
	dB <- rcdk::get.desc.names(dc[2])
	dC <- rcdk::get.desc.names(dc[3])
	dD <- rcdk::get.desc.names(dc[4])
	dE <- rcdk::get.desc.names(dc[5])
	allDescsA <- rcdk::eval.desc(mols, dA)
	allDescsB <- rcdk::eval.desc(mols, dB)
	allDescsC <- rcdk::eval.desc(mols, dC)
	allDescsD <- rcdk::eval.desc(mols, dD)
	allDescsE <- rcdk::eval.desc(mols, dE)
	total <- c(allDescsA,allDescsB,allDescsC,allDescsD,allDescsE)
	if (length (dc) == 6)#in some versions rcdk generate 6 descriptor categeory
	{
		dF <- rcdk::get.desc.names(dc[6])
		allDescsF <- rcdk::eval.desc(mols, dF)
		total <- c(allDescsA,allDescsB,allDescsC,allDescsD,allDescsE,allDescsF)
	
	}
	utils::write.table(total, filename , sep=",", col.names=NA, row.names=TRUE )#Descriptors output
	return (mols)
	
}
#This is a function to fix a bug in fp.sim.matrix function of fingerprint package in version 3.5.4
sim.matrix.ver2<- function(fp1,fp2,met)
{
	h=0
	w=0
	matrix_sim<-mat.or.vec(length(fp1),length(fp2))
	for (h in 1:length(fp1)) 
	{
		for (w in 1:length(fp2)) 
		{
			matrix_sim[h,w]<-fingerprint::distance(fp1[[h]],fp2[[w]],method=met) 
		}
	}
	return (matrix_sim)
}	 
mm<-Nofdesc
filename2<-activityfile########
IC50 <- utils::read.csv(file=filename2, header=TRUE, sep=",")#Reading IC50 file
if (activity!=0)
{
	#IC50$IC50<- 10^(-1*IC50$IC50)
	IC50$IC50<- -1*(log10(IC50$IC50))
}
des1<-des_generation (,SDFfile)#descriptor generation for the Main compound set 
if (newdataset!=0)
{	
	des2<-des_generation (filename="descriptors_newset.csv",newdataset)#descriptor generation for the vs compound set 
	#molsnew <- rcdk::load.molecules(newdataset)#reading input molecules
	if (propertyfield_newset == "title")
	{
		gnew<- lapply(des2, rcdk::get.title)
	}
	if (propertyfield_newset != "title")
	{
		gnew<- lapply(des2, rcdk::get.property, propertyfield_newset)#extraction of molecules name in s_m_Source_File field
	}
	MyDatanew <- utils::read.csv(file="descriptors_newset.csv", header=TRUE, sep=",")#Reading output file as an input
	MyDatanew$name <- gnew
	dfnew <- MyDatanew
	dfnew <- dfnew[,colSums(is.na(dfnew))<nrow(dfnew)]#Removing NA columns
	newnew2<-dfnew[-1]#for omiting a useless coloumn
	newnew2<-newnew2[-(ncol(newnew2))]#for omiting last column that is name column
	prep_tablenew <- newnew2	
	Trainnew <- c()
	Testnew  <- prep_tablenew
	if (centerscale=="TRUE")
	{
		#Data pre process center and scale
		preProcValues2new <- caret::preProcess(Testnew, method = c("center", "scale"))
		Test_CSnew <- stats::predict(preProcValues2new, Testnew)
	}

}
if (newdataset2!=0)
{	
	des3<-des_generation (filename="descriptors_newset2.csv",newdataset2)#descriptor generation for the vs compound set 
	#molsnew <- rcdk::load.molecules(newdataset)#reading input molecules
	if (propertyfield_newset2 == "title")
	{
		gnew2<- lapply(des3, rcdk::get.title)
	}
	if (propertyfield_newset2 != "title")
	{
		gnew2<- lapply(des3, rcdk::get.property, propertyfield_newset2)#extraction of molecules name in s_m_Source_File field
	}
	MyDatanew2 <- utils::read.csv(file="descriptors_newset2.csv", header=TRUE, sep=",")#Reading output file as an input
	MyDatanew2$name <- gnew2
	dfnew2 <- MyDatanew2
	dfnew2 <- dfnew2[,colSums(is.na(dfnew2))<nrow(dfnew2)]#Removing NA columns
	newnew3<-dfnew2[-1]#for omiting a useless coloumn
	newnew3<-newnew3[-(ncol(newnew3))]#for omiting last column that is name column
	prep_tablenew2 <- newnew3	
	Trainnew2 <- c()
	Testnew2  <- prep_tablenew2
	if (centerscale=="TRUE")
	{
		#Data pre process center and scale
		preProcValues2new2 <- caret::preProcess(Testnew2, method = c("center", "scale"))
		Test_CSnew2 <- stats::predict(preProcValues2new2, Testnew2)
	}

}
if (propertyfield == "title")
{
	g<- lapply(des1, rcdk::get.title)
}
if (propertyfield != "title")
{
	g<- lapply(des1, rcdk::get.property, propertyfield)#extraction of molecules name in propertyfield variable
}
MyData <- utils::read.csv(file="descriptors.csv", header=TRUE, sep=",")#Reading output file as an input
MyData$name <- g
merged <- merge (IC50,MyData,by="name", sort=F)#Merge two tables by names, dont sort by names
df <- merged
df <- df[,colSums(is.na(df))<nrow(df)]#Removing NA columns
#removing near zero variance
nzv <- caret::nearZeroVar(df, freqCut = 95/5, uniqueCut = 10)#If you get Near zero variance warning for test set you can avoid it by adjusting these two parameters
new<- df[, -nzv]
new2<-new[-1]#remove name coloumn
new2<-new2[-2]#remove a useless coloumn called X
#removing highly correlated coloumns
descrCor <-  stats::cor(new2[2:ncol(new2)])#excepton for IC50 column (first column)
highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = correlation)
if (length(highlyCorDescr)== 0)
{
	prep_table <- new2
}
if (length(highlyCorDescr)!= 0)
{
	prep_table <- new2[,-(highlyCorDescr+1)]#final table without correlated coloumns
}
#Data splitting for both descriptor table and fragments
#set.seed(SEED)
if (testset==0)
{
	if (partition != 1)
	{
		trainIndex <- caret::createDataPartition(prep_table$IC50, p = partition, list = FALSE, times = 1)#Data partitioning according to the IC50
		Train <- prep_table[ trainIndex,]
		Test  <- prep_table[-trainIndex,]
		train.mols <- des1[trainIndex]
		test.mols <- des1[-trainIndex]
		fps_test <- lapply(test.mols, rcdk::get.fingerprint, type='maccs')
		
	} 
}
if (testset!=0)
{
	Train <- prep_table[-testset,]
	Test  <- prep_table[testset,]
	train.mols <- des1[-testset]
	test.mols <- des1[testset]
	fps_test <- lapply(test.mols, rcdk::get.fingerprint, type='maccs')
}
if (partition == 1)
{
	Train <- prep_table
	Test  <- c()
	train.mols <- des1
	test.mols <- c()
}
fps_train <- lapply(train.mols, rcdk::get.fingerprint, type='maccs')
#fp.train <- fingerprint::fp.sim.matrix(fps_train, method='tanimoto')#A complete table of Tanimoto similarity index for train set
fp.train<-sim.matrix.ver2(fps_train, fps_train, "tanimoto")
#fp.test <- fingerprint::fp.sim.matrix(fps_test, fps_train, method='tanimoto')#A complete table of Tanimoto similarity index for test set against train set
fp.test<-sim.matrix.ver2(fps_test, fps_train, "tanimoto")
#set the y (outcome) and x (predictors) coloumns
MATRIX<- as.matrix(Train)
if (partition != 1)
{
	MATRIX_T<- as.matrix(Test)
}
if (centerscale=="TRUE")
{
	#Data pre process center and scale
	preProcValues <- caret::preProcess(Train[,-1], method = c("center", "scale"))#excepton for IC50 column (first column, very important indeed
	Train_CS <- stats::predict(preProcValues, Train)
	MATRIX<- as.matrix(Train_CS)
	if (partition != 1)
	{
		preProcValues3 <- caret::preProcess(Test[,-1], method = c("center", "scale"))#excepton for IC50 column (first column, very important indeed
		Test_CS <- stats::predict(preProcValues3, Test)
		MATRIX_T<- as.matrix(Test_CS)
	}
}
nn<-ncol(MATRIX)#number of descriptors after pre processing steps
predictors <-MATRIX[,2:nn]
outcome <-MATRIX[,1]
#Variable selection
#method=c("exhaustive","backward", "forward", "seqrep"), really.big=FALSE,...)
#Use exhaustive search, forward selection, backward selection or sequential replacement to search.
Selected_variables<-leaps::regsubsets(y=outcome,x=predictors,nbest=1, nvmax=Nofdesc-1,method=(des_sel_meth),all.best=FALSE, really.big=TRUE)
t<-summary(Selected_variables)
sortr<-sort(t$adjr2)#Sorted adjusted R2 vector
l<-length(sortr)#Length of the vector
v<-sortr[l]#Largest adjusted R2
p<-match(c(v),t$adjr2)#Place of the largest adjusted R2 in the vector
stats::coef(Selected_variables, p)
D<-stats::coef(Selected_variables, p)#Intercept and coefficients for the selected descriptors
d<-c()#making an empty vector for selected variable names
ph<-c()#making an empty vector for selected variable places
#A loop to find the place (column) of the selected variables in the MATRIX
mm<-(length(D)-1)#because in forward selection number of variables are unpredictable
for (i in 1:mm+1)
{
	d<-c(d,names(D[i]))
	ph<-c(ph,match(c(d[i-1]),names(Train))) 
}
#After variable seletion
objlm <- caret::train(x = MATRIX[,ph[1:mm]], y = outcome , method = "lm", trControl = caret::trainControl(method = "LOOCV"), verbose = TRUE)#cross validation methods: "LOOCV", 
s_objlm<-summary(objlm)
if (partition != 1)
{
	#Prediction on TEST set 
	objtest<-stats::predict(objlm, newdata=MATRIX_T[,ph[1:mm]])
}
#R2 pred
M<-mean(Train$IC50)
SD=0
PRESS=0
if (partition != 1)
{
	for (i in 1:length(objtest))
	{
		PRESS<-(Test$IC50[i]-objtest[i])^2+PRESS
		SD<-(Test$IC50[i]-M)^2+SD
	}
	RPred<-1-(PRESS/SD)
	Test$IC50#Observed Y for test set
	RPred#R2 pred
	objtest#Predicted Y for test set
}
#Plotting
allpointsx<-c(outcome, Test$IC50)
allpointsy<- c(objlm$pred$pred, objtest)
grDevices::pdf("Plot-MLR.pdf")
graphics::plot(outcome, objlm$pred$pred,  pch=19, col="blue", xlab="Observed", ylab="Predicted", xlim=c(min(allpointsx)-0.5, max(allpointsx)+0.5), ylim=c(min(allpointsy)-0.5, max(allpointsy)+1))#Plot observed vs predicted for train set
if (partition != 0.5)
{
	graphics::points(Test$IC50, objtest, pch=2, col= "red")#Adding Test compounds to the plot
}
grDevices::dev.off()
#Aplicabality domain based on descriptor values and fingerprints
dimM<-dimnames(Train)[[1]]#Vector for Train set nrow
CI<-mat.or.vec(nrow(MATRIX),mm)#Empty matrix for complete table of indexes for train set
Out<-mat.or.vec((nrow(MATRIX))*mm,5)#Empty matrix for train set outliers 
colnames(Out)<- c("Row number in train set", "Row number in entire set", "Molecule name", 
"Out of range descriptor" ,"Standardized value")# 
dimM_name<-c()#Vector for Train set names
Tanimoto<-mat.or.vec(nrow(MATRIX),2)#Empty matrix to store Tanimoto similarity indexes for train set
colnames(Tanimoto)<- c("Minimum Tanimoto similarity index", "Average Tanimoto similarity index")
MEAN<-c()#Parameter needed for calculating CI
SD2<-c()#Parameter needed for calculating CI
for (i in 1:mm)#i for descriptor column after variable selection
{
	MEAN[i]<- mean (MATRIX[,ph[i]])
	SD2[i]<- stats::sd (MATRIX[,ph[i]])
	for (j in 1:nrow(MATRIX))#j for a compound in train set
	{
		CI[j,i]<- abs((MATRIX[j,ph[i]]- MEAN[i])/SD2[i])
		temp_yy<-strtoi(dimM[j])
		dimM_name[j]<-g[[temp_yy]]
		Tanimoto[j,1]<- min (fp.train[j,])
		Tanimoto[j,2]<- mean (fp.train[j,])
		if (CI[j,i]>Cutoff)#only train set outliers will be listed in Out
		{
			Out[j+((i-1)*nrow(MATRIX)),1]<- j
			Out[j+((i-1)*nrow(MATRIX)),2]<- dimM[j]
			Out[j+((i-1)*nrow(MATRIX)),3]<- g[[temp_yy]]
			Out[j+((i-1)*nrow(MATRIX)),4]<- d[i]
			Out[j+((i-1)*nrow(MATRIX)),5]<- CI[j,i]
		}
	}
	if (partition != 1)
	{
		if (i == 1)
		{
			dimM_T<-dimnames(MATRIX_T)[[1]]#Vector for Test set nrow
			CI_Test<-mat.or.vec(nrow(MATRIX_T),mm)#Empty matrix for complete table of indexes for test set
			Out_Test<-mat.or.vec((nrow(MATRIX_T))*mm,5)#Empty matrix for test set outliers 
			colnames(Out_Test)<- c("Row number in test set", "Row number in entire set", "Molecule name", "Out of range descriptor" ,"Standardized value")
			dimM_T_name<-c()#Vector for Test set names
			Tanimoto_T<-mat.or.vec(nrow(MATRIX_T),3)#Empty matrix to store Tanimoto similarity indexes for test set
			colnames(Tanimoto_T)<-c("Maximum Tanimoto similarity index","Minimum Tanimoto similarity index", "Average Tanimoto similarity index")
		}
		for (n in 1:nrow(MATRIX_T))#n for a compound in train set
		{
			CI_Test[n,i]<- abs((MATRIX_T[n,ph[i]]- MEAN[i])/SD2[i])
			temp_zz<-strtoi(dimM_T[n])
			dimM_T_name[n]<-g[[temp_zz]]
			Tanimoto_T[n,1]<- max (fp.test[n,])
			Tanimoto_T[n,2]<- min (fp.test[n,])
			Tanimoto_T[n,3]<- mean (fp.test[n,])
			if (CI_Test[n,i]>Cutoff)#only test set outliers will be listed in OUT
			{
				Out_Test[n+((i-1)*nrow(MATRIX_T)),1]<- n
				Out_Test[n,2]<- dimM_T[n]
				Out_Test[n+((i-1)*nrow(MATRIX_T)),3]<- g[[temp_zz]]
				Out_Test[n+((i-1)*nrow(MATRIX_T)),4]<- d[i]
				Out_Test[n+((i-1)*nrow(MATRIX_T)),5]<- CI_Test[n,i]
			}
		}
		if (i==mm)
		{
			Out_Test2 <- Out_Test[Out_Test[,1]>0,]#removing zeros from Out_Test
			colnames(CI_Test)<- d
			if (length (Out_Test2)== 0)
			{
				Out_Test2<-0
			}
			rownames(fp.test)<- dimM_T_name
			colnames(fp.test)<- dimM_name#dimM
			rownames(Tanimoto_T)<-dimM_T_name#dimM
		}
	}
	if (newdataset!=0)
	{
		if (i==1)
		{
			Testnewd<-Testnew[,d]#Extraction of selecet variables column from vs table
			MATRIX_Tnew<- as.matrix(Testnewd)#convert it to matrix
			dimM_Tnew<- rownames(Testnew)#Vector for Train set nrows
			CI_Testnew<-mat.or.vec(nrow(MATRIX_Tnew),mm)#Empty matrix for complete table of indexes for vs set
			Out_Testnew<-mat.or.vec((nrow(MATRIX_Tnew))*mm,4)#Empty matrix for vs set outliers 
			colnames(Out_Testnew)<- c("Row number in the new test set", "Molecule name", "Out of range descriptor" ,"Standardized value")
			dimM_Tnew_name<-c()#Vector for Train set names
			Tanimoto_Tnew<-mat.or.vec(nrow(MATRIX_Tnew),3)#Empty matrix to store Tanimoto similarity indexes for new set
			colnames(Tanimoto_Tnew)<-colnames(Tanimoto_T)
			fps_vs <- lapply(des2, rcdk::get.fingerprint, type='maccs')
			#fp.vs <- fingerprint::fp.sim.matrix(fps_vs, fps_train, method='tanimoto')#Tanimoto table agains train set
			fp.vs<-sim.matrix.ver2(fps_vs, fps_train, "tanimoto")
		}	
		if (centerscale=="TRUE")
		{	
			Test_CSnewd<-Test_CSnew[,d]#Extraction of selecet variables column from vs table
			MATRIX_Tnew<- as.matrix(Test_CSnewd)#convert it to matrix
		}
		for (l in 1:nrow(MATRIX_Tnew))#l for a compound in train set
		{
			CI_Testnew[l,i]<- abs((MATRIX_Tnew[l,d[i]]- MEAN[i])/SD2[i])
			temp_zznew<-strtoi(dimM_Tnew[l])
			dimM_Tnew_name[l]<-gnew[[temp_zznew]]
			Tanimoto_Tnew[l,1]<- max (fp.vs[l,])
			Tanimoto_Tnew[l,2]<- min (fp.vs[l,])
			Tanimoto_Tnew[l,3]<- mean (fp.vs[l,])

			if (CI_Testnew[l,i]>Cutoff)#only newset outliers will be listed in OUT_Testnew
			{
				Out_Testnew[l+((i-1)*nrow(MATRIX_Tnew)),1]<- l
				Out_Testnew[l+((i-1)*nrow(MATRIX_Tnew)),2]<- gnew[[temp_zznew]]
				Out_Testnew[l+((i-1)*nrow(MATRIX_Tnew)),3]<- d[i]
				Out_Testnew[l+((i-1)*nrow(MATRIX_Tnew)),4]<- CI_Testnew[l,i]

			}
		}
		if (i==mm)
		{
			Out_Testnew2 <- Out_Testnew[Out_Testnew[,1]>0,]#removing zeros from Out_Testnew
			colnames(CI_Testnew)<- d
			if (length (Out_Testnew2)== 0)
			{
				Out_Testnew2<-0
			}
			objtestnew<-stats::predict(objlm, newdata=MATRIX_Tnew)
			rownames(fp.vs)<- dimM_Tnew_name#rownames(Testnewd)[1]#~dimM_Tnew
			colnames(fp.vs)<- dimM_name#dimM
			rownames(Tanimoto_Tnew)<-dimM_Tnew_name#dimM
	
		}
	}
	if (newdataset2!=0)
	{
		if (i==1)
		{
			Testnewd2<-Testnew2[,d]#Extraction of selecet variables column from vs table
			MATRIX_Tnew2<- as.matrix(Testnewd2)#convert it to matrix
			dimM_Tnew2<- rownames(Testnew2)#Vector for Train set nrows
			CI_Testnew2<-mat.or.vec(nrow(MATRIX_Tnew2),mm)#Empty matrix for complete table of indexes for vs set
			Out_Testnew3<-mat.or.vec((nrow(MATRIX_Tnew2))*mm,4)#Empty matrix for vs set outliers 
			colnames(Out_Testnew3)<- c("Row number in the new test set 2", "Molecule name", "Out of range descriptor" ,"Standardized value")
			dimM_Tnew_name2<-c()#Vector for Train set names
			Tanimoto_Tnew2<-mat.or.vec(nrow(MATRIX_Tnew2),3)#Empty matrix to store Tanimoto similarity indexes for new set
			colnames(Tanimoto_Tnew2)<-colnames(Tanimoto_T)
			fps_vs2 <- lapply(des3, rcdk::get.fingerprint, type='maccs')
			#fp.vs <- fingerprint::fp.sim.matrix(fps_vs2, fps_train, method='tanimoto')#Tanimoto table agains train set
			fp.vs2<-sim.matrix.ver2(fps_vs2, fps_train, "tanimoto")
		}	
		if (centerscale=="TRUE")
		{	
			Test_CSnewd2<-Test_CSnew2[,d]#Extraction of selecet variables column from vs table
			MATRIX_Tnew2<- as.matrix(Test_CSnewd2)#convert it to matrix
		}
		for (l in 1:nrow(MATRIX_Tnew2))#l for a compound in train set
		{
			CI_Testnew2[l,i]<- abs((MATRIX_Tnew2[l,d[i]]- MEAN[i])/SD2[i])
			temp_zznew2<-strtoi(dimM_Tnew2[l])
			dimM_Tnew_name2[l]<-gnew2[[temp_zznew2]]
			Tanimoto_Tnew2[l,1]<- max (fp.vs2[l,])
			Tanimoto_Tnew2[l,2]<- min (fp.vs2[l,])
			Tanimoto_Tnew2[l,3]<- mean (fp.vs2[l,])

			if (CI_Testnew2[l,i]>Cutoff)#only newset outliers will be listed in OUT_Testnew
			{
				Out_Testnew3[l+((i-1)*nrow(MATRIX_Tnew2)),1]<- l
				Out_Testnew3[l+((i-1)*nrow(MATRIX_Tnew2)),2]<- gnew2[[temp_zznew2]]
				Out_Testnew3[l+((i-1)*nrow(MATRIX_Tnew2)),3]<- d[i]
				Out_Testnew3[l+((i-1)*nrow(MATRIX_Tnew2)),4]<- CI_Testnew2[l,i]

			}
		}
		if (i==mm)
		{
			Out_Testnew22 <- Out_Testnew3[Out_Testnew3[,1]>0,]#removing zeros from Out_Testnew
			colnames(CI_Testnew2)<- d
			if (length (Out_Testnew22)== 0)
			{
				Out_Testnew22<-0
			}
			objtestnew2<-stats::predict(objlm, newdata=MATRIX_Tnew2)
			rownames(fp.vs2)<- dimM_Tnew_name2#rownames(Testnewd)[1]#~dimM_Tnew
			colnames(fp.vs2)<- dimM_name#dimM
			rownames(Tanimoto_Tnew2)<-dimM_Tnew_name2#dimM
	
		}
	}	

}
Out2 <- Out[Out[,1]>0,]#removing zeros from Out
colnames(CI)<- d
if (length (Out2)== 0)
{
	Out2<-0
}
rownames(fp.train)<- dimM_name#dimM
colnames(fp.train)<- dimM_name#dimM
rownames(Tanimoto)<- dimM_name#dimM
#Setting not used variables to zero to bypass getting strange outputs
if (partition == 1)
{
	Test$IC50=0
	objtest=0
	RPred=0
	CI_Test=0
	Out_Test=0
	fp.test=0
}
if (newdataset==0)
{
	CI_Testnew=0
	Out_Testnew2=0
	objtestnew=0
	fp.vs=0
	Tanimoto_Tnew=0
	#dimM_Tnew_name=0
	MATRIX_objtestnew=0
}
if (newdataset2==0)
{
	CI_Testnew2=0
	Out_Testnew22=0
	objtestnew2=0
	fp.vs2=0
	Tanimoto_Tnew2=0
	#dimM_Tnew_name2=0
	MATRIX_objtestnew2=0
}

#Making rows and column names for the final tables
MATRIX_objtrain<-mat.or.vec(length(objlm$pred$pred),2)
rownames(MATRIX_objtrain)<-dimM_name 
colnames(MATRIX_objtrain)<- c("train_observed", "train_predicted")
MATRIX_objtrain[,1]<-outcome
MATRIX_objtrain[,2]<-objlm$pred$pred

MATRIX_objtest<-mat.or.vec(length(objtest),2)
rownames(MATRIX_objtest)<-dimM_T_name 
colnames(MATRIX_objtest)<- c("test_observed", "test_predicted")
MATRIX_objtest[,1]<-Test$IC50
MATRIX_objtest[,2]<-objtest

rownames(CI)<-dimM_name 
rownames(CI_Test)<-dimM_T_name

if (newdataset!=0)
{
MATRIX_objtestnew<-mat.or.vec(length(objtestnew),2)
colnames(MATRIX_objtestnew)<- c("name","newset_predicted")
MATRIX_objtestnew[,1]<-dimM_Tnew_name
MATRIX_objtestnew[,2]<-objtestnew
rownames(CI_Testnew)<-dimM_Tnew_name 
} 
if (newdataset2!=0)
{
MATRIX_objtestnew2<-mat.or.vec(length(objtestnew2),2)
colnames(MATRIX_objtestnew2)<- c("name","newset2_predicted")
MATRIX_objtestnew2[,2]<-objtestnew2
MATRIX_objtestnew2[,1]<-dimM_Tnew_name2 
rownames(CI_Testnew2)<-dimM_Tnew_name2 
} 


result <- list(Q2=objlm$results$Rsquared,R2=s_objlm$r.squared, adjusted_R2=s_objlm$adj.r.square, RMSE=objlm$results$RMSE, 
F_statistics=s_objlm$fstatistic, Descriptors=d, Model=objlm$finalModel, train=MATRIX_objtrain, 
test=MATRIX_objtest, newset=MATRIX_objtestnew ,newset2=MATRIX_objtestnew2 ,R2_pred=RPred,
standardized_des_train=CI, standardized_des_test=CI_Test,  standardized_des_newset=CI_Testnew, standardized_des_newset2=CI_Testnew2,
AD_outlier_train=Out2, AD_outlier_test=Out_Test2, AD_outlier_newset=Out_Testnew2, AD_outlier_newset2=Out_Testnew22,
Tanimoto_train=fp.train, Tanimoto_test=fp.test, Tanimoto_newset=fp.vs, Tanimoto_newset2=fp.vs2,
Tanimoto_train_sum=Tanimoto, Tanimoto_test_sum=Tanimoto_T, Tanimoto_newset_sum=Tanimoto_Tnew, Tanimoto_newset2_sum=Tanimoto_Tnew2)
}
