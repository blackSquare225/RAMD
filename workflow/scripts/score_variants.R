library(data.table)

args = commandArgs(trailingOnly=TRUE)

fi <- args[1]
fo <- args[2]

df <- fread(fi, header=TRUE, sep="\t")

name_ <- unlist(strsplit(basename(fi[1]),".varscan.vep.tsv"))
dir <- dirname(fi[1])

wanted <- c("SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred","MutationAssessor_pred","FATHMM_pred","fathmm-XF_coding_pred","fathmm-MKL_coding_pred","PROVEAN_pred",
	"MetaSVM_pred","MetaLR_pred","ClinPred_pred","DEOGEN2_pred","LIST-S2_pred","MetaRNN_pred","SIFT4G_pred","ESM1b_pred","M-CAP_pred")

SCORE <- c()

for (i in 1:nrow(df)){

	if (all(df[i,..wanted] == ".")){

		resultNULL <- 99999
		SCORE[i] <- resultNULL

	} else {

		if (df[i,"SIFT_pred"]!="."){

			if (grepl("&",df[i,"SIFT_pred"])){

				SIFT <- unlist(strsplit(as.character(df[i,"SIFT_pred"]),"&", fixed=TRUE))
				SIFT <- sapply(SIFT, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultSIFT <- sum(SIFT, na.rm=TRUE)/length(which(!is.na(SIFT)))


			} else {

				SIFT <- sapply(df[i,"SIFT_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultSIFT <- sum(SIFT)/length(which(!is.na(SIFT)))
			}

		} else {

			resultSIFT <- NA

		} 

		if (df[i,"Polyphen2_HDIV_pred"]!="."){

			if (grepl("&",df[i,"Polyphen2_HDIV_pred"])){

				Polyphen2_HDIV_pred <- unlist(strsplit(as.character(df[i,"Polyphen2_HDIV_pred"]),"&", fixed=TRUE))
				Polyphen2_HDIV_pred <- sapply(Polyphen2_HDIV_pred, function(x) ifelse(x=="D",2,ifelse(x=="P",1,ifelse(x=="B",0,NA))))
				resultHDIV <- sum(Polyphen2_HDIV_pred, na.rm=TRUE)/length(which(!is.na(Polyphen2_HDIV_pred)))


			} else {

				Polyphen2_HDIV_pred <- sapply(df[i,"Polyphen2_HDIV_pred"], function(x) ifelse(x=="D",2,ifelse(x=="P",1,ifelse(x=="B",0,NA))))
				resultHDIV <- sum(Polyphen2_HDIV_pred)/length(which(!is.na(Polyphen2_HDIV_pred)))

			}

		} else {

			resultHDIV <- NA

		}

		if (df[i,"Polyphen2_HVAR_pred"]!="."){

			if (grepl("&",df[i,"Polyphen2_HVAR_pred"])){

				Polyphen2_HVAR_pred <- unlist(strsplit(as.character(df[i,"Polyphen2_HVAR_pred"]),"&", fixed=TRUE))
				Polyphen2_HVAR_pred <- sapply(Polyphen2_HVAR_pred, function(x) ifelse(x=="D",2,ifelse(x=="P",1,ifelse(x=="B",0,NA))))
				resultHVAR <- sum(Polyphen2_HVAR_pred, na.rm=TRUE)/length(which(!is.na(Polyphen2_HVAR_pred)))


			} else {

				Polyphen2_HVAR_pred <- sapply(df[i,"Polyphen2_HVAR_pred"], function(x) ifelse(x=="D",2,ifelse(x=="P",1,ifelse(x=="B",0,NA))))
				resultHVAR <- sum(Polyphen2_HVAR_pred)/length(which(!is.na(Polyphen2_HVAR_pred)))
				
			}		

		} else {

			resultHVAR <- NA

		}

		if (df[i,"LRT_pred"]!="."){

			if (grepl("&",df[i,"LRT_pred"])){

				LRT_pred <- unlist(strsplit(as.character(df[i,"LRT_pred"]),"&", fixed=TRUE))
				LRT_pred <- sapply(LRT_pred, function(x) ifelse(x=="D",2,ifelse(x=="U",1,ifelse(x=="N",0,NA))))
				resultLRT <- sum(LRT_pred, na.rm=TRUE)/length(which(!is.na(LRT_pred)))


			} else {

				LRT_pred <- sapply(df[i,"LRT_pred"], function(x) ifelse(x=="D",2,ifelse(x=="U",1,ifelse(x=="N",0,NA))))
				resultLRT <- sum(LRT_pred)/length(which(!is.na(LRT_pred)))
				
			}

		} else {

			resultLRT <- NA

		}

		if (df[i,"MutationTaster_pred"]!="."){

			if (grepl("&",df[i,"MutationTaster_pred"])){

				MutationTaster_pred <- unlist(strsplit(as.character(df[i,"MutationTaster_pred"]),"&", fixed=TRUE))
				MutationTaster_pred <- sapply(MutationTaster_pred, function(x) ifelse(x=="A",2,ifelse(x=="D",1,ifelse(x=="N",0.5,ifelse(x=="P",0,NA)))))
				resultMT <- sum(MutationTaster_pred, na.rm=TRUE)/length(which(!is.na(MutationTaster_pred)))


			} else {

				MutationTaster_pred <- sapply(df[i,"MutationTaster_pred"], function(x) ifelse(x=="A",2,ifelse(x=="D",1,ifelse(x=="N",0.5,ifelse(x=="P",0,NA)))))
				resultMT <- sum(MutationTaster_pred)/length(which(!is.na(MutationTaster_pred)))
				
			}		

		} else {

			resultMT <- NA

		}

		if (df[i,"MutationAssessor_pred"]!="."){

			if (grepl("&",df[i,"MutationAssessor_pred"])){

				MutationAssessor_pred <- unlist(strsplit(as.character(df[i,"MutationAssessor_pred"]),"&", fixed=TRUE))
				MutationAssessor_pred <- sapply(MutationAssessor_pred, function(x) ifelse(x=="H",2,ifelse(x=="M",1,ifelse(x=="L",0.5,ifelse(x=="N",0,NA)))))
				resultMA <- sum(MutationAssessor_pred, na.rm=TRUE)/length(which(!is.na(MutationAssessor_pred)))


			} else {

				MutationAssessor_pred <- sapply(df[i,"MutationAssessor_pred"], function(x) ifelse(x=="H",2,ifelse(x=="M",1,ifelse(x=="L",0.5,ifelse(x=="N",0,NA)))))
				resultMA <- sum(MutationAssessor_pred)/length(which(!is.na(MutationAssessor_pred)))
				
			}

		} else {

			resultMA <- NA

		}

		if (df[i,"FATHMM_pred"]!="."){

			if (grepl("&",df[i,"FATHMM_pred"])){

				FATHMM_pred <- sapply(df[i,"FATHMM_pred"], function(x) ifelse(x=="D",2,0))
				FATHMM_pred <- sapply(FATHMM_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultFATHMM <- sum(FATHMM_pred, na.rm=TRUE)/length(which(!is.na(FATHMM_pred)))


			} else {

				FATHMM_pred <- sapply(df[i,"FATHMM_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultFATHMM <- sum(FATHMM_pred)/length(which(!is.na(FATHMM_pred)))
			}

		} else {

			resultFATHMM <- NA

		}

		if (df[i,"fathmm-XF_coding_pred"]!="."){

			if (grepl("&",df[i,"fathmm-XF_coding_pred"])){

				fathmm_XF_coding_pred <- unlist(strsplit(as.character(df[i,"fathmm-XF_coding_pred"]),"&", fixed=TRUE))
				fathmm_XF_coding_pred <- sapply(fathmm_XF_coding_pred, function(x) ifelse(x=="D",2,ifelse(x=="N",0,NA)))
				resultXF <- sum(fathmm_XF_coding_pred, na.rm=TRUE)/length(which(!is.na(fathmm_XF_coding_pred)))


			} else {

				fathmm_XF_coding_pred <- sapply(df[i,"fathmm-XF_coding_pred"], function(x) ifelse(x=="D",2,ifelse(x=="N",0,NA)))
				resultXF <- sum(fathmm_XF_coding_pred)/length(which(!is.na(fathmm_XF_coding_pred)))

			}

		} else {

			resultXF <- NA

		}

		if (df[i,"fathmm-MKL_coding_pred"]!="."){

			if (grepl("&",df[i,"fathmm-MKL_coding_pred"])){

				fathmm_MKL_coding_pred <- unlist(strsplit(as.character(df[i,"fathmm-MKL_coding_pred"]),"&", fixed=TRUE))
				fathmm_MKL_coding_pred <- sapply(fathmm_MKL_coding_pred, function(x) ifelse(x=="D",2,ifelse(x=="N",0,NA)))
				resultMKL <- sum(fathmm_MKL_coding_pred, na.rm=TRUE)/length(which(!is.na(fathmm_MKL_coding_pred)))


			} else {

				fathmm_MKL_coding_pred <- sapply(df[i,"fathmm-MKL_coding_pred"], function(x) ifelse(x=="D",2,ifelse(x=="N",0,NA)))
				resultMKL <- sum(fathmm_MKL_coding_pred)/length(which(!is.na(fathmm_MKL_coding_pred)))

			}

		} else {

			resultMKL <- NA

		}

		if (df[i,"PROVEAN_pred"]!="."){

			if (grepl("&",df[i,"PROVEAN_pred"])){

				PROVEAN_pred <- unlist(strsplit(as.character(df[i,"PROVEAN_pred"]),"&", fixed=TRUE))
				PROVEAN_pred <- sapply(PROVEAN_pred, function(x) ifelse(x=="D",2,ifelse(x=="N",0,NA)))
				resultPROV <- sum(PROVEAN_pred, na.rm=TRUE)/length(which(!is.na(PROVEAN_pred)))


			} else {

				PROVEAN_pred <- sapply(df[i,"PROVEAN_pred"], function(x) ifelse(x=="D",2,ifelse(x=="N",0,NA)))
				resultPROV <- sum(PROVEAN_pred)/length(which(!is.na(PROVEAN_pred)))

			}
		
		} else {

			resultPROV <- NA

		}

		if (df[i,"MetaSVM_pred"]!="."){

			if (grepl("&",df[i,"MetaSVM_pred"])){

				MetaSVM_pred <- unlist(strsplit(as.character(df[i,"MetaSVM_pred"]),"&", fixed=TRUE))
				MetaSVM_pred <- sapply(MetaSVM_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultSVM <- sum(MetaSVM_pred, na.rm=TRUE)/length(which(!is.na(MetaSVM_pred)))


			} else {

				MetaSVM_pred <- sapply(df[i,"MetaSVM_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultSVM <- sum(MetaSVM_pred)/length(which(!is.na(MetaSVM_pred)))

			}

		} else {

			resultSVM <- NA

		}

		if (df[i,"MetaLR_pred"]!="."){

			if (grepl("&",df[i,"MetaLR_pred"])){

				MetaLR_pred <- unlist(strsplit(as.character(df[i,"MetaLR_pred"]),"&", fixed=TRUE))
				MetaLR_pred <- sapply(MetaLR_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultLR <- sum(MetaLR_pred, na.rm=TRUE)/length(which(!is.na(MetaLR_pred)))

			} else {

				MetaLR_pred <- sapply(df[i,"MetaLR_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultLR <- sum(MetaLR_pred)/length(which(!is.na(MetaLR_pred)))

			}		

		} else {

			resultLR <- NA

		}

		if (df[i,"ClinPred_pred"]!="."){


			if (grepl("&",df[i,"ClinPred_pred"])){

				ClinPred_pred <- unlist(strsplit(as.character(df[i,"ClinPred_pred"]),"&", fixed=TRUE))
				ClinPred_pred <- sapply(ClinPred_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultClin <- sum(ClinPred_pred, na.rm=TRUE)/length(which(!is.na(ClinPred_pred)))


			} else {

				ClinPred_pred <- sapply(df[i,"ClinPred_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultClin <- sum(ClinPred_pred)/length(which(!is.na(ClinPred_pred)))

			}

		} else {

			resultClin <- NA

		}

		if (df[i,"DEOGEN2_pred"]!="."){

			if (grepl("&",df[i,"DEOGEN2_pred"])){

				DEOGEN2_pred <- unlist(strsplit(as.character(df[i,"DEOGEN2_pred"]),"&", fixed=TRUE))
				DEOGEN2_pred <- sapply(DEOGEN2_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultDEO <- sum(DEOGEN2_pred, na.rm=TRUE)/length(which(!is.na(DEOGEN2_pred)))


			} else {

				DEOGEN2_pred <- sapply(df[i,"DEOGEN2_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultDEO <- sum(DEOGEN2_pred)/length(which(!is.na(DEOGEN2_pred)))

			}

		} else {

			resultDEO <- NA

		}

		if (df[i,"LIST-S2_pred"]!="."){

			if (grepl("&",df[i,"LIST-S2_pred"])){

				LIST_S2_pred <- unlist(strsplit(as.character(df[i,"LIST-S2_pred"]),"&", fixed=TRUE))
				LIST_S2_pred <- sapply(LIST_S2_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultLIST <- sum(LIST_S2_pred, na.rm=TRUE)/length(which(!is.na(LIST_S2_pred)))


			} else {

				LIST_S2_pred <- sapply(df[i,"LIST-S2_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultLIST <- sum(LIST_S2_pred)/length(which(!is.na(LIST_S2_pred)))

			}

		} else {

			resultLIST <- NA

		}

		if (df[i,"MetaRNN_pred"]!="."){

			if (grepl("&",df[i,"MetaRNN_pred"])){

				MetaRNN_pred <- unlist(strsplit(as.character(df[i,"MetaRNN_pred"]),"&", fixed=TRUE))
				MetaRNN_pred <- sapply(MetaRNN_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultRNN <- sum(MetaRNN_pred, na.rm=TRUE)/length(which(!is.na(MetaRNN_pred)))


			} else {

				MetaRNN_pred <- sapply(df[i,"MetaRNN_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultRNN <- sum(MetaRNN_pred)/length(which(!is.na(MetaRNN_pred)))

			}

		} else {

			resultRNN <- NA

		}

		if (df[i,"SIFT4G_pred"]!="."){

			if (grepl("&",df[i,"SIFT4G_pred"])){

				SIFT4G_pred <- unlist(strsplit(as.character(df[i,"SIFT4G_pred"]),"&", fixed=TRUE))
				SIFT4G_pred <- sapply(SIFT4G_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				result4G <- sum(SIFT4G_pred, na.rm=TRUE)/length(which(!is.na(SIFT4G_pred)))


			} else {

				SIFT4G_pred <- sapply(df[i,"SIFT4G_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				result4G <- sum(SIFT4G_pred)/length(which(!is.na(SIFT4G_pred)))

			}

		} else {

			result4G <- NA

		}

		if (df[i,"ESM1b_pred"]!="."){

			if (grepl("&",df[i,"ESM1b_pred"])){

				ESM1b_pred <- unlist(strsplit(as.character(df[i,"ESM1b_pred"]),"&", fixed=TRUE))
				ESM1b_pred <- sapply(ESM1b_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultES <- sum(ESM1b_pred, na.rm=TRUE)/length(which(!is.na(ESM1b_pred)))


			} else {

				ESM1b_pred <- sapply(df[i,"ESM1b_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultES <- sum(ESM1b_pred)/length(which(!is.na(ESM1b_pred)))

			}

		} else {

			resultES <- NA

		}

		if (df[i,"M-CAP_pred"]!="."){

			if (grepl("&",df[i,"M-CAP_pred"])){

				M_CAP_pred <- unlist(strsplit(as.character(df[i,"M-CAP_pred"]),"&", fixed=TRUE))
				M_CAP_pred <- sapply(M_CAP_pred, function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultMCAP <- sum(M_CAP_pred, na.rm=TRUE)/length(which(!is.na(M_CAP_pred)))


			} else {

				M_CAP_pred <- sapply(df[i,"M-CAP_pred"], function(x) ifelse(x=="D",2,ifelse(x=="T",0,NA)))
				resultMCAP <- sum(M_CAP_pred)/length(which(!is.na(M_CAP_pred)))

			}		

		} else {

			resultMCAP <- NA

		}

	resultTOT <- c(resultSIFT,resultHDIV,resultHVAR,resultLRT,resultMT,resultMA,resultFATHMM,resultXF,resultMKL,resultPROV,resultSVM,
	resultLR,resultClin,resultDEO,resultLIST,resultRNN,result4G,resultES,resultMCAP)

	resultTOT <- resultTOT[!is.na(resultTOT)]

	SCORE[i] <- round(sum(resultTOT)/length(resultTOT),2)

	}
}

df$SCORE <- SCORE

write.table(df, fo, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)