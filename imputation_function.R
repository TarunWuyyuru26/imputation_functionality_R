
rm(list = ls(all = TRUE))

Nexus_imputer <- function(data, num_entity = 'mean'){
  
  # Custom Functionalities
    # Multi Mode Function
      multi_mode <- function(ipvec, multimodal_op = FALSE, verbose = FALSE){
        
        ipvec = as.character(ipvec)
        freq_of_levels = table(ipvec)
        maxfreq_count = max(freq_of_levels)
        
        freq_of_levels = ifelse(freq_of_levels == maxfreq_count, TRUE, FALSE)
        
        multimodal_instances = rownames(which(x = freq_of_levels, arr.ind = TRUE))
        count_of_multimodal_instances = length(multimodal_instances)
        
        if(verbose)
          cat("No. of Instances for Mode is :", count_of_multimodal_instances, "\n")
        
        if(multimodal_op == TRUE){
          FinalMode = multimodal_instances
          
          if(verbose)
            cat("Unique Mode/s :", FinalMode, "\n")
          
          return(FinalMode)
          
        }else if(multimodal_op == FALSE){
          rndm_index = sample(x = 1:length(multimodal_instances), size = 1)
          FinalMode = multimodal_instances[rndm_index]
          
          if(verbose)
            cat("Multi Mode/s :", FinalMode, "\n")
          
          return(FinalMode)
        }
      }

    # Checking if atleast one NA is present in a Column or a Vector or complete DF
      atlstOneNA <- function(ip){
        NA_logical_Op = is.na(ip)
        if((TRUE %in% NA_logical_Op) == TRUE){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
  
  # Load basic Libraries
    
    
  # Main Process
    data = ipdata
    num_entity = 'mean'
    
    logical_op = sapply(data, is.numeric)
    num_cols = names(which(logical_op == TRUE))
    non_num_cols = names(which(logical_op != TRUE))
    
    num_df = data[,num_cols]
    non_num_df = data[,non_num_cols]
    
    if(num_entity == "mean"){
      num_cols_data = lapply(num_df, 
                             function(x){
                               temp = na.omit(x)
                               return(mean(temp))
                               })
      
      if(is.null(ncol(non_num_df)) != TRUE){
        non_num_data = lapply(non_num_df, 
                              function(x){
                                temp = na.omit(x)
                                multi_mode(ipvec = temp, 
                                           multimodal_op = FALSE)
                                })
      }else if(is.null(ncol(non_num_df)) == TRUE){
        temp_df = data.frame(non_num_df)
        colnames(temp_df) = non_num_cols
        non_num_data = lapply(na.omit(temp_df), function(x){multi_mode(ipvec = x, multimodal_op = FALSE)})
      }
      
      # Object
        object = list("numdata" = num_cols_data, "non_numdata" = non_num_data)
        
      # Imputing Process
        for (col in colnames(data)){
          if(atlstOneNA(data[,col]) == TRUE){
            if((class(data[,col]) == "numeric") | (class(data[,col]) == "integer")){
              NA_logical = match(x = data[,col], table = NA)
              NA_locs = which(NA_logical == 1)
              data[NA_locs,col] <- object$numdata[[col]]
            }else if((class(data[,col]) != "numeric") | (class(data[,col]) != "integer")){
              NA_logical = match(x = data[,col], table = NA)
              NA_locs = which(NA_logical == 1)
              data[NA_locs,col] <- object$non_numdata[[col]]
            }
          }else{
            next
          }
        }
        
    }else if(num_entity == "median"){
      num_cols_data = lapply(num_df, 
                             function(x){
                               temp = na.omit(x)
                               return(median(temp))
                             })
      
      if(is.null(ncol(non_num_df)) != TRUE){
        non_num_data = lapply(non_num_df, 
                              function(x){
                                temp = na.omit(x)
                                multi_mode(ipvec = temp, 
                                           multimodal_op = FALSE)
                              })
      }else if(is.null(ncol(non_num_df)) == TRUE){
        temp_df = data.frame(non_num_df)
        colnames(temp_df) = non_num_cols
        non_num_data = lapply(na.omit(temp_df), function(x){multi_mode(ipvec = x, multimodal_op = FALSE)})
      }
      
      # Object
        object = list("numdata" = num_cols_data, "non_numdata" = non_num_data)
      
      # Imputing Process
        for (col in colnames(data)){
          if(atlstOneNA(data[,col]) == TRUE){
            if((class(data[,col]) == "numeric") | (class(data[,col]) == "integer")){
              NA_logical = match(x = data[,col], table = NA)
              NA_locs = which(NA_logical == 1)
              data[NA_locs,col] <- object$numdata[[col]]
            }else if((class(data[,col]) != "numeric") | (class(data[,col]) != "integer")){
              NA_logical = match(x = data[,col], table = NA)
              NA_locs = which(NA_logical == 1)
              data[NA_locs,col] <- object$non_numdata[[col]]
            }
          }else{
            next
          }
        }
      
    }else if(num_entity == 'mice'){
      
      # loading required library
        suppressPackageStartupMessages(library(mice))
      
      # Imputation Process
        imputed_Data = mice(data = data, m=5, maxit = 50, method = 'pmm', seed = 500)
        data = mice::complete(data = imputed_Data)
        
    }else if(num_entity == 'amelia'){
      
      # loading required library
        suppressPackageStartupMessages(library(Amelia))
      
      # Imputation Process
        amelia_fit = amelia(data, m=5, parallel = "multicore", noms = non_num_cols[1])
        data = amelia_fit$imputations$imp1
      
    }else if(num_entity == 'missforest'){
      
      # loading required library
        suppressPackageStartupMessages(library(missForest))
      
      # Imputation Process
        missforest_impdf <- missForest(xmis = data, maxiter = 5, ntree = 10, 
                                     mtry = 3, replace = TRUE)
        data = missforest_impdf$ximp
      
    }else if(num_entity == 'hmisc'){
      
      # loading required library
        suppressPackageStartupMessages(library(Hmisc))
      
      # Imputation Process
        form_obj = as.formula(paste('', paste(colnames(data), collapse=" + "), sep=" ~ "))
        data <- aregImpute(formula = form_obj, data = data, 
                                 n.impute = 5, group = na.omit(non_num_df[,1]))
        
      
    }else if(num_entity == 'mi'){
      
      # loading required library
        suppressPackageStartupMessages(library(mi))
      
      # Imputation Process
        mi_data <- mi(data, seed = 335, parallel = TRUE)
        data = mi_data
    }
    
  # Returning the Final Imputed Dataframe
    return(data)
}

# Implementation

ipdata = read.csv('D:/iFiles/iris.csv', na.strings = c(''))

colSums(is.na(ipdata))

imp_data = Nexus_imputer(data = ipdata, num_entity = 'mean')
# imp_data = Nexus_imputer(data = ipdata, num_entity = 'median')
# imp_data = Nexus_imputer(data = ipdata, num_entity = 'mice')
# imp_data = Nexus_imputer(data = ipdata, num_entity = 'amelia')
# imp_data = Nexus_imputer(data = ipdata, num_entity = 'missforest')
# imp_data = Nexus_imputer(data = ipdata, num_entity = 'hmisc')
# imp_data = Nexus_imputer(data = ipdata, num_entity = 'mi')

colSums(is.na(imp_data))

View(imp_data)

