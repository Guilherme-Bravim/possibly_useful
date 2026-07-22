#' @param predictors A character vector containing the names of the predictor variables used by the model. These variables must be available in both `data.training` and `data.prediction`. An exemple: predictors<-"~ LAT + LON"
#' @param data.training A data frame containing the training data used to fit the Random Forest model. It must include complete observations for the response variable and all predictor variables specified in `formula`. An example dataset is available through `get_data.training_example()`.
#' @param data.prediction A data frame containing the data used for prediction. It must include the predictor variables required by the model, with properly named columns. Columns preceding `from` may contain auxiliary information (e.g., identifiers or coordinates) and are ignored during prediction. The `from` argument specifies the first column containing the variables to be predicted, and all columns from that position onward are treated as target variables. The names of all predictor and target variable columns must match the corresponding column names in `data.training`.
#' @param from argument specifies the first column containing the predictor variables, and all columns from that position onward are used for prediction. By structure, start in 4.
#' @param lat Character string specifying the name of the latitude column in `data.prediction` and `param data.training`. This variable is required and must be present in the input data to perform spatial predictions. The default column name is `"LAT"`.
#' @param lon Character string specifying the name of the longitude column in `data.prediction` and `param data.training`. This variable is required and must be present in the input data to perform spatial predictions. The default column name is `"LON"`.
#' @param crs Coordinate Reference System (CRS) used to define the spatial reference of the input coordinates. The CRS specifies how geographic coordinates are represented and allows spatial data to be correctly interpreted and combined with other spatial datasets. It can be provided as an EPSG code or as a CRS object. The default is EPSG:31982 (SIRGAS 2000 / UTM zone 22S).
#' @param group.cv Character string specifying the grouping variable used for cross-validation to estimate the predictive performance of the model. The grouping variable defines the units that will be used to separate training and validation datasets. By default, the variable `"YYYMMDD"` is used, representing the sampling date. Users may specify any other variable of interest, such as day, environment, location, experimental group or combination between factors.
#' 

message("Function developed for interpolation of environmental data, based on the RFSI (Random Forest Spatial Interpolation) methodology (Sekulic et al. 2020, https://github.com/AleksandarSekulic/RFSI). 
        
The function automatically selects the best number of trees and neighbors for the interpolation of each variable.")

library(meteo)
library(sp)
library(sf)
library(sftime)
library(terra)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)
library(ranger)


get_data.training_example <- function() {
  
  mat_w_exemplo <- matrix(
    c(1, -20.75, -42.87, 25.4, 35.8,
      2, -21.10, -43.20, 24.1, 42.1,
      3, -20.95, -42.95, 26.0, 39.7),
    nrow = 3, byrow = TRUE,
    dimnames = list(NULL, c("env", "lat", "lon", "temp", "temp.max"))
  )
  
  message("Exemplo de estrutura esperada para o argumento 'data' de get_envcov():")
  print(mat_w_exemplo)
  
  return(invisible(mat_w_exemplo))
}


get_envcov <- function(predictors, data.training, data.prediction, lat="LAT", lon="LON",
                       n.trees=seq(300,1000,by=200), n.neighbors=c(3:10),
                       from=4, crs = 31982, group.cv="env") {

  vizinhos<- n.neighbors 
  arvores<-n.trees
  vars<-colnames(data.prediction)[-c(1:(from-1))]
  w<-as.data.frame(data.training)
  
  Z<-data.prediction
  
  w$group_cv<- w[[group.cv]]
  cv.factor<-unique(w$group_cv)
  
  results<-list()
  for(h in 1 : length(vars)){
    
    fm.RFSI <- as.formula(paste(vars[h], predictors))
    variavel<-vars[h]
    
    for(i in 1:length(vizinhos)){
      for(j in 1:length(arvores)){
        
        predito<-list()
        
        for(k in 1:length(cv.factor)) { # preparando o arquivo de dados
          
          # Preparando dados treinamento ------------    
          dt.training <- w %>%
            filter(group_cv != cv.factor[k]) %>%
            st_as_sf(coords = c(lon, lat), crs = crs)
            
          # recriando coordenadas  
          coords_train <- st_coordinates(dt.training)
      
          dt.training[[lon]] <- coords_train[, 1]
          dt.training[[lat]] <- coords_train[, 2]

          # Preparando dados de predicao ------------
          newdata <- w %>% st_as_sf(coords = c(lon, lat), crs = crs)
          
          coords_pred <- st_coordinates(newdata)
          newdata[[lon]] <- coords_pred[, 1]
          newdata[[lat]] <- coords_pred[, 2]
          
          newdata[newdata$group_cv == cv.factor[k], variavel] <- NA

          # Modelo de treinamento -----------
          rfsi_model <- rfsi(formula = fm.RFSI,
                             data = dt.training,
                             n.obs = vizinhos[i],
                             s.crs = st_crs(dt.training),
                             p.crs = st_crs(importance = "impurity",
                                            seed = 42,dt.training),
                             cpus = parallel::detectCores() - 1,
                             progress = TRUE,
                             importance = "impurity",
                             seed = 42,
                             num.trees = arvores[j],
                             mtry = 5,
                             splitrule = "variance",
                             min.node.size = 10,
                             sample.fraction = 0.95)
          
          # Predicao ---------------------------------
          rfsi_prediction <- pred.rfsi(model = rfsi_model,
                                       data = dt.training,
                                       obs.col = vars[h],
                                       newdata = newdata,
                                       output.format = "sf",     # pontos, nao raster
                                       zero.tol = 0,
                                       s.crs = st_crs(dt.training),
                                       newdata.s.crs = st_crs(dt.training),
                                       p.crs = st_crs(dt.training),
                                       cpus = 1,
                                       progress = TRUE
          )
          

          Local<-w%>%filter(group_cv == cv.factor[k])
          
          pred.value <- rfsi_prediction %>%
            st_drop_geometry() %>%
            bind_cols(st_coordinates(rfsi_prediction)) %>%   
            filter(X %in% unique(Local$LON), Y %in% unique(Local$LAT))
          
          pred.value$real<-Local[[variavel]]

          predito[[length(predito)+1]]<-data.frame(pred.value)  
          
      }# CP 
        
        cp<-as.data.frame(do.call("rbind",predito))
        tab<-data.frame(N_vizinho=vizinhos[i],
                        N_arvores=arvores[j],
                        cp=cor(cp$real,cp$pred))
        
        results[[length(results)+1]]<-tab  
        
      }  # n arvores
    }# n vizinhos  
    
    resutados<-do.call("rbind",results)%>%arrange(desc(cp))
    
    #---------------------------------------------------  
    # FAZENDO AS PREDICOES <----------------------------
    #---------------------------------------------------  
    
    # Preparando dados treinamento ------------    
    data.training2 <- w %>%
      st_as_sf(coords = c(lon, lat), crs = crs)
    
    coords_train2 <- st_coordinates(data.training2)
    data.training2[[lon]] <- coords_train2[,1]
    data.training2[[lat]] <- coords_train2[,2]
    
    model.train<-rfsi(formula = fm.RFSI,
                      data = data.training2,
                      n.obs = resutados[1,1],
                      s.crs = st_crs(data.training2),
                      p.crs = st_crs(data.training2),
                      cpus = parallel::detectCores() - 1,
                      progress = TRUE,
                      importance = "impurity",
                      seed = 42,
                      num.trees = resutados[1,2],
                      mtry = 5,
                      splitrule = "variance",
                      min.node.size = 10,
                      sample.fraction = 0.95)

    newdata <- Z%>%
      st_as_sf(coords = c(lon, lat), crs = crs)
    
    # recriando coordenadas  
    coords<- st_coordinates(newdata)
    newdata[[lon]] <- coords[,1]
    newdata[[lat]] <- coords[,2]

    rfsi_prediction2 <- pred.rfsi(model = model.train,
                                 data = data.training2,
                                 obs.col = vars[h],
                                 newdata = newdata,
                                 output.format = "sf",   # <- pontos, nao SpatRaster
                                 zero.tol = 0,
                                 s.crs = st_crs(data.training2),
                                 newdata.s.crs = st_crs(data.training2),
                                 p.crs = st_crs(data.training2),
                                 cpus = 1,
                                 progress = TRUE)
    
    env.pred<-as.data.frame(rfsi_prediction,xy=TRUE)
    
    env.pred <- rfsi_prediction2 %>%
      st_drop_geometry() %>%
      bind_cols(st_coordinates(rfsi_prediction2))%>%
      dplyr::select(!staid)
    
    colnames(env.pred)<-c( vars[h], lon, lat)
    env.pred<-as.data.frame(env.pred)
    
    Z <- Z %>%
      dplyr::select(-all_of(vars[h]))
    
    Z<-merge(Z,env.pred, by = c(lon, lat))
    
  }
  
  saveRDS(Z, file="new_data_untested_env.RDS")
  return(Z)
  
  saveRDS(Z, file= paste0("interpolated_data_", format(Sys.time(), "%Y-%m-%d_%H-%M"), ".RDS"))
  
} # Fecha a funcao
























head(rfsi_prediction$staid)
head(w$id)
