
#' @title get_envcov_otimizada
#' @description Versao otimizada em memoria de get_envcov(). Interpola
#'   covariaveis ambientais para pontos nao amostrados via RFSI (Random
#'   Forest Spatial Interpolation, Sekulic et al. 2020), evitando o
#'   colapso de memoria em bancos de predicao grandes (ex.: >100 mil
#'   linhas). Em vez de unir (merge) a tabela completa de predicao a
#'   cada variavel interpolada, os resultados sao acumulados em lista e
#'   unidos uma unica vez ao final, alem de liberar memoria (rm/gc) a
#'   cada variavel e salvar checkpoints periodicos em disco.
#'
#' @param predictors Character string com a formula dos preditores
#'   espaciais do modelo, no formato aceito por rfsi() (ex.:
#'   "~ LAT + LON"). Essas variaveis devem existir tanto em
#'   `data.training` quanto em `data.prediction`.
#' @param data.training Data frame com os dados de treinamento do
#'   modelo Random Forest. Deve conter observacoes completas da
#'   variavel resposta e dos preditores especificados em `predictors`,
#'   alem das colunas de coordenadas (`lat`, `lon`) e da variavel de
#'   agrupamento (`group.cv`).
#' @param data.prediction Data frame com os pontos onde a interpolacao
#'   sera realizada (ex.: grade de pixels do Brasil a 1 km2). As
#'   colunas anteriores a `from` sao tratadas como auxiliares
#'   (identificadores, coordenadas) e ignoradas na predicao; a partir
#'   de `from`, todas as colunas sao tratadas como variaveis alvo a
#'   serem interpoladas.
#' @param lat Character string com o nome da coluna de latitude em
#'   `data.training` e `data.prediction`. Padrao "LAT".
#' @param lon Character string com o nome da coluna de longitude em
#'   `data.training` e `data.prediction`. Padrao "LON".
#' @param n.trees Vetor numerico com os valores de numero de arvores
#'   (num.trees do ranger) testados na busca de hiperparametros via
#'   validacao cruzada. Padrao seq(300, 1000, by = 200).
#' @param n.neighbors Vetor numerico com os valores de numero de
#'   vizinhos (n.obs do rfsi) testados na busca de hiperparametros via
#'   validacao cruzada. Padrao c(3:10).
#' @param from Posicao (indice de coluna) a partir da qual
#'   `data.prediction` contem as variaveis alvo a serem interpoladas.
#'   Colunas antes desse indice sao tratadas como auxiliares. Padrao 4.
#' @param crs Sistema de referencia de coordenadas (Coordinate
#'   Reference System) usado para interpretar `lat`/`lon` como dados
#'   espaciais. Aceita codigo EPSG ou objeto CRS. Padrao EPSG:31982
#'   (SIRGAS 2000 / UTM zona 22S).
#' @param group.cv Character string com o nome da variavel usada para
#'   definir os grupos de validacao cruzada (leave-one-group-out), que
#'   separa dados de treino e validacao a cada rodada. Padrao "env".
#' @param salvar.a.cada Numero inteiro: a cada quantas variaveis
#'   processadas a funcao salva um checkpoint em disco
#'   (`arquivo.checkpoint`), permitindo retomar o processamento sem
#'   perder o progresso caso o R trave por falta de memoria. Padrao 20.
#' @param arquivo.checkpoint Character string com o caminho/nome do
#'   arquivo RDS usado para salvar e retomar o checkpoint de
#'   progresso. Padrao "checkpoint_envcov.RDS".
#'
#' @return Data frame `Z.final`, correspondente a `data.prediction`
#'   com todas as variaveis alvo (colunas a partir de `from`)
#'   completamente interpoladas (sem NA). O resultado tambem e salvo
#'   automaticamente em "new_data_untested_env.RDS".

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
library(dplyr)


get_envcov_otimizada <- function(predictors, data.training, data.prediction,
                                 lat = "LAT", lon = "LON",
                                 n.trees = seq(300, 1000, by = 200),
                                 n.neighbors = c(3:10),
                                 from = 4, crs = 31982, group.cv = "env",
                                 salvar.a.cada = 20,
                                 arquivo.checkpoint = "checkpoint_envcov.RDS") {
  
  vizinhos <- n.neighbors
  arvores  <- n.trees
  vars <- colnames(data.prediction)[-c(1:(from - 1))]
  w <- as.data.frame(data.training)
  
  w$group_cv <- w[[group.cv]]
  cv.factor  <- unique(w$group_cv)
  
  # ── conversao para sf feita UMA UNICA VEZ (nao a cada variavel) ──
  newdata.base <- data.prediction %>%
    st_as_sf(coords = c(lon, lat), crs = crs, remove = FALSE)
  coords.base <- st_coordinates(newdata.base)
  newdata.base[[lon]] <- coords.base[, 1]
  newdata.base[[lat]] <- coords.base[, 2]
  
  # guarda so lon/lat (para o join final) + lista de resultados por variavel
  chave.pred <- data.prediction[, c(lon, lat)]
  resultados.vars <- vector("list", length(vars))
  names(resultados.vars) <- vars
  
  # retomar de um checkpoint, se existir (ver Passo extra mais abaixo)
  inicio <- 1
  if (file.exists(arquivo.checkpoint)) {
    chk <- readRDS(arquivo.checkpoint)
    resultados.vars[names(chk$resultados.vars)] <- chk$resultados.vars
    inicio <- chk$proxima.variavel
    message("Retomando checkpoint a partir da variavel ", inicio, "/", length(vars))
  }
  
  for (h in inicio:length(vars)) {
    
    fm.RFSI  <- as.formula(paste(vars[h], predictors))
    variavel <- vars[h]
    message("Processando variavel ", h, "/", length(vars), ": ", variavel)
    
    results <- list()
    
    for (i in seq_along(vizinhos)) {
      for (j in seq_along(arvores)) {
        
        predito <- list()
        
        for (k in seq_along(cv.factor)) {
          
          dt.training <- w %>%
            filter(group_cv != cv.factor[k]) %>%
            st_as_sf(coords = c(lon, lat), crs = crs)
          
          coords_train <- st_coordinates(dt.training)
          dt.training[[lon]] <- coords_train[, 1]
          dt.training[[lat]] <- coords_train[, 2]
          
          newdata <- w %>% st_as_sf(coords = c(lon, lat), crs = crs)
          coords_pred <- st_coordinates(newdata)
          newdata[[lon]] <- coords_pred[, 1]
          newdata[[lat]] <- coords_pred[, 2]
          newdata[newdata$group_cv == cv.factor[k], variavel] <- NA
          
          rfsi_model <- rfsi(formula = fm.RFSI,
                             data = dt.training,
                             n.obs = vizinhos[i],
                             s.crs = st_crs(dt.training),
                             p.crs = st_crs(dt.training),
                             cpus = parallel::detectCores() - 1,
                             progress = FALSE,
                             importance = "impurity",
                             seed = 42,
                             num.trees = arvores[j],
                             mtry = 5,
                             splitrule = "variance",
                             min.node.size = 10,
                             sample.fraction = 0.95)
          
          rfsi_prediction <- pred.rfsi(model = rfsi_model,
                                       data = dt.training,
                                       obs.col = vars[h],
                                       newdata = newdata,
                                       output.format = "sf",
                                       zero.tol = 0,
                                       s.crs = st_crs(dt.training),
                                       newdata.s.crs = st_crs(dt.training),
                                       p.crs = st_crs(dt.training),
                                       cpus = 1,
                                       progress = FALSE)
          
          Local <- w %>% filter(group_cv == cv.factor[k])
          
          pred.value <- rfsi_prediction %>%
            st_drop_geometry() %>%
            bind_cols(st_coordinates(rfsi_prediction)) %>%
            filter(X %in% unique(Local$LON), Y %in% unique(Local$LAT))
          
          pred.value$real <- Local[[variavel]]
          predito[[length(predito) + 1]] <- data.frame(pred.value)
          
          # libera objetos pesados desta rodada de CV
          rm(rfsi_model, rfsi_prediction, dt.training, newdata, Local, pred.value)
        }
        
        cp  <- as.data.frame(do.call("rbind", predito))
        tab <- data.frame(N_vizinho = vizinhos[i], N_arvores = arvores[j],
                          cp = cor(cp$real, cp$pred))
        results[[length(results) + 1]] <- tab
        
        rm(predito, cp); gc(verbose = FALSE)
      }
    }
    
    resutados <- do.call("rbind", results) %>% arrange(desc(cp))
    
    # ── predicao final para esta variavel, no banco grande ──
    data.training2 <- w %>% st_as_sf(coords = c(lon, lat), crs = crs)
    coords_train2 <- st_coordinates(data.training2)
    data.training2[[lon]] <- coords_train2[, 1]
    data.training2[[lat]] <- coords_train2[, 2]
    
    model.train <- rfsi(formula = fm.RFSI,
                        data = data.training2,
                        n.obs = resutados[1, 1],
                        s.crs = st_crs(data.training2),
                        p.crs = st_crs(data.training2),
                        cpus = parallel::detectCores() - 1,
                        progress = FALSE,
                        importance = "impurity",
                        seed = 42,
                        num.trees = resutados[1, 2],
                        mtry = 5,
                        splitrule = "variance",
                        min.node.size = 10,
                        sample.fraction = 0.95)
    
    rfsi_prediction2 <- pred.rfsi(model = model.train,
                                  data = data.training2,
                                  obs.col = vars[h],
                                  newdata = newdata.base,   # reusa sf ja pronto
                                  output.format = "sf",
                                  zero.tol = 0,
                                  s.crs = st_crs(data.training2),
                                  newdata.s.crs = st_crs(data.training2),
                                  p.crs = st_crs(data.training2),
                                  cpus = 1,
                                  progress = FALSE)
    
    env.pred <- rfsi_prediction2 %>%
      st_drop_geometry() %>%
      bind_cols(st_coordinates(rfsi_prediction2)) %>%
      dplyr::select(!staid)
    colnames(env.pred) <- c(vars[h], lon, lat)
    
    # guarda so o resultado desta variavel (leve) -> NAO faz merge aqui
    resultados.vars[[variavel]] <- env.pred
    
    # libera tudo que foi usado nesta variavel
    rm(results, resutados, data.training2, model.train,
       rfsi_prediction2, env.pred)
    gc(verbose = FALSE)
    
    # checkpoint periodico -> se o R cair, voce retoma sem perder tudo
    if (h %% salvar.a.cada == 0 || h == length(vars)) {
      saveRDS(list(resultados.vars = resultados.vars, proxima.variavel = h + 1),
              file = arquivo.checkpoint)
      message("  Checkpoint salvo (", h, "/", length(vars), " variaveis)")
    }
  }
  
  # ── join final: UMA UNICA VEZ, no fim de tudo ──
  message("Unindo todas as variaveis interpoladas...")
  Z.final <- chave.pred
  for (variavel in vars) {
    if (is.null(resultados.vars[[variavel]])) next
    Z.final <- dplyr::left_join(Z.final, resultados.vars[[variavel]], by = c(lon, lat))
  }
  
  saveRDS(Z.final, file = "new_data_untested_env.RDS")
  
  # limpa o checkpoint depois de concluir com sucesso
  if (file.exists(arquivo.checkpoint)) file.remove(arquivo.checkpoint)
  
  return(Z.final)
}
