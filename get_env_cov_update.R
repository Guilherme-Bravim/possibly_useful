
#' @title get_envcov_otimizada
#' @description Versao otimizada em memoria de get_envcov(). Interpola
#'   covariaveis ambientais para pontos nao amostrados via RFSI (Random
#'   Forest Spatial Interpolation, Sekulic et al. 2020), evitando o
#'   colapso de memoria em bancos de predicao grandes (ex.: >100 mil
#'   linhas). Os resultados de cada variavel sao gravados diretamente
#'   em uma matriz pre-alocada de tamanho fixo (nunca em uma lista
#'   crescente nem via merge/join repetido), alem de liberar memoria
#'   (rm/gc) a cada variavel e salvar checkpoints periodicos em disco.
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
#' @param usar.inteiro.escalado Logico. Se TRUE, os valores
#'   interpolados sao arredondados para `casas.decimais` casas e
#'   armazenados como inteiro multiplicado por 10^casas.decimais (ex.:
#'   25.437 com 2 casas vira 2544L), ocupando 4 bytes/valor em vez dos
#'   8 bytes/valor de um double -- reduzindo pela metade a memoria e o
#'   tamanho em disco do resultado final. Para reverter aos valores
#'   originais, dividir por 10^casas.decimais (o fator usado fica
#'   salvo no atributo "fator.escala" do objeto retornado). Padrao
#'   FALSE (mantem double, comportamento identico ao original).
#' @param casas.decimais Numero inteiro de casas decimais de precisao
#'   desejadas. Usado para arredondar o resultado final quando
#'   `usar.inteiro.escalado = FALSE`, e para definir o fator de escala
#'   quando `usar.inteiro.escalado = TRUE`. Padrao 2.
#'
#' @return Data frame `Z.final`, correspondente a `data.prediction`
#'   com todas as variaveis alvo (colunas a partir de `from`)
#'   completamente interpoladas (sem NA). O resultado tambem e salvo
#'   automaticamente em "new_data_untested_env.RDS". Se
#'   `usar.inteiro.escalado = TRUE`, os valores estao em inteiro
#'   escalado (ver `usar.inteiro.escalado` acima).

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
                                 arquivo.checkpoint = "checkpoint_envcov.RDS",
                                 usar.inteiro.escalado = FALSE,
                                 casas.decimais = 2) {
  
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
  
  # chave de juncao ESTAVEL: nao confiar na ordem devolvida por pred.rfsi;
  # cada variavel eh casada de volta por lon/lat via match(), nunca por
  # posicao de linha
  chave.pred <- data.prediction[, c(lon, lat)]
  chave.id   <- paste(chave.pred[[lon]], chave.pred[[lat]], sep = "_")
  n.pred     <- nrow(chave.pred)
  
  # ── UNICA estrutura que "cresce": uma matriz pre-alocada de tamanho
  # fixo. Substitui o acumulo de ~869 data.frames (491465 linhas cada)
  # que antes eram somente unidos (left_join) no final -> aquele padrao
  # ainda mantinha tudo em RAM simultaneamente antes do join, e por isso
  # continuava estourando memoria em bancos com muitas variaveis.
  resultado.matriz <- matrix(
    if (usar.inteiro.escalado) NA_integer_ else NA_real_,
    nrow = n.pred, ncol = length(vars),
    dimnames = list(NULL, vars)
  )
  fator.escala <- 10^casas.decimais  # usado somente se usar.inteiro.escalado = TRUE
  
  # retomar de um checkpoint, se existir
  inicio <- 1
  if (file.exists(arquivo.checkpoint)) {
    chk <- readRDS(arquivo.checkpoint)
    resultado.matriz[, colnames(chk$resultado.matriz)] <- chk$resultado.matriz
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
    
    # ── grava direto na coluna da matriz, casando por lon/lat ──
    # match() garante a linha certa mesmo se pred.rfsi devolver em
    # ordem diferente da original. Nada e acumulado em lista.
    id.pred    <- paste(env.pred[[lon]], env.pred[[lat]], sep = "_")
    pos.destino <- match(id.pred, chave.id)
    valores <- env.pred[[vars[h]]]
    if (usar.inteiro.escalado) {
      valores <- as.integer(round(valores * fator.escala))
    }
    resultado.matriz[pos.destino, variavel] <- valores
    
    # libera tudo que foi usado nesta variavel
    rm(results, resutados, data.training2, model.train,
       rfsi_prediction2, env.pred, id.pred, pos.destino, valores)
    gc(verbose = FALSE)
    
    # checkpoint periodico -> se o R cair, voce retoma sem perder tudo
    # (salva so a matriz, nao mais uma lista de data.frames)
    if (h %% salvar.a.cada == 0 || h == length(vars)) {
      saveRDS(list(resultado.matriz = resultado.matriz, proxima.variavel = h + 1),
              file = arquivo.checkpoint)
      message("  Checkpoint salvo (", h, "/", length(vars), " variaveis)")
    }
  }
  
  # ── monta o resultado final SEM nenhum merge/join ──
  # cbind direto: chave.pred e resultado.matriz ja estao na mesma ordem
  # de linhas (n.pred), entao nao ha custo extra de juncao aqui.
  message("Montando resultado final...")
  if (usar.inteiro.escalado) {
    # mantem como inteiro escalado no RDS (4 bytes/valor, metade do double)
    # para reverter: as.data.frame(resultado.matriz) / 10^casas.decimais
    Z.final <- cbind(chave.pred, as.data.frame(resultado.matriz))
    attr(Z.final, "fator.escala") <- fator.escala
    message("  Valores salvos como inteiro escalado (fator = ", fator.escala,
            "). Para reverter: valores / ", fator.escala)
  } else {
    resultado.matriz <- round(resultado.matriz, casas.decimais)
    Z.final <- cbind(chave.pred, as.data.frame(resultado.matriz))
  }
  
  saveRDS(Z.final, file = "new_data_untested_env.RDS")
  
  # limpa o checkpoint depois de concluir com sucesso
  if (file.exists(arquivo.checkpoint)) file.remove(arquivo.checkpoint)
  
  return(Z.final)
}

