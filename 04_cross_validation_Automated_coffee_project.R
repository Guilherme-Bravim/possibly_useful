#' @factors corresponde aos fatores a serem testados, pode ser informado 1 ou um vetor de numero diferentes de fatores. Os valores informados devem ser maior do que 1 
#' @phen.data Dados fenotipicos, ambientes devem ser codificados como env e genotipos como gen
#' @env.data Data frame contendo apenas as informacoes ambientais, inserindo os codigos dos ambientes como nome de linha
#' @n.comp.pls  corresponde ao numero de componentes da PLS a serem testados, pode ser informado 1 ou um vetor com diferentes valores. 
#' @description indica uma descricao para adicionar ao titulo do resultado a ser dalvo

cross.val <- function(factors=1, phen.data, env.data, n.comp.pls=1,description="Cross_validation"){
  
  for(k in factors){
    
    num.env<- nlevels(phen.data$env)
    num.gen<- nlevels(phen.data$gen)
    num.year<- nlevels(phen.data$ano)
    name.env<- levels(phen.data$env)
    name.gen<- levels(phen.data$gen)
    
# Modelando efeito aleatorio no FA    
    random_formula <- as.formula(
      paste("~ gen:fa(env,", k, ")")
    )
    
# Modelo 
    mod <- asreml(fixed = prod ~  env,
                  random = random_formula,
                  data = phen.data,
                  maxit = 1000,
                  na.action = na.method(x="include", y = "include"))
    mod<-update(mod)
  
# Extraindo info do modelo completo 
    library(tidyr)
    library(tibble)
    source("fa_outs_mod.R")

  fa.models = list(mod)
  fa.res = lapply(fa.models, function(x) fa.outs(x, name.env = "env",
                                                   name.gen = "gen"))
    
  # Extraindo os scores geneticos e cargas ambientais ----------------------------
  all <- fa.res[[1]]$blups
  fa = fa.res[[1]]
    
  #-------------------------------------------------------------------------------
  ######### CV using PLS
  #-------------------------------------------------------------------------------
  env.info<-env.data
  
  for(j in n.comp.pls){
    #EBLUPs_matrix <- matrix(NA, nrow = num.gen, ncol = num.env, dimnames = list(name.gen, name.env))
    EBLUPs_matrix<-data.frame(gen=unique(phen.data$gen))
    
    for (excluded_env in name.env) {
      
      # fenotypic data
      train_data <- subset(phen.data, env != excluded_env) %>%
        droplevels() %>%
        arrange(env)
      
      # enviromental data
      coffe <-as.data.frame(env.info) 
      train_covamb <- coffe[rownames(coffe) != excluded_env, ]
      
      mod1 <- asreml(fixed = prod ~ env,
                    random = random_formula,
                    data = train_data,
                    maxit = 1000,
                    na.action = na.method(x="include", y = "include"))
      mod1 <- update(mod1)
  
      fa.res1 <- fa.outs(mod1, name.env = "env", name.gen = "gen")
      
      data.pls.lamb <- data.frame(
        lambda = I(fa.res1$rot.loads),
        CovAmb = I(scale(as.matrix(train_covamb)))
      )
      
      pls_model <- pls::plsr(lambda ~ CovAmb, validation = "none", data = data.pls.lamb, ncomp = j)
      coef_pls <- coef(pls_model, intercept = TRUE)[,,1]
      
      coffe <- scale(coffe)
      matcovamb_pred <- coffe[rownames(coffe) == excluded_env, , drop = FALSE]
      

      fa_names <- paste0("fa", 1:k)
      pred_lambda <- matcovamb_pred %*% coef_pls[-1, fa_names, drop = FALSE] # Predição dos loadings
      pred_lambda <- sweep( pred_lambda, 2, coef_pls["(Intercept)", fa_names], "+" ) # Adiciona o intercepto
      
      
      EBLUPs_pred <- pred_lambda %*% t(fa.res1$rot.scores[, fa_names, drop = FALSE])
      
      df <- data.frame(gen=colnames(EBLUPs_pred), blups = as.vector(EBLUPs_pred))
      colnames(df)<-c("gen",excluded_env)
      EBLUPs_matrix<-merge(EBLUPs_matrix,df,by="gen",all.x = TRUE)
      
    }
    
    EBLUPs_matrix_df <- EBLUPs_matrix[, colnames(EBLUPs_matrix) != "gen"]
    rownames(EBLUPs_matrix_df)<-EBLUPs_matrix$gen
    
# Correlacao BLUP vs BLUP ------------------------------------------------------
    correlations <- list()
    ambientes <- unique(all$env)
    
    #ambientes<-unique(dt.test$env)
    for (env in ambientes) {
      subset_all <- all %>% filter(env == !!env)
      if (env %in% colnames(EBLUPs_matrix_df)) {
        common_gen <- intersect(subset_all$gen, rownames(EBLUPs_matrix_df))
        if (length(common_gen) > 1) {
          valores_all <- subset_all %>% filter(gen %in% common_gen) %>% arrange(gen) %>% pull(conditional)
          valores_eblup <- EBLUPs_matrix_df[common_gen, env, drop = FALSE] %>% arrange(rownames(.)) %>% pull()
          cor_value <- cor(valores_all, valores_eblup, method = "spearman",  use = "complete.obs")
          correlations[[env]] <- cor_value
        }
      }
    }
    
    correlation_df <- data.frame(
      Ambiente = names(correlations),
      Spearman_Correlation = unlist(correlations)
    )
    
    meanCV <- mean(correlation_df$Spearman_Correlation, na.rm = TRUE)
    
# Correlacao BLUP vs fenotipo --------------------------------------------------
    correlations2 <- list()

    #ambientes<-unique(dt.test$env)
    for (env in ambientes) {
      subset_all2 <- phen.data %>% data.frame()%>%filter(env == !!env)%>%select(env,gen,prod)
      
      if (env %in% colnames(EBLUPs_matrix_df)) {
        common_gen2 <- intersect(subset_all2$gen, rownames(EBLUPs_matrix_df))
        if (length(common_gen2) > 1) {
          valores.fen <- subset_all2 %>% filter(gen %in% common_gen2) %>% arrange(gen)
          valores_eblup2<-EBLUPs_matrix_df[common_gen2, env, drop = FALSE]
          valores_eblup2$gen<-rownames(valores_eblup2)
          valores.cor<-merge(valores_eblup2,valores.fen, by="gen")

          cor_value2 <- cor(valores.cor$prod,valores.cor[,env], method = "spearman",  use = "complete.obs")
          correlations2[[env]] <- cor_value2
        }
      }
    }
    
    correlation_df2 <- data.frame(
      Ambiente = names(correlations2),
      Spearman_Correlation = unlist(correlations2)
    )
    
    meanCV2 <- mean(correlation_df2$Spearman_Correlation, na.rm = TRUE)
    
    cross.validation<-data.frame(method=c("BLUP_BLUP","BLUP_fen"),cv.spearman=c(meanCV,meanCV2))
    
    write.table(cross.validation, file = paste0(description,"_FA",k,"_ncomp", j, ".txt"))
  }
    }
      }
  






