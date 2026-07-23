## Funcao 
fa.outs = function(model, name.env, name.gen){
  
  require(asreml)
  
  data = as.data.frame(model$mf)
  num.env = nlevels(data[, name.env])
  num.gen = nlevels(data[, name.gen])
  sum = summary(model)
  load = sum$varcomp[grep('fa\\d+', rownames(sum$varcomp)),]
  var = sum$varcomp[grep('var', rownames(sum$varcomp)), 1]
  
  # Loadings' rotation via SVD
  load$fa = regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings = do.call(cbind, lapply(split(load[,c('fa','component')], f = load$fa), 
                                       function(x) x[,'component']))
  rownames(mat.loadings) = levels(data[, name.env])
  
  svd.lambda = svd(mat.loadings)
  D = diag(svd.lambda$d^2, nrow = length(svd.lambda$d))
  if(sum(svd.lambda$u[,1] < 0)/nrow(svd.lambda$u) > 0.5){
    svd.lambda$u = -1 * svd.lambda$u
    svd.lambda$v = -1 * svd.lambda$v
  }
  
  mat.loadings.star = svd.lambda$u
  colnames(mat.loadings.star) = colnames(mat.loadings)
  rownames(mat.loadings.star) = rownames(mat.loadings)
  
  # Scores rotation via SVD
  scores = data.frame(summary(model, coef = T)$coef.random)[
    grep('Comp', rownames(data.frame(summary(model, coef = T)$coef.random))),
  ]
  scores$fa = regmatches(rownames(scores), regexpr("Comp\\d+", rownames(scores)))
  scor.vec = do.call(c, lapply(split(scores, f = scores$fa), function(x) x[,'solution']))
  scor.vec.star = kronecker(sqrt(D) %*% t(svd.lambda$v), diag(num.gen)) %*% scor.vec
  scor.mat.star = matrix(scor.vec.star, nrow = num.gen, ncol = length(unique(scores$fa)),
                         dimnames = list(levels(data[, name.gen]), unique(load$fa)))
  
  
  # Full variance-covariance genetic matrix and genetic correlation matrix
  Gvcov = mat.loadings.star %*% tcrossprod(D, mat.loadings.star) + diag(var)
  Gcor = cov2cor(Gvcov)
  
  # Percentage of explained variance
  expvar = sum(diag(mat.loadings.star %*% tcrossprod(D, mat.loadings.star)))/
    sum(diag(Gvcov))
  expvar.j = matrix(nrow = ncol(mat.loadings), ncol = num.env, 
                    dimnames = list(colnames(mat.loadings), levels(data[, name.env])))
  for (i in 1:nrow(expvar.j)) {
    for (j in 1:ncol(expvar.j)) {
      expvar.j[i,j] = 100 * mat.loadings.star[j,i]^2 * diag(D)[i]/
        (sum(mat.loadings.star[j,]^2 * diag(D)) + var[j])
    }
  }
  
  # Average semivariances ratio
  lambdacross = mat.loadings.star %*% tcrossprod(D, mat.loadings.star)
  fullsv = list()
  for (i in levels(data[,name.env])) {
    semivar = matrix(nrow = num.env, ncol = 3,
                     dimnames = list(levels(data[,name.env]),
                                     c('i', 'j', 'semivar')))
    for (j in levels(data[,name.env])) {
      semivar[j,] = c(i, j, .5 * (lambdacross[i,i] + lambdacross[j,j]) -
                        lambdacross[i,j])
    }
    fullsv[[i]] = semivar
    rm(semivar)
  }
  fullsv = data.frame(do.call(rbind, fullsv))
  fullsv$semivar = as.numeric(fullsv$semivar)
  fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[,-1]
  colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv)- 1)) * sum(fullsv[upper.tri(fullsv)])
  rm(fullsv)
  fullsv = list()
  for (i in levels(data[,name.env])) {
    semivar = matrix(nrow = num.env, ncol = 3,
                     dimnames = list(levels(data[,name.env]),
                                     c('i', 'j', 'semivar')))
    for (j in levels(data[,name.env])) {
      semivar[j,] = c(i,j,
                      0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
    }
    fullsv[[i]] = semivar
  }
  fullsv = data.frame(do.call(rbind, fullsv))
  fullsv$semivar = as.numeric(fullsv$semivar)
  fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[,-1]
  colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  G_ASV = 2/(nrow(fullsv) * (nrow(fullsv)- 1)) * sum(fullsv[upper.tri(fullsv)])
  ASVR = lambda_ASV/G_ASV
  
  # eBLUPs
  num.factors<-as.numeric(sub(".*?,\\s*([^)]*)\\).*", "\\1", names(model$vparameters)[1]))
  
  comb <- expand.grid(gen = unique(data[,name.gen]), env = unique(data[,name.env]))
  selecionar<-paste0(name.gen,"_",comb$gen,":fa(",name.env,", ", num.factors ,")_",comb$env)
  
  blups<-as.data.frame(coef(model)$random)  
  blups$cod<-rownames(blups)
  
  blups<-blups%>%filter(cod%in%selecionar)
  blups[[name.gen]]<-sub(paste0("^",name.gen,"_"), "", sub(":.*", "", blups$cod))
  blups[[name.env]]<-sub(".*_", "", blups$cod)
  
  temp = data.frame(
    name.env = rep(levels(data[, name.env]), each = nlevels(data[, name.gen])),
    name.gen = rep(levels(data[, name.gen]), times = nlevels(data[, name.env])), 
    marginal = kronecker(mat.loadings.star, diag(num.gen)) %*% scor.vec.star
  )
  
  colnames(temp)<-c(name.env, name.gen,"marginal")
  blups<-merge(blups,temp,by=c(name.env, name.gen))
  
  colnames(blups)[which(colnames(blups) == 'effect')] = 'conditional'
  rownames(blups)<-NULL
  
  blups<-blups[,c(name.env, name.gen,"conditional","marginal")]
  
  # Environment-wise generalized heritabilities
  
  results = list('rot.loads' = mat.loadings.star, 
                 'Gvcov' = Gvcov, 
                 'Gcor' = Gcor, 
                 'diagnostics' = c(
                   expvar = round(expvar*100, 3),
                   ASVR = round(ASVR * 100, 3),
                   aic = round(sum$aic, 3), 
                   bic = round(sum$bic, 3),
                   logl = round(sum$loglik, 3)
                 ),
                 'expvar_j' = expvar.j,
                 "rot.scores" = scor.mat.star, 
                 'blups' = blups)
  
  return(results)
}