# inspired by JM:::survfitJMCR.jointModel
predict.jointModel.CR <- 
  function (object, newData.long, newData.surv, idVar = "id", survTimes = NULL, 
            last.time = NULL, simulate = F, M = 1000, CI.levels = c(0.025, 0.975), formT, 
            GHk = NULL, GKk = NULL, ...) {
    # Now, we have to specify the last.time for each subject (as vector)
    # formT is the formula in the survival part without the Surv part and without the stratification,
    # i.e. without '+ strata(.)'.
    # Dans le code, on suppose que les ID sont triés de la même manière dans newData.long et newData.surv
    # il faut également vérifier que les id dans les deux tables correspondent
    # vérifier que le tps landmark ne dépasse pas le tps de dernière mesure ds newData.long
    # "GHk" argument: - number of GH quad points to compute the integral over the RE
    #                 -> added: if unspecified, it will used by default object$control$GHk
    # Même chose avec GKk.
    # Il faudrait que la longueur de survTimes puisse être > 1 , pour l'instant c'est impossible
    if (!inherits(object, "jointModel")) 
      stop("Use only with 'jointModel' objects.\n")
    if (!is.data.frame(newData.long) || nrow(newData.long) == 0 ||
          !is.data.frame(newData.surv) || nrow(newData.surv) == 0) 
      stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newData.long[[idVar]]) || is.null(newData.surv[[idVar]])) 
      stop("'idVar' not in 'newData.long.\n'")
    if (is.null(survTimes) || !is.numeric(survTimes)) 
      stop("''survTimes' has to be numeric.\n") # added
    if (object$CompRisk != TRUE)
      stop("The object must be a joint model with competing risks.\n'") # voir si on est obligé d'avoir cette condition
    if (object$method != "spline-PH-GH")
      stop("How did you do to fit a joint model with competing risks and method != 'spline-PH-GH' ?\n")
    if (!identical(unique(newData.long[[idVar]]), unique(newData.surv[[idVar]])))
      stop("The id in the longitudinal and the survival datasets are not the same.\n")
    if (object$LongFormat)
      stop("This function is not implemented for Longformat datasets.\n")
    if (is.null(GHk)) # added (nb of quad pts for the integrals over the RE)
      GHk <- object$control$GHk # added
    if (is.null(GKk)) # added (nb of quad pts for the integrals over the time)
      GKk <- object$control$GKk # added
    method <- object$method
    timeVar <- object$timeVar
    interFact <- object$interFact
    parameterization <- object$parameterization
    derivForm <- object$derivForm
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
    TermsX <- object$termsYx
    TermsZ <- object$termsYz
    TermsX.deriv <- object$termsYx.deriv
    TermsZ.deriv <- object$termsYz.deriv
    mfX <- model.frame(TermsX, data = newData.long)
    mfZ <- model.frame(TermsZ, data = newData.long)
    formYx <- formula(delete.response(TermsX))
    formYz <- object$formYz
    na.ind <- as.vector(attr(mfX, "na.action"))
    na.ind <- if (is.null(na.ind)) {
      rep(TRUE, nrow(newData.long))
    }
    else {
      !seq_len(nrow(newData.long)) %in% na.ind
    }
    id <- as.numeric(unclass(newData.long[[idVar]]))
    id <- id. <- match(id, unique(id))
    id <- id[na.ind]
    y <- model.response(mfX)
    X <- model.matrix(formYx, mfX)
    Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
    TermsT <- object$termsT
    data.id <- newData.long[tapply(row.names(newData.long), id, tail, n = 1L), ]
    idT <- newData.surv[[idVar]]
    idT <- match(idT, unique(idT))
    kk <- attr(TermsT, "specials")$strata
    strt <- eval(attr(TermsT, "variables"), newData.surv)[[kk]] # vérifier quand on inverse l'ordre des colonnes (t_Death avant t_Rec)
    W <- model.matrix(formT, data = newData.surv)
    WintF.vl <- WintF.sl <- as.matrix(rep(1, nrow(newData.surv))) # Idem: à construire moi-même
    if (!is.null(interFact)) {
      if (!is.null(interFact$value)) 
        WintF.vl <- model.matrix(interFact$value, data = newData.surv)
      if (!is.null(interFact$slope)) 
        WintF.sl <- model.matrix(interFact$slope, data = newData.surv)
    }
    obs.times <- split(newData.long[[timeVar]][na.ind], id)
    if (!is.numeric(last.time) || length(last.time) != nrow(data.id))
      stop("\nnot appropriate value for 'last.time' argument.")
    times.to.pred <- lapply(last.time, function(t) survTimes[survTimes > t])
    n <- object$n
    n.tp <- length(last.time) # nb of ID to predict
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    lag <- object$y$lag
    betas <- object$coefficients[["betas"]]
    sigma <- object$coefficients[["sigma"]]
    D <- object$coefficients[["D"]]
    diag.D <- ncol(D) == 1 & nrow(D) > 1
    if (diag.D) 
      D <- diag(c(D))
    gammas <- object$coefficients[["gammas"]]
    alpha <- object$coefficients[["alpha"]]
    Dalpha <- object$coefficients[["Dalpha"]]
    gammas.bs <- object$coefficients[["gammas.bs"]]
    list.thetas <- list(betas = betas, log.sigma = log(sigma), 
                        gammas = gammas, alpha = alpha, Dalpha = Dalpha, 
                        gammas.bs = gammas.bs, 
                        D = if (diag.D) log(diag(D)) else JM:::chol.transf(D))
    if (ncww == 1) {
      W <- NULL
      ncww <- 0
    }
    else {
      W <- W[, -1, drop = FALSE]
      ncww <- ncww - 1
    }
    Q <- object$x$Q
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    Var.thetas <- vcov(object)
    obs.times.surv <- split(data.id[[timeVar]], unique(idT))
    K <- nrow(W) / n.tp # number of causes of event
    survMats <- survMats.last <- vector("list", n.tp)
    
    environment(JMCR.ModelMats) <- environment(JMCR.log.posterior.b) <- environment(JMCR.log.posterior2.b) <- 
      environment(JMCR.ModelMats) <- environment(JMCR.ModelMats.st) <- 
      environment(JMCR.ModelMats.hazard.st) <- environment(JMCR.S.b.st) <- 
      environment(JMCR.hazard.b.st) <- environment(S.pred.punct) <- 
      environment(S.last.punct) <- environment()
    
    for (i in seq_len(n.tp)) {
      # Complete data (and design matrices) at the prediction time
      # survMats[[i]] <- lapply(times.to.pred[[i]], JMCR.ModelMats, ii = i, K = K)
      # Complete data (and design matrices) at the landmark time
      survMats.last[[i]] <- JMCR.ModelMats(last.time[i], ii = i, K = K)
    }
    
    # Computation of the posterior mode of the RE and the associated variance
    #  useful for the computation of the integral over the RE
    modes.b <- matrix(0, n.tp, ncz) # pour l'intégrale du dessous
    Vars.b <- vector("list", n.tp)
    ff <- function(b, y, tt, mm, i) -JMCR.log.posterior.b(b, y, Mats = tt, ii = i)
    for (i in seq_len(n.tp)) {
      opt <- try(optim(rep(0, ncz), ff, y = y, tt = survMats.last, 
                       i = i, method = "BFGS", hessian = TRUE), 
                 TRUE)
      if (inherits(opt, "try-error")) {
        gg <- function(b, y, tt, mm, i) JM:::cd(b, ff, y = y, tt = tt, i = i)
        opt <- optim(rep(0, ncz), ff, gg, y = y, tt = survMats.last, 
                     i = i, method = "BFGS", hessian = TRUE)
      }
      modes.b[i, ] <- opt$par
      Vars.b[[i]] <- solve(opt$hessian)
    } # environ 50ms par indiv
    
    # ok et vérifié jusqu'ici (04/03/16)
    
    modes2.b <- replicate(K, matrix(0, n.tp, ncz), simplify = FALSE) # pour l'intégrale du dessus
    Vars2.b <- replicate(K, vector("list", n.tp), simplify = FALSE)
    ff2 <- function(b, y, tt, mm, i, k) { environment(); -JMCR.log.posterior2.b(b, y, ii = i, k)}
    for (i in seq_len(n.tp)) { ## IL FAUDRAIT PARALLELISER ICI!
      for (k in seq_len(K)) { # on met une boucle sur les causes d'évt
        opt2 <- try(optim(modes.b[i, ], ff2, y = y,
                          i = i, k = k, method = "BFGS", hessian = TRUE), 
                    TRUE) # environ 0.5s par indiv et par cause d'évt
        if (inherits(opt2, "try-error")) {
          gg2 <- function(b, y, tt, mm, i) JM:::cd(b, ff, y = y, i = i)
          opt2 <- optim(modes.b[i, ], ff2, gg2, y = y,
                        i = i, k = k, method = "BFGS", hessian = TRUE)
        }
        modes2.b[[k]][i, ] <- opt2$par
        Vars2.b[[k]][[i]] <- solve(opt2$hessian)
      }
    }
    
    res <- replicate(K, vector("list", n.tp), simplify = FALSE)
    GH <- JM:::gauher(GHk)
    b <- as.matrix(expand.grid(rep(list(GH$x), ncz)))
    wGH <- as.matrix(expand.grid(rep(list(GH$w), ncz)))
    wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
    for (i in seq_len(n.tp)) {
      VCdets.last <- Vars.b[[i]]
      b.last <- rep(modes.b[i, ], each = nrow(b)) + sqrt(2) * t(solve(chol(solve(VCdets.last))) %*% t(b))
      S.last <- c(exp(apply(b.last, MARGIN = 1, 
                            FUN = function(x) JMCR.log.posterior.b(x, y, survMats.last, ii = i))) * 
                    1/det(chol(solve(VCdets.last)))) %*% wGH
      for (k in seq_len(K)) { # on met une boucle sur les causes d'évt
        VCdets.pred <- Vars2.b[[k]][[i]]
        b.pred <- rep(modes2.b[[k]][i, ], each = nrow(b)) + sqrt(2) * t(solve(chol(solve(VCdets.pred))) %*% t(b))
        S.pred <- c(exp(apply(b.pred, MARGIN = 1,
                              FUN = function(x) JMCR.log.posterior2.b(x, y, ii = i, k))) *
                      1/det(chol(solve(VCdets.pred)))) %*% wGH
        res[[k]][[i]] <- cbind(times = times.to.pred[[i]], predSurv = S.pred/S.last)
        rownames(res[[k]][[i]]) <- seq_along(S.pred)
      }
    }
    
    res.b <- replicate(K, vector("list", n.tp), simplify = FALSE)
    for (i in seq_len(n.tp)) {
      S.last.b <- S.last.punct(modes.b[i, ], survMats.last, ii = i)
      for (k in seq_len(K)) {
        S.pred.b <- S.pred.punct(modes2.b[[k]][i, ], ii = i, k) #modifié le 13/10/2016 (modes.b en modes2.b)
        res.b[[k]][[i]] <- cbind(times = times.to.pred[[i]], predSurv = S.pred.b/S.last.b)
        rownames(res.b[[k]][[i]]) <- seq_along(S.pred.b)
      }
    }
    
    if (simulate) {
      thetas.new <- mvrnorm(n = M, mu = thetas, Sigma = Var.thetas)
      res.m <- res.m.b <- replicate(n.tp, matrix(, nrow = M, ncol = K), simplify = F)
      for(m in 1:M){
        thetas.m <- relist(thetas.new[m,], skeleton = list.thetas)
        betas <- thetas.m$betas
        sigma <- exp(thetas.m$log.sigma)
        gammas <- thetas.m$gammas
        alpha <- thetas.m$alpha
        Dalpha <- thetas.m$Dalpha
        D <- thetas.m$D
        D <- if (diag.D) 
          exp(D)
        else JM:::chol.transf(D)
        gammas.bs <- thetas.m$gammas.bs
        
        modes.b <- matrix(0, n.tp, ncz)
        Vars.b <- vector("list", n.tp)
        ff <- function(b, y, tt, mm, i) -JMCR.log.posterior.b(b, y, Mats = tt, ii = i)
        for (i in seq_len(n.tp)) {
          opt <- try(optim(rep(0, ncz), ff, y = y, tt = survMats.last, 
                           i = i, method = "BFGS", hessian = TRUE), 
                     TRUE)
          if (inherits(opt, "try-error")) {
            gg <- function(b, y, tt, mm, i) JM:::cd(b, ff, y = y, tt = tt, i = i)
            opt <- optim(rep(0, ncz), ff, gg, y = y, tt = survMats.last, 
                         i = i, method = "BFGS", hessian = TRUE)
          }
          modes.b[i, ] <- opt$par
          Vars.b[[i]] <- solve(opt$hessian)
        }
        
        modes2.b <- replicate(K, matrix(0, n.tp, ncz), simplify = FALSE)
        Vars2.b <- replicate(K, vector("list", n.tp), simplify = FALSE)
        ff2 <- function(b, y, tt, mm, i, k) { environment(); -JMCR.log.posterior2.b(b, y, ii = i, k)}
        for (i in seq_len(n.tp)) {
          for (k in seq_len(K)) {
            opt2 <- try(optim(modes.b[i, ], ff2, y = y,
                              i = i, k = k, method = "BFGS", hessian = TRUE), 
                        TRUE)
            if (inherits(opt2, "try-error")) {
              gg2 <- function(b, y, tt, mm, i) JM:::cd(b, ff, y = y, i = i)
              opt2 <- optim(modes.b[i, ], ff2, gg2, y = y,
                            i = i, k = k, method = "BFGS", hessian = TRUE)
            }
            modes2.b[[k]][i, ] <- opt2$par
            Vars2.b[[k]][[i]] <- solve(opt2$hessian)
          }
        }
        
        for (i in seq_len(n.tp)) {
          VCdets.last <- Vars.b[[i]]
          b.last <- rep(modes.b[i, ], each = nrow(b)) + sqrt(2) * t(solve(chol(solve(VCdets.last))) %*% t(b))
          S.last <- c(exp(apply(b.last, MARGIN = 1, 
                                FUN = function(x) JMCR.log.posterior.b(x, y, survMats.last, ii = i))) * 
                        1/det(chol(solve(VCdets.last)))) %*% wGH
          for (k in seq_len(K)) {
            VCdets.pred <- Vars2.b[[k]][[i]]
            b.pred <- rep(modes2.b[[k]][i, ], each = nrow(b)) + sqrt(2) * t(solve(chol(solve(VCdets.pred))) %*% t(b))
            S.pred <- c(exp(apply(b.pred, MARGIN = 1,
                                  FUN = function(x) JMCR.log.posterior2.b(x, y, ii = i, k))) *
                          1/det(chol(solve(VCdets.pred)))) %*% wGH
            res.m[[i]][m,k] <- S.pred/S.last
          }
        }
        
        for (i in seq_len(n.tp)) {
          S.last.b <- S.last.punct(modes.b[i, ], survMats.last, ii = i)
          for (k in seq_len(K)) {
            S.pred.b <- S.pred.punct(modes2.b[[k]][i, ], ii = i, k)
            res.m.b[[i]][m,k] <- S.pred.b/S.last.b
          }
        }
        cat(paste("Monte Carlo sample: ", m, "/", M, sep = ""), "\n")
      }
      res.m.old <- res.m
      res.m.b.old <- res.m.b
      res.m <- res.m.b <- replicate(K, matrix(, nrow = n.tp, ncol = 5), simplify = F)
      for(k in 1:K) {
        for (i in seq_len(n.tp)) {
          res.m[[k]][i,] <- c(times.to.pred[[i]],
                              colMeans(res.m.old[[i]], na.rm = T)[k],
                              apply(res.m.old[[i]], 2, median, na.rm = T)[k],
                              apply(res.m.old[[i]], 2, quantile, probs = CI.levels[1], na.rm = T)[k],
                              apply(res.m.old[[i]], 2, quantile, probs = CI.levels[2], na.rm = T)[k])
          res.m.b[[k]][i,] <- c(times.to.pred[[i]],
                                colMeans(res.m.b.old[[i]], na.rm = T)[k],
                                apply(res.m.b.old[[i]], 2, median, na.rm = T)[k],
                                apply(res.m.b.old[[i]], 2, quantile, probs = CI.levels[1], na.rm = T)[k],
                                apply(res.m.b.old[[i]], 2, quantile, probs = CI.levels[2], na.rm = T)[k])
        }
        colnames(res.m[[k]]) <- colnames(res.m.b[[k]]) <- c("times", "Mean", "Median", "Lower", "Upper")
      }
    }
    names(res) <- names(res.b) <- levels(strt) # modif (22/03/2017)
    if (simulate) {
      names(res.m) <- names(res.m.b) <- levels(strt)
      return(list(res.marg = res, res.punct = res.b,
                  res.MC.marg = res.m, res.MC.punct = res.m.b))
    }
    else {
      return(list(res.marg = res, res.punct = res.b))
    }
  }


# inspired by:
# edit(JM:::log.posterior.b)
# edit(JM:::S.b)
# edit(JM:::ModelMats)

JMCR.ModelMats <- function(time, ii, K)  # OK
{
  id.GK <- rep(ii, each = GKk)
  idT.GK <- rep(which(idT == ii), each = GKk)
  wk <- JM:::gaussKronrod(GKk)$wk
  sk <- JM:::gaussKronrod(GKk)$sk
  P <- time/2
  st <- P * (sk + 1)
  P <- rep(P, K)
  st <- rep(st, K)
  data.id2 <- data.id[id.GK, ]
  data.id2[[timeVar]] <- pmax(st - lag, 0)[seq_along(id.GK)]
  out <- list(st = st, wk = rep(wk, length(P)), P = P)
  if (parameterization %in% c("value", "both")) {
    mfX <- model.frame(delete.response(TermsX), data = data.id2)
    mfZ <- model.frame(TermsZ, data = data.id2)
    out$Xs <- model.matrix(formYx, mfX)
    out$Zs <- model.matrix(formYz, mfZ)
    out$Ws.intF.vl <- WintF.vl[idT.GK, , drop = FALSE]
  }
  if (parameterization %in% c("slope", "both")) {
    mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
    mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
    out$Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
    out$Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
    out$Ws.intF.sl <- WintF.sl[idT.GK, , drop = FALSE]
  }
  out
}

JMCR.log.posterior.b <- function(b, y, Mats, ii) # OK
{
  id.i <- id %in% ii
  idT.i <- idT %in% ii
  X.i <- X[id.i, , drop = FALSE]
  Z.i <- Z[id.i, , drop = FALSE]
  mu.y <- as.vector(X.i %*% betas) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
  logNorm <- dnorm(y[id.i], mu.y, sigma, TRUE)
  log.p.yb <- sum(logNorm)
  log.p.b <- JM:::dmvnorm(b, rep(0, ncz), D, TRUE)
  st <- Mats[[ii]]$st
  wk <- Mats[[ii]]$wk
  P <- Mats[[ii]]$P
  Xs <- Mats[[ii]]$Xs
  Zs <- Mats[[ii]]$Zs
  Xs.deriv <- Mats[[ii]]$Xs.deriv
  Zs.deriv <- Mats[[ii]]$Zs.deriv
  Ws.intF.vl <- Mats[[ii]]$Ws.intF.vl
  Ws.intF.sl <- Mats[[ii]]$Ws.intF.sl
  if (parameterization %in% c("value", "both")) 
    Ys <- as.vector(Xs %*% betas + rowSums(Zs * rep(b, each = nrow(Zs))))
  if (parameterization %in% c("slope", "both")) 
    Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
    rowSums(Zs.deriv * rep(b[indRandom], each = nrow(Zs)))
  tt <- switch(parameterization, 
               value = c(Ws.intF.vl %*% alpha) * Ys,
               slope = c(Ws.intF.sl %*% Dalpha) * Ys.deriv,
               both = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
  eta.tw <- if (!is.null(W)) {
    as.vector(W[idT.i, , drop = FALSE] %*% gammas) # modified
  }
  else 0
  log.survival <- {
    kn <- object$control$knots
    strt.i <- strt[idT.i]
    W2s <- mapply(function(k, t) splineDesign(k, t, ord = object$control$ord, outer.ok = TRUE), # outer.ok added
                  kn, split(st, rep(strt.i, each = GKk)), SIMPLIFY = FALSE) # OK
    W2s <- mapply(function(w2s, ind) {
      out <- matrix(0, nrow(Xs) * K, ncol(w2s))
      strt.i.s <- rep(strt.i, each = GKk)
      out[strt.i.s == ind, ] <- w2s
      out
    }, W2s, levels(strt.i), SIMPLIFY = FALSE)
    W2s <- do.call(cbind, W2s)
    W2s.gammas.bs <- c(W2s %*% gammas.bs)
    W2s.gammas.bs[which(W2s.gammas.bs == 0)] <- -Inf
    Vi <- exp(W2s.gammas.bs + tt)
    idT.GK.i <- rep(seq_along(P), each = GKk)
    -sum(exp(eta.tw) * P * tapply(wk * Vi, idT.GK.i, sum))
  }
  log.p.yb + log.survival + log.p.b
}

JMCR.log.posterior2.b <- function (b, y, ii, k)
{  
  id.i <- id %in% ii
  X.i <- X[id.i, , drop = FALSE]
  Z.i <- Z[id.i, , drop = FALSE]
  mu.y <- as.vector(X.i %*% betas) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
  logNorm <- dnorm(y[id.i], mu.y, sigma, TRUE)
  log.p.yb <- sum(logNorm)
  log.p.b <- JM:::dmvnorm(b, rep(0, ncz), D, TRUE)
  
  id.GK <- rep(ii, each = GKk)
  wk <- JM:::gaussKronrod(GKk)$wk
  sk <- JM:::gaussKronrod(GKk)$sk
  low <- last.time[ii]
  upp <- times.to.pred[[ii]]
  P <- (upp - low)/2
  P1 <- (upp + low)/2
  st <- c(outer(P, sk) + P1)
  
  Mats.survival.st <- JMCR.ModelMats.st(time = st, ii = ii, K = K) # vérifié (04/03/16)
  Mats.hazard.st <- JMCR.ModelMats.hazard.st(time = st, ii = ii) # vérifié (04/03/16)
  survival.b.st <- JMCR.S.b.st(b = b, ii = ii, Mats = Mats.survival.st) # vérifié (04/03/16)
  hazard.b.st <- JMCR.hazard.b.st(b = b, ii = ii, Mats = Mats.hazard.st) # vérifié (04/03/16)
  
  integrand.k.st <- survival.b.st * hazard.b.st[(GKk * k - (GKk - 1)) : (GKk * k)]  
  
  log(P * sum(wk * integrand.k.st)) + log.p.yb + log.p.b
}

JMCR.ModelMats.hazard.st <- function (time, ii) # ok (04/03/16)
{
  id.GK <- rep(ii, GKk)
  idT.GK <- rep(which(idT == ii), each = GKk)
  data.id <- data.id[id.GK, ]
  data.id[[timeVar]] <- pmax(time - lag, 0)
  out <- list(time = time)
  if (parameterization %in% c("value", "both")) {
    mfX.id <- model.frame(delete.response(TermsX), data = data.id)
    mfZ.id <- model.frame(TermsZ, data = data.id)
    out$Xtime <- model.matrix(formYx, mfX.id)
    out$Ztime <- model.matrix(formYz, mfZ.id)
    out$W.intF.vl <- WintF.vl[idT.GK, ]
  }
  if (parameterization %in% c("slope", "both")) {
    mfX.deriv.id <- model.frame(TermsX.deriv, data = data.id)
    mfZ.deriv.id <- model.frame(TermsZ.deriv, data = data.id)
    out$Xtime.deriv <- model.matrix(derivForm$fixed, mfX.deriv.id)
    out$Ztime.deriv <- model.matrix(derivForm$random, mfZ.deriv.id)
    out$W.intF.sl <- WintF.sl[idT.GK, ]
  }
  out
}

JMCR.ModelMats.st <- function (time, ii, K) # modif le 30/03/2017
{
  id.GK <- rep(ii, each = GKk^2)
  idT.GK <- rep(which(idT == ii), each = GKk^2)
  wk <- JM:::gaussKronrod(GKk)$wk
  sk <- JM:::gaussKronrod(GKk)$sk
  P <- time/2
  st <- c(t(outer(P, sk + 1))) # modif (30/03/2017)
  P <- rep(P, K)
  st <- rep(st, K)
  data.id2 <- data.id[id.GK, ]
  data.id2[[timeVar]] <- pmax(st - lag, 0)[seq_along(id.GK)]
  out <- list(st = st, wk = wk, P = P)
  if (parameterization %in% c("value", "both")) {
    mfX <- model.frame(delete.response(TermsX), data = data.id2)
    mfZ <- model.frame(TermsZ, data = data.id2)
    out$Xs <- model.matrix(formYx, mfX)
    out$Zs <- model.matrix(formYz, mfZ)
    out$Ws.intF.vl <- WintF.vl[idT.GK, , drop = FALSE]
  }
  if (parameterization %in% c("slope", "both")) {
    mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
    mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
    out$Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
    out$Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
    out$Ws.intF.sl <- WintF.sl[idT.GK, , drop = FALSE]
  }
  out
}

JMCR.S.b.st <- function (b, ii, Mats) # modif le 30/03/2017
{
  idT.i <- idT %in% ii
  st <- Mats$st
  wk <- Mats$wk
  P <- Mats$P
  Xs <- Mats$Xs
  Zs <- Mats$Zs
  Xs.deriv <- Mats$Xs.deriv
  Zs.deriv <- Mats$Zs.deriv
  Ws.intF.vl <- Mats$Ws.intF.vl
  Ws.intF.sl <- Mats$Ws.intF.sl
  if (parameterization %in% c("value", "both")) 
    Ys <- as.vector(Xs %*% betas + rowSums(Zs * rep(b, each = nrow(Zs))))
  if (parameterization %in% c("slope", "both")) 
    Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
    rowSums(Zs.deriv * rep(b[indRandom], each = nrow(Zs)))
  tt <- switch(parameterization,
               value = c(Ws.intF.vl %*% alpha) * Ys,
               slope = c(Ws.intF.sl %*% Dalpha) * Ys.deriv,
               both = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
  eta.tw <- if (!is.null(W)) {
    rep(as.vector(W[idT.i, , drop = FALSE] %*% gammas), each = GKk)
  }
  else 0
  log.survival.i <- {
    kn <- object$control$knots
    strt.i <- strt[idT.i]
    W2s <- mapply(function(k, t) splineDesign(k, t, ord = object$control$ord, outer.ok = TRUE),
                  kn, split(st, rep(strt.i, each = GKk^2)), SIMPLIFY = FALSE) # modif (30/03/2017)
    W2s <- mapply(function(w2s, ind) {
      out <- matrix(0, nrow(Xs) * K, ncol(w2s))
      strt.i.s <- rep(strt.i, each = GKk^2)
      out[strt.i.s == ind, ] <- w2s
      out
    }, W2s, levels(strt.i), SIMPLIFY = FALSE)
    W2s <- do.call(cbind, W2s)
    W2s.gammas.bs <- c(W2s %*% gammas.bs)
    W2s.gammas.bs[which(W2s.gammas.bs == 0)] <- -Inf
    Vi <- exp(W2s.gammas.bs + tt)
    idT.GK.i <- rep(seq_along(P), each = GKk)
    -exp(eta.tw) * P * tapply(wk * Vi, idT.GK.i, sum)
  }
  exp(rowSums(matrix(log.survival.i, nrow = GKk)))
}

JMCR.hazard.b.st <- function (b, ii, Mats) # modif le 30/03/2017
{
  idT.i <- idT %in% ii
  time <- Mats$time
  Xtime <- Mats$Xtime
  Ztime <- Mats$Ztime
  Xtime.deriv <- Mats$Xtime.deriv
  Ztime.deriv <- Mats$Ztime.deriv
  W.intF.vl <- Mats$W.intF.vl
  W.intF.sl <- Mats$W.intF.sl
  if (parameterization %in% c("value", "both")) 
    Ytime <- as.vector(Xtime %*% betas + rowSums(Ztime * rep(b, each = nrow(Ztime))))
  if (parameterization %in% c("slope", "both")) 
    Ytime.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + 
    rowSums(Ztime.deriv * rep(b[indRandom], each = nrow(Ztime)))
  tt <- switch(parameterization,
               value = c(W.intF.vl %*% alpha) * Ytime,
               slope = c(W.intF.sl %*% Dalpha) * Ytime.deriv,
               both = c(W.intF.vl %*% alpha) * Ytime + 
                 c(W.intF.sl %*% Dalpha) * Ytime.deriv)
  eta.tw <- if (!is.null(W)) {
    rep(as.vector(W[idT.i, , drop = FALSE] %*% gammas), each = GKk)
  }
  else 0
  kn <- object$control$knots
  strt.i <- strt[idT.i]
  W2 <- mapply(function(k, t) splineDesign(k, t, ord = object$control$ord, outer.ok = TRUE),
               kn, split(rep(time, each = K), strt.i), SIMPLIFY = FALSE) # modif (30/03/2017)
  W2 <- mapply(function(w2, ind) {
    out <- matrix(0, nrow(Xtime) * K, ncol(w2))
    strt.i <- rep(strt.i, each = GKk)
    out[strt.i == ind, ] <- w2
    out
  }, W2, levels(strt.i), SIMPLIFY = FALSE)
  W2 <- do.call(cbind, W2)
  W2.gammas.bs <- c(W2 %*% gammas.bs)
  W2.gammas.bs[which(W2.gammas.bs == 0)] <- -Inf
  Vi <- exp(W2.gammas.bs + tt)
  exp(eta.tw) * Vi
}

## Pour calculer les probas punctuelles en mode(p(b|...))
S.pred.punct <- function (b, ii, k)
{
  id.GK <- rep(ii, each = GKk)
  wk <- JM:::gaussKronrod(GKk)$wk
  sk <- JM:::gaussKronrod(GKk)$sk
  low <- last.time[ii]
  upp <- times.to.pred[[ii]]
  P <- (upp - low)/2
  P1 <- (upp + low)/2
  st <- c(outer(P, sk) + P1)
  
  Mats.survival.st <- JMCR.ModelMats.st(time = st, ii = ii, K = K)
  Mats.hazard.st <- JMCR.ModelMats.hazard.st(time = st, ii = ii)
  survival.b.st <- JMCR.S.b.st(b = b, ii = ii, Mats = Mats.survival.st)
  hazard.b.st <- JMCR.hazard.b.st(b = b, ii = ii, Mats = Mats.hazard.st)
  
  integrand.k.st <- survival.b.st * hazard.b.st[(GKk * k - (GKk - 1)) : (GKk * k)]  
  
  P * sum(wk * integrand.k.st)
}

S.last.punct <- function (b, Mats, ii)
{
  idT.i <- idT %in% ii
  st <- Mats[[ii]]$st
  wk <- Mats[[ii]]$wk
  P <- Mats[[ii]]$P
  Xs <- Mats[[ii]]$Xs
  Zs <- Mats[[ii]]$Zs
  Xs.deriv <- Mats[[ii]]$Xs.deriv
  Zs.deriv <- Mats[[ii]]$Zs.deriv
  Ws.intF.vl <- Mats[[ii]]$Ws.intF.vl
  Ws.intF.sl <- Mats[[ii]]$Ws.intF.sl
  if (parameterization %in% c("value", "both")) 
    Ys <- as.vector(Xs %*% betas + rowSums(Zs * rep(b, each = nrow(Zs))))
  if (parameterization %in% c("slope", "both")) 
    Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
    rowSums(Zs.deriv * rep(b[indRandom], each = nrow(Zs)))
  tt <- switch(parameterization, 
               value = c(Ws.intF.vl %*% alpha) * Ys,
               slope = c(Ws.intF.sl %*% Dalpha) * Ys.deriv,
               both = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
  eta.tw <- if (!is.null(W)) {
    as.vector(W[idT.i, , drop = FALSE] %*% gammas)
  }
  else 0
  log.survival.b <- {
    kn <- object$control$knots
    strt.i <- strt[idT.i]
    W2s <- mapply(function(k, t) splineDesign(k, t, ord = object$control$ord, outer.ok = TRUE),
                  kn, split(st, rep(strt.i, each = GKk)), SIMPLIFY = FALSE)
    W2s <- mapply(function(w2s, ind) {
      out <- matrix(0, nrow(Xs) * K, ncol(w2s))
      strt.i.s <- rep(strt.i, each = GKk)
      out[strt.i.s == ind, ] <- w2s
      out
    }, W2s, levels(strt.i), SIMPLIFY = FALSE)
    W2s <- do.call(cbind, W2s)
    W2s.gammas.bs <- c(W2s %*% gammas.bs)
    W2s.gammas.bs[which(W2s.gammas.bs == 0)] <- -Inf
    Vi <- exp(W2s.gammas.bs + tt)
    idT.GK.i <- rep(seq_along(P), each = GKk)
    -sum(exp(eta.tw) * P * tapply(wk * Vi, idT.GK.i, sum))
  }
  exp(log.survival.b)
}