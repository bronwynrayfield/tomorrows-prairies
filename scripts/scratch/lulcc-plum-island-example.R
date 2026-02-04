library(lulcc)

## Plum Island Ecosystems

## HACK: temporarily override compareRaster to ignore CRS checks
.compareRaster_orig <- raster::compareRaster
unlockBinding("compareRaster", asNamespace("raster"))
assign("compareRaster",
       function(..., crs = FALSE, stopiffalse = TRUE) .compareRaster_orig(..., crs = crs, stopiffalse = stopiffalse),
       envir = asNamespace("raster"))
lockBinding("compareRaster", asNamespace("raster"))


## load data
data(pie)

## observed maps
obs <- ObsLulcRasterStack(x=pie,
                          pattern="lu",
                          categories=c(1,2,3),
                          labels=c("Forest","Built","Other"),
                          t=c(0,6,14))
obs
plot(obs)  #not working but ok
crossTabulate(obs, times=c(0,14))

## explanatory variables
ef <- ExpVarRasterList(x=pie, pattern="ef")
ef
part <- partition(x=obs[[1]], size=0.1, spatial=TRUE)
train.data <- getPredictiveModelInputData(obs=obs, ef=ef, cells=part[["train"]])
forms <- list(Built ~ ef_001+ef_002+ef_003,
              Forest ~ ef_001+ef_002,
              Other ~ ef_001+ef_002)
glm.models <- glmModels(formula=forms, family=binomial, data=train.data, obs=obs)
rpart.models <- rpartModels(formula=forms, data=train.data, obs=obs)
rf.models <- randomForestModels(formula=forms, data=train.data, obs=obs)

## test ability of models to predict allocation of forest, built and other
## land uses in testing partition
test.data <- getPredictiveModelInputData(obs=obs, ef=ef, cells=part[["test"]])
glm.pred <- PredictionList(models=glm.models, newdata=test.data)
glm.perf <- PerformanceList(pred=glm.pred, measure="rch")
rpart.pred <- PredictionList(models=rpart.models, newdata=test.data)
rpart.perf <- PerformanceList(pred=rpart.pred, measure="rch")
rf.pred <- PredictionList(models=rf.models, newdata=test.data)
rf.perf <- PerformanceList(pred=rf.pred, measure="rch")
plot(list(glm=glm.perf, rpart=rpart.perf, rf=rf.perf))

## test ability of models to predict location of urban gain 1985 to 1991
part <- rasterToPoints(obs[[1]], fun=function(x) x != 2, spatial=TRUE)
test.data <- getPredictiveModelInputData(obs=obs, ef=ef, cells=part, t=6)
glm.pred <- PredictionList(models=glm.models[[2]], newdata=test.data)
glm.perf <- PerformanceList(pred=glm.pred, measure="rch")
plot(list(glm=glm.perf))

## obtain demand scenario
dmd <- approxExtrapDemand(obs=obs, tout=0:14)
matplot(dmd, type="l", ylab="Demand (no. of cells)", xlab="Time point",
        lty=1, col=c("Green","Red","Blue"))
legend("topleft", legend=obs@labels, col=c("Green","Red","Blue"), lty=1)

## get neighbourhood values
w <- matrix(data=1, nrow=3, ncol=3)
nb <- NeighbRasterStack(x=obs[[1]], weights=w, categories=2)

## create CLUE-S model object
clues.rules <- matrix(data=1, nrow=3, ncol=3, byrow=TRUE)
clues.parms <- list(jitter.f=0.0002,
                    scale.f=0.000001,
                    max.iter=1000,
                    max.diff=50,
                    ave.diff=50)
clues.model <- CluesModel(obs=obs,
                          ef=ef,
                          models=glm.models,
                          time=0:14,
                          demand=dmd,
                          elas=c(0.2,0.2,0.2),
                          rules=clues.rules,
                          params=clues.parms)

## Create Ordered model
ordered.model <- OrderedModel(obs=obs,
                              ef=ef,
                              models=glm.models,
                              time=0:14,
                              demand=dmd,
                              order=c(2,1,3))

## perform allocation
clues.model <- allocate(clues.model)
ordered.model <- allocate(ordered.model, stochastic=TRUE)

## pattern validation
## CLUE-S
clues.tabs <- ThreeMapComparison(x=clues.model,
                                 factors=2^(1:8),
                                 timestep=14)
plot(clues.tabs)
plot(clues.tabs, category=1, factors=2^(1:8)[c(1,3,5,7)])

## Ordered
ordered.tabs <- ThreeMapComparison(x=ordered.model,
                                   factors=2^(1:8),
                                   timestep=14)
plot(ordered.tabs)
plot(ordered.tabs, category=1, factors=2^(1:8)[c(1,3,5,7)])

## calculate agreement budget and plot
## CLUE-S
clues.agr <- AgreementBudget(x=clues.tabs)
plot(clues.agr, from=1, to=2)
## Ordered
ordered.agr <- AgreementBudget(x=ordered.tabs)
plot(ordered.agr, from=1, to=2)

## calculate Figure of Merit and plot
## CLUE-S
clues.fom <- FigureOfMerit(x=clues.tabs)
p1 <- plot(clues.fom, from=1, to=2)
## Ordered
ordered.fom <- FigureOfMerit(x=ordered.tabs)
p2 <- plot(ordered.fom, from=1, to=2)


## Plot observed and predicted maps for t6
## land-use codes
lu_vals   <- c(1, 2, 3)
lu_labels <- c("Forest", "Built", "Other")

## choose colours (example)
lu_cols <- c("darkgreen", "red", "grey70")

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(obs[[2]],
     col = lu_cols,
     breaks = c(0.5, 1.5, 2.5, 3.5),
     legend = FALSE,
     main = "Observed t6 (1991)")

plot(pred_t6,
     col = lu_cols,
     breaks = c(0.5, 1.5, 2.5, 3.5),
     legend = FALSE,
     main = "Predicted t6")

par(mfrow = c(1, 1))
