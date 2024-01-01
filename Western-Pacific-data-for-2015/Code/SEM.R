#Create artboard
dev.new(title="Partial RDA", width=10,height=8,
        noRStudioGD = TRUE)
install.packages("lavaan")

install.packages("rlang")
install.packages("haven")
install.packages("Hmisc")
install.packages("semPolt")
install.packages("devtools")
library(devtools)
install_github("jslefche/piecewiseSEM@devel",build_vignette=F)
library(lavaan)
library(haven)
library(Hmisc)

install_github("SachaEpskamp/semPlot")
library(semPlot)
library(sem)
data <- read.csv(file.choose(),header = T, row.names = 1)#OTU
df <- scale(data,center = T, scale = TRUE)

model <- '
Temperature ~ Latitude +SiO3
Prokaryotic_diversity ~ NH4 + NO2 + NO3+PO4+Temperature +Salinity +SSTA +SiO3
Eukaryotic_diversity~ PO4 + Temperature + NO3 + NH4 + NO2 + Salinity +SSTA+SiO3'#The model needs to be adjusted by multiple tests

df <- as.data.frame(df)
fit_model <- lavaan::sem(model,data = df, se = 'bootstrap')

lavaan::summary(fit_model, standardize = TRUE, rsq = TRUE, modindices = TRUE)

lavaan::fitmeasures(fit_model, c("chisq", "df", "pvalue", "cfi", "gfi", "rmsea", 
                         "AIC"))#The smaller the chisq, the better. The pvalue is required > 0.05, the closer cfi\gfi is to 1, the better rmsea is to 0, and the smaller AIC is, the better
semPaths(fit_model,"std",edge.label.cex = 0.8,fade = FALSE, layout = "spring",
         optimizeLatRes = FALSE, residuals = FALSE)
