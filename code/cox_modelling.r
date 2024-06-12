library(survival)
library(survcomp)
library(readxl)
library(tools)
library(yaml)

##### VARIABLES #####
# Main directory - directory containing labelled feature csvs
# Negative controls - list of negative controls to run CPH on
# Model features - list of features to use for CPH model
# Model weights - matrix of weights to apply for trained CPH model


###### STEPS #######
# Pretrained CPH Model
# 1. Load csv file for testing patient set
# 2. Drop the patient ID and survival labels from the dataframe
# 3. Confirm that features in model features list match the radiomics
# 3. Convert the features dataframe to a data.matrix
# 4. Multiple the weight matrix by the features
# 5. Set up time and event label
# 6. Calculate concordance index 


# From scratch CPH Model training
# 1. Load csv file of training patient set 
# 2. Drop the patient ID and survival labels from the dataframe
# 3. Confirm that features in model features list match the radiomics
# 4. Fit cph model with the select features with the training set
# 5. Get weights for model
# 6. Return the trained CPH model

loadDataFile <- function(dataFilePath) { #nolint
    dataFileType = file_ext(dataFilePath)
    if (length(dataFileType) == 0) {
        # No file passed, return 0
        return(NULL)
    }
    if (dataFileType == "csv") {
        loadedDataframe <- read.csv(dataFilePath, header = TRUE, sep = ",", check.names = FALSE)
    } else if (dataFileType == "xlsx") {
        loadedDataframe <- read_excel(dataFilePath)
    } else {
        stop("Radiogenomic data file must be a .csv or .xlsx.")
    }

    return(loadedDataframe)
}


trainCoxModel <- function(csvTrainingFeatures,
                          patientIDLabel,
                          survTimeLabel,
                          survEventLabel,
                          modelFeatureList,
                          modelFeatureWeights){ #nolint

    trainRadiomicsData <- loadDataFile(csvTrainingFeatures)

    timeData <- trainRadiomicsData[, survTimeLabel]
    eventData <- trainRadiomicsData[, survEventLabel]

    model.fit <- coxph(Surv(timeData, eventData) ~ modelFeatureWeights,
                       x = TRUE,
                       y = TRUE,
                       method = "breslow",
                       data = trainRadiomicsData)

    # Get weights from model and return those two

    return(model.fit)
}

testCoxModel <- function(csvTestingFeatures,
                         patientIDLabel,
                         survTimeLabel,
                         survEventLabel,
                         modelFeatureList,
                         modelFeatureWeights){ #nolint

    testRadiomicsData <- loadDataFile(csvTestingFeatures)

    # INSERT CHECK FOR RADIOMIC COLUMNS AND MODELFEATURELIST

    testRadiomicFeatures <- testRadiomicsData[, !(colnames(testRadiomicsData) %in% c(patientIDLabel, survTimeLabel, survEventLabel))]

    # testRadiomicFeatures <- subset(testRadiomicsData, select = -c(patientIDLabel, survTimeLabel, survEventLabel))

    testRadMatrix <- data.matrix(testRadiomicFeatures)

    radiomicHazards <- testRadMatrix %*% modelFeatureWeights

    timeLabel <- testRadiomicsData[, survTimeLabel]
    eventLabel <- testRadiomicsData[, survEventLabel]

    performanceResults <- concordance.index(x = radiomicHazards,
                                            surv.time = timeLabel,
                                            surv.event = eventLabel,
                                            method = "noether",
                                            alpha = 0.5,
                                            alternative = "two.sided")

    return(performanceResults)
}


##### VARIABLES #####
# Main directory - directory containing labelled feature csvs
# Negative controls - list of negative controls to run CPH on
# Model features - list of features to use for CPH model
# Model weights - matrix of weights to apply for trained CPH model
datasetConfig <- read_yaml("config/HEAD-NECK-RADIOMICS-HN1_config.yaml")

featureDirPath <- datasetConfig$output_dir_path
datasetName <- datasetConfig$dataset_name

signatureName <- "aerts"
signatureConfig <- read_yaml(paste("signatures/", signatureName, ".yaml", sep = ""))
sigFeatures <- names(signatureConfig$signature)
sigWeights <- matrix(unlist(signatureConfig$signature))

testDirPath <- paste(featureDirPath, "/", signatureName, "_signature/", signatureName, "_radiomic_features_", datasetName, ".csv", sep = "")

# trainDirPath <- paste(featureDirPath, "/training/", signatureName, "_signature/training_", signatureName, "_radiomics_features_", datasetName, ".csv", sep = "")
# testDirPath <- paste(featureDirPath, "/test/", signatureName, "_signature/test_", signatureName, "_radiomics_features_", datasetName, ".csv", sep = "")

if (datasetConfig$need_bool_event == 1) {
    eventLabel <- paste(datasetConfig$outcome_status$event_label, "_bool", sep = "")
} else {
    eventLabel <- datasetConfig$outcome_status$event_label
}


testResults <- testCoxModel(csvTestingFeatures = testDirPath,
                            patientIDLabel = "patientID",
                            survTimeLabel = datasetConfig$outcome_status$time_label,
                            survEventLabel = eventLabel,
                            modelFeatureList = sigFeatures,
                            modelFeatureWeights = sigWeights)