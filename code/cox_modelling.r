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
                          survTimeLabel,
                          survEventLabel,
                          modelFeatureList){ #nolint

    trainRadiomicsData <- loadDataFile(csvTrainingFeatures)

    timeData <- trainRadiomicsData[, survTimeLabel]
    eventData <- trainRadiomicsData[, survEventLabel]

    trainModelFeatures <- trainRadiomicsData[, modelFeatureList]

    model.fit <- coxph(Surv(timeData, eventData) ~ .,
                       x = TRUE,
                       y = TRUE,
                       method = "breslow",
                       data = trainModelFeatures)

    # Get weights from model and return those two
    return(model.fit$coefficients)
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


saveSignature <- function(modelFeatureNames, modelFeatureWeights, outputDir, signatureName){ #nolint 
    # Save out model weights for CPH model

    # Setupt output file name
    outputFile <- file(paste(outputDir, "/", signatureName, ".yaml", sep = ""), "w")

    # Convert model weights to list and set up names
    modelFeatureList <- as.list(modelFeatureWeights)
    names(modelFeatureList) <- modelFeatureNames

    # Name the signature so ouptut is correctly formatted
    finalSignature <- list(signature = modelFeatureList)

    # Write out the signature
    write_yaml(finalSignature, outputFile)
    close(outputFile)
}

signatureYAMLSetup <- function(signatureName) { #nolint 
    signatureConfig <- read_yaml(paste("signatures/", signatureName, ".yaml", sep = ""))
    # Names of the features in the signature
    sigFeatureNames <- names(signatureConfig$signature)
    # Weights for the features in the signature
    sigWeights <- matrix(unlist(signatureConfig$signature))

    signature <- list(sigFeatureNames, sigWeights)
    names(signature) <- c("names", "weights")

    return(signature)
}


applySignature <- function(configFilePath, signatureName) { #nolint

    datasetConfig <- read_yaml(configFilePath)

    # Path to directory containing radiomics features
    featureDirPath <- datasetConfig$output_dir_path
    # Name of the dataset to run CPH on
    datasetName <- datasetConfig$dataset_name

    # Signature setup - get the signature features and weights
    signature <- signatureYAMLSetup(signatureName)
    sigFeatureNames <- signature$names
    sigWeights <- signature$weights

    # Determine whether to load data as train test split or not
    if (datasetConfig$train_test_split == 1) {
        trainDirPath <- paste(featureDirPath, "/training/", signatureName, "_signature/training_", signatureName, "_radiomic_features_", datasetName, ".csv", sep = "")
        testDirPath <- paste(featureDirPath, "/test/", signatureName, "_signature/test_", signatureName, "_radiomic_features_", datasetName, ".csv", sep = "")
    } else {
        testDirPath <- paste(featureDirPath, "/", signatureName, "_signature/", signatureName, "_radiomic_features_", datasetName, ".csv", sep = "")
    }

    # Outcome status setup
    # Time to event label
    timeLabel <- datasetConfig$outcome_status$time_label
    # If the event label is boolean, convert it to a string
    if (datasetConfig$need_bool_event == 1) {
        eventLabel <- paste(datasetConfig$outcome_status$event_label, "_bool", sep = "")
    } else { # Otherwise, use the event label as is
        eventLabel <- datasetConfig$outcome_status$event_label
    }

    # Run the CPH model with signature features and weights on test set
    testResults <- testCoxModel(csvTestingFeatures = testDirPath,
                                patientIDLabel = "patientID",
                                survTimeLabel = datasetConfig$outcome_status$time_label,
                                survEventLabel = eventLabel,
                                modelFeatureList = sigFeatures,
                                modelFeatureWeights = sigWeights)

    print("Original radiomics features")
    print(paste("c-index:", testResults$c.index, sep = " "))
    # print(paste("se:", testResults$se, sep = " "))
    print(paste("lower:", testResults$lower, sep = " "))
    print(paste("upper:", testResults$upper, sep = " "))
    print(paste("p-value:", testResults$p.value, sep = " "))
    print("------------------------------------------------------------")

    # Repeat the above for each negative control
    for (negControl in datasetConfig$negative_control_names) {
        # Add the negative control name to the dataset name
        datasetName <- paste(negControl, datasetConfig$dataset_name, sep = "_")

        if (datasetConfig$train_test_split == 1) {
            trainDirPath <- paste(featureDirPath, "/training/", signatureName, "_signature/training_", signatureName, "_radiomic_features_", datasetName, ".csv", sep = "")
            testDirPath <- paste(featureDirPath, "/test/", signatureName, "_signature/test_", signatureName, "_radiomic_features_", datasetName, ".csv", sep = "")
        } else {
            testDirPath <- paste(featureDirPath, "/", signatureName, "_signature/", signatureName, "_radiomic_features_", datasetName, ".csv", sep = "")
        }

        if (datasetConfig$need_bool_event == 1) {
            eventLabel <- paste(datasetConfig$outcome_status$event_label, "_bool", sep = "")
        } else {
            eventLabel <- datasetConfig$outcome_status$event_label
        }

        ncTestResults <- testCoxModel(csvTestingFeatures = testDirPath,
                                      patientIDLabel = "patientID",
                                      survTimeLabel = datasetConfig$outcome_status$time_label,
                                      survEventLabel = eventLabel,
                                      modelFeatureList = sigFeatures,
                                      modelFeatureWeights = sigWeights)

        print("")
        print(negControl)
        print(paste("c-index:", ncTestResults$c.index, sep = " "))
        # print(paste("se:", ncTestResults$se, sep = " "))
        print(paste("lower:", ncTestResults$lower, sep = " "))
        print(paste("upper:", ncTestResults$upper, sep = " "))
        print(paste("p-value:", ncTestResults$p.value, sep = " "))
        print("")
        print("------------------------------------------------------------")
    }
}



##### VARIABLES #####
# Main directory - directory containing labelled feature csvs
# Negative controls - list of negative controls to run CPH on
# Model features - list of features to use for CPH model
# Model weights - matrix of weights to apply for trained CPH model
datasetConfigPath <- "config/RADCURE_config.yaml"
signatureName <- "RADCURE_radiomics"

applySignature(datasetConfigPath, signatureName)
