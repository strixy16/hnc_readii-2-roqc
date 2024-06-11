import pandas as pd 
import numpy as np 

from pandas import Series, DataFrame
from collections.abc import Sequence, Mapping

def subsetDataframe(dataframe:DataFrame,
                    columnName:str,
                    includeValues:Series | DataFrame | Sequence | Mapping = None,
                    excludeValues:Series | DataFrame | Sequence | Mapping = None):
    """Get rows of a pandas DataFrame based on values in a specific column.

    Parameters
    ----------
    dataframe : pd.DataFrame
        Dataframe to get a subset of rows from
    
    columnName : str
        Column to use to subset the dataframe based on includeValues and excludeValues
    
    includeValues: iterable, Series, DataFrame, list
        Values marking rows to include from the columName column in the subset dataframe
    
    excludeValues: iterable, Series, DataFrame, list
        Values marking rows to exclude from the columnName column in the subset dataframe


    Returns
    -------
    pd.DataFrame
        Dataframe with subset of rows based on the provided values for the columnName column.
    
    """
    if columnName in ["Index", "index"]:
        subsetByIndex = True
    else:
        subsetByIndex = False
    
    if (includeValues is None) and (excludeValues is None):
        raise ValueError("Must provide one of includeValues or excludeValues.")

    if (not subsetByIndex) and (columnName not in dataframe.columns):
        raise ValueError((columnName + " does not exist in this dataframe."))
    
    if includeValues is not None:
        if subsetByIndex:
            dataframe = dataframe[dataframe.index.isin(includeValues)]
        else:
            # Get rows of dataframe that have values in columnName column that are in the includeValues
            dataframe = dataframe[dataframe[columnName].isin(includeValues)]

    if excludeValues is not None:
        if subsetByIndex:
            dataframe = dataframe[~dataframe.index.isin(excludeValues)]
        else:
            # Get rows of dataframe that don't are not in excludeValues in the columnName column
            dataframe = dataframe[~dataframe[columnName].isin(excludeValues)]

    return dataframe


def outcomeLabelSetup(dfClinical:DataFrame,
                      statusLabel:str = "Status",
                      statusValues:list = [],
                      followupLabel:str = "Length FU"):
    """Function to setup status and follow up labels for predictive models. 

    Parameters
    ----------
    dfClinical : DataFrame
        Dataframe containing clinical data. Must include a column with statusLabel and followupLabel.
    
    statusLabel : str
        Name of outcome status label in dfClinical (e.g. Status, Event)
    
    statusValues : list
        List of status values for boolean conversion, in the order they should be numbered. 
        If empty list passed, will get list of unique values from statusLabel column and they will
        be numbered in alphabetical order.
    
    followupLabel : str
        Name of follow-up label column in dfClinical. Should point to a numeric column of something like
        days to last follow-up or death. 
    
    Returns
    -------
    outcomeLabels : DataFrame
        Dataframe with a boolean version of the status label and the follow up label
    """

    # If no status values passed, get all of them from the status column
    if statusValues == []:
        statusValues = list(dfClinical[statusLabel].unique())

    # Convert status to boolean
    statusNumeric = list(range(0, len(statusValues)))
    statusBoolMap = dict(zip(statusValues, statusNumeric))

    statusBool = dfClinical[statusLabel].map(statusBoolMap)
    statusBoolLabel = statusLabel + "_bool"

    outcomeLabels = dfClinical[followupLabel]
    outcomeLabels = outcomeLabels.to_frame()
    outcomeLabels.insert(0, statusBoolLabel, statusBool)

    return outcomeLabels


def splitDataSetup(dfClinical:DataFrame,
                   dfFeatures:DataFrame,
                   splitVariables:dict,
                   ):
    """Function to split clinical and feature data tables into subgroups (e.g. train and test) based on clinical variable.
       

    Parameters
    ----------
    dfClinical : DataFrame
        Dataframe containing clinical data. Must include a column with labels from splitVariables.
        The index of this dataframe must correspond with dfFeatures' index.

    dfFeatures : DataFrame
        Dataframe containing feature data. This table will be split based on the patient index of the clinical table.
        The index of this dataframe must correspond with dfClinical's index.

    splitVariables : dict
        Dictionary where keys are labels of columns in dfClinical to base the split upon and values is a list of values to select
        in each data split.

    Returns
    -------
    splitClinicalDataframes : dict
        Dictionary of split up clinical dataframes. Keys indicate the value used to select the subset, value is the Dataframe
    
    splitFeatureDataframes : dict
        Dictionary of split up feature dataframes. Keys indicate the value used to select the clinical subset, value is the Dataframe

    Examples
    --------
    >>> splitDataSetup(clinicalData, featureData, {'data_split': ['train', 'test']})
    """

    # Initialize dictionaries to contain the split up dataframes
    splitClinicalDataframes = {}
    splitFeatureDataframes = {}

    # Get the clinical variable to split by and the list of values to select in each subgroup
    for variable, values in splitVariables.items():
        print("Getting split for ", variable)

        for value in values:
            # Get the clinical subset for rows with value in the variable column
            splitClinical = subsetDataframe(dataframe = dfClinical,
                                           columnName = variable,
                                           includeValues = [value])
            # Save this subset dataframe into the output dictionary
            splitClinicalDataframes[value] = splitClinical

            # Get the feature subset based on the indices from splitClinical
            splitFeature = subsetDataframe(dataframe = dfFeatures,
                                           columnName = "index",
                                           includeValues = splitClinical.index)
            # Save this subset dataframe into the output dictionary
            splitFeatureDataframes[value] = splitFeature
    
    return splitClinicalDataframes, splitFeatureDataframes



def filterDataSetup(dfClinical:DataFrame, 
                    dfFeatures:DataFrame,
                    subsetIncludeVariables: dict = None,
                    subsetExcludeVariables: dict = None):
    """Filter clinical and feature datasets by patients existing in both and any specified clinical inclusion/exclusion variables

    Parameters
    ----------
    dfClinical : DataFrame
        Dataframe containing clinical data. Must include a column with labels from subsetInclude/ExcludeVariables.
        The index of this dataframe must correspond with dfFeatures' index.

    dfFeatures : DataFrame
        Dataframe containing feature data. Must include a column with labels from subsetInclude/ExcludeVariables.
        The index of this dataframe must correspond with dfClinical's index.

    subsetIncludeVariables : dict
        Dictionary where keys are names of columns to filter for inclusion from the clinical dataframe as keys and 
        the values is a list of values from that column to select included rows.
    
    subsetExcludeVariables : dict
        Dictionary where keys are names of columns to filter for exclusion from the clinical dataframe as keys and 
        the values is a list of values from that column to select excluded rows.
        Example:
        subsetExcludeVariables = {'RADCURE-challenge': [0],
                                'Ds Site': ['Sarcoma', 'Unknown', 'Paraganglioma', 'Salivary Glands', 'Other', 'benign tumor', 'Orbit', 'Lacrimal gland'] }
    
    Returns
    -------
    filteredClinicalData : DataFrame
        Dataframe of clinical data, with the patient identifier as the index, filtered by the patients with features 
        available and any subset variables provided.
    filteredFeatureData : DataFrame
        Dataframe of features data, with the patient identifier as the index, filtered by any subset variables provided.
    """

    # Get list of patients with features
    patientListWithFeats = dfFeatures.index

    # Get clinical data with features
    filteredClinicalData = dfClinical[dfClinical.index.isin(patientListWithFeats)]
    # subsetDataframe(dfClinical, columnName=dfClinical.index, includeValues=patientListWithFeats)

    # Filter the clinical data by variables with values specified to include
    if subsetIncludeVariables is not None:
        for variable, values in subsetIncludeVariables.items():
            filteredClinicalData = subsetDataframe(dataframe = filteredClinicalData,
                                                   columnName = variable,
                                                   includeValues = values)

    # Filter the clinical data by variables with values specified to exclude
    if subsetExcludeVariables is not None:
        for variable, values in subsetExcludeVariables.items():
            filteredClinicalData = subsetDataframe(dataframe = filteredClinicalData,
                                                   columnName = variable,
                                                   excludeValues = values)

    filteredPatientList = filteredClinicalData.index

    # Update the features data based on the subset variables
    filteredFeatureData = dfFeatures[dfFeatures.index.isin(filteredPatientList)]

    return filteredClinicalData, filteredFeatureData


def getPatientIdentifierLabel(dataframeToSearch:DataFrame):
    """Function to find a column in a dataframe that contains some form of patient ID or case ID (case-insensitive). 
       If multiple found, will return the first match.

    Parameters
    ----------
    dataframeToSearch : DataFrame
        Dataframe to look for a patient ID column in.

    Returns
    -------
    str
        Label for patient identifier column from the dataframe.
    """

    # regex to get patient identifier column name in the dataframes
    # catches case-insensitive variations of patient_id, patid, pat id, case_id, case id, caseid
    regexSearchTerm = r"(?i)pat(ient)?(_| )?id|(?i)case(_| )?id|id"

    patIdentifier = dataframeToSearch.filter(regex=regexSearchTerm).columns.to_list()

    if len(patIdentifier) > 1:
        print("Multiple patient identifier labels found. Using the first one.")
    
    elif len(patIdentifier) == 0:
        raise ValueError("Dataframe doesn't have a recognizeable patient ID column. Must contain patient or case ID.")

    return patIdentifier[0]


def dropPyradiomicsDiagnostics(dfPyradiomicsFeatures:DataFrame):
    """ Function to get out just the features from a Pyradiomics output that includes diagnostics columns before the features.
        Assumes the features start after the last diagnostics column.
    Parameters
    ----------
    dfPyradiomicsFeatures : DataFrame
        Dataframe of Pyradiomics features with diagnostics and other columns before the features
    
    Returns
    -------
    featsOnlyRadiomics : DataFrame
        Dataframe with just the radiomic features
    
    """
    # Find all the columns that begin with diagnostics
    diagnosticRadiomics = dfPyradiomicsFeatures.filter(regex=r"diagnostics_*")

    if not diagnosticRadiomics.empty:
        # Get the last diagnostics column index - the features begin in the next column
        lastDiagnosticIdx = dfPyradiomicsFeatures.columns.get_loc(diagnosticRadiomics.columns[-1])
        # Drop all the columns before the features start
        featsOnlyRadiomics = dfPyradiomicsFeatures.iloc[:, lastDiagnosticIdx+1:]

    else:
        originalRadiomics = dfPyradiomicsFeatures.filter(regex=r'^original_*')
        if not originalRadiomics.empty:
            # Get the first original feature column index - the features begin in this column
            firstOriginalIdx = dfPyradiomicsFeatures.columns.get_loc(originalRadiomics.columns[0])
            # Drop all the columns before the features start
            featsOnlyRadiomics = dfPyradiomicsFeatures.iloc[:, firstOriginalIdx:]
        else:
            raise ValueError("PyRadiomics file doesn't contain any diagnostics or original feature columns, so can't find beginning of features.")

    return featsOnlyRadiomics


def modelDataSetup(clinicalDataPath,
                   featureDataPath,
                   clinicalPatID = None,
                   featurePatID = None,
                   subsetIncludeVariables = None,
                   subsetExcludeVariables = None,):
    """
    Parameters
    ----------
    clinicalDataPath : str
        Path to the clinical data file, expect xlsx file
    featureDataPath : str
        Path to the features data file, expect csv file
    clinicalPatID : str
        Name of the column in the clinical data that has patient identifiers. If not passed, will look for 
        a variable that resembles patient_id, pat_id, or case id (case-insensitive)
    subsetIncludeVariables : dict
        Dictionary where keys are names of columns to filter for inclusion from the clinical dataframe as keys and 
        the values is a list of values from that column to select included rows.
    subsetExcludeVariables : dict
        Dictionary where keys are names of columns to filter for exclusion from the clinical dataframe as keys and 
        the values is a list of values from that column to select excluded rows.
        Example:
        subsetExcludeVariables = {'RADCURE-challenge': [0],
                                'Ds Site': ['Sarcoma', 'Unknown', 'Paraganglioma', 'Salivary Glands', 'Other', 'benign tumor', 'Orbit', 'Lacrimal gland'] }
    """

    # Load in the clinical data into a DataFrame
    completeClinicalData = pd.read_excel(clinicalDataPath)
    # Load in the features data into a DataFrame
    featureData = pd.read_csv(featureDataPath)

    if clinicalPatID is None:
        clinicalPatID = getPatientIdentifierLabel(completeClinicalData)

    if featurePatID is None: 
        featurePatID = getPatientIdentifierLabel(featureData)

    # Set patient IDs as index in both dataframes
    completeClinicalData = completeClinicalData.set_index(clinicalPatID)
    featureData = featureData.set_index(featurePatID) 

    # Filter patients by those with features, and then any variables marked with exclusion conditions
    filteredClinical, filteredRadiomics = filterDataSetup(completeClinicalData, featureData, subsetIncludeVariables, subsetExcludeVariables)

    
    

    return None


if __name__ == "__main__":

    # Step 1: Load clinical and radiomic data into pandas dataframes
    # Step 2: get intersection of patients from clinical and radiomics
    # Step 2a (optional): Get subset of patients to work with - e.g. RADCURE challenge, disease site
    #                  note - can be an input arg, column in clinical data

    # Step 3: separate clinical into train and test based on column or with random seed
    # Step 4: separate radiomics into train and test based on labels in clinical data
    # Step 5: set patient_id as index in both dataframes
    # Step 5: save out train and test 

    clinicalDataPath = "/Users/katyscott/Documents/RADCURE/RADCURE_TCIA_Clinical June 13 2023.xlsx"
    radiomicDataPath = "/Users/katyscott/Documents/RADCURE/RADCURE_updated_data/uhn_radcure_plus_aerts/snakemake_RADCURE_radiomic_features.csv"
    subsetExcludeVariables = {'RADCURE-challenge': [0],
                              'Ds Site': ['Sarcoma', 'Unknown', 'Paraganglioma', 'Salivary Glands', 'Other', 'benign tumor', 'Orbit', 'Lacrimal gland'] }

    filteredClinical, filteredRadiomics = modelDataSetup(clinicalDataPath, radiomicDataPath, subsetExcludeVariables=subsetExcludeVariables)

    
