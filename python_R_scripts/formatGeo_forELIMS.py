import pandas as pd
import argparse
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--geoPrediction")
parser.add_argument("--metadata")

args = parser.parse_args()
timeStamp = datetime.today().strftime('%Y-%m-%d_%H%M')

allGeoPredictions = pd.read_csv(args.geoPrediction, sep = "\t")
metadataIn = pd.read_excel(args.metadata)
# metadataSheet = metadataIn.active
print(allGeoPredictions)
combine = pd.merge(allGeoPredictions, metadataIn, left_on = "Sample", right_on = "Seq_ID_final", how = "left")

combine = combine[combine['trueCUID'].notna()]

# print(combine.columns)
columnList = ['Sample', 'Region_Prediction_1', 'Region_Prediction_2', 'Region_Prediction_3', 'CSID_y','trueCUID']
selectedColumns = combine[combine.columns.intersection(columnList)]


selectedColumns[['prediction_1', 'probability_1']] = selectedColumns['Region_Prediction_1'].str.split(" ", expand = True)
selectedColumns[['prediction_2', 'probability_2']] = selectedColumns['Region_Prediction_2'].str.split(" ", expand = True)
selectedColumns[['prediction_3', 'probability_3']] = selectedColumns['Region_Prediction_3'].str.split(" ", expand = True)

selectedColumns['probability_1'] = selectedColumns['probability_1'].map(lambda x: x.lstrip('(').rstrip(')'))
selectedColumns['probability_2'] = selectedColumns['probability_2'].map(lambda x: x.lstrip('(').rstrip(')'))
selectedColumns['probability_3'] = selectedColumns['probability_3'].map(lambda x: x.lstrip('(').rstrip(')'))

finalColumnList = ['Sample','CSID_y', 'trueCUID', 'prediction_1', 'probability_1', 'prediction_2', 'probability_2', 'prediction_3', 'probability_3']
finalOut = selectedColumns[selectedColumns.columns.intersection(finalColumnList)]

finalOut = finalOut.rename(columns={'CSID_y': 'CSID', 'trueCUID': 'CUID'})

fileName = timeStamp + '_pvivax_geoPredictions_forELIMS_beforeClean.csv'
finalOut.to_csv(fileName, index =False)