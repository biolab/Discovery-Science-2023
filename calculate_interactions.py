import pandas as pd
import numpy as np
import os

from joblib import Parallel, delayed

from method.utils import Logp1Scaler, StandardScaler, load_TCGA_survival_data, load_METABRIC_survival_data
from method.rmst_diff import rmst_diff


l1000 = set()
with open('data/L1000.txt', 'r') as file:
    for line in file:
        l1000.add(line.rstrip())


def run(dataset: str):
    # load data
    if dataset == 'METABRIC':
        data = load_METABRIC_survival_data()
        data = data.dropna()
        data = data.loc[:,~data.columns.duplicated()].copy()
    else:
        data = load_TCGA_survival_data(dataset)

    time_col = data.pop('time')
    event_col = data.pop('event')

    # Remove low expressed genes with 75th percentile lower than 10
    if not dataset == 'METABRIC':
        to_keep = np.percentile(data, 75, axis=0) >= 10
        data = data.loc[:, to_keep]


    # transform data
    data = Logp1Scaler().fit_transform(data)
    data = StandardScaler().fit_transform(data)
    data = pd.merge(data, pd.concat([time_col, event_col], axis=1), left_index=True, right_index=True)

    # compute interactions
    y = data[['time', 'event']].copy()
    X = data[data.columns.difference(['time', 'event'])].copy()
    X = X[X.columns.intersection(list(l1000)[:10])]
    print(f'DONE: {dataset}, data size {X.shape}')

    interaction = rmst_diff(X, y)
    
    df = pd.DataFrame(interaction, columns=['interaction', 'rmst_f1', 'rmst_f2', 'rmst_f1+f2', 'rmst_f1-f2', 'rmst_f1*f2','limit_is_lower'])
    df['rmst_diff_f1+f2'] = (df[['rmst_f1', 'rmst_f2']].max(axis=1) -  df['rmst_f1+f2']).abs()
    df['rmst_diff_f1-f2'] = (df[['rmst_f1', 'rmst_f2']].max(axis=1) -  df['rmst_f1-f2']).abs()
    df['rmst_diff_f1*f2'] = (df[['rmst_f1', 'rmst_f2']].max(axis=1) -  df['rmst_f1*f2']).abs()

    df.to_csv(f'{dataset}.csv', index=False)

    return f'DONE: {dataset}, data size {X.shape}'


DATASETS = [
    "BLCA",
    "BRCA",
    "CESC",
    "COAD",
    "GBM",
    "HNSC",
    "KIRC",
    "KIRP",
    "LAML",
    "LGG",
    "LIHC",
    "LUAD",
    "LUSC",
    "OV",
    "PRAD",
    "READ",
    "SKCM",
    "STAD",
    "THCA",
    "UCEC",
    "METABRIC"
]

Parallel(n_jobs=1)(delayed(run)(data) for data in DATASETS)

