import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler as SKStandardScaler
from sklearn.preprocessing import MinMaxScaler as SKMinMaxScaler


class StandardScaler(SKStandardScaler):
    """
    Wraper for sklearn's StandardScaler so that it returns
    pandas DataFrame when transforming.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def fit_transform(self, X, y=None):
        if type(X) == pd.DataFrame:
            return pd.DataFrame(
                super().fit_transform(X), index=X.index, columns=X.columns
            )
        return super().fit_transform(X)

    def transform(self, X, y=None):
        if type(X) == pd.DataFrame:
            return pd.DataFrame(super().transform(X), index=X.index, columns=X.columns)
        return super().transform(X)


class Logp1Scaler(BaseEstimator, TransformerMixin):
    """
    A standard log + 1 transformation for expression data
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def log_function(self, X):
        return np.log(X + 1)

    def fit_transform(self, X, y=None):
        return self.log_function(X)

    def transform(self, X, y=None):
        return self.log_function(X)


class MinMaxScaler(SKMinMaxScaler):
    """
    Wraper for sklearn's MinMaxScaler so that it returns
    pandas DataFrame when transforming.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def fit_transform(self, X, y=None):
        if type(X) == pd.DataFrame:
            return pd.DataFrame(
                super().fit_transform(X), index=X.index, columns=X.columns
            )
        return super().fit_transform(X)

    def transform(self, X, y=None):
        if type(X) == pd.DataFrame:
            return pd.DataFrame(super().transform(X), index=X.index, columns=X.columns)
        return super().transform(X)
