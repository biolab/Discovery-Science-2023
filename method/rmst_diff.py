import numpy as np
import pandas as pd

from itertools import combinations

from lifelines import KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time


def fit_KM_sequence(sample_indicator: np.ndarray, event_indicator: np.ndarray):
    """
    Fits a sequence of survival probabilities where each step is unique timepoint from the data.
    The method assumes that samples are ordered by their survival time.

    Arguments
    ---------
    sample_indicator: np.ndarray (n_samples,)
        Indicator with 1 if sample is in cohort and 0 if it is not.
    event_indicator: np.ndarray (n_samples,)
        Indicator with 1 if event happened and 0 if it is censored.

    Returns
    -------
    np.ndarray (n_samples,)
        Survival probability of the cohort where indices correspond to the time points in the data.
    """
    n_series = np.cumsum(sample_indicator[::-1])[::-1]
    return np.append(
        [1], np.cumprod((n_series - event_indicator * sample_indicator) / n_series)
    )


def difference_RMST(km1: np.ndarray, km2: np.ndarray, time_values: np.ndarray):
    """
    Calculates the difference between the KM1 and KM1 restricted mean survival time.
    Basically integrates over the curves.
    KM curves are of equal length by design.

    Arguments
    ---------
    km1: np.ndarray (n_samples, )
        KM 1 curve with survival probabilities for each time point.
    km2: np.ndarray (n_samples, )
        KM 2 curve with survival probabilities for each time point.
    time_values: np.ndarray (n_samples, )
        Array of survival time points in the dataset

    Returns
    -------
    float
        The difference in RMST (AOC of difference in RMST)
    """
    LEN = len(time_values)  # the last value is the maximal (resticted at) time
    dt = np.ediff1d(time_values, to_begin=time_values[0])
    return np.abs(np.sum((km1[:LEN] - km2[:LEN]) * dt))




# def split_by_median(feature, feature_values, durations, events, TIME_LIMIT):
#     cutoff = feature_values.median()
#     strata = (feature_values > cutoff).astype(bool)

#     time_limit = min(durations[strata].max(), durations[~strata].max())

#     limit_is_lower = True
#     if TIME_LIMIT <= time_limit:
#         time_limit = TIME_LIMIT
#     else:
#         print("TIME_LIMIT is lower than maximal time in one of the cohorts.")
#         limit_is_lower = False

#     kmf1 = KaplanMeierFitter().fit(durations[strata], events[strata])
#     rmst_1 = restricted_mean_survival_time(kmf1, t=time_limit)

#     kmf2 = KaplanMeierFitter().fit(durations[~strata], events[~strata])
#     rmst_2 = restricted_mean_survival_time(kmf2, t=time_limit)

#     total_rmst = abs(rmst_1 - rmst_2)

#     return (feature, total_rmst, time_limit, limit_is_lower)


def fast_split_by_median(feature, feature_values, sorted_time, sorted_events, TIME_LIMIT):
    cutoff = feature_values.median()
    strata = (feature_values > cutoff).astype(bool)
    time_limit = min(sorted_time[strata].max(), sorted_time[~strata].max())

    limit_is_lower = True
    if TIME_LIMIT <= time_limit:
        time_limit = TIME_LIMIT
    else:
        print("TIME_LIMIT is lower than maximal time in one of the cohorts.")
        limit_is_lower = False

    km1 = fit_KM_sequence(strata, sorted_events)
    km2 = fit_KM_sequence(~strata, sorted_events)
    dif_rstm = difference_RMST(km1, km2, list(sorted_time[sorted_time <= time_limit]) + [time_limit])

    return (feature, dif_rstm, limit_is_lower)


def rmst_diff(X: pd.DataFrame, y=None, return_all:bool=False):
        interaction = []

        y.loc[:, 'time'] = y['time'] + np.arange(len(y))*1e-5
        y.sort_values(by=['time'], inplace=True)
        X = X.loc[y.index] # sort X by y index

        time_col = y['time'].values 
        event_col = y['event'].values

        TIME_LIMIT = np.percentile(time_col, 75)


        for feature_comb in combinations(X.columns, 2):
            feature1, feature2 = feature_comb

            _, feature1_rmst, lil_f1 = fast_split_by_median(f'{feature1}', X[feature1], time_col, event_col, TIME_LIMIT)
            _, feature2_rmst, lil_f2 = fast_split_by_median(f'{feature2}', X[feature2], time_col, event_col, TIME_LIMIT)
            
            _, interaction_add_rmst, lil_f1p2 = fast_split_by_median(f'{feature1}+{feature2}', X[feature1] + X[feature2], time_col, event_col, TIME_LIMIT)
            _, interaction_diff_rmst, lil_f1m2 = fast_split_by_median(f'{feature1}-{feature2}', X[feature1] - X[feature2], time_col, event_col, TIME_LIMIT)
            _, interaction_mult_rmst, lil_f1i2 = fast_split_by_median(f'{feature1}*{feature2}', X[feature1] * X[feature2], time_col, event_col, TIME_LIMIT)

            limit_is_lower = True in [lil_f1, lil_f2, lil_f1p2, lil_f1m2, lil_f1i2]

            temp = np.array([feature1_rmst, feature2_rmst])
            # self.interaction.append((interaction_term, interaction_mult_rmst, np.abs(interaction_mult_rmst - np.max(temp))))

            hits = [
                interaction_add_rmst if return_all or np.all(interaction_add_rmst > temp) else np.nan,
                interaction_diff_rmst if return_all or np.all(interaction_diff_rmst > temp) else np.nan,
                interaction_mult_rmst if return_all or np.all(interaction_mult_rmst > temp) else np.nan,
            ]

            if np.isnan(hits).all():
                continue
            else:
                interaction.append([f'{feature1}*{feature2}'] + [feature1_rmst, feature2_rmst] + hits + [limit_is_lower])
            
            # if np.all(interaction_mult_rmst > temp):
            #   self.interaction.append((interaction_term, interaction_mult_rmst, np.abs(interaction_mult_rmst - np.max(temp)), limit_is_lower))

        
        # self.interaction = sorted(self.interaction, key=lambda x: x[1], reverse=True)
        return interaction