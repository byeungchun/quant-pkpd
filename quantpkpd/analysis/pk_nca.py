""" PK NCA (Non Compartment Analysis) """

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from dataclasses import asdict

from quantpkpd.model.record import (
    REQ_INPUT_PARAMS,
    DoseAdmMtd,
    DoseRecord,
    ObservationRecord,
    ObsDownSlopeRegression
)
from quantpkpd.util.nonmen_file_handler import (
    load_nonmen_format_file
)


def exec_pk_nca(datafile: str, params: dict):
    """ main function to execute PK NCA analysis 

        Args:
            datafile (str): dosing, observation data
            params (dict): user parameters
             - datafile_type: datafile format (nonmen)
             - dose_adm_method: dose administration method

    """

    if not isinstance(datafile, str):
        raise TypeError("Nonmen format CSV file path should be string")

    if not os.path.exists(datafile):
        raise RuntimeError(f"File ({datafile}) is not exist")

    if not all([_key in params.keys() for _key in REQ_INPUT_PARAMS]):
        raise ValueError(
            f'params({params.keys()}) not contains required values({REQ_INPUT_PARAMS})')

    if params['datafile_type'] == 'nonmen':
        dataset = load_nonmen_format_file(datafile, params)
    else:
        raise RuntimeError(
            f'Only Nonmen data format available - {params["datafile_type"]}')

    # PK calculation
    analysis_res = {}
    for _id in dataset['patient']:
        dose_rec = dataset['dose'][_id]
        obs_rec = dataset['obs'][_id]
        # STEP 1 - find best down slope in observation
        best_down_regression = _obs_find_best_slope(dose_rec, obs_rec)

    return analysis_res


def _obs_find_best_slope(
    dose_rec: DoseRecord,
    obs_rec: ObservationRecord,
    rsqrt_tolerance: float = 1e-4
):
    """ find best slope in observation

    Args:
        dose_rec
        obs_rec
        rsqrt_tolerance: Tolerance, see Phoneix WinNonlin 6.4 User's Guide p33

    """

    obs_ts = pd.DataFrame(obs_rec.event_concentration_ts)
    obs_ts = obs_ts.loc[obs_ts.DV.dropna().index]  # remove None values
    obs_ts = obs_ts[obs_ts.DV.apply(lambda x: x != 0)]
    obs_ts = obs_ts.sort_values('TIME').reset_index(drop=True)

    res = ObsDownSlopeRegression()

    if len(obs_ts.DV.drop_duplicates()) == 1:
        res.intercept_b = obs_ts.DV.iloc[0]
    else:
        idx_begin = (obs_ts.DV).idxmax()
        if dose_rec.adm_mtd != DoseAdmMtd.BOLUS:
            idx_begin += 1
        idx_end = obs_ts.index[-1]

        # Too small to calculate regresssion
        if len(obs_ts[idx_begin:idx_end+1]) < 2:
            return res

        ols_res = {}
        for i in np.arange(idx_begin, idx_end-1):
            ols_res[i] = _exec_regression(obs_ts.loc[i:idx_end+1])
        best_res = _find_best_rsqrt_more_x_values(ols_res, rsqrt_tolerance)

        if len(best_res) != 0:
            raise RuntimeError(
                'cannot find best down slope regression. check observation event record data')
        res = ObsDownSlopeRegression(**best_res)

    return res


def _find_best_rsqrt_more_x_values(ols_res: dict, rsqrt_tolerance: float):
    """ find best rsqrt and more x values """

    df = pd.DataFrame(ols_res).T
    # STEP1: find maximum r square value
    max_rsqrt = df[df.lambda_z_cnt > 2].r_sqrt_adj.max()
    # STEP2: longest x values from smaller error tolerance
    df2 = df[(max_rsqrt - df.r_sqrt_adj) < rsqrt_tolerance]

    return ols_res[df2.lambda_z_cnt.idxmax()]


def _exec_regression(obs_ts: pd.DataFrame):
    """ calculate a regression of observation event record """
    x = np.array(obs_ts.TIME)
    X = sm.add_constant(x)
    y = np.array(np.log(obs_ts.DV))
    model = sm.OLS(y, X).fit()

    res = ObsDownSlopeRegression(
        r_sqrt=model.rsquared,
        r_sqrt_adj=model.rsquared_adj,
        lambda_z_cnt=len(obs_ts.TIME),
        lambda_z=abs(model.params[1]),
        intercept_b=model.params[0],
        corr_xy=np.corrcoef(x, y)[0, 1],
        lambda_z_lower_t=x[0],
        lambda_z_upper_t=x[-1],
        last_predict_y=np.exp(model.predict()[-1])
    )

    return asdict(res)
