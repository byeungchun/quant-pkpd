""" PK/PD record file handler """
import os
import pandas as pd

from quantpkpd.model.record import (
    Patient,
    DoseRecord,
    DoseAdmMtd,
    ObservationRecord,
    NONMEN_REQ_COLS,
    NONMEN_PATIENT_MAPPER
)


def load_nonmen_format_file(csv_input: str, rec_params: dict) -> dict:
    """
    load NONMEN datset format file. Two formats support - PRED, PREDPP

    Args:
        csv_input (str): csv dataset filepath
        rec_params (dict): user input parameter set

    Return:
        dict: patient, dosing, observation
    """

    df_all = pd.read_csv(csv_input)

    if not all([col in df_all.columns for col in NONMEN_REQ_COLS]):
        raise ValueError(
            f"CSV file does not contain required columns {NONMEN_REQ_COLS}")

    rec_patient = _extract_patient_info(df_all)

    rec_dosing, rec_observation = _extract_dosing_obs_record(
        df_all, rec_params)

    return {'patient': rec_patient, 'dose': rec_dosing, 'obs': rec_observation}


def _extract_patient_info(df_all: pd.DataFrame) -> dict:
    """ extract patient information from NONMEN formatted CSV file"""

    patient_col = {x: v for x, v in NONMEN_PATIENT_MAPPER.items()
                   if x in df_all.columns}
    # extract patient information
    df_patient = df_all[patient_col.keys()].drop_duplicates()
    # rename column name with Patient class variables
    df_patient = df_patient.rename(columns=patient_col)

    if len(df_patient) != len(df_all['#ID'].drop_duplicates()):
        raise ValueError('#ID columns is not unique')

    patients = {}
    for patient in df_patient.to_dict(orient='records'):
        patients[patient['pid']] = Patient(**patient)

    return patients


def _extract_dosing_obs_record(df_all: pd.DataFrame, rec_params: dict):
    """ extract dosing record from NONMEM formatted CSV file"""

    if len(df_all[df_all['DV'] == '.']) == 0:
        df_dose = df_all[NONMEN_REQ_COLS]
        df_obs = df_all[NONMEN_REQ_COLS]
    else:
        df_dose = df_all[df_all['DV'] == '.'][NONMEN_REQ_COLS]
        df_obs = df_all[df_all['DV'] != '.'][NONMEN_REQ_COLS]

    df_dose = df_dose.drop(columns=['DV']).drop_duplicates()
    df_obs = df_obs.drop(columns=['AMT']).drop_duplicates()

    if len(df_dose) == 0:
        raise ValueError('Dose record is not exist in the CSV file')

    if len(df_obs) == 0:
        raise ValueError('Obs record is not exist in the CSV file')

    try:
        df_dose = df_dose.astype(float)
        df_obs = df_obs.astype(float)
    except ValueError as err:
        print('Values in DV or AMT column is not number ', err)
        raise

    dose_records = {}
    obs_records = {}
    for _id in set(df_dose['#ID']):
        _df_dose_pid = df_dose[df_dose['#ID'] == _id]
        _df_obs_pid = df_obs[df_obs['#ID'] == _id]

        if sum(_df_dose_pid['AMT']) <= 0 or sum(_df_obs_pid['DV']) <= 0:
            raise ValueError(
                'Dose, observation value is not proper dosing value')

        _event_dose_ts = _df_dose_pid[['TIME', 'AMT']].to_dict(orient='list')
        dose_records[_id] = DoseRecord(
            pid=_id,
            dose_amt=sum(_df_dose_pid['AMT']),
            adm_mtd=DoseAdmMtd(rec_params['dose_admin_method']),
            duration=sum(_df_dose_pid['TIME']),
            event_dose_ts=_event_dose_ts
            # TODO: add other parameter (ie. time_unit, dose_unit, mol_weight)
        )

        _event_obs_ts = _df_obs_pid[['TIME', 'DV']].to_dict(orient='list')
        obs_records[_id] = ObservationRecord(
            pid=_id,
            event_concentration_ts=_event_obs_ts
            # TODO: add other parameter (ie. event_obs_df, concentration_unit, time_unit)
        )

    return dose_records, obs_records
