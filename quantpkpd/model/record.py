""" Model class for PK-PD experiment record """
import numpy as np
from enum import Enum
from pydantic.dataclasses import dataclass

REQ_INPUT_PARAMS = ['datafile_type']

NONMEN_REQ_COLS = ['#ID', 'TIME', 'DV', 'AMT']
NONMEN_PATIENT_MAPPER = {'#ID': 'pid', 'SEX': 'gender',
                         'AGE': 'age', 'HT': 'height', 'WT': 'weight'}


class DoseAdmMtd(Enum):
    """ Enumeration for dosing administration methods """
    EXTRAVASCULAR = 1
    BOLUS = 2
    INFUSION = 3


class Unit(Enum):
    """ variable units in dosing and obseration records """
    MILLIGRAM = 1
    HOUR = 2
    MICROGRAM_LITER = 3


@dataclass
class Patient:
    """ PK/PD patient information """
    pid: int
    gender: str = 'X'
    age: float = 0
    weight: float = 0
    height: float = 0


@dataclass
class DoseRecord:
    """ Dosing record """
    pid: int
    dose_amt: float
    adm_mtd: DoseAdmMtd
    duration: float
    event_dose_ts: dict
    time_unit: Unit = Unit.HOUR
    dose_unit: Unit = Unit.MILLIGRAM
    mol_weight: float = 0


@dataclass
class ObservationRecord:
    """ Obersation Record """
    pid: int
    event_concentration_ts: dict
    event_obs_df: dict = None
    concentration_unit: Unit = Unit.MICROGRAM_LITER
    time_unit: Unit = Unit.HOUR


@dataclass
class ObsDownSlopeRegression:
    """ Regresssion result for observation down slope """
    r_sqrt: float = np.nan  # R2, R-squared
    r_sqrt_adj: float = np.nan  # R2ADJ, Adj. R-squared
    lambda_z_cnt: int = 0  # LAMZPT, number of X points
    lambda_z: float = np.nan  # LAMZ, absolute value of X coefficient
    intercept_b: float = np.nan  # b0, constant coefficient
    corr_xy: float = np.nan  # CORRXY, correlation value of X and Y
    lambda_z_lower_t: float = np.nan  # LANZLL, most lower bound of X
    lambda_z_upper_t: float = np.nan  # LAMZUL, most upper bound of X
    last_predict_y: float = np.nan  # CLSTP, exp(predicted_y[-1])
