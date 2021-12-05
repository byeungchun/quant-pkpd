""" record Model test """

from pydantic import ValidationError
from quantpkpd.model.record import Patient


def test_record():
    """test patient class"""
    patient = Patient(1, 'M', 43, 70.1, 171.1)

    assert patient.pid == 1
    assert patient.gender == 'M'
    assert patient.age == 43
    assert patient.weight == 70.1

    patient = Patient('1', 'M', 43, 70.1, 171.1)

    assert isinstance(patient.pid, int)

    is_exception = False
    try:
        patient = Patient('aaa', 'M', 43, 70, 127)
    except ValidationError:
        is_exception = True

    assert is_exception

    patient = Patient(1)

    assert patient.gender == 'X'
