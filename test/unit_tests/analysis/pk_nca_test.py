""" test pk_nca """

from quantpkpd.analysis.pk_nca import exec_pk_nca


def test_exec_pk_nca():
    """ test pk_nca executor """

    is_exception = False
    try:
        exec_pk_nca(1, {})
    except TypeError:
        is_exception = True

    assert is_exception

    is_exception = False
    try:
        exec_pk_nca('aaa', {})
    except RuntimeError:
        is_exception = True
    assert is_exception

    is_exception = False
    try:
        exec_pk_nca('test/unit_tests/files/predpp_sample1.csv',
                    {'dose_admin_method': 1})
    except ValueError:
        is_exception = True
    assert is_exception

    exec_pk_nca('test/unit_tests/files/theoph.csv',
                {'dose_admin_method': 1, 'datafile_type': 'nonmen', 'r_sqrt_adj': 0.9})

    assert True
