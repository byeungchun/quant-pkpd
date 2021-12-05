""" Recard file handler test """

from quantpkpd.util.nonmen_file_handler import load_nonmen_format_file


def test_load_nonmen_format_file():
    """ test load_nonmen_format_file function"""

    records = load_nonmen_format_file(
        'test/unit_tests/files/predpp_sample1.csv',
        {'dose_admin_method': 2})

    assert isinstance(records, dict)
    assert len(records.keys()) == 3
