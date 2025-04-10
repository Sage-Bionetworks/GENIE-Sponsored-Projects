import logging
import pytest
from unittest import mock

import pandas as pd
from pandas.testing import assert_frame_equal
import synapseclient

from geniesp import bpc_redcap_export_mapping as bpc_export

LOGGER = logging.getLogger(__name__)


@pytest.fixture
def mock_syn():
    yield mock.Mock(spec=synapseclient.Synapse)


def test_that_get_drug_variable_names_gets_expected_list():
    var_names = bpc_export.get_drug_variable_names()
    assert var_names == [
        "drugs_drug_1",
        "drugs_drug_oth1",
        "drugs_drug_2",
        "drugs_drug_oth2",
        "drugs_drug_3",
        "drugs_drug_oth3",
        "drugs_drug_4",
        "drugs_drug_oth4",
        "drugs_drug_5",
        "drugs_drug_oth5",
    ]


def test_get_mapping_data_calls_grs_if_use_grs_is_true(mock_syn):
    with mock.patch.object(mock_syn, "get") as mock_get, mock.patch.object(
        pd, "read_csv"
    ):
        bpc_export.get_mapping_data(
            syn=mock_syn, synid_file_grs="synGRS", synid_file_dd="synDD", use_grs=True
        )
        mock_get.assert_called_with("synGRS")


def test_get_mapping_data_calls_dd_if_use_grs_is_false(mock_syn):
    with mock.patch.object(mock_syn, "get") as mock_get, mock.patch.object(
        pd, "read_csv"
    ):
        bpc_export.get_mapping_data(
            syn=mock_syn, synid_file_grs="synGRS", synid_file_dd="synDD", use_grs=False
        )
        mock_get.assert_called_with("synDD")


@pytest.mark.parametrize(
    "input_mapping, var_names, output_mapping",
    [
        (
            pd.DataFrame(
                {
                    "Variable / Field Name": ["drugs_drug_1", "drugs_drug_2"],
                    "Choices, Calculations, OR Slider Labels": [
                        "D001, Aspirin | D002, Ibuprofen | D003, Paracetamol",
                        "D004, Tylenol |",
                    ],
                }
            ),
            ["drugs_drug_1", "drugs_drug_2"],
            {
                "Aspirin": "D001",
                "Ibuprofen": "D002",
                "Paracetamol": "D003",
                "Tylenol": "D004",
            },
        ),
        (
            pd.DataFrame(
                {
                    "Variable / Field Name": ["drugs_drug_1"],
                    "Choices, Calculations, OR Slider Labels": ["D001, Aspirin|"],
                }
            ),
            ["drugs_drug_1"],
            {"Aspirin": "D001"},
        ),
        (
            pd.DataFrame(
                {
                    "Variable / Field Name": ["ethnicity"],
                    "Choices, Calculations, OR Slider Labels": ["1"],
                }
            ),
            ["drugs_drug_1"],
            {},
        ),
        (
            pd.DataFrame(
                {
                    "Variable / Field Name": ["drugs_drug_1"],
                    "Choices, Calculations, OR Slider Labels": [
                        "D001, Aspirin (alternative) | D002, Ibuprofen | D003, Paracetamol"
                    ],
                    "extra_column": ["test1"],
                }
            ),
            ["drugs_drug_1"],
            {"Aspirin": "D001", "Ibuprofen": "D002", "Paracetamol": "D003"},
        ),
    ],
    ids=[
        "multiple_drug_vars",
        "empty_split",
        "nothing_to_parse",
        "parenthesis_split",
    ],
)
def test_that_parse_drug_mappings(input_mapping, var_names, output_mapping):
    result = bpc_export.parse_drug_mappings(mapping=input_mapping, var_names=var_names)
    assert result == output_mapping


@pytest.mark.parametrize(
    "input_data, oncotree_dict, expected_warning",
    [
        (
            pd.DataFrame(
                dict(
                    ONCOTREE_CODE=["Renal Cell Carcinoma", "Renal Clear Cell Carcinoma"]
                )
            ),
            {
                "RCC": {"CANCER_TYPE": "Renal Cell Carcinoma"},
                "OVARY": {"CANCER_TYPE": "Ovarian Cancer"},
            },
            "There are invalid values in ONCOTREE_CODE column in the clinical df: ['Renal Clear Cell Carcinoma', 'Renal Cell Carcinoma']",
        ),
        (
            pd.DataFrame(dict(ONCOTREE_CODE=["Renal Cell Carcinoma", "RCC"])),
            {
                "RCC": {"CANCER_TYPE": "Renal Cell Carcinoma"},
                "OVARY": {"CANCER_TYPE": "Ovarian Cancer"},
            },
            "There are invalid values in ONCOTREE_CODE column in the clinical df: ['Renal Cell Carcinoma']",
        ),
    ],
    ids=["all_invalid", "some_invalid"],
)
def test_that_check_oncotree_codes_gives_expected_warning_when_invalid_codes(
    caplog, input_data, oncotree_dict, expected_warning
):
    with caplog.at_level(logging.WARNING):
        bpc_export.check_oncotree_codes(df=input_data, oncotree_dict=oncotree_dict)
    assert expected_warning in caplog.text


def test_that_check_oncotree_codes_gives_no_warning_when_all_codes_valid(caplog):
    input_data = pd.DataFrame(dict(ONCOTREE_CODE=["RCC", "OVARY"]))
    oncotree_dict = {
        "RCC": {"CANCER_TYPE": "Renal Cell Carcinoma"},
        "OVARY": {"CANCER_TYPE": "Ovarian Cancer"},
    }
    with caplog.at_level(logging.WARNING):
        bpc_export.check_oncotree_codes(df=input_data, oncotree_dict=oncotree_dict)
    assert (
        "There are invalid values in ONCOTREE_CODE column in the clinical df"
        not in caplog.text
    )


def test_that_get_derived_variable_file_gets_file_correctly(mock_syn):
    test_df = pd.DataFrame(
        dict(
            cohort = ["BLADDER", "BrCa", "BLADDER"],
            record_id = ["GENIE-SAGE-1", "GENIE-SAGE-2", "GENIE-SAGE-3"]
            )
    )
    with mock.patch.object(mock_syn, "get") as mock_syn_get, mock.patch.object(
        pd, "read_csv", return_value = test_df
        ) as mock_read_csv:
            output = bpc_export.get_derived_variable_file(
                mock_syn, 
                derived_var_synid = "synZZZZ", 
                cohort = "BLADDER"
                )
            assert_frame_equal(
                output.reset_index(drop=True), pd.DataFrame(
                    dict(
                        cohort = ["BLADDER", "BLADDER"],
                        record_id = ["GENIE-SAGE-1", "GENIE-SAGE-3"]
                        )
                ).reset_index(drop=True),
                check_index_type=False
            )
    

@pytest.mark.parametrize(
    "input, clinical, expected",
    [
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "record_id": ["GENIE-1", "GENIE-1"],
                    "cpt_genie_sample_id": ["GENIE-1-3", "GENIE-1-4"],
                    "cpt_seq_date": ["2017", "2018"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": [None, None],
                }
            )),
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "record_id": ["GENIE-1", "GENIE-1"],
                    "cpt_genie_sample_id": ["GENIE-1-1", "GENIE-1-3"],
                    "cpt_seq_date": ["2017", "2018"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2017", None],
                }
            )),
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "record_id": ["GENIE-1", "GENIE-1"],
                    "cpt_genie_sample_id": ["GENIE-1-1", "GENIE-1-2"],
                    "cpt_seq_date": ["2017", "2018"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2017", "2018"],
                }
            )),
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "record_id": ["GENIE-1", "GENIE-1"],
                    "cpt_genie_sample_id": ["GENIE-1-1", "GENIE-1-2"],
                    "cpt_seq_date": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            )),
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2012", "2013"],
                }
            ),
         pd.DataFrame(
                {
                    "record_id": ["GENIE-1", "GENIE-1", "GENIE-1"],
                    "cpt_genie_sample_id": ["GENIE-1-1", "GENIE-1-2", "GENIE-1-2"],
                    "cpt_seq_date": ["2014", "2015", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ))
        ],
    ids = [
        "none_replaced", 
        "some_replaced", 
        "all_replaced", 
        "the_same", 
        "replacement_has_dups"
        ]
)
def test_that_replace_cpt_seq_date_replaces_correctly_with_derived_variable_replacement_type(input, clinical, expected):
    output = bpc_export.replace_cpt_seq_date(
        input_data = input, 
        replacement_data= clinical,
        cpt_seq_date_replacement_type = "derived_variable"
        )
    assert_frame_equal(output, expected)


@pytest.mark.parametrize(
    "input, clinical, expected",
    [
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "PATIENT_ID": ["GENIE-1", "GENIE-1"],
                    "SAMPLE_ID": ["GENIE-1-3", "GENIE-1-4"],
                    "SEQ_YEAR": ["2017", "2018"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": [None, None],
                }
            )),
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-3"],
                    "SEQ_YEAR": ["2017", "2018"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2017", None],
                }
            )),
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "SEQ_YEAR": ["2017", "2018"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2017", "2018"],
                }
            )),
        (pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "SEQ_YEAR": ["2014", "2015"],
                }
            ),
         pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                    "CPT_SEQ_DATE": ["2014", "2015"],
                }
            ))
        ],
    ids = ["none_replaced_with_extra_cols", "some_replaced", "all_replaced", "the_same"]
)
def test_that_replace_cpt_seq_date_replaces_correctly_with_main_genie_replacement_type(input, clinical, expected):
    output = bpc_export.replace_cpt_seq_date(
        input_data = input, 
        replacement_data= clinical, 
        cpt_seq_date_replacement_type = "main_genie"
        )
    assert_frame_equal(output, expected)
    
    
def test_that_replace_cpt_seq_date_raises_value_error():
    with pytest.raises(
        ValueError, 
        match = "cpt_seq_date_replacement_type: invalid_cpt_seq_date_replacement_type invalid!"
        ):
        output = bpc_export.replace_cpt_seq_date(
            input_data = pd.DataFrame(
                    {
                        "SAMPLE_ID": ["GENIE-1-1", "GENIE-1-2"],
                        "CPT_SEQ_DATE": ["2014", "2015"],
                    }
                ), 
            replacement_data= pd.DataFrame(), 
            cpt_seq_date_replacement_type = "invalid_cpt_seq_date_replacement_type"
            )
