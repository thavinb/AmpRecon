import pytest
import json
from grc_process_mccoil_io import __getSNPnumbersFrom, loadBarcodeDef, loadSamplesBarcode, writeMcCOILat
from WriteCOIgrc import parseCOIout, writeCOIgrc
import os

# ---- BARCODE2McCOILin unit tests -----

### test __getSNPNumbersFrom

@pytest.fixture
def barcode_def_dct():
    return {
        1: {"Chromosome": "Pf3D7_02_v3", "Locus": 376222, "RefAllele": "A"},
        2: {"Chromosome": "Pf3D7_02_v3", "Locus": 376223, "RefAllele": "T"},
        3: {"Chromosome": "Pf3D7_02_v3", "Locus": 376224, "RefAllele": "C"},
        4: {"Chromosome": "Pf3D7_02_v3", "Locus": 376225, "RefAllele": "G"},
        5: {"Chromosome": "Pf3D7_02_v3", "Locus": 376226, "RefAllele": "T"},
    }

def test_valid_barcode_returns_list_of_numbers(barcode_def_dct):
    """
    test if function returns a list of numbers
    """
    barcode = "ATCGT"
    result = __getSNPnumbersFrom(barcode, barcode_def_dct)
    assert isinstance(result, list)
    assert len(result) == 5
    assert all((isinstance(x, float) or isinstance(x,int)) for x in result)


def test_invalid_characters_raise_error(barcode_def_dct):
    """
    if non ATCGXN char in the barcode, then it should raise an error
    """
    barcode = "ATC2G"
    with pytest.raises(Exception):
        __getSNPnumbersFrom(barcode, barcode_def_dct)

def test_reference_allele_set_as_N_raises_error(barcode_def_dct):
    """
    if reference allele is set as 'N', it should raise an error
    """
    barcode_def_dct[1]["RefAllele"] = "N"
    barcode = "ATCGT"
    with pytest.raises(Exception):
        __getSNPnumbersFrom(barcode, barcode_def_dct)

def test_missing_positions_are_assigned_negative_one(barcode_def_dct):
    """
    If barcode base is 'X' (missing), it should be assign as -1
    """
    barcode = "ATCGX"
    result = __getSNPnumbersFrom(barcode, barcode_def_dct)
    assert result[4] == -1

def test_heterozygous_positions_are_assigned_point_five(barcode_def_dct):
    """
    If barcode base is 'N' (heterozygoues), it should be assign as 0.5
    """
    barcode = "ATCGN"
    result = __getSNPnumbersFrom(barcode, barcode_def_dct)
    assert result[4] == 0.5

def test_homozygous_minor_positions_are_assigned_zero(barcode_def_dct):
    """
    If barcode base is equal to its ref allele, it should be assign as 0
    """
    barcode = "ATCGT"
    result = __getSNPnumbersFrom(barcode, barcode_def_dct)
    assert result == [0,0,0,0,0]

def test_homozygous_major_positions_are_assigned_one(barcode_def_dct):
    """
    If barcode base is not equal to its ref allele, it should be assign as 1
    """
    barcode = "GTCAT"
    result = __getSNPnumbersFrom(barcode, barcode_def_dct)
    assert result == [1,0,0,1,0]


### test loading json barcode definition

@pytest.fixture
def barcode_def_json(tmp_path):
    # Set up a fixture for the barcode definition dictionary
    barcode_def = {
        "barcode_ref": {
            "001": {"Chromosome": "Pf3D7_02_v3", "Locus": 376222, "RefAllele": "A"},
            "002": {"Chromosome": "Pf3D7_02_v3", "Locus": 470013, "RefAllele": "G"},
            "003": {"Chromosome": "Pf3D7_03_v3", "Locus": 656861, "RefAllele": "T"}
        }
    }
    barcode_def_file = tmp_path / "barcode_def_test.json"
    barcode_def_file.write_text(json.dumps(barcode_def))

    return barcode_def["barcode_ref"]

def test_loadBarcodeDef(barcode_def_json, tmp_path):
    '''
    Test if loading a barcode definition dictionary from a JSON file
    returns the expected result
    '''
    barcode_def_file = tmp_path / "barcode_def_test.json"
    barcode_def_file.write_text(json.dumps({"barcode_ref": barcode_def_json}))
    assert loadBarcodeDef(str(barcode_def_file)) == barcode_def_json

def test_loadBarcodeDef(barcode_def_json, tmp_path):
    barcode_def_file = tmp_path  / "barcode_def_test.json"
    barcode_def_file.write_text(json.dumps({"barcode_ref":barcode_def_json}))

    barcode_def = loadBarcodeDef(str(barcode_def_file))

    assert set(barcode_def.keys()) == set([1, 2, 3]), "Unexpected keys in barcode definition"
    assert all(isinstance(k, int) for k in barcode_def.keys()), "Keys should be integers"
    assert all(isinstance(v, dict) for v in barcode_def.values()), "Values should be dictionaries"
    assert all(isinstance(v['Chromosome'], str) for v in barcode_def.values()), "Chromosome values should be strings"
    assert all(isinstance(v['Locus'], int) for v in barcode_def.values()), "Locus values should be integers"
    assert all(isinstance(v['RefAllele'], str) for v in barcode_def.values()), "RefAllele values should be strings"


### test loadSamplesBarcode()

@pytest.fixture
def samples_brcd_path(tmpdir):
    # Create temporary directory and file for testing
    test_file = tmpdir.join("test_samples.tsv")
    test_file.write("Sample\tBarcode\nPT59269\tAGTCAGTGCAAACGTCTAA\nSPT4416\tXXXXAXNXXAXXXXXXXXX\n")
    return str(test_file)

def test_loadSamplesBarcode_output_type(samples_brcd_path):
    # Test if the output is a dictionary
    output = loadSamplesBarcode(samples_brcd_path)
    assert type(output) == dict

def test_loadSamplesBarcode_keys_exist(samples_brcd_path):
    # Test if expected keys exist in the output
    output = loadSamplesBarcode(samples_brcd_path)
    assert all([key in output for key in ["PT59269", "SPT4416"]])

def test_loadSamplesBarcode_barcode_exist(samples_brcd_path):
    # Test if barcode key exists in output
    output = loadSamplesBarcode(samples_brcd_path)
    assert all(["barcode" in value for value in output.values()])

def test_loadSamplesBarcode_barcode_type(samples_brcd_path):
    # Test if barcode value is a string
    output = loadSamplesBarcode(samples_brcd_path)
    assert all([type(value["barcode"]) == str for value in output.values()])

def test_loadSamplesBarcode_error_on_wrong_columns(samples_brcd_path):
    # Test if an AssertionError is raised if the input file contains more than 2 columns
    with open(samples_brcd_path, "a") as barcode_infl:
        barcode_infl.write("SPT9999\tNNNNN\tExtraColumn\n")
    with pytest.raises(AssertionError):
        loadSamplesBarcode(samples_brcd_path)

### test writeMcCOILat()
@pytest.fixture
def barcodes_samples_dct():
    return {
        'PT59269': {'barcode':'ATCGT', 'coi_in': [0,0,0,0,0]},
        'SPT4416': {'barcode':'XXXGT', 'coi_in': [-1,-1,-1,0,0]},
        'SPT4417': {'coi_in':'ATTGN', 'coi_in':[0,0,1,0,0.5]}
    }

@pytest.fixture
def out_flpath(tmpdir):
    return os.path.join(tmpdir, "test_output.txt")

def test_writeMcCOILat(barcode_def_dct, barcodes_samples_dct, out_flpath):
    writeMcCOILat(out_flpath, barcode_def_dct, barcodes_samples_dct)
    assert os.path.exists(out_flpath)
    with open(out_flpath, "r") as f:
        # check header
        header = f.readline().strip().split("\t")
        assert header == ["sample_id", "Pf3D7_02_v3_376222", "Pf3D7_02_v3_376223", 
                        "Pf3D7_02_v3_376224", "Pf3D7_02_v3_376225", "Pf3D7_02_v3_376226"]
        # check rows
        rows = f.readlines()
        assert len(rows) == 3
        assert rows[0].strip().split("\t") == ["PT59269", "0", "0", "0","0","0"]
        assert rows[1].strip().split("\t") == ["SPT4416", "-1", "-1", "-1","0","0"]
        assert rows[2].strip().split("\t") == ["SPT4417", "0", "0", "1","0","0.5"]
# --------------------------------------

def test_parseCOIout(tmp_path):
    # create temporary file with test data
    coiout_path = os.path.join(tmp_path, "test_coiout.tsv")
    with open(coiout_path, "w") as f:
        f.write("file\tCorP\tname\tmean\tmedian\tsd\tquantile0.025\tquantile0.975\n")
        f.write("COIout\tC\tSPT59269\t1\t1\t0\t1\t1\n")
        f.write("COIout\tC\tSPT42677\t2\t2\t0.48162\t1\t2\n")
        f.write("COIout\tC\tSPT31022\t3\t3\t0.13761\t3\t3\n")
        f.write("COIout\tC\tSPT32861\t2\t2\t0.00636\t2\t2\n")
    
    # call the function and check the output
    expected_output = {"SPT59269": {"CorP": "C", "name": "SPT59269", "mean": "1", 
                                    "median": "1", "sd": "0", "quantile0.025": "1",
                                    "quantile0.975": "1"},
                       "SPT42677": {"CorP": "C", "name": "SPT42677", "mean": "2",
                                    "median": "2", "sd": "0.48162", "quantile0.025": "1",
                                    "quantile0.975": "2"},
                       "SPT31022": {"CorP": "C", "name": "SPT31022", "mean": "3",
                                     "median": "3", "sd": "0.13761", "quantile0.025": "3",
                                    "quantile0.975": "3"},
                       "SPT32861": {"CorP": "C", "name": "SPT32861", "mean": "2",
                                    "median": "2", "sd": "0.00636", "quantile0.025": "2",
                                    "quantile0.975": "2"}
                    }
    assert parseCOIout(coiout_path) == expected_output