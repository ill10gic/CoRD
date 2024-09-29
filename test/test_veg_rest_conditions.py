import pytest
from cord.ripcas_dflow import determine_veg_reset_value
from pandas import read_excel

@pytest.mark.parametrize("input_case", [
    (100, 100, 1000, 1000),
    (100, 100, 1001, 1000),
    (100, 100, 3999, 1000),
    (100, 200, 1000, 1000),
    (100, 300, 1000, 1000),
    (200, 100, 1000, 1000),
    (200, 100, 3000, 1000),
    (200, 200, 3000, 1000),
    (200, 200, 2500, 2000),
    (400, 300, 3999, 2000)
])
def test_determine_veg_reset_value(input_case):
    df = read_excel("test/data/veg_roughness_shearres.xlsx")
    
    #q, hbfl, veg_min, veg_max, expected_reset_val    
    assert input_case[3] == determine_veg_reset_value(input_case[0], input_case[1], input_case[2], df)