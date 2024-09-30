import pytest
import copy
import six
# from cord import ESRIAsc
from cord.ripcas_dflow import ripcas, ESRIAsc, determine_veg_reset_value, veg2n, Pol
from pandas import read_excel

@pytest.mark.parametrize("input_case", [
    # (100, 100, 1000, 1000),
    # (100, 100, 1001, 1000),
    # (100, 100, 3999, 1000),
    # (100, 200, 1000, 1000),
    # (100, 300, 1000, 1000),
    # (200, 100, 1000, 1000),
    # (200, 100, 3005, 1000),
    # (200, 200, 3000, 1000),
    # (200, 200, 2500, 2000),
    # (400, 300, 3999, 2000),
    # (-9999, -9999, 1200, 2000),
    # (300, 200, 1007, 1000),
    (200, 200, 3007, 1000),
    # (200,100,3007, 1000),
    # (300, 200, 3007, 3000)
    #(qzone,hbfl,veg-map-value, expected reset value)
])
def test_determine_veg_reset_value(input_case):
    pandas_dataframe = read_excel("test/data/veg_roughness_shearres_update_veg_rules_colormap.xlsx")
    
    for idx in range(0,200000):
        #q, hbfl, veg_min, veg_max, expected_reset_val    
        return_vale=  determine_veg_reset_value(
            input_case[0], #qzone 
            input_case[1], #hbfl
            input_case[2], #veg-map-value
            pandas_dataframe) #excel data
        print("return value is: {}".format(return_vale))
        assert input_case[3] == return_vale
    
def test_veg2n_and_polygon_code():
    # convert the vegetation .asc to roughness .pol, write to veg_path
    Pol.from_ascii(veg2n(ESRIAsc("test/data/Vegmap_1000s_2024jun.asc"), "test/data/veg_roughness_shearres_update_veg_rules_colormap.xlsx", 0.035)).write('test/data/test_veg2n_and_polygon_code.asc')
    
def test_table():
    pandas_dataframe = read_excel("test/data/veg_roughness_shearres_update_veg_rules_colormap.xlsx")
    
    q_zone_map_column = pandas_dataframe["Q_zone_map"]
    hbfl_zone_map_column = pandas_dataframe["HBFL_zone_map"]
    veg_input_map_min_column = pandas_dataframe["veg_input_map_min"]
    veg_input_map_max_column = pandas_dataframe["veg_input_map_max"]
    vegetation_reset_value_column = pandas_dataframe["vegetation_reset_value"]
    print("hello")
    # this is actually zero indexed - it skips the header
    for i in range(0, 36):  # TODO: set this to figure out the number of conditions dynamically
        print("hello")
        print("row {}".format(i))
        print("{},{},{},{},{}".format(q_zone_map_column.iloc[i], hbfl_zone_map_column.iloc[i],veg_input_map_min_column.iloc[i], veg_input_map_max_column.iloc[i], vegetation_reset_value_column.iloc[i]))
    assert False
    
def test_ripcas():
    #shear asc
    shear_asc_path = "test/data/test-ripcas-real-dflow-3-output/shear_out.asc"
    previous_veg_asc_path = "test/data/test-ripcas-real-dflow-2-output/vegetation.asc"
    current_veg3_ACTUAL_from_wheeler_path = "test/data/test-ripcas-real-dflow-3-output/vegetation.asc"
    ripcas_expected = "test/data/test-ripcas-real-dflow-3-output/ripcas_expected.asc"
    
    actual_asc = ESRIAsc(current_veg3_ACTUAL_from_wheeler_path) 
    expected_asc = ESRIAsc(ripcas_expected)
    are_equal = True
    for idx in range(len(actual_asc.data)):
        if (actual_asc.data[idx] != expected_asc.data[idx]):
            are_equal = False
            print("not equal at idx: {}".format(idx))
            print("actual != expected: {} != {}".format(actual_asc.data[idx], expected_asc.data[idx]))
    assert are_equals
    
def test_ripcas_with_hardcoded_version():
    #shear asc
    shear_asc = ESRIAsc("test/data/test-ripcas-real-dflow-3-output/shear_out.asc")
    previous_veg_asc = ESRIAsc("test/data/test-ripcas-real-dflow-2-output/vegetation.asc")
    q_zone_map = ESRIAsc("test/data/rastert_rastert5_Q_zone.asc")
    hbfl_map = ESRIAsc("test/data/rastert_boundar1_HBFL_zone.asc")
    excel_rules_dataframe = "test/data/veg_roughness_shearres_update_veg_rules_colormap.xlsx"
    hardcoded_ripcas_output = ripcas_with_hardcoded(previous_veg_asc, q_zone_map, hbfl_map, shear_asc, excel_rules_dataframe)
    current_ripcas_output = ripcas(previous_veg_asc, q_zone_map, hbfl_map, shear_asc, excel_rules_dataframe)
    are_equals = True
    for idx in range(len(hardcoded_ripcas_output.data)):
        if (hardcoded_ripcas_output.data[idx] != current_ripcas_output.data[idx]):
            are_equal = False
            break
            #print("not equal at idx: {}".format(idx))
            #print("actual != expected: {} != {}".format(hardcoded_ripcas_output.data[idx], current_ripcas_output.data[idx]))
    assert are_equals
    hardcoded_ripcas_output.write_unflattened_asc("test/data/test-ripcas-real-dflow-3-output/ripcas_expected.asc")
    read_hardcoded_ripcas = ESRIAsc("test/data/test-ripcas-real-dflow-3-output/ripcas_expected.asc")
    assert read_hardcoded_ripcas == hardcoded_ripcas_output

def test_make_ACTUAL_ripcas_output_byhand():
    shear_asc = ESRIAsc("test/data/test-ripcas-real-dflow-3-output/shear_out.asc")
    previous_veg_asc = ESRIAsc("test/data/test-ripcas-real-dflow-2-output/vegetation.asc")
    q_zone_map = ESRIAsc("test/data/rastert_rastert5_Q_zone.asc")
    hbfl_map = ESRIAsc("test/data/rastert_boundar1_HBFL_zone.asc")
    excel_rules_dataframe = read_excel("test/data/veg_roughness_shearres_update_veg_rules_colormap.xlsx")
    accurate_ripcas_output = ripcas_with_hardcoded(shear_asc, previous_veg_asc, q_zone_map, hbfl_map, excel_rules_dataframe)
    accurate_ripcas_output.write_unflattened_asc("test/data/test-ripcas-real-dflow-3-output/ripcas_expected.asc")

def determine_veg_reset_value_test(zone_map_value, hbfl_map_value, vegetation_map_value, ripcas_excel_dataframe):
    """ Uses recruitment rules to determine reset value for vegatation code

    Args:
        zone_map_value (number): Q Zone map value of the current location being processed in recruitment rules. 
        hbfl_map_value (number): HBFL zone map value of the current location being processed in recruitment rules
        vegetation_map_value (number): previous iteration's vegetation value of the current location being processed in recruitment rules
        ripcas_excel_dataframe (dataframe): all the recruitment rules and variables from the excel spreadsheet. Values are associate with vegetation succession model 
    Returns:
        (number) the reset vegetation code determined by the rule set
    """
    
    if (zone_map_value == 100 and hbfl_map_value == 100 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 1000
    elif (zone_map_value == 100 and hbfl_map_value == 200 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 1000
    elif (zone_map_value == 100 and hbfl_map_value == 300 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 1000
    elif (zone_map_value == 100 and hbfl_map_value == -9999 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 1000
    elif (zone_map_value == 200 and hbfl_map_value == 100 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 1000
    elif (zone_map_value == 200 and hbfl_map_value == 100 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 1000
    elif (zone_map_value == 200 and hbfl_map_value == 100 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    elif (zone_map_value == 200 and hbfl_map_value == 200 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 1000
    elif (zone_map_value == 200 and hbfl_map_value == 200 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 1000
    elif (zone_map_value == 200 and hbfl_map_value == 200 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    elif (zone_map_value == 200 and hbfl_map_value == 300 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 1000
    elif (zone_map_value == 200 and hbfl_map_value == -9999 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 1000
    #13
    elif (zone_map_value == 200 and hbfl_map_value == 300 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 3000
    elif (zone_map_value == 200 and hbfl_map_value == -9999 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 3000
    elif (zone_map_value == 200 and hbfl_map_value == 300 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    elif (zone_map_value == 200 and hbfl_map_value == -9999 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    elif (zone_map_value == 300 and hbfl_map_value == 100 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 1000
    elif (zone_map_value == 300 and hbfl_map_value == 100 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    elif (zone_map_value == 300 and hbfl_map_value == 100 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 3000
    #20
    elif (zone_map_value == 300 and hbfl_map_value == 200 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 1000
    elif (zone_map_value == 300 and hbfl_map_value == 200 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    elif (zone_map_value == 300 and hbfl_map_value == 200 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 3000
    elif (zone_map_value == 300 and hbfl_map_value == 300 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 3000
    elif (zone_map_value == 300 and hbfl_map_value == -9999 and vegetation_map_value >= 1000 and vegetation_map_value <= 1999):
        return 3000
    elif (zone_map_value == 300 and hbfl_map_value == 300 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    elif (zone_map_value == 300 and hbfl_map_value == -9999 and vegetation_map_value >= 2000 and vegetation_map_value <= 2999):
        return 2000
    #27
    elif (zone_map_value == 300 and hbfl_map_value == 300 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 2000
    elif (zone_map_value == 300 and hbfl_map_value == -9999 and vegetation_map_value >= 3000 and vegetation_map_value <= 3999):
        return 2000
    elif (zone_map_value == 400 and hbfl_map_value == 100 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000
    elif (zone_map_value == -9999 and hbfl_map_value == 100 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000
    elif (zone_map_value == 400 and hbfl_map_value == 200 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000
    elif (zone_map_value == -9999 and hbfl_map_value == 200 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000
    #33
    elif (zone_map_value == 400 and hbfl_map_value == 300 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000
    elif (zone_map_value == 400 and hbfl_map_value == -9999 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000
    elif (zone_map_value == -9999 and hbfl_map_value == 300 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000
    #36
    elif (zone_map_value == -9999 and hbfl_map_value == -9999 and vegetation_map_value >= 1000 and vegetation_map_value <= 3999):
        return 2000    
    raise RuntimeError('Veg reset value not found with given rules')

def ripcas_with_hardcoded(vegetation_map, zone_map, hbfl_map, shear_map, ripcas_required_data):
    """
    Simple version of the CASiMiR model for vegetation succession. Before the
    model is run, we check that all the unique values from vegetation_map are
    present in the shear_resistance_dict. Otherwise the process will fail
    wherever the vegetation map value is not present in the dictionary
    on lookup.

    Arguments:
        vegetation_map (str or ESRIAsc): location on disk or ESRIAsc
            representation of the vegetation map
        zone_map (str or ESRIAsc): location on disk or ESRIAsc representation
            of the zone map.
        hbfl_map (str or ESRIAsc): location on disk or ESRIAsc representation
            of the hbfl map.
        shear_map (str or ESRIAsc): location on disk or ESRIAsc representation
            of the shear stress map
        ripcas_required_data (str): Excel spreadsheet of data needed for
            ripcas run. Encompasses the landscape model for the watershed.
            Can have one or two 'Code' columns and must have exactly one
            'shear_resis' column and exactly one 'n_val' column

    Returns:
        (ESRIAsc) vegetation map updated with new values corresponding to
            succession rules
    """
    print('in ripcas')
    print('vegetation_map')
    print(vegetation_map)
    print('zone_map')
    print(zone_map)
    print('hbfl_map')
    print(hbfl_map)
    # Check that veg, shear, and zone maps are the right type (either string or ESRIAsc)
    if isinstance(vegetation_map, six.string_types):
        vegetation_map = ESRIAsc(vegetation_map)
    elif not isinstance(vegetation_map, ESRIAsc):
        raise TypeError('vegetation_map must be type str or ESRIAsc')

    if isinstance(shear_map, six.string_types):
        shear_map = ESRIAsc(shear_map)
    elif not isinstance(shear_map, ESRIAsc):
        raise TypeError('shear_map must be type str or ESRIAsc')

    if isinstance(zone_map, six.string_types):
        zone_map = ESRIAsc(zone_map)

    if isinstance(hbfl_map, six.string_types):
        hbfl_map = ESRIAsc(hbfl_map)

    elif not isinstance(zone_map, ESRIAsc):
        raise TypeError('zone_map must be type str of ESRIAsc, not ' +
                        str(type(zone_map)))

    # Check that all the parameters from Excel are imported correctly into pandas data frame
    if isinstance(ripcas_required_data, six.string_types):
        cas_df = read_excel(ripcas_required_data)

        # sanity check to make sure our lookup is correct
        #vegetation and shear stress columns
        assert 'Code.1' in cas_df  # TODO generalize later
        assert 'shear_resis' in cas_df
        
        # succession ruleset columns
        assert "Q_zone_map" in cas_df
        assert "HBFL_zone_map" in cas_df
        assert "veg_input_map_min" in cas_df
        assert "veg_input_map_max" in cas_df
        assert "vegetation_reset_value" in cas_df

        shear_resistance_dict = dict(
            # TODO verify this is the "2nd instance of Code column" I.E. zero indexed
            zip(cas_df['Code.1'], cas_df['shear_resis'])
        )
        

    else:
        raise TypeError('shear_resistance_dict must be type str')

    # initialize the vegetation map that will be returned
    ret_veg_map = copy.deepcopy(vegetation_map)

    # run through each cell of the ESRI ascii file
    for idx in range(len(shear_map.data)):
        # determine whether or not the vegetation should be reset to age zero
        veg_val = int(vegetation_map.data[idx])
        if veg_val != 0:

            # Make sure shear map and veg map at this cell have data
            # To Do: Enable vegetation to age, even if shear at that cell is nodata
            is_not_nodata = (
                shear_map.data[idx] != shear_map.NODATA_value and
                vegetation_map.data[idx] != vegetation_map.NODATA_value and
                vegetation_map.data[idx] > 0
            )

            # set reset to true if shear is over threshold
            veg_needs_reset = ( # TODO no change needed here
                is_not_nodata and
                shear_map.data[idx] > shear_resistance_dict[veg_val] # TODO where do these originate or are read from (files)?
                # TODO HBFL is "old zone map", will use a new Q Zone in addition
                
            )

            if veg_needs_reset:
                # reset vegetation to age zero with veg type appropriate to zone
                # ret_veg_map.data[idx] = zone_map.data[idx] # TODO - this line is a problem 
                ret_veg_map.data[idx] = determine_veg_reset_value_test(zone_map.data[idx], hbfl_map.data[idx], vegetation_map.data[idx], cas_df)

            # whether or not the vegetation was destroyed, age by one
            if is_not_nodata:
                ret_veg_map.data[idx] += 1

    return ret_veg_map