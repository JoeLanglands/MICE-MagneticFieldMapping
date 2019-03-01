


def get_polarity(data_field):
    min_field, max_field = data_field[0].Bz, data_field[0].Bz
    for point in data_field:
        if point.Bz < min_field:
            min_field = point.Bz
        elif point.Bz > max_field:
            max_field = point.Bz
            
    if abs(min_field) > abs(max_field):
        return -1.0
    elif abs(min_field) < abs(max_field):
        return 1.0



def cut_field_for_fit(data_field, min_bz=0.005):
    """Function to cut the data field appropriately for the geofit and coilfit to save run time.
    This function takes a field as a list of measurements and cuts parts of it off so that it is 
    smaller.  This will save on run time for the fitting process.  It cuts off the low field regions
    whilst also removing data from the central probe since it can't be trusted!
    """
    if data_field == None:
        return None
    new_field = []
    for f in data_field:
        if f.Bz > min_bz and f.sensorNumber != 0:
            new_field.append(f)

    return new_field
