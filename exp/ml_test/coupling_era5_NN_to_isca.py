import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import tensorflow as tf
from netCDF4 import Dataset
from keras.models import Model, Sequential
from keras.layers import Input, Dense, LeakyReLU, Dropout
from keras.layers import LeakyReLU
from keras.optimizers import Adam
from keras.models import model_from_json
import csv
import xarray as xr

# The ERA5 data is gridded to match Isca's grid (lon:128, lat:64) so the ML output should also match this  
# Which is why I'm largely ignoring what Joe did

# I'm assuming you would only do this for one month at a time so dataset is at a fixed point in time
# I have code to loop through times but I'm assuming that's unnecessary because it should be a monthly averaged atmos_monthly.nc file I think?

def create_input_file_from_ANN(exp_name, month_num):
    data_dir = '/emmy-noether/home/sit204/isca_data'
    file_list = [f'{data_dir}/{exp_name}/run{month_num:04d}/atmos_half_monthly.nc']
    dataset = xr.open_mfdataset(file_list, decode_times = False)
    
    sdor_dataset = xr.open_dataset('/home/links/sit204/Isca/exp/ml_test/input/2021_t42_std.nc')

    # get rid of time dimension so the reshaping works later (there's only one time anyway)
    dataset = dataset.sel(time = dataset.time.values[0])

    
    # Need to select vertical Level Closest to the surface? by taking dataset and selecting the surface using model height rather than pfull hopefully
    data = dataset.sel(pfull = 1000., method = 'nearest')  
    # !!SELECT lowest /surface height!! assuming that you get that to be outputted in atmos_monthly.nc
    
    # All of this assumes you're working with data in the plane / at the surface and the output of this function is t2m_std for the
    # surface plane (which is a 64 x 128 array that should match Isca's lat lon grid)

    # So then you end up with a dataset that is for a fixed time and height that is only 64 x 128 with Isca's grid spacing - I
    # assume it's called data

    # Can use min, max imported from training to normalise or calculate min/max from Isca data (decide later)
    def normalise_with_training(datain, training_min , training_max):
        dataout=(datain-training_min)/(training_max-training_min)
        return dataout;

    # Load data about min and max feature values from training (for normalising)

    feature_min = {}
    # Open the CSV file for reading
    with open('training_min.csv', 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        feature_min = next(reader)

    feature_max = {}
    # Open the CSV file for reading
    with open('training_max.csv', 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        feature_max = next(reader)


    def unnormalise(normalised_data, training_min , training_max):
        return ((training_max - training_min) *(normalised_data)) + training_min


    # this function will reshape the features
    # Here df is a surface map of a single variable at a single time at Isca's resolution (should be 64 x 128)
    def reshape(df): 
        # Reshape to 1d vector
        var = np.reshape(df.values,(64*128, 1))
        return var 

    # list the input features (changed from era5 names to Isca names)

    
    # sdor needs to be changed to isca equivalent, but the rest should be ok
    features = ['u_10m', 'v_10m', 'temp_2m', 'sdor', 'rh_2m', 'zsurf']

    isca_to_era5_translator = {
        'u_10m':'u10',
        'v_10m':'v10',
        'temp_2m':'t2m',
        'sdor':'sdor',
        'rh_2m':'rh2',
        'zsurf':'zsurf'
    }

    normalised_features = {}


    # assuming data is a dataset with all the required features at surface (64 x 128) at fixed time (i.e. from atmos_monthly.nc
    # with surface selected)
    # this loops through the features listed and normalises them, and then concatenates them all into big_data which is the input
    # array for prediction
    for f in features:
        if f=='sdor':
            df = sdor_dataset['sdor'][0,...]
        else:
            df = data[f]

        df = reshape(df) # this should be reshaped

        # normalises the feature
        normalised_features[f'{f}'] = normalise_with_training(df, training_min = float(feature_min[isca_to_era5_translator[f'{f}']]), training_max = float(feature_max[isca_to_era5_translator[f'{f}']]))

    features = list(normalised_features.keys())

    big_data = np.append(normalised_features[features[0]], normalised_features[features[1]], axis = 1) # join 0th and 1st feature into big data
    for i in range(2, len(normalised_features)): # join remaining features
        big_data=np.append(big_data, normalised_features[features[i]], axis = 1)
    
    # big data should be normalised and reshaped array with all the input features from isca

    # import model weights
    fileout='SR_resampled_minusd2m4L_16N'
    print('About to read in',fileout,'.json')
    json_file=open(fileout+'.json', 'r')
    loaded_model_json=json_file.read()
    json_file.close()
    model=model_from_json(loaded_model_json)
    # read in the most recent weights saved when training that model.
    model.load_weights(fileout+'.weights.h5')
    # Write out the model architecture of the model you have read in.
    model.summary()

    # use model to predict t2m_std
    pred_t2m_std = model.predict(big_data)

    # unnormalise using training min and max
    unnormalised_pred_t2m_std = unnormalise(pred_t2m_std, training_min = float(feature_min['t2m_std']), training_max = float(feature_max['t2m_std']))

    # reshape to 2d prediction array (which should match isca's grid spacing)
    reshaped_pred_t2m_std = np.reshape(unnormalised_pred_t2m_std, (64, 128))

    # return the 2d array which should match isca's grid
    return reshaped_pred_t2m_std


