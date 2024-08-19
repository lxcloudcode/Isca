
import os

import numpy as np

from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress
import create_std_input_sara as csi

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb= SocratesCodeBase.from_directory(GFDL_BASE)
#cb.enable_fftw3()

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp_name = 'ml_test_with_ml_full_sd_mean_clim_half_month'
exp = Experiment(exp_name, codebase=cb)
exp.clear_rundir()

do_simple_global_value = False


#Tell model how to write diagnostics
diag = DiagTable()
#diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_file('atmos_half_monthly', 15, 'days', time_units='days')
# diag.add_file('atmos_pentad', 5, 'days', time_units='days')
#diag.add_file('atmos_daily', 1, 'days', time_units='days')

# diag.add_file('atmos_timestep', 400, 'seconds', time_units='days')


#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'temp_2m', time_avg=True)
diag.add_field('atmosphere', 'u_10m', time_avg=True)
diag.add_field('atmosphere', 'v_10m', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)
diag.add_field('atmosphere', 'sphum_2m', time_avg=True)
diag.add_field('atmosphere', 'rh_2m', time_avg=True)


diag.add_field('atmosphere', 'pert_t', time_avg=True)
# diag.add_field('atmosphere', 'pert_q', time_avg=True)

diag.add_field('atmosphere', 'bucket_depth', time_avg=True)
#diag.add_field('atmosphere', 'bucket_depth_cond', time_avg=True)
#diag.add_field('atmosphere', 'bucket_depth_conv', time_avg=True)
#diag.add_field('atmosphere', 'bucket_depth_lh', time_avg=True)

diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'albedo', time_avg=True)
diag.add_field('mixed_layer', 'ice_conc', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)

diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'zsurf', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)

diag.add_field('vert_turb', 'z_pbl', time_avg=True)

#radiative tendencies
#diag.add_field('socrates', 'soc_tdt_lw', time_avg=True)
#diag.add_field('socrates', 'soc_tdt_sw', time_avg=True)
#diag.add_field('socrates', 'soc_tdt_rad', time_avg=True)

#net (up) and down surface fluxes
#diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
#diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
#diag.add_field('socrates', 'soc_surf_flux_lw_down', time_avg=True)
#diag.add_field('socrates', 'soc_surf_flux_sw_down', time_avg=True)
#net (up) TOA and downard fluxes
#diag.add_field('socrates', 'soc_olr', time_avg=True)
#diag.add_field('socrates', 'soc_toa_sw', time_avg=True) 
#diag.add_field('socrates', 'soc_toa_sw_down', time_avg=True)

# additional output options commented out 
#diag.add_field('socrates', 'soc_flux_lw', time_avg=True)
#diag.add_field('socrates', 'soc_flux_sw', time_avg=True)
#diag.add_field('socrates', 'soc_co2', time_avg=True)
diag.add_field('socrates', 'soc_ozone', time_avg=True) 
#diag.add_field('socrates', 'soc_coszen', time_avg=True) 
#diag.add_field('socrates', 'soc_spectral_olr', time_avg=True)

exp.diag_table = diag

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 15,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    'socrates_rad_nml': {
        'stellar_constant':1370.,
        'lw_spectral_filename':'INPUT/sp_lw_ga3_1.txt',
        'sw_spectral_filename':'INPUT/sp_sw_ga3_0.txt',
        'do_read_ozone': True,
        'ozone_file_name':'ozone_1990_cmip5',
        'ozone_field_name':'ozone_1990_cmip5',
        'dt_rad':3600,
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':False,
        'co2_ppmv':400.,
        #'solday': 90
    }, 
    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': do_simple_global_value,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,            
        'two_stream_gray': False,     #Use the grey radiation scheme
        'do_socrates_radiation': True,
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use simple Betts miller convection        
        'land_option' : 'input',                     # !Use land mask from input file
        'land_file_name' : 'INPUT/era-spectral7_T42_64x128.out.nc',  # !Tell model where to find input file
        'land_roughness_prefactor' :10.0,           # !How much rougher to make land than ocean
        'roughness_mom'   : 2.e-04,                  # !Ocean roughness lengths  
        'roughness_heat'  : 2.e-04,                  # !Ocean roughness lengths  
       'roughness_moist' : 2.e-04,                  # !Ocean roughness lengths              
        'bucket':True,
        'init_bucket_depth_land':0.15, #Set initial bucket depth over land, default = 20
    },

    'ml_interface_nml': {
          'conv_input_file': 'ml_generated_std',
          'tstd_field_name': 'tstd',
          'qstd_field_name': 'qstd',
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': do_simple_global_value,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': do_simple_global_value,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': do_simple_global_value,
        'old_dtaudv': True,
        'land_evap_prefactor': 0.7,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 200.,
        'prescribe_initial_dist':True,
        'evaporation':True,  
        'depth': 20.,                          #Depth of mixed layer used
        'albedo_value': 0.20,                  #Albedo value used      
        'do_qflux' : False, #Don't use analytic formula for q-fluxes 
        'load_qflux' : True, #Do load q-flux field from an input file
        'time_varying_qflux' : True, #q-flux will be time-varying
        'qflux_file_name' : 'soc_ga3_bucket_jra_55_ice_temps', #Name of q-flux input file
        'do_read_sst' : False, #Read in sst values from input file
        'do_sc_sst' : False, #Do specified ssts (need both to be true)
        # 'sst_file' : 'sst_clim_amip', #Set name of sst input file
        'specify_sst_over_ocean_only' : True, #Make sure sst only specified in regions of ocean.       
        'land_albedo_prefactor': (0.325/0.20),
        'land_option':'input', #             !Tell mixed layer to get land mask from input file
        'land_h_capacity_prefactor': 0.1,    #!What factor to multiply mixed-layer depth by over land. 
        'update_albedo_from_ice':True,
        'ice_albedo_value':0.7, 
        'binary_ice_albedo':False,
        'ice_file_name':'siconc_clim_amip',
        'specify_sst_over_sea_ice':False, 
        # 'linearly_interpolate_sea_ice_temp_and_sst':True,
        # 'ice_sst_file': 'temp_2m_input',
        'update_land_mask_from_ice':False,

    },


    'qe_moist_convection_nml': {
        'rhbm':0.5,
        'Tmin':100.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':do_simple_global_value,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':do_simple_global_value,
        'tcmax':150,
        'tcmin':-200
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,      
    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':40,      #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03,
        'ocean_topog_smoothing': 0.0,
    },

    'spectral_init_cond_nml' : {
        'topog_file_name': 'era-spectral7_T42_64x128.out.nc',
        'topography_option' : 'input',
        'initial_temperature': 240.,
    }

})

#Lets do a run!
if __name__=="__main__":

        inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'), os.path.join(base_dir,'input/era-spectral7_T42_64x128.out.nc'), os.path.join(base_dir,'input/soc_ga3_bucket_jra_55_ice_temps.nc'), os.path.join(base_dir,'input/sp_lw_ga3_1.txt'), os.path.join(base_dir,'input/sp_sw_ga3_0.txt'), os.path.join(base_dir,'input/siconc_clim_amip.nc'), os.path.join(base_dir,'input/ozone_1990_cmip5.nc'), os.path.join(base_dir,'input/ml_generated_std.nc')]
        exp.inputfiles = inputfiles

        cb.compile()
        exp.set_resolution('T42', 40)
        #Set up the experiment object, with the first argument being the experiment name.
        #This will be the name of the folder that the data will appear in.
        PATH_TO_FILE = '/home/links/jc1420/isca_data/ml_test_with_ml_1/restarts/res0120.tar.gz'
        #PATH_TO_FILE = '/home/links/jc1420/isca_data/ml_test_with_ml_full_sd_mean_clim_month/res0127.tar'
        exp.run(121,  num_cores=NCORES, multi_node=False, use_restart = False, restart_file=PATH_TO_FILE, overwrite_data = False)
                # print('about to create new file')
        input_file_name, full_file_path = csi.create_input_file_final(exp_name, 121)
        # print('bingpot')

        exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'), os.path.join(base_dir,'input/era-spectral7_T42_64x128.out.nc'), os.path.join(base_dir,'input/soc_ga3_bucket_jra_55_ice_temps.nc'), os.path.join(base_dir,'input/sp_lw_ga3_1.txt'), os.path.join(base_dir,'input/sp_sw_ga3_0.txt'), os.path.join(base_dir,'input/siconc_clim_amip.nc'), os.path.join(base_dir,'input/ozone_1990_cmip5.nc'), full_file_path]
        exp.update_namelist(
            {'ml_interface_nml':{
                'conv_input_file': input_file_name[:-3]
            },
            'idealized_moist_phys_nml':{
                'read_conv_perturb_input_file':True,
                'perturb_conv_with_ml':True,
            }
            }
        )
        for i in range(122,361):
            exp.run(i, num_cores=NCORES, multi_node=False)
            # print('about to create new file')
            input_file_name, full_file_path = csi.create_input_file_final(exp_name, i)
            # print('bingpot')
            # makes it so the ML aspect of ISCA is now running
            exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'), os.path.join(base_dir,'input/era-spectral7_T42_64x128.out.nc'), os.path.join(base_dir,'input/soc_ga3_bucket_jra_55_ice_temps.nc'), os.path.join(base_dir,'input/sp_lw_ga3_1.txt'), os.path.join(base_dir,'input/sp_sw_ga3_0.txt'), os.path.join(base_dir,'input/siconc_clim_amip.nc'), os.path.join(base_dir,'input/ozone_1990_cmip5.nc'), full_file_path]
            exp.update_namelist(
                {'ml_interface_nml':{
                    'conv_input_file': input_file_name[:-3]
                },
                'idealized_moist_phys_nml':{
                    'read_conv_perturb_input_file':True,
                    'perturb_conv_with_ml':True,
                }
                }
            )
