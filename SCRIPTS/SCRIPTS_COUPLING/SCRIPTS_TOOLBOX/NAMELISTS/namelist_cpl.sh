#-------------------------------------------------------------------------------
# CPL
#-------------------------------------------------------------------------------

# Namelist
#---------
# Note: namelist example files are provided in OASIS_IN/
# if you want to use a pre-built weight file for grid interpolations, point to 
# e.g. namcouple.base.oa.smtho2a 
export namcouplename=namcouple.base.${RUNtype} 

# Coupling frequency
#-------------------
export CPL_FREQ=3600

# Restart files for OASIS
#------------------------
# If TRUE: create OASIS restart files from pre-existing atm/oce/wav outputs.
# If FALSE: create OASIS restart files from calm conditions (need to read at least the grid for each model)
export CPL_restart="FALSE"
export oce_rst_file="${OCE_FILES_DIR}/croco_grd.nc" 
export oce_rst_timeind=-1 # time index (-1 is last) in the file to extract as restart
export atm_rst_file="${ATM_FILES_DIR}/wrfinput_dXX_2005_01_01_00" # the domain dXX will be automatically replaced
export atm_rst_timeind=-1 # time index (-1 is last) in the file to extract as restart
export wav_rst_file="${WAV_FILES_DIR}/ww3.200501.nc" 
export wav_rst_timeind=-1 # time index (-1 is last) in the file to extract as restart

# Weight files for grid interpolations
#-------------------------------------
# If FALSE files for grid interpolations will automatically be generated during the run
# If TRUE use pre-built weight file for grid interpolations
export WEIGHT_FLAG="FALSE" 
# If TRUE : filenames
# They should be placed in OASIS_FILES
# Note: for the moment only works with 1-atm and 1-oce domain
export weight_o2a=""
export weight_a2o=""

