# %%
####################################################################################################
# 2019 
# ----
from pathlib import Path
from echopop.workflows.nwfsc_feat import cli_utils
####################################################################################################
# PARAMETER ENTRY
# ---------------
# ** ENTER FILE INFORMATION FOR ALL INGESTED DATASETS.
# ** ADDITIONAL PARAMETERIZATIONS THROUGHOUT THE SCRIPT SHOULD BE EDITED BASED ON SPECIFIC NEEDS. 
# ** MAKE SURE TO EDIT WITH CARE.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRINT CONSOLE LOGGING MESSAGES
# ---- When set to `True`, logging information will be printed in the terminal/console as the 
# ---- script progresses
try: 
    # ---- FOR CLI USE
    VERBOSE = cli_utils.get_verbose()
except Exception:
    # ---- FOR INTERACTIVE REPL USE
    VERBOSE = True
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA ROOT DIRECTORY
#DATA_ROOT = Path("C:/Data/EchopopData/echopop_2019")
DATA_ROOT = Path("C:/rthomas/Projects/EchoPop_validation/Data")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REPORTS SAVE DIRECTORY
#REPORTS_DIR = DATA_ROOT / "reports"
REPORTS_DIR = Path("C:/rthomas/Projects/EchoPop_validation/Reports")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ALREADY PROCESSED NASC FILE ? 
# ---- When False, the raw NASC exports will be processed. When True, the pre-formatted NASC 
# ---- spreadsheet will be read in. This also requires defining `NASC_EXPORTS_SHEET`
NASC_PREPROCESSED = False
# NASC EXPORTS FILE(S)
Year=2019 
runyearstr=str(Year) #added by RT

NASC_year_path="Exports/"+runyearstr+"/"  
#NASC_EXPORTS_FILES = DATA_ROOT / "raw_nasc/"
NASC_EXPORTS_FILES = DATA_ROOT / NASC_year_path
# NASC EXPORTS SHEET
NASC_EXPORTS_SHEET = "Sheet1"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BIODATA FILE
#BIODATA_FILE = DATA_ROOT / "Biological/1995-2023_biodata_redo.xlsx"
#BIODATA_FILE = DATA_ROOT / "Biological/Updated/1977-2023_Survey_Biodata_old_not_updated.xlsx"
BIODATA_FILE = DATA_ROOT / "Biological/Updated/old master biodata files/1995-2023_biodata_redo.xlsx"
# BIODATA SHEETS
# ---- Assign the sheetnames to 'catch', 'length', 'specimen'
BIODATA_SHEETS = {
    "catch": "biodata_catch",
    "length": "biodata_length",
    "specimen": "biodata_specimen",
}
# BIODATA PROCESSING
# ---- This is used to parse the biodata master spreadsheet, which is required for aligning the 
# ---- biodata with ancillary files such as stratification and transect-haul mappings. This should 
# ---- define the "ships" based on their IDs with the associated survey IDs. If an offset should be 
# ---- added to the haul numbers, that must also be defined here. The target species should also be 
# ---- defined here. 
BIODATA_SHIP_SPECIES = {
    # "ships": {
    #     160: {
    #         "survey": 202106
    #     },
    #     584: {
    #         "survey": 202113,
    #         "haul_offset": 200
    #     }
    # },
    "ships": {
        160: {
            "survey": 201906
        },
        584: {
            "survey": 2019097,
            "haul_offset": 200
        }
    },    
    "species_code": [22500]
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HAUL STRATIFICATION FILE
Strata_year_path="Stratification/"+runyearstr+"/"  
#HAUL_STRATA_FILE = DATA_ROOT / "Stratification/US_CAN strata 2019_final.xlsx"
HAUL_STRATA_FILE = DATA_ROOT / Strata_year_path / "US&CAN strata 2019_EP_08-Jan-2025.xlsx"
# HAUL STRATIFICATION SHEET MAP
# ---- Valid keys are limited to "ks" and "inpfc"
HAUL_STRATA_SHEETS = {
    "inpfc": "INPFC",
    "ks": "Base KS",
}
# GEOGRAPHIC STRATIFICATION FILE
#GEOSTRATA_FILE = DATA_ROOT / "Stratification/Stratification_geographic_Lat_2019_final.xlsx"
GEOSTRATA_FILE = DATA_ROOT / Strata_year_path / "Stratification_geographic_Lat_2019_EP_08-Jan-2025.xlsx"
# GEOGRAPHIC STRATIFICATION SHEET MAP
# ---- Valid keys are limited to "ks" and "inpfc"
GEOSTRATA_SHEETS = {
    "inpfc": "INPFC",
    #"ks": "stratification1",
    "ks": "Base KS",
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KRIGING MESH FILE 
KRIGING_MESH_FILE = (
    #DATA_ROOT / "Kriging_files/Kriging_grid_files/krig_grid2_5nm_cut_centroids_2013.xlsx"
    DATA_ROOT / "Kriging files & parameters/Kriging grid files/krig_grid2_5nm_cut_centroids_2013.xlsx"
)
# KRIGING MESH SHEET
KRIGING_MESH_SHEET = "krigedgrid2_5nm_forChu"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KRIGING AND VARIOGRAM PARAMETERS FILE
KRIGING_VARIOGRAM_PARAMETERS_FILE = (
    #DATA_ROOT / "Kriging_files/default_vario_krig_settings_2019_US_CAN.xlsx"
    DATA_ROOT / "Kriging files & parameters/2019/default_vario_krig_settings_2019_US&CAN.xlsx"
)
# KRIGING AND VARIOGRAM PARAMETERS SHEET
KRIGING_VARIGORAM_PARAMETERS_SHEET = "Sheet1"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 200m ISOBATH FILE
ISOBATH_FILE = (
    #DATA_ROOT / "Kriging_files/Kriging_grid_files/transformation_isobath_coordinates.xlsx"
    DATA_ROOT / "Kriging files & parameters/Kriging grid files/transformation_isobath_coordinates.xlsx"
)
# 200m ISOBATH SHEET
ISOBATH_SHEET = "Smoothing_EasyKrig"

# %%
####################################################################################################
####################################################################################################
# !!! START OF PROCESSING SCRIPT !!
# !!! EDIT CODE BELOW WITH CARE !!
####################################################################################################
####################################################################################################
import logging
import numpy as np
import pandas as pd
from lmfit import Parameters
from echopop.workflows.nwfsc_feat import functions as feat, parameters as feat_parameters, Reporter
import echopop.ingest as ingestion
from echopop import geostatistics, inversion, utils
from echopop.survey import biology, proportions, stratified, transect
from echopop.workflows.nwfsc_feat import apportionment as feat_apportion, biology as feat_biology
####################################################################################################
# FORMAT LOGGER
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
    
logging.basicConfig(
    level=logging.INFO if VERBOSE else logging.WARNING,
    format="%(message)s")


# %%
# ==================================================================================================
# DATA INGESTION 
# ==================================================================================================
# INGEST NASC DATA 
if NASC_PREPROCESSED:
    logging.info(f"Reading pre-generated NASC export file: '{NASC_EXPORTS_FILES.as_posix()}'.")

    # DEFINE COLUMN MAPPING
    FEAT_TO_ECHOPOP_COLUMNS = {
        "transect": "transect_num",
        "region id": "region_id",
        "vessel_log_start": "distance_s",
        "vessel_log_end": "distance_e",
        "spacing": "transect_spacing",
        "layer mean depth": "layer_mean_depth",
        "layer height": "layer_height",
        "bottom depth": "bottom_depth",
        "assigned haul": "haul_num",
    }

    # Read file
    df_nasc = ingestion.nasc.read_nasc_file(
        filename=NASC_EXPORTS_FILES,
        sheetname=NASC_EXPORTS_SHEET,
        column_name_map=FEAT_TO_ECHOPOP_COLUMNS
    )
else:
    logging.info(
        f"Beginning NASC export ingestion for files in: '{NASC_EXPORTS_FILES.as_posix()}'."
    )

    # MERGE EXPORTS
    logging.info(
        "---- Merging NASC exports...\n"
        "     Filename transect pattern: 'T(\\d+)'\n"
        "     Default transect spacing: 10.0 nmi\n"
        "     Default latitude threshold: 60.0 deg."
    )    
    df_intervals, df_exports = ingestion.nasc.merge_echoview_nasc(
        nasc_path = NASC_EXPORTS_FILES,
        filename_transect_pattern = r"T(\d+)",
        default_transect_spacing = 10.0,
        default_latitude_threshold = 60.0,
    )

    # EXPORT REGION NAME MAPPING
    REGION_NAME_EXPR_DICT = {
        "REGION_CLASS": {
            "Age-1 Hake": "^(?:h1a(?![a-z]|m))",
            "Age-1 Hake Mix": "^(?:h1am(?![a-z]|1a))",
            "Hake": "^(?:h(?![a-z]|1a)|hake(?![_]))",
            "Hake Mix": "^(?:hm(?![a-z]|1a)|hake_mix(?![_]))",
        },
        "HAUL_NUM": {
            "[0-9]+",
        },
        "COUNTRY": {
            "CAN": "^[cC]",
            "US": "^[uU]",
        },
    }    
    
    # PROCESS REGION NAMES 
    logging.info(
        "---- Processing export region names\n"
        "     Applying CAN haul number offset: 200"
    )
    df_exports_with_regions = ingestion.nasc.process_region_names(
        df=df_exports,
        region_name_expr_dict=REGION_NAME_EXPR_DICT,
        can_haul_offset=200,
    )

    # GENERATE TRANSECT-REGION-HAUL KEY
    logging.info(
        "---- Generating transect-region-haul key mapping\n"
        "     Searching for the export regions: 'Age-1 Hake', 'Age-1 Hake Mix', 'Hake', 'Hake Mix'"
    )
    df_transect_region_haul_key = ingestion.nasc.generate_transect_region_haul_key(
        df=df_exports_with_regions,
        filter_list=["Age-1 Hake", "Age-1 Hake Mix", "Hake", "Hake Mix"]
    )

    # CONSOLIDATE THE EXPORTS WITH TRANSECT-REGION-HAUL MAPPINGS
    logging.info(
        "---- Finalizing NASC export ingestion\n"
        "     Searching for the export regions: 'Age-1 Hake', 'Age-1 Hake Mix', 'Hake', 'Hake Mix'"
        "     Imputing overlapping region IDs within each interval: True"
    )
    df_nasc = ingestion.nasc.consolidate_echvoiew_nasc(
        df_merged=df_exports_with_regions,
        interval_df=df_intervals,
        region_class_names=["Age-1 Hake", "Age-1", "Age-1 Hake Mix", "Hake", "Hake Mix"],
        impute_region_ids=True,
        transect_region_haul_key_df=df_transect_region_haul_key
    )
logging.info(
    "NASC ingestion complete\n"
    "'df_nasc' created."
)

# %%
BIODATA_FILE

# %%
# ==================================================================================================
# INGEST BIODATA
logging.info(
    f"Beginning biodata ingestion for: '{BIODATA_FILE.as_posix()}'."
)

# BIODATA DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_BIODATA_COLUMNS = {
    "frequency": "length_count",
    "haul": "haul_num",
    "weight_in_haul": "weight",
}

# BIODATA LABEL MAPPING
BIODATA_SEX = {
    "sex": {
        1: "male",
        2: "female",
        3: "unsexed"
    }
}

# READ IN DATA
dict_df_bio = ingestion.load_biological_data(
    biodata_filepath=BIODATA_FILE, 
    #BIODATA_SHEETS=BIODATA_SHEETS, 
    biodata_sheet_map=BIODATA_SHEETS,
    column_name_map=FEAT_TO_ECHOPOP_BIODATA_COLUMNS, 
    subset_dict=BIODATA_SHIP_SPECIES, 
   biodata_label_map=BIODATA_SEX
)
# ---- Remove specimen hauls
feat_biology.remove_specimen_hauls(dict_df_bio)
logging.info(
    "Biodata ingestion complete\n"
    "'dict_df_bio' created."
)

# %%
# ==================================================================================================
# INGEST STRATIFICATION DATA
logging.info(
    "Loading stratification files..."
)

# HAUL-BASED STRATIFICATION DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_STRATA_COLUMNS = {
    "fraction_hake": "nasc_proportion",
    "haul": "haul_num",
    "stratum": "stratum_num",
}

# READ IN STRATA FILE 
logging.info(
    f"Load in haul-based stratification: '{HAUL_STRATA_FILE.as_posix()}'."
)
df_dict_strata = ingestion.load_strata(
    strata_filepath=HAUL_STRATA_FILE, 
    strata_sheet_map=HAUL_STRATA_SHEETS, 
    column_name_map=FEAT_TO_ECHOPOP_STRATA_COLUMNS
)
logging.info(
    "Haul-based stratification loading complete\n"
    "'df_dict_strata' created."
)

# GEOGRAPHIC-BASED STRATIFICATION DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS = {
    "latitude (upper limit)": "northlimit_latitude",
    "stratum": "stratum_num",
}

# READ IN GEOSTRATA FILE
logging.info(
    f"Load in geographic-based stratification: '{GEOSTRATA_FILE.as_posix()}'."
)
df_dict_geostrata = ingestion.load_geostrata(
    geostrata_filepath=GEOSTRATA_FILE, 
    geostrata_sheet_map=GEOSTRATA_SHEETS, 
    column_name_map=FEAT_TO_ECHOPOP_GEOSTRATA_COLUMNS
)
logging.info(
    "Geographic-based stratification loading complete\n"
    "'df_dict_geostrata' created."
)

# %%
# ==================================================================================================
# LOAD KRIGING MESH FILE
logging.info(
    f"Loading kriging mesh file: '{KRIGING_MESH_FILE.as_posix()}'."
)

# KRIGING MESH DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_MESH_COLUMNS = {
    "centroid_latitude": "latitude",
    "centroid_longitude": "longitude",
    "fraction_cell_in_polygon": "fraction",
}

# LOAD MESH
df_mesh = ingestion.load_mesh_data(
    mesh_filepath=KRIGING_MESH_FILE, 
    sheet_name=KRIGING_MESH_SHEET, 
    column_name_map=FEAT_TO_ECHOPOP_MESH_COLUMNS
)
logging.info(
    "Kriging mesh loading complete\n"
    "'df_mesh' created."
)

# %%
# ==================================================================================================
# LOAD ISOBATH FILE
logging.info(
    f"Loading isobath file: '{ISOBATH_FILE}'."
)
df_isobath = ingestion.load_isobath_data(
    isobath_filepath=ISOBATH_FILE,
    sheet_name=ISOBATH_SHEET
)
logging.info(
    f"Iosbath loading complete\n"
    "'df_isobath' created."
)

# %%
# ==================================================================================================
# LOAD KRIGING AND VARIOGRAM PARAMETERS
logging.info(
    f"Loading variogram and kriging parameters: '{KRIGING_VARIOGRAM_PARAMETERS_FILE.as_posix()}'."
)

# PARAMETERS DATAFRAME COLUMN NAME MAPPING
FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS = {
    "hole": "hole_effect_range",
    "lscl": "correlation_range",
    "nugt": "nugget",
    "powr": "decay_power",
    "ratio": "aspect_ratio",
    "res": "lag_resolution",
    "srad": "search_radius",
}

# LOAD IN PARAMETERS
dict_kriging_params, dict_variogram_params = ingestion.load_kriging_variogram_params(
    geostatistic_params_filepath=KRIGING_VARIOGRAM_PARAMETERS_FILE,
    sheet_name=KRIGING_VARIGORAM_PARAMETERS_SHEET,
    column_name_map=FEAT_TO_ECHOPOP_GEOSTATS_PARAMS_COLUMNS
)
logging.info(
    "Variogram and kriging parameter loading complete\n"
    "---- 'dict_variogram_params' created [variogram]\n"
    "---- 'dict_kriging_params' created [kriging]"
)

# %%
# ==================================================================================================
# INITIAL DATA PROCESSING
# ==================================================================================================
# APPLY STRATIFICATION DEFINITIONS TO BIODATA AND NASC
logging.info("Applying strata to datasets...")

# HAUL-BASED STRATA
logging.info(
    "Applying haul-based strata to 'dict_df_bio' and 'df_nasc'.\n"
    "     Default stratum: 0\n"
    "     New columns:\n"
    "         INPFC: 'stratum_inpfc'\n"
    "         KS: 'stratum_ks'"
)
# ---- BIODATA [INPFC]
dict_df_bio = ingestion.join_strata_by_haul(data=dict_df_bio,
                                            strata_df=df_dict_strata["inpfc"],
                                            default_stratum=0,
                                            stratum_name="stratum_inpfc")
# ---- BIODATA [KS]
dict_df_bio = ingestion.join_strata_by_haul(data=dict_df_bio,
                                            strata_df=df_dict_strata["ks"],
                                            default_stratum=0,
                                            stratum_name="stratum_ks")
# ---- NASC [INPFC]
df_nasc = ingestion.join_strata_by_haul(data=df_nasc,
                                        strata_df=df_dict_strata["inpfc"],
                                        default_stratum=0,
                                        stratum_name="stratum_inpfc")
# ---- NASC [KS]
df_nasc = ingestion.join_strata_by_haul(data=df_nasc,
                                        strata_df=df_dict_strata["ks"],
                                        default_stratum=0,
                                        stratum_name="stratum_ks")

# GEOGRAPHIC-BASED STRATA
logging.info(
    "Applying geographic-based strata to 'df_nasc' and 'df_mesh.\n"
    "     New columns:\n"
    "         INPFC: 'geostratum_inpfc'\n"
    "         KS: 'geostratum_ks'"
)
# ---- NASC [INPFC]
df_nasc = ingestion.join_geostrata_by_latitude(data=df_nasc,
                                               geostrata_df=df_dict_geostrata["inpfc"],
                                               stratum_name="geostratum_inpfc")
# ---- NASC [KS]
df_nasc = ingestion.join_geostrata_by_latitude(data=df_nasc,
                                               geostrata_df=df_dict_geostrata["ks"],
                                               stratum_name="geostratum_ks")
# ---- MESH [INPFC]
df_mesh = ingestion.join_geostrata_by_latitude(data=df_mesh, 
                                               geostrata_df=df_dict_geostrata["inpfc"], 
                                               stratum_name="geostratum_inpfc")
# ---- MESH [KS]
df_mesh = ingestion.join_geostrata_by_latitude(data=df_mesh, 
                                               geostrata_df=df_dict_geostrata["ks"], 
                                               stratum_name="geostratum_ks")
logging.info("Strata application complete!")
# ==================================================================================================
# BINIFY DATA 
logging.info(
    "Binning biodata ['dict_df_bio'] into discrete age and length bins\n"
    "     Age bins: [1, 2, 3, ..., 20, 21, 22]\n"
    "     Length bins: [2.0, 4.0, 6.0, ... 76.0, 78.0, 80.0]\n"
    "     New columns:\n"
    "         Age: 'age_bin'\n"
    "         Length: 'length_bin'"
)

# AGE-BINS
AGE_BINS = np.linspace(start=1., stop=22, num=22)
utils.binify(
    data=dict_df_bio, bins=AGE_BINS, bin_column="age",
)

# LENGTH-BINS
LENGTH_BINS = np.linspace(start=2., stop=80., num=40)
utils.binify(
    data=dict_df_bio, bins=LENGTH_BINS, bin_column="length", 
)
logging.info("Age and length binning complete!")

# %%
# ==================================================================================================
# FIT LENGTH-WEIGHT REGRESSION
logging.info("Fitting length-weight regression for each sex and for all fish.")

# CREATE DICTIONARY CONTAINER
dict_length_weight_coefs = {}

# ALL FISH
dict_length_weight_coefs["all"] = dict_df_bio["specimen"].assign(sex="all").groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)

# SEX-SPECIFIC
dict_length_weight_coefs["sex"] = dict_df_bio["specimen"].groupby(["sex"]).apply(
    biology.fit_length_weight_regression,
    include_groups=False
)
logging.info("Fitting length-weight regression for each sex and for all fish complete!")

# %%
# ==================================================================================================
# COMPUTE MEAN WEIGHTS PER LENGTH BIN
logging.info(
    "Computing the mean weight per length bin for each sex and for all fish.\n"
    "     Impute missing length bins using modeled weights: True"
    "     Minimum specimen count per bin: 5"
    )

# SEX-SPECIFIC
df_binned_weights_sex = feat_biology.length_binned_weights(
    data=dict_df_bio["specimen"],
    length_bins=LENGTH_BINS,
    regression_coefficients=dict_length_weight_coefs["sex"],
    impute_bins=True,
    minimum_count_threshold=5
)

# ALL FISH
df_binned_weights_all = feat_biology.length_binned_weights(
    data=dict_df_bio["specimen"].assign(sex="all"),
    length_bins=LENGTH_BINS,
    regression_coefficients=dict_length_weight_coefs["all"],
    impute_bins=True,
    minimum_count_threshold=5,
)

# COMBINE
binned_weight_table = pd.concat([df_binned_weights_sex, df_binned_weights_all], axis=1)
logging.info(
    "Length-binned mean weight calculations complete\n"
    "'binned_weight_table' created."
    )

# %%
# ==================================================================================================
# COMPUTE COUNT DISTRIBUTIONS PER AGE- AND LENGTH-BINS
logging.info(
    "Computing the counts per age- and length-bins across sex.\n"
    "     Stratifying by: 'stratum_ks'"
    "     Grouping by: 'sex'"
    )

# DICTIONARY CONTAINER
dict_df_counts = {}

# AGED
dict_df_counts["aged"] = proportions.compute_binned_counts(
    data=dict_df_bio["specimen"].dropna(subset=["age", "length", "weight"]), 
    groupby_cols=["stratum_ks", "length_bin", "age_bin", "sex"], 
    count_col="length",
    agg_func="size"
)

# UNAGED
dict_df_counts["unaged"] = proportions.compute_binned_counts(
    data=dict_df_bio["length"].copy().dropna(subset=["length"]), 
    groupby_cols=["stratum_ks", "length_bin", "sex"], 
    count_col="length_count",
    agg_func="sum"
)
logging.info(
    "Count distributions across age, length, and sex complete\n"
    "'dict_df_counts' created."
    )

# %%
# ==================================================================================================
# COMPUTE NUMBER PROPORTIONS
logging.info(
    "Computing number proportions across age and length bins\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Excluding: 'sex'='unsexed' from 'dict_df_counts['aged']'"
    )
dict_df_number_proportions = proportions.number_proportions(
    data=dict_df_counts, 
    group_columns=["stratum_ks"],
    exclude_filters={"aged": {"sex": "unsexed"}},
)
logging.info(
    "Number proportions calculation complete\n"
    "'dict_df_number_proportions' created\n"
    )
# ==================================================================================================
# COMPUTE BINNED WEIGHTS
logging.info(
    "Computing the summed weights per age- and length-bins across sex.\n"
    "     Stratifying by: 'stratum_ks'"
    "     Grouping by: 'sex'"
    "     Excluding: 'sex'='unsexed'"
    )

# DICTIONARY CONTAINER
dict_df_weight_distr = {}

# AGED
dict_df_weight_distr["aged"] = proportions.binned_weights(
    length_dataset=dict_df_bio["specimen"],
    include_filter = {"sex": ["female", "male"]},
    interpolate_regression=False,
    contrast_vars="sex",
    table_cols=["stratum_ks", "sex", "age_bin"]
)

# UNAGED
logging.info(
    "Unaged binned weights require additional processing steps.\n"
    "     Interpolating binned length-weight regression estimates: True"
    )
dict_df_weight_distr["unaged"] = proportions.binned_weights(
    length_dataset=dict_df_bio["length"],
    length_weight_dataset=binned_weight_table,
    include_filter = {"sex": ["female", "male"]},
    interpolate_regression=True,
    contrast_vars="sex",
    table_cols=["stratum_ks", "sex"]
)
logging.info(
    "Summed weights per age- and length-bins across sex computation complete\n"
    "'dict_df_weight_distr' created."
    )

# %%
# ==================================================================================================
# COMPUTE WEIGHT PROPORTIONS
logging.info(
    "Computing weight proportions across age and length bins\n"
    "     Stratifying by: 'stratum_ks'"
    "     Grouping by: 'sex'"
    )

# DICTIONARY CONTAINER
dict_df_weight_proportions = {}

# AGED WEIGHT PROPORTIONS
logging.info("Computing aged weight proportions...")
dict_df_weight_proportions["aged"] = proportions.weight_proportions(
    weight_data=dict_df_weight_distr, 
    catch_data=dict_df_bio["catch"], 
    group="aged",
    stratum_col="stratum_ks"
)

# UNAGED SCALING
logging.info("Scaling unaged binned weights...")
standardized_sexed_unaged_weights_df = proportions.scale_weights_by_stratum(
    weights_df=dict_df_weight_distr["unaged"], 
    reference_weights_df=dict_df_bio["catch"].groupby(["stratum_ks"])["weight"].sum(),
    stratum_col="stratum_ks",
)

# UNAGED WEIGHT PROPORTIONS
logging.info(
    "Computing unaged weight proportions\n"
    "     Scaling weight proportions in reference to the aged estimates"
    )
dict_df_weight_proportions["unaged"] = proportions.scale_weight_proportions(
    weight_data=standardized_sexed_unaged_weights_df, 
    reference_weight_proportions=dict_df_weight_proportions["aged"], 
    catch_data=dict_df_bio["catch"], 
    number_proportions=dict_df_number_proportions,
    binned_weights=binned_weight_table["all"],
    group="unaged",
    group_columns = ["sex"],
    stratum_col = "stratum_ks"
)
logging.info(
    "Weight proportions calculation complete\n"
    "'dict_df_weight_proportions' created."
    )

# %%
# ==================================================================================================
# NASC TO BIOMASS CONVERSION
# ==================================================================================================
# INVERSION
logging.info(
    "Beginning inversion based on hake-specific TS-length regression coefficients\n"
    "     Model: 20.0 x log[10](L) + -68.0\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Imputing missing strata: True\n"
    "     Treating hauls as replicates: True"
    )

# DEFINE INVERSION MODEL PARAMETERS
MODEL_PARAMETERS = {
    "ts_length_regression": {
        "slope": 20.,
        "intercept": -68.
    },
    "stratify_by": ["stratum_ks"],
    "expected_strata": df_dict_strata["ks"].stratum_num.unique(),
    "impute_missing_strata": True,
    "haul_replicates": True,
}

# INITIALIZE INVERSION OBJECT
invert_hake = inversion.InversionLengthTS(MODEL_PARAMETERS)
logging.info("Inversion-class object 'invert_hake' created...")

# INVERT NUMBER DENSITY
df_nasc = invert_hake.invert(
    df_nasc=df_nasc, df_length=[dict_df_bio["length"], dict_df_bio["specimen"]]
)
logging.info(
    "Number density inversion complete\n"
    "     New column in 'df_nasc':\n"
    "         Number density (animals nm^-2): 'number_density'"
    )

# %%
# ==================================================================================================
# CONVERT TO BIOMASS
logging.info("Converting number density estimates into biomass estimates")

# SET TRANSECT INTERVAL DISTANCES
logging.info(
    "Defining transect interval distances...\n"
    "     Along-transect interval distance (nmi) threshold: 0.5 nmi"
)
transect.compute_interval_distance(df_nasc=df_nasc, interval_threshold=0.05)

# SET TRANSECT INTERVAL AREAS
logging.info("Defining transect interval areas...")
df_nasc["area_interval"] = (
    df_nasc["transect_spacing"] * df_nasc["distance_interval"]
)

# COMPUTE ABUNDANCE
logging.info(
    "Compute interval abundances...\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Grouping by: 'sex'\n"
    "     Excluding: 'sex'='unsexed' from 'dict_df_number_proportions'"    
)
feat_biology.compute_abundance(
    dataset=df_nasc,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
    exclude_filter={"sex": "unsexed"},
    number_proportions=dict_df_number_proportions
)

# COMPUTE STRATUM-AVERAGED WEIGHTS
df_averaged_weight = proportions.stratum_averaged_weight(
    proportions_dict=dict_df_number_proportions, 
    binned_weight_table=binned_weight_table,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
)

# COMPUTE BIOMASS
logging.info(
    "Compute interval biomass...\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Grouping by: 'sex'\n"  
)
feat_biology.compute_biomass(
    dataset=df_nasc,
    stratify_by=["stratum_ks"],
    group_by=["sex"],
    df_average_weight=df_averaged_weight,
)
logging.info(
    "NASC to biomass conversion complete\n"
    "     New columns in 'df_nasc':\n"
    "         Sex-specific number densities (animals nmi^-2): "
    "'number_density_female'/'number_density_male'\n"
    "         Abundance (animals): 'abundance'/'abundance_female'/'abundance_male'\n"
    "         Biomass density (kg nmi^-2): 'biomass_density'/'biomass_density_female'/"
    "'biomass_density_male'\n"
    "         Biomass (kg): 'biomass'/'biomass_female'/'biomass_male'"
    )

# %%
# ==================================================================================================
# DISTRIBUTE POPULATION ESTIMATES ACROSS AGE AND LENGTH BINS
logging.info(
    "Distribute population estimates across age- and length-bins\n"
    "     Stratifying by: 'stratum_ks'\n"
    "     Grouping by: 'sex'"
)

# ABUNDANCE
logging.info("Distributing abundances...")
dict_transect_abundance_table = feat_apportion.distribute_population_estimates(
    data=df_nasc,
    proportions=dict_df_number_proportions,
    variable="abundance",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
)
logging.info("Abundance distributions complete\n'dict_transect_abundance_table' created.")
# BIOMASS [ALL]
logging.info("Distributing biomass...")
dict_transect_biomass_table = feat_apportion.distribute_population_estimates(
    data=df_nasc,
    proportions=dict_df_weight_proportions,
    variable="biomass",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
)
logging.info("Biomass distribution complete\n'dict_transect_biomass_table' created.")
# BIOMASS [AGED-ONLY]
logging.info("Distributing biomass...\n     Aged-only weight proportions: True")
df_transect_aged_biomass_table = feat_apportion.distribute_population_estimates(
    data=df_nasc,
    proportions=dict_df_weight_proportions["aged"],
    variable="biomass",
    group_by=["sex", "age_bin", "length_bin"],
    stratify_by=["stratum_ks"],
)
logging.info("Aged-biomass distribution complete\n'df_transect_aged_biomass_table' created.")

# %%


# %%
# ==================================================================================================
# GEOSTATISTICS
# ==================================================================================================
# INITIAL
logging.info("Beginning geostatistical analysis...")
# ==================================================================================================
# COORDINATE TRANSFORMATION
logging.info(
    "Transform spatial coordinates for 'df_nasc' and 'df_mesh'\n"
    "     Reference coordinates: 'df_isobath'\n"
    "     Longitudinal offset: -124.78338 deg.E\n"
    "     Latitudinal offset: 45.0 deg.N"
)

# NASC
df_nasc, delta_longitude, delta_latitude = geostatistics.transform_coordinates(
    data = df_nasc,
    reference = df_isobath,
    x_offset = -124.78338,
    y_offset = 45.,   
)

# MESH
df_mesh, _, _ = geostatistics.transform_coordinates(
    data = df_mesh,
    reference = df_isobath,
    x_offset = -124.78338,
    y_offset = 45.,   
    delta_x=delta_longitude,
    delta_y=delta_latitude
)
logging.info(
    "Coordinate transformation complete\n"
    "     New columns:\n"
    "          Transformed longitude: 'x'\n"
    "          Transformed latitude: 'y'\n"
)

# %% [markdown]
# 


