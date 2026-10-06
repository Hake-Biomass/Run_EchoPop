"""
hake.py

Analysis-style hake biomass and abundance workflow for scientist users.

This script keeps the same core calculations as the standard workflow but
organizes execution as a linear analysis pipeline:
1. Load input data
2. Prepare shared stratified tables
3. Run model-specific estimation
4. Run geostatistics and kriging
5. Export reports
"""

# ---------------------------------------------------------
# IMPORTS
# ---------------------------------------------------------
import os
import shutil
from dataclasses import dataclass
from datetime import datetime

import copy
from dotenv import load_dotenv
from lmfit import Parameters
import logging
import numpy as np
import pandas as pd
import xarray as xr

from src.echopop_workflows.workflows.plots import (
    plot_hauls_in_chunks,
    plot_selectivity_sum_to_one_density_by_haul,
)
from src.echopop_workflows.workflows.utils import Configure
from echopop import geostatistics, ingest, inversion, utils
from echopop.reports import Reporter
from echopop.survey import apportionment, biology, proportions, selectivity, stratified, transect

# ---------------------------------------------------------
# CONSTANTS & CONFIGURATION
# ---------------------------------------------------------
load_dotenv()
LOGGING_LEVEL = os.getenv('LOGGING_LEVEL') or 'INFO' # Default to INFO if not set
logger = logging.getLogger(__name__)
logging.basicConfig(level=LOGGING_LEVEL, format="%(message)s")

logger.info("Logging level set to {}".format(LOGGING_LEVEL))
logger.info("Reading .env variables...")

logger.info("Reading configuration yaml files...")
config = Configure()

# ---------------------------------------------------------
# DATA INPUT
# ---------------------------------------------------------

# NASC
def _ingest_nasc():
    logger.info("Starting NASC data ingestion...")
    if config.nasc.is_preprocessed:
        logger.info(f"Reading pre-generated NASC export file: '{config.nasc.file.as_posix()}'.")

        # Read file
        nasc = ingest.nasc.read_nasc_file(
            filename=config.nasc.file,
            sheetname=config.nasc.sheets["nasc"],
            column_name_map=config.nasc.column_mapping,
            haul_uid_config=config.biodata.haul_uid_config,
        )

        logger.info("NASC data read successfully from pre-generated file.")
    else:
        logger.info("Ingesting and preprocessing NASC data...")

        if config.nasc.remove_age1:
            CLASS_REGIONS = ["Hake", "Hake Mix"]
        else:
            CLASS_REGIONS = ["Age-1 Hake", "Age-1 Hake Mix", "Hake", "Hake Mix"]

        # MERGE EXPORTS
        df_intervals, df_exports = ingest.nasc.merge_echoview_nasc(
            file_directory = config.nasc.file,
            filename_transect_pattern = r"T(\d+)",
            default_transect_spacing = 10.0,
            default_latitude_threshold = 60.0,
        )

        ## >>>> Filter out any class "unknown" from df_exports (RT added)
        # In 2015 on x49 there was a region that was named as a hake mix region, but later changed to unknown.
        # Regions in EchoPop are filtered by region name, not by class
        rows = df_exports[df_exports["region_class"] == "unknown"].index
        df_exports.drop(rows, inplace=True)

        # PROCESS REGION NAMES
        df_exports_with_regions = ingest.nasc.process_region_names(
            nasc_cells=df_exports,
            region_name_expr=config.nasc.region_name_expr,
            can_haul_offset=config.biodata.subset["can_haul_offset"],
        )

        # GENERATE TRANSECT-REGION-HAUL KEY
        df_transect_region_haul_key = ingest.nasc.generate_transect_region_haul_key(
            region_data=df_exports_with_regions,
            filter_list=CLASS_REGIONS
        )

        # AGE-1 DOMINATED HAUL REMOVAL
        if config.nasc.remove_age1:
            df_transect_region_haul_key = utils.apply_filters(
                df_transect_region_haul_key, exclude_filter={"haul_num": config.biodata.age1_dominated_hauls}
            )

        # CONSOLIDATE THE EXPORTS WITH TRANSECT-REGION-HAUL MAPPINGS
        nasc = ingest.nasc.consolidate_echvoiew_nasc(
            nasc_data=df_exports_with_regions,
            interval_data=df_intervals,
            region_class_names=CLASS_REGIONS,
            impute_region_ids=True,
            transect_region_haul_key=df_transect_region_haul_key,
            haul_uid_config=config.biodata.haul_uid_config,
        )

        logger.info("NASC data ingested and preprocessed successfully.")


    if config.nasc.drop_transects["start"] and config.nasc.drop_transects["stop"]:
        logger.info(
            f"Dropping transects in range: {config.nasc.drop_transects['start']}-{config.nasc.drop_transects['stop']}"
        )
        nasc = utils.apply_filters(
            nasc,
            exclude_filter={
                "transect_num":
                    np.arange(config.nasc.drop_transects["start"],
                              config.nasc.drop_transects["stop"])
            }
        )

    return nasc

# BIODATA
def _ingest_biodata():
    logger.info("Starting biological data ingestion...")

    if config.biodata.format == "excel":
        biodata = ingest.load_biological_data(
            biodata_filepath=config.biodata.file,
            biodata_sheet_map=config.biodata.sheets,
            column_name_map=config.biodata.column_mapping,
            survey_subset=config.biodata.ship_species,
            biodata_label_map=config.biodata.label_mapping,
            haul_uid_config=config.biodata.haul_uid_config,
        )
    elif config.biodata.format == "db":
        biodata = ingest.load_biodata_db_views(
            db_credentials=config.db_credentials,
            biodata_table_map=config.biodata.sheets,
            subset_dict=config.biodata.ship_species,
            column_name_map=config.biodata.column_mapping
        )


    if "catch" not in biodata:
        raise ValueError("catch not found in biodata sheets. Check biodata file and column mapping.")
    if "length" not in biodata:
        raise ValueError("length not found in biodata sheets. Check biodata file and column mapping.")
    if "specimen" not in biodata:
        raise ValueError("specimen not found in biodata sheets. Check biodata file and column mapping.")


    if config.biodata.age1_dominated_hauls:
        logger.info(f"Removing age-1 dominated hauls: {config.biodata.age1_dominated_hauls}")
        biodata = {
            key: utils.apply_filters(dataset, exclude_filter={"haul_num": config.biodata.age1_dominated_hauls})
            for key, dataset in biodata.items()
        }

    logger.info("Biological data ingested successfully.")
    return biodata

# STRATIFICATION
def _ingest_stratification():
    logger.info("Starting stratification data ingestion...")

    strata = ingest.load_strata(
        strata_filepath=config.strata.file,
        strata_sheet_map=config.strata.sheets,
        column_name_map=config.strata.column_mapping,
        haul_uid_config=config.biodata.haul_uid_config,
    )

    if "inpfc" not in strata:
        raise ValueError("INPFC not found in geostrata sheets. Check strata file and column mapping.")
    if "ks" not in strata:
        raise ValueError("KS not found in geostrata sheets. Check strata file and column mapping.")
    if "stratum_num" not in strata["inpfc"].columns:
        raise ValueError("stratum_num column not found in INPFC strata. Check strata file and column mapping.")
    if "stratum_num" not in strata["ks"].columns:
        raise ValueError("stratum_num column not found in KS strata. Check strata file and column mapping.")

    geostrata = ingest.load_geostrata(
        geostrata_filepath=config.strata.geostrata_file,
        geostrata_sheet_map=config.strata.geostrata_sheets,
        column_name_map=config.strata.geostrata_column_mapping,
    )

    if "inpfc" not in geostrata:
        raise ValueError("INPFC not found in geostrata sheets. Check geostrata file and column mapping.")
    if "ks" not in geostrata:
        raise ValueError("KS not found in geostrata sheets. Check geostrata file and column mapping.")
    if "stratum_num" not in geostrata["inpfc"].columns:
        raise ValueError("stratum_num column not found in INPFC strata. Check geostrata file and column mapping.")
    if "stratum_num" not in geostrata["ks"].columns:
        raise ValueError("stratum_num column not found in KS strata. Check geostrata file and column mapping.")

    logger.info("Stratification data ingested successfully.")
    return strata, geostrata

# Kriging
def _ingest_kriging():
    logger.info("Starting kriging data ingestion...")

    mesh = ingest.load_mesh_data(
        mesh_filepath=config.kriging.mesh_file,
        sheet_name=config.kriging.mesh_sheets["mesh"],
        column_name_map=config.kriging.mesh_column_mapping
    )

    isobath = ingest.load_isobath_data(
        isobath_filepath=config.kriging.isobath_file,
        sheet_name=config.kriging.isobath_sheets["isobath"],
    )

    kriging_params, variogram_params = ingest.load_kriging_variogram_params(
        geostatistic_params_filepath=config.kriging.variogram_file,
        sheet_name=config.kriging.variogram_sheets["variogram_parameters"],
        column_name_map=config.kriging.variogram_column_mapping
    )

    logger.info("Kriging data ingested successfully.")
    return mesh, isobath, kriging_params, variogram_params

# ---------------------------------------------------------
# DATA PROCESSING
# ---------------------------------------------------------
def _apply_stratification(dict_df_bio, df_nasc, df_mesh, df_dict_strata, df_dict_geostrata):
    """Apply haul-based and geographic strata to biodata, nasc and mesh frames."""
    # haul-based
    dict_df_bio = ingest.join_strata_by_uid(data=dict_df_bio, strata=df_dict_strata["inpfc"],
                                               default_stratum=0, stratum_name="stratum_inpfc")
    dict_df_bio = ingest.join_strata_by_uid(data=dict_df_bio, strata=df_dict_strata["ks"],
                                               default_stratum=0, stratum_name="stratum_ks")
    df_nasc = ingest.join_strata_by_uid(data=df_nasc, strata=df_dict_strata["inpfc"],
                                           default_stratum=0, stratum_name="stratum_inpfc")
    df_nasc = ingest.join_strata_by_uid(data=df_nasc, strata=df_dict_strata["ks"],
                                           default_stratum=0, stratum_name="stratum_ks")
    # geographic-based
    df_dict_geostrata["inpfc"].drop_duplicates("northlimit_latitude", inplace=True)
    df_nasc = ingest.join_geostrata_by_latitude(data=df_nasc, geostrata=df_dict_geostrata["inpfc"],
                                                   stratum_name="geostratum_inpfc")
    df_nasc = ingest.join_geostrata_by_latitude(data=df_nasc, geostrata=df_dict_geostrata["ks"],
                                                   stratum_name="geostratum_ks")
    df_mesh = ingest.join_geostrata_by_latitude(data=df_mesh, geostrata=df_dict_geostrata["inpfc"],
                                                   stratum_name="geostratum_inpfc")
    df_mesh = ingest.join_geostrata_by_latitude(data=df_mesh, geostrata=df_dict_geostrata["ks"],
                                                   stratum_name="geostratum_ks")
    logger.info("Strata applied to biodata, nasc and mesh.")
    return dict_df_bio, df_nasc, df_mesh


def _binify_biodata(dict_df_bio, age_bins=None, length_bins=None):
    """Binify ages and lengths in biodata dict in-place and return the dict."""
    age_bins = age_bins if age_bins is not None else np.linspace(start=1., stop=22, num=22)
    length_bins = length_bins if length_bins is not None else np.linspace(start=2., stop=80., num=40)
    utils.binify(data=dict_df_bio, bins=age_bins, bin_column="age")
    utils.binify(data=dict_df_bio, bins=length_bins, bin_column="length")
    logger.info("Binning complete.")
    return dict_df_bio, age_bins, length_bins


def _fit_length_weight_regressions(dict_df_bio):
    """Fit length-weight regressions and return dictionary of coefficients for 'sex' and 'all'."""
    coefs = {}
    coefs["all"] = dict_df_bio["specimen"].assign(sex="all").groupby(["sex"]).apply(
        biology.fit_length_weight_regression, include_groups=False
    )
    coefs["sex"] = dict_df_bio["specimen"].groupby(["sex"]).apply(
        biology.fit_length_weight_regression, include_groups=False
    )
    logger.info("Length-weight regressions fitted.")
    return coefs


def _compute_binned_weights(dict_df_bio, length_bins, length_weight_coefs):
    """Compute mean weight per length bin for sexes and 'all' and return concatenated xarray DataArray."""
    da_binned_weights_sex = biology.length_binned_weights(
        data=dict_df_bio["specimen"],
        length_bins=length_bins,
        regression_coefficients=length_weight_coefs["sex"],
        impute_bins=True,
        minimum_count_threshold=5,
    )
    da_binned_weights_all = biology.length_binned_weights(
        data=dict_df_bio["specimen"].assign(sex="all"),
        length_bins=length_bins,
        regression_coefficients=length_weight_coefs["all"],
        impute_bins=True,
        minimum_count_threshold=5,
    )
    da_binned_weight_table = xr.concat([da_binned_weights_sex, da_binned_weights_all], dim="sex")
    logger.info("Binned weights computed.")
    return da_binned_weight_table


def _apply_net_selectivity(dict_df_bio, da_binned_weight_table):
    """Apply net selectivity expansion factors to specimen-level data."""
    specimen_data = dict_df_bio["specimen"].merge(
        dict_df_bio["catch"][["uid", "net_num"]].drop_duplicates("uid"),
        on="uid",
        how="left",
    )

    specimen_data_selectivity = selectivity.assign_selectivity_expansion(
        specimen_data,
        config.selectivity.parameters,
        net_column="net_num",
    )
    if config.selectivity.plot_during_apply:
        plot_hauls_in_chunks(
            specimen_data_selectivity,
            save_pdf_path=None,
            skip_show=config.selectivity.skip_show_plots,
            logger=logger,
        )
    logger.info("Net selectivity expansion factors applied.")
    return specimen_data_selectivity


def _compute_counts_and_proportions(dict_df_bio, aged_data=None, aged_count_col="length", aged_agg_func="size"):
    """Compute aged/unaged counts and number proportions stratified by stratum name."""
    aged_data = (
        aged_data
        if aged_data is not None
        else dict_df_bio["specimen"].dropna(subset=["age", "length", "weight"])
    )
    ds_counts = xr.Dataset()
    ds_counts["aged"] = proportions.compute_binned_counts(
        data=aged_data,
        groupby_cols=[config.strata.stratum_name, "length_bin", "age_bin", "sex"],
        count_col=aged_count_col,
        agg_func=aged_agg_func,
    )
    ds_counts["unaged"] = proportions.compute_binned_counts(
        data=dict_df_bio["length"].copy().dropna(subset=["length"]),
        groupby_cols=[config.strata.stratum_name, "length_bin", "sex"],
        count_col="length_count",
        agg_func="sum",
    )
    dict_ds_number_proportion = proportions.number_proportions(
        data=ds_counts,
        stratum_dim=config.strata.stratum_name,
        exclude_filters={"aged": {"sex": "unsexed"}},
    )
    logger.info("Counts and number proportions computed.")
    return dict_ds_number_proportion


def _compute_binned_weight_distributions(dict_df_bio, da_binned_weight_table, aged_length_data=None):
    """Compute summed binned weights (aged/unaged) and return weight proportions dict."""
    aged_length_data = aged_length_data if aged_length_data is not None else dict_df_bio["specimen"]
    ds_da_weight_dist = xr.Dataset()
    ds_da_weight_dist["aged"] = proportions.binned_weights(
        length_data=aged_length_data,
        include_filter={"sex": ["male", "female"]},
        interpolate_regression=False,
        group_columns=[config.strata.stratum_name, "sex", "age_bin"],
    )
    ds_da_weight_dist["unaged"] = proportions.binned_weights(
        length_data=dict_df_bio["length"],
        include_filter={"sex": ["male", "female"]},
        interpolate_regression=True,
        length_weight_data=da_binned_weight_table,
        group_columns=[config.strata.stratum_name, "sex"],
    )
    logger.info("Binned weight distributions computed.")
    return ds_da_weight_dist

def _compute_binned_weight_proportions(dict_df_bio, da_binned_weight_table, dict_ds_number_proportion, ds_da_weight_dist):
    """Compute weight proportions for aged and unaged groups and return as dict."""
    dict_da_weight_proportion = {}
    dict_da_weight_proportion["aged"] = proportions.weight_proportions(
        weight_data=ds_da_weight_dist["aged"],
        catch_data=dict_df_bio["catch"],
        stratum_dim=config.strata.stratum_name,
    )
    dict_da_weight_proportion["unaged"] = proportions.fitted_weight_proportions(
        weight_data=ds_da_weight_dist["unaged"],
        aged_weight_proportions=dict_da_weight_proportion["aged"],
        number_proportions=dict_ds_number_proportion["unaged"],
        binned_weights=da_binned_weight_table.sel(sex="all"),
        stratum_dim=config.strata.stratum_name,
    )
    logger.info("Binned weight proportions computed.")
    return dict_da_weight_proportion


def _compute_model_proportions(
    dict_df_bio,
    da_binned_weight_table,
    aged_data=None,
    aged_count_col="length",
    aged_agg_func="size",
    weight_distributions=None,
):
    """Build number and weight proportions with a configurable aged-data source."""
    number_proportions = _compute_counts_and_proportions(
        dict_df_bio,
        aged_data=aged_data,
        aged_count_col=aged_count_col,
        aged_agg_func=aged_agg_func,
    )
    if weight_distributions is None:
        weight_distributions = _compute_binned_weight_distributions(
            dict_df_bio,
            da_binned_weight_table,
            aged_length_data=aged_data,
        )
    weight_proportions = _compute_binned_weight_proportions(
        dict_df_bio,
        da_binned_weight_table,
        number_proportions,
        weight_distributions,
    )
    return number_proportions, weight_proportions, weight_distributions


def _setup_inversion_and_invert(df_nasc, dict_df_bio, df_dict_strata):
    """Configure inversion object, run inversion and return (invert_hake, df_nasc_inverted)."""

    model_parameters = config.other.inversion_parameters
    model_parameters["expected_strata"] = df_dict_strata["ks"].stratum_num.unique()

    invert_hake = inversion.InversionLengthTS(model_parameters)
    df_nasc_inverted = invert_hake.invert(nasc_data=df_nasc, length_data=[dict_df_bio["length"], dict_df_bio["specimen"]])
    logger.info("Inversion complete.")
    return invert_hake, df_nasc_inverted

def _compute_transect_density_estimates(df_nasc, number_proportions, da_binned_weight_table):
    """Compute transect-level area, abundance, and biomass from proportions."""
    logger.info("Converting number density estimates into biomass estimates")
    transect.compute_interval_distance(
        nasc_data=df_nasc, interval_threshold=config.nasc.transect_interval_distance_threshold
    )
    df_nasc["area_interval"] = df_nasc["transect_spacing"] * df_nasc["distance_interval"]

    biology.compute_abundance(
        transect_data=df_nasc,
        exclude_filter={"sex": "unsexed"},
        number_proportions=number_proportions,
    )
    da_averaged_weight = proportions.stratum_averaged_weight(
        number_proportions=number_proportions,
        length_weight_data=da_binned_weight_table,
        stratum_dim=config.strata.stratum_name,
    )
    biology.compute_biomass(transect_data=df_nasc, stratum_weights=da_averaged_weight)
    return da_averaged_weight


def _remove_age1_contribution(
    df_nasc,
    number_proportions,
    weight_proportions,
    age_source_key=None,
):
    """Optionally remove age-1 contribution from transect-level estimates."""
    if not config.nasc.remove_age1:
        return df_nasc.copy()

    number_slice_source = (
        number_proportions[age_source_key] if age_source_key is not None else number_proportions
    )
    weight_slice_source = (
        weight_proportions[age_source_key] if age_source_key is not None else weight_proportions
    )
    age1_filter = {"age_bin": [1]}
    age1_nasc_proportions = proportions.get_nasc_proportions_slice(
        number_proportions=number_slice_source,
        stratum_dim=config.strata.stratum_name,
        ts_length_regression_parameters=config.other.inversion_parameters["ts_length_regression"],
        include_filter=age1_filter,
    )
    age1_number_proportions = proportions.get_number_proportions_slice(
        number_proportions=number_slice_source,
        stratum_dim=config.strata.stratum_name,
        include_filter=age1_filter,
    )
    age1_weight_proportions = proportions.get_weight_proportions_slice(
        weight_proportions=weight_slice_source,
        stratum_dim=config.strata.stratum_name,
        include_filter=age1_filter,
        number_proportions=number_proportions,
        length_threshold_min=config.biodata.age1_length_threshold,
        weight_proportion_threshold=config.biodata.age1_weight_proportion_threshold,
    )
    df_nasc_proc = apportionment.remove_group_from_estimates(
        transect_data=df_nasc,
        group_proportions=xr.Dataset(
            {
                "nasc": age1_nasc_proportions,
                "abundance": age1_number_proportions,
                "biomass": age1_weight_proportions,
            }
        ),
    )
    logger.info("Age-1 contribution removal complete.")
    return df_nasc_proc


def _log_transect_summary(df_nasc_proc):
    """Log key transect-level totals used by scientists to verify model outputs."""
    logger.info("NASC to biomass conversion complete.")
    logger.info(
        f"----------------------\n"
        f"Transect-based results\n"
        f"     Total NASC: {df_nasc_proc['nasc'].sum():.1f} m²nmi⁻²\n"
        f"     Mean number density: {df_nasc_proc['number_density'].mean():.1f} animals nmi⁻²\n"
        f"     Total abundance: {df_nasc_proc['abundance'].sum():.0f} fish\n"
        f"     Mean biomass density: {df_nasc_proc['biomass_density'].mean() * 1e-6:.3f} kmt nmi⁻²\n"
        f"     Total biomass: {df_nasc_proc['biomass'].sum() * 1e-6:.1f} kmt"
    )

def _compute_biomass_and_abundance_estimates(
    df_nasc, dict_ds_number_proportion, da_binned_weight_table, dict_da_weight_proportion
):
    """Compute transect biomass/abundance and distributed tables."""
    da_averaged_weight = _compute_transect_density_estimates(
        df_nasc, dict_ds_number_proportion, da_binned_weight_table
    )
    df_nasc_proc = _remove_age1_contribution(
        df_nasc=df_nasc,
        number_proportions=dict_ds_number_proportion,
        weight_proportions=dict_da_weight_proportion,
        age_source_key="aged",
    )
    _log_transect_summary(df_nasc_proc)

    logger.info("Distributing abundances...")
    dict_ds_transect_abundance_table = apportionment.distribute_population_estimates(
        data=df_nasc,
        proportions=dict_ds_number_proportion,
        variable="abundance",
        group_columns=["sex", "age_bin", "length_bin", config.strata.stratum_name],
    )

    logger.info("Distributing biomass...")
    dict_ds_transect_biomass_table = apportionment.distribute_population_estimates(
        data=df_nasc,
        proportions=dict_da_weight_proportion,
        variable="biomass",
        group_columns=["sex", "age_bin", "length_bin", config.strata.stratum_name],
    )
    dict_ds_transect_biomass_table["standardized_unaged"] = (
        apportionment.distribute_unaged_from_aged(
            population_table=dict_ds_transect_biomass_table["unaged"],
            reference_table=dict_ds_transect_biomass_table["aged"],
            stratum_dim=config.strata.stratum_name,
            impute=True,
            impute_variable=['age_bin'],
        )
    )
    da_transect_biomass_table = apportionment.sum_population_tables(
        population_tables={
            "aged": dict_ds_transect_biomass_table["aged"],
            "unaged": dict_ds_transect_biomass_table["standardized_unaged"],
        },
    )
    logger.info("Biomass and abundance estimates computed and distributed to transect level.")
    return (
        df_nasc_proc,
        dict_ds_transect_abundance_table,
        da_transect_biomass_table,
        da_averaged_weight,
    )

def _geostatistical_analysis(df_nasc_proc, df_mesh, df_isobath, dict_variogram_params):
    """Perform variogram analysis and fit variogram model to NASC-based biomass density estimates."""
    # COORDINATE TRANSFORMATION
    # NASC
    df_nasc_proc, delta_longitude, delta_latitude = geostatistics.transform_coordinates(
        data=df_nasc_proc,
        reference=df_isobath,
        x_offset=config.kriging.geostat_offset_x,
        y_offset=config.kriging.geostat_offset_y,
    )

    # MESH
    df_mesh, _, _ = geostatistics.transform_coordinates(
        data=df_mesh,
        reference=df_isobath,
        x_offset=config.kriging.geostat_offset_x,
        y_offset=config.kriging.geostat_offset_y,
        delta_x=delta_longitude,
        delta_y=delta_latitude
    )

    logger.info("Coordinates transformed for geostatistical analysis.")


    if config.kriging.variogram["optimize"]:
        # INITIALIZE VARIOGRAM-CLASS OBJECT
        vgm = geostatistics.Variogram(
            lag_resolution=config.kriging.variogram["lag_resolution"],
            n_lags=config.kriging.variogram["n_lags"],
            coordinate_names=("x", "y"),
        )

        # EMPIRICAL VARIOGRAM
        vgm.calculate_empirical_variogram(
            data=df_nasc_proc,
            variable="biomass_density",
            azimuth_filter=True,
            azimuth_angle_threshold=config.kriging.variogram["azimuth_angle_threshold"],
        )

        # SET UP FITTING PARAMETERS
        # ----- lmfit.Parameters tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
        logger.info(
            f"Optimizing variogram parameters using non-linear least-squares\n"
            f"     Model: Exponential-Bessel (['exponential', 'bessel'])\n"
            f"     Initial values:\n"
            f"          Nugget: {dict_variogram_params['nugget']}\n"
            f"          Sill: {dict_variogram_params['sill']}\n"
            f"          Correlation range: {dict_variogram_params['correlation_range']}\n"
            f"          Hole effect range: {dict_variogram_params['hole_effect_range']}\n"
            f"          Decay power exponent: {dict_variogram_params['decay_power']}"
        )
        variogram_parameters_lmfit = Parameters()
        variogram_parameters_lmfit.add_many(
            ("nugget", dict_variogram_params["nugget"], True, 0.),
            ("sill", dict_variogram_params["sill"], True, 0.),
            ("correlation_range", dict_variogram_params["correlation_range"], True, 0.),
            ("hole_effect_range", dict_variogram_params["hole_effect_range"], True, 0.),
            ("decay_power", dict_variogram_params["decay_power"], True, 1.25, 1.75),
        )

        # OPTIMIZATION PARAMETERS
        OPTIM_ARGS = {
            "max_nfev": None, "ftol": 1e-08, "gtol": 1e-8, "xtol": 1e-8, "diff_step": None,
            "tr_solver": "exact", "x_scale": 1., "jac": "2-point"
        }
        logger.info(
            f"Optimization arguments:\n"
            f"{OPTIM_ARGS}"
        )

        # RUN MINIMIZER
        best_fit_parameters = vgm.fit_variogram_model(
            model=["exponential", "bessel"],
            model_parameters=variogram_parameters_lmfit,
            optimizer_kwargs=OPTIM_ARGS,
        )
        logger.info(
            f"Variogram parameter fitting complete\n"
            f"     Best-fit parameters:\n"
            f"     {best_fit_parameters}"
        )
    else:
        best_fit_parameters = {
            "nugget": dict_variogram_params["nugget"],
            "sill": dict_variogram_params["sill"],
            "hole_effect_range": dict_variogram_params["hole_effect_range"],
            "correlation_range": dict_variogram_params["correlation_range"],
            "decay_power": dict_variogram_params["decay_power"]
        }

    logger.info("Variogram analysis complete.")
    return best_fit_parameters, df_nasc_proc, df_mesh

def _kriging_analysis(best_fit_parameters, df_nasc_proc, df_mesh, dict_da_weight_proportion, da_averaged_weight, invert_hake):
    """Perform ordinary kriging to interpolate biomass density estimates across the full mesh and convert to NASC."""
    def _get_transect_mesh_region_map():
        configured_year = str(getattr(config, "year", "2019"))
        mapping_name = f"transect_mesh_region_{configured_year}"
        mapping_fn = getattr(utils.feat_parameters, mapping_name, None)
        if mapping_fn is None:
            fallback_name = "transect_mesh_region_2019"
            mapping_fn = getattr(utils.feat_parameters, fallback_name)
            logger.warning(
                "No FEAT transect-region map '%s'; falling back to '%s'.",
                mapping_name,
                fallback_name,
            )
        return mapping_fn

    def _can_apply_feat_crop(transects: pd.DataFrame, mapping_fn) -> tuple[bool, str]:
        transect_numbers = set(transects["transect_num"].unique())

        for region in (1, 2, 3):
            region_start, region_end, _, _ = mapping_fn(region)
            has_region_data = any(
                region_start <= transect_num <= region_end
                for transect_num in transect_numbers
            )
            if not has_region_data:
                return False, (
                    f"no transects available for FEAT region {region} "
                    f"({region_start}-{region_end})"
                )

            if region == 3:
                missing_endpoints = [
                    transect_num
                    for transect_num in (region_start, region_end)
                    if transect_num not in transect_numbers
                ]
                if missing_endpoints:
                    missing_label = ", ".join(str(value) for value in missing_endpoints)
                    return False, (
                        "missing region-3 endpoint transects required by FEAT crop: "
                        f"{missing_label}"
                    )

        return True, ""

    KRIGING_PARAMETERS = {
        "search_radius": best_fit_parameters["correlation_range"] * config.kriging.search_radius_multiplier,
        "aspect_ratio": config.kriging.aspect_ratio,
        "k_min": config.kriging.k_min,
        "k_max": config.kriging.k_max,
    }

    # VARIOGRAM PARAMETERS CONTAINER
    VARIOGRAM_PARAMETERS = {
        "model": ["exponential", "bessel"],
        **best_fit_parameters
    }

    # INITIALIZE CLASS OBJECT
    krg = geostatistics.Kriging(
        mesh=df_mesh,
        kriging_params=KRIGING_PARAMETERS,
        variogram_params=VARIOGRAM_PARAMETERS,
        coordinate_names=("x", "y"),
    )

    # REGISTER KRIGING METHOD
    krg.register_search_strategy("FEAT_strategy", utils.feat_functions.western_boundary_search_strategy)
    # ---- Parameterize
    transect_western_extents = utils.feat_functions.get_survey_western_extents(
        transects=df_nasc_proc,
        coordinate_names=("x", "y"),
        latitude_threshold=config.kriging.feat_latitude_threshold,
    )
    FEAT_STRATEGY_KWARGS = {
        "western_extent": transect_western_extents,
    }

    if not config.kriging.extrapolate:
        transect_mesh_region_map = _get_transect_mesh_region_map()
        can_use_feat_crop, reason = _can_apply_feat_crop(df_nasc_proc, transect_mesh_region_map)
        if can_use_feat_crop:
            krg.crop_mesh(
                crop_function=utils.feat_functions.transect_ends_crop,
                transects=df_nasc_proc,
                latitude_resolution=1.25 / 60.0,
                transect_mesh_region_function=transect_mesh_region_map,
            )
        else:
            logger.warning(
                "Skipping FEAT transect-ends crop (%s); using default hull crop instead.",
                reason,
            )
            krg.crop_mesh(transects=df_nasc_proc)

        krg.register_search_strategy("FEAT_strategy", utils.feat_functions.western_boundary_search_strategy)

    # RUN KRIGING
    df_kriged_results = krg.krige(
        transects=df_nasc_proc,
        variable="biomass_density",
        extrapolate=config.kriging.extrapolate,
        default_mesh_cell_area=config.kriging.default_mesh_cell_area,
        adaptive_search_strategy="FEAT_strategy",
        custom_search_kwargs=FEAT_STRATEGY_KWARGS
    )
    logger.info("Ordinary kriging complete.")


    # CONVERT TO BIOMASS
    df_kriged_results["biomass"] = df_kriged_results["biomass_density"] * df_kriged_results["area"]

    # BIOMASS TO NASC
    apportionment.mesh_biomass_to_nasc(
        mesh_data=df_kriged_results,
        biodata=dict_da_weight_proportion,
        group_columns=["sex", config.strata.stratum_name],
        mesh_biodata_link={config.strata.geostratum_name: config.strata.stratum_name},
        stratum_weights=da_averaged_weight.sel(sex="all"),
        stratum_sigma_bs=invert_hake.sigma_bs_strata,
    )

    # SUMMARIZE KRIGING RESULTS
    logger.info(
        f"----------------------\n"
        f"Kriging-based results\n"
        f"     Total derived NASC: {df_kriged_results['nasc'].sum():.1f} m²nmi⁻²\n"
        f"     Total derived abundance: {df_kriged_results['abundance'].sum():.0f} fish\n"
        f"     Total biomass: {df_kriged_results['biomass'].sum() * 1e-6:.1f} kmt"
    )

    logger.info("Biomass to NASC conversion complete.")
    return df_kriged_results

def _distribute_population_estimates(df_kriged_results, dict_ds_number_proportion, dict_da_weight_proportion):
    """Distribute kriged population estimates to age-length bins and return distributed abundance and biomass tables."""
    dict_ds_kriged_abundance_table = apportionment.distribute_population_estimates(
        data=df_kriged_results,
        proportions=dict_ds_number_proportion,
        variable="abundance",
        group_columns=["sex", "age_bin", "length_bin", config.strata.stratum_name],
        data_proportions_link={config.strata.geostratum_name: config.strata.stratum_name}
    )

    # SCALE UNAGED ABUNDANCE
    dict_ds_kriged_abundance_table[
        "standardized_unaged"] = apportionment.distribute_unaged_from_aged(
        population_table=dict_ds_kriged_abundance_table["unaged"],
        reference_table=dict_ds_kriged_abundance_table["aged"],
        stratum_dim=config.strata.stratum_name,
        impute=True,
        impute_variable=['age_bin'],
    )

    # BIOMASS [ALL]
    dict_ds_kriged_biomass_table = apportionment.distribute_population_estimates(
        data=df_kriged_results,
        proportions=dict_da_weight_proportion,
        variable="biomass",
        group_columns=["sex", "age_bin", "length_bin", config.strata.stratum_name],
        data_proportions_link={config.strata.geostratum_name: config.strata.stratum_name}
    )

    # SCALE UNAGED BIOMASS
    dict_ds_kriged_biomass_table["standardized_unaged"] = apportionment.distribute_unaged_from_aged(
        population_table=dict_ds_kriged_biomass_table["unaged"],
        reference_table=dict_ds_kriged_biomass_table["aged"],
        stratum_dim=config.strata.stratum_name,
        impute=True,
        impute_variable=["age_bin"],
    )

    # CONSOLIDATE
    # ---- ABUNDANCE
    logger.info("Consolidating abundance tables...")
    da_kriged_abundance_table = apportionment.sum_population_tables(
        population_tables={
            "aged": dict_ds_kriged_abundance_table["aged"],
            "unaged": dict_ds_kriged_abundance_table["standardized_unaged"]
        },
    )
    # ---- Biomass
    logger.info("Consolidating biomass tables...")
    da_kriged_biomass_table = apportionment.sum_population_tables(
        population_tables={
            "aged": dict_ds_kriged_biomass_table["aged"],
            "unaged": dict_ds_kriged_biomass_table["standardized_unaged"]
        },
    )

    if config.nasc.remove_age1:
        da_kriged_abundance_table = apportionment.reallocate_excluded_estimates(
            population_table=da_kriged_abundance_table,
            exclusion_filter={"age_bin": [1]},
            group_columns=["sex"],
        )
        da_kriged_biomass_table = apportionment.reallocate_excluded_estimates(
            population_table=da_kriged_biomass_table,
            exclusion_filter={"age_bin": [1]},
            group_columns=["sex"],
        )
        da_kriged_abundance_table_aged = apportionment.reallocate_excluded_estimates(
            population_table=dict_ds_kriged_abundance_table["aged"],
            exclusion_filter={"age_bin": [1]},
            group_columns=["sex"],
        )
        # ---- Construct dictionary of tables
        dict_ds_kriged_abundance_table = {
            "aged": da_kriged_abundance_table_aged,
            "unaged": dict_ds_kriged_abundance_table["unaged"],
        }
    else:
        da_kriged_abundance_table = da_kriged_abundance_table
        da_kriged_biomass_table = da_kriged_biomass_table
        dict_ds_kriged_abundance_table = copy.deepcopy(dict_ds_kriged_abundance_table)

    logger.info("Population estimate distribution complete.")
    return da_kriged_abundance_table, da_kriged_biomass_table, dict_ds_kriged_abundance_table

def _jolly_and_hampton_analysis(df_nasc_proc, df_kriged_results, df_dict_geostrata):
    "Stratified analysis to estimate uncertainties (Jolly and Hampton, 1990)"

    # INITIALIZE JOLLYHAMPTON CLASS OBJECT
    jh = stratified.JollyHampton(config.kriging.jh_params)

    jh.stratified_bootstrap(data=df_nasc_proc,
                            stratum_dim="geostratum_inpfc",
                            variable="biomass")

    df_jh_transect_results = jh.summarize(ci_percentile=0.95, ci_method="t-jackknife")

    # REPORT
    logger.info(
        f"Mean transect CV [95% CI]: "
        f"{df_jh_transect_results.loc[('survey', 'cv')]['mean']:.3f} "
        f"[{df_jh_transect_results.loc[('survey', 'cv')]['low']:.3f}, "
        f"{df_jh_transect_results.loc[('survey', 'cv')]['high']:.3f}]\n"
        f"    Resampling/bootstrapping bias: "
        f"{df_jh_transect_results.loc[('survey', 'cv')]['bias']:.3f}"
    )

    # RUN ON KRIGED DATA
    # ---- Create virtual transects
    kriged_transects = jh.create_virtual_transects(
        mesh_data=df_kriged_results,
        geostrata=df_dict_geostrata["inpfc"],
        stratum_dim="geostratum_inpfc",
        variable="biomass",
    )
    # ---- Run rest of flow

    jh.stratified_bootstrap(data=kriged_transects,
                            stratum_dim="geostratum_inpfc",
                            variable="biomass")

    df_jh_kriged_results = jh.summarize(ci_percentile=0.95, ci_method="t-jackknife")

    # REPORT
    logger.info(
        f"Mean kriging CV [95% CI]: "
        f"{df_jh_kriged_results.loc[('survey', 'cv')]['mean']:.3f} "
        f"[{df_jh_kriged_results.loc[('survey', 'cv')]['low']:.3f}, "
        f"{df_jh_kriged_results.loc[('survey', 'cv')]['high']:.3f}]\n"
        f"    Resampling/bootstrapping bias: "
        f"{df_jh_kriged_results.loc[('survey', 'cv')]['bias']:.3f}"
    )

    return df_jh_transect_results, df_jh_kriged_results

# ---------------------------------------------------------
# DATA EXPORT/REPORTING
# ---------------------------------------------------------

def generate_summary_report(nasc_proccessed, kriged_results, target_dir, transect_length_age_biomass=None):
    transect_biomass_text = f"     Total age2+ biomass: {nasc_proccessed['biomass'].sum() * 1e-6:.1f} kmt\n"
    if transect_length_age_biomass is not None:
        age_bin_coord = transect_length_age_biomass["age_bin"]

        def _age_bin_midpoint(value):
            if isinstance(value, pd.Interval):
                return float(value.mid)
            return float(value)

        age1_mask_values = np.array([np.isclose(_age_bin_midpoint(value), 1.0) for value in age_bin_coord.values])
        age1_mask = xr.DataArray(
            age1_mask_values,
            coords={age_bin_coord.dims[0]: age_bin_coord.values},
            dims=age_bin_coord.dims,
        )
        age1_leaks = float(transect_length_age_biomass.where(age1_mask, other=0.0).sum().item())
        age1_plus_total = float(transect_length_age_biomass.sum().item())
        age2_plus_total = max(age1_plus_total - age1_leaks, 0.0)

        transect_biomass_text = (
            f"     Total age1+ biomass: {age1_plus_total * 1e-6:.1f} kmt\n"
            f"     age1 biomass (age1 leaks): {age1_leaks * 1e-6:.1f} kmt\n"
            f"     Total age2+ biomass: {age2_plus_total * 1e-6:.1f} kmt\n"
        )

    output_text = (
        f"----------------------\n"
        f"Transect-based results\n"
        f"     Total NASC: {nasc_proccessed['nasc'].sum():.1f} m²nmi⁻²\n"
        f"     Mean number density: {nasc_proccessed['number_density'].mean():.1f} animals nmi⁻²\n"
        f"     Total abundance: {nasc_proccessed['abundance'].sum():.0f} fish\n"
        f"     Mean biomass density: {nasc_proccessed['biomass_density'].mean() * 1e-6:.3f} kmt nmi⁻²\n"
        f"{transect_biomass_text}"
        f"----------------------\n"
        f"Kriging-based results\n"
        f"     Total derived NASC: {kriged_results['nasc'].sum():.1f} m²nmi⁻²\n"
        f"     Total derived abundance: {kriged_results['abundance'].sum():.0f} fish\n"
        f"     Total age2+ biomass: {kriged_results['biomass'].sum() * 1e-6:.1f} kmt\n"
    )

    with open(target_dir / "summary_report.txt", "a", encoding="utf-8") as file:
        file.write(output_text)
    #logger.info(output_text)

def _create_reports(dict_df_bio, df_nasc_proc, df_kriged_results, ds_da_weight_dist,
                    dict_ds_kriged_abundance_table, da_kriged_biomass_table,
                    invert_hake, da_averaged_weight, output_dir,
                    len_age_abundance=None, len_age_biomass=None, selectivity_plot_data=None):
    """Create reports and export to output directory."""
    if config.reports.export_mode.lower() == "postgres":
        raise NotImplementedError("export_mode='postgres' is not implemented yet.")

    output_dir.mkdir(parents=True, exist_ok=True)

    reporter = Reporter(output_dir, verbose=True)

    if selectivity_plot_data is not None:
        plot_hauls_in_chunks(
            selectivity_plot_data,
            save_pdf_path=output_dir,
            skip_show=config.selectivity.skip_show_plots,
            output_stem="expansion_counts_comparison",
            logger=logger,
        )
        plot_selectivity_sum_to_one_density_by_haul(
            selectivity_plot_data,
            save_output_path=output_dir,
            skip_show=config.selectivity.skip_show_plots,
            output_stem="expansion_sum_to_one_density_comparison",
            logger=logger,
        )

    generate_summary_report(
        df_nasc_proc,
        df_kriged_results,
        output_dir,
        transect_length_age_biomass=len_age_biomass,
    )

    target_dir = output_dir / "config_files"
    target_dir.mkdir(parents=True, exist_ok=True)

    for config_file in config.yaml_configs:
        try:
            shutil.copy(config_file, target_dir)
        except FileNotFoundError:
            logger.warning("Config source file was not found: '%s'", config_file)
        except PermissionError:
            logger.error("Permission denied while copying config file: '%s'", config_file)

    # AGED-LENGTH HAUL
    reporter.aged_length_haul_counts_report(
        filename=config.output.aged_length_haul_counts,
        sheetnames=config.output.sheets,
        bio_data=dict_df_bio["specimen"].dropna(subset=["age", "length", "weight"])
    )

    # TOTAL LENGTH HAUL COUNTS
    reporter.total_length_haul_counts_report(
        filename=config.output.total_length_haul_counts,
        sheetnames=config.output.sheets,
        bio_data=dict_df_bio
    )

    # KRIGED AGED BIOMASS

    # All values
    reporter.kriged_aged_biomass_mesh_report(
        filename=config.output.kriged_aged_biomass_mesh_full,
        sheetnames=config.output.sheets_all,
        kriged_data=df_kriged_results,
        weight_data=ds_da_weight_dist["aged"],
        kriged_stratum_link={config.strata.geostratum_name: config.strata.stratum_name},
    )

    # Nonzero values
    reporter.kriged_aged_biomass_mesh_report(
        filename=config.output.kriged_aged_biomass_mesh_nonzero,
        sheetnames=config.output.sheets_all,
        kriged_data=df_kriged_results[df_kriged_results["biomass"] > 0.],
        weight_data=ds_da_weight_dist["aged"],
        kriged_stratum_link={config.strata.geostratum_name: config.strata.stratum_name},
    )

    # KRIGERD MESH RESULTS

    # All values
    reporter.kriged_mesh_results_report(
        filename=config.output.kriged_biomass_mesh_full,
        sheetname=config.output.sheet,
        kriged_data=df_kriged_results,
        kriged_stratum=config.strata.geostratum_name,
        kriged_variable="biomass",
        sigma_bs_data=invert_hake.sigma_bs_strata,
        sigma_bs_stratum=config.strata.stratum_name,
    )

    # Nonzero values
    reporter.kriged_mesh_results_report(
        filename=config.output.kriged_biomass_mesh_nonzero,
        sheetname=config.output.sheet,
        kriged_data=df_kriged_results[df_kriged_results["abundance"] > 0.],
        kriged_stratum=config.strata.geostratum_name,
        kriged_variable="biomass",
        sigma_bs_data=invert_hake.sigma_bs_strata,
        sigma_bs_stratum=config.strata.stratum_name,
    )

    # KRIGING INPUT
    reporter.kriging_input_report(
        filename=config.output.kriging_input_report,
        sheetname=config.output.sheet,
        transect_data=df_nasc_proc,
    )

    if len_age_abundance is not None:
        # TRANSECT LENGTH-AGE ABUNDANCES
        reporter.transect_length_age_abundance_report(
            filename=config.output.transect_length_age_abundance_report,
            sheetnames=config.output.sheets,
            datatables=len_age_abundance,
        )
        # KRIGED LENGTH-AGE ABUNDANCES
        reporter.kriged_length_age_abundance_report(
            filename=config.output.kriged_length_age_abundance_report,
            sheetnames=config.output.sheets,
            datatables=dict_ds_kriged_abundance_table,
        )

    if len_age_biomass is not None:
        # TRANSECT LENGTH-AGE BIOMASS
        reporter.transect_length_age_biomass_report(
            filename=config.output.transect_length_age_biomass_report,
            sheetnames=config.output.sheets,
            datatable=len_age_biomass,
        )
        # KRIGED LENGTH-AGE BIOMASS
        reporter.kriged_length_age_biomass_report(
            filename=config.output.kriged_length_age_biomass_report,
            sheetnames=config.output.sheets,
            datatable=da_kriged_biomass_table,
        )

        # TRANSECT AGED BIOMASS

        # Full values
        reporter.transect_aged_biomass_report(
            filename=config.output.transect_aged_biomass_report_full,
            sheetnames=config.output.sheets_all,
            transect_data=df_nasc_proc,
            weight_data=ds_da_weight_dist["aged"],
        )

        # Nonzero values
        reporter.transect_aged_biomass_report(
            filename=config.output.transect_aged_biomass_report_nonzero,
            sheetnames=config.output.sheets_all,
            transect_data=df_nasc_proc[df_nasc_proc["biomass"] > 0.],
            weight_data=ds_da_weight_dist["aged"],
        )

    # TRANSECT RESULTS

    # Full values
    reporter.transect_population_results_report(
        filename=config.output.transect_population_results_full,
        sheetname=config.output.sheet,
        transect_data=df_nasc_proc,
        weight_strata_data=da_averaged_weight,
        sigma_bs_stratum=invert_hake.sigma_bs_strata,
        stratum_name=config.strata.stratum_name,
    )

    # Nonzero values
    reporter.transect_population_results_report(
        filename=config.output.transect_population_results_nonzero,
        sheetname=config.output.sheet,
        transect_data=df_nasc_proc[df_nasc_proc["biomass"] > 0.],
        weight_strata_data=da_averaged_weight,
        sigma_bs_stratum=invert_hake.sigma_bs_strata,
        stratum_name=config.strata.stratum_name,
    )

    logger.info(f"Data export completed successfully. Output available in directory: '{output_dir}'.")

# ---------------------------------------------------------
# ORCHESTRATION
# ---------------------------------------------------------
@dataclass
class LoadedInputs:
    nasc_raw: pd.DataFrame
    biodata_raw: dict
    strata_tables: dict
    geostrata_tables: dict
    kriging_mesh_raw: pd.DataFrame
    isobath_table: pd.DataFrame
    variogram_parameters: dict


@dataclass
class PreparedInputs:
    biodata_stratified: dict
    nasc_stratified: pd.DataFrame
    kriging_mesh_stratified: pd.DataFrame
    binned_weight_table: xr.DataArray
    weight_distributions: xr.Dataset
    inversion_model: inversion.InversionLengthTS
    nasc_inverted: pd.DataFrame


def _stage_load_inputs() -> LoadedInputs:
    """Stage 1: Load all source tables needed by the full analysis."""
    logger.info("----------Stage 1: Loading input data----------")
    nasc_raw = _ingest_nasc()
    biodata_raw = _ingest_biodata()
    strata_tables, geostrata_tables = _ingest_stratification()
    kriging_mesh_raw, isobath_table, _kriging_params, variogram_parameters = _ingest_kriging()
    return LoadedInputs(
        nasc_raw=nasc_raw,
        biodata_raw=biodata_raw,
        strata_tables=strata_tables,
        geostrata_tables=geostrata_tables,
        kriging_mesh_raw=kriging_mesh_raw,
        isobath_table=isobath_table,
        variogram_parameters=variogram_parameters,
    )


def _stage_prepare_shared_inputs(loaded_inputs: LoadedInputs) -> PreparedInputs:
    """
    Stage 2: Build stratified and binned inputs reused by every model run.

    Why: this keeps model loops focused on model differences only.
    """
    logger.info("----------Stage 2: Preparing shared analysis inputs----------")
    biodata_stratified, nasc_stratified, kriging_mesh_stratified = _apply_stratification(
        loaded_inputs.biodata_raw,
        loaded_inputs.nasc_raw,
        loaded_inputs.kriging_mesh_raw,
        loaded_inputs.strata_tables,
        loaded_inputs.geostrata_tables,
    )
    biodata_stratified, _age_bins, length_bins = _binify_biodata(biodata_stratified)
    length_weight_coefs = _fit_length_weight_regressions(biodata_stratified)
    binned_weight_table = _compute_binned_weights(biodata_stratified, length_bins, length_weight_coefs)
    weight_distributions = _compute_binned_weight_distributions(biodata_stratified, binned_weight_table)
    inversion_model, nasc_inverted = _setup_inversion_and_invert(
        nasc_stratified, biodata_stratified, loaded_inputs.strata_tables
    )
    return PreparedInputs(
        biodata_stratified=biodata_stratified,
        nasc_stratified=nasc_stratified,
        kriging_mesh_stratified=kriging_mesh_stratified,
        binned_weight_table=binned_weight_table,
        weight_distributions=weight_distributions,
        inversion_model=inversion_model,
        nasc_inverted=nasc_inverted,
    )


def _run_single_model(
    model_name: str,
    run_timestamp: str,
    loaded_inputs: LoadedInputs,
    prepared_inputs: PreparedInputs,
):
    """Stages 3-5 for one output model."""
    logger.info(f"----------Processing model: {model_name}----------")
    output_dir = config.output_dir / f"{model_name}_{run_timestamp}"
    selectivity_plot_data = None

    if model_name == "net_selectivity":
        selectivity_plot_data = _apply_net_selectivity(
            prepared_inputs.biodata_stratified,
            prepared_inputs.binned_weight_table,
        )
        number_proportions, weight_proportions, model_weight_distributions = _compute_model_proportions(
            prepared_inputs.biodata_stratified,
            prepared_inputs.binned_weight_table,
            aged_data=selectivity_plot_data.dropna(subset=["age", "length", "weight"]),
            aged_count_col="selectivity_expansion",
            aged_agg_func="sum",
        )
    else:
        number_proportions, weight_proportions, model_weight_distributions = _compute_model_proportions(
            prepared_inputs.biodata_stratified,
            prepared_inputs.binned_weight_table,
            weight_distributions=prepared_inputs.weight_distributions,
        )

    (
        nasc_processed,
        transect_abundance,
        transect_biomass,
        average_weight,
    ) = _compute_biomass_and_abundance_estimates(
        prepared_inputs.nasc_inverted,
        number_proportions,
        prepared_inputs.binned_weight_table,
        weight_proportions,
    )

    logger.info("----------Stage 4: Geostatistics and kriging----------")
    best_fit_params, nasc_processed, kriging_mesh = _geostatistical_analysis(
        nasc_processed,
        prepared_inputs.kriging_mesh_stratified,
        loaded_inputs.isobath_table,
        loaded_inputs.variogram_parameters,
    )
    kriged_results = _kriging_analysis(
        best_fit_params,
        nasc_processed,
        kriging_mesh,
        weight_proportions,
        average_weight,
        prepared_inputs.inversion_model,
    )
    _jolly_and_hampton_analysis(nasc_processed, kriged_results, loaded_inputs.geostrata_tables)

    logger.info(
        "----------Stage 5: Exporting outputs----------"
    )
    _, kriged_biomass, kriged_abundance_tables = _distribute_population_estimates(
        kriged_results,
        number_proportions,
        weight_proportions,
    )
    _create_reports(
        prepared_inputs.biodata_stratified,
        nasc_processed,
        kriged_results,
        model_weight_distributions,
        kriged_abundance_tables,
        kriged_biomass,
        prepared_inputs.inversion_model,
        average_weight,
        output_dir,
        selectivity_plot_data=selectivity_plot_data,
        len_age_abundance=transect_abundance,
        len_age_biomass=transect_biomass,
    )

    logger.info("Done")
    return {
        "model_name": model_name,
        "output_dir": output_dir,
        "transect_abundance": transect_abundance,
        "transect_biomass": transect_biomass,
        "kriged_abundance": kriged_abundance_tables,
        "kriged_biomass": kriged_biomass,
    }


def main():
    run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    loaded_inputs = _stage_load_inputs()
    prepared_inputs = _stage_prepare_shared_inputs(loaded_inputs)

    logger.info("----------Stage 3: Running requested models----------")
    model_results_by_name = {}
    for model_name in config.models:
        model_results_by_name[model_name] = _run_single_model(
            model_name,
            run_timestamp,
            loaded_inputs,
            prepared_inputs,
        )

if __name__ == "__main__":
    main()
