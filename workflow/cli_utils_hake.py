import argparse

def get_verbose():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args, unknown = parser.parse_known_args()
    return args.verbose

def get_compare():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compare", action="store_false", help="Compare generated reports to EchoPro"
    )
    args, unknown = parser.parse_known_args()
    return args.compare

def get_year():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, default=2015, help="Set year")
    args, unknown = parser.parse_known_args()
    return args.year

def get_extrap_flag():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--extrap",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable or disable extrapolation",
    )
    args, unknown = parser.parse_known_args()
    return args.extrap

def get_strata_type():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--strata",
        choices=["ks", "inpfc"],
        default="ks",
        help="Set stratification type",
    )
    args, unknown = parser.parse_known_args()
    return args.strata