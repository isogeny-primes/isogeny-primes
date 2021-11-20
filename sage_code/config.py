from pathlib import Path

SOURCE_DIR = Path(__file__).parent
PROJECT_DIR = SOURCE_DIR.parent
DATA_DIR = SOURCE_DIR.joinpath("data_files")

FORMAL_IMMERSION_DATA_AT_2_PATH = DATA_DIR.joinpath("formal_immersion_at_2.json")
BAD_FORMAL_IMMERSION_DATA_PATH = DATA_DIR.joinpath("bad_formal_immersion_data.json")
QUADRATIC_POINTS_DATA_PATH = DATA_DIR.joinpath("quadratic_points_catalogue.json")

SMALL_GONALITIES = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 59, 71}
