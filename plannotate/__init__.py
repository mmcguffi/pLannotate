__version__ = "1.2.2"

import warnings

import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
pd.set_option("future.no_silent_downcasting", True)
