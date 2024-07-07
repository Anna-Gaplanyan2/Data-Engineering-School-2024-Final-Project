import pandas as pd
import os
from sqlalchemy import create_engine
import logging
import time

logging.basicConfig(level=logging.INFO, format='%Y-%m-%d %H:%M:%S - %(levelname)s - %(message)s')

properties_file_path = "molecule_properties.csv"
dictionary_file_path = "input_files/molecule_dictionary.csv"
similarities_file_path = "combined_molecule_similarities.parquet"
dimension_table_csv_path = "input_files/dimension_table.csv"

logging.info("Starting the process to create the dimension table.")

logging.info("Loading CSV and Parquet files.")
properties_df = pd.read_csv(properties_file_path)
dictionary_df = pd.read_csv(dictionary_file_path)
similarities_df = pd.read_parquet(similarities_file_path)
logging.info("Files loaded successfully.")

dictionary_selected_df = dictionary_df[["chembl_id", "molecule_type"]]
properties_selected_df = properties_df[
    [
        "mw_freebase",
        "alogp",
        "psa",
        "cx_logp",
        "molecular_species",
        "full_mwt",
        "aromatic_rings",
        "heavy_atoms",
    ]
]

logging.info("Concatenating the columns from dictionary and properties data.")
concatenated_df = pd.concat(
    [dictionary_selected_df.reset_index(drop=True), properties_selected_df.reset_index(drop=True)],
    axis=1,
)
logging.info(f"Concatenated DataFrame contains {len(concatenated_df)} rows.")

logging.info("Filtering concatenated DataFrame based on source_chembl_id in the Parquet file.")
valid_chembl_ids = set(similarities_df["source_chembl_id"])
filtered_df = concatenated_df[concatenated_df["chembl_id"].isin(valid_chembl_ids)]
logging.info(f"Filtered DataFrame contains {len(filtered_df)} rows.")

logging.info("Validating the number of rows in the dimension table.")
valid_dictionary_count = len(dictionary_selected_df[dictionary_selected_df["chembl_id"].isin(valid_chembl_ids)])
if len(filtered_df) == valid_dictionary_count:
    logging.info("The number of rows in the filtered DataFrame matches the number of valid chembl_ids.")
else:
    logging.warning(
        f"The number of rows in the filtered DataFrame ({len(filtered_df)}) does not match the number of valid chembl_ids ({valid_dictionary_count})."
    )

logging.info(f"Saving the dimension table to {dimension_table_csv_path}.")
filtered_df.to_csv(dimension_table_csv_path, index=False)
logging.info("Dimension table saved locally.")

db_params = {
    "dbname": os.environ["dbname"],
    "user": os.environ["user"],
    "password": os.environ["password"],
    "host": os.environ["host"]
}

logging.info("Creating database engine.")
engine = create_engine(
    f"postgresql://{db_params['user']}:{db_params['password']}@{db_params['host']}/{db_params['dbname']}"
)

try:
    logging.info("Uploading the dimension table to PostgreSQL.")
    start_time = time.time()
    filtered_df.to_sql("molecule_dimension", engine, if_exists="append", index=False)
    end_time = time.time()
    logging.info(f"Dimension table uploaded to PostgreSQL successfully in {end_time - start_time:.2f} seconds.")
except Exception as e:
    logging.error(f"Error occurred while uploading to PostgreSQL: {e}")
