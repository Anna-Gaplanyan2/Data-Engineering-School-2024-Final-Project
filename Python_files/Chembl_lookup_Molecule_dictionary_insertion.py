import pandas as pd
import os
import multiprocessing as mp
from sqlalchemy import create_engine, Table, Column, Integer, String, MetaData, Text, Float, exc
import logging
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

db_params = {
    "dbname": os.environ["dbname"],
    "user": os.environ["user"],
    "password": os.environ["password"],
    "host": os.environ["host"]
}

db_url = (
    f"postgresql+psycopg2://{db_params['user']}:{db_params['password']}@"
    f"{db_params['host']}/{db_params['dbname']}"
)

csv_files = {
    "chembl_id_lookup": "input_files/chembl_lookup.csv",
    "molecule_dictionary": "input_files/molecule_dictionary.csv"
}

schema_name = "agaplanyan"

def create_schema_and_tables(engine):
    metadata = MetaData(schema=schema_name)

    Table(
        "chembl_id_lookup", metadata,
        Column("chembl_id", String, primary_key=True),
        Column("entity_type", String),
        Column("status", String),
        Column("last_active", Integer),
        Column("resource_url", Text)
    )

    Table(
        "molecule_dictionary", metadata,
        Column("chembl_id", String, primary_key=True),
        Column("pref_name", String),
        Column("max_phase", Integer),
        Column("structure_type", String),
        Column("molecule_type", String),
        Column("first_approval", Integer),
        Column("oral", String),
        Column("parenteral", String),
        Column("topical", String),
        Column("black_box_warning", String),
        Column("natural_product", String),
        Column("prodrug", Float),
        Column("dosed_ingredient", String),
        Column("therapeutic_flag", String),
        Column("chirality", Integer),
        Column("usan_year", Integer),
        Column("usan_stem", String),
        Column("availability_type", Integer),
        Column("first_in_class", String),
        Column("inorganic_flag", String),
        Column("chebi_par_id", Integer),
        Column("polymer_flag", String),
        Column("usan_substem", String),
        Column("usan_stem_definition", String),
        Column("indication_class", String)
    )

    metadata.create_all(engine)
    logger.info(f"Tables created successfully in schema {schema_name}.")

def clean_data(df, table_name):
    boolean_columns = [
        "oral", "parenteral", "topical", "black_box_warning", "natural_product",
        "dosed_ingredient", "therapeutic_flag", "first_in_class", "inorganic_flag",
        "polymer_flag"
    ]
    integer_columns = [
        "max_phase", "chirality", "usan_year", "availability_type", "chebi_par_id"
    ]

    for col in boolean_columns:
        if col in df.columns:
            df[col] = df[col].astype(str).apply(lambda x: x if x in ["0", "1", "nan"] else "nan")

    for col in integer_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.replace({pd.NA: None, float("nan"): None})

    if "chembl_id" in df.columns:
        df = df[df["chembl_id"].notnull()]

    return df

def load_csv_to_db_chunk(db_url, table_name, df_chunk, retry_attempts=3):
    attempt = 0
    while attempt < retry_attempts:
        try:
            engine = create_engine(db_url)
            logger.info(f"Loading chunk into {schema_name}.{table_name} (Attempt {attempt + 1})")
            df_chunk = clean_data(df_chunk, table_name)
            if df_chunk.empty:
                logger.error(f"DataFrame is empty after cleaning. Skipping chunk for {table_name}.")
                return
            df_chunk.to_sql(table_name, engine, if_exists="append", index=False, schema=schema_name)
            logger.info(f"Chunk loaded into {schema_name}.{table_name}")
            return
        except exc.SQLAlchemyError as e:
            logger.error(f"Error loading chunk into {schema_name}.{table_name}: {e}")
            attempt += 1
            time.sleep(2)
        except Exception as e:
            logger.error(f"Unexpected error loading chunk into {schema_name}.{table_name}: {e}")
            return
    logger.error(f"Failed to load chunk into {schema_name}.{table_name} after {retry_attempts} attempts")

def load_data_multiprocessing(db_url, table_name, csv_path):
    num_processes = min(mp.cpu_count(), 4)
    chunk_size = 100000

    try:
        logger.info(f"Reading CSV file {csv_path}")
        df_iter = pd.read_csv(csv_path, chunksize=chunk_size, dtype=str)
    except Exception as e:
        logger.error(f"Error reading CSV file {csv_path}: {e}")
        return

    try:
        with mp.Pool(processes=num_processes) as pool:
            logger.info(f"Starting multiprocessing for {table_name}")
            for df_chunk in df_iter:
                logger.info(f"Submitting chunk for {table_name}")
                pool.apply_async(load_csv_to_db_chunk, args=(db_url, table_name, df_chunk))

            pool.close()
            pool.join()
            logger.info(f"Finished multiprocessing for {table_name}")
    except Exception as e:
        logger.error(f"Error during multiprocessing: {e}")
        pool.terminate()
    finally:
        pool.close()
        pool.join()

if __name__ == "__main__":
    try:
        engine = create_engine(db_url)
        create_schema_and_tables(engine)
        for table_name, csv_path in csv_files.items():
            logger.info(f"Loading data from {csv_path} into {schema_name}.{table_name}")
            load_data_multiprocessing(db_url, table_name, csv_path)
            logger.info(f"Finished loading data from {csv_path} into {schema_name}.{table_name}")
    except Exception as e:
        logger.error(f"Error: {e}")
