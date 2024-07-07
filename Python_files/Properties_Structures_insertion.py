import pandas as pd
import os
from sqlalchemy import create_engine, Table, Column, Integer, String, MetaData, Text, Float, exc
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

db_params = {
    "dbname": os.environ["dbname"],
    "user": os.environ["user"],
    "password": os.environ["password"],
    "host": os.environ["host"]
}

db_url = f"postgresql+psycopg2://{db_params['user']}:{db_params['password']}@{db_params['host']}/{db_params['dbname']}"

csv_files = {
    "compound_properties": "molecule_properties.csv",
    "compound_structures": "input_files/molecule_structures_with_chembl_id.csv"
}

schema_name = "agaplanyan"

def create_schema_and_tables(engine):
    metadata = MetaData(schema=schema_name)

    compound_properties = Table(
        "compound_properties", metadata,
        Column("molregno", Integer, primary_key=True),
        Column("alogp", Float),
        Column("aromatic_rings", Integer),
        Column("cx_logd", Float),
        Column("cx_logp", Float),
        Column("cx_most_apka", Float),
        Column("full_molformula", String),
        Column("full_mwt", Float),
        Column("hba", Integer),
        Column("hba_lipinski", Integer),
        Column("hbd", Integer),
        Column("hbd_lipinski", Integer),
        Column("heavy_atoms", Integer),
        Column("molecular_species", String),
        Column("mw_freebase", Float),
        Column("mw_monoisotopic", Float),
        Column("np_likeness_score", Float),
        Column("psa", Float),
        Column("qed_weighted", Float),
        Column("ro3_pass", String),
        Column("rtb", Integer),
        Column("num_lipinski_ro5_violations", Integer),
        Column("num_ro5_violations", Integer),
        Column("cx_most_bpka", Float)
    )

    compound_structures = Table(
        "compound_structures", metadata,
        Column("chembl_id", String, primary_key=True),
        Column("canonical_smiles", Text),
        Column("standard_inchi", Text),
        Column("standard_inchi_key", Text),
        Column("molfile", Text)
    )

    metadata.create_all(engine)
    logger.info(f"Tables created successfully in schema {schema_name}.")

def clean_data(df, table_name):
    boolean_columns = ["oral", "parenteral", "topical", "black_box_warning", "natural_product",
                       "dosed_ingredient", "therapeutic_flag", "first_in_class", "inorganic_flag",
                       "polymer_flag", "ro3_pass"]
    integer_columns = ["max_phase", "chirality", "usan_year", "availability_type", "chebi_par_id"]

    if table_name == "compound_properties":
        boolean_columns = []
        integer_columns = []

    for col in boolean_columns:
        if col in df.columns:
            df[col] = df[col].astype(str).apply(lambda x: x if x in ["0", "1", "nan"] else "nan")

    for col in integer_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.replace({pd.NA: None, float("nan"): None})

    if table_name != "compound_properties" and "chembl_id" in df.columns:
        df = df[df["chembl_id"].notnull()]

    return df

def load_csv_to_db(db_url, table_name, csv_path):
    try:
        engine = create_engine(db_url)
        logger.info(f"Loading data from {csv_path} into {schema_name}.{table_name}")
        df = pd.read_csv(csv_path, dtype=str)
        df = clean_data(df, table_name)
        if df.empty:
            logger.error(f"DataFrame is empty after cleaning. Skipping load for {table_name}.")
            return
        df.to_sql(table_name, engine, if_exists="append", index=False, schema=schema_name)
        logger.info(f"Data loaded into {schema_name}.{table_name} from {csv_path}")
    except exc.SQLAlchemyError as e:
        logger.error(f"Error loading data into {schema_name}.{table_name} from {csv_path}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error loading data into {schema_name}.{table_name} from {csv_path}: {e}")

if __name__ == "__main__":
    try:
        engine = create_engine(db_url)
        create_schema_and_tables(engine)
        for table_name, csv_path in csv_files.items():
            load_csv_to_db(db_url, table_name, csv_path)
    except Exception as e:
        logger.error(f"Error: {e}")
