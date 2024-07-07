import pandas as pd
import os
import re
import logging
from sqlalchemy import create_engine, Table, Column, String, MetaData, Text, PrimaryKeyConstraint
from sqlalchemy.dialects.postgresql import VARCHAR
from sqlalchemy.exc import IntegrityError

logging.basicConfig(level=logging.INFO)

db_params = {
    "dbname": os.environ["dbname"],
    "user": os.environ["user"],
    "password": os.environ["password"],
    "host": os.environ["host"]
}

db_url = f"postgresql+psycopg2://{db_params['user']}:{db_params['password']}@{db_params['host']}/{db_params['dbname']}"

csv_files = [
    "input_files/data_03_2024.csv",
    "input_files/data_04_2024.csv",
    "input_files/data_05_2024.csv",
    "input_files/data_06_2024.csv"
]

def clean_chembl_id(chembl_id):
    match = re.match(r"(CHEMBL\d+)", chembl_id)
    if match:
        return match.group(1)
    digits = "".join(filter(str.isdigit, chembl_id))
    if digits:
        return f"CHEMBL{digits}"
    return chembl_id

def load_and_clean_data(file_path, has_smiles=True, chunksize=10000):
    if os.path.getsize(file_path) == 0:
        logging.warning(f"File {file_path} is empty. Skipping.")
        return pd.DataFrame()

    for chunk in pd.read_csv(file_path, on_bad_lines="skip", encoding="latin1", chunksize=chunksize):
        chunk.columns = [col.strip().lower() for col in chunk.columns]
        if "molecule name" in chunk.columns:
            chunk.rename(columns={"molecule name": "chembl_id", "smiles": "canonical_smiles"}, inplace=True)
        if "chembl_id" in chunk.columns:
            chunk["chembl_id"] = chunk["chembl_id"].apply(lambda x: clean_chembl_id(x.strip()))
        if has_smiles and "canonical_smiles" in chunk.columns:
            chunk["canonical_smiles"] = chunk["canonical_smiles"].apply(lambda x: x.strip() if pd.notna(x) else None)
            chunk.dropna(subset=["canonical_smiles"], inplace=True)
        yield chunk

def create_schema_and_table(engine):
    metadata = MetaData()
    staging_new_compounds = Table(
        "staging_new_compounds", metadata,
        Column("chembl_id", VARCHAR, nullable=False),
        Column("canonical_smiles", Text),
        Column("date_added", String, nullable=False),
        PrimaryKeyConstraint("chembl_id", "date_added")
    )
    metadata.create_all(engine)
    logging.info("Schema and table created.")

def load_csv_to_db_chunk(engine, df_chunk):
    logging.info(f"Loading chunk with {len(df_chunk)} records into the database.")
    df_chunk.to_sql("staging_new_compounds", engine, if_exists="append", index=False)

if __name__ == "__main__":
    try:
        engine = create_engine(db_url)
        create_schema_and_table(engine)
        for file_path in csv_files:
            if not os.path.isfile(file_path):
                logging.error(f"File {file_path} does not exist.")
                continue
            file_name = os.path.basename(file_path)
            logging.info(f"Processing file: {file_name}")
            match = re.search(r"data_(\d{2})_(\d{4})\.csv", file_name)
            if match:
                month, year = match.groups()
                date_added = f"{year}-{month}"
                logging.info(f"Extracted date_added: {date_added}")
                for chunk in load_and_clean_data(file_path):
                    chunk["date_added"] = date_added
                    logging.info(f"Processed chunk with {len(chunk)} records.")
                    load_csv_to_db_chunk(engine, chunk)
            else:
                logging.error(f"Filename {file_name} does not match expected pattern.")
    except Exception as e:
        logging.error(f"Error: {e}")
