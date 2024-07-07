import pandas as pd
import boto3
import os
import re
import logging
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import gc
from concurrent.futures import ThreadPoolExecutor

logging.basicConfig(level=logging.INFO)

required_env_vars = ["AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"]
for var in required_env_vars:
    if var not in os.environ:
        raise EnvironmentError(f"Environment variable {var} is not set")


def clean_chembl_id(chembl_id):
    match = re.match(r"(CHEMBL\d+)", chembl_id)
    if match:
        return match.group(1)
    digits = "".join(filter(str.isdigit, chembl_id))
    if digits:
        return f"CHEMBL{digits}"
    return chembl_id


def load_and_clean_data(file_path, has_smiles=True, chunksize=100000):
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


def download_files_from_s3(bucket_name, s3_folder, local_folder):
    session = boto3.Session(
        aws_access_key_id=os.environ["AWS_ACCESS_KEY_ID"],
        aws_secret_access_key=os.environ["AWS_SECRET_ACCESS_KEY"],
        aws_session_token=os.environ["AWS_SESSION_TOKEN"]
    )
    s3_client = session.client("s3")

    if not os.path.exists(local_folder):
        os.makedirs(local_folder)

    logging.info("Listing files in S3 bucket...")
    response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=s3_folder)
    file_paths = []
    for obj in response.get("Contents", []):
        key = obj["Key"]
        if key.endswith(".csv"):
            local_path = os.path.join(local_folder, os.path.basename(key))
            logging.info(f"Downloading {key} to {local_path}...")
            s3_client.download_file(bucket_name, key, local_path)
            logging.info(f"Downloaded {key} to {local_path}")
            file_paths.append(local_path)
    return file_paths


def upload_file_to_s3(local_path, s3_path, s3_bucket_name):
    s3_client = boto3.client("s3")
    try:
        s3_client.upload_file(local_path, s3_bucket_name, s3_path)
        logging.info(f"Uploaded {local_path} to s3://{s3_bucket_name}/{s3_path}")
    except Exception as e:
        logging.error(f"Failed to upload {local_path} to s3://{s3_bucket_name}/{s3_path}: {e}")


bucket_name = "de-school-2024-aws"
s3_folder = "final_task/input_files/"
local_folder = "input_files"
file_paths = download_files_from_s3(bucket_name, s3_folder, local_folder)

molecule_structures = pd.read_csv(os.path.join(local_folder, "molecule_structures.csv"))
logging.info(f"Loaded molecule_structures.csv with columns: {molecule_structures.columns.tolist()}")

molecule_dictionary = pd.read_csv(os.path.join(local_folder, "molecule_dictionary.csv"))
logging.info(f"Loaded molecule_dictionary.csv with columns: {molecule_dictionary.columns.tolist()}")

molecule_structures["chembl_id"] = molecule_dictionary["chembl_id"]
logging.info("Assigned chembl_ids to molecule_structures")

additional_data_list = []
for file_path in file_paths:
    additional_data_chunks = load_and_clean_data(file_path, has_smiles=True)
    for chunk in additional_data_chunks:
        additional_data_list.append(chunk)

if additional_data_list:
    additional_data = pd.concat(additional_data_list, ignore_index=True)
    molecule_structures = pd.concat([molecule_structures, additional_data], ignore_index=True)

molecule_structures = molecule_structures[molecule_structures["canonical_smiles"].notna()]
molecule_structures = molecule_structures[molecule_structures["canonical_smiles"] != "nan"]

valid_smiles = molecule_structures["canonical_smiles"].apply(lambda x: Chem.MolFromSmiles(x) is not None)
molecule_structures = molecule_structures[valid_smiles]


def compute_similarity_scores(selected_mol, all_mols, source_chembl_id, all_chembl_ids):
    selected_fp = AllChem.GetMorganFingerprintAsBitVect(selected_mol, 2, nBits=2048)
    all_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in all_mols]

    similarity_scores = [DataStructs.TanimotoSimilarity(selected_fp, all_fp) for all_fp in all_fps]

    similarity_df = pd.DataFrame({
        "source_chembl_id": source_chembl_id,
        "target_chembl_id": all_chembl_ids,
        "similarity_score": similarity_scores
    })

    return similarity_df


def find_top_10_similar(similarity_df):
    top_10_similar_df = similarity_df.nlargest(10, "similarity_score")
    has_duplicates = (top_10_similar_df["similarity_score"] == top_10_similar_df["similarity_score"].iloc[-1]).sum() > 1
    top_10_similar_df["has_duplicates_of_last_largest_score"] = has_duplicates
    return top_10_similar_df


chosen_molecules = molecule_structures["chembl_id"].sample(100, random_state=42).tolist()

s3_client = boto3.client("s3")
s3_bucket_name = "de-school-2024-aws"

output_folder = "output"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

selected_molecule_structures = molecule_structures[molecule_structures["chembl_id"].isin(chosen_molecules)]
all_mols = [Chem.MolFromSmiles(smile) for smile in molecule_structures["canonical_smiles"]]
all_chembl_ids = molecule_structures["chembl_id"].tolist()


def process_molecule(mol_id):
    selected_mol = Chem.MolFromSmiles(
        selected_molecule_structures[selected_molecule_structures["chembl_id"] == mol_id]["canonical_smiles"].values[0])
    if selected_mol is None:
        logging.warning(f"Skipping molecule {mol_id} due to invalid SMILES string.")
        return
    similarity_df = compute_similarity_scores(selected_mol, all_mols, mol_id, all_chembl_ids)

    similarity_scores_file = os.path.join(output_folder, f"similarity_scores_{mol_id}.parquet")
    similarity_df.to_parquet(similarity_scores_file, index=False)

    s3_similarity_scores_path = f"final_task/gaplanyan_anna/similarity_scores_{mol_id}.parquet"
    upload_file_to_s3(similarity_scores_file, s3_similarity_scores_path, s3_bucket_name)

    top_10_similar_df = find_top_10_similar(similarity_df)
    top_10_similar_file = os.path.join(output_folder, f"top_10_similar_{mol_id}.parquet")
    top_10_similar_df.to_parquet(top_10_similar_file, index=False)

    s3_top_10_similar_path = f"final_task/gaplanyan_anna/top_10_similar_{mol_id}.parquet"
    upload_file_to_s3(top_10_similar_file, s3_top_10_similar_path, s3_bucket_name)

    del selected_mol, similarity_df, top_10_similar_df
    gc.collect()


with ThreadPoolExecutor(max_workers=8) as executor:
    executor.map(process_molecule, chosen_molecules)

logging.info("Uploaded all results to S3.")
