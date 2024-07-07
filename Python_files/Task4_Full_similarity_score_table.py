# Not done because of my laptop limitations
import boto3
import os
import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import logging
import pyarrow.parquet as pq
import pyarrow as pa
from multiprocessing import Pool, cpu_count
import faiss
import numpy as np
import asyncio
import aiofiles
import aioboto3

logging.basicConfig(level=logging.INFO, format="%Y-%m-%d %H:%M:%S - %(levelname)s - %(message)s")

required_env_vars = ["AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"]
for var in required_env_vars:
    if var not in os.environ:
        raise EnvironmentError(f"Environment variable {var} is not set")

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

def clean_chembl_id(chembl_id):
    match = re.match(r"(CHEMBL\d+)", chembl_id)
    if match:
        return match.group(1)
    digits = "".join(filter(str.isdigit, chembl_id))
    if digits:
        return f"CHEMBL{digits}"
    return chembl_id

def load_and_clean_data(file_path, has_smiles=True, chunksize=182000):
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

def compute_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    else:
        logging.warning(f"Invalid SMILES string: {smiles}")
        return None

def compute_fingerprints(df, smiles_column="canonical_smiles"):
    logging.info("Computing fingerprints for dataframe...")
    with Pool(cpu_count()) as pool:
        df["fingerprint"] = pool.map(compute_fingerprint, df[smiles_column])
    df.dropna(subset=["fingerprint"], inplace=True)
    logging.info("Computed fingerprints for dataframe")
    return df

def fingerprint_to_array(fp):
    arr = np.zeros((1,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def build_index(chunk):
    f = len(chunk[0])
    index = faiss.IndexFlatL2(f)
    faiss.normalize_L2(chunk)
    index.add(chunk)
    return index

def build_faiss_index_parallel(fingerprint_chunks):
    num_cores = cpu_count()
    logging.info(f"Building FAISS index using {num_cores} cores...")
    with Pool(num_cores) as pool:
        indices = pool.map(build_index, fingerprint_chunks)

    combined_index = faiss.IndexFlatL2(indices[0].d)
    for index in indices:
        combined_index.add(index.reconstruct_n(0, index.ntotal))

    return combined_index

def query_faiss_index(index, fingerprints, num_neighbors=10):
    distances, indices = index.search(fingerprints, num_neighbors)
    return distances, indices

def process_similarity_query(input_df, molecule_structures, faiss_index, num_neighbors=10):
    input_fingerprints = np.array([fingerprint_to_array(fp) for fp in input_df["fingerprint"]]).astype("float32")
    if input_fingerprints.shape[0] == 0:
        logging.warning("No valid fingerprints found in input data. Skipping similarity calculation.")
        return pd.DataFrame()

    distances, indices = query_faiss_index(faiss_index, input_fingerprints, num_neighbors=num_neighbors)

    results = []
    for i, (d, idx) in enumerate(zip(distances, indices)):
        if isinstance(d, (list, np.ndarray)) and isinstance(idx, (list, np.ndarray)):
            for distance, index in zip(d, idx):
                results.append({
                    "source_chembl_id": input_df.iloc[i]["chembl_id"],
                    "target_chembl_id": molecule_structures.iloc[index]["chembl_id"],
                    "similarity_score": 1 - distance
                })
    return pd.DataFrame(results)

def save_to_parquet(df, file_path):
    logging.info(f"Saving dataframe to Parquet file {file_path}...")
    table = pa.Table.from_pandas(df)
    pq.write_table(table, file_path, compression="snappy")
    logging.info(f"Saved dataframe to Parquet file {file_path}")

async def upload_file_to_s3(file_path, bucket_name, s3_key):
    logging.info(f"Uploading {file_path} to s3://{bucket_name}/{s3_key}...")
    session = aioboto3.Session(
        aws_access_key_id=os.environ["AWS_ACCESS_KEY_ID"],
        aws_secret_access_key=os.environ["AWS_SECRET_ACCESS_KEY"],
        aws_session_token=os.environ["AWS_SESSION_TOKEN"]
    )
    async with session.client("s3") as s3_client:
        try:
            await s3_client.head_object(Bucket=bucket_name, Key=s3_key)
            logging.info(f"File {s3_key} exists in bucket {bucket_name}. It will be overwritten.")
        except:
            logging.info(f"File {s3_key} does not exist in bucket {bucket_name}. It will be created.")
        async with aiofiles.open(file_path, "rb") as f:
            await s3_client.upload_fileobj(f, bucket_name, s3_key)
    logging.info(f"Uploaded {file_path} to s3://{bucket_name}/{s3_key}")

async def save_and_upload_batch(similarity_scores_df, batch_index, bucket_name):
    result_file = f"similarity_scores_batch_{batch_index}.parquet"
    save_to_parquet(similarity_scores_df, result_file)
    await upload_file_to_s3(result_file, bucket_name, f"final_task/gaplanyan_anna/{os.path.basename(result_file)}")

async def process_and_upload_similarity_scores(input_df, molecule_structures, faiss_index, bucket_name, batch_size=10000):
    tasks = []
    batch_index = 0
    results = []

    for chembl_id in input_df["chembl_id"].unique():
        logging.info(f"Computing similarity scores for {chembl_id}...")
        molecule_df = input_df[input_df["chembl_id"] == chembl_id]
        similarity_scores_df = process_similarity_query(molecule_df, molecule_structures, faiss_index)

        if similarity_scores_df.empty:
            logging.warning(f"No similarity scores computed for {chembl_id}. Skipping.")
            continue

        results.append(similarity_scores_df)

        if len(results) >= batch_size:
            combined_results = pd.concat(results)
            tasks.append(save_and_upload_batch(combined_results, batch_index, bucket_name))
            batch_index += 1
            results = []

    if results:
        combined_results = pd.concat(results)
        tasks.append(save_and_upload_batch(combined_results, batch_index, bucket_name))

    await asyncio.gather(*tasks)

if __name__ == "__main__":
    bucket_name = "de-school-2024-aws"
    s3_folder = "final_task/input_files/"
    local_folder = "input_files/"

    logging.info("Starting download of files from S3...")
    file_paths = download_files_from_s3(bucket_name, s3_folder, local_folder)
    logging.info("Completed download of files from S3")

    molecule_dictionary = pd.read_csv("input_files/molecule_dictionary.csv")
    molecule_structures = pd.concat(load_and_clean_data("input_files/molecule_structures.csv", has_smiles=True))
    molecule_structures = molecule_structures.merge(molecule_dictionary[["chembl_id"]], left_index=True, right_index=True, how="left")
    molecule_structures = compute_fingerprints(molecule_structures)
    fingerprints = np.array([fingerprint_to_array(fp) for fp in molecule_structures["fingerprint"]]).astype("float32")
    chunk_size = 182000
    fingerprint_chunks = [fingerprints[i:i + chunk_size] for i in range(0, len(fingerprints), chunk_size)]
    logging.info("Building FAISS index in parallel...")
    faiss_index = build_faiss_index_parallel(fingerprint_chunks)
    all_molecules_df = molecule_structures.copy()

    for input_file in file_paths:
        logging.info(f"Processing file {input_file}...")
        for input_chunk in load_and_clean_data(input_file, has_smiles=True):
            input_chunk = compute_fingerprints(input_chunk)
            all_molecules_df = pd.concat([all_molecules_df, input_chunk], ignore_index=True)

    all_fingerprints = np.array([fingerprint_to_array(fp) for fp in all_molecules_df["fingerprint"]]).astype("float32")
    logging.info("Updating FAISS index with combined fingerprints...")
    all_fingerprint_chunks = [all_fingerprints[i:i + chunk_size] for i in range(0, len(all_fingerprints), chunk_size)]
    faiss_index = build_faiss_index_parallel(all_fingerprint_chunks)

    asyncio.run(process_and_upload_similarity_scores(all_molecules_df, molecule_structures, faiss_index, bucket_name, batch_size=10000))

    logging.info("All files processed and uploaded to S3")
