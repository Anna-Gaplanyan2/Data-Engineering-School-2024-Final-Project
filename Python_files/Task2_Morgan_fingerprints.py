import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool, cpu_count
import boto3
import os
import time
import logging
import pyarrow.parquet as pq
import pyarrow as pa

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

required_env_vars = ['AWS_ACCESS_KEY_ID', 'AWS_SECRET_ACCESS_KEY', 'AWS_SESSION_TOKEN']
for var in required_env_vars:
    if var not in os.environ:
        raise EnvironmentError(f'Environment variable {var} is not set')

def compute_fingerprint(smiles):
    if pd.isna(smiles):
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048).ToBitString()
    else:
        return None

def process_chunk(chunk):
    logging.info(f'Starting to compute fingerprints for a chunk of size {len(chunk)}')
    chunk['canonical_smiles'] = chunk['canonical_smiles'].fillna('')
    chunk['fingerprint'] = chunk['canonical_smiles'].apply(compute_fingerprint)
    logging.info('Finished computing fingerprints for chunk')
    return chunk

def upload_file_to_s3(file_path, retries=3):
    session = boto3.Session(
        aws_access_key_id=os.environ['AWS_ACCESS_KEY_ID'],
        aws_secret_access_key=os.environ['AWS_SECRET_ACCESS_KEY'],
        aws_session_token=os.environ['AWS_SESSION_TOKEN']
    )
    s3_client = session.client('s3')

    bucket_name = 'de-school-2024-aws'
    folder = 'final_task/gaplanyan_anna/'
    s3_key = os.path.join(folder, os.path.basename(file_path))

    attempt = 0
    while attempt < retries:
        try:
            logging.info(f'Uploading {file_path} to S3 (Attempt {attempt + 1})')
            s3_client.upload_file(file_path, bucket_name, s3_key)
            logging.info(f'Uploaded {file_path} to s3://{bucket_name}/{s3_key}')
            return
        except Exception as e:
            logging.error(f'Failed to upload {file_path} to S3 on attempt {attempt + 1}: {str(e)}')
            attempt += 1
            time.sleep(2 ** attempt)
    logging.error(f'Failed to upload {file_path} to S3 after {retries} attempts')

def process_file_in_chunks(file_path, chunk_size=182000):
    chunk_files = []
    for i, chunk in enumerate(pd.read_csv(file_path, chunksize=chunk_size)):
        logging.info(f'Processing chunk {i}')
        start_time = time.time()
        with Pool(cpu_count()) as pool:
            processed_chunk = pool.map(process_chunk, [chunk])[0]
        computation_time = time.time() - start_time
        logging.info(f'Chunk {i} computation time: {computation_time:.2f} seconds')

        local_chunk_file = f'processed_chunk_{i}.parquet'
        logging.info(f'Saving chunk {i} with fingerprints to {local_chunk_file}')
        table = pa.Table.from_pandas(processed_chunk)
        pq.write_table(table, local_chunk_file, compression='snappy')
        chunk_files.append(local_chunk_file)

        upload_start_time = time.time()
        upload_file_to_s3(local_chunk_file)
        upload_time = time.time() - upload_start_time
        logging.info(f'Chunk {i} upload time: {upload_time:.2f} seconds')

    return chunk_files

def merge_chunks(chunk_files):
    df_list = [pd.read_parquet(chunk_file) for chunk_file in chunk_files]
    return pd.concat(df_list, ignore_index=True)

if __name__ == '__main__':
    logging.info('Starting processing')

    file_path = 'molecule_structures.csv'
    processed_chunk_files = process_file_in_chunks(file_path)
    logging.info('Processing completed')
