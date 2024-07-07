import boto3
import pandas as pd
import os
import logging
from io import BytesIO
from sqlalchemy import create_engine

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

bucket_name = "de-school-2024-aws"
prefix = "final_task/gaplanyan_anna/top_10_similar"

db_params = {
    "dbname": os.environ["dbname"],
    "user": os.environ["user"],
    "password": os.environ["password"],
    "host": os.environ["host"]
}

logging.info("Initializing boto3 session")
session = boto3.Session(
    aws_access_key_id=os.environ["AWS_ACCESS_KEY_ID"],
    aws_secret_access_key=os.environ["AWS_SECRET_ACCESS_KEY"],
    aws_session_token=os.environ["AWS_SESSION_TOKEN"]
)

logging.info("Initializing S3 client")
s3 = session.client("s3")

logging.info("Listing objects in the S3 bucket")
response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
files = [item["Key"] for item in response["Contents"] if item["Key"].startswith("final_task/gaplanyan_anna/top_10_similar") and item["Key"].endswith(".parquet")]

logging.info(f"Found {len(files)} files to process")

dataframes = []
for file in files:
    logging.info(f"Downloading file {file}")
    obj = s3.get_object(Bucket=bucket_name, Key=file)
    data = obj["Body"].read()
    df = pd.read_parquet(BytesIO(data))
    dataframes.append(df)

logging.info("Combining all dataframes")
combined_df = pd.concat(dataframes, ignore_index=True)

logging.info("Creating database connection string")
connection_string = f"postgresql+psycopg2://{db_params['user']}:{db_params['password']}@{db_params['host']}/{db_params['dbname']}"
engine = create_engine(connection_string)

table_name = "molecule_similarities"

logging.info("Uploading DataFrame to PostgreSQL")
combined_df.to_sql(table_name, engine, if_exists="replace", index=False)

logging.info("Data successfully uploaded to PostgreSQL.")
