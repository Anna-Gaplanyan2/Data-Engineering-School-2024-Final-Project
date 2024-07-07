import requests
import xml.etree.ElementTree as ET
import pandas as pd
import logging
from multiprocessing import Pool, cpu_count
from time import sleep

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

def fetch_data(url):
    results = []
    logging.info(f'Fetching data from: {url}')
    retries = 3
    for _ in range(retries):
        try:
            response = requests.get(url)
            if response.status_code == 200:
                root = ET.fromstring(response.content)
                for item in root.findall('.//chembl_id_lookup'):
                    record = {
                        'chembl_id': item.find('chembl_id').text,
                        'entity_type': item.find('entity_type').text,
                        'status': item.find('status').text,
                        'last_active': item.find('last_active').text if item.find('last_active') is not None else 'N/A',
                        'resource_url': item.find('resource_url').text
                    }
                    results.append(record)
                break
            else:
                logging.error(f'Failed to fetch data from: {url} with status code {response.status_code}')
        except requests.RequestException as e:
            logging.error(f'Request failed: {e}')
        sleep(2)
    return results

def fetch_all_data(endpoint, total_records, start_page=1, limit=1000):
    all_results = []
    total_pages = (total_records // limit) + 1
    urls = [f'{BASE_URL}/{endpoint}?limit={limit}&offset={i * limit}' for i in range(start_page - 1, total_pages)]

    with Pool(cpu_count()) as pool:
        for i, data in enumerate(pool.imap_unordered(fetch_data, urls)):
            all_results.extend(data)
            logging.info(f'Fetched {len(data)} records from page {i + 1}')

    return all_results

if __name__ == '__main__':
    response = requests.get(f'{BASE_URL}/chembl_id_lookup?limit=1')
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        total_records = int(root.find('.//total_count').text)
        logging.info(f'Total records to fetch: {total_records}')
    else:
        logging.error('Failed to fetch the total number of records.')
        total_records = 0

    all_data = fetch_all_data('chembl_id_lookup', total_records, start_page=1, limit=1000)

    df = pd.DataFrame(all_data)
    df.to_csv('chembl_id_lookup_data_with_retires.csv', index=False)
