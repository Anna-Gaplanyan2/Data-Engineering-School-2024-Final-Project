import requests
import xml.etree.ElementTree as ET
import pandas as pd
from multiprocessing import Pool, cpu_count
import logging


logging.basicConfig(level=logging.INFO, format="%Y-%m-%d %H:%M:%S - %(levelname)s - %(message)s")


BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"


def fetch_data(url, params=None):
    if params is None:
        params = {}
    results = []
    logging.info(f"Fetching data from: {url}")
    response = requests.get(url, params=params)
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        for molecule in root.findall(".//molecule"):
            record = {}
            properties = molecule.find(".//molecule_properties")
            structures = molecule.find(".//molecule_structures")

            record.update({
                "chembl_id": molecule.find("molecule_chembl_id").text if molecule.find(
                    "molecule_chembl_id") is not None else "N/A",
                "pref_name": molecule.find("pref_name").text if molecule.find("pref_name") is not None else "N/A",
                "max_phase": molecule.find("max_phase").text if molecule.find("max_phase") is not None else "N/A",
                "structure_type": molecule.find("structure_type").text if molecule.find(
                    "structure_type") is not None else "N/A",
                "molecule_type": molecule.find("molecule_type").text if molecule.find(
                    "molecule_type") is not None else "N/A",
                "first_approval": molecule.find("first_approval").text if molecule.find(
                    "first_approval") is not None else "N/A",
                "oral": molecule.find("oral").text if molecule.find("oral") is not None else "N/A",
                "parenteral": molecule.find("parenteral").text if molecule.find("parenteral") is not None else "N/A",
                "topical": molecule.find("topical").text if molecule.find("topical") is not None else "N/A",
                "black_box_warning": molecule.find("black_box_warning").text if molecule.find(
                    "black_box_warning") is not None else "N/A",
                "natural_product": molecule.find("natural_product").text if molecule.find(
                    "natural_product") is not None else "N/A",
                "prodrug": molecule.find("prodrug").text if molecule.find("prodrug") is not None else "N/A",
                "dosed_ingredient": molecule.find("dosed_ingredient").text if molecule.find(
                    "dosed_ingredient") is not None else "N/A",
                "therapeutic_flag": molecule.find("therapeutic_flag").text if molecule.find(
                    "therapeutic_flag") is not None else "N/A",
                "chirality": molecule.find("chirality").text if molecule.find("chirality") is not None else "N/A",
                "usan_year": molecule.find("usan_year").text if molecule.find("usan_year") is not None else "N/A",
                "usan_stem": molecule.find("usan_stem").text if molecule.find("usan_stem") is not None else "N/A",
                "availability_type": molecule.find("availability_type").text if molecule.find(
                    "availability_type") is not None else "N/A",
                "first_in_class": molecule.find("first_in_class").text if molecule.find(
                    "first_in_class") is not None else "N/A",
                "inorganic_flag": molecule.find("inorganic_flag").text if molecule.find(
                    "inorganic_flag") is not None else "N/A",
                "chebi_par_id": molecule.find("chebi_par_id").text if molecule.find(
                    "chebi_par_id") is not None else "N/A",
                "polymer_flag": molecule.find("polymer_flag").text if molecule.find(
                    "polymer_flag") is not None else "N/A",
                "usan_substem": molecule.find("usan_substem").text if molecule.find(
                    "usan_substem") is not None else "N/A",
                "usan_stem_definition": molecule.find("usan_stem_definition").text if molecule.find(
                    "usan_stem_definition") is not None else "N/A",
                "indication_class": molecule.find("indication_class").text if molecule.find(
                    "indication_class") is not None else "N/A"
            })

            if properties is not None:
                record.update({
                    "alogp": properties.find("alogp").text if properties.find("alogp") is not None else "N/A",
                    "aromatic_rings": properties.find("aromatic_rings").text if properties.find(
                        "aromatic_rings") is not None else "N/A",
                    "cx_logd": properties.find("cx_logd").text if properties.find("cx_logd") is not None else "N/A",
                    "cx_logp": properties.find("cx_logp").text if properties.find("cx_logp") is not None else "N/A",
                    "cx_most_apka": properties.find("cx_most_apka").text if properties.find(
                        "cx_most_apka") is not None else "N/A",
                    "full_molformula": properties.find("full_molformula").text if properties.find(
                        "full_molformula") is not None else "N/A",
                    "full_mwt": properties.find("full_mwt").text if properties.find("full_mwt") is not None else "N/A",
                    "hba": properties.find("hba").text if properties.find("hba") is not None else "N/A",
                    "hba_lipinski": properties.find("hba_lipinski").text if properties.find(
                        "hba_lipinski") is not None else "N/A",
                    "hbd": properties.find("hbd").text if properties.find("hbd") is not None else "N/A",
                    "hbd_lipinski": properties.find("hbd_lipinski").text if properties.find(
                        "hbd_lipinski") is not None else "N/A",
                    "heavy_atoms": properties.find("heavy_atoms").text if properties.find(
                        "heavy_atoms") is not None else "N/A",
                    "molecular_species": properties.find("molecular_species").text if properties.find(
                        "molecular_species") is not None else "N/A",
                    "mw_freebase": properties.find("mw_freebase").text if properties.find(
                        "mw_freebase") is not None else "N/A",
                    "mw_monoisotopic": properties.find("mw_monoisotopic").text if properties.find(
                        "mw_monoisotopic") is not None else "N/A",
                    "np_likeness_score": properties.find("np_likeness_score").text if properties.find(
                        "np_likeness_score") is not None else "N/A",
                    "psa": properties.find("psa").text if properties.find("psa") is not None else "N/A",
                    "qed_weighted": properties.find("qed_weighted").text if properties.find(
                        "qed_weighted") is not None else "N/A",
                    "ro3_pass": properties.find("ro3_pass").text if properties.find("ro3_pass") is not None else "N/A",
                    "rtb": properties.find("rtb").text if properties.find("rtb") is not None else "N/A",
                    "num_lipinski_ro5_violations": properties.find(
                        "num_lipinski_ro5_violations").text if properties.find(
                        "num_lipinski_ro5_violations") is not None else "N/A",
                    "num_ro5_violations": properties.find("num_ro5_violations").text if properties.find(
                        "num_ro5_violations") is not None else "N/A",
                    "cx_most_bpka": properties.find("cx_most_bpka").text if properties.find(
                        "cx_most_bpka") is not None else "N/A"
                })

            if structures is not None:
                record.update({
                    "canonical_smiles": structures.find("canonical_smiles").text if structures.find(
                        "canonical_smiles") is not None else "N/A",
                    "standard_inchi": structures.find("standard_inchi").text if structures.find(
                        "standard_inchi") is not None else "N/A",
                    "standard_inchi_key": structures.find("standard_inchi_key").text if structures.find(
                        "standard_inchi_key") is not None else "N/A",
                    "molfile": structures.find("molfile").text if structures.find("molfile") is not None else "N/A"
                })

            results.append(record)
    else:
        logging.error(f"Failed to fetch data from: {url}")
    return results


def fetch_all_data(endpoint, total_records, start_page=1, limit=500):
    all_results = []
    total_pages = (total_records // limit) + 1
    urls = [f"{BASE_URL}/{endpoint}?limit={limit}&offset={i * limit}" for i in range(start_page - 1, total_pages)]

    with Pool(cpu_count()) as pool:
        for i, data in enumerate(pool.imap_unordered(fetch_data, urls)):
            all_results.extend(data)
            logging.info(f"Fetched {len(data)} records from page {i + 1}")

    return all_results


if __name__ == "__main__":
    response = requests.get(f"{BASE_URL}/molecule?limit=1")
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        total_records = int(root.find(".//total_count").text)
        logging.info(f"Total records to fetch: {total_records}")
    else:
        logging.error("Failed to fetch the total number of records.")
        total_records = 0

    all_data = fetch_all_data("molecule", total_records, start_page=1, limit=500)

    molecule_dictionary = []
    molecule_properties = []
    molecule_structures = []

    for molecule in all_data:
        molecule_dictionary.append({
            "chembl_id": molecule.get("chembl_id"),
            "pref_name": molecule.get("pref_name"),
            "max_phase": molecule.get("max_phase"),
            "structure_type": molecule.get("structure_type"),
            "molecule_type": molecule.get("molecule_type"),
            "first_approval": molecule.get("first_approval"),
            "oral": molecule.get("oral"),
            "parenteral": molecule.get("parenteral"),
            "topical": molecule.get("topical"),
            "black_box_warning": molecule.get("black_box_warning"),
            "natural_product": molecule.get("natural_product"),
            "prodrug": molecule.get("prodrug"),
            "dosed_ingredient": molecule.get("dosed_ingredient"),
            "therapeutic_flag": molecule.get("therapeutic_flag"),
            "chirality": molecule.get("chirality"),
            "usan_year": molecule.get("usan_year"),
            "usan_stem": molecule.get("usan_stem"),
            "availability_type": molecule.get("availability_type"),
            "first_in_class": molecule.get("first_in_class"),
            "inorganic_flag": molecule.get("inorganic_flag"),
            "chebi_par_id": molecule.get("chebi_par_id"),
            "polymer_flag": molecule.get("polymer_flag"),
            "usan_substem": molecule.get("usan_substem"),
            "usan_stem_definition": molecule.get("usan_stem_definition"),
            "indication_class": molecule.get("indication_class"),
        })

        molecule_properties.append({
            "alogp": molecule.get("alogp"),
            "aromatic_rings": molecule.get("aromatic_rings"),
            "cx_logd": molecule.get("cx_logd"),
            "cx_logp": molecule.get("cx_logp"),
            "cx_most_apka": molecule.get("cx_most_apka"),
            "full_molformula": molecule.get("full_molformula"),
            "full_mwt": molecule.get("full_mwt"),
            "hba": molecule.get("hba"),
            "hba_lipinski": molecule.get("hba_lipinski"),
            "hbd": molecule.get("hbd"),
            "hbd_lipinski": molecule.get("hbd_lipinski"),
            "heavy_atoms": molecule.get("heavy_atoms"),
            "molecular_species": molecule.get("molecular_species"),
            "mw_freebase": molecule.get("mw_freebase"),
            "mw_monoisotopic": molecule.get("mw_monoisotopic"),
            "np_likeness_score": molecule.get("np_likeness_score"),
            "psa": molecule.get("psa"),
            "qed_weighted": molecule.get("qed_weighted"),
            "ro3_pass": molecule.get("ro3_pass"),
            "rtb": molecule.get("rtb"),
            "num_lipinski_ro5_violations": molecule.get("num_lipinski_ro5_violations"),
            "num_ro5_violations": molecule.get("num_ro5_violations"),
            "cx_most_bpka": molecule.get("cx_most_bpka")
        })

        molecule_structures.append({
            "canonical_smiles": molecule.get("canonical_smiles"),
            "standard_inchi": molecule.get("standard_inchi"),
            "standard_inchi_key": molecule.get("standard_inchi_key"),
            "molfile": molecule.get("molfile")
        })

    df_molecule_dictionary = pd.DataFrame(molecule_dictionary)
    df_molecule_properties = pd.DataFrame(molecule_properties)
    df_molecule_structures = pd.DataFrame(molecule_structures)

    logging.info("Molecule Dictionary")
    logging.info(df_molecule_dictionary.head())

    logging.info("Molecule Properties")
    logging.info(df_molecule_properties.head())

    logging.info("Molecule Structures")
    logging.info(df_molecule_structures.head())

    df_molecule_dictionary.to_csv("molecule_dictionary.csv", index=False)
    df_molecule_properties.to_csv("molecule_properties.csv", index=False)
    df_molecule_structures.to_csv("molecule_structures.csv", index=False)

    a = pd.read_csv("input_files/molecule_dictionary.csv")
    print(a.columns)

    