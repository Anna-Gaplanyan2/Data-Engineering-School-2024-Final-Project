# Data Engineering School 2024 Final Project
## Overview
The objective of this project is to identify the ten most similar molecules for each set of input molecules. The project includes obtaining data from ChemBL, managing it to generate molecular fingerprints and determining similarity scores.
## Prerequisites
- Python 3.x
- Required Python packages (`requirements.txt`)
- AWS S3 Bucket access
- PostgreSQL Database
## Usage
### Ingest ChemBL Data
- python Task1_Chembl_lookup.py
- python Task1_Molecule.py
- python Chembl_lookup_Molecule_dictionary_insertion.py
- python Properties_Structures_insertion.py
- python New_files_insertion.py
### Compute Morgan Fingerprints
- python Task2_Morgan_fingerprints.py
### Compute Tanimoto Similarity Scores
- python Task3_Tanimoto_similarity_scores.py
### Full Similarity Score Table
- python Task4_Full_similarity_score_table.py
### Extract Top-10 Similar Molecules
- python Task5.py
### Create Data Marts
- python Task6_Dimension_Table.py
- python Task6_Fact_Table.py
### Database Views
- 7_A.sql
- 8_A.sql
- 8_B.sql
- 8_C.sql
## Results
- **Example 1:**
  ```json
  {
    "source_molecule": "CHEMBL1081575",
    "top_10_similar_molecules": [
      {"chembl_id": "CHEMBL1081575", "similarity_score": 1.00},
      {"chembl_id": "CHEMBL1079438", "similarity_score": 0.85},
      {"chembl_id": "CHEMBL1079439", "similarity_score": 0.82},
      {"chembl_id": "CHEMBL1076160", "similarity_score": 0.82},
      {"chembl_id": "CHEMBL1080259", "similarity_score": 0.77},
      {"chembl_id": "CHEMBL1075643", "similarity_score": 0.73},
      {"chembl_id": "CHEMBL1076158", "similarity_score": 0.72},
      {"chembl_id": "CHEMBL1081931", "similarity_score": 0.65},
      {"chembl_id": "CHEMBL1075668", "similarity_score": 0.65},
      {"chembl_id": "CHEMBL3425597", "similarity_score": 0.64},  
     
    ]
  }
## Author
Anna Gaplanyan
