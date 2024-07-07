CREATE OR REPLACE VIEW avg_similarity_scores AS
SELECT
  CASE
    WHEN GROUPING(source_chembl_id) = 1 AND GROUPING(aromatic_rings) = 1 AND GROUPING(heavy_atoms) = 1 THEN 'TOTAL'
    ELSE source_chembl_id
  END AS source_chembl_id,
  aromatic_rings,
  heavy_atoms,
  AVG(similarity_score) AS avg_similarity_score
FROM
  molecule_similarities
JOIN
  molecule_dimension ON molecule_similarities.source_chembl_id = molecule_dimension.chembl_id
GROUP BY
  GROUPING SETS (
    (source_chembl_id),
    (aromatic_rings, heavy_atoms),
    (heavy_atoms),
    ()
  );
