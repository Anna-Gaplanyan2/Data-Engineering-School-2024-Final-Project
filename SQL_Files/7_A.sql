CREATE VIEW avg_similarity_score_per_source_id AS
SELECT
    source_chembl_id,
    AVG(similarity_score) AS avg_similarity_score
FROM
    molecule_similarities
GROUP BY
    source_chembl_id;
