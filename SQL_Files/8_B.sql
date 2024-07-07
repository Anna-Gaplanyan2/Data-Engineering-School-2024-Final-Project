CREATE VIEW molecule_similarity_with_top_targets AS
WITH ranked_similarities AS (
    SELECT
        source_chembl_id,
        target_chembl_id,
        similarity_score,
        ROW_NUMBER() OVER (PARTITION BY source_chembl_id ORDER BY similarity_score DESC) AS rn
    FROM
        molecule_similarities
)
SELECT
    r1.source_chembl_id,
    r1.target_chembl_id,
    r1.similarity_score,
    r2.target_chembl_id AS next_most_similar_target_chembl_id,
    r3.target_chembl_id AS second_most_similar_target_chembl_id
FROM
    ranked_similarities r1
LEFT JOIN
    ranked_similarities r2 ON r1.source_chembl_id = r2.source_chembl_id AND r2.rn = 2
LEFT JOIN
    ranked_similarities r3 ON r1.source_chembl_id = r3.source_chembl_id AND r3.rn = 3
WHERE
    r1.rn = 1;
