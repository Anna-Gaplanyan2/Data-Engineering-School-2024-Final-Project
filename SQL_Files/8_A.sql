CREATE VIEW molecule_similarity_pivot AS
WITH random_sources AS (
    SELECT source_chembl_id
    FROM (
        SELECT DISTINCT source_chembl_id, ROW_NUMBER() OVER (ORDER BY RANDOM()) AS rn
        FROM molecule_similarities
    ) sub
    WHERE rn <= 10
),
pivot_data AS (
    SELECT
        target_chembl_id,
        source_chembl_id,
        similarity_score
    FROM
        molecule_similarities
    WHERE
        source_chembl_id IN (SELECT source_chembl_id FROM random_sources)
)
SELECT
    target_chembl_id,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 0) THEN similarity_score END) AS source_1,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 1) THEN similarity_score END) AS source_2,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 2) THEN similarity_score END) AS source_3,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 3) THEN similarity_score END) AS source_4,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 4) THEN similarity_score END) AS source_5,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 5) THEN similarity_score END) AS source_6,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 6) THEN similarity_score END) AS source_7,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 7) THEN similarity_score END) AS source_8,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 8) THEN similarity_score END) AS source_9,
    MAX(CASE WHEN source_chembl_id = (SELECT source_chembl_id FROM random_sources LIMIT 1 OFFSET 9) THEN similarity_score END) AS source_10
FROM
    pivot_data
GROUP BY
    target_chembl_id
ORDER BY
    target_chembl_id;
