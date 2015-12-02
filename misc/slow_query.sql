SELECT * FROM 
(
	SELECT 
		nms.*,
		structure.transcript_id AS structure_transcript_id
		FROM (
			SELECT nucleotide_measurement_set.*
			FROM nucleotide_measurement_set 
			ORDER BY nucleotide_measurement_set.coverage DESC
			LIMIT 20 OFFSET 0
		) AS nms LEFT OUTER JOIN structure ON 
			structure.transcript_id = nms.transcript_id AND
			structure.structure_prediction_run_id = 2 
) AS jnms, 
	transcript,
	gene_location

WHERE
	jnms.nucleotide_measurement_run_id = 1 AND 
	transcript.id = jnms.transcript_id AND 
	transcript.gene_id = gene_location.gene_id AND 
	gene_location.strain_id = 'Col_0' 

GROUP BY jnms.transcript_id



SELECT 
	nucleotide_measurement_set.*, 
	transcript.*
	-- gene_location.*
	-- structure.*
FROM transcript, gene_location, nucleotide_measurement_set 
-- LEFT OUTER JOIN structure ON 
-- 	structure.transcript_id = nucleotide_measurement_set.transcript_id AND
-- 	structure.structure_prediction_run_id = 2 
WHERE 
	nucleotide_measurement_set.nucleotide_measurement_run_id = 1 AND 
	transcript.id = nucleotide_measurement_set.transcript_id AND 
	transcript.gene_id = gene_location.gene_id 
	-- gene_location.strain_id = 'Col_0' 

GROUP BY transcript.id 
ORDER BY nucleotide_measurement_set.coverage DESC
LIMIT 2 OFFSET 0 # 1.8702




SELECT 
	nucleotide_measurement_set.*, 
	transcript.*, 
	gene_location.*, 
	structure.*
FROM transcript, gene_location, nucleotide_measurement_set 
LEFT OUTER JOIN structure ON 
	structure.transcript_id = nucleotide_measurement_set.transcript_id AND
	structure.structure_prediction_run_id = 2 
WHERE 
	nucleotide_measurement_set.nucleotide_measurement_run_id = 1 AND 
	transcript.id = nucleotide_measurement_set.transcript_id AND 
	transcript.gene_id = gene_location.gene_id AND 
	gene_location.strain_id = 'Col_0' 

GROUP BY transcript.id 
ORDER BY nucleotide_measurement_set.coverage DESC
LIMIT 1 OFFSET 0
