CREATE VIEW `experiments` AS
    SELECT 
        experiment_name,
        count(*) n_files,
        sum(scans_n) n_scans
    FROM
        jvln.experiments,
        jvln.files
    WHERE
        jvln.experiments.file_pk = jvln.files.pk
        group by experiment_name