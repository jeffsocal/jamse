<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Database;

use ProxyMySQL\Transactions\Prepared;
use ProxyMySQL\Transactions\Simple;

class Methods
{

    protected $Prepared;

    protected $Simple;

    public function __construct(string $server, int $port = 3306)
    {
        $this->Prepared = new Prepared($server, $port);
        $this->Prepared->setSchema('jvln');

        $this->Simple = new Simple($server, $port);
        $this->Simple->setSchema('jvln');
    }

    public function putMethods(array $method)
    {
        /*
         * place the methods into the db
         */
        $this->Prepared->setStatement("INSERT INTO jvln.methods
                             (parameters)
                             VALUES (?)");
        $this->Prepared->addVarible('parameters', json_encode($method), 's');
        $this->Prepared->paramPut();
        return $this->getLastInsertID();
    }

    public function getLastInsertID()
    {
        return $this->Prepared->getLastInsertID();
    }

    public function getFileIDPK(string $file_id)
    {
        $query = "SELECT pk from jvln.files WHERE file_id = ?;";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;

        return $data['pk'][0];
    }

    public function listExperiments()
    {
        $query = "SELECT *  
                    FROM jvln.summary_experiments;";

        $data = $this->Simple->sqlGet($query);

        if (is_false($data))
            return FALSE;

        return $data;
    }

    public function getFilesByMatch($str)
    {
        $query = "SELECT 
                        experiment_name,
                        file_name,
                        file_id,
                        scans_n,
                        CONCAT_WS(' ', experiment_name,
                                file_name,
                                file_id,
                                scans_n) as text
                    FROM
                        jvln.experiments,
                        jvln.files
                    WHERE
                        jvln.experiments.file_pk = jvln.files.pk
                    HAVING
                        text like ?";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('text', "%" . $str . "%");

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;

        return table_dropcols($data, 'text');
    }

    public function getFilesByExperiment($str)
    {
        $query = "SELECT 
                        experiment_name, file_name, file_id, scans_n
                    FROM
                        jvln.experiments,
                        jvln.files
                    WHERE
                        jvln.experiments.file_pk = jvln.files.pk
                        AND experiment_name like ?;";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('experiment_name', "%" . $str . "%");

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;

        return $data;
    }

    public function fileScansExist($file_name)
    {
        $query = "SELECT file_name, file_id, count(*) count
                    FROM
                        jvln.files,
                        jvln.scans
                    WHERE
                        jvln.scans.file_pk = jvln.files.pk
                        AND file_name = ?
                    GROUP BY file_name, file_id;";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_name', $file_name);

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;
        else
            return table_invert($data)[0];
    }

    public function filePsmsExist($file_id)
    {
        $query = "SELECT file_name, file_id, count(distinct(scan_pk)) count
                    FROM
                        jvln.files,
                        jvln.scans,
                        jvln.psms
                    WHERE
                        jvln.scans.file_pk = jvln.files.pk AND
                        jvln.scans.pk = jvln.psms.scan_pk
                        AND file_id = ?
                    GROUP BY file_name, file_id;";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;
        else
            return table_invert($data)[0];
    }

    public function fileFDRExist($file_id, $fdr_metric = "psm_fcorr")
    {
        $query = "SELECT file_name, file_id, count(distinct(scan_pk)) count
                    FROM
                        jvln.files, jvln.scans,
                        jvln.psms, jvln.validations
                    WHERE
                        jvln.scans.file_pk = jvln.files.pk AND
                        jvln.scans.pk = jvln.psms.scan_pk AND
                        jvln.psms.pk = jvln.validations.psm_pk AND
                        file_id = ? AND
                        fdr_metric = ?
                    GROUP BY file_name, file_id;";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);
        $this->Prepared->addVarible('fdr_metric', $fdr_metric);

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;
        else
            return table_invert($data)[0];
    }

    public function fileMetricsExist($file_id)
    {
        $query = "SELECT file_name, file_id, count(distinct(metrics.pk)) count
                    FROM
                        jvln.files, jvln.metrics
                    WHERE
                        jvln.metrics.file_pk = jvln.files.pk 
                        AND file_id = ?
                    GROUP BY file_name, file_id;";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);

        $data = $this->Prepared->paramGet();

        if (is_false($data))
            return FALSE;
        else
            return table_invert($data)[0];
    }
}

?>