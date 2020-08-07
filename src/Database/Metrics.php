<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Database;

class Metrics extends Methods
{

    private $Spectra;

    private $PeptideMatches;

    public function __construct(string $server, int $port = 3306)
    {
        parent::__construct($server, $port);

        $this->Spectra = new Spectra($server, $port);
        $this->PeptideMatches = new PeptideMatches($server, $port);
    }

    function getScanMetrics(string $file_id, int $pmatch = 0)
    {
        $query = "SELECT
                      scans.pk scan_pk, file_id, title,
                      rt_sec, pre_mz, pre_z, int_rsd, int_rad,
                      int_iqr, int_ntp, t_peaks, n_peaks,
                      spec_prob, spec_fcorr
                  FROM
                      jvln.scans,
                      jvln.files
                  WHERE
                      jvln.scans.file_pk = jvln.files.pk
                      AND file_id = ?";

        if ($pmatch > 0)
            $query .= " pmatch <= ?";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);

        if ($pmatch > 0)
            $this->Prepared->addVarible('pmatch', $pmatch);

        return $this->Prepared->paramGet();
    }

    function getPsmMetrics(string $file_id, float $qvalue = 0.95)
    {
        $query = "SELECT
                  pre_mma, psm_z, psm_i, psm_o, psm_mma
                  FROM jvln.psms, jvln.scans, jvln.files
                  WHERE jvln.scans.file_pk = jvln.files.pk
                  AND jvln.psms.scan_pk = jvln.scans.pk
                  AND hit_rank = 1
                  AND psm_prob >= ?
                  AND file_id = ?";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('psm_qvalue', $qvalue);
        $this->Prepared->addVarible('file_id', $file_id);

        return $this->Prepared->paramGet();
    }

    function updatePsms(array $table)
    {
        $this->Simple->setTable('psms');
        return $this->Simple->sqlMod($table);
    }

    function putMetricData(array $table, string $file_pk)
    {
        $table['file_pk'] = array_fill(0, table_length($table), $file_pk);
        $this->Simple->setTable('metrics');
        return $this->Simple->sqlPut($table);
    }
}

?>