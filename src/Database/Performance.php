<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Database;

class Performance extends Methods
{

    public function __construct(string $server = '127.0.0.1', int $port = 3306)
    {
        parent::__construct($server, $port);
    }

    function putPerfData(array $table, array $method)
    {
        $method_pk = $this->putMethods($method);
        if (is_false($method_pk))
            return FALSE;

        $table['method_pk'] = array_fill(0, table_length($table), $method_pk);

        $this->Simple->setTable('validations');
        return $this->Simple->sqlPut($table, 5000);
    }

    function getPsmPerf(string $file_id, string $metric = 'psm_fcorr', float $value = 0.05)
    {
        if (strstr("psm_fcorr|psm_score", $metric))
            return $this->getPsmFDR($file_id, $metric, $value);

        if (strstr("psm_prob|hit_prob", $metric))
            return $this->getPsmProb($file_id, $metric, $value);

        return FALSE;
    }

    private function getPsmFDR(string $file_id, string $metric = 'psm_fcorr', float $value = 0.05)
    {
        if (! strstr("psm_fcorr|psm_score", $metric))
            return FALSE;

        $query = "SELECT
                      file_id, scan_id,
                      psms.pk pk,    rt_sec,    pre_mz,
                      pre_z,    int_rsd,    int_rad,    int_iqr,
                      int_ntp,    spec_fcorr,   spec_prob,
                      peptide,	pre_mma,    psm_z,
                      psm_i,    psm_o,       psm_mma,
                      psm_score,   psm_fcorr, psm_prob, hit_prob,    hit_rank,
                      fdr_prob, fdr_eprob, fdr_qvalue, fdr_metric
                  FROM
                      jvln.psms,
                      jvln.scans,
                      jvln.files,
                      jvln.validations
                  WHERE
                      jvln.scans.file_pk = jvln.files.pk
                  	  AND jvln.psms.scan_pk = jvln.scans.pk
                      AND jvln.validations.psm_pk = jvln.psms.pk
                  AND file_id = ?
                  AND fdr_qvalue <= ?
                  AND fdr_metric = ?
                  ORDER BY fdr_qvalue";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);
        $this->Prepared->addVarible('fdr_qvalue', $value);
        $this->Prepared->addVarible('fdr_metric', $metric);

        return $this->Prepared->paramGet();
    }

    private function getPsmProb(string $file_id, string $metric = 'psm_prob', float $value = 0.95)
    {
        if (! strstr("psm_prob|hit_prob", $metric))
            return FALSE;

        $query = "SELECT
                      file_id, scan_id,    psms.pk pk,    rt_sec,    pre_mz,
                      pre_z,    int_rsd,    int_rad,    int_iqr,
                      int_ntp,   spec_fcorr,   spec_prob,
                      peptide,	pre_mma,    psm_z,
                      psm_i,    psm_o,       psm_mma,
                      psm_score,   psm_fcorr, psm_prob, hit_prob,    hit_rank
                  FROM
                      jvln.psms,
                      jvln.scans,
                      jvln.files
                  WHERE
                      jvln.scans.file_pk = jvln.files.pk
                  	  AND jvln.psms.scan_pk = jvln.scans.pk
                  AND file_id = ?
                  AND " . $metric . " >= ?
                  AND hit_rank = 1 ;";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);
        // $this->Prepared->addVarible('fdr_metric', $metric);
        $this->Prepared->addVarible('fdr_qvalue', $value);

        return $this->Prepared->paramGet();
    }
}

?>