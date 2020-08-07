<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace Jvln\Database;

class PeptideMatches extends Methods
{

    public function __construct(string $server, int $port = 3306)
    {
        parent::__construct($server, $port);
    }

    function putPsmData(array $table, array $method)
    {
        $method_pk = $this->putMethods($method);
        if (is_false($method_pk))
            return FALSE;

        $table['method_pk'] = array_fill(0, table_length($table), $method_pk);

        $this->Simple->setTable('psms');
        return $this->Simple->sqlPut($table);
    }

    function putPsmAltData(array $table)
    {
        $this->Simple->setTable('psm_alts');
        return $this->Simple->sqlPut($table);
    }

    function getPsms(string $file_id, int $rank = 0)
    {
        $query = "SELECT psms.pk pk, 
                  scan_id, rt_sec, pre_mz, pre_z, 
                  int_rsd, int_rad, int_iqr, int_ntp, 
                  spec_fcorr, spec_prob, peptide,
                  pre_mma, psm_z, psm_i, psm_ovr_y, psm_ovr_b, 
                  psm_mma, 
                  psm_score, psm_fcorr, psm_fcorr_y, psm_fcorr_b, 
                  psm_prob, hit_rank
                  FROM jvln.psms, jvln.scans, jvln.files
                  WHERE jvln.scans.file_pk = jvln.files.pk 
                  AND jvln.psms.scan_pk = jvln.scans.pk 
                  AND file_id = ?";

        if ($rank > 0)
            $query .= " AND hit_rank <= ?";

        $this->Prepared->setStatement($query);
        $this->Prepared->addVarible('file_id', $file_id);

        if ($rank > 0)
            $this->Prepared->addVarible('hit_rank', $rank);

        return $this->Prepared->paramGet();
    }
}

?>