<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function json_peptides(array $json, string $form = 'table')
{
    if (! key_exists('jaevalen', $json))
        return FALSE;

    if ($form == 'json')
        return $json['jaevalen'];

    $jvln_stat = $json['jaevalen']['stats'];
    $jvln_data = $json['jaevalen']['psms'];

    $fill_out = [
        'scan_id',
        'hit_rank',
        'time_sec',
        'psm_n',
        'peptide',
        'psm_score',
        'pre_mma',
        'psm_z',
        'psm_i',
        'psm_o',
        'psm_r',
        'psm_mma',
        'psm_pvalue',
        'psm_evalue'
    ];

    $fill_out = array_combine($fill_out, array_fill(0, count($fill_out), ''));

    $data_out = [];
    foreach ($jvln_stat as $i => $v) {

        $data_out[] = $fill_out;
        $n = count($data_out) - 1;

        $data_out[$n]['scan_id'] = $i;
        $data_out[$n]['psm_n'] = $v['psm_n'];
        $data_out[$n]['time_sec'] = $v['time_sec'];
        
        $this_set = $data_out[$n];

        if (! key_exists($i, $jvln_data))
            continue;

        foreach ($jvln_data[$i] as $hit_rank => $meta) {

            if (key_exists('hit_rank', $meta)) {
                if ($meta['hit_rank'] != 1)
                    continue;
                $data_out[$n] = array_merge($this_set, $meta);
            } else {
                $n += $hit_rank;
                $data_out[$n] = array_merge($this_set, $meta);
            }
        }
    }

    return table_invert($data_out);
}
?>