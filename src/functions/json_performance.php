<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function json_performance($json)
{
    if (! key_exists('performance', $json))
        return FALSE;
    
        $row_data = [
            'scan_id',
            'psm_eprob',
            'psm_prob',
            'psm_qvalue'
        ];

    $perf_data = $json['performance']['psms'];
    $data_out = [];
    foreach ($perf_data as $i => $scan_values) {
        $data_out[$i] = [];
        $scan_values['scan_id'] = $i;

        foreach ($row_data as $var) {
            $data_out[$i][$var] = $scan_values[$var];
        }
    }

    return table_invert($data_out);
}
?>