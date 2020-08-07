<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function json_scans($json)
{
    $row_data = [
        'file',
        'title',
        'scan_id',
        'rt_sec',
        'pre_mz',
        'pre_z',
        'int_rsd',
        'int_rad',
        'int_iqr',
        'int_ntp',
        'a_bins',
        't_peaks',
        'n_peaks',
        'pmatch',
        'ematch'
    ];

    $scan_data = $json['scans'];
    $file_name = basename($json['file']['name']);

    $data_out = [];
    foreach ($scan_data as $i => $scan_values) {
        $data_out[$i] = [];
        $scan_values['file'] = $file_name;
        $scan_values['scan_id'] = $i;

        foreach ($row_data as $var) {
            $data_out[$i][$var] = $scan_values[$var];
        }
    }

    return table_invert($data_out);
}
?>