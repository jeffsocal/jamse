<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
function json_proteins($json)
{
    if (! key_exists('proteins', $json))
        return FALSE;

    $data = $json['proteins'];
    $file_name = basename($json['file']['name']);

    foreach ($data as $i => $values) {
        unset($data[$i]['peptides']);
        $data[$i]['file'] = $file_name;
    }

    return table_invert($data);
}

function json_proteingroups($json)
{
    if (! key_exists('protein_groups', $json))
        return FALSE;

    $data = $json['protein_groups'];

    return $data;
}
?>