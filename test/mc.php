<?php

use ProxySci\Omics\Proteomics\Fragments;

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
$frg = new Fragments();

ini_set("memory_limit", "4096M");

$m = opts('m');
$zs = explode(",", opts('z'));

if ($zs[0] == '')
    $zs = [
        0,
        1,
        2,
        3
    ];

print_message("mass", $m);

$table = [];
foreach ($zs as $z) {
    
    $table['z'][] = $z;
    $table['n mass'][] = neutralMass($m, $z);
    $table['z mass'][] = chargeMass($m, $z);
    
    $table['n hash'][] = $frg->massToHash(neutralMass($m, $z));
    $table['z hash'][] = $frg->massToHash(chargeMass($m, $z));
}

print_table($table);
?>