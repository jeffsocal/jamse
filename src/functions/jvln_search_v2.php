<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use function BenTools\CartesianProduct\cartesian_product;
use Jvln\Engines\TandemHash;

function jvln_search_v2($pre_m, $peaks, $n = 4)
{
    $frg = new TandemHash();
    $hash = array_values(array_map(array(
        $frg,
        'massToHash'
    ), ($peaks)));
    
    $zs = range(1, 5);
    $as = range(0, - 1);
    // $as = [
    // 0
    // ];
    $p_hs = [];
    foreach ($zs as $z) {
        foreach ($as as $a) {
            $p_mz = round(neutralMass($pre_m + $a / $z, $z), 3);
            $p_hs[] = "pmv" . $frg->massToHash($p_mz, false);
        }
    }
    
    $p_hs = "(" . array_tostring($p_hs, " | ", "") . ")";
    
    $must = ' ' . array_tostring(array_slice($hash, 0, $n), ' ', '');
    $mayb = ' "' . array_tostring(array_slice($hash, $n), ' ', '') . '"/' . $n;
    
    return $p_hs . $must . $mayb;
}

?>