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

function getSphinxQueryV1($pre_m, $peaks, $n = 4, $k = 7)
{
    $frg = new TandemHash();
    $hash = array_values(array_map(array(
        $frg,
        'massToHash'
    ), ($peaks)));
    
    $zs = range(1, 3);
    $as = range(0, - 2);
    $as = [
        0
    ];
    $p_hs = [];
    foreach ($zs as $z) {
        foreach ($as as $a) {
            $p_mz = round(neutralMass($pre_m + $a, $z), 3);
            $p_hs[] = "pmv" . $frg->massToHash($p_mz);
        }
    }
    
    $p_hs = "(" . array_tostring($p_hs, " | ", "") . ") ";
    
    $must = array_slice($hash, 0, $k);
    $mayb = array_slice($hash, $k);
    $mayb = 'MAYBE ' . array_tostring($mayb, ' MAYBE ', '');
    
    $queries = [];
    if ($n == 0 and $k == 0) {
        
        $queries[] = $p_hs . $mayb;
        return $queries;
    }
    
    $cases = array_fill(1, $n, array_keys($must));
    $cp = cartesian_product($cases);
    
    $kcn = [];
    foreach ($cp->asArray() as $p) {
        if (count(array_unique($p)) < $n)
            continue;
        
        sort($p);
        $kcn[array_tostring($p, '', '')] = $p;
    }
    
    sort($must);
    foreach ($kcn as $p) {
        $must_q = '';
        foreach ($must as $n => $v) {
            if (! key_exists($n, array_flip($p)))
                $must_q .= 'MAYBE ';
            $must_q .= $v . ' ';
        }
        $queries[] = $p_hs . $must_q . $mayb;
    }
    
    return $queries;
}

?>