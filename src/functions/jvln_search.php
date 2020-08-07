<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use Jvln\Engines\TandemHash;

function jvln_search(array $peaks, float $premass = NULL, array $precharge = NULL, array $preisotope = NULL, int $topn = 4, float $quorum = 0.5)
{
    $frg = new TandemHash();
    $hash = array_values(array_map(array(
        $frg,
        'massToHash'
    ), ($peaks)));
    
    $p_hs = '';
    if ($premass != NULL) {
        
        if ($precharge == NULL)
            $precharge = [
                1
            ];
        
        if ($preisotope == NULL)
            $preisotope = [
                0
            ];
        
        $p_hs = [];
        foreach ($precharge as $z) {
            foreach ($preisotope as $a) {
                $p_mz = round(neutralMass($premass + $a / $z, $z), 3);
                $hash_bin = $frg->precursorMassSpread($p_mz);
                foreach ($hash_bin as $hash_mass) {
                    $p_hs[] = $frg->precursorMassToHash($hash_mass);
                }
            }
        }
        
        $p_hs = "(" . array_tostring($p_hs, " | ", "") . ")";
    }
    $must = ' ' . array_tostring(array_slice($hash, 0, $topn), ' ', '');
    sort($hash);
    $mayb = ' "' . array_tostring($hash, ' ', '') . '"/' . $quorum;
    
    return $p_hs . $must . $mayb;
}

?>