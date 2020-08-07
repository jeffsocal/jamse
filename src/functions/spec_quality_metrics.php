<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function spec_quality_metrics($table, $pre_mz, $pre_z = 2, $difpeaks = 17.5)
{
    $spec = $table;
    
    $int_sum = array_sum($spec['int']);
    $int_max = array_max($spec['int']);
    $int_mean = array_mean($spec['int']);
    $int_median = array_median($spec['int']);
    $int_sdev = array_stdev($spec['int']);
    $int_adev = array_absdev($spec['int']);
    
    $out['int_sum_log'] = round(log($int_sum), 2);
    $out['int_max_log'] = round(log($int_max), 2);
    $out['int_mean_log'] = round(log($int_mean), 3);
    $out['int_median_log'] = round(log($int_median), 3);
    $out['int_sdev_log'] = round(log($int_sdev), 4);
    $out['int_adev_log'] = round(log($int_adev), 4);
    
    /*
     * cut out precursor peaks
     * deleting peaks that are within 2 Da of precuror
     */
    
    $spec = table_sort($spec, 'int', 'desc', 'number');
    $spec['peak'] = range(0, table_length($spec) - 1, 1);
    $spec = $ref = array_rowtocol($spec);
    
    foreach ($spec as $i => $v) {
        if (abs($v['mz'] - $pre_mz) <= 2)
            unset($spec[$i]);
    }
    $spec = array_values($spec);
    
    /*
     * iterate through array removing peaks around
     * the tallest peaks by intensity
     */
    $toss = [];
    $keep = [];
    while (count($keep) < 50) {
        
        if (count($spec) == 0)
            break;
        
        $mz = $spec[0]['mz'];
        $in = $spec[0]['int'];
        $pk = $spec[0]['peak'];
        $keep[] = $pk;
        
        unset($spec[0]);
        
        foreach ($spec as $i => $v) {
            $diff = $mz - $v['mz'];
            if (abs($diff) <= $difpeaks) {
                unset($spec[$i]);
                if (abs($diff - 1) <= 1.05 and $v['int'] > $in / 2) {
                    
                    $toss[] = $pk;
                    $keep[] = $v['peak'];
                }
            }
        }
        $spec = array_values($spec);
    }
    
    $keep = array_diff($keep, $toss);
    
    $keep = array_intersect_key($ref, array_flip($keep));
    
    if (count($keep) == 0)
        return array_values($table['mz']);
    
    $bins = array();
    foreach ($keep as $i => $vals) {
        
        if ($vals['int'] < $int_median)
            continue;
        
        $bin = ceil($vals['mz'] / 114);
        
        $bins[] = $bin;
    }
    
    $out['avg_bin_occupation'] = count(array_unique($bins)) / ceil(neutralMass($pre_mz, $pre_z) / 114);
    
    return $out;
   
}

?>