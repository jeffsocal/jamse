<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function spec_get_peaks_by_int($table, $pre_mz, $npeaks = 40, $difpeaks = 17.5)
{
    $spec = $table;
    $spec_iso = [];
    $spec_rm = [];
    
    /*
     * SCBI cuts peaks below mz 224 (for indexing reasons)
     */
    $func = function ($x) {
        return $x <= 224;
    };
    $spec_rm = array_merge($spec_rm, array_keys(array_filter($spec['mz'], $func)));
    
    /*
     * OMSSA cuts out precursor peaks
     * OMSSA deleting peaks that are within 2 Da of precuror
     */
    $spec_len = count($spec['int']);
    $a_pmz = array_fill(1, $spec_len, $pre_mz);
    $spec_rm = array_merge($spec_rm, array_keys(array_filter(array_map("dropByPreMz", $spec['int'], $a_pmz))));
    $spec = table_reindex(table_droprows($spec, array_unique($spec_rm)));
    
    // /*
    // * OMSSA cuts peaks below 2.5% of the maximum intensity
    // */
    // $spec_rm = [];
    // $spec_len = count($spec['int']);
    // $a_int = array_fill(1, $spec_len, array_max($spec['int']) * .005);
    // $spec_rm = array_merge($spec_rm, array_keys(array_filter(array_map("dropByInt", $spec['int'], $a_int))));
    
    // $spec = table_reindex(table_sort(table_droprows($spec, array_unique($spec_rm)), 'int', 'desc', 'number'));
    
    $spec = table_sort($spec, 'int', 'desc', 'number');
    $spec_keep = [];
    $spec_rm = [];
    
    /*
     * OMSSA deleting peaks that are within 2 Da of mz ranked by int
     * OMSSA deleting peaks that are within 27 Da of mz ranked by int
     */
    $spec_a = $spec_b = $spec['mz'];
    foreach ($spec_a as $ia => $a_m) {
        
        if (key_exists($ia, $spec_rm))
            continue;
        
        $spec_keep[$ia] = $ia;
        
        foreach ($spec_b as $ib => $b_m) {
            
            if (key_exists($ib, $spec_keep))
                continue;
            
            if (abs($a_m - $b_m) < $difpeaks)
                $spec_rm[$ib] = $ib;
        }
    }
    $spec = table_sort(table_reindex(table_droprows($spec, array_unique($spec_rm))), 'int', 'desc', 'number');
    
    $spec_a = $spec_b = $spec['mz'];
    foreach ($spec_a as $ia => $a_m) {
        
        if (key_exists($ia, $spec_rm))
            continue;
        
        foreach ($spec_b as $ib => $b_m) {
            
            if (abs($a_m - $b_m - 18) < 1)
                $spec_rm[$ib] = $ib;
        }
    }
    $spec = table_sort(table_reindex(table_droprows($spec, array_unique($spec_rm))), 'int', 'desc', 'number');
    
    /*
     * take the top 20 peaks
     */
    $spec = table_head($spec, $npeaks);
    
    $spec = table_sort(table_reindex(table_droprows($spec, array_unique($spec_rm))), 'mz', 'ascn', 'number');
    
    return $spec['mz'];
}

?>