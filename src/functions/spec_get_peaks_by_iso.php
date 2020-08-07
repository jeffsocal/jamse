<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use function BenTools\CartesianProduct\cartesian_product;

function spec_get_peaks_by_iso($table, $npeaks = 40)
{
    $spec = $table;
    $spec_iso = [];
    $spec_rm = [];
    
    $spec_int_desc = table_sort($spec, 'int', 'desc', 'number');
    
    /*
     * determine the isotope cluster and charge state of fragment ions
     * iterate over several charge states
     */
    for ($z = 2; $z >= 1; $z --) {
        /*
         * run through each pairwise combination
         */
        for ($i = 0; $i < table_length($spec_int_desc); $i ++) {
            
            if (array_search($i, $spec_rm))
                continue;
            
            $this_mz = $spec_int_desc['mz'][$i];
            
            /*
             * pull out the peaks within a suspected isotopic cluster
             */
            $this_sp = getIsoCluster($spec_int_desc, $this_mz);
            
            /*
             * find the peaks that fit the current z profile
             */
            $this_iso = getIsoPeaks($this_sp, $this_mz, $z);
            
            if (is_false($this_iso))
                continue;
            
            $spec_rm = array_merge($spec_rm, array_keys(array_intersect($spec_int_desc['mz'], $this_iso)));
            
            $spec_iso['z' . $z][$this_iso[0]] = $this_iso;
        }
    }
    
    $spec['z'] = array_fill(0, count($spec['int']), '');
    
    if (key_exists('z1', $spec_iso)) {
        foreach (array_keys($spec_iso['z1']) as $n => $mz) {
            if ($key = array_search($mz, $spec['mz']))
                $spec['z'][$key] = 1;
        }
    }
    
    if (key_exists('z2', $spec_iso)) {
        foreach (array_keys($spec_iso['z2']) as $n => $mz) {
            if ($key = array_search($mz, $spec['mz']))
                $spec['z'][$key] = 2;
        }
    }
    
    $spec_rm = [];
    $max_int = array_max($spec['int']);
    $min_int = array_min($spec['int']);
    
    /*
     * SCBI cuts peaks below mz 224 (for indexing reasons)
     */
    $func = function ($x) {
        return $x <= 224;
    };
    $spec_rm = array_merge($spec_rm, array_keys(array_filter($spec['mz'], $func)));
    
    /*
     * keep just the mono-isotopes
     */
    $func = function ($x) {
        return $x == '';
    };
    $spec_rm = array_merge($spec_rm, array_keys(array_filter($spec['z'], $func)));
    
    /*
     * OMSSA cuts peaks below 2.5% of the maximum intensity
     */
    foreach ($spec['int'] as $n => $v) {
        if ($v <= $max_int * 0.025)
            $spec_rm[] = $n;
    }
    
    $spec['i'] = range(0, count($spec['int']) - 1);
    
    /*
     * OMSSA cuts out precursor peaks
     * OMSSA deleting peaks that are within 2 Da of precuror
     */
    
    $spec = table_reindex(table_droprows($spec, array_unique($spec_rm)));
    
    if (table_length($spec) < 3)
        return false;
    
    $spec = table_sort($spec, 'mz', 'ascn', 'number');
    $spec = table_head($spec, $npeaks);
    foreach ($spec['mz'] as $n => $mz) {
        $spec['nm'][$n] = chargeMass(neutralMass($spec['mz'][$n], $spec['z'][$n]), 1);
    }
    
    return $spec['nm'];
}

?>