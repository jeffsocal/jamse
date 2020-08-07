<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use function BenTools\CartesianProduct\cartesian_product;

function averagineMass()
{
    return 114.2446;
}

function averagineBins($mz, $z)
{
    return ceil(neutralMass($mz, $z) / averagineMass() * 2);
}

function neutralMass($mz, $z)
{
    if ($z == 0)
        $z = 1;
    
    return ($mz - 1.00727646688) * $z;
}

function chargeMass($nm, $z)
{
    if ($z == 0)
        return $nm;
    
    return ($nm / $z) + 1.00727646688;
}

function is_wholenumber($num)
{
    if (strpos($num, '.') !== false)
        return FALSE;
    
    return TRUE;
}

function getIsoCluster($arr, $mz)
{
    $this_rm = [];
    foreach ($arr['mz'] as $key => $val) {
        if ($val < $mz | $val > ($mz + 4.1))
            $this_rm[] = $key;
    }
    
    return table_reindex(table_droprows($arr, array_unique($this_rm)));
}

function getIsoPeaks($arr, $mz, $z)
{
    $wtf = array(
        'a' => $arr['mz'],
        'b' => $arr['mz']
    );
    
    $cp = cartesian_product($wtf);
    
    $this_rm = [];
    foreach ($cp->asArray() as $key => $val) {
        if (abs(abs($val['a'] - $val['b']) - (1 / $z)) < 0.05)
            $this_rm = array_unique(array_merge($this_rm, array_values($val)));
    }
    
    if (count($this_rm) == 0)
        return false;
    
    $this_iso = array(
        'a' => [
            $mz
        ],
        'b' => $this_rm
    );
    $cp = cartesian_product($this_iso);
    $this_rm = [];
    foreach ($cp->asArray() as $key => $val) {
        if (is_true(is_wholenumber(round(abs(abs($val['b'] - $val['a'])), 1) * ($z))))
            $this_rm = array_unique(array_merge($this_rm, array_values($val)));
    }
    
    if (count($this_rm) == 0)
        return false;
    
    sort($this_rm);
    
    return $this_rm;
}

?>