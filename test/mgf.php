<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use ProxySci\MassSpec\Fragmentation\ReadMGF;
use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Fragments;
use function BenTools\CartesianProduct\cartesian_product;
use ProxyIO\File\Delim\WriteDelim;

$tmr = new Timer();
// $mgf = new ReadMGF("/data/SCBI/20171201_DBS_20Samples/mgf/20171201_SoCal_Samples18-9.mgf");
$mgf = new ReadMGF("/data/OpenAccess/MSdata/yeast_gs/data/mgf/01_061220_000.mgf");
$frg = new Fragments();
$csv = new WriteDelim();

$seq = opts('s');

echo PHP_EOL;

$data = $mgf->getData();

// $gets = array_map(function ($x) {
// return $x + 4486;
// }, range(1, 20, 1));

// print_r($gets);

// $gets = range(0,count($data),1);

$gets = [
    932,
    2330
];

// function neutralMass($mz, $z)
// {
//     return $mz * $z - $z * 1.00727646688;
// }

// function chargeMass($mz, $z)
// {
//     return neutralMass($mz, $z) + 1.00727646688;
// }

// function is_wholenumber($num)
// {
//     if (strpos($num, '.') !== false)
//         return FALSE;
    
//     return TRUE;
// }

// function getIsoCluster($arr, $mz)
// {
//     $this_rm = [];
//     foreach ($arr['mz'] as $key => $val) {
//         if ($val < $mz | $val > ($mz + 4.1))
//             $this_rm[] = $key;
//     }
    
//     return table_reindex(table_droprows($arr, array_unique($this_rm)));
// }

// function getIsoPeaks($arr, $mz, $z)
// {
//     $wtf = array(
//         'a' => $arr['mz'],
//         'b' => $arr['mz']
//     );
    
//     $cp = cartesian_product($wtf);
    
//     $this_rm = [];
//     foreach ($cp->asArray() as $key => $val) {
//         if (abs(abs($val['a'] - $val['b']) - (1 / $z)) < 0.05)
//             $this_rm = array_unique(array_merge($this_rm, array_values($val)));
//     }
    
//     if (count($this_rm) == 0)
//         return false;
    
//     $this_iso = array(
//         'a' => [
//             $mz
//         ],
//         'b' => $this_rm
//     );
//     $cp = cartesian_product($this_iso);
//     $this_rm = [];
//     foreach ($cp->asArray() as $key => $val) {
//         if (is_true(is_wholenumber(round(abs(abs($val['b'] - $val['a'])), 1) * ($z))))
//             $this_rm = array_unique(array_merge($this_rm, array_values($val)));
//     }
    
//     if (count($this_rm) == 0)
//         return false;
    
//     sort($this_rm);
    
//     return $this_rm;
// }

foreach ($gets as $get) {
    
    $this_data = $data[$get];
    
    $spec = $this_data['spec'];
    $prec = preg_replace("/\s+.*/", "", $this_data['pepmass']);
    
    // print_table($spec, 300);
    
    print_message('scan', $get);
    print_message('  title', $this_data['title']);
    print_message('  rtinseconds', $this_data['rtinseconds']);
    print_message('  max int', max($spec['int']));
    print_message("precursor mass", $prec);
    print_message(" z1", round(neutralMass($prec, 1), 3));
    print_message(" z2", round(neutralMass($prec, 2), 3));
    print_message(" z3", round(neutralMass($prec, 3), 3));
    
    $spec_int_desc = table_sort($spec, 'int', 'desc', 'number');
    
    $spec_iso = [];
    $spec_rm = [];
    
    /*
     * determine the isotope cluster and charge state of fragment ions
     * iterate over several charge states
     */
    // for ($z = 1; $z <= 2; $z ++) {
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
            else
                echo $mz . PHP_EOL;
        }
    }
    
    if (key_exists('z2', $spec_iso)) {
        foreach (array_keys($spec_iso['z2']) as $n => $mz) {
            
            if ($key = array_search($mz, $spec['mz']))
                $spec['z'][$key] = 2;
            else
                echo $mz . PHP_EOL;
        }
    }
    
    // print_r($spec_iso);
    // exit;
    /*
     * inspect the spectrum here
     */
    $csv->writeCsvFile("/data/tmp_spec.csv", $spec);
    // echo PHP_EOL;
    // exit();
    
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
    
    // echo $prec . "\t" . $max_int . "\t" . $max_int * 0.025 . PHP_EOL;
    /*
     * OMSSA cuts peaks below 2.5% of the maximum intensity
     */
    foreach ($spec['int'] as $n => $v) {
        if ($v <= $max_int * 0.025)
            $spec_rm[] = $n;
    }
    
    $spec['i'] = range(0,count($spec['int'])-1);
    
    // eval('$func = function ($x) { return $x <= ' . $max_int * 0.025 . '; };');
    // $spec_rm = array_merge($spec_rm, array_keys(array_filter($spec['int'], $func)));
    // eval('$func = function ($x) { return $x == ' . $min_int . '; };');
    // $spec_rm = array_merge($spec_rm, array_keys(array_filter($spec['int'], $func)));
    
    // /*
    // * OMSSA cuts out precursor peaks
    // */
    // eval('$func = function ($x) { return ($x > ' . ($prec - 2) . ' & $x < ' . ($prec + 2) . '); };');
    // $spec_rm = array_merge($spec_rm, array_keys(array_filter($spec['mz'], $func)));
    
    // /*
    // * OMSSA deleting peaks that are within 2 Da of precuror
    // */
    // eval('$func = function ($x) { return ($x > ' . ($prec - 2) . ' & $x < ' . ($prec + 2) . '); };');
    // $spec_rm = array_merge($spec_rm, array_keys(array_filter($spec['mz'], $func)));
    
    
    $spec = table_reindex(table_droprows($spec, array_unique($spec_rm)));
//     print_table($spec, 300);
    
    
    $spec = table_sort($spec, 'int', 'desc', 'number');
    $spec = table_head($spec, 20);
    foreach ($spec['mz'] as $n => $mz) {
        $spec['nm'][$n] = chargeMass(neutralMass($spec['mz'][$n], $spec['z'][$n]), 1);
    }

    print_table($spec, 20);
    
    // /*
    // * OMSSA deleting peaks that are within 2 Da of mz ranked by mz
    // * OMSSA deleting peaks that are within 27 Da of mz ranked by mz
    // */
    // $spec_rm = [];
    // foreach ($spec['mz'] as $n => $mz) {
    // if (array_key_exists($n, array_flip($spec_rm)))
    // continue;
    
    // eval('$func = function ($x) { return ($x > ' . ($mz - 27) . ' & $x < ' . ($mz + 27) . '); };');
    // $keys = array_diff(array_keys(array_filter($spec['mz'], $func)), array(
    // $n
    // ));
    // $spec_rm = array_unique(array_merge($spec_rm, $keys));
    // }
    
    // $spec = table_reindex(table_droprows($spec, array_unique($spec_rm)));
    // $spec = table_sort($spec, 'mz', 'ascn', 'number');
    
    $spec['hash'] = array_map(array(
        $frg,
        'massToHash'
    ), ($spec['nm']));
    
    // $spec = table_sort($spec, 'mz', 'ascn', 'number');
    
    
    $zs = range(1, 3);
    $as = range(0, - 2);
    $as = [0];
    $p_hs = [];
//     foreach ($zs as $z) {
//         foreach ($as as $a) {
//             $p_mz = round(neutralMass(($prec + $a / $z), $z), 3);
//             $p_hs[] = "pmv" . $frg->massToHash($p_mz);
//         }
//     }
    
    foreach ($zs as $z) {
        foreach ($as as $a) {
            $p_mz = round(neutralMass($prec, $z), 3);
            $p_hs[] = "pmv" . $frg->massToHash($p_mz);
        }
    }
    
    $p_hs = "(" . array_tostring($p_hs, " | ", "") . ")";
    
    // echo "select *, WEIGHT() from i_fpep_ms WHERE MATCH('" . $p_hs . " " . array_tostring($spec['hash'], " ", "") . "') OPTION ranker=expr('bm25 - abs($p_mz - peptide_mass)/$p_mz * 1e6');";
    echo "select *, WEIGHT() from i_fpep_yeast WHERE MATCH('" . $p_hs . " " . array_tostring($spec['hash'], " MAYBE ", "") . "') OPTION ranker=expr('bm15*sum(hit_count)');";
    echo PHP_EOL . PHP_EOL;
    exit();
}

echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

?>