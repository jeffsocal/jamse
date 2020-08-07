#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * jvlnmgf - jaevalen search tandem ms file
 *
 * SYNOPSIS
 * jvlnmgf -p <path> -f <regex>
 *
 * DESCRIPTION
 * searches the jaevalen engine for matching peptides based on the
 * precursor and fragments masses
 *
 *
 * COMMAND LINE
 *
 * -t : flag for testing, input the scan N to test
 *
 * -p : path to data file
 *
 * -f : file regex
 *
 * -i : require n number of top in search (default: '4,3,2,1,0')
 *
 * -n : max number of input peaks (default: 40)
 *
 * -d : deltas for intensity noise reduction(default: '17.5,27,35,70')
 *
 * -a : noise reduction algorithm (default: 'int,iso')
 *
 * -s : jscore minimum (default: 9000)
 *
 *
 * EXAMPLE
 *
 * jvlnmgf -p /home/scbi/data -f '^00[1-3].+\.mgf'
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxyIO\File\Delim\WriteDelim;
use ProxySci\MassSpec\Fragmentation\ReadMGF;
use ProxySci\Omics\Proteomics\Fragments;
use ProxyTime\ProgressTimer;
use ProxyTime\Timer;
use function BenTools\CartesianProduct\cartesian_product;
use Jvln\Engines\PIMS;

// $mgf = new ReadMGF("/data/SCBI/20171201_DBS_20Samples/mgf/20171201_SoCal_Samples18-9.mgf");

$cli = new Cli();
$frg = new Fragments('yb', 1, '');
$pims = new PIMS();
$tmr = $tmi = new ProgressTimer();

$test = $cli->getVar('t');
$test_table = false;
$path = $cli->getVar('p');
// $fpeak = $cli->getVar('m');
$npeak = $cli->getVar('n');
$ipeak = $cli->getVar('i');
$dpeak = $cli->getVar('d');
$fpeak = $cli->getVar('a');
$score_min = $cli->getVar('s');

if (is_null($npeak))
    $npeak = 40;

if (is_null($ipeak))
    $ipeak = '4,3,2,1,0';

if (is_null($dpeak))
    $dpeak = '17.5,27,35,70';

if (is_null($fpeak))
    $fpeak = 'int,iso';

if (is_null($score_min))
    $score_min = 9000;

$ipeaks = explode(",", $ipeak);
$dpeaks = explode(",", $dpeak);

$cp_n = cartesian_product([
    $ipeaks,
    $dpeaks,
    [
        'int'
    ]
]);
$cp_n = $cp_n->asArray();
if (! strstr($fpeak, "int"))
    $cp_n = [];

$cp_i = cartesian_product([
    $ipeaks,
    [
        0
    ],
    [
        'iso'
    ]
]);
$cp_i = $cp_i->asArray();
if (! strstr($fpeak, "iso"))
    $cp_i = [];

$cp_ui = array_merge($cp_n, $cp_i);
// print_r($cp_ui);

$file_regex = $cli->getVar('f');
$files = array_diff(scandir($path), array(
    '..',
    '.'
));

foreach ($files as $n => $file) {
    if (! preg_match("/" . $file_regex . "/", $file))
        unset($files[$n]);
}

$tmr->start();

foreach ($files as $file) {
    
    print_message('file', $file, ".", 40);
    
    $csv = new WriteDelim();
    $mgf = new ReadMGF($path . $file);
    $pims->setIndex('i_fpep_yeast');
    
    print_message('  time', $tmr->timeinstr(), ".", 40);
    
    $data = $mgf->getData();
    $tmr->start('  ' . count($data) . ' scans', count($data));
    $table = array();
    
    $pims = new PIMS();
    $pims->setIndex('i_fpep_yeast');
    
    $file = preg_replace("/\.mgf/", "", $file);
    $file_out = $path . "../pims/" . $file;
    $file_out .= "_pims_" . str_replace(",", "", $fpeak);
    $file_out .= "_d" . str_replace(",", "-", $dpeak);
    $file_out .= "_i" . str_replace(",", "-", $ipeak);
    $file_out .= "_n" . $npeak;
    $file_out .= ".csv";
    
    foreach ($data as $scan => $this_data) {
        
        $tmr->barPercent();
        $tmi->getTime();
        
        if ($test != NULL) {
            if ($scan != $test)
                continue;
            
            echo PHP_EOL;
        }
        
        $t_spec = new Timer();
        
        $spec = $this_data['spec'];
        $prec = round(preg_replace("/\s+.*/", "", $this_data['pepmass']), 4);
        
        $title = preg_replace("/\"|\,/", "", $this_data['title']);
        $rtsec = round(preg_replace("/\"|\,/", "", $this_data['rtinseconds']));
        
        $this_table = false;
        /*
         * int / iso peak picking
         */
        foreach ($cp_ui as $pks_ui) {
            
            $pks_c = $pks_ui[0];
            $pks_d = $pks_ui[1];
            $pks_f = $pks_ui[2];
            
            if ($pks_f == 'int' & $pks_c == 0 & $pks_d > 18)
                continue;
            
            // if ($pks_f == 'int' & $pks_c > 1 & $pks_d < 18)
            // continue;
            
            if ($pks_f == 'iso')
                $peaks = getPeaksByIsotope($spec, $npeak);
            
            if ($pks_f == 'int')
                $peaks = getPeaksByInt($spec, $prec, $npeak, $pks_d);
            
            if (is_false($peaks) or count($peaks) < 2)
                continue;
            
            $t_search = new Timer();
            
            $sphinxq = getSphinxQuery($prec, $peaks, $pks_c);
            
            if ($test != NULL) {
                print_message('config ', $pks_f . " n" . $npeak . " c" . $pks_c . " d" . $pks_d);
                // echo '$pks <- c(' . array_tostring($peaks, ",", "") . ")" . PHP_EOL;
                echo array_tostring($peaks, ",", "") . PHP_EOL;
                // echo PHP_EOL;
                // echo $sphinxq;
                // echo PHP_EOL;
            }
            
            $this_table = $pims->search($sphinxq);
            
            if (is_false($this_table))
                continue;
            
            $row_n = count($this_table['score']);
            $this_table['pks_c'] = array_fill(0, $row_n, $pks_c);
            $this_table['pks_f'] = array_fill(0, $row_n, $pks_f);
            $this_table['pks_n'] = array_fill(0, $row_n, count($peaks));
            $this_table['pks_d'] = array_fill(0, $row_n, $pks_d);
            
            $query_hash = array_values(array_map(array(
                $frg,
                'massToHash'
            ), ($peaks)));
            
            for ($i = 0; $i < $row_n; $i ++) {
                $this_peptide = $this_table['peptide'][$i];
                $this_hash = $frg->getFragmentWordArray($this_peptide);
                $this_overlap = count(array_intersect($this_hash, $query_hash));
                $this_table['pks_i'][$i] = count($this_hash);
                $this_table['pks_o'][$i] = $this_overlap;
            }
            
            $search_time = round($t_search->timeinmsec(), 3);
            $this_table['title'] = array_fill(0, $row_n, $title);
            $this_table['scan'] = array_fill(0, $row_n, $scan);
            $this_table['file'] = array_fill(0, $row_n, $file);
            $this_table['pre_mass'] = array_fill(0, $row_n, $prec);
            $this_table['pre_rt_s'] = array_fill(0, $row_n, $rtsec);
            $this_table['stime_sec'] = array_fill(0, $row_n, $search_time);
            
            if ($test != NULL) {
                $test_table = table_bind($test_table, $this_table);
            }
            
            if (! file_exists($file_out))
                $csv->writeCsvFile($file_out, $this_table, 'w');
            else
                $csv->writeCsvFile($file_out, $this_table, 'a');
            
            if (! is_false($this_table) and array_max($this_table['score']) > $score_min)
                break;
        }
        
        /*
         *
         */
        if ($test != NULL) {
            print_message('search time msec', time_toString($search_time, TRUE));
            echo PHP_EOL;
            echo $title . PHP_EOL;
            
            if (! is_false($test_table)) {
                unset($test_table['title']);
                $test_table = table_sort($test_table, 'pks_o', 'desc', 'number');
                print_table($test_table, 100);
            }
            exit();
        }
    }
    
    echo PHP_EOL;
    if ($test != NULL)
        continue;
}

?>