#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * jvlnnt - jaevalen parallel search tandem ms file
 *
 * SYNOPSIS
 * jvlnnt -p <path> -f <regex>
 *
 * DESCRIPTION
 * searches the jaevalen engine for matching peptides based on the
 * precursor and fragments masses
 *
 * COMMAND LINE
 *
 * -t : flag for testing, input the scan N to test
 *
 * -nt : number of threads (default: 1)
 *
 * -p : path to data file
 *
 * -f : file regex as a means to filter
 *
 * -s : jscore minimum (default: 9000)
 *
 * -l : return limit (default: 10)
 *
 * -r : quorum match ratio (default: 0.50)
 *
 * -n : n top must have (default: 1)
 *
 * -db : data base
 *
 * -m : precursor search mode
 * -- a -- asis ........ considers neutral mass i.e. (mz-1)*z
 * -- i -- isotope ..... considers the isotopes a0 a1 a2 were acquired
 * -- z -- m/z ......... considers charged mass (mz) at charge states 1 2 3
 *
 * EXAMPLE
 *
 * jvlnnt -p /home/scbi/data -f 00 -m aiz -db i_fpep_uniprot
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxyIO\File\Write;
use Jvln\System\Passthru;

$cli = new Cli();
$wrt = new Write();
$exe = new Passthru();

$ui['file_regex'] = $cli->getVar('f', "\.csv$");
// $u['t'] = $ui['test'] = $cli->getVar('t', NULL);
$u['nt'] = $ui['number of threads'] = $cli->getVar('nt', 1);
$u['p'] = $ui['path'] = preg_replace("/\/+$/", "", $cli->getVar('p')) . "/";
$u['db'] = $ui['search database'] = $cli->getVar('db');
$u['t'] = $ui['search mode'] = $cli->getVar('m', "a");
$u['s'] = $ui['search score min'] = $cli->getVar('s', 100);
$u['m'] = $ui['search return n max'] = $cli->getVar('l', 10);
$u['r'] = $ui['peak quorum ratio'] = $cli->getVar('r', 0.5);
$u['nmh'] = $ui['peak n must have'] = $cli->getVar('nmh', 1);
$u['npm'] = $ui['peak min n'] = $cli->getVar('npm', 4);

echo PHP_EOL;
$exe->cli->header("START JAEVALEN SEARCH");
$exe->cli->header("CONFIG VALUES");
foreach ($ui as $param => $val) {
    $exe->cli->message($param, $val);
}

$files = list_files($ui['path'], $ui['file_regex']);

foreach ($files as $file) {
    /*
     * shard file
     */
    $exe->cli->header($file);
    $file_tmp = $ui['path'] . $file . '_body';
    $clis = [
        'rm ' . $ui['path'] . 'x*_' . $file,
        'tail -n +2 ' . $ui['path'] . $file . ' > ' . $file_tmp,
        'split -n r/' . $u['nt'] . ' ' . $file_tmp . ' --additional-suffix=_' . $file,
        'head -1 ' . $ui['path'] . $file . ' > ' . $ui['path'] . $file . '_head'
    ];
    foreach ($clis as $cli) {
        if (is_false($exe->dohup($cli, true)))
            $exe->cli->header("ERROR - EARLY TERMINATION", true);
    }
    /*
     * locate sharded files
     */
    $shard_files = list_files($ui['path'], "^x.*" . $file);
    
    /*
     * add header back into sharded files
     */
    $clis = [];
    foreach ($shard_files as $shard_file) {
        $clis = array_merge($clis, [
            'cat ' . $ui['path'] . $file . '_head ' . $ui['path'] . $shard_file . ' > ' . $ui['path'] . $shard_file . '_tmp',
            'mv ' . $ui['path'] . $shard_file . '_tmp ' . $ui['path'] . $shard_file
        ]);
    }
    
    foreach ($clis as $cli) {
        if (is_false($exe->dohup($cli, false)))
            $exe->cli->header("ERROR - EARLY TERMINATION", true);
    }
    
    /*
     * parallel search
     */
    foreach ($shard_files as $shard_file) {
        $u['f'] = $shard_file;
        $cli = 'jvlncsv';
        foreach ($u as $p => $v) {
            $cli .= ' -' . $p . ' ' . $v;
        }
        
        if (is_false($exe->nohup($cli, TRUE)))
            $exe->cli->header("ERROR - EARLY TERMINATION", TRUE);
        $pids[] = $exe->getPID();
    }
    $exe->monitorPIDs($pids, 'parallel jvlncsv');
    
    /*
     * clean up shard
     */
    $clis = [
        'rm ' . $ui['path'] . 'x*_' . $file,
        'rm ' . $ui['path'] . $file . "_*"
    ];
    foreach ($clis as $cli) {
        if (is_false($exe->dohup($cli, true))) {
            $exe->cli->header("ERROR - EARLY TERMINATION", true);
        }
    }
}

$exe->cli->header("COMPLETED");

echo PHP_EOL;

?>