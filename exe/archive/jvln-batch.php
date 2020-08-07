#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * jvln - jaevalen search tandem ms file
 *
 * SYNOPSIS
 * jvln -p <path> -f <regex>
 *
 * DESCRIPTION
 * searches the jaevalen engine for matching peptides based on the
 * precursor and fragments masses
 *
 * COMMAND LINE
 *
 * -t : flag for testing, input the scan N to test
 *
 * -p : path to data file
 *
 * -f : file regex as a means to filter
 *
 * -n : scan limit n (optional)
 *
 * -s : jscore minimum (default: 0.1)
 *
 * -l : hits per spectrum return limit (default: 3)
 *
 * -r : quorum match ratio (default: 0.50)
 *
 * -nmh : n top must have at least 1 (default: 5)
 *
 * -npn : n peaks min (default: 4)
 *
 * -npx : n peaks max (default: 21)
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
 * jvln -p /home/scbi/data -f 00 -m aiz -db i_fpep_uniprot
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxyTime\ProgressTimer;
use ProxyTime\Timer;
use Jvln\Engines\JvlnQuery;

$cli = new Cli();

$test = $cli->getVar('t', false);
$path = $cli->getVar('p', "./");
$file_regex = $cli->getVar('f', '');

$jvln_i = 'i_fpep_' . $cli->getVar('db', 'uniprot');
$jvln_sl = $cli->getVar('n');
$jvln_sm = $cli->getVar('s', 0.1);
$jvln_lh = $cli->getVar('l', 3);
$jvln_quorum = $cli->getVar('npq', 1);
$jvln_nmusthave = $cli->getVar('nph', 5);
$jvln_npeaksmax = $cli->getVar('npx', 21);
$jvln_npeaksmin = $cli->getVar('npn', 4);

$jvln_s = '127.0.0.1';
$jvln_p = 9307;
$jvln_db = NULL;
$jvln_un = 'root';
$jvln_pw = '';

/*
 * set up the JVLN query
 */
$jvn = new JvlnQuery($jvln_i);
$jvn->setQuorum($jvln_quorum);
$jvn->setNMustHavePeaks($jvln_nmusthave);
$jvn->setReturnLimit($jvln_lh);
$db_stats = $jvn->getDBStats($jvln_s);

$files = list_files($path, $file_regex . '.*\.json$', TRUE);

$pgt = new ProgressTimer();
$pgt->header("JVLN START");
$pgt->header("CONFIG");
$pgt->message(" file path", $path);
$pgt->message(" file regex", $file_regex);
$pgt->message(" file n", count($files));
$pgt->message(" jvln database", $jvln_i);
$pgt->message(" jvln score-min", $jvln_sm);
$pgt->message(" pks  limit", $jvln_lh);
$pgt->message(" pks  quorum", $jvln_quorum);
$pgt->message(" pks  n top must have", $jvln_nmusthave);
$pgt->message(" pks  n minimum", $jvln_npeaksmin);
$pgt->message(" pks  n maximum", $jvln_npeaksmax);

$pgt->header("BEGIN");
foreach ($files as $file) {
    
    $jvln_out = [];
    
    $jvln_out = [
        'parameters' => [
            'server' => $jvln_s,
            'database' => $db_stats,
            'psm_score_min' => $jvln_sm,
            'psm_limit' => $jvln_lh,
            'peak_quorum' => $jvln_quorum,
            'peak_n_top_must_have' => $jvln_nmusthave,
            'peak_n_minimum' => $jvln_npeaksmin
        ]
    ];
    
    $pgt = new ProgressTimer();
    
    $data = json_decode(file_get_contents($file), TRUE);
    
    $data_to_remove = array_diff(array_keys($data), [
        'file',
        'scans'
    ]);
    if (count($data_to_remove) > 0)
        foreach ($data_to_remove as $rm) {
            unset($data[$rm]);
        }
    
    $data_pre_mass = $data['scans']['pre_mass'];
    $data_pre_z = $data['scans']['pre_z'];
    $data_peaks = $data['scans']['peaks'];
    
    $scans = array_keys($data_pre_mass);
    
    if (! is_null($jvln_sl))
        $scans = array_slice($scans, 0, $jvln_sl);
    
    $pgt->start("n:" . count($scans) . " f:" . basename($file, '.json'), count($scans));
    
    /*
     * SQL SERVER CONN
     */
    
    $jvln_c = @new mysqli($jvln_s, $jvln_un, $jvln_pw, $jvln_db, $jvln_p);
    
    if ($jvln_c->connect_errno) {
        $pgt->message(" - connect error: " . $jvln_c->connect_errno, $jvln_c->connect_error);
        $pgt->message(" - thread id", $jvln_c_id);
        $pgt->header("EARLY TERMINATION", TRUE);
    }
    
    $jvln_c_id = $jvln_c->thread_id;
    
    /*
     * FOREACH SPECTRA -- START
     */
    $jvln_q = '';
    $scan_batch = [];
    foreach ($scans as $scan_i => $scan) {
        
        $tmr = new Timer();
        
        $pgt->barPercent();
        // $pgt->barCount();
        
        $this_pre_mass = $data_pre_mass[$scan];
        $this_pre_z = $data_pre_z[$scan];
        $this_peaks = $data_peaks[$scan];
        
        $jvln_data = [];
        /*
         * don't bother if we don't have enough peaks
         */
        if (count($this_peaks) >= $jvln_npeaksmin) {
            
            /*
             * slice down to the number of peaks requested
             */
            $this_peaks = array_slice($this_peaks, 0, $jvln_npeaksmax);
            
            $scan_batch[] = [
                'scan' => $scan,
                'peaks' => $this_peaks,
                'pre_mass' => $this_pre_mass,
                'pre_z' => $this_pre_z
            ];
            
            $this_pre_za = [
                $this_pre_z
            ];
            
            if ($this_pre_z == 0 | $this_pre_z > 5) {
                $this_pre_z = 2;
                $this_pre_za = range(1, 4);
            }
            
            $ia = range(0, 1);
            
            if ($scan_i == 5) {
                $jvln_q .= $jvn->getQuery($this_peaks, 100, $this_pre_za, $ia);
            } else {
                $jvln_q .= $jvn->getQuery($this_peaks, $this_pre_mass, $this_pre_za, $ia);
            }
            
            if (($scan_i + 1) % 10 == 0) {
                // if ($jvln_r = $jvln_c->query($jvln_q)) {
                // $jvln_data = $jvln_r->fetch_all(MYSQLI_ASSOC);
                // $jvln_r->free();
                // }
                
                if (! mysqli_multi_query($jvln_c, $jvln_q))
                    die("query failed");
                
                $time_sec = $tmr->timeinmsec();
                
                $scan_i_sub = - 1;
                do {
                    // fetch and print result set
                    if ($jvln_r = mysqli_store_result($jvln_c)) {
                        
                        $scan_i_sub ++;
                        
                        $this_scan = $scan_batch[$scan_i_sub]['scan'];
                        
                        // while ($row = mysqli_fetch_row($result))
                        // printf("id=%s\n", $row[0]);
                        $jvln_data = $jvln_r->fetch_all(MYSQLI_ASSOC);
                        $jvln_out['time_sec'][$this_scan] = round($time_sec / 10, 4);
                        
                        if (count($jvln_data) > 0) {
                            
                            $this_peaks = $scan_batch[$scan_i_sub]['peaks'];
                            $this_pre_mass = $scan_batch[$scan_i_sub]['pre_mass'];
                            $this_pre_z = $scan_batch[$scan_i_sub]['pre_z'];
                            
                            $jvln_data = jvln_hits_munge($jvln_data, $this_peaks, $this_pre_mass, $this_pre_z);
                            $jvln_out['psms'][$this_scan] = array_rowtocol($jvln_data);
                        }
                        
                        mysqli_free_result($jvln_r);
                    }
                    
                    // // print divider
                    // if (mysqli_more_results($jvln_c))
                    // printf("------\n");
                } while (mysqli_next_result($jvln_c));
                
                // print_r($jvln_out);
                
                // echo PHP_EOL;
                // echo $jvln_q;
                // echo PHP_EOL;
                // echo $scan_i;
                // echo PHP_EOL;
                
                $scan_batch = [];
                $jvln_q = '';
            }
        }
    }
    
    $jvln_err = $jvln_c->error;
    $jvln_c->close();
    
    $data['jaevalen'] = $jvln_out;
    
    file_put_contents($file, json_encode($data));
    
    /*
     * summary
     */
    // if (! is_false($stats_table))
    // print_table(table_summary($stats_table));
}

$pgt->header("COMPLETE");

?>