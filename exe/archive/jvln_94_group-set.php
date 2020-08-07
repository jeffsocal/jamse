#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * jvlnpro - jaevalen search tandem ms file
 *
 * SYNOPSIS
 * jvlnpro -p <path> -f <regex>
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
 * -l : return limit (default: 10)
 *
 * -db : database
 *
 * EXAMPLE
 *
 * jvlnpro -p /home/scbi/data -f 00 -m aiz -db i_fpep_uniprot
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxyTime\ProgressTimer;
use ProxyIO\File\Delim\WriteDelim;
use ProxyTime\Timer;

// $mgf = new ReadMGF("/data/SCBI/20171201_DBS_20Samples/mgf/20171201_SoCal_Samples18-9.mgf");

$cli = new Cli();
$csv = new WriteDelim();
$pgt = new ProgressTimer();

$test = $cli->getVar('t', false);
$path = preg_replace("/\/+$/", "", $cli->getVar('p')) . "/";
$file_regex = $cli->getVar('f', "\.results.prot.json$");

$files = list_files($path, $file_regex, TRUE);

$pgt->header("JVLN START");
$pgt->header("CONFIG");
$pgt->message(" file path", $path);
$pgt->message(" file regex", $file_regex);
$pgt->message(" file n", count($files));

$file_out = $path . "protein_group.json";

$jvln_table = [];
$pgt->header("BEGIN");
$prot_accounting_final = [];
$pept_data = [];

$pept_data = [];
$scan_data = [];
$prot_data = [];

foreach ($files as $file_path) {
    
    $file = preg_replace("/.*\//", "", $file_path);
    $cli->message(" ==> read ", $file);
    
    $json_cont = json_decode(file_get_contents($file_path), TRUE);
    
    $file = preg_replace("/.*\//", "", $file_path);
    
    $prots = $json_cont['proteins'];
    $scans = $json_cont['scans'];
    /*
     * SCAN DATA
     */
    $scan_i = 0;
    foreach ($scans as $scan) {
        if (count($scan['psms']) == 0)
            continue;
        $scan_i ++;
        
        $peptide_bare = preg_replace("/I/", "L", getRootSequence($scan['psms'][0]['peptide']));
        
        $scan_data[$scan_i] = $scan['psms'][0];
        $scan_data[$scan_i]['peptide_bare'] = $peptide_bare;
        
        if (! key_exists($peptide_bare, $pept_data))
            $pept_data[$peptide_bare] = 0;
        
        $pept_data[$peptide_bare] += $scan['psms'][0]['score'];
    }
    /*
     * PROTEIN DATA
     */
    $scan_i = 0;
    foreach ($prots as $prot_id => $prot) {
        if (! key_exists($prot_id, $prot_data)) {
            $prot_data[$prot_id] = $prot;
        } else {
            $prot_data[$prot_id]['peptides'] = array_unique(array_merge($prot_data[$prot_id]['peptides'], $prot['peptides']));
        }
    }
    unset($prots);
    unset($scans);
    unset($json_cont);
}

$total_peptides = count($pept_data);

$cli->message("total peptides", $total_peptides);

// $pept_data = [];
// $scan_data = [];
// $prot_data = [];

while (count($pept_data) > 0) {
    
    /*
     * PROTEIN DATA
     */
    
    $cli->flush(" ==> decomposing proteins ", "proteins: " . count($prot_data) . ' peptides: ' . count($pept_data) . ' of ' . $total_peptides);
    
    $keep_peptides = [];
    $prot_accounting = [];
    foreach ($prot_data as $protein_id => $protein_data) {
        
        $this_scores = array_intersect_key($pept_data, array_flip($protein_data['peptides']));
        $this_peptides = array_keys($this_scores);
        
        if (count($this_peptides) == 0) {
            unset($prot_data[$protein_id]);
            continue;
        }
        
        /*
         * calculate the protein group score
         */
        sort($this_peptides);
        $this_prot['id'] = $protein_id;
        $this_prot['hits'] = count($this_peptides);
        $this_prot['score'] = array_sum($this_scores);
        $this_prot['protein_group'] = md5(array_tostring($this_peptides, ' ', ''));
        
        /*
         * keep a lookup table of peptides (as keys) that belong to the group
         */
        $keep_peptides[$this_prot['protein_group']] = $this_scores;
        $prot_accounting[] = $this_prot;
    }
    
    if (count($prot_accounting) == 0) {
        break;
    }
    
    // reshape the data
    $prot_accounting = array_rowtocol($prot_accounting);
    
    // filter by score max
    $score_max = array_max($prot_accounting['score']);
    $prot_accounting = table_droprows_lt($prot_accounting, 'score', $score_max);
    $prot_group = $prot_accounting['protein_group'][0];
    
    $rows = array_diff($prot_accounting['protein_group'], [
        $prot_group
    ]);
    
    $prot_accounting = table_droprows($prot_accounting, $rows);
    
    $pept_data = array_diff_key($pept_data, $keep_peptides[$prot_group]);
    
    $prot_accounting_final[] = [
        'protein_group' => $prot_accounting['protein_group'][0],
        'protein_score' => $score_max,
        'proteins' => $prot_accounting['id'],
        'peptides' => $keep_peptides[$prot_accounting['protein_group'][0]]
    ];
}

$prot_accounting_final_r = array_rowtocol($prot_accounting_final);
$prot_accounting_final_grps = $prot_accounting_final_r['protein_group'];
unset($prot_accounting_final_r['protein_group']);

$prot_accounting_final_r = table_sort($prot_accounting_final_r, 'protein_score', 'desc', 'number');

$prot_accounting_final = array_rowtocol($prot_accounting_final_r);

$prot_accounting_final = array_combine($prot_accounting_final_grps, $prot_accounting_final);

$json_cont['protein_groups'] = $prot_accounting_final;

file_put_contents($file_out, json_encode($json_cont));

$pgt->header("COMPLETE");

?>