<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "2048M");

use ProxySci\Omics\Proteomics\Digest;
use ProxySci\Omics\Proteomics\Fasta;
use ProxySci\Omics\Proteomics\Peptide;
use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Fragments;
use ProxyIO\File\Delim\WriteDelim;
use ProxyTime\ProgressTimer;

$tmr = new Timer();
$prg = new ProgressTimer();
$dig = new Digest();
$fst = new Fasta();
$amu = new Peptide();
$frg = new Fragments('yb', '1', '');
$csv = new WriteDelim();

echo PHP_EOL;
$prots = $fst->getAllEntries(TRUE);
echo PHP_EOL;
$prg->proTimerSize(count($prots));
$prg->proTimerStart("fraging " . count($prots) . " proteins");
foreach ($prots as $prot) {
    
    $prg->proTimerPrint();
    // $prot = $fst->getProtein('ALBU_HUMAN');
    $peps = $dig->getPeptides($prot['seq']);
    
    // $table['peptide'] = $peps;
    // $table['mass'] = array();
    foreach ($peps as $n => $seq) {
        $wtf = $frg->getFragmentMassArray($seq);
        
        $table = [];
        
        $table['frag'] = array_keys($wtf);
        $table['mass'] = array_values($wtf);
        
        $csv->writeCsvFile("./dat/frag_masses.csv", $table, 'a');
    }
}

echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

?>