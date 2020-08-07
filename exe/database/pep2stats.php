#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * pepstats
 *
 * SYNOPSIS
 * pepstats
 *
 * DESCRIPTION
 * generates the sphinx database given config json file
 *
 * -f: <file_name>
 * 
 * EXAMPLE
 *
 * pepstats -f ./data/peptides
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxySci\Omics\Proteomics\Peptide;

$cli = new Cli();

$file = $cli->getVar('f');

$amu = new Peptide();

$handle = fopen($file, "r");
$n = 0;
$id = 0;
if ($handle) {
    echo "peptide,length,tryp_mc,molwt,homology" . PHP_EOL;
    while (($line = fgets($handle)) !== false) {
        $n ++;
        
        $line = explode(',', strtoupper($line));
        $seq = $line[0];
        $hom = $line[2];
        
        $mw = $amu->getMolecularWeight($seq);
        
        $mc = max(0, count(preg_grep("/K|R/", str_split($seq))) - 1);
        
        echo $seq . "," . strlen($seq) . "," . $mc . "," . $mw . "," . $hom . PHP_EOL;
    }
    
    fclose($handle);
} else {
    echo "FAILED TO READ LINE";
}

?>