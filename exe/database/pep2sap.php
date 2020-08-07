#!/usr/bin/env php
<?php
use ProxyIO\Cli;
use ProxySci\Omics\Proteomics\Peptide;
use ProxySci\Omics\Proteomics\Fragments;

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * pep2ptm - apply ptms to digested peptides
 *
 * SYNOPSIS
 * pep2sap -f <ptm_file>
 *
 * DESCRIPTION
 *
 * COMMAND LINE
 *
 * -f : <full_path_to_file>
 *
 * EXAMPLE
 *
 * ptm2sap -f /var/sphinx/yeast/data/peptides
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$amu = new Peptide();
$frg = new Fragments();
$file = $cli->getVar('f');
$min_seq_leg = $cli->getVar('n', 5);

echo '<?xml version="1.0" encoding="utf-8"?>
        <sphinx:docset>
        <sphinx:schema>
        <sphinx:field name="content" />
        <sphinx:attr name="peptide" type="string" />
        <sphinx:attr name="length" type="int" />
        </sphinx:schema>
        ' . PHP_EOL;

$handle = fopen($file, "r");
$n = 0;

if ($handle) {
    while (($line = fgets($handle)) !== false) {
        $n ++;
        
        $seq = trim(strtoupper(preg_replace("/\,.+/", "", $line)));
        
        if (strlen($seq) < $min_seq_leg)
            continue;
        
        $mc = max(0, count(preg_grep("/K|R/", str_split($seq))) - 1);
        if ($mc > 2)
            continue;
        
        $seqs = str_split('__' . $seq . '__');
        $content = "";
        for ($i = 2; $i < count($seqs); $i ++) {
            $content .= $seqs[$i - 2] . $seqs[$i - 1] . $seqs[$i] . " ";
        }
        echo '<sphinx:document id="' . hexdec(uniqid()) . '">' . PHP_EOL;
        echo '<peptide>' . $seq . '</peptide>' . PHP_EOL;
        echo '<length>' . strlen($seq) . '</length>' . PHP_EOL;
        echo '<content>' . $content . "</content>" . PHP_EOL;
        echo '</sphinx:document>' . PHP_EOL . PHP_EOL;
    }
    
    fclose($handle);
} else {
    echo "FAILED TO READ LINE";
}
echo '</sphinx:docset>' . PHP_EOL;

?>