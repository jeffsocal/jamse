#!/usr/bin/env php
<?php

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
 * pep2ptm -i <index_name> -f <ptm_file>
 *
 * DESCRIPTION
 *
 * COMMAND LINE
 *
 * -f : <full_path_to_file>
 *
 * -i : <path_to_ini>
 *
 * EXAMPLE
 *
 * ptm2frag -f /var/sphinx/yeast/data/xaa_peptides_ptm_acccd
 *
 */
use ProxyIO\Cli;
use ProxySci\Omics\Proteomics\Peptide;
use Jvln\Engines\TandemHash;

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$amu = new Peptide();
$ui_file = $cli->getVar('f');

$ui_path_json = list_files('./', ".json", TRUE);
if (count($ui_path_json) == 0)
    $cli->message('no config file found', "", TRUE);
$ui_path_json = $ui_path_json[0];

$ion_series = 'yb';
$ion_charge = '12';
$ion_decay = '';
$ion_mass_max = 1800;

if (! is_null($ui_path_json)) {
    $ini = parse_json_file($ui_path_json);
    
    if (key_exists('fragmentation', $ini)) {
        
        if (key_exists('ion_series', $ini['fragmentation']))
            $ion_series = $ini['fragmentation']['ion_series'];
        
        if (key_exists('ion_charge', $ini['fragmentation']))
            $ion_charge = $ini['fragmentation']['ion_charge'];
        
        if (key_exists('ion_decay', $ini['fragmentation']))
            $ion_decay = $ini['fragmentation']['ion_decay'];
        
        if (key_exists('ion_mass_max', $ini['fragmentation']))
            $ion_mass_max = $ini['fragmentation']['ion_mass_max'];
    }
}

$frg = new TandemHash($ion_series, $ion_charge, $ion_decay);
$frg->setMassCutoff($ion_mass_max);

echo '<?xml version="1.0" encoding="utf-8"?>
        <sphinx:docset>
        <sphinx:schema>
        <sphinx:field name="content" />
        <sphinx:attr name="peptide" type="string" />
        </sphinx:schema>
        ' . PHP_EOL;

$handle = fopen($ui_file, "r");
if ($handle) {
    while (($line = fgets($handle)) !== false) {
        
        $pep = trim(strtoupper(preg_replace("/\,.+/", "", $line)));
        $mw = $amu->getMolecularWeight($pep);
        
        $this_id = hexdec(uniqid());
        
        $w = $frg->getFragmentWordArray($pep);
        
        /*
         * 0.3 da precursor mass tolerance
         */
        $mwh = $frg->precursorMassToHash($mw);
        
        /*
         * 1 da precursor mass tolerance
         */
        // $mwh = 'pmv' . $frg->massToHash($mw);
        
        echo '<sphinx:document id="' . $this_id . '">' . PHP_EOL;
        echo '<peptide>' . $pep . '</peptide>' . PHP_EOL;
        echo '<content>' . $mwh . ' ' . array_tostring($w, ' ', '') . "</content>" . PHP_EOL;
        echo '</sphinx:document>' . PHP_EOL . PHP_EOL;
    }
    
    fclose($handle);
} else {
    echo "FAILED TO READ LINE";
}
echo '</sphinx:docset>' . PHP_EOL;

?>