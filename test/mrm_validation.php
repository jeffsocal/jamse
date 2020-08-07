<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * fragment
 *
 * SYNOPSIS
 * fragment -p <peptide_sequence>
 *
 * DESCRIPTION
 * display a table of fragment peaks for unit testing
 *
 * COMMAND LINE
 *
 * -p : <peptide_sequence> -- e.g. SAMPLER
 * -k : <keep_percent> -- default 0.33
 * -s : <ion_series> -- <ybcza> -- (default: yb)
 * -z : <ion_series_charge> -- <1234> -- (default: 1)
 * -d : <ion_decay_series> -- <aw> ---- (default: NA)
 * -m : <ion_mass_max> ------ (default: 1800)
 *
 * EXAMPLE
 *
 * fragment -p SAMPLE
 */
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Peptide;
use ProxyIO\File\Delim\WriteDelim;
use ProxyIO\Cli;
use Jvln\Engines\TandemHash;
use Jvln\Engines\ReactionMonitoring;

$cli = new Cli();
$tmr = new Timer();
$amu = new Peptide();
$csv = new WriteDelim();
$jvln = new ReactionMonitoring('192.168.11.4');
$jvln->setIndex('i_fpep_uniprot');

$pep = $cli->getVariable('p');
$keep = min(1, max(0, 1 - $cli->getVariable('k', 0.33)));
$ui_series = $cli->getVariable('s', 'yb');
$ui_charge = $cli->getVariable('z', '1');
$ui_decays = $cli->getVariable('d', '');
$ui_massmax = $cli->getVariable('m', 1800);

$frg = new TandemHash($ui_series, $ui_charge, $ui_decays);
$frg->setMassRange(76, $ui_massmax);
// $this->m_int = 584.9754;
// $this->m_slp = 0.9995063;

$wtf_m = $frg->getFragmentMassArray($pep);

// print_r($wtf_m);
// print_r($wtf_h);

echo PHP_EOL;
print_message("peptide", $pep, '.', 30);
$nmass = $amu->getMolecularWeight($pep);
preg_match_all("/\d/", $ui_charge, $ui_charges);

$pre_z = 1;
$pre_mz = chargeMass($nmass, $pre_z);

$hit = [
    'insilico_fragments' => $wtf_m
];
foreach ($wtf_m as $fragment => $mass1) {
    if ($mass1 > 1800)
        continue;
    
    foreach ($wtf_m as $fragment2 => $mass2) {
        if ($fragment == $fragment2)
            continue;
        
        if ($mass2 > 1800)
            continue;
        
        if ($mass2 > $mass1)
            continue;
        
        $query = array(
            $fragment => $mass1,
            $fragment2 => $mass2
        );
        $results = $jvln->search($query, $pre_mz, $pre_z);
        
        echo $fragment . " + " . $fragment2 . " => " . $results . PHP_EOL;
    }
}
// print_r($hit);
echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

// $csv->writeCsvFile("/data/OpenAccess/MSdata/yeast_gs/peptides/" . $pep . ".fragments.csv", $new_table);

?>
