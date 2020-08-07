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

$cli = new Cli();
$tmr = new Timer();
$amu = new Peptide();
$csv = new WriteDelim();

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
$wtf_h = $frg->getFragmentWordArray($pep);

// print_r($wtf_m);
// print_r($wtf_h);

$wtf_keys = array_intersect(array_keys($wtf_m), array_keys($wtf_h));

$new_table['frag'] = array_values($wtf_keys);
$new_table['mass'] = array_values(array_intersect_key($wtf_m, array_flip($wtf_keys)));
$new_table['hash'] = array_values(array_intersect_key($wtf_h, array_flip($wtf_keys)));
$new_table['indx'] = array_map([
    $frg,
    'massToIndex'
], $new_table['mass']);

echo PHP_EOL;
print_table($new_table, 100);

$nmass = $amu->getMolecularWeight($pep);

$pre_table = [];
$pre_table['charge'] = [];
$pre_table['mz'] = [];
$pre_table['hash'] = [];
$pre_table['index'] = [];
for ($z = 0; $z < 4; $z ++) {
    
    $zm = chargeMass($nmass, $z);
    $masses = $frg->precursorMassSpread($zm);
    
    foreach ($masses as $mass) {
        $pre_table['charge'][] = $z;
        $pre_table['mz'][] = $mass;
        $pre_table['hash'][] = $frg->massToHash($mass, FALSE);
        $pre_table['index'][] = $frg->massToIndex($mass, FALSE);
        $pre_table['index_pre'][] = $frg->precursorMassToHash($mass);
    }
}

print_table($pre_table, 100);

echo PHP_EOL;
echo PHP_EOL;

function rnd2(float $float)
{
    return round($float, 2);
}

$z = 2;
shuffle($wtf_m);
array_splice($wtf_m, 0, floor(count($wtf_m) * $keep));
sort($wtf_m);
$wtf_m = array_map('rnd2', $wtf_m);

echo '?pmz=' . round(chargeMass($nmass, $z), 2) . '&pz=' . $z . '&pks=' . array_tostring($wtf_m, '+', '');
echo PHP_EOL;
echo PHP_EOL;
echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);
echo PHP_EOL;
print_message("peptide", $pep, '.', 30);
print_message(" length", strlen($pep), '.', 30);
echo PHP_EOL;
// $csv->writeCsvFile("/data/OpenAccess/MSdata/yeast_gs/peptides/" . $pep . ".fragments.csv", $new_table);

?>
