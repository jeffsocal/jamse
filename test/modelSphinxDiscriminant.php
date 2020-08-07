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
 * -p : <peptide_sequence> e.g. SAMPLER
 *
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
use Jvln\Engines\TandemHash;
use ProxySci\Omics\Proteomics\Fasta;
use ProxySci\Omics\Proteomics\Digest;
use Jvln\Engines\Spectrum;
use function BenTools\CartesianProduct\cartesian_product;
use ProxyTime\ProgressTimer;
use ProxyIO\Cli;

$prg = new ProgressTimer();
$cli = new Cli();
$amu = new Peptide();
$csv = new WriteDelim();

$dig = new Digest(7, 20);
$fst = new Fasta();
$jvln = new Spectrum('192.168.11.4');
$jvln->setIndex('i_fpep_uniprot');

$folder = $cli->getVariable('f', 'modeling');
$ui_series = $cli->getVariable('s', 'yb');
$ui_charge = $cli->getVariable('z', '12');
$ui_decays = $cli->getVariable('d', '');
$ui_massmax = $cli->getVariable('m', 1800);

$frg = new TandemHash($ui_series, $ui_charge, $ui_decays);
// echo PHP_EOL;

$prot = $fst->getProtein('CO3_HUMAN');
$peps = $dig->getPeptides($prot['seq']);
$peps = array_diff($peps, preg_grep("/[KR].+/", $peps));

// $cases['factors'] = '1';
// $cases['jvln01'] = '(doc_word_count*0.086 + sum(tf_idf)*13.5 - 0.95) * 5000';
// $cases['jvln02'] = '(doc_word_count*0.3883 + sum(sum_idf)*4.9 + bm25a*-4.354 + 1.0033e-13) * 5000';
$cases['jvln01'] = '(doc_word_count*0.085 + bm15*-0.029 + sum(sum_idf)*25.12 + 14.168) * (100000 / 2.67)';

$iteration = date("Ymdhis");

// $cases = expand_grid($cases);
// print_r($cases);
// exit;
function addMMA(float $val, $delta = 0.25)
{
    $delta *= 100;
    $val += (rand(0, $delta * 2) - $delta) / 100;
    return $val;
}

foreach ($peps as $pep) {
    
    $new_table = [];
    $prg->start($pep);
    
    $pre_z = 2;
    $pre_mz = addMMA(chargeMass($frg->getMolecularWeight($pep), $pre_z));
    
    $masses = $frg->getFragmentMassArray($pep);
    $masses = array_map('addMMA', $masses);
    
    foreach (range(12,2) as $frac) {
        
        for ($n = 0; $n < 100; $n ++) {
            
            shuffle($masses);
            
            $peaks = array_slice($masses, 0, $frac);
            $query_id = uniqueId();
            
            foreach ($cases as $case_it => $case) {
                // $jvln->setRankerFunction(array_tostring($case, '*', ''));
                $prg->barCount(5);
                
                if ($case == '')
                    continue;
                
                $jvln->setRankerFunction($case);
                
                $tmr = new Timer();
                $hits = $jvln->search($peaks, $pre_mz, $pre_z);
                
                print_r($hits);
                exit;
                
                if ($hits != FALSE) {
                    
                    $hits_n = count($hits);
                    
                    if (key_exists('factors', $hits[0])) {
                        foreach ($hits as $nh => $array) {
                            
                            $factors = json_decode($array['factors'], JSON_OBJECT_AS_ARRAY);
                            unset($array['factors']);
                            unset($hits[$nh]);
                            
                            $fields = $factors['fields'][0];
                            unset($factors['fields']);
                            unset($factors['words']);
                            
                            $hits[$nh] = array_merge($array, $factors, $fields);
                        }
                    }
                    
                    $hits = array_rowtocol($hits);
                    
                    $hits['query_frag_n'] = array_fill(0, $hits_n, count($peaks));
                    $hits['query_frag_r'] = array_fill(0, $hits_n, $frac / count($masses));
                    $hits['query_peptide'] = array_fill(0, $hits_n, $pep);
                    $hits['query_id'] = array_fill(0, $hits_n, $query_id);
                    $hits['query_function'] = array_fill(0, $hits_n, $case_it);
                    $hits['timer_msec'] = array_fill(0, $hits_n, $tmr->timeinmsec());
                    
                    $new_table = table_bind($new_table, $hits);
                }
            }
        }
    }
    echo PHP_EOL;
    $file_name = get_include_path() . "./dat/" . $folder . "/" . $pep . "_jvln_modeling_" . $iteration . ".csv";
    $csv->writeCsvFile($file_name, $new_table);
}

?>
