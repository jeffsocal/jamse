<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * digest
 *
 * SYNOPSIS
 * fragment -p <peptide_sequence>
 *
 * DESCRIPTION
 * display a table of fragment peaks for unit testing
 *
 * COMMAND LINE
 *
 * -p : <protein_name> ---- e.g. ALBU_HUMAN
 * -ul : <upper_length_limit> -- default: 36
 * -ll : <lower_length_limit> -- default: 6
 * -um : <upper_mass_limit> -- default: 6000
 * -lm : <lower_mass_limit> -- default: 600
 * -mc : <missed_clevage_max> -- default: 2
 * -rx : <enzyme_regex>
 *
 * EXAMPLE
 *
 * digest -p ALBU_HUMAN
 */
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxySci\Omics\Proteomics\Digest;
use ProxySci\Omics\Proteomics\Fasta;
use ProxyTime\Timer;
use ProxyIO\Cli;
use function BenTools\CartesianProduct\cartesian_product;

$cli = new Cli();
$tmr = new Timer();
$dig = new Digest();
$fst = new Fasta();

$ui_p = $cli->getVar("p", "ALBU_HUMAN");

// echo PHP_EOL;

$prots = explode("+", $ui_p);
$seqs = [];
foreach ($prots as $prot) {
    $prot_data = $fst->getProtein($prot);
    $seqs[$prot] = $prot_data['seq'];
}

$prots = array_keys($seqs);
for ($i = 1; $i < count($prots); $i ++) {
    echo $prots[$i - 1] . "\t";
    echo $prots[$i] . "\t";
    echo levenshtein($seqs[$prots[$i - 1]], $seqs[$prots[$i]]);
    echo PHP_EOL;
}

// echo PHP_EOL;
// print_table($table, table_length($table));
echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

?>