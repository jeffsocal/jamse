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
use ProxySci\Omics\Proteomics\Peptide;
use ProxyTime\Timer;
use ProxyIO\Cli;

$cli = new Cli();
$tmr = new Timer();
$dig = new Digest();
$fst = new Fasta();
$amu = new Peptide();

$ui_p = $cli->getVar("p", "ALBU_HUMAN");
$ui_ul = $cli->getVar("ul", 53);
$ui_ll = $cli->getVar("ll", 6);
$ui_um = $cli->getVar("um", 6000);
$ui_lm = $cli->getVar("lm", 400);
$ui_mc = $cli->getVar("mc", 3);
$ui_rx = $cli->getVar("rx", "/.*?[KR\#]/");

// echo PHP_EOL;

$dig->setLengthLimits($ui_ll, $ui_ul);
$dig->setMassLimits($ui_lm, $ui_um);
$dig->setMissedClevageMax($ui_mc);
$dig->setEnzymeRegex($ui_rx);

$prot = $fst->getProtein($ui_p);
$peps = $dig->getPeptides($prot['seq']);


$table['n'] = range(1, count($peps), 1);
$table['peptide'] = $peps;

echo PHP_EOL;
print_table($table, 25);
// print_table(table_summary($table));
print_r($prot);
echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

?>