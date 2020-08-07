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
 * -rx : <enzyme_regex> -- default: "/.*?[K|R]/"
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

$fst->getFastaPath("/data/resources/uniprot/saccharomyces_db.fasta");
$ui_regex = $cli->getVar("r", "ALBU_HUMAN");

// echo PHP_EOL;


$prot = $fst->getHomology($ui_regex);
echo PHP_EOL;
print_r($prot);
echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

?>