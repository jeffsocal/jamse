#!/usr/bin/env php
<?php
use ProxyIO\Cli;
use ProxySci\Omics\Proteomics\Peptide;
use ProxySci\Omics\Proteomics\PeptideModifications;
use Jvln\System\JsonConfig;

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
 * -i : <index_name> found as ini/index_name.ini
 *
 * EXAMPLE
 *
 * pep2ptm -i yeast -f peptides.csv
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$jsc = new JsonConfig();

$file = $cli->getVar('f');
$test = $cli->getVar('t');

/*
 * get INI definitions
 */
if ($jsc->getConfigFile() == false)
    $cli->message('no config file found', "", TRUE);

$unimod = $jsc->getVariable('unimod');
$ui_ul = $jsc->getVariable('peptide-ptm', 'length_upper_limit', 48);
$ui_ll = $jsc->getVariable('peptide-ptm', 'length_lower_limit', 6);
$ui_um = $jsc->getVariable('peptide-ptm', 'mass_upper_limit', 4500);
$ui_lm = $jsc->getVariable('peptide-ptm', 'mass_lower_limit', 600);
$ui_mc = $jsc->getVariable('peptide-ptm', 'missed_clevage_max', 2);

$amu = new Peptide();
$mod = new PeptideModifications($unimod);

$handle = fopen($file, "r");
$n = 0;
$id = 0;
if ($handle) {
    while (($line = fgets($handle)) !== false) {
        $n ++;
        // if ($n < 8)
        // continue;
        
        $seq = trim(strtoupper(preg_replace("/\,.+/", "", $line)));
        $seq_len = strlen($seq);
        if ($seq_len < $ui_ll || $seq_len > $ui_ul)
            continue;
        
        /*
         * remove over weight peptides
         */
        $mw = $amu->getMolecularWeight($seq);
        if ($mw < $ui_lm || $mw > $ui_um)
            continue;
        
        /*
         * skip over large miscleavages
         */
        $mc = max(0, count(preg_grep("/K|R/", str_split($seq))) - 1);
        if ($mc > $ui_mc)
            continue;
        
        $peps = $mod->getPTMPeptides($seq);
        
        foreach ($peps as $pep) {
            
            $mw = $amu->getMolecularWeight($pep);
            if ($mw > $ui_um)
                continue;
            
            if ($test == 'debug')
                continue;
            echo $pep . PHP_EOL;
        }
        
        if ($test == 'debug')
            print_row_neat([
                $seq,
                $mc,
                $mw,
                $id,
                memory_get_peak_usage() / (1e6) . " MB",
                memory_get_usage() / (1e6) . " MB"
            ], [
                37,
                4,
                17,
                8,
                20,
                20
            ]);
    }
    
    fclose($handle);
} else {
    echo "FAILED TO READ LINE";
}

?>