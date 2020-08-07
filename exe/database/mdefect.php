#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * mdefect - expose fragmentation values for stats calc
 *
 * SYNOPSIS
 * mdefect -f <ptm_file>
 *
 * DESCRIPTION
 *
 * COMMAND LINE
 *
 * -f : <full_path_to_file>
 * -s : <ion series> abycz                  default = yb
 * -z : <charge state> 1,2,3,4              default = 1
 * -d : <meta decay> aw (amonia / water)    default = 
 * -l : <min peptide length>                default = 6
 * -m : <molwt max>                         default = 6000
 * -c : <mis-clevage max>                   default = 2
 *
 * EXAMPLE
 *
 * mdefect -f /var/sphinx/yeast/data/xaa_peptides_ptm_acccd
 *
 */
use ProxyIO\Cli;
use ProxySci\Omics\Proteomics\Fragments;
use ProxySci\Omics\Proteomics\Peptide;

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$file = $cli->getVar('f');
$series = $cli->getVar('s', 'yb');
$charge = $cli->getVar('z', 1);
$decay = $cli->getVar('d', '');
$len_min = $cli->getVar('l', 6);
$mwt_max = $cli->getVar('m', 6000);
$mcl_max = $cli->getVar('c', 2);

$frg = new Fragments($series, $charge, $decay);
$amu = new Peptide();

$handle = fopen($file, "r");
if ($handle) {
    while (($line = fgets($handle)) !== false) {
        $pep = trim(strtoupper(preg_replace("/\,.+/", "", $line)));

        /*
         * remove short sequences
         */
        if (strlen(getRootSequence($pep)) <= $len_min)
            continue;

        /*
         * remove over weight peptides
         */
        $mw = $amu->getMolecularWeight($pep);
        if ($mw > $mwt_max)
            continue;

        /*
         * skip over large miscleavages
         */
        if (max(0, count(preg_grep("/K|R/", str_split($pep))) - 1) > $mcl_max)
            continue;

        $m = $frg->getFragmentMassArray($pep);

        echo array_tostring($m, PHP_EOL, '');
        echo PHP_EOL;
    }

    fclose($handle);
} else {
    echo "FAILED TO READ LINE";
}

?>