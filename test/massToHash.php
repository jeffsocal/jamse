<?php
use ProxySci\Omics\Proteomics\Fragments;
use ProxyIO\Cli;

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
 * EXAMPLE
 *
 * fragment -p SAMPLE
 */

set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$frg = new Fragments();

ini_set("memory_limit", "4096M");

$zs = $cli->getVariable('n');
$zs = explode(",", $zs);

if ($zs[0] == '')
    $zs = [
        922.4038696,
        461.9765015,
        583.6765747,
        850.295166,
        335.2226563
    ];

foreach ($zs as $z) {
    echo $frg->massToHash($z).PHP_EOL;
//     print_message($z, $frg->massToHash($z, TRUE));
}

?>