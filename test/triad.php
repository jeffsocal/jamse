#!/usr/bin/env php
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
 * triad -p <peptide_sequence>
 *
 * DESCRIPTION
 * create the peptide triad
 *
 * COMMAND LINE
 *
 * -p : <peptide_sequence>
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

$seq = $cli->getVariable('p');

$seqs = str_split('__' . $seq . '__');
$content = [];
for ($i = 2; $i < count($seqs); $i ++) {
    $content[] = $seqs[$i - 2] . $seqs[$i - 1] . $seqs[$i];
}

echo PHP_EOL;

echo array_tostring($content, ' ', '') . "     " . count($content);
echo PHP_EOL;

echo strtolower(array_tostring($content, ' ', '')) . "     " . count($content);
echo PHP_EOL;

echo PHP_EOL;

?>