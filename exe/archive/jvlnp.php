#!/usr/bin/env php
<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * jvlnp - jaevalen search by peaks
 *
 * SYNOPSIS
 * jvlnp -p <precursor_mass> -f '<fragment_masses>'
 *
 * DESCRIPTION
 * searches the jaevalen engine for matching peptides based on the
 * precursor and fragments masses
 *
 *
 * COMMAND LINE
 *
 * -p : <float> precursor mass value
 *
 * -f : '<float, array>' fragment mass values must be enclosed in single or double quotes
 *
 *
 * EXAMPLE
 *
 * akert jvlnp -p 936.493 -f '365.23 436.38 549.37'
 */

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxyIO\File\Delim\WriteDelim;
use ProxySci\Omics\Proteomics\Fragments;
use ProxyTime\ProgressTimer;
use ProxyTime\Timer;
use Jvln\Engines\PIMS;

$cli = new Cli();
$pims = new PIMS();
$frg = new Fragments('yb', 1, '');
$tmr = $tmi = new ProgressTimer();
$csv = new WriteDelim();

$peaks = $cli->getVar('f');
$prec = $cli->getVar('p');

$tmr->start();

$pims->setIndex('i_fpep_yeast');

if (is_false($peaks))
    exit();

$peaks = explode(" ", $peaks);
$pks_c = 1;
$pks_d = count($peaks);

$t_search = new Timer();

$sphinxq = getSphinxQuery($prec, $peaks, 0);

$this_table = $pims->search($sphinxq);

if (is_false($this_table))
    exit();

$row_n = count($this_table['score']);
$this_table['pks_c'] = array_fill(0, $row_n, $pks_c);
$this_table['pks_d'] = array_fill(0, $row_n, $pks_d);

$query_hash = array_values(array_map(array(
    $frg,
    'massToHash'
), ($peaks)));

for ($i = 0; $i < $row_n; $i ++) {
    $this_peptide = $this_table['peptide'][$i];
    $this_hash = $frg->getFragmentWordArray($this_peptide);
    $this_overlap = count(array_intersect($this_hash, $query_hash));
    $this_table['pks_i'][$i] = count($this_hash);
    $this_table['pks_o'][$i] = $this_overlap;
}

$search_time = round($t_search->timeinmsec(), 3);
$this_table['pre_mass'] = array_fill(0, $row_n, $prec);
$this_table['stime_sec'] = array_fill(0, $row_n, $search_time);

print_message('search time msec', time_toString($search_time, TRUE));
echo PHP_EOL;

$test_table = table_sort($this_table, 'pks_o', 'desc', 'number');
print_table($this_table, 100);

?>
