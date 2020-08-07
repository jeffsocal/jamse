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
 * metabolites_hmdbxml -p <xml_path>
 *
 * DESCRIPTION
 * parse the hmdb-xml database file - to csv
 *
 * xml-dataset is about 4.5G which is too big for
 * SimpleXML, therefor we will just read line-by-line
 * thankfully the datafile is CRLF for each element
 *
 * COMMAND LINE
 *
 * -p : <xml_path>
 *
 * EXAMPLE
 *
 */
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Peptide;
use ProxyIO\File\Delim\WriteDelim;
use ProxyIO\Cli;
use ProxyTime\ProgressTimer;

$cli = new Cli();
$tmr = new Timer();
$prg = new ProgressTimer();
$amu = new Peptide();
$csv = new WriteDelim();

$keys = $these_keys = [
    'id',
    'accession',
    'mode',
    'peaks'
];

$path = $cli->getVariable('p');
$fout = $cli->getVariable('o', "../hmdb_metabolite_fragmasses");

$files = list_files($path, "xml", TRUE);
$prg->start("parse msms files", count($files));

foreach ($files as $n => $file) {
    $prg->barPercent();
    
    $fxml = new SimpleXMLElement(file_get_contents($file));
    
    $peaks = [];
    foreach ($fxml->{"ms-ms-peaks"}->{"ms-ms-peak"} as $peak) {
        $peaks[] = $peak->{"mass-charge"}->__toString();
        // $peaks['int'][] = $peak->{"intensity"}->__toString();
    }
    
    $table['peaks'] = $peaks;
    $table['id'] = array_fill(0, count($peaks), $fxml->{"id"}->__toString());
    $table['accession'] = array_fill(0, count($peaks), $fxml->{"database-id"}->__toString());
    $table['mode'] = array_fill(0, count($peaks), strtolower($fxml->{"ionization-mode"}->__toString()));
    
    if (file_exists($path . $fout . ".csv"))
        $csv->writeCsvFile($path . $fout . ".csv", $table, 'a');
    else
        $csv->writeCsvFile($path . $fout . ".csv", $table);
    
}

// $csv = new WriteDelim();
// $tmr = new Timer();

// $table = array_rowtocol($data);

// print_message("TOTAL", $tmr->timeinstr(TRUE), '.', 30);

?>