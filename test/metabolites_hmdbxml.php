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
 * metabolites_hmdbxml -f <xml_file>
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
 * -f : <xml_file>
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

$file = $cli->getVariable('f');

$csv = new WriteDelim();
$tmr = new Timer();

if (! preg_match("/.XML$/i", $file))
    $cli->message("error", "bad file", true);

print_message("file", $file, '.', 30);

$keys = $these_keys = [
    'accession',
    'name',
    'chemical_formula',
    'monisotopic_molecular_weight',
    'iupac_name',
    'traditional_iupac',
    'cas_registry_number',
    'smiles',
    'inchi',
    'inchikey'
];

$n = 0;
$data = [];

$prg->start("parsing xml");

$handle = fopen($file, "r");
if ($handle) {
    while (($line = fgets($handle)) !== false) {
        // process the line read.
        
        if (preg_match('/\<metabolite\>/', $line)) {
            $n ++;
            
            $prg->barCount();
            
            // if ($n == 100)
            // break;
            
            $data[$n] = array_combine($keys, array_fill(0, count($keys), ''));
            $these_keys = $keys;
            continue;
        }
        
        foreach ($these_keys as $i => $key) {
            
            if (preg_match('/\<' . $key . '\>/', $line)) {
                $xml = new SimpleXMLElement($line);
                $data[$n][$key] = $xml->__toString();
                
                if (preg_grep('/' . $key . '/', [
                    'accession',
                    'name'
                ]))
                    unset($these_keys[$i]);
                
                if ($key == 'traditional_iupac')
                    $data[$n][$key] = clean_nonascii($data[$n][$key]);
                
                break;
            }
        }
    }
    
    fclose($handle);
} else {
    // error opening the file.
}

$table = array_rowtocol($data);
$csv->writeCsvFile($file . ".csv", $table);

print_message("TOTAL", $tmr->timeinstr(TRUE), '.', 30);

?>