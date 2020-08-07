#!/usr/bin/env php
<?php
use ProxyIO\Cli;
use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Fasta;

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * prot2sap - digest proteins for indexing homology
 *
 * SYNOPSIS
 * prot2sap
 *
 * DESCRIPTION
 * generates the sphinx xml idexing file for peptides based
 * on the input ini file
 *
 * COMMAND LINE
 *
 * <none>
 *
 * EXAMPLE
 *
 * prot2sap
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$tmr = new Timer();

$rm = "\!|\@|\#|\$|\%|\^|\&|\*|\(|\)|\[|\]|\{|\}|\<|\>";

/*
 * get INI definitions
 */
$ui_path_json = list_files('./', ".json", TRUE);
if (count($ui_path_json) == 0)
    $cli->message('no config file found', "", TRUE);
$ui_path_json = $ui_path_json[0];

$ini = parse_json_file($ui_path_json);

/*
 * SET THE DATABASE
 */
$path = $ini['fasta']['path'];
$files = $ini['fasta']['files'];

function ipc($p, $n = 0)
{
    $seq = $p['seq'];
    $out = '';
    
    if (strlen($seq) <= 5)
        return false;
    
    $out .= '<sphinx:document id="' . $n . '">' . PHP_EOL;
    $out .= '<protein_name>' . $p['name'] . '</protein_name>' . PHP_EOL;
    $out .= '<protein_id>' . $p['id'] . '</protein_id>' . PHP_EOL;
    $out .= '<protein_desc>' . $p['desc'] . '</protein_desc>' . PHP_EOL;
    $out .= '<organism>' . $p['org'] . "</organism>" . PHP_EOL;
    $out .= '<sequence>' . $seq . "</sequence>" . PHP_EOL;
    $out .= '<content>';
    
    $seq = preg_replace("/K/", "K_", $seq);
    $seq = preg_replace("/R/", "R_", $seq);
    
    $seqs = str_split('_' . $seq . '_');
    for ($i = 2; $i < count($seqs); $i ++) {
        $out .= strtolower($seqs[$i - 2] . $seqs[$i - 1] . $seqs[$i] . " ");
    }
    
    $out .= "</content>" . PHP_EOL;
    $out .= '</sphinx:document>' . PHP_EOL . PHP_EOL;
    return $out;
}

$head = '<?xml version="1.0" encoding="utf-8"?>
        <sphinx:docset>
        <sphinx:schema>
        <sphinx:field name="content" />
        <sphinx:attr name="protein_name" type="string" />
        <sphinx:attr name="protein_id" type="string" />
        <sphinx:attr name="protein_desc" type="string" />
        <sphinx:attr name="organism" type="string" />
        <sphinx:attr name="sequence" type="string" />
        </sphinx:schema>
        ' . PHP_EOL;

$n = 0;
$entry = FALSE;

echo $head;

foreach ($files as $file) {
    
    $handle = fopen($path . $file, "r");
    $fst = new Fasta($file);
    
    if ($handle) {
        while (($line = fgets($handle)) !== false) {
            /*
             * process the line read.
             */
            if (preg_match('/^>/', $line)) {
                
                if (! is_false($entry)) {
                    $n ++;
                    
                    $entry = preg_replace("/" . $rm . "/", "", $entry);
                    $entry = preg_replace("/\//", "-", $entry);
                    echo ipc($fst->parseFastaEntry($entry), $n);
                    
                    // if ($n % 100000 == 0)
                    // print_message($n, $tmr->timeinstr());
                }
                
                $entry = trim($line, '>') . " ";
            } else {
                $entry .= $line;
            }
        }
        /*
         * get the last one
         */
        $n ++;
        $entry = preg_replace("/" . $rm . "/", "", $entry);
        $entry = preg_replace("/\//", "-", $entry);
        echo ipc($fst->parseFastaEntry($entry), $n);
        
        fclose($handle);
    } else {
        // error opening the file.
    }
}

echo '</sphinx:docset>' . PHP_EOL;

// UniProtKB 125,355,233
// - UniProtKB/Swiss-Prot 558,125 ->
// - UniProtKB/TrEMBL 124,797,108
// UniRef100 154,752,138
// UniRef90 78,915,455
// UniRef50 32,365,088
// UniParc 225,932,239

?>
