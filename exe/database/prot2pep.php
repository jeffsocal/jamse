#!/usr/bin/env php
<?php
use ProxyIO\Cli;
use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Fasta;
use ProxySci\Omics\Proteomics\Digest;
use Jvln\System\JsonConfig;

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
$dig = new Digest();
$jsc = new JsonConfig();

$rm = "\!|\@|\#|\$|\%|\^|\&|\*|\(|\)|\[|\]|\{|\}|\<|\>";

/*
 * get INI definitions
 */
if ($jsc->getConfigFile() == false)
    $cli->message('no config file found', "", TRUE);


/*
 * SET THE DATABASE
 */

$fasta_path = $jsc->getVariable('fasta', 'path');
$fasta_files = $jsc->getVariable('fasta', 'files');

$ui_ul = $jsc->getVariable('peptide', 'length_upper_limit', 53);
$ui_ll = $jsc->getVariable('peptide', 'length_lower_limit', 6);
$ui_um = $jsc->getVariable('peptide', 'mass_upper_limit', 6000);
$ui_lm = $jsc->getVariable('peptide', 'mass_lower_limit', 400);
$ui_mc = $jsc->getVariable('peptide', 'missed_clevage_max', 2);
$ui_rx = $jsc->getVariable('peptide', 'enzyme_regex', '/.*?[KR\#]/');

/*
 * cli flag for decoy database
 */
$is_decoy = $jsc->getVariable('peptide', 'decoy', false);

$dig->setLengthLimits($ui_ll, $ui_ul);
$dig->setMassLimits($ui_lm, $ui_um);
$dig->setMissedClevageMax($ui_mc);
$dig->setEnzymeRegex($ui_rx);

function ipc($p, $n = 0, $dig, $is_decoy)
{
    $seq = $p['seq'];
    $out = '';
    
    // $pep = new Digest();
    
    if (strlen($seq) <= 5)
        return false;
    
    if (! is_false($is_decoy)) {
        $p['name'] .= "_decoy";
        preg_match_all('/\w/', $seq, $rseq);
        shuffle($rseq[0]);
        $seq = array_tostring($rseq[0], '', '');
        unset($rseq);
    }
    
    $out .= '<sphinx:document id="' . $n . '">' . PHP_EOL;
    $out .= '<protein_source>' . $p['source'] . '</protein_source>' . PHP_EOL;
    $out .= '<protein_name>' . $p['name'] . '</protein_name>' . PHP_EOL;
    $out .= '<protein_id>' . $p['id'] . '</protein_id>' . PHP_EOL;
    $out .= '<protein_desc>' . $p['desc'] . '</protein_desc>' . PHP_EOL;
    $out .= '<organism>' . $p['org'] . "</organism>" . PHP_EOL;
    $out .= '<sequence>' . $seq . "</sequence>" . PHP_EOL;
    $out .= '<content>';
    
    $peptides = $dig->getPeptides($seq);
    if (is_false($peptides))
        return false;
    
    $peptides = preg_replace("/I/", "L", $peptides);
    
    $out .= array_tostring($peptides, ' ', '');
    
    $out .= "</content>" . PHP_EOL;
    $out .= '</sphinx:document>' . PHP_EOL . PHP_EOL;
    return $out;
}

$head = '<?xml version="1.0" encoding="utf-8"?>
        <sphinx:docset>
        <sphinx:schema>
        <sphinx:field name="content" />
        <sphinx:attr name="protein_source" type="string" />
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

foreach ($fasta_files as $file) {
    
    $handle = fopen($fasta_path . $file, "r");
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
                    $parsed = $fst->parseFastaEntry($entry);
                    $parsed['source'] = str_replace(".fasta", "", $file);
                    echo ipc($parsed, $n, $dig, $is_decoy);
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
        $parsed = $fst->parseFastaEntry($entry);
        $parsed['source'] = str_replace(".fasta", "", $file);
        echo ipc($parsed, $n, $dig, $is_decoy);
        
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
