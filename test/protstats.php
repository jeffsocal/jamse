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
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$tmr = new Timer();
$dig = new Digest();
$jsc = new JsonConfig("/var/jvlnse/uniprot/");

$rm = "\!|\@|\#|\$|\%|\^|\&|\*|\(|\)|\[|\]|\{|\}|\<|\>";

/*
 * get INI definitions
 */
if ($jsc->getConfigFile() == false)
    $cli->message('no config file found', "", TRUE);

/*
 * SET THE DATABASE
 */

$fasta_path = $jsc->getVariable('fasta', 'path', '/data/UniProt/');
$fasta_files = $jsc->getVariable('fasta', 'files', 'uniprot_sprot.fasta');

$ui_ul = $jsc->getVariable('peptide', 'length_upper_limit', 53);
$ui_ll = $jsc->getVariable('peptide', 'length_lower_limit', 6);
$ui_um = $jsc->getVariable('peptide', 'mass_upper_limit', 6000);
$ui_lm = $jsc->getVariable('peptide', 'mass_lower_limit', 400);
$ui_mc = $jsc->getVariable('peptide', 'missed_clevage_max', 2);
$ui_rx = $jsc->getVariable('peptide', 'enzyme_regex', '/.*?[KR\#]/');

/*
 * cli flag for decoy database
 */

$dig->setLengthLimits($ui_ll, $ui_ul);
$dig->setMassLimits($ui_lm, $ui_um);
$dig->setMissedClevageMax($ui_mc);
$dig->setEnzymeRegex($ui_rx);

$n = 0;
$entry = FALSE;

$out = [];
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
                    $entry = $fst->parseFastaEntry($entry);
                    
                    $seqs = preg_split("//", $entry['seq']);
                    
                    if (! key_exists($entry['org'], $out))
                        $out[$entry['org']] = array_fill_keys(range('A', 'Z'), 0);
                    
                    for ($i = 1; $i < count($seqs)-1; $i ++) {
                        $out[$entry['org']][$seqs[$i]] ++;
                    }
                    
                }
                
                $entry = trim($line, '>') . " ";
            } else {
                $entry .= $line;
            }
        }
        
        $seqs = preg_split("//", $entry['seq']);
        
        if (! key_exists($entry['org'], $out))
            $out[$entry['org']] = array_fill_keys(range('A', 'Z'), 0);
            
            for ($i = 1; $i < count($seqs)-1; $i ++) {
                $out[$entry['org']][$seqs[$i]] ++;
            }
        
        /*
         * get the last one
         */
        $n ++;
        $entry = preg_replace("/" . $rm . "/", "", $entry);
        $entry = preg_replace("/\//", "-", $entry);
        $entry = $fst->parseFastaEntry($entry);
        
        fclose($handle);
    } else {
        // error opening the file.
    }
}

file_put_contents("/var/jvlnse/uniprot_aa_accounting.json", json_encode($out));

?>
