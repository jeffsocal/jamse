#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
generates a fasta file composed of a subset of proteins
based on the regular expression provided. The path to
the fasta file is determined from the molecular.ini file.";

use ProxyIO\Cli;
use ProxyIO\File\Write;
use ProxyTime\ProgressTimer;
use ProxySci\Omics\Proteomics\Fasta;

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$ptm = new ProgressTimer();
$doc = new Write();
$fst = new Fasta();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('RELEASE v1.0.19');
$cli->Help->setUsage('subfasta -f <fasta_files> -o <file_out_name> -r <regex>');
$cli->Help->setExample('-r saccharomyces -o saccharomyces_uniprot.fasta');

/*
 * load the variables
 * flag : default-value : description
 */
$files = $cli->getVar('f', 'uniprot_sprot.fasta', '<fasta_files> comma delimited');
$dest = $cli->getVar('o', 'sub.fasta', '<file_out_name>');
$regex = $cli->getVar('r', NULL, '<regex> the expression used to pull the subset');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

if (is_null($regex))
    $regex = '';

$regex = preg_replace("/\"|\'/", "", $regex);

$files = explode(",", $files);

$n = 0;
$entry = FALSE;

$ptm->cli->header('START');
$ptm->cli->message('fasta', array_tostring($files, ', ', ''));
$ptm->cli->message('out', $dest);
$ptm->cli->message('regex', $regex);

// $ptm->start('parsing fasta');
$ptm->start('parsing fasta');

foreach ($files as $file) {
    $handle = fopen($file, "r");
    if ($handle) {
        while (($line = fgets($handle)) !== false) {
            /*
             * process the line read.
             */
            if (preg_match('/^>/', $line)) {

                $ptm->barCount(1000, $n);
                // $ptm->barPercent();

                if (! is_false($entry)) {

                    if (preg_match('/' . $regex . '/i', $entry)) {
                        $n ++;
                        // $ptm->barCount(1000, $n);

                        $p = $fst->parseFastaEntry($entry);

                        $entry = ">sp|" . $p['id'] . "|";
                        $entry .= $p['name'] . " " . $p['desc'];
                        $entry .= " OS=" . $p['org'];
                        $entry .= PHP_EOL . chunk_split(str_shuffle($p['seq']), 60, PHP_EOL);

                        $doc->writeFile($dest, $entry, 'a');
                    }

                    // if ($n % 1000 == 0)
                    // $ptm->printMessage($n, $tmr->timeinstr());
                }
                $entry = $line;
            } else {
                $entry .= $line;
            }
        }
        /*
         * get the last one
         */
        $n ++;
        if (preg_match('/' . $regex . '/', $entry))
            $doc->writeFile($dest, $entry, 'a');

        fclose($handle);
    } else {
        // error opening the file.
    }
}
$ptm->barCount(1000, $n);
echo PHP_EOL;
$ptm->cli->header('FINISHED');

// UniProtKB 125,355,233
// - UniProtKB/Swiss-Prot 558,125 ->
// - UniProtKB/TrEMBL 124,797,108
// UniRef100 154,752,138
// UniRef90 78,915,455
// UniRef50 32,365,088
// UniParc 225,932,239

?>