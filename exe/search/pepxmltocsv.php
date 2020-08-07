#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
converts a proteomics pepxml to csv";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxyIO\File\Delim\WriteDelim;
use ProxyTime\ProgressTimer;

$cli = new Cli();
$csv = new WriteDelim();
$prg = new ProgressTimer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('RC v1.9.191014');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('...');

/*
 * load the variables
 * flag : default-value : description
 */
$path = $cli->getVar('p', "./", 'path to data file');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

$files = list_files($path, $file_regex . '.*\.xml$', TRUE);

foreach ($files as $file) {

    $data = [];
    $cli->message('read', basename($file));

    $xml = new SimpleXMLElement(file_get_contents($file));

    $spectrum_queries = $xml->msms_run_summary->spectrum_query;
    $prg->start("  parse", count($spectrum_queries));

    foreach ($spectrum_queries as $spectrum_query) {

        $prg->barPercent();
        /*
         * get the spectrum values
         */
        $spectrum = [];
        foreach ($spectrum_query->attributes() as $a => $v) {
            $spectrum[$a] = preg_replace("/\,/", '', preg_replace("/\"|\'/", ' ', $v->__toString()));
        }

        foreach ($spectrum_query->search_result as $search_result) {
            foreach ($search_result->search_hit as $search_hit) {

                /*
                 * get the peptide hit values
                 */
                $hit = [];
                foreach ($search_hit->attributes() as $a => $v) {
                    $hit[$a] = $v->__toString();
                }

                /*
                 * get the ptm modified peptide
                 */
                $hit['modified_peptide'] = "";
                if (!empty($search_hit->modification_info)) {
                    $hit['modified_peptide'] = $search_hit->modification_info->attributes()->__toString();
                }

                /*
                 * get the peptide hit score values
                 */
                $score = [];
                foreach ($search_hit->search_score as $search_score) {

                    $score[(string) $search_score['name']] = (string) $search_score['value'];
                }
                /*
                 * assemble the table
                 */
                $data[] = array_merge($spectrum, $hit, $score);
            }
        }
    }

    $table = array_rowtocol($data);
    $csv->writeCsvFile($path . basename($file, ".xml") . ".csv", $table);
}

?>