#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
converts a proteomics featureXML to csv";

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

$files = list_files($path, $file_regex . '.*\.featureXML$', TRUE);

$data = [];
foreach ($files as $file) {

    $data = [];
    $cli->message('read', basename($file));

    $xml = new SimpleXMLElement(file_get_contents($file));

    $feature_list = $xml->featureList;

    $prg->start("  parse", count($feature_list->feature));

    foreach ($feature_list->feature as $feature) {

        $prg->barPercent();
        /*
         * get the spectrum values
         */
        $feature_data = [];

        foreach ($feature->position as $feature_value) {
            $name = 'elution';
            if (sizeof($feature_data) == 1)
                $name = 'masscharge';
            $feature_data[$name] = (string) $feature_value->__toString();
        }

        $feature_data['charge'] = $feature->charge->__toString();
        $feature_data['intensity'] = $feature->intensity->__toString();

        $data[] = $feature_data;
    }

    $table = table_invert($data);

    $csv->writeCsvFile($path . basename($file, ".xml") . ".csv", $table);
}

?>