#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
convert json result files to CSV";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\File\Delim\WriteDelim;
use ProxyTime\Timer;
use ProxyIO\Cli;

$cli = new Cli();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v1.1.190916');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('-p /home/scbi/data -f som.*json -db yeast');

/*
 * load the variables
 * flag : default-value : description
 */
$path = $cli->getVar('p', "./", 'path to data file');
$file_new = $cli->getVar('o', "jvln", 'output file name');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

$files = list_files($path, $file_regex . '.*\.json$', TRUE);

$cli->header("JVLN START");
$cli->header("CONFIG");
$cli->message(" file path", $path);
$cli->message(" file regex", $file_regex);
$cli->message(" file n", count($files));

$cli->header("BEGIN");

$json_out = [];
$file_n = 0;
foreach ($files as $file) {
    $file_n++;
    $tmr = new Timer();
    $cli->flush("f:" . basename($file, '.json'), 'read');
    $data = json_decode(file_get_contents($file), TRUE);

    $file_id = 'f' . $file_n;
    $data['file']['file_id'] = $file_id;

    $json_out['file'][] = $data['file'];
    foreach ($data['scans'] as $scan_id => $scan_data) {
        $json_out['scans'][$file_id . $scan_id] = $scan_data;
    }
    $json_out['jaevalen']['parameters'][$file_id] = $data['jaevalen']['parameters'];
    foreach ($data['jaevalen']['stats'] as $scan_id => $scan_data) {
        $json_out['jaevalen']['stats'][$file_id . $scan_id] = $scan_data;
    }
    foreach ($data['jaevalen']['psms'] as $scan_id => $scan_data) {
        $json_out['jaevalen']['psms'][$file_id . $scan_id] = $scan_data;
    }

    $cli->flusheol("f:" . basename($file, '.json'), $tmr->timeinstr(TRUE));
}

file_put_contents($path."/../".$file_new."_combined_".date("ymdhms").".json", json_encode($json_out));

$cli->header("COMPLETE");

?>