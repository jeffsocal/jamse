#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
searches the jaevalen engine for matching peptides based on the
precursor and fragments masses";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\Timer;
use ProxyIO\Cli;
use Jvln\Database\Metrics;

$cli = new Cli();
$tmr = new Timer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.9.190917');
$cli->Help->setUsage('-f <file_id>');

/*
 * load the variables
 * flag : default-value : description
 */
$ui_test = $cli->getVar('t', false, 'flag for testing, input the scan N to test');
$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');
$file_id = $cli->getVar('f', NULL, 'file id in database');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

if (is_null($file_id))
    $cli->message("ERROR", "file id (-f) not set", TRUE);
/*
 * establish the connection to the database
 */
$cli->header("START");
/*
 * get scan data
 */
$mtrx = new Metrics($ui_server);

if (! is_false($data = $mtrx->fileMetricsExist($file_id)))
    $cli->message("file has been assessed", $data['file_id'] . " n:" . $data['count'], TRUE);

$tmr->start();
$cli->flush("get data", '...');

$file_pk = $mtrx->getFileIDPK($file_id);
$data = $mtrx->getScanMetrics($file_id);

if (is_false($data))
    $cli->flusheol("ERROR", "no data", TRUE);

$data = table_dropcols($data, [
    'scan_pk',
    'file_id',
    'title'
]);
$table = table_summary($data);

$table_z = [];
foreach (table_invert($data) as $n => $row) {
    $z = min(3, $row['pre_z']);
    $table_z[$z][] = $row;
}

foreach ($table_z as $z => $zdata) {
    $zdata = table_invert($zdata);
    $zdata = table_dropcols($zdata, [
        'pre_z',
        'int_rsd',
        'int_rad',
        'n_peaks',
        't_peaks'
    ]);
    $zdata = table_summary($zdata);
    $z = array_fill(0, table_length($zdata), $z);
    $zdata['column'] = array_map(function ($x, $z) {
        return $x . "_z" . $z;
    }, $zdata['column'], $z);
    $table = table_bind($table, $zdata);
}

/*
 * get psm data
 */
$data = $mtrx->getPsmMetrics($file_id);
$data = table_dropcols($data, [
    'scan_pk',
    'file_id',
    'title',
    'pre_z',
    'psm_score'
]);
$table = table_bind($table, table_summary($data));

$table = table_dropcols($table, [
    'unique'
]);
$table = table_renamecol($table, 'column', 'metric');
$mtrx->putMetricData($table, $file_pk);
$cli->flusheol("get data", "finished " . $tmr->timeinstr(TRUE));
$cli->header("COMPLETE");

?>