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
use Jvln\Database\PeptideMatches;
use ProxyIO\File\Delim\Delim;
use ProxyIO\Exec;

$cli = new Cli();
$tmr = new Timer();
$csv = new Delim();
$exe = new Exec();

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

if (is_null($file_id))
    $cli->message("ERROR", "file id (-f) not set", TRUE);

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

/*
 * establish the connection to the database
 */
$cli->header("START");
$psms = new PeptideMatches($ui_server);

$tmr->start();
$cli->flush("get data");
$data = $psms->getPsms($file_id);
$file_dscore = "/var/jvlnse/local/dscore/" . $file_id . ".csv";
$csv->writeCsvFile($file_dscore, $data);
$cli->flusheol("get data", $tmr->timeinstr(TRUE));

$exe->setExitOnFail(TRUE);
$exe->run("Rscript /var/jvlnse/bin/jvln-dscore-csv.R --file=" . $file_dscore);
$cli->message("END Rscript", $tmr->timeinstr(TRUE));

$table_dscore = $csv->read($file_dscore);
$table_dscore = table_keepcols($table_dscore, [
    'pk',
    'psm_dscore',
    'hit_rank'
]);
$psms->updatePsms($table_dscore);
$exe->deleteFile($file_dscore);
$cli->header("COMPLETE");

?>