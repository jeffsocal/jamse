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
use Jvln\Database\Performance;

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
$ui_variable = $cli->getVar('v', 'psm_fcorr', 'variable to run TD FDR on');
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
// $perf = new Performance($ui_server);

if (! is_false($data = $psms->fileFDRExist($file_id, $ui_variable)))
    $cli->message("file has been FDR'd", $data['file_id'] . " n:" . $data['count'], TRUE);

$tmr->start();
$cli->flush("get data");
/*
 * need rank 2 to determine the FDR
 */
$data = $psms->getPsms($file_id, 4);
$file_fdr = "/var/jvlnse/local/fdr/" . $file_id . ".csv";
$csv->writeCsvFile($file_fdr, $data);
$cli->flusheol("get data", $tmr->timeinstr(TRUE));
/*
 * run the external R script
 */
$exe->setExitOnFail(TRUE);
$exe->run("Rscript /var/jvlnse/bin/jvln-perf-mm-csv.R --metric=" . $ui_variable . " --file=" . $file_fdr);
$cli->message("END Rscript", $tmr->timeinstr(TRUE));

/*
 * munge and save out the associated tables
 */
$file_sum = preg_replace("/\.csv/", "_summary.json", $file_fdr);
$summary_fdr = json_decode(file_get_contents($file_sum), TRUE);
$table_fdr = $csv->read($file_fdr);
$table_fdr = table_dropcols($table_fdr, [
    'fdr_metric'
]);

$table_fdr = table_renamecol($table_fdr, 'pk', 'psm_pk');

/*
 * array_map fails here for some reason
 */
for ($i = 0; $i < table_length($table_fdr); $i ++) {
    $table_fdr['fdr_metric'][$i] = $ui_variable;
    $table_fdr['fdr_qvalue'][$i] = sigfigs($table_fdr['fdr_qvalue'][$i], 4);
    $table_fdr['fdr_prob'][$i] = sigfigs($table_fdr['fdr_prob'][$i], 4);
    $table_fdr['fdr_eprob'][$i] = sigfigs($table_fdr['fdr_eprob'][$i], 4);
}

$perf = new Performance();
$perf->putPerfData($table_fdr, $summary_fdr);
$exe->deleteFile($file_fdr);
$exe->deleteFile($file_sum);
$cli->header("COMPLETE");

?>