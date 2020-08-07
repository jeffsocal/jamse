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
use Jvln\System\Passthru;

$cli = new Cli();
$tmr = new Timer();
$exe = new Passthru();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.9.190917');
$cli->Help->setUsage('-f <file_regex>');

/*
 * load the variables
 * flag : default-value : description
 */
$ui_expr = $cli->getVar('x', '', 'experiment name');
$ui_mthd = $cli->getVar('m', 'mzml', 'data format type [mgf|mzml]');
$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');
$ui_index = $cli->getVar('db', 'uniprot', 'proteomic data base');
$ui_path = $cli->getVar('p', "./", 'path to data file');
$ui_regex = $cli->getVar('f', '', 'file regex as a means to filter');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

if ($ui_expr == '')
    $ui_expr = uniqid("exp");

if (! file_exists("/var/jvlnse/" . $ui_index))
    $cli->message("ERROR", "database (-db) does not exists on this machine", TRUE);

if ($ui_mthd == 'mgf') {
    $exe->cli->header("MGF-2-SQL");
    $cmd = "mgf2sql -p " . $ui_path . " -f " . $ui_regex . " -x " . $ui_expr . " -ip " . $ui_server;
    if (is_false($exe->fedup($cmd, true)))
        $exe->cli->header("ERROR - EARLY TERMINATION", true);
} elseif ($ui_mthd == 'mzml') {
    $exe->cli->header("MZML-2-SQL");
    $cmd = "mzml2sql -p " . $ui_path . " -f " . $ui_regex . " -x " . $ui_expr . " -ip " . $ui_server;
    if (is_false($exe->fedup($cmd, true)))
        $exe->cli->header("ERROR - EARLY TERMINATION", true);
} else {
    $exe->cli->message("ERROR", "datafile (-m) format type not recognized", true);
}
/*
 * grab all the file ids echo'd out by the mgf parser
 */
$stdout = $exe->getStdout();
preg_match_all("/mz[a-z0-9]{13}/", array_tostring($stdout, " ", ''), $file_ids);

/*
 * run through each unique file id searching and sql'n the data
 */
$exe->cli->header("START");
foreach ($file_ids[0] as $file_id) {

    $cmd = "jvlnsql-search -db " . $ui_index . " -f " . $file_id . " -ip " . $ui_server;
    if (is_false($exe->fedup($cmd, true)))
        $exe->cli->header("ERROR - EARLY TERMINATION", true);

    // $cmd = "jvlnsql-dscore -f " . $file_id;
    // if (is_false($exe->fedup($cmd, true)))
    // $exe->cli->header("ERROR - EARLY TERMINATION", true);

    foreach ([
        'psm_fcorr',
        'psm_score'
    ] as $fdr_metric) {
        $cmd = "jvlnsql-fdr -f " . $file_id . " -v " . $fdr_metric . " -ip " . $ui_server;
        if (is_false($exe->fedup($cmd, true)))
            $exe->cli->header("ERROR - EARLY TERMINATION", true);
    }

    $cmd = "jvlnsql-metrics -f " . $file_id . " -ip " . $ui_server;
    if (is_false($exe->fedup($cmd, true)))
        $exe->cli->header("ERROR - EARLY TERMINATION", true);
}
$exe->cli->message("time", $tmr->timeinstr());
$exe->cli->header("COMPLETE");

?>