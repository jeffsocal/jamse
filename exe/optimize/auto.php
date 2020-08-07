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
use Jvln\Database\Experiments;

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

/*
 * grab all the file ids echo'd out by the mgf parser
 */
$exp = new Experiments();
$files = $exp->getExpr($ui_expr);

$file_ids = $files['file_id'];

/*
 * run through each unique file id searching and sql'n the data
 */
$exe->cli->header("START");
foreach ($file_ids as $file_id) {

    $cmd = "~/git/jvln-psm/exe/optimize/search.php -db " . $ui_index . " -f " . $file_id;
    if (is_false($exe->fedup($cmd, true)))
        $exe->cli->header("ERROR - EARLY TERMINATION", true);
}
$exe->cli->message("time", $tmr->timeinstr());
$exe->cli->header("COMPLETE");

?>