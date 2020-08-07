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
$cli->Help->setUsage('-f <file_id>');

/*
 * load the variables
 * flag : default-value : description
 */
$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');
$ui_index = $cli->getVar('db', 'uniprot', 'proteomic data base');
$file_id = $cli->getVar('f', NULL, 'file id in database');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

if (is_null($file_id))
    $cli->message("ERROR", "file id (-f) not set", TRUE);

if (! file_exists("/var/jvlnse/" . $ui_index))
    $cli->message("ERROR", "database (-db) does not exists on this machine", TRUE);

$exe->cli->header("JVLN SEARCH");
$cmd = "php ~/git/jvln-psm/exe/server/jvlnsql-search.php -db " . $ui_index . " -f " . $file_id;
if (is_false($exe->dohup($cmd, true)))
    $exe->cli->header("ERROR - EARLY TERMINATION", true);

$exe->cli->header("JVLN DSCORE");
$cmd = "php ~/git/jvln-psm/exe/server/jvlnsql-dscore.php -f " . $file_id;
if (is_false($exe->dohup($cmd, true)))
    $exe->cli->header("ERROR - EARLY TERMINATION", true);

$exe->cli->header("JVLN FDR");
$cmd = "php ~/git/jvln-psm/exe/server/jvlnsql-fdr.php -f " . $file_id;
if (is_false($exe->dohup($cmd, true)))
    $exe->cli->header("ERROR - EARLY TERMINATION", true);

?>