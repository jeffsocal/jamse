#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
identifies and associates all proteins from which a given peptide
originated";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use Jvln\Database\Experiments;
use ProxyTime\ProgressTimer;
use ProxyTime\Timer;
use ProxyIO\Cli;

// $mgf = new ReadMGF("/data/SCBI/20171201_DBS_20Samples/mgf/20171201_SoCal_Samples18-9.mgf");

$cli = new Cli();
$pgt = new ProgressTimer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v1.1.190521');
$cli->Help->setUsage('-f <file_regex> -e <experiment_name>');
$cli->Help->setUsage('-f <file_regex> -e <experiment_name> -a <action>');

/*
 * load the variables
 * flag : default-value : description
 */
$ui_expr = $cli->getVar('w', '', 'experiment name or file name or file id');
$ui_actn = $cli->getVar('a', 'lookup', 'action [lookup, remove, create]');
$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');

$exps = new Experiments($ui_server);

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

if (! stristr("lookup, remove, create", $ui_actn))
    $cli->message("user input error", "action string not a valid choice [lookup, remove, create]", TRUE);

if ($ui_actn == 'lookup') {
    $table = $exps->getTableByMatch($ui_expr);
    if (is_false($table))
        $pgt->cli->message("ERROR:", "no files found for that group", TRUE);

    echo PHP_EOL;
    print_table($table, table_length($table));
}

if ($ui_actn == 'remove') {
    $exps->deletePSMSByMatch($ui_expr);
}

$pgt->header("COMPLETE");

?>