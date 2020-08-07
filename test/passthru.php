<?php
use ProxyIO\Cli;
use ProxyTime\Timer;
use Jvln\System\Passthru;

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * fragment
 *
 * SYNOPSIS
 * fragment -p <peptide_sequence>
 *
 * DESCRIPTION
 * display a table of fragment peaks for unit testing
 *
 * COMMAND LINE
 *
 * -p : <peptide_sequence> e.g. SAMPLER
 *
 * EXAMPLE
 *
 * fragment -p SAMPLE
 */

set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$tmr = new Timer();
$exe = new Passthru();

$pids = $cli->getVariable('p');

$table = $exe->monitorPIDs(explode(" ", $pids));

print_table($table);

?>
