<?php

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

set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\Timer;
use ProxyIO\Cli;

$cli = new Cli();

$n = $cli->getVariable('n', 123);

echo PHP_EOL;

# new test ########################################
$tmr = new Timer();
print_message("number", $n, '.', 30);

print_message("base26", base26_encode($n), '.', 30);

print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);


?>