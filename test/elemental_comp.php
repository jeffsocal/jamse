<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxySci\Atomic\Mass;
use ProxyTime\Timer;

$tmr = new Timer();
$atom = new Mass();
echo PHP_EOL;


echo $atom->getMolecularWeight("C16H34NO6PH");

echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

?>