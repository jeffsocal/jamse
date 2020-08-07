<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
use ProxyTime\Timer;


$set = range(1e10,2e10,10000);

echo PHP_EOL;

# new test ########################################
$tmr = new Timer();
foreach($set as $val){
    $a = numberToAlphaHash($val);
}
echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);

echo PHP_EOL;

# new test ########################################
$tmr = new Timer();
foreach($set as $val){
    $a = crc32($val);
}
echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);


?>