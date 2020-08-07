<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Fragments;

$frg = new Fragments();

if( is_null($h = opts('h')) )
    $h = 'axgn';

echo PHP_EOL;

# new test ########################################
$tmr = new Timer();
print_message("hash", $h, '.', 30);



print_message("integer", base26_decode($h), '.', 30);

print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);


?>