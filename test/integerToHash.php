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

if( is_null($n = opts('n')) )
    $n = 1234;

echo PHP_EOL;

# new test ########################################
$tmr = new Timer();
print_message("number", $n, '.', 30);


$hash = $frg->massToHash($n);
$hint = $frg->hashToMass($hash);

print_message("hash", $hash, '.', 30);
print_message("b26d", base26_decode($hash), '.', 30);
print_message("hash-to-mass", $hint, '.', 30);

print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);


?>