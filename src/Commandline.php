#!/usr/bin/env php
<?php

// ###############################################################################
// Copyright (2016) SoCal Bioinformatics Inc. All rights reserved.
// This script is the confidential and proprietary product of SoCal
// Bioinformatics Inc. Any unAuthenticated reproduction or transfer of
// the contents herein is strictly prohibited.
//
// ###############################################################################
// AUTH: Jeff Jones | SoCal Bioinformatics Inc
// DATE: 2017.01.10
// OWNER: SoCal Bioinformatics, Inc., Glendale CA
// PROJECT:
// DESC: =
// ###############################################################################
//
/*
 * IMPORT MAIN CLASSES
 */
set_include_path(__DIR__ . '/../');

require_once get_include_path() . 'vendor/autoload.php';
require_once get_include_path() . 'src/Autoload/Functions.php';

function opts($var = '')
{
    $args = $_SERVER['argv'];
    if ($key = array_search('-' . $var, $args)) {
        return $args[$key + 1];
    } else {
        return NULL;
    }
}

/*
 * MAIN
 */
$cmd = $argv;
$exe = $cmd[1];

/*
 * CONTROLLER
 */
$_ENV['PID'] = uniqueId();

if (! file_exists(get_include_path() . 'exe/' . $exe . '.php')) {
    systemError($exe . " is an unknown program");
}

if (preg_match("/-h|--help/", array_tostring($_SERVER['argv']))) {
    
    $con = file_get_contents(get_include_path() . 'exe/' . $exe . '.php');
    
    if (preg_match("/(?<=MANUAL).+\n[\s\*\n\w\-\.\>\<\'\"\:\,]+/", $con, $man) != FALSE) {
        $man[0] = preg_replace("/(?<=\n)\s+\*/", " ", $man[0]);
        echo $man[0] . PHP_EOL . PHP_EOL;
    }
    
    exit();
}

require_once get_include_path() . 'exe/' . $exe . '.php';

echo PHP_EOL;

?>
