<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function array_probability($array)
{
    $probability = 1;
    foreach ($array as $prob) {
        $probability *= (1 - $prob);
    }
    return 1 - $probability;
}

?>