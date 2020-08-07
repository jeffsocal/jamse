<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use function BenTools\CartesianProduct\cartesian_product;
use Jvln\Engines\TandemHash;

function dropByInt($x, $int)
{
    return $x <= $int;
}

function dropByPreMz($x, $pre)
{
    return $x >= $pre - 2 & $x <= $pre + 2;
}

?>