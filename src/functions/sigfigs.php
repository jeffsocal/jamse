<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function sigfigs($val, $n = 3)
{
    if ($val == 0)
        return $val;

    $l = log10(abs($val));

    if ($l >= 0)
        return truncate($val, $n);

    $n += ceil(- $l);

    return truncate($val, $n);
}

?>