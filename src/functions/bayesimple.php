<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function bayesimple(float $prior = 0.5, float $likelyhood = 0.5)
{
    $top = ($prior * $likelyhood);
    $btm = $top + ((1 - $prior) * (1 - $likelyhood));

    if ($btm == 0)
        return 0;

    return $posterior = $top / $btm;
}

?>