<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

function digest_protein($sequence, $sub_il = FALSE, $lengthUpperLimit = 36, $lengthLowerLimit = 6)
{
    $sequence = preg_replace("/[\n\r\t\s\.\-\*]/", '', $sequence);
    
    if (is_true($sub_il))
        $sequence = preg_replace("/I/", "L", $sequence);
    
    /*
     * get all regular peptides
     */
    preg_match_all("/.*?[K|R]/", $sequence, $peptides);
    
    /*
     * get the last non-specific tryptic peptide and add it to the end of the list
     */
    preg_match_all("/.*?[K|R]/", strrev($sequence), $peptides_rev);
    
    //
    if (key_exists(0, $peptides_rev[0])) {
        $peptides_rev = strrev(preg_replace("/[K|R]/", "", $peptides_rev[0][0]));
        if ($peptides_rev != "")
            $peptides[0][] = $peptides_rev;
    }
    
    if (sizeof($peptides[0]) == 0) {
        /*
         * systemError("Low Peptide Count: $i");
         */
        $peptides = array(
            $sequence
        );
    } else {
        $peptides = $peptides[0];
    }
    
    //
    $i = sizeof($peptides);
    $concat_peptide = "";
    $array_peptides = array();
    
    $s = 0;
    for ($n = 0; $n < $i; $n ++) {
        $concat_peptide .= $peptides[$n];
        $stl = strlen($concat_peptide);
        if ($stl <= $lengthUpperLimit) {
            if ($stl >= $lengthLowerLimit) {
                $array_peptides[] = $concat_peptide;
                
                if ($n < 5) {
                    $array_peptides[] = preg_replace("/^M/", '', $concat_peptide);
                }
            }
            $s ++;
        } else {
            $concat_peptide = "";
            $n -= $s;
            $s = 0;
        }
        if ($n == $i - 1 and $s != 0) {
            $concat_peptide = "";
            $n -= ($s - 1);
            $s = 0;
        }
    }
    return $array_peptides;
}

?>