<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use ProxySci\Omics\Proteomics\Fasta;

function parse_fasta($e, $f, $n = 0)
{
    $fst = new Fasta($f);
    
    $rm = "\!|\@|\#|\$|\%|\^|\&|\*|\(|\)|\[|\]|\{|\}|\<|\>";
    $e = preg_replace("/" . $rm . "/", "", $e);
    $e = preg_replace("/\//", "-", $e);
    
    $p = $fst->parseFastaEntry($e);
    
    $p['id'] = str_replace(".fasta", "", $f) . '_' . $p['id'];
    $out = indexer_protein_content($p, $n);
    
    if (is_false($out))
        return;
    
    // echo $e . PHP_EOL;
    echo $out;
}

?>