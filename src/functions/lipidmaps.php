<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use ProxySci\Omics\Proteomics\Fasta;

function parse_sdf($e, $f, $n = 0, $r = FALSE)
{
    // if ($n != 11027705)
    // return;
    $continue = FALSE;
    if (preg_match("/" . $r . "/", $e))
        $continue = TRUE;
    
    if ($continue == FALSE)
        return;
    
    $fst = new Fasta(preg_replace("/[\_0-9].+/", "", $f));
    
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