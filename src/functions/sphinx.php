<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");


function indexer_protein_content($p, $n = 0)
{
    $out = '';
    $peps = digest_protein($p['seq'], TRUE);
    
    if (count($peps) == 0)
        return false;
    
    $out .= '<sphinx:document id="' . $n . '">' . PHP_EOL;
    $out .= '<protein_name>' . $p['name'] . '</protein_name>' . PHP_EOL;
    $out .= '<protein_id>' . $p['id'] . '</protein_id>' . PHP_EOL;
    $out .= '<protein_desc>' . $p['desc'] . '</protein_desc>' . PHP_EOL;
    $out .= '<organism>' . $p['org'] . "</organism>" . PHP_EOL;
    $out .= '<content>' . array_tostring($peps, ' ', '') . "</content>" . PHP_EOL;
    $out .= '</sphinx:document>' . PHP_EOL . PHP_EOL;
    
    return $out;
}


?>