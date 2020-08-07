<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "4096M");

use ProxyIO\File\Delim\WriteDelim;
use ProxyTime\Timer;

$csv = new WriteDelim();
$tmr = new Timer();

$path = opts('p');
$files = explode(",", opts('f'));

if ($files[0] == '') {
    $files = array_diff(scandir($path), array(
        '..',
        '.'
    ));
}

foreach ($files as $file) {
    
    if (! preg_match("/.protXML$/", $file))
        continue;
    
    $data_peptides = [];
    $data_proteins = [];
    print_message("file", $file, '.', 30);
    
    $xml = new SimpleXMLElement(file_get_contents($path . $file));
    
    foreach ($xml->protein_group as $protein_group) {
        
        $data_prot_group = [];
        foreach ($protein_group->attributes() as $a => $v) {
            $data_prot_group[$a] = $v->__toString();
        }
        unset($data_prot_group['probability ']);
        
        foreach ($protein_group->protein as $protein) {
            $data_prot = array_fill_keys([
                'protein_name',
                'probability',
                'n_indistinguishable_proteins',
                'percent_coverage',
                'group_sibling_id',
                'total_number_peptides',
                'total_number_distinct_peptides',
                'pct_spectrum_ids',
                'confidence',
                'pseudo_name',
                'subsuming_protein'
            ], '');
            
            foreach ($protein->attributes() as $a => $v) {
                $data_prot[$a] = $v->__toString();
            }
            unset($data_prot['unique_stripped_peptides']);
            
            $data_proteins[] = array_merge($data_prot_group, $data_prot);
            
            $data_pep = [];
            foreach ($protein->peptide as $peptide) {
                foreach ($peptide->attributes() as $a => $v) {
                    $data_pep[$a] = $v->__toString();
                }
                
                $data_peptides[] = array_merge($data_prot_group, $data_prot, $data_pep);
            }
        }
    }
    
    $data_peptides = array_rowtocol($data_peptides);
    $data_proteins = array_rowtocol($data_proteins);
    
    print_table($data_proteins, 1000);
    exit();
    
    $csv->writeCsvFile($path . $file . ".csv", $table);
}
print_message("TOTAL", $tmr->timeinstr(TRUE), '.', 30);

?>