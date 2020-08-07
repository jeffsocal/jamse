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
    
    if (! preg_match("/.pepXML$/", $file))
        continue;
    
    $data = [];
    print_message("file", $file, '.', 30);
    
    $xml = new SimpleXMLElement(file_get_contents($path . $file));
    
    foreach ($xml->msms_run_summary->spectrum_query as $spectrum_query) {
        
        /*
         * get the spectrum values
         */
        $spectrum = [];
        foreach ($spectrum_query->attributes() as $a => $v) {
            $spectrum[$a] = $v->__toString();
        }
        
        foreach ($spectrum_query->search_result as $search_result) {
            foreach ($search_result->search_hit as $search_hit) {
                
                /*
                 * get the peptide hit values
                 */
                $hit = [];
                foreach ($search_hit->attributes() as $a => $v) {
                    $hit[$a] = $v->__toString();
                }
                
                /*
                 * get the peptide hit score values
                 */
                $score = [];
                foreach ($search_hit->search_score as $search_score) {
                    
                    $score[(string) $search_score['name']] = (string) $search_score['value'];
                }
                
                /*
                 * get the peptide prophet values
                 */
                foreach ($search_hit->analysis_result->peptideprophet_result as $peptideprophet_result) {
                    
                    $score['peptide_prophet_p'] = (string) $peptideprophet_result['probability'];
                    
                    foreach ($peptideprophet_result->search_score_summary->parameter as $peptideprophet_parameter) {
                        
                        $score[(string) $peptideprophet_parameter['name']] = (string) $peptideprophet_parameter['value'];
                        
                    }
                }
                /*
                 * assemble the table
                 */
                $data[] = array_merge($spectrum, $hit, $score);
            }
        }
        
    }
    
    $table = array_rowtocol($data);
    $csv->writeCsvFile($path . $file . ".csv", $table);
}
print_message("TOTAL", $tmr->timeinstr(TRUE), '.', 30);

?>