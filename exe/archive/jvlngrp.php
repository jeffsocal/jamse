#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
assembles peptides into protein groups by iteratively assigning 
peptides to the most representative protein, removing those 
peptides, and running the process again until all peptides have 
been assigned to a protein or group";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\Timer;
use ProxyIO\Cli;

$cli = new Cli();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.3.190521');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('-p /home/scbi/data -f som.*json');

/*
 * load the variables
 * flag : default-value : description
 */
$path = $cli->getVar('p', "./", 'path to data file');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

$files = list_files($path, $file_regex . '.*\.json$', TRUE);

$cli->header("JVLN START");
$cli->header("CONFIG");
$cli->message(" file path", $path);
$cli->message(" file regex", $file_regex);
$cli->message(" file n", count($files));

$jvln_table = [];
$cli->header("BEGIN");

foreach ($files as $file) {
    
    $prot_accounting_final = [];
    
    $data = json_decode(file_get_contents($file), TRUE);
    
    if (! key_exists('proteins', $data)) {
        $cli->message("f:" . basename($file, '.json'), "no jaevalen data found");
        continue;
    }
    
    if (! key_exists('performance', $data)) {
        $pgt->message("f:" . basename($file, '.json'), "no performance data found");
        continue;
    } else {
        $id_cutoff = $data['performance']['summary']['cutoff_value'];
        $id_metric = $data['performance']['summary']['cutoff_metric'];
    }
    
    $data_perf = $data['performance']['psms'];
    $data_pept = $data['homology'];
    $data_prot = $data['proteins'];
    
    $id_pep_probmax = [];
    $id_out_grp = [];
    
    foreach ($data_perf as $scan_id => $scan) {
        
        if (count($scan) > 0) {
            $peptide_score = $scan[$id_metric];
            $peptide_prob = $scan['psm_prob'];
            
            if ($peptide_score < $id_cutoff)
                continue;
            
            $peptide = $scan['peptide'];
            $peptide_bare = preg_replace("/I/", "L", getRootSequence($peptide));
            
            if (! key_exists($peptide_bare, $id_pep_probmax))
                $id_pep_probmax[$peptide_bare] = 0;
            
            $id_pep_probmax[$peptide_bare] = max($id_pep_probmax[$peptide_bare], $peptide_prob);
        }
    }
    
    $id_pep_probmax = array_intersect_key($id_pep_probmax, $data_pept);
    
    $total_peptides = count($id_pep_probmax);
    
    $timer = new Timer();
    
    while (count($id_pep_probmax) > 0) {
        
        /*
         * PROTEIN DATA
         */
        
        $stdout = ' proteins: ' . count($data_prot);
        $stdout .= ' peptides: ' . count($id_pep_probmax);
        $stdout .= ' of ' . $total_peptides;
        $stdout .= ' -- time: ' . $timer->timeinstr(true);
        
        $cli->flush("f:" . basename($file, '.json'), $stdout);
        
        $keep_peptides = [];
        $prot_accounting = [];
        
        foreach ($data_prot as $protein_id => $protein_data) {
            
            $this_probs = array_intersect_key($id_pep_probmax, array_flip($protein_data['peptides']));
            $this_peptides = array_keys($this_probs);
            
            if (count($this_peptides) == 0) {
                unset($data_prot[$protein_id]);
                continue;
            }
            
            /*
             * calculate the protein group score
             */
            sort($this_peptides);
            $this_prot['id'] = $protein_id;
            $this_prot['hits'] = count($this_peptides);
            $this_prot['score'] = array_probability($this_probs) * count($this_peptides);
            $this_prot['protein_group'] = 'gid' . substr(md5(array_tostring($this_peptides, ' ', '')), 3, 8);
            
            // if (count($this_probs) > 5) {
            // print_r($this_probs);
            // echo "\n";
            // print_r(array_sum($this_probs));
            // echo "\n";
            // print_r(array_product($this_probs));
            // echo "\n";
            // exit();
            // }
            /*
             * keep a lookup table of peptides (as keys) that belong to the group
             */
            $keep_peptides[$this_prot['protein_group']] = $this_probs;
            $prot_accounting[] = $this_prot;
        }
        
        if (count($prot_accounting) == 0) {
            break;
        }
        
        // reshape the data
        $prot_accounting = array_rowtocol($prot_accounting);
        
        // filter by score max
        $score_max = array_max($prot_accounting['score']);
        $prot_accounting = table_droprows_lt($prot_accounting, 'score', $score_max);
        $prot_group = $prot_accounting['protein_group'][0];
        
        $prot_accounting = array_rowtocol($prot_accounting);
        
        /*
         * pull out non matching group_ids
         */
        foreach ($prot_accounting as $i => $i_data) {
            if ($i_data['protein_group'] != $prot_group)
                unset($prot_accounting[$i]);
        }
        $prot_accounting = array_rowtocol($prot_accounting);
        
        /*
         * for the results
         */
        $prot_peptides = $keep_peptides[$prot_group];
        $prot_prob = array_probability($prot_peptides);
        $prot_accounting_final[] = [
            'protein_group' => $prot_group,
            'protein_score' => $prot_prob * count($prot_peptides),
            'protein_prob' => $prot_prob,
            'proteins' => $prot_accounting['id'],
            'peptides' => $prot_peptides
        ];
        
        /*
         * remove peptides from the pool for the next iteration
         */
        $id_pep_probmax = array_diff_key($id_pep_probmax, $keep_peptides[$prot_group]);
    }
    
    /*
     * one last print message
     */
    $stdout = ' proteins: ' . count($data_prot);
    $stdout .= ' peptides: ' . count($id_pep_probmax);
    $stdout .= ' of ' . $total_peptides;
    $stdout .= ' -- time: ' . $timer->timeinstr(true);
    
    $cli->flusheol("f:" . basename($file, '.json'), $stdout);
    
    $prot_accounting_final = array_rowtocol($prot_accounting_final);
    $prot_accounting_final = table_sort($prot_accounting_final, 'protein_score', 'desc', 'number');
    $prot_accounting_final = table_dropcols($prot_accounting_final, 'protein_score');
    $prot_accounting_ids = $prot_accounting_final['protein_group'];
    unset($prot_accounting_final['protein_group']);
    $prot_accounting_final = array_rowtocol($prot_accounting_final);
    
    $prot_accounting_final = array_combine($prot_accounting_ids, $prot_accounting_final);
    
    $data['protein_groups'] = $prot_accounting_final;
    
    file_put_contents($file, json_encode($data));
}

$cli->header("COMPLETE");

?>