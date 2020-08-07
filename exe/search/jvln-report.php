#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
convert json result files to CSV";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\File\Delim\WriteDelim;
use ProxyTime\Timer;
use ProxyIO\Cli;

$cli = new Cli();
$csv = new WriteDelim();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v1.1.190916');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('-p /home/scbi/data -f som.*json -db yeast');
$cli->Help->setExample('-db yeast -nph 5');

/*
 * load the variables
 * flag : default-value : description
 */
$path = $cli->getVar('p', "./", 'path to data file');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');

$ui_methods = $cli->getVar('m', 'pep+pro', '<methods> default = scans+peptides+proteins
 -- this notation strings together computational modes
 -- with partial spelling accepted (i.e. sca+pep+pro)
 -- 
 -- methods:
 -- scans       stats on every scan collected
 -- peptides    stats on every scan collected with peptide hits
 -- ranker      defines the possible protein-group ranking options');

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

$cli->header("BEGIN");

foreach ($files as $file) {

    $tmr = new Timer();
    $data = json_decode(file_get_contents($file), TRUE);

    $cli->flush("f:" . basename($file, '.json'), 'read');

    /*
     * display the ranker
     */
    if (preg_match("/rank[er]*/", $ui_methods)) {

        $out_ranker = [];

        $data_prot = $data['proteins'];
        foreach ($data_prot as $data_prot_data) {
            $out_ranker_id = preg_replace("/.*_/", "", $data_prot_data['protein_name']);
            $out_ranker[$out_ranker_id] = $out_ranker_id;
        }

        print_r(array_values($out_ranker));

        exit();
    }


    /*
     * run spectrum report
     */
    if (preg_match("/sc[ans]*|pep[tides]*/", $ui_methods)) {

        $report_type = 'scans';
        $data_scans = $data_out = json_scans($data);
        if (is_false($data_scans))
            $cli->message("error", "no scan data", TRUE);

        if (preg_match("/pep[tides]*/", $ui_methods)) {
            $report_type = 'peptides';
            $data_peptides = json_peptides($data);
            $data_perf = json_performance($data);
            
            if (! is_false($data_peptides)) {
                $data_out = table_merge($data_scans, $data_peptides, 'scan_id');
                if (! is_false($data_perf)) 
                    $data_out = table_merge($data_out, $data_perf, 'scan_id');
            }
        }

        echo PHP_EOL;
        print_table(table_summary($data_out), 100);
        
        $file_out = "./" . basename($file, '.json') . ".report-" . $report_type . ".csv";
        $csv->writeCsvFile($file_out, $data_out);
        $cli->flusheol("f:" . basename($file, '.json'), $tmr->timeinstr(TRUE));
    }

    /*
     * run spectrum report
     */
    if (preg_match("/mod[eling]*/", $ui_methods)) {

        $report_type = 'modeling';

//         if (! key_exists('jaevalen', $data)) {
//             $cli->flusheol("     fail", "no protein grouping");
//             continue;
//         }

//         $scan_data = $data['scans'];
//         $jvln_data = $data['jaevalen']['psms'];
//         $jvln_stat = $data['jaevalen']['stats'];

//         $data_out = [];

//         $n = 0;
//         foreach ($jvln_data as $i => $jvln_table) {
//             $jvln_table = table_invert($jvln_table);
//             foreach ($jvln_table as $r => $jvln_values) {
//                 $n ++;
//                 $data_out[$n] = [];

//                 $data_out[$n]['file'] = basename($data['file']['name']);
//                 $data_out[$n]['scan_id'] = $i;

//                 $scan_values = $scan_data[$i];
//                 unset($scan_values['peaks']);

//                 $data_out[$n] = array_merge($data_out[$n], $scan_values, $jvln_values);
//                 $data_out[$n]['hit_rank'] = $r + 1;
//             }
//         }

//         $data_out = table_invert($data_out);

//         $file_out = "./" . basename($file, '.json') . ".report-" . $report_type . ".csv";
//         $csv->writeCsvFile($file_out, $data_out);
//         $cli->flusheol("f:" . basename($file, '.json'), $tmr->timeinstr(TRUE));
    }
}

$cli->header("COMPLETE");

?>