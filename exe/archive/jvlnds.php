#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
Applies a modeled discriminant score to jvln search results
with an estimated AUC of 0.96.

-- [dscore] the value of the discriminat function applied to the PSM";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\ProgressTimer;
use MathPHP\Probability\Distribution\Discrete\Poisson;
use ProxyIO\Cli;

$cli = new Cli();
$pgt = new ProgressTimer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.1');
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

$pgt->header("JVLN START");
$pgt->header("CONFIG");
$pgt->message(" file path", $path);
$pgt->message(" file regex", $file_regex);
$pgt->message(" file n", count($files));

function discriminant_model($array)
{
    /*
     * linear regression discriminat score
     *
     * estimated AUC 0.855 10x10 FCV
     * /* score 1.190338
     * /* int_iqr 3.750421
     * /* a_bins -0.410630
     * /* abs_pre_mma -0.500655
     * /* n_match -0.913862
     * /* r_overlap 19.979971
     * /* abs_mma_mean -0.546755
     */
    $dscore = $array['score'] * 1.190338;
    $dscore += $array['int_iqr'] * 3.750421;
    $dscore += $array['a_bins'] * - 0.410630;
    $dscore += abs($array['pre_mma']) * - 0.500655;
    $dscore += $array['psm_o'] * - 0.913862;
    $dscore += $array['psm_r'] * 19.979971;
    $dscore += abs($array['psm_mma']) * - 0.546755;
    $dscore = $dscore + 0.4;
    
    return $dscore;
}

$pgt->header("BEGIN");
foreach ($files as $file) {
    
    $data = json_decode(file_get_contents($file), TRUE);
    
    if (! key_exists('jaevalen', $data)) {
        $pgt->message("f:" . basename($file, '.json'), "no jaevalen data found");
        continue;
    }
    
    $data_jvln = $data['jaevalen']['psms'];
    $stat_jvln = $data['jaevalen']['stats'];
    $data_scan = $data['scans'];
    
    $pgt->start(" f: " . basename($file), count($data_jvln));
    foreach ($data_jvln as $scan_id => $hit_data) {
        
        $pgt->barPercent();
        $hit_data = array_rowtocol($hit_data);
        $scan_data = $data_scan[$scan_id];
        $stat_data = $stat_jvln[$scan_id];
        foreach ($hit_data as $i => $hit) {
            $ds_array = discriminant_model(array_merge($hit, $scan_data, $stat_data));
            $hit_data[$i] = array_merge($hit, $ds_array);
        }
        $hit_data = array_rowtocol($hit_data);
        $hit_data = table_sort($hit_data, 'dscore', 'desc', 'number');
        
        $data['jaevalen']['psms'][$scan_id] = $hit_data;
    }
    
//     file_put_contents($file, json_encode($data));
}

$pgt->header("COMPLETE");

?>