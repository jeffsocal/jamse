#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
converts a proteomics MGF to a JSON file that has been filtered
to just the peaks of interest, also calculates the spectral p-value
(an estimate for the likelyhood the spectra contains enough information
to generate a credable PSM) and the N (estimated) corrected e-value,
based on the same Possion probability distribution proposed by Geer, L. Y. et al. 
J. Proteome Res. 3, 958–964 (2004).";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyIO\Cli;
use ProxySci\MassSpec\Fragmentation\ReadMGF;
use ProxyTime\ProgressTimer;
use MathPHP\Probability\Distribution\Discrete\Poisson;

$cli = new Cli();
$prg = new ProgressTimer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.9.190917');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('-p ./ -m int -d 4 -n 30');

/*
 * load the variables
 * flag : default-value : description
 */
$path = $cli->getVar('p', "./", 'path to data file');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');

$peaks_max_n = $cli->getVar('n', 30, 'max number of peaks');
$peaks_filter_da = $cli->getVar('d', 4, 'dalton plus-minus to filter from top intensity');
$scan_limit = $cli->getVar('l', NULL, 'limit number of scans to process');
$method = $cli->getVar('m', 'int', '<methods> default = intensity
 -- this notation strings together computational modes
 -- with partial spelling accepted (i.e. int or iso)
 -- 
 -- methods:
 -- intensity  filter based on area subtraction from top peaks
 -- isotopic   determine the spectral monoisotopic peaks');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

if (preg_match("/iso[topic]*/", $method))
    $method = 'isotopic determination';

if (preg_match("/int[ensity]*/", $method))
    $method = 'peak intensity filtering';

$files = list_files($path, $file_regex . '.*\.mgf$', TRUE);

$prg->cli->header("MGF-JSON START");
foreach ($files as $file) {

    /*
     * read in the data
     */
    $mgf = new ReadMGF($file);
    $data = $mgf->getData();

    $mgf_n = count($data);
    $pad_n = ceil(log10($mgf_n));

    if (! is_null($scan_limit)) {
        $data = array_slice($data, 0, $scan_limit);
        $mgf_n = count($data);
        $prg->cli->message(' limiting to', $mgf_n);
    }

    /*
     * form the methods output
     */
    $file_atts = [
        'name' => $file,
        'scans_n' => $mgf_n,
        'method' => $method,
        'peaks_max_n' => $peaks_max_n,
        'peaks_filter_da' => $peaks_filter_da
    ];

    if (preg_match("/iso[topic]*/", $method))
        unset($file_atts['peaks_filter_da']);

    /*
     * construct the base name
     */
    $file = basename($file, '.mgf');
    $file_out = "./" . $file . ".json";

    $json_out = [];
    $prg->start($file . ':' . $mgf_n, $mgf_n);

    $n_break = 0;
    foreach ($data as $scan => $this_data) {

        $n_break ++;
        $prg->barPercent();

        $scan_id = 'sn' . str_pad($scan, $pad_n, '0', STR_PAD_LEFT);

        $pre_m = 0;
        $pre_z = 0;
        $title = 0;
        $rtsec = 0;

        if (key_exists('pepmass', $this_data))
            $pre_m = round(preg_replace("/\s+.*/", "", $this_data['pepmass']), 4);

        if (key_exists('charge', $this_data))
            $pre_z = preg_replace("/[\+\s+].*/", "", $this_data['charge']) * 1;

        /*
         * override the zeros
         */
        $pre_z = max(1, $pre_z) * 1;

        if (key_exists('title', $this_data))
            $title = preg_replace("/[\"\,]+/", "", $this_data['title']);

        if (key_exists('rtinseconds', $this_data))
            $rtsec = round($this_data['rtinseconds'], 2);

        $this_spectra = $this_data['spec'];
        /*
         * round the spectral values
         */ 
        $this_spectra['mz'] = array_map(function ($x) {
            return round($x, 3);
        }, $this_spectra['mz']);

        $this_spectra['int'] = array_map(function ($x) {
            return round($x, 2);
        }, $this_spectra['int']);

        /*
         * pull out the peaks
         */ 
        if (preg_match("/iso[topic]*/", $method))
            $peaks = spec_get_peaks_by_iso($this_spectra, $peaks_max_n);

        if (preg_match("/int[ensity]*/", $method))
            $peaks = spec_get_peaks_by_int2($this_spectra, $pre_m, $pre_z, $peaks_max_n, $peaks_filter_da);

        if (is_false($peaks))
            continue;

        $json_out[$scan_id]['title'] = $title;
        $json_out[$scan_id]['rt_sec'] = $rtsec;
        $json_out[$scan_id]['pre_mz'] = $pre_m;
        $json_out[$scan_id]['pre_z'] = $pre_z;

        foreach ($peaks as $var => $val) {
            $json_out[$scan_id][$var] = $val;
        }

        $json_out[$scan_id]['peaks'] = mysql_escape_string(gzdeflate($json_out[$scan_id]['peaks']));
        $json_out[$scan_id]['spectra'] = mysql_escape_string(gzdeflate($this_spectra));

        print_r($json_out[$scan_id]);
        exit;
        
        $n_mass = neutralMass($pre_m, $pre_z);
        $n_bins = averagineBins($pre_m, $pre_z);
        /*
         * expected number of peptides in the search space give a
         * database of ~8000 non-redundant peptides
         */
        $n_search = max(250, 2800 + $n_mass * - 0.8);
        /*
         * Poisson distribution mean
         */
        $p_mean = (2 * 1 * $peaks['n_peaks'] * $n_bins) / $n_mass;
        /*
         * Poisson distribution
         */
        $poisson = new Poisson($p_mean);
        $pmatch = abs(1 - $poisson->cdf($n_bins * $peaks['a_bins']));
        $json_out[$scan_id]['pmatch'] = sigfigs($pmatch, 4);
        $json_out[$scan_id]['ematch'] = sigfigs($pmatch * $n_search, 4);

        /*
         * if there is a non-float p-value in match fail out
         */
        if ($json_out[$scan_id]['pmatch'] . "" == 'NAN')
            $pgt->cli->flusheol("ERROR", "exception to spectral quality scan:" . sizeof($json_out), TRUE);
    }

    // /*
    // * pull out the scan numbers and use as array keys
    // */
    // $scan_keys = array_values($json_out['scan']);
    // $scan_keys = array_map(function ($x) {
    // return 'scan_' . str_pad($x, 6, 0, STR_PAD_LEFT);
    // }, $scan_keys);

    // unset($json_out['scan']);

    // foreach ($json_out as $col => $data) {
    // $json_out[$col] = array_combine($scan_keys, $data);
    // }

    /*
     * format the output
     */
    $json_out = [
        'file' => $file_atts,
        'scans' => $json_out
    ];

    /*
     * save the output
     * this is stored as long-data, i.e. as column-based arrays
     */
    file_put_contents($file_out, json_encode($json_out));
    $prg->message("save file", $file_out);
}

$prg->header("COMPLETE");

?>