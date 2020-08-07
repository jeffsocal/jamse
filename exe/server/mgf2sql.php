#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
converts a proteomics MGF to a SQL table that has been filtered
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
use Jvln\Database\Spectra;
use ProxyTime\Timer;
use MathPHP\Probability\Distribution\Continuous\LogNormal;

$cli = new Cli();
$prg = new ProgressTimer();
$tmr = new Timer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('RC v1.9.191014');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('-p ./ -m int -d 4 -n 30');

/*
 * load the variables
 * flag : default-value : description
 */
$path = $cli->getVar('p', "./", 'path to data file');
$ui_expr = $cli->getVar('x', '', 'experiment name');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');

$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');

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

/*
 * establish the connection to the database
 */
$spec = new Spectra($ui_server);

if (preg_match("/iso[topic]*/", $method))
    $method = 'isotopic determination';

if (preg_match("/int[ensity]*/", $method))
    $method = 'peak intensity filtering';

if ($ui_expr == '')
    $ui_expr = uniqid("exp");

$files = list_files($path, $file_regex . '.*\.mgf$', TRUE);

$prg->cli->header("MGF-JSON START");
foreach ($files as $file) {

    if (! is_false($data = $spec->fileScansExist(basename($file)))) {
        $prg->cli->message("file exists", $data['file_id'] . " n:" . $data['count']);
        continue;
    }

    /*
     * read in the data
     */
    $mgf_id = uniqid("mz");
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
    $file_data = [
        'file' => basename($file),
        'file_id' => $mgf_id,
        'path' => $path,
        'scans_n' => $mgf_n,
        'experiment_name' => $ui_expr
    ];
    $peaks_method = [
        'file_id' => $mgf_id,
        'peaks_max_n' => $peaks_max_n,
        'peak_filtering' => $method,
        'peak_filter_da' => $peaks_filter_da
    ];
    if (preg_match("/iso[topic]*/", $method))
        unset($method['peaks_filter_da']);

    /*
     * construct the base name
     */
    $file = basename($file, '.mgf');

    $json_out = [];
    $prg->start($file . ':' . $mgf_n, $mgf_n);

    $n_break = 0;
    foreach ($data as $scan => $this_data) {

        $n_break ++;
        $prg->barPercent();

        $scan_id = $scan;
        $pre_m = 0;
        $pre_z = 0;
        $pre_int = 0;
        $title = 0;
        $rtsec = 0;

        if (key_exists('pepmass', $this_data)) {
            $pre_m = round(preg_replace("/\s+.*/", "", $this_data['pepmass']), 4);
            $pre_int = round(preg_replace("/.+\s/", "", $this_data['pepmass']) * 1, 4);
        }

        if (key_exists('charge', $this_data))
            $pre_z = preg_replace("/[\+\s+].*/", "", $this_data['charge']) * 1;

        /*
         * override the zeros
         */
        $pre_z = max(1, $pre_z) * 1;

        if (key_exists('title', $this_data))
            $title = preg_replace("/[\'\\\"\,]+/", "", $this_data['title']);

        if (key_exists('rtinseconds', $this_data))
            $rtsec = round($this_data['rtinseconds'], 2);

        $this_spectra = $this_data['spec'];
        if ($pre_int == 0)
            $pre_int = array_sum($this_spectra);

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
            $peaks = [];

        $json_out[$scan_id]['scan_id'] = $scan_id;
        $json_out[$scan_id]['title'] = $title;
        $json_out[$scan_id]['rt_sec'] = $rtsec;
        $json_out[$scan_id]['pre_mz'] = $pre_m;
        $json_out[$scan_id]['pre_z'] = $pre_z;
        $json_out[$scan_id]['pre_int'] = round(log10($pre_int), 2);

        foreach ($peaks as $var => $val) {
            $json_out[$scan_id][$var] = $val;
        }

        $json_out[$scan_id]['peaks'] = $json_out[$scan_id]['peaks'];
        $json_out[$scan_id]['spectrum'] = $this_spectra;

        // $n_mass = neutralMass($pre_m, max(2, $pre_z));
        // $n_bins = averagineBins($pre_m, max(2,$pre_z));
        // $a_bins = max(0.01, $peaks['a_bins']);
        // $a_mass = averagineMass();
        // /*
        // * expected number of peptides in the search space give a
        // * database of ~8000 non-redundant peptides
        // */
        // $n_search = max(250, 2800 + $n_mass * - 0.8);
        // /*
        // * log-normal distribution
        // * fit from poly^2 :: overlap ~ pre_nmass from large dataset
        // */
        // $p_mean = - 0.5577285 * pow(log($n_mass), 2) + 8.6588187 * log($n_mass) - 31.8397488;
        // $p_sd = 0.06632917 * pow(log($n_mass), 2) - 1.17749888 * log($n_mass) + 5.43862010;
        // $psm_lnorm = new LogNormal($p_mean, $p_sd);

        // $pmatch = $psm_lnorm->cdf($n_mass * $a_bins / $a_mass);
        // if (is_nan($pmatch))
        // $pmatch = 0;

        // $json_out[$scan_id]['spec_prob'] = sigfigs($pmatch, 5);

        // /*
        // * if there is a non-float p-value in match fail out
        // */
        // if ($json_out[$scan_id]['spec_prob'] . "" == 'NAN')
        // $prg->cli->flusheol("ERROR", "exception to spectral quality scan:" . sizeof($json_out), TRUE);
    }

    $tmr->start();
    /*
     * load file attributes
     */
    $file_pk = $spec->putFileData($file_data);

    if (is_false($file_pk))
        $cli->message("ERROR", 'failed to load file attributes to db', TRUE);
    /*
     * save the output
     * this is stored as long-data, i.e. as column-based arrays
     */

    $table = table_invert($json_out);

    /*
     * not all scans made it .. some have no peaks
     */

    $table = table_reindex($table);
    $table['peaks'] = array_map('json_encode', $table['peaks']);

    $spec->putScanData($table, $peaks_method, $file_pk);
    $prg->cli->flusheol('push to sql db', $tmr->timeinstr() . " | " . $mgf_id);
}

$prg->header("COMPLETE");

?>