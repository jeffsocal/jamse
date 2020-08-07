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

ini_set("memory_limit", "8192M");
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

// ###################################################################################################
// XML
// id: MS:1000133 name: collision-induced dissociation .................................... synonym: "CID"
// id: MS:1000134 name: plasma desorption ................................................. synonym: "PD"
// id: MS:1000135 name: post-source decay ................................................. synonym: "PSD"
// id: MS:1000136 name: surface-induced dissociation ...................................... synonym: "SID"
// id: MS:1000242 name: blackbody infrared radiative dissociation ......................... synonym: "BIRD"
// id: MS:1000250 name: electron capture dissociation ..................................... synonym: "ECD"
// id: MS:1000262 name: infrared multiphoton dissociation ................................. synonym: "IRMPD"
// id: MS:1000282 name: sustained off-resonance irradiation ............................... synonym: "SORI"
// id: MS:1000433 name: low-energy collision-induced dissociation.......................... synonym: "LECID"
// id: MS:1000435 name: photodissociation ................................................. synonym: "MPD"
// id: MS:1000598 name: electron transfer dissociation .................................... synonym: "ETD"
// id: MS:1000599 name: pulsed q dissociation ............................................. synonym: "PQD"
// id: MS:1001880 name: in-source collision-induced dissociation .......................... synonym: "ISCID"
// id: MS:1002000 name: LIFT .............................................................. synonym: "LIFT"
// id: MS:1002631 name: Electron-Transfer/Higher-Energy Collision Dissociation (EThcD) .... synonym: "EThcD"
// id: MS:1000422 name: beam-type collision-induced dissociation .......................... synonym: "HCD"
// id: MS:1002472 name: trap-type collision-induced dissociation .......................... synonym: "HCD"
// id: MS:1002679 name: supplemental collision-induced dissociation ....................... synonym: "HCD"

$cid['MS:1000133'] = "CID";
$cid['MS:1000134'] = "PD";
$cid['MS:1000135'] = "PSD";
$cid['MS:1000136'] = "SID";
$cid['MS:1000242'] = "BIRD";
$cid['MS:1000250'] = "ECD";
$cid['MS:1000262'] = "IRMPD";
$cid['MS:1000282'] = "SORI";
$cid['MS:1000433'] = "LECID";
$cid['MS:1000435'] = "MPD";
$cid['MS:1000598'] = "ETD";
$cid['MS:1000599'] = "PQD";
$cid['MS:1001880'] = "ISCID";
$cid['MS:1002000'] = "LIFT";
$cid['MS:1002631'] = "EThcD";
$cid['MS:1000422'] = "HCD";
$cid['MS:1002472'] = "HCD";
$cid['MS:1002679'] = "HCD";

$time['minute'] = 60;
$time['second'] = 1;
$time['millisecond'] = 1 / 1000;
// ###################################################################################################

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

$files = list_files($path, $file_regex . '.*\.mz(ML|ml)$', TRUE);

$prg->cli->header("PARSE MZML START");
foreach ($files as $file) {

    if (! is_false($data = $spec->fileScansExist(basename($file)))) {
        $prg->cli->message("file exists", $data['file_id'] . " n:" . $data['count']);
        continue;
    }

    /*
     * read in the data
     */
    $mgf_id = uniqid("mz");
    $mzml = new XMLReader();

    if (! is_null($scan_limit)) {
        $data = array_slice($data, 0, $scan_limit);
        $mgf_n = count($data);
        $prg->cli->message(' limiting to', $mgf_n);
    }

    $scans = [];
    $mzml->open($file);

    while ($mzml->read() && $mzml->name != 'spectrumList');

    $n_spec = $mzml->getAttribute('count') * 1;

    /*
     * ======================================== file info
     */
    $peaks_method = [
        'file_id' => $mgf_id,
        'peaks_max_n' => $peaks_max_n,
        'peak_filtering' => $method,
        'peak_filter_da' => $peaks_filter_da
    ];
    /*
     * form the methods output
     */
    $file_data = [
        'file' => basename($file),
        'file_id' => $mgf_id,
        'path' => $path,
        'scans_n' => $n_spec,
        'experiment_name' => $ui_expr
    ];
    /*
     * load file attributes
     */
    $file_pk = $spec->putFileData($file_data);

    if (is_false($file_pk))
        $cli->message("ERROR", 'failed to load file attributes to db', TRUE);

    $n = 0;
    $prg->start("reading mzml", $n_spec * 2);
    while ($mzml->read()) { // start reading.

        if ($mzml->name != 'spectrum')
            continue;

        $scan = [];

        $scan['scan_id'] = $mzml->getAttribute('index') + 1;

        $n ++;

        $prg->barPercent();

        try {
            $xml = simplexml_load_string($mzml->readOuterXML());
            $is_ms2 = FALSE;
        } catch (Exception $e) {
            continue;
        }

        /*
         * get: TITLE, MS LEVEL
         */
        foreach ($xml->cvParam as $param) {
            $accession = $param->attributes()->accession->__toString();
            $name = $param->attributes()->name->__toString();
            $value = $param->attributes()->value->__toString();

            if ($name == "spectrum title" || $accession == 'MS:1000796')
                $scan['title'] = substr(preg_replace("/\"|\'/", "", $value), 0, 254);

            if ($name == "ms level" and $value == 2)
                $is_ms2 = TRUE;
        }

        if (is_false($is_ms2))
            continue;

        /*
         * get: RET TIME SEC
         */
        foreach ($xml->scanList->scan->cvParam as $param) {
            $accession = $param->attributes()->accession->__toString();
            $name = $param->attributes()->name->__toString();
            $value = $param->attributes()->value->__toString();
            $unitName = $param->attributes()->unitName;

            if (! is_null($unitName))
                $unitName = $unitName->__toString();

            if ($name == "scan start time" || $accession == 'MS:1000016')
                $scan['rt_sec'] = round($value * $time[$unitName], 1);

            if ($name == "ion injection time" || $accession == 'MS:1000927')
                $scan['ms_it_sec'] = sigfigs($value * $time[$unitName], 3);
        }

        $prec = $xml->precursorList->precursor;
        foreach ($prec->activation->cvParam as $param) {
            $accession = $param->attributes()->accession->__toString();
            $value = $param->attributes()->value->__toString();

            if (key_exists($accession, $cid))
                $scan['ms_ct'] = $cid[$accession];

            if ($accession == "MS:1000045")
                $scan['ms_ce'] = $value;
        }

        /*
         * get: PRE MZ, Z, INT
         */
        $scan['pre_mz'] = 1;
        $scan['pre_z'] = 1;
        $scan['pre_int'] = 1;
        foreach ($prec->selectedIonList->selectedIon->cvParam as $param) {
            $name = $param->attributes()->name->__toString();
            $value = $param->attributes()->value->__toString();

            if ($name == "selected ion m/z")
                $scan['pre_mz'] = round($value, 4);

            if ($name == "charge state")
                $scan['pre_z'] = $value;

            if ($name == "peak intensity")
                $scan['pre_int'] = $value;
        }

        /*
         * get: MS2 SPECTRUM
         */
        $got_spectrum = TRUE;
        try {
            $this_spectra = [];
            foreach ($xml->binaryDataArrayList->binaryDataArray as $bina) {
                $is_64bit = FALSE;
                $is_zlibc = FALSE;
                $fl_array = 'mz';
                foreach ($bina->cvParam as $param) {

                    $attr = $param->attributes()->name->__toString();

                    if (preg_match("/64[-bit]*/i", $attr))
                        $is_64bit = TRUE;

                    if (preg_match("/zlib/i", $attr))
                        $is_zlibc = TRUE;

                    if (preg_match("/intensity/i", $attr))
                        $fl_array = 'int';
                }

                $spec_vals = $bina->binary->__toString();
                if ($is_64bit == TRUE)
                    try {
                        $spec_vals = base64_decode($spec_vals);
                    } catch (Exception $e) {
                        $got_spectrum = FALSE;
                        continue;
                    }
                if ($is_zlibc == TRUE)
                    try {
                        $spec_vals = zlib_decode($spec_vals);
                    } catch (Exception $e) {
                        $got_spectrum = FALSE;
                        continue;
                    }

                $this_spectra[$fl_array] = unpack('d*', $spec_vals);
            }
        } catch (Exception $e) {
            /*
             * fail silently
             */
            continue;
        }
        /*
         * error in spectrum decoding
         * fail silently
         */
        if ($got_spectrum == FALSE)
            continue;

        if ($scan['pre_int'] == 1)
            $scan['pre_int'] = array_sum($this_spectra);

        $pre_m = $scan['pre_mz'];
        $pre_z = $scan['pre_z'];
        /*
         * pull out the peaks
         */
        if (preg_match("/iso[topic]*/", $method))
            $peaks = spec_get_peaks_by_iso($this_spectra, $peaks_max_n);

        if (preg_match("/int[ensity]*/", $method))
            $peaks = spec_get_peaks_by_int2($this_spectra, $pre_m, $pre_z, $peaks_max_n, $peaks_filter_da);

        if (is_false($peaks))
            $peaks = [];

        foreach ($peaks as $var => $val) {
            $scan[$var] = $val;
        }

        $scan['pre_int'] = round(log10($scan['pre_int'] + 1), 2);

        /*
         * truncate and round the spectral values
         */
        $this_spectra = table_sort($this_spectra, 'int', 'desc', 'number');
        $this_spectra = table_head($this_spectra, 360);
        $this_spectra = table_sort($this_spectra, 'mz', 'ascn', 'number');

        $this_spectra['mz'] = array_map(function ($x) {
            return round($x, 3);
        }, $this_spectra['mz']);

        $this_spectra['int'] = array_map(function ($x) {
            return round($x, 2);
        }, $this_spectra['int']);

        $scan['peaks'] = $scan['peaks'];
        $scan['spectrum'] = $this_spectra;

        $scans[] = $scan;

        if (sizeof($scans) >= 1000) {
            /*
             * save the output
             * this is stored as long-data, i.e. as column-based arrays
             */
            $table = table_invert($scans);

            /*
             * not all scans made it .. some have no peaks
             */
            $table = table_reindex($table);
            $table['peaks'] = array_map('json_encode', $table['peaks']);

            $spec->putScanData($table, $peaks_method, $file_pk);

            $scans = [];
            unset($table);
        }
    }
    $mzml->close();

    // ###################################################################################################

    if (sizeof($scans) > 0) {
        /*
         * save the output
         * this is stored as long-data, i.e. as column-based arrays
         */
        $table = table_invert($scans);

        /*
         * not all scans made it .. some have no peaks
         */
        $table = table_reindex($table);
        $table['peaks'] = array_map('json_encode', $table['peaks']);

        $spec->putScanData($table, $peaks_method, $file_pk);
    }

    $prg->cli->flusheol(' done', $mgf_id);
}

$prg->header("COMPLETE");

?>