#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
searches the jaevalen engine for matching peptides based on the
precursor and fragments masses";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\ProgressTimer;
use ProxyTime\Timer;
use Jvln\Engines\JvlnQuery;
use Jvln\System\JsonConfig;
use ProxyIO\Cli;
use Jvln\Database\Spectra;
use Jvln\Database\PeptideMatches;
use ProxyIO\File\Delim\Delim;
use ProxyIO\File\Delim\WriteDelim;

$cli = new Cli();
$pgt = new ProgressTimer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.9.190917');
$cli->Help->setUsage('-f <file_id>');
$cli->Help->setExample('-db yeast -nph 5');

/*
 * load the variables
 * flag : default-value : description
 */
$ui_test = $cli->getVar('t', false, 'flag for testing, input the scan N to test');

$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');

$file_id = $cli->getVar('f', '', 'file id in database');
$jvln_sl = $cli->getVar('n', 0, 'scan limit n');
$jvln_index = $cli->getVar('db', 'uniprot', 'proteomic data base');
$jvln_sm = $cli->getVar('s', 0.1, 'minimum score to report');
$jvln_lh = $cli->getVar('l', 10, 'hits per spectrum return limit');
$jvln_quorum = $cli->getVar('nq', 3, 'quorum match number ');
$jvln_nmusthave = $cli->getVar('nm', 3, 'n top must have at least 1');
$jvln_npeaksmin = $cli->getVar('nl', 4, 'n peaks min');
$jvln_npeaksmax = $cli->getVar('nu', 36, 'n peaks max');
$jvln_range_iso = $cli->getVar('ir', "0:1", 'precursor isotopic range to test, use 0 for no range
-- ');

$jvln_range_chr = $cli->getVar('zr', "2:3", 'precursor charge state range to test, including as stated
use 0 to use charge state as is');

$ui_optimize = $cli->getVar('o', FALSE, 'optimize the search by sending to csv');
/*
 * display the help if needed
 */
$cli->Help->displayHelp();

/*
 * establish the connection to the database
 */
$spec = new Spectra($ui_server);
$psms = new PeptideMatches($ui_server);

$jvln_s = '127.0.0.1';
$jvln_db = NULL;
$jvln_un = 'root';
$jvln_pw = '';

$pgt->header("JVLN START");

/*
 * abort if the file has already been searched
 */
if (! is_false($data = $spec->filePsmsExist($file_id)) & is_false($ui_optimize))
    $pgt->cli->message("file has been searched", $data['file_id'] . " n:" . $data['count'], TRUE);

$jsc = new JsonConfig("/var/jvlnse/" . $jvln_index);
$jvln_port = $jsc->getVariable('index', 'port') + 1;

if (is_null($jvln_port))
    $cli->message("not a viable database", $jvln_index, TRUE);

$pgt->message(" db port", $jvln_port);

/*
 * set up the JVLN query
 */
$jvn = new JvlnQuery($jvln_index, $jvln_port);
$jvn->setQuorum($jvln_quorum);
$jvn->setNMustHavePeaks($jvln_nmusthave);
$jvn->setReturnLimit(max($jvln_lh, 100));
$db_stats = $jvn->getDBStats($jvln_s);

$pgt->header("CONFIG");
$pgt->message(" jvln database", $jvln_index);
$pgt->message(" jvln score-min", $jvln_sm);
$pgt->message(" pks  limit", $jvln_lh);
$pgt->message(" pks  quorum", $jvln_quorum);
$pgt->message(" pks  n top must have", $jvln_nmusthave);
$pgt->message(" pks  n minimum", $jvln_npeaksmin);
$pgt->message(" pks  n maximum", $jvln_npeaksmax);

$pgt->message(" val  range iso", $jvln_range_iso);
$pgt->message(" val  range charges", $jvln_range_chr);

/*
 * format the range for isotopic considerations
 */
if ($jvln_range_iso == 0) {
    $jvln_range_iso = [
        0
    ];
} elseif (! strstr($jvln_range_iso, ":")) {
    $jvln_range_iso = [
        0
    ];
} else {
    $jvln_range_iso = array_unique(explode(":", $jvln_range_iso));
}
/*
 * format the range for isotopic considerations
 */
if ($jvln_range_chr == 0) {
    $jvln_range_chr = FALSE;
} elseif (! strstr($jvln_range_chr, ":")) {
    $jvln_range_chr = [
        2
    ];
} else {
    $jvln_range_chr = array_unique(explode(":", $jvln_range_chr));
    $jvln_range_chr = range($jvln_range_chr[0], $jvln_range_chr[1]);
}

$pgt = new ProgressTimer();
$pgt->header("BEGIN");
/*
 * get the data to process
 */
if (! is_false($ui_test)) {
    $scans = $spec->getScan($file_id, $ui_test);
} else {
    $scans = $spec->getScans($file_id, $jvln_sl);
}
if ($scans == FALSE)
    $pgt->cli->message("ERROR", "..wtf..", TRUE);

$scans = table_invert($scans);

$jvln_params = [
    'server' => $jvln_s,
    'database' => $db_stats,
    'psm_score_min' => $jvln_sm,
    'psm_limit' => $jvln_lh,
    'peak_quorum' => $jvln_quorum,
    'peak_n_top_must_have' => $jvln_nmusthave,
    'peak_n_minimum' => $jvln_npeaksmin
];

$jvln_out = [];
$jvln_alt = [];

/*
 * SQL SERVER CONN for JVLN
 */

$jvln_c = @new mysqli($jvln_s, $jvln_un, $jvln_pw, $jvln_db, $jvln_port);

if ($jvln_c->connect_errno) {
    $pgt->message(" - connect error: " . $jvln_c->connect_errno, $jvln_c->connect_error);
    $pgt->header("EARLY TERMINATION", TRUE);
}
$jvln_c_id = $jvln_c->thread_id;

// $pgt->start("n:" . count($scans) . " f:" . $file_id, count($scans));
$pgt->cli->flush("start", "---");

/*
 * FOREACH SPECTRA -- START
 */

$tmr = new Timer();
$tmr->start();
foreach ($scans as $scan_n => $scan) {

    // $pgt->barPercent();
    // $pgt->barCount();

    $this_scan_pk = $scan['scan_pk'];
    $this_pre_mass = $scan['pre_mz'];
    $this_pre_z = $scan['pre_z'];
    $this_peaks = json_decode(trim($scan['peaks'], '"'));

    $jvln_data = [];
    /*
     * don't bother if we don't have enough peaks
     */
    if (count($this_peaks) >= $jvln_npeaksmin) {

        /*
         * slice down to the number of peaks requested
         */
        $this_peaks = array_slice($this_peaks, 0, $jvln_npeaksmax);

        /*
         * propigate the charge states to test
         */
        $this_pre_za = [
            $this_pre_z
        ];
        if ($jvln_range_chr != FALSE) {
            $this_pre_za = array_merge($jvln_range_chr, $this_pre_za);
            $this_pre_za = range(min($this_pre_za), max($this_pre_za));
        }
        /*
         * get a count of all peptides by precursor considered
         */
        // $jvln_pre_n = 0;
        // $jvln_pn = $jvn->getPrecursorCount($this_pre_mass, $this_pre_za, $jvln_range_iso);
        // if ($jvln_r = $jvln_c->query($jvln_pn)) {
        // $jvln_pre_n = $jvln_r->fetch_all(MYSQLI_ASSOC);
        // $jvln_pre_n = $jvln_pre_n[0]['n'];
        // $jvln_r->free();
        // }

        /*
         * get the list of matching peptides
         */

        if (is_false($ui_optimize)) {
            $jvln_q = $jvn->getQuery($this_peaks, $this_pre_mass, $this_pre_za, $jvln_range_iso);
        } else {
            $jvln_q = $jvn->getQueryModel($this_peaks, $this_pre_mass, $this_pre_za, $jvln_range_iso);
        }

        if ($jvln_r = $jvln_c->query($jvln_q)) {
            $jvln_data = $jvln_r->fetch_all(MYSQLI_ASSOC);
            $jvln_r->free();
        }

        // $jvln_time_sec = number_format($tmr->timeinmsec(), 4);

        if (count($jvln_data) > 0) {

//             print_r($jvln_data);
//             exit();

            $jvln_data = jvln_hits_munge_v2($jvln_data, $this_peaks, $this_pre_mass, $this_pre_z, $jvln_lh);

            $jvln_data['hits'] = table_fillCol($jvln_data['hits'], 'scan_pk', $this_scan_pk);
            $jvln_data['alts'] = table_fillCol($jvln_data['alts'], 'scan_pk', $this_scan_pk);

            $jvln_alt = table_bind($jvln_alt, $jvln_data['alts']);
            foreach (table_invert($jvln_data['hits']) as $row) {
                $jvln_out[] = $row;
            }
        }
    }
}
$pgt->cli->flusheol("finished", $tmr->timeinstr(TRUE));
$jvln_err = $jvln_c->error;
$jvln_c->close();

$jvln_out = table_invert($jvln_out);
$jvln_out = table_dropcols($jvln_out, 'alt_prob');
$jvln_alt = table_dropcols($jvln_alt, 'rank_prob');

if (! is_false($ui_test)) {
    print_table($jvln_out, min(50, $jvln_lh));

    if (table_length($jvln_alt) > 0)
        print_table($jvln_alt, min(50, $jvln_lh));

    exit();
}

if (! is_false($ui_optimize)) {

    $csv = new WriteDelim();
    $csv->writeCsvFile("/var/jvlnse/local/optimize/" . $file_id . ".csv", $jvln_out);

    $pgt->header("COMPLETE");
    exit();
}

$psms->putPsmData($jvln_out, $jvln_params);

if (table_length($jvln_alt) > 0)
    $psms->putPsmAltData($jvln_alt);

$pgt->header("COMPLETE");

?>