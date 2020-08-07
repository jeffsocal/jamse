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

$cli = new Cli();
$pgt = new ProgressTimer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.9.190917');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('-p /home/scbi/data -f som.*json -db yeast');
$cli->Help->setExample('-db yeast -nph 5');

/*
 * load the variables
 * flag : default-value : description
 */
$test = $cli->getVar('t', false, 'flag for testing, input the scan N to test');
$path = $cli->getVar('p', "./", 'path to data file');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');
$jvln_sl = $cli->getVar('n', NULL, 'scan limit n');
$jvln_index = $cli->getVar('db', 'uniprot', 'proteomic data base');
$jvln_sm = $cli->getVar('s', 0.1, 'minimum score to report');
$jvln_lh = $cli->getVar('l', 2, 'hits per spectrum return limit');
$jvln_quorum = $cli->getVar('nq', 1, 'quorum match number ');
$jvln_nmusthave = $cli->getVar('nm', 0, 'n top must have at least 1');
$jvln_npeaksmin = $cli->getVar('nl', 4, 'n peaks min');
$jvln_npeaksmax = $cli->getVar('nu', 36, 'n peaks max');
$jvln_range_iso = $cli->getVar('ir', "0:1", 'precursor isotopic range to test, use 0 for no range
-- ');

$jvln_range_chr = $cli->getVar('zr', "2:3", 'precursor charge state range to test, including as stated
use 0 to use charge state as is');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

$jvln_s = '127.0.0.1';
$jvln_db = NULL;
$jvln_un = 'root';
$jvln_pw = '';

$pgt->header("JVLN START");
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
$jvn->setReturnLimit($jvln_lh);
$db_stats = $jvn->getDBStats($jvln_s);

$files = list_files($path, $file_regex . '.*\.json$', TRUE);

$pgt->header("CONFIG");
$pgt->message(" file path", $path);
$pgt->message(" file regex", $file_regex);
$pgt->message(" file n", count($files));
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


$pgt->header("BEGIN");
foreach ($files as $file) {

    $jvln_out = [];

    $jvln_out = [
        'parameters' => [
            'server' => $jvln_s,
            'database' => $db_stats,
            'psm_score_min' => $jvln_sm,
            'psm_limit' => $jvln_lh,
            'peak_quorum' => $jvln_quorum,
            'peak_n_top_must_have' => $jvln_nmusthave,
            'peak_n_minimum' => $jvln_npeaksmin
        ]
    ];

    $pgt = new ProgressTimer();

    /*
     * read the file
     */
    $data = json_decode(file_get_contents($file), TRUE);

    /*
     * remove any downstream data, i.e we are going to overwrite
     */
    $data_to_remove = array_diff(array_keys($data), [
        'file',
        'scans'
    ]);
    if (count($data_to_remove) > 0)
        foreach ($data_to_remove as $rm) {
            unset($data[$rm]);
        }

    /*
     * get the scans
     */
    $scans = $data['scans'];

    if (! is_false($test)) {
        $scans = [];
        $scans[] = $data['scans'][$test];
    }
    /*
     * only run a portion of the scans
     */
    if (! is_null($jvln_sl))
        $scans = array_slice($scans, 0, $jvln_sl);

    /*
     * SQL SERVER CONN
     */

    $jvln_c = @new mysqli($jvln_s, $jvln_un, $jvln_pw, $jvln_db, $jvln_port);

    if ($jvln_c->connect_errno) {
        $pgt->message(" - connect error: " . $jvln_c->connect_errno, $jvln_c->connect_error);
        $pgt->header("EARLY TERMINATION", TRUE);
    }
    $jvln_c_id = $jvln_c->thread_id;

    $pgt->start("n:" . count($scans) . " f:" . basename($file, '.json'), count($scans));

    /*
     * FOREACH SPECTRA -- START
     */

    $tmr = new Timer();
    foreach ($scans as $scan_n => $scan) {

        $pgt->barPercent();
        // $pgt->barCount();

        $this_pre_mass = $scan['pre_mz'];
        $this_pre_z = $scan['pre_z'];
        $this_peaks = $scan['peaks'];

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
            $jvln_pre_n = 0;
            $jvln_pn = $jvn->getPrecursorCount($this_pre_mass, $this_pre_za, $jvln_range_iso);
            if ($jvln_r = $jvln_c->query($jvln_pn)) {
                $jvln_pre_n = $jvln_r->fetch_all(MYSQLI_ASSOC);
                $jvln_pre_n = $jvln_pre_n[0]['n'];
                $jvln_r->free();
            }

            /*
             * get the list of matching peptides
             */
            $tmr->start();
            $jvln_q = $jvn->getQuery($this_peaks, $this_pre_mass, $this_pre_za, $jvln_range_iso);

            if ($jvln_r = $jvln_c->query($jvln_q)) {
                $jvln_data = $jvln_r->fetch_all(MYSQLI_ASSOC);
                $jvln_r->free();
            }

            $jvln_out['stats'][$scan_n]['time_sec'] = number_format($tmr->timeinmsec(), 4);
            $jvln_out['stats'][$scan_n]['psm_n'] = $jvln_pre_n;
            if (count($jvln_data) > 0) {

                $jvln_data = jvln_hits_munge($jvln_data, $this_peaks, $this_pre_mass, $this_pre_z, $jvln_pre_n);
                $jvln_out['psms'][$scan_n] = $jvln_data;
                if (! is_false($test))
                    $jvln_out['query'][$scan_n] = $jvln_q;
            }
        }
    }

    $jvln_err = $jvln_c->error;
    $jvln_c->close();

    if (! is_false($test)) {
        print_table($jvln_out['psms'][0], $jvln_lh);

        echo $jvln_out['query'][0] . PHP_EOL;

        $pgt->header("TEST COMPLETE");
        exit();
    }

    $data['jaevalen'] = $jvln_out;

    file_put_contents($file, json_encode($data));

    /*
     * summary
     */
    // if (! is_false($stats_table))
    // print_table(table_summary($stats_table));
}

$pgt->header("COMPLETE");

?>