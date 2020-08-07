#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
calculates the q-value [qvalue], [eprob] and [psm_prob] for each peptide-spectrum-match.

-- [psm_qvalue] the specific FDR associated with a give score value computed from a population of suspected positives (target) and known negatives (decoy)
-- [psm_eprob] the posterior error probability proposed by Kall, L. et. al J. Proteome Res. 7, 40–44 (2008).
-- [psm_prob] the peptide PSM probability calculated using the dscore value according to Keller, A., et. al. Anal. Chem. 74, 5383–5392 (2002).";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\ProgressTimer;
use MathPHP\Statistics\KernelDensityEstimation;
use ProxyIO\Cli;

// $mgf = new ReadMGF("/data/SCBI/20171201_DBS_20Samples/mgf/20171201_SoCal_Samples18-9.mgf");

$cli = new Cli();
$pgt = new ProgressTimer();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v0.2.190917');
$cli->Help->setUsage('-p <path> -f <regex>');
$cli->Help->setExample('-p /home/scbi/data -f som.*json -fdr 0.05');

/*
 * load the variables
 * flag : default-value : description
 */
$path = $cli->getVar('p', "./", 'path to data file');
$file_regex = $cli->getVar('f', '', 'file regex as a means to filter');
$valu = $cli->getVar('v', "dscore", 'variable to use for fdr analysis');
$jvln_fdr = $cli->getVar('fdr', 0.1, 'fdr cutoff for print analysis');

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
$pgt->message(" user fdr", $jvln_fdr);

$jvln_table = [];
$pgt->header("BEGIN");
foreach ($files as $file) {

    $data = json_decode(file_get_contents($file), TRUE);

    if (! key_exists('jaevalen', $data)) {
        $pgt->message("f:" . basename($file, '.json'), "no jaevalen data found");
        continue;
    }

    $data_jvln = $data['jaevalen']['psms'];

    // store the scores raw
    $stats_scores_t = [];
    $stats_scores_f = [];

    // bin the scan_ids by score
    $stats_ids_t = [];
    $stats_ids_f = [];

    $data_out = [];

    foreach ($data_jvln as $scan_id => $scan) {

        if (! key_exists($valu, $scan))
            $cli->message("ERROR", "$valu not found", TRUE);

        // store the scores raw
        $stats_scores_t[$scan_id] = $scan[$valu][0];
        if (count($scan[$valu]) > 1)
            $stats_scores_f[$scan_id] = $scan[$valu][1];

        // bin the scan_ids by score
        $stats_ids_t[floor($scan[$valu][0] * 1000)][] = $scan_id;
        if (count($scan[$valu]) > 1)
            $stats_ids_f[floor($scan[$valu][1] * 1000)][] = $scan_id;
    }

    $min = min(array_merge(array_keys($stats_ids_t), array_keys($stats_ids_f)));
    $max = max(array_merge(array_keys($stats_ids_t), array_keys($stats_ids_f)));

    $count_t = 0;
    $count_f = 0;
    // $score = 0;
    $target_fdr = 0;
    $target_score = 0;
    $peptide_true = 0;

    $kde_t = new KernelDensityEstimation($stats_scores_t);
    $kde_f = new KernelDensityEstimation($stats_scores_f);

    $pgt->start("f:" . basename($file, '.json'), ($max - $min));

    for ($s = $max; $s >= $min; $s --) {

        $pgt->barPercent();

        $scan_ids = [];

        if (key_exists($s, $stats_ids_t)) {
            $count_t += count($stats_ids_t[$s]);
            $scan_ids = $stats_ids_t[$s];
        }

        if (key_exists($s, $stats_ids_f))
            $count_f += count($stats_ids_f[$s]);

        if (count($scan_ids) == 0)
            continue;

        $score = $s / 1000;

        $pvalue_f = $kde_f->evaluate($score);
        $pvalue_t = $kde_t->evaluate($score);

        $pvalue = $pvalue_f / ($pvalue_f + $pvalue_t);
        $qvalue = $count_f / ($count_f + $count_t);

        if ($qvalue <= $jvln_fdr) {
            $peptide_true += count($scan_ids);
            $target_fdr = max($qvalue, $target_fdr);
            $target_score = max($score, $target_score);
        }

        foreach ($scan_ids as $scan_id) {
            $data_out[$scan_id]['eprob'] = number_format($pvalue, 5);
            $data_out[$scan_id]['qvalue'] = number_format($qvalue, 5);
        }
    }

    $n_pass = $peptide_true;
    $p_pass = number_format($n_pass / count($data_jvln) * 100, 1) . "%";

    $out_message = "pass: " . $n_pass;
    $out_message .= " (" . $p_pass . ")";
    $out_message .= "  score cutoff: " . str_pad($target_score, 6);
    $out_message .= "  fdr est.: " . number_format($target_fdr, 4);
    $pgt->message("f:" . basename($file, '.json'), $out_message);

    $data['performance'] = [
        'variable' => $valu,
        'cutoff' => $score,
        'estimated fdr' => $target_fdr,
        'n pass' => $n_pass,
        'percent pass' => $p_pass,
        'psms' => $data_out
    ];

//     file_put_contents($file, json_encode($data));
}

$pgt->header("COMPLETE");

?>