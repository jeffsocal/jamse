#!/usr/bin/env php
<?php
$Copyright = "
Written by Jeff Jones (jeff@socalbioinformatics.com)
Copyright (2016) SoCal Bioinformatics Inc.
See LICENSE.txt for the license.";

$Description = "
identifies and associates all proteins from which a given peptide
originated";

ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

use ProxyTime\ProgressTimer;
use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\Protein;
use Jvln\System\JsonConfig;
use ProxyIO\Cli;
use ProxyIO\File\Delim\WriteDelim;
use Jvln\Database\Performance;
use ProxyMySQL\Transact;

// $mgf = new ReadMGF("/data/SCBI/20171201_DBS_20Samples/mgf/20171201_SoCal_Samples18-9.mgf");

$cli = new Cli();
$tmr = new Timer();
$pgt = new ProgressTimer();
$pro = new Protein();
$csv = new WriteDelim();

/*
 * setup the help display
 */
$cli->Help->setCopyright($Copyright);
$cli->Help->setDescription($Description);
$cli->Help->setVersion('BETA v1.1.190521');
$cli->Help->setUsage('-f <file_id>');
$cli->Help->setExample('-f <file_id> -db yeast -m psm_score');

/*
 * load the variables
 * flag : default-value : description
 */
$ui_file = $cli->getVar('f', '', 'file id in database');
$ui_expr = $cli->getVar('x', '', 'experiment name');

$ui_expx = $cli->getVar('xi', 'indv', 'instructions for experiment [show|individual|group]');
$ui_report = $cli->getVar('r', 'protein', 'report type [protein|peptide]');

if (preg_match("/pro[tein]*/", $ui_report)) {
    $ui_report = 'protein';
} else {
    $ui_report = 'peptide';
}

$jvln_fdr = $cli->getVar('fdr', 0.05, 'target FDR');
$jvln_m = $cli->getVar('m', 'psm_score', 'target FDR metric [psm_fcorr, psm_score, psm_prob]');
$jvln_index = $cli->getVar('db', 'uniprot', 'proteomic data base');
$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

$pgt->header("JVLN START");
$psms = new Performance($ui_server);

$jsc = new JsonConfig("/var/jvlnse/" . $jvln_index);

if (is_false($jsc->getConfigFile()))
    $cli->message("user input error", "choose a valid database", TRUE);

$jvln_port = $jsc->getVariable('index', 'port') + 1;

$jvln_s = '127.0.0.1';
$jvln_db = NULL;
$jvln_un = 'root';
$jvln_pw = '';

if ($ui_expx == 'show') {
    $table = $psms->listExperiments();
    echo PHP_EOL;
    print_table($table, table_length($table));
    exit();
}

if ($ui_file == '' and $ui_expr == '')
    $cli->message("user input error", "choose a file or experiment to run", TRUE);

if (! preg_match("/psm_fcorr|psm_score|psm_prob|hit_prob/", $jvln_m))
    $cli->message("user input error", "metric not a valid choice [psm_fcorr, psm_score, psm_prob]", TRUE);

if (preg_match("/psm_prob|hit_prob/", $jvln_m))
    $jvln_fdr = 1 - $jvln_fdr;

$jvln_port = $jsc->getVariable('index', 'port') + 1;
if (is_null($jvln_port))
    $cli->message("not a viable database", $jvln_index, TRUE);

$pgt->header("CONFIG");
$pgt->message(" file id", $ui_file);
$pgt->message(" target FDR", $jvln_fdr);
$pgt->message(" metric FDR", $jvln_m);
$pgt->message(" jvln database", $jvln_index);
$pgt->message(" db port", $jvln_port);

function topProtein($v_pro, $v_pep, $v_pep_org)
{
    $v_pro_stat = [];
    foreach ($v_pro as $protein => $peptides) {
        $stats = proteinProbabilities($peptides, $v_pep, $v_pep_org);
        $v_pro_stat[$protein] = $stats['prot_prob_score'];
    }

    return array_keys(array_intersect($v_pro_stat, [
        array_max($v_pro_stat)
    ]));
}

function proteinProbabilities($peptides, $v_pep, $v_pep_org)
{
    $v_pep_t = 0; // total observational frequency
    $v_pep_u = 0; // unique observational frequency
    $v_pep_a = []; // probabilities of all associated peptides
    $v_pep_i = []; // probabilities of all identity peptides
    $v_pep_r = []; // probabilities of round robin remaining peptides

    $v_pep = array_intersect_key($v_pep, array_flip($peptides));
    $v_pep_org = array_intersect_key($v_pep_org, array_flip($peptides));

    foreach ($v_pep_org as $vars) {
        $v_pep_t += $vars['obs_total'];
        $v_pep_u += $vars['obs_unique'];
        $v_pep_a[] = $vars['max_prob'];
        if (count($vars['source_protein_id']) == 1)
            $v_pep_i[] = $vars['max_prob'];
    }

    foreach ($v_pep as $vars) {
        $v_pep_r[] = $vars['max_prob'];
    }

    $out['n_pep_t'] = $v_pep_t; // total observational frequency
    $out['n_pep_u'] = $v_pep_u; // unique observational frequency
    $out['n_pep_a'] = count($peptides); // unique sequence observational frequency
    $out['n_pep_i'] = count($v_pep_i); // identity peptides
    $out['n_pep_r'] = count($v_pep_r); // round robin remaining peptides

    $out['prot_prob_a'] = array_probability($v_pep_a);
    $out['prot_prob_i'] = array_probability($v_pep_i);
    $out['prot_prob_r'] = array_probability($v_pep_r);

    $score = $out['prot_prob_a'] * $out['n_pep_a'];
    $score += $out['prot_prob_i'] * $out['n_pep_i'];
    $score += $out['prot_prob_r'] * $out['n_pep_r'];

    $score *= $out['n_pep_t'];

    $out['pep_r'] = $out['n_pep_i'] / $out['n_pep_a']; // identity / unique
    $out['prot_prob_score'] = log($score);

    return $out;
}

$pgt->header("BEGIN");

$tmr->start();

if ($ui_file != '') {
    $ui_files = [
        $ui_file
    ];
} elseif ($ui_expr != '') {
    $exp_files = $psms->getFilesByExperiment($ui_expr);
    if (is_false($exp_files))
        $pgt->cli->message("ERROR:", "no files found for that group", TRUE);

    $ui_files = $exp_files['file_id'];
}

$data_tables = [];
foreach ($ui_files as $ui_file) {
    $this_table = $psms->getPsmPerf($ui_file, $jvln_m, $jvln_fdr);
    $pgt->cli->message($ui_file . " n peptides", table_length($this_table));
    if (is_false($this_table))
        $pgt->cli->message("ERROR:", "no data", TRUE);

    $data_tables[] = $this_table;
}

if ($ui_expr != '' & $ui_expx == 'group') {
    $new_table = [];
    foreach ($data_tables as $data_table) {
        $new_table = table_bind($new_table, $data_table);
    }
    $data_tables = [];
    $data_tables[] = $new_table;
}

foreach ($data_tables as $data_table) {
    $table_prot = []; // master table of protein meta data
    $v_pep_prob = []; // bare peptides by max probabilities
    $v_pep_obst = []; // number of observed peptides
    $v_pep_obsu = []; // number of observed peptides
    $v_pep = []; // array of peptides by key proteins by value
    $v_pro = []; // array of proteins by key peptides by value
    foreach (table_invert($data_table) as $vals) {
        $peptide = $peptide_ptm = $vals['peptide'];

        if (strstr("psm_prob|hit_prob", $jvln_m)) {
            $prob = $vals[$jvln_m];
        } else {
            $prob = $vals['fdr_prob'];
        }

        $peptide = getRootSequence($peptide);
        if (! key_exists($peptide, $v_pep_prob)) {
            $v_pep_prob[$peptide] = $prob;
            $v_pep_obst[$peptide] = 0;
            $v_pep_obsu[$peptide] = [];
        }

        $v_pep_obst[$peptide] ++;
        $v_pep_obsu[$peptide][$peptide_ptm] = 1;

        $v_pep_prob[$peptide] = max($v_pep_prob[$peptide], $prob);
    }

    /*
     * master array of peptide stats
     */
    $v_pep = table_invert([
        'max_prob' => $v_pep_prob,
        'obs_total' => $v_pep_obst,
        'obs_unique' => array_map('array_sum', $v_pep_obsu)
    ]);

    /*
     * cleanup
     */
    unset($v_pep_prob);
    unset($v_pep_obst);
    unset($v_pep_obsu);

    /*
     * SQL SERVER CONN
     */
    $jvln_c = @new mysqli($jvln_s, $jvln_un, $jvln_pw, $jvln_db, $jvln_port);
    if ($jvln_c->connect_errno) {
        $pgt->message(" - connect error: " . $jvln_c->connect_errno, $jvln_c->connect_error);
        $pgt->header("EARLY TERMINATION", TRUE);
    }

    /*
     * FOREACH SPECTRA -- START
     */
    $pgt->start(" total peptides:" . count($v_pep), count($v_pep));

    foreach ($v_pep as $peptide => $stats) {

        $pgt->barPercent();

        $jvln_q = "SELECT 
                    id source_protein_id,
                    protein_source, protein_name, protein_id, 
                    protein_desc, organism, sequence 
                    FROM peptides 
                    WHERE 
                    MATCH('" . strtolower($peptide) . "') limit 100;";

        $jvln_data = FALSE;
        if ($jvln_r = $jvln_c->query($jvln_q)) {
            $jvln_data = $jvln_r->fetch_all(MYSQLI_ASSOC);
            $jvln_r->free();
        }

        $peptide = strtoupper($peptide);

        /*
         * munge a protid onto the db index for uniqueness
         */
        $jvln_data = table_invert($jvln_data);
        $jvln_data['source_protein_id'] = array_map(function ($x) {
            return "protid_" . $x;
        }, $jvln_data['source_protein_id']);

        /*
         * map the peptide to proteins
         */
        $v_pep[$peptide]['source_protein_id'] = array_unique($jvln_data['source_protein_id']);

        /*
         * return same dataframe
         * -- silly sphinx does not allow CONCAT in query
         */
        $table_prot = table_bind($table_prot, $jvln_data);
    }
    /*
     * close out the sql connection
     */
    $jvln_err = $jvln_c->error;
    $jvln_c->close();

    /*
     * construct the protein array
     */
    foreach ($v_pep as $peptide => $stats) {
        foreach ($stats['source_protein_id'] as $protein) {
            if (! key_exists($protein, $v_pro))
                $v_pro[$protein] = [];

            $v_pro[$protein][] = $peptide;
        }
    }

    /*
     * create the protein lookup array
     */
    $col_prot = $table_prot['source_protein_id'];
    $table_prot = array_combine($col_prot, table_invert($table_prot));

    /*
     * table_prot // master table of protein meta data
     * v_pep_prob // bare peptides by max probabilities
     * v_pep_obsv // number of observed peptides
     * v_pep // array of peptides by key proteins by value
     * v_pro // array of proteins by key peptides by value
     */

    if (preg_match("/pro[tein]*/", $ui_report)) {

        $v_pep_full = $v_pep;

        $total_peptides = count($v_pep_full);
        $timer = new Timer();
        $prot_accounting_final = [];
        $last_n_pep = 0;
        while (count($v_pep) > 0 & count($v_pro) > 0) {

            $last_n_pep = count($v_pep);

            $stdout = 'prot: ' . count($v_pro);
            $stdout .= '  pept: ' . count($v_pep);
            $stdout .= ' of ' . $total_peptides;
            $stdout .= ' -- time: ' . $timer->timeinstr(true);

            $cli->flush(" building proteins", $stdout);

            $this_peptides_final = [];
            $this_protein_final = '';

            $top_hit_stored = FALSE;
            $this_proteins = topProtein($v_pro, $v_pep, $v_pep_full);
            foreach ($this_proteins as $this_protein) {
                $this_peptides = $v_pro[$this_protein];
                $this_probabilities = proteinProbabilities($this_peptides, $v_pep, $v_pep_full);

                if ($this_probabilities['pep_r'] == 0) {
                    $v_pro = array_diff_key($v_pro, array_flip([
                        $this_protein
                    ]));
                } elseif ($top_hit_stored == FALSE) {
                    $this_peptides_final = $this_peptides;
                    $this_protein_final = $this_protein;

                    $this_protein_acc = $table_prot[$this_protein];
                    $this_protein_acc['sequence_coverage'] = $pro->getCoverage($this_protein_acc['sequence'], $this_peptides)['coverage'];
                    $this_probabilities['prot_prob_score'] = round($this_probabilities['prot_prob_score'], 2);
                    $prot_accounting_final[] = array_merge($this_protein_acc, $this_probabilities);
                    $top_hit_stored = TRUE;
                }
            }

            /*
             * remove peptides from the pool for the next iteration
             */
            $v_pep = array_diff_key($v_pep, array_flip($this_peptides_final));
            $v_pro = array_diff_key($v_pro, array_flip([
                $this_protein_final
            ]));

            if ($last_n_pep == count($v_pep)) {
                /*
                 * remove all proteins with no residual peptides
                 */
                foreach (array_keys($v_pro) as $this_protein) {
                    $this_peptides = $v_pro[$this_protein];
                    $this_probabilities = proteinProbabilities($this_peptides, $v_pep, $v_pep_full);

                    // print_r($v_pro[$this_protein]);
                    // print_r($this_probabilities);
                    // exit;

                    if ($this_probabilities['n_pep_r'] == 0) {
                        $v_pro = array_diff_key($v_pro, array_flip([
                            $this_protein
                        ]));
                    }
                }
            }
        }

        /*
         * one last print message
         */
        $stdout = 'prot: ' . count($v_pro);
        $stdout .= '  pept: ' . count($v_pep);
        $stdout .= ' of ' . $total_peptides;
        $stdout .= ' -- time: ' . $timer->timeinstr(true);

        $cli->flusheol(" complete", $stdout);
        $prot_accounting_final = table_invert($prot_accounting_final);
    }

    if (preg_match("/pep[tide]*/", $ui_report)) {

        $prot_accounting_final = [];
        foreach (table_invert($data_table) as $peptide_vals) {
            $peptide = getRootSequence($peptide_vals['peptide']);

            foreach ($v_pep[$peptide]['source_protein_id'] as $protein_key) {
                $protein_vals = $table_prot[$protein_key];
                unset($protein_vals['source_protein_id']);
                unset($protein_vals['sequence']);
            }

            $prot_accounting_final[] = array_merge($peptide_vals, $protein_vals);
        }
        $prot_accounting_final = table_invert($prot_accounting_final);
    }

    $ui_file_name = $ui_file;
    if ($ui_expr != '' & $ui_expx == 'group')
        $ui_file_name = strtolower(preg_replace("/\s+/", "_", trim($ui_expr)));

    $file_name = "/var/jvlnse/local/" . $ui_report . "s/" . $ui_file_name . "_" . $jvln_m . "-" . $jvln_fdr . "_" . $ui_report . ".csv";
    $csv->writeCsvFile($file_name, $prot_accounting_final);
}
$pgt->header("COMPLETE");

?>