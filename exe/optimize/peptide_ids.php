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
use Jvln\System\JsonConfig;
use ProxyIO\Cli;
use ProxyMySQL\Transact;
use Jvln\Database\PeptideMatches;
use ProxyIO\File\Delim\Delim;

// $mgf = new ReadMGF("/data/SCBI/20171201_DBS_20Samples/mgf/20171201_SoCal_Samples18-9.mgf");

$cli = new Cli();
$tmr = new Timer();
$pgt = new ProgressTimer();
$csv = new Delim();

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
$jvln_index = $cli->getVar('db', 'uniprot', 'proteomic data base');
$ui_server = $cli->getVar('ip', '127.0.0.1', 'server ip location');

/*
 * display the help if needed
 */
$cli->Help->displayHelp();

$pgt->header("JVLN START");
$psms = new PeptideMatches($ui_server);

$jsc = new JsonConfig("/var/jvlnse/" . $jvln_index);

if (is_false($jsc->getConfigFile()))
    $cli->message("user input error", "choose a valid database", TRUE);

$jvln_port = $jsc->getVariable('index', 'port') + 1;

$jvln_s = '127.0.0.1';
$jvln_db = NULL;
$jvln_un = 'root';
$jvln_pw = '';

if ($ui_expx == 'show') {
    if ($ui_expr != '') {
        $table = $psms->getFilesByExperiment($ui_expr);
        if (is_false($table))
            $pgt->cli->message("ERROR:", "no files found for that group", TRUE);
    } else {
        $table = $psms->listExperiments();
    }
    echo PHP_EOL;
    print_table($table, table_length($table));
    exit();
}

if ($ui_file == '' and $ui_expr == '')
    $cli->message("user input error", "choose a file or experiment to run", TRUE);

$jvln_port = $jsc->getVariable('index', 'port') + 1;
if (is_null($jvln_port))
    $cli->message("not a viable database", $jvln_index, TRUE);

$pgt->header("CONFIG");
$pgt->message(" file", $ui_file);
$pgt->message(" db-index", $jvln_index);
$pgt->message(" db-port", $jvln_port);

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
    $file_name = "/var/jvlnse/local/optimize/" . $ui_file . ".csv";
    $this_table = $csv->read($file_name);

    $pgt->cli->message($ui_file . " n peptides", table_length($this_table));
    if (is_false($this_table))
        $pgt->cli->message("ERROR:", "no data", TRUE);

    $data_tables[$ui_file] = $this_table;
}

if ($ui_expr != '' & $ui_expx == 'group') {
    $new_table = [];
    foreach ($data_tables as $data_table) {
        $new_table = table_bind($new_table, $data_table);
    }
    $data_tables = [];
    $data_tables[$ui_expr] = $new_table;
}

foreach ($data_tables as $ui_file => $data_table) {
    $table_prot = []; // master table of protein meta data
    $v_pep_prob = []; // bare peptides by max probabilities
    $v_pep_obst = []; // number of observed peptides
    $v_pep_obsu = []; // number of observed peptides
    $v_pep = []; // array of peptides by key proteins by value
    $v_pro = []; // array of proteins by key peptides by value
    foreach (table_invert($data_table) as $vals) {

        $peptide = $peptide_ptm = $vals['peptide'];
        $prob = $vals['psm_prob'];

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

    $file_name = "/var/jvlnse/local/optimize/" . $ui_file . ".csv";
    $csv->writeCsvFile($file_name, $prot_accounting_final);
}
$pgt->header("COMPLETE");

?>