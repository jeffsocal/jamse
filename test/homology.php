#!/usr/bin/env php
<?php
use ProxyIO\Cli;
use ProxySci\Omics\Proteomics\Peptide;
use ProxyTime\Timer;
use ProxyTime\ProgressTimer;

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 *
 * NAME
 * homology - find homologous peptides
 *
 * SYNOPSIS
 * homology -p <peptide_sequence>
 *
 * DESCRIPTION
 *
 * COMMAND LINE
 *
 * -p : <peptide_sequence>
 *
 * EXAMPLE
 *
 * homology -p TIPNLVNGLFK
 * homology -p QLICKSS[K42.01]AAFSTATMASYPHL[K188.03]
 *
 */
ini_set("memory_limit", "4096M");
set_include_path(__DIR__ . '/../');
require_once get_include_path() . 'src/Autoload/Functions.php';
require_once get_include_path() . 'vendor/autoload.php';

$cli = new Cli();
$amu = new Peptide();
$pgt = new ProgressTimer();

$peptide = $cli->getVar('p');
$jvln_lh = $cli->getVar('l', 10);
$jvln_i = $cli->getVar('db', 'yeast');

$jvln_s = '127.0.0.1';
$jvln_p = 55101;
$jvln_db = NULL;
$jvln_un = 'root';
$jvln_pw = '';

$peptide = preg_replace("/i/", "l", strtolower(getRootSequence($peptide)));

/*
 * SQL SERVER CONN
 */

$jvln_c = @new mysqli($jvln_s, $jvln_un, $jvln_pw, $jvln_db, $jvln_p);
$jvln_c_id = $jvln_c->thread_id;

if ($jvln_c->connect_errno) {
    $pgt->message(" - connect error: " . $jvln_c->connect_errno, $jvln_c->connect_error);
    $pgt->message(" - thread id", $jvln_c_id);
    $pgt->header("EARLY TERMINATION", TRUE);
}

/*
 * FOREACH SPECTRA -- START
 */
$jvln_q = "SELECT protein_name, protein_id, protein_desc
                        FROM peptides
                        WHERE
                        MATCH('" . $peptide . "') limit 100;";

$jvln_data = FALSE;
if ($jvln_r = $jvln_c->query($jvln_q)) {
    $jvln_data = $jvln_r->fetch_all(MYSQLI_ASSOC);
    $jvln_r->free();
}

$jvln_err = $jvln_c->error;
$jvln_c->close();

if (count($jvln_data) > 0)
    print_table(array_rowtocol($jvln_data), count($jvln_data));
?>