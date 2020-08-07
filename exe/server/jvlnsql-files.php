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

// DELETE from jvln.validations where psm_pk in( select pk psm_pk FROM jvln.psms where scan_pk in (SELECT pk scan_pk from jvln.scans where file_pk > 55));
// DELETE FROM jvln.psms where scan_pk in (SELECT pk scan_pk from jvln.scans where file_pk > 55);
// DELETE FROM jvln.psm_alts where scan_pk in (SELECT pk scan_pk from jvln.scans where file_pk > 55);
// DELETE from jvln.scans where file_pk > 55;
// DELETE from jvln.experiments where file_pk > 55;
// DELETE from jvln.files where pk > 55;

$prg->header("COMPLETE");

?>