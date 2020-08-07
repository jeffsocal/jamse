<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
ini_set("memory_limit", "16384M");
use ProxyTime\Timer;
use ProxySci\Omics\Proteomics\PeptideModifications;
use ProxySci\Omics\Proteomics\Peptide;

$tmr = new Timer();
$mod = new PeptideModifications();
$amu = new Peptide();

$n = opts('n');
if (is_null($n))
    $n = 1;

$seq = opts('s');
if (is_null($seq))
    $seq = "AAACAIACHKNHKHTSSCKRVCSSTATRAAACAIACHKNHKHTSSCKRVCSSTATR";

// $seq = "IVRPEVDVMCTAFHDNEETFIK";
// $seq = "VTTDTWPRRAAQEPIIIIIIR";
// $seq = "AISQVIQIPIEVIQADSPCITIGEEYDKFIK";

$seq = "";

echo PHP_EOL;
print_message('peptide', $seq, '.', 30);

if (strlen($seq) <= 6)
    exit();

print_message('  mass', $amu->getMolecularWeight($seq), '.', 30);

$mc = max(0, count(preg_grep("/K|R/", str_split($seq))) - 1);
print_message('  mc', $mc, '.', 30);

$aas = $mod->getAminosThatHaveMods();

$aas = array_tostring($aas, '|', '');
preg_match_all('/' . $aas . '/', $seq, $wtf);
print_message('  mod-aa', count($wtf[0]), '.', 30);
print_message('  mods', $n, '.', 30);

$peps = $mod->getPTMPeptides($seq, $n);
print_message('  peps', count($peps), '.', 30);

// print_r($peps);

echo PHP_EOL;
print_message("TIMER", $tmr->timeinstr(TRUE), '.', 30);
print_message("MEMORY", memory_get_peak_usage() / (1e6) . " MB", ".", 30)?>
